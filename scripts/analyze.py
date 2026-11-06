#!/usr/bin/env python3
"""Validate GSE165081 metadata/matrix; optionally test processed M-values."""

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import multipletests


def canonical_group(value):
    """Map recognized PD and Control/CN labels; otherwise return None."""
    if pd.isna(value):
        return None
    text = re.sub(r'\s+', ' ', str(value).strip().lower()).replace('_', ' ').replace('-', ' ')
    if text in {"pd", 'parkinson', 'parkinson disease', "parkinson's disease", 'parkinsons disease'}:
        return 'PD'
    if text in {'cn', 'control', 'healthy control', 'control normal', 'normal control'}:
        return 'Control'
    return None


def read_csv_header(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        try:
            header = next(csv.reader(stream))
        except StopIteration as error:
            raise ValueError(f'CSV file is empty: {path}') from error
    if len(header) < 2:
        raise ValueError(f'Expected an identifier and data columns in {path}.')
    duplicates = sorted(name for name, count in Counter(header).items() if count > 1)
    if duplicates:
        raise ValueError(f'Duplicate CSV headers in {path}: {duplicates}')
    return header


def read_first_column(path):
    frame = pd.read_csv(path, compression='infer', usecols=[0], dtype=str, keep_default_na=False)
    if frame.empty:
        raise ValueError(f'CSV contains no data rows: {path}')
    return frame.iloc[:, 0].astype(str).str.strip().tolist()


def matrix_candidate_ids(matrix_path):
    """Inspect sample-axis candidates without presuming the CSV orientation."""
    header = read_csv_header(matrix_path)
    first_column = read_first_column(matrix_path)
    column_ids = [str(value).strip() for value in header[1:]]
    column_ids = [value for value in column_ids if value]
    row_ids = [value for value in first_column if value]
    return header, column_ids, row_ids


def find_sample_column(samples, matrix_ids, requested=None):
    """Find metadata column by overlap with IDs actually found in the matrix."""
    if requested:
        if requested not in samples.columns:
            raise ValueError(f'Sample-ID column {requested!r} not found. Headers: {list(samples.columns)}')
        return requested
    matrix_ids = set(matrix_ids)
    candidates = []
    scores = {}
    for column in samples.columns:
        values = samples[column].astype(str).str.strip()
        nonempty = values[values != '']
        overlap = set(nonempty) & matrix_ids
        if overlap:
            scores[column] = len(overlap)
    if not scores:
        raise ValueError('No metadata column values match IDs in matrix headers or first column. Inspect files.')
    best_score = max(scores.values())
    candidates = [column for column, score in scores.items() if score == best_score]
    if len(candidates) != 1:
        raise ValueError(f'Sample-ID column detection is ambiguous: {candidates}. Use --sample-column with actual header.')
    return candidates[0]


def detect_group_column(samples, requested=None):
    if requested:
        if requested not in samples.columns:
            raise ValueError(f'Group column {requested!r} not found. Headers: {list(samples.columns)}')
        return requested
    candidates = [column for column in samples.columns
                  if set(samples[column].map(canonical_group).dropna()) == {'PD', 'Control'}]
    if len(candidates) != 1:
        raise ValueError(f'Could not uniquely detect PD/Control column; candidates={candidates}. Use --group-column.')
    return candidates[0]


def labeled_metadata(samples, sample_column, group_column):
    result = samples.copy()
    result[sample_column] = result[sample_column].astype(str).str.strip()
    labels = result[group_column].astype(str).str.strip()
    result['_group'] = labels.map(canonical_group)
    bad = result['_group'].isna()
    if bad.any():
        raise ValueError(f'Missing or unrecognized group labels: {sorted(set(labels[bad]))}')
    return result


def detect_matrix_orientation(matrix_path, labeled_ids):
    header, column_ids, row_ids = matrix_candidate_ids(matrix_path)
    labeled_ids = [str(value).strip() for value in labeled_ids]
    col_counts, row_counts = Counter(column_ids), Counter(row_ids)
    col_matches = sum(col_counts.get(value, 0) == 1 for value in labeled_ids)
    row_matches = sum(row_counts.get(value, 0) == 1 for value in labeled_ids)
    if col_matches == len(labeled_ids) and row_matches != len(labeled_ids):
        orientation, axis_ids = 'samples_in_columns', column_ids
    elif row_matches == len(labeled_ids) and col_matches != len(labeled_ids):
        orientation, axis_ids = 'samples_in_rows', row_ids
    else:
        raise ValueError(f'Cannot unambiguously locate all labeled IDs: total={len(labeled_ids)}, column matches={col_matches}, row matches={row_matches}.')
    duplicates = sorted(value for value, count in Counter(axis_ids).items() if count > 1)
    if duplicates:
        raise ValueError(f'Duplicate IDs on matrix sample axis: {duplicates[:20]}')
    return orientation, axis_ids, header


def validate_sample_matches(labeled, sample_column, axis_ids):
    metadata_counts = Counter(labeled[sample_column].astype(str).str.strip())
    matrix_counts = Counter(map(str, axis_ids))
    duplicate_metadata = sorted(value for value, count in metadata_counts.items() if count > 1)
    duplicate_matrix = sorted(value for value, count in matrix_counts.items() if count > 1)
    unmatched = sorted(value for value in metadata_counts if matrix_counts.get(value, 0) != 1)
    extras = sorted(value for value in matrix_counts if value not in metadata_counts)
    print(f'Duplicate metadata IDs: {duplicate_metadata}')
    print(f'Duplicate matrix IDs: {duplicate_matrix}')
    print(f'Labeled IDs absent/non-unique in matrix: {unmatched}')
    print(f'Matrix IDs without labeled metadata: {extras[:20]}')
    if duplicate_metadata or duplicate_matrix or unmatched:
        raise ValueError('Sample validation failed: every labeled PD/Control sample must match exactly one matrix sample.')
    counts = Counter(labeled['_group'])
    print(f'Final validated group counts: PD={counts["PD"]}, Control={counts["Control"]}')
    return dict(counts), extras


def load_matrix(matrix_path, orientation, sample_ids):
    matrix = pd.read_csv(matrix_path, compression='infer', index_col=0, low_memory=False)
    if orientation == 'samples_in_rows':
        matrix = matrix.T
    elif orientation != 'samples_in_columns':
        raise ValueError(f'Unknown orientation: {orientation}')
    matrix.columns = [str(value).strip() for value in matrix.columns]
    matrix.index = [str(value).strip() for value in matrix.index]
    if not matrix.columns.is_unique or not matrix.index.is_unique:
        raise ValueError('Loaded matrix contains duplicate sample or feature IDs.')
    missing = sorted(set(map(str, sample_ids)) - set(matrix.columns))
    if missing:
        raise ValueError(f'Samples absent after loading matrix: {missing[:20]}')
    return matrix.loc[:, list(map(str, sample_ids))]


def validate_pairs(labeled, sample_column, pair_column):
    if pair_column not in labeled.columns:
        raise ValueError(f'Pair-ID column {pair_column!r} not found.')
    pairs = labeled[[sample_column, '_group', pair_column]].copy()
    pairs[pair_column] = pairs[pair_column].astype(str).str.strip()
    missing = pairs[pair_column].eq('') | pairs[pair_column].str.lower().isin({'nan', 'none', 'na'})
    if missing.any():
        raise ValueError(f'Missing pair IDs: {pairs.loc[missing, sample_column].tolist()[:20]}')
    invalid = []
    for pair_id, rows in pairs.groupby(pair_column, sort=False):
        if len(rows) != 2 or Counter(rows['_group']) != Counter(['PD', 'Control']):
            invalid.append(pair_id)
    if invalid:
        raise ValueError(f'Each pair must contain exactly one PD and one Control; invalid IDs={invalid[:20]}')
    pd_map = pairs[pairs['_group'] == 'PD'].set_index(pair_column)[sample_column].to_dict()
    control_map = pairs[pairs['_group'] == 'Control'].set_index(pair_column)[sample_column].to_dict()
    pair_ids = sorted(pd_map)
    return pair_ids, pd_map, control_map


def calculate_results(matrix, labeled, sample_column, design, pair_column=None):
    """Vectorized per-feature tests; no Python loop over matrix features."""
    counts = Counter(labeled['_group'])
    if counts['PD'] < 2 or counts['Control'] < 2:
        raise ValueError(f'At least two samples per group are required: {dict(counts)}')
    if design == 'paired':
        if not pair_column:
            raise ValueError('Paired analysis requires --pair-column.')
        pairs, pd_map, control_map = validate_pairs(labeled, sample_column, pair_column)
        pd_ids = [pd_map[pair] for pair in pairs]
        control_ids = [control_map[pair] for pair in pairs]
        pd_array = matrix[pd_ids].apply(pd.to_numeric, errors='coerce').to_numpy(float)
        control_array = matrix[control_ids].apply(pd.to_numeric, errors='coerce').to_numpy(float)
        valid = np.isfinite(pd_array) & np.isfinite(control_array)
        n = valid.sum(axis=1)
        pd_safe = np.where(valid, pd_array, np.nan)
        control_safe = np.where(valid, control_array, np.nan)
        pd_mean = np.nanmean(pd_safe, axis=1)
        control_mean = np.nanmean(control_safe, axis=1)
        test = ttest_rel(pd_safe, control_safe, axis=1, nan_policy='omit')
        n_pd = n_control = n
        test_name = 'paired t-test'
    elif design == 'independent':
        pd_ids = labeled.loc[labeled['_group'] == 'PD', sample_column].tolist()
        control_ids = labeled.loc[labeled['_group'] == 'Control', sample_column].tolist()
        pd_array = matrix[pd_ids].apply(pd.to_numeric, errors='coerce').to_numpy(float)
        control_array = matrix[control_ids].apply(pd.to_numeric, errors='coerce').to_numpy(float)
        n_pd = np.isfinite(pd_array).sum(axis=1)
        n_control = np.isfinite(control_array).sum(axis=1)
        pd_mean = np.nanmean(pd_array, axis=1)
        control_mean = np.nanmean(control_array, axis=1)
        test = ttest_ind(pd_array, control_array, axis=1, equal_var=False, nan_policy='omit')
        test_name = 'Welch two-sample t-test (unpaired; matching not modeled)'
    else:
        raise ValueError('Auto mode must not call calculate_results.')

    statistic = np.asarray(test.statistic, dtype=float)
    pvalue = np.asarray(test.pvalue, dtype=float)
    valid_feature = (n_pd >= 2) & (n_control >= 2) & np.isfinite(pvalue)
    if not valid_feature.any():
        raise ValueError('No features have enough non-missing values to test.')
    result = pd.DataFrame({
        'feature_id': matrix.index.astype(str)[valid_feature],
        'n_PD': n_pd[valid_feature],
        'n_Control': n_control[valid_feature],
        'mean_PD': pd_mean[valid_feature],
        'mean_Control': control_mean[valid_feature],
        'mean_difference_PD_minus_Control': pd_mean[valid_feature] - control_mean[valid_feature],
        'test_statistic': statistic[valid_feature],
        'p_value': pvalue[valid_feature],
    })
    result['q_value_BH'] = multipletests(result['p_value'].to_numpy(), method='fdr_bh')[1]
    return result.sort_values('p_value'), test_name


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')


def write_statistical_outputs(results, summary, output_dir):
    out = Path(output_dir)
    figure_dir = out / 'figures'
    figure_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(out / 'differential_methylation.csv', index=False)
    write_json(out / 'analysis_summary.json', summary)
    p = results['p_value'].clip(lower=np.finfo(float).tiny)
    plt.figure(figsize=(7, 5))
    plt.scatter(results['mean_difference_PD_minus_Control'], -np.log10(p), s=6, alpha=0.5)
    plt.xlabel('Mean M-value difference (PD − Control)')
    plt.ylabel('−log10(p-value)')
    plt.title('GSE165081: PD versus Control')
    plt.tight_layout()
    plt.savefig(figure_dir / 'methylation_volcano.png', dpi=150)
    plt.close()


def run_analysis(args):
    samples = pd.read_csv(args.samples, dtype=str, keep_default_na=False)
    if samples.empty or list(samples.columns) != read_csv_header(args.samples):
        raise ValueError('Sample sheet is empty or parsed headers differ from the actual CSV.')
    # Detect sample field by intersection with actual matrix header/first-column IDs.
    _, column_candidates, row_candidates = matrix_candidate_ids(args.matrix)
    sample_column = find_sample_column(samples, column_candidates + row_candidates, args.sample_column)
    group_column = detect_group_column(samples, args.group_column)
    labeled = labeled_metadata(samples, sample_column, group_column)
    ids = labeled[sample_column].astype(str).str.strip().tolist()
    orientation, axis_ids, _ = detect_matrix_orientation(args.matrix, ids)
    counts, extras = validate_sample_matches(labeled, sample_column, axis_ids)
    # Load and verify numeric matrix even in validation-only auto mode.
    matrix = load_matrix(args.matrix, orientation, ids)
    numeric = matrix.apply(pd.to_numeric, errors='coerce')
    if numeric.empty or not np.isfinite(numeric.to_numpy(dtype=float)).any():
        raise ValueError('Matrix contains no finite numeric M-value measurements.')

    validation = {
        'accession': 'GSE165081',
        'matrix_file': str(args.matrix),
        'sample_sheet': str(args.samples),
        'matrix_orientation': orientation,
        'sample_id_column': sample_column,
        'group_column': group_column,
        'validated_group_counts': counts,
        'validated_sample_ids': ids,
        'unlabeled_matrix_sample_ids': extras,
        'validation_status': 'passed',
        'design_requested': args.design,
        'statistical_analysis': 'not run' if args.design == 'auto' else 'requested',
    }
    print(f'Matrix orientation: {orientation}')
    print(f'Sample-ID column: {sample_column}')
    print(f'Group column: {group_column}')
    print(f'Final validated group counts: {counts}')
    if args.design == 'auto':
        write_json(Path(args.output_dir) / 'validation_summary.json', validation)
        print('Validation completed successfully. No statistical tests were run.')
        return

    results, test_name = calculate_results(matrix, labeled, sample_column, args.design, args.pair_column)
    validation.update({
        'design': args.design,
        'pair_column': args.pair_column,
        'statistical_test': test_name,
        'multiple_testing_correction': 'Benjamini-Hochberg FDR',
        'features_tested': int(len(results)),
        'limitations': [
            'No covariates are modeled.',
            'Associations do not establish causality or diagnostic performance.',
            'Processed M-values are analyzed; raw IDAT files are not processed.',
        ],
    })
    write_statistical_outputs(results, validation, args.output_dir)


def build_parser():
    parser = argparse.ArgumentParser(description='Validate or analyze processed GSE165081 M-values.')
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--sample-column')
    parser.add_argument('--group-column')
    parser.add_argument('--design', choices=['auto', 'paired', 'independent'], default='auto')
    parser.add_argument('--pair-column')
    return parser


def main():
    args = build_parser().parse_args()
    if args.design == 'paired' and not args.pair_column:
        print('ERROR: paired design requires --pair-column.', file=sys.stderr)
        return 2
    if args.design == 'independent' and args.pair_column:
        print('ERROR: --pair-column cannot be used with independent design.', file=sys.stderr)
        return 2
    try:
        run_analysis(args)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
