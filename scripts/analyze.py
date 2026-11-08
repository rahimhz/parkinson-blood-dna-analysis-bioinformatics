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


def validate_numeric_matrix(matrix_path, chunksize=2000):
    """Check that the matrix contains numeric measurements without loading it all at once."""
    has_finite_value = False
    with pd.read_csv(
        matrix_path, compression='infer', index_col=0,
        low_memory=False, chunksize=chunksize
    ) as reader:
        for chunk in reader:
            numeric = chunk.apply(pd.to_numeric, errors='coerce')
            if np.isfinite(numeric.to_numpy(dtype=float)).any():
                has_finite_value = True
                break
    if not has_finite_value:
        raise ValueError('Matrix contains no finite numeric M-value measurements.')


def finite_row_mean(values):
    """Calculate row means without warnings for rows containing only missing values."""
    finite = np.isfinite(values)
    counts = finite.sum(axis=1)
    totals = np.where(finite, values, 0.0).sum(axis=1)
    return np.divide(
        totals, counts,
        out=np.full(values.shape[0], np.nan, dtype=float),
        where=counts > 0
    )


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
        pd_mean = finite_row_mean(pd_safe)
        control_mean = finite_row_mean(control_safe)
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
        pd_mean = finite_row_mean(pd_array)
        control_mean = finite_row_mean(control_array)
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


def calculate_adjusted_results(matrix, labeled, sample_column):
    """Fit vectorized per-feature OLS: M-value ~ PD status + age + sex."""
    required = {'Age', 'Sex'}
    missing_columns = required - set(labeled.columns)
    if missing_columns:
        raise ValueError(f'Adjusted analysis requires metadata columns: {sorted(missing_columns)}')

    age = pd.to_numeric(labeled['Age'], errors='coerce').to_numpy(float)
    if not np.isfinite(age).all():
        raise ValueError('Age contains missing or non-numeric values; adjusted analysis cannot proceed.')

    sex_text = labeled['Sex'].astype(str).str.strip().str.lower()
    sex_map = {'m': 1.0, 'male': 1.0, 'f': 0.0, 'female': 0.0}
    sex = sex_text.map(sex_map)
    if sex.isna().any() or sex.nunique() != 2:
        raise ValueError('Sex must contain both recognized categories (M/F or Male/Female) with no missing values.')

    group = (labeled['_group'] == 'PD').to_numpy(dtype=float)
    n = len(labeled)
    design_matrix = np.column_stack([
        np.ones(n),
        group,
        age - age.mean(),
        sex.to_numpy(dtype=float),
    ])
    rank = int(np.linalg.matrix_rank(design_matrix))
    if rank != design_matrix.shape[1]:
        raise ValueError('Adjusted model is rank-deficient: PD status, age, and sex cannot be separated in this sample.')

    residual_df = n - rank
    if residual_df <= 0:
        raise ValueError('Adjusted model has no residual degrees of freedom.')

    values = matrix.apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
    complete = np.isfinite(values).all(axis=1)
    if not complete.any():
        raise ValueError('No methylation features have complete values across all samples for adjusted OLS.')
    y = values[complete].T
    coefficients = np.linalg.lstsq(design_matrix, y, rcond=None)[0]
    residuals = y - design_matrix @ coefficients
    residual_variance = np.sum(residuals ** 2, axis=0) / residual_df
    covariance_base = np.linalg.inv(design_matrix.T @ design_matrix)
    standard_error = np.sqrt(np.maximum(residual_variance * covariance_base[1, 1], 0.0))
    effect = coefficients[1]
    statistic = np.divide(
        effect, standard_error,
        out=np.full_like(effect, np.nan),
        where=standard_error > 0
    )
    from scipy.stats import t as student_t
    pvalue = 2.0 * student_t.sf(np.abs(statistic), df=residual_df)
    valid = np.isfinite(pvalue) & np.isfinite(effect)
    if not valid.any():
        raise ValueError('Adjusted OLS produced no valid feature-level p-values.')

    pd_ids = labeled.loc[labeled['_group'] == 'PD', sample_column].tolist()
    control_ids = labeled.loc[labeled['_group'] == 'Control', sample_column].tolist()
    pd_values = values[:, [matrix.columns.get_loc(sample) for sample in pd_ids]]
    control_values = values[:, [matrix.columns.get_loc(sample) for sample in control_ids]]
    mean_pd = finite_row_mean(pd_values)
    mean_control = finite_row_mean(control_values)
    feature_ids = matrix.index.astype(str).to_numpy()[complete][valid]
    result = pd.DataFrame({
        'feature_id': feature_ids,
        'n_PD': len(pd_ids),
        'n_Control': len(control_ids),
        'mean_PD': mean_pd[complete][valid],
        'mean_Control': mean_control[complete][valid],
        'mean_difference_PD_minus_Control': mean_pd[complete][valid] - mean_control[complete][valid],
        'adjusted_PD_coefficient': effect[valid],
        'test_statistic': statistic[valid],
        'p_value': pvalue[valid],
    })
    result['q_value_BH'] = multipletests(result['p_value'].to_numpy(), method='fdr_bh')[1]
    summary = {
        'model': 'M-value ~ PD_status + age + sex',
        'n_samples': n,
        'n_PD': int(group.sum()),
        'n_Control': int(n - group.sum()),
        'age_mean': float(age.mean()),
        'age_sd': float(age.std(ddof=1)),
        'sex_counts': {str(k): int(v) for k, v in sex_text.value_counts().to_dict().items()},
        'features_in_matrix': int(len(values)),
        'features_tested_complete_cases': int(complete.sum()),
        'features_excluded_for_missing_or_nonfinite_values': int((~complete).sum()),
        'residual_degrees_of_freedom': int(residual_df),
        'model_rank': rank,
        'interpretation': 'Adjusted PD coefficient is the estimated PD-Control difference in M-values conditional on age and sex.',
    }
    return result.sort_values('p_value'), 'OLS t-test for PD coefficient adjusted for age and sex', summary


def write_json(path, payload):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')


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
    # Validation-only mode scans the matrix in chunks to avoid a second full-size
    # in-memory copy. Statistical modes load the selected samples for testing.
    if args.design == 'auto':
        validate_numeric_matrix(args.matrix)
        matrix = None
    else:
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

    if args.design == 'adjusted':
        results, test_name, model_summary = calculate_adjusted_results(matrix, labeled, sample_column)
    else:
        results, test_name = calculate_results(matrix, labeled, sample_column, args.design, args.pair_column)
        model_summary = {}
    validation.update({
        'design': args.design,
        'pair_column': args.pair_column,
        'model_summary': model_summary,
        'statistical_test': test_name,
        'multiple_testing_correction': 'Benjamini-Hochberg FDR',
        'features_tested': int(len(results)),
        'limitations': [
            'Adjusted mode includes age and sex; independent and paired t-test modes do not model covariates.',
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
    parser.add_argument('--design', choices=['auto', 'paired', 'independent', 'adjusted'], default='auto')
    parser.add_argument('--pair-column')
    return parser


def main():
    args = build_parser().parse_args()
    if args.design == 'paired' and not args.pair_column:
        print('ERROR: paired design requires --pair-column.', file=sys.stderr)
        return 2
    if args.design in {'independent', 'adjusted'} and args.pair_column:
        print('ERROR: --pair-column cannot be combined with independent or adjusted design.', file=sys.stderr)
        return 2
    try:
        run_analysis(args)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
