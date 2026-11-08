"""Offline unit tests with synthetic data; never downloads GEO files."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
import numpy as np
from scripts.analyze import (
    calculate_adjusted_results, calculate_results, detect_group_column, detect_matrix_orientation,
    find_sample_column, labeled_metadata, load_matrix, run_analysis,
    validate_numeric_matrix,
    validate_pairs, validate_sample_matches,
)


class AnalysisTests(unittest.TestCase):
    def write_csv(self, directory, name, content):
        path = Path(directory) / name
        path.write_text(content, encoding='utf-8')
        return str(path)

    def test_sample_and_group_column_detection_by_matrix_overlap(self):
        samples = pd.DataFrame({'odd_header': ['S1', 'S2'], 'condition': ['PD', 'Control']})
        self.assertEqual(find_sample_column(samples, ['S1', 'S2']), 'odd_header')
        self.assertEqual(detect_group_column(samples), 'condition')

    def test_missing_and_invalid_group_labels_fail(self):
        for label in ['', 'unknown']:
            samples = pd.DataFrame({'sample': ['S1', 'S2'], 'group': ['PD', label]})
            with self.assertRaisesRegex(ValueError, 'Missing or unrecognized'):
                labeled_metadata(samples, 'sample', 'group')

    def test_duplicate_and_unmatched_ids_fail(self):
        duplicate = pd.DataFrame({'sample': ['S1', 'S1'], '_group': ['PD', 'Control']})
        with self.assertRaisesRegex(ValueError, 'Sample validation failed'):
            validate_sample_matches(duplicate, 'sample', ['S1', 'S2'])
        unmatched = pd.DataFrame({'sample': ['S1', 'S3'], '_group': ['PD', 'Control']})
        with self.assertRaisesRegex(ValueError, 'Sample validation failed'):
            validate_sample_matches(unmatched, 'sample', ['S1', 'S2'])

    def test_load_columns_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_csv(directory, 'cols.csv', 'feature,S1,S2\ncgA,1,2\ncgB,3,4\n')
            orientation, axis, _ = detect_matrix_orientation(path, ['S1', 'S2'])
            matrix = load_matrix(path, orientation, ['S1', 'S2'])
            self.assertEqual(orientation, 'samples_in_columns')
            self.assertEqual(axis, ['S1', 'S2'])
            self.assertEqual(matrix.loc['cgA', 'S2'], 2)

    def test_load_rows_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_csv(directory, 'rows.csv', 'sample,cgA,cgB\nS1,1,3\nS2,2,4\n')
            orientation, axis, _ = detect_matrix_orientation(path, ['S1', 'S2'])
            matrix = load_matrix(path, orientation, ['S1', 'S2'])
            self.assertEqual(orientation, 'samples_in_rows')
            self.assertEqual(axis, ['S1', 'S2'])
            self.assertEqual(matrix.loc['cgB', 'S1'], 3)

    def test_paired_analysis_synthetic(self):
        samples = pd.DataFrame({
            'sample': ['P1', 'C1', 'P2', 'C2'],
            '_group': ['PD', 'Control', 'PD', 'Control'],
            'pair': ['A', 'A', 'B', 'B'],
        })
        pairs, pd_map, control_map = validate_pairs(samples, 'sample', 'pair')
        self.assertEqual(pairs, ['A', 'B'])
        matrix = pd.DataFrame({pd_map['A']: [2, 5], control_map['A']: [1, 4],
                               pd_map['B']: [3, 7], control_map['B']: [1, 4]},
                              index=['cgA', 'cgB'])
        result, name = calculate_results(matrix, samples, 'sample', 'paired', 'pair')
        self.assertEqual(name, 'paired t-test')
        self.assertEqual(len(result), 2)
        self.assertTrue(result['q_value_BH'].notna().all())

    def test_invalid_pair_fails(self):
        samples = pd.DataFrame({'sample': ['P1', 'C1', 'P2'],
                                '_group': ['PD', 'Control', 'PD'],
                                'pair': ['A', 'A', 'B']})
        with self.assertRaisesRegex(ValueError, 'Each pair must contain'):
            validate_pairs(samples, 'sample', 'pair')

    def test_chunked_numeric_validation_accepts_numeric_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_csv(directory, 'numeric.csv',
                'feature,S1,S2\ncgA,1.2,2.3\ncgB,3.4,4.5\n')
            validate_numeric_matrix(path, chunksize=1)

    def test_chunked_numeric_validation_rejects_non_numeric_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_csv(directory, 'text.csv',
                'feature,S1,S2\ncgA,missing,unknown\ncgB,none,NA\n')
            with self.assertRaisesRegex(ValueError, 'no finite numeric'):
                validate_numeric_matrix(path, chunksize=1)

    def test_age_sex_adjusted_ols(self):
        samples = pd.DataFrame({
            'sample': [f'S{i}' for i in range(8)],
            '_group': ['PD', 'PD', 'PD', 'PD', 'Control', 'Control', 'Control', 'Control'],
            'Age': [60, 62, 64, 66, 61, 63, 65, 67],
            'Sex': ['M', 'F', 'M', 'F', 'M', 'F', 'M', 'F'],
        })
        age = samples['Age'].to_numpy(float)
        pd_flag = (samples['_group'] == 'PD').to_numpy(float)
        sex_flag = (samples['Sex'] == 'M').to_numpy(float)
        base = 0.15 * (age - age.mean()) + 0.4 * sex_flag
        noise = np.array([0.12, -0.08, 0.04, -0.03, -0.05, 0.09, -0.11, 0.02])
        values = pd_flag * 1.8 + base + noise
        matrix = pd.DataFrame([values, values * 0.5 + np.array([0.1, -0.1, 0.0, 0.05, -0.04, 0.08, -0.02, 0.03])],
                              index=['cg_effect', 'cg_second'],
                              columns=samples['sample'])
        results, name, summary = calculate_adjusted_results(matrix, samples, 'sample')
        self.assertIn('adjusted for age and sex', name)
        self.assertEqual(summary['n_samples'], 8)
        self.assertEqual(summary['features_tested_complete_cases'], 2)
        self.assertIn('adjusted_PD_coefficient', results.columns)
        self.assertTrue(results['q_value_BH'].between(0, 1).all())

    def test_adjusted_ols_requires_age_and_sex(self):
        samples = pd.DataFrame({
            'sample': ['P1', 'P2', 'C1', 'C2'],
            '_group': ['PD', 'PD', 'Control', 'Control'],
        })
        matrix = pd.DataFrame([[1.0, 2.0, 1.5, 1.7]], columns=samples['sample'])
        with self.assertRaisesRegex(ValueError, 'requires metadata columns'):
            calculate_adjusted_results(matrix, samples, 'sample')

    def test_auto_validation_does_not_create_statistical_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            matrix_path = self.write_csv(root, 'matrix.csv',
                'feature,S1,S2\ncgA,1,2\ncgB,3,4\n')
            sheet_path = self.write_csv(root, 'samples.csv',
                'unusual_sample_header,condition\nS1,PD\nS2,Control\n')
            out = root / 'out'
            from argparse import Namespace
            args = Namespace(matrix=matrix_path, samples=sheet_path,
                             output_dir=str(out), sample_column=None,
                             group_column=None, design='auto', pair_column=None)
            with patch('scripts.analyze.load_matrix', side_effect=AssertionError(
                    'Validation-only mode should not load the full matrix')):
                run_analysis(args)
            self.assertTrue((out / 'validation_summary.json').is_file())
            self.assertFalse((out / 'differential_methylation.csv').exists())
            self.assertFalse((out / 'analysis_summary.json').exists())


if __name__ == '__main__':
    unittest.main()
