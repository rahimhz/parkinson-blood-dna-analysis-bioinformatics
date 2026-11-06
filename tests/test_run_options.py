"""Test Bash design-option handling without GEO downloads."""
import shlex
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def resolve_in_bash(requested_design, pair_column):
    helper = ROOT / 'scripts' / 'run_options.sh'
    command = (f'source {shlex.quote(str(helper))}; '
               f'resolve_design {shlex.quote(requested_design)} {shlex.quote(pair_column)}')
    return subprocess.run(['bash', '-c', command], capture_output=True, text=True, check=False)


class RunOptionTests(unittest.TestCase):
    def test_auto_default(self):
        result = resolve_in_bash('auto', '')
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), 'auto')

    def test_pair_option_selects_paired(self):
        result = resolve_in_bash('auto', 'Pair_ID')
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), 'paired')

    def test_pair_and_independent_conflict(self):
        result = resolve_in_bash('independent', 'Pair_ID')
        self.assertEqual(result.returncode, 2)
        self.assertIn('cannot be combined', result.stderr)

    def test_help(self):
        result = subprocess.run([str(ROOT / 'run.sh'), '-h'], cwd=ROOT,
                                capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0)
        self.assertIn('Usage:', result.stdout)

    def test_conflicting_cli_options_fail_before_download(self):
        result = subprocess.run([str(ROOT / 'run.sh'), '-p', 'Pair_ID', '-d', 'independent'],
                                cwd=ROOT, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn('cannot be combined', result.stderr)


if __name__ == '__main__':
    unittest.main()
