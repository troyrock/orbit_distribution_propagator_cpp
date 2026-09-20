"""Checkpoint integrity and real C++ runner workflow regression tests."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('coverage_grid', REPO / 'tools/coverage_grid.py')
grid = importlib.util.module_from_spec(spec)
spec.loader.exec_module(grid)
EXE = Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else None


class GridTests(unittest.TestCase):
    def test_axes_and_cartesian_product(self):
        self.assertEqual(grid.numbers('1,2,4'), [1, 2, 4])
        for text in ('2,1', '1,1', 'NaN', 'inf'):
            with self.assertRaises(argparse.ArgumentTypeError):
                grid.numbers(text)
        cases = grid.make_cases([[1000, 2000], [0, 90], [0], [1], [0.01]])
        self.assertEqual(len(cases), 4)
        self.assertEqual(len({row['run_id'] for row in cases}), 4)
        self.assertEqual(cases[-1]['inclination_deg'], 90)

    def test_checkpoint_requires_unique_matching_parameters(self):
        cases = grid.make_cases([[1000], [0], [0], [100], [1]])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            path = directory / 'batch_one.jsonl'
            path.write_text(json.dumps(cases[0]) + '\n')
            self.assertEqual(len(grid.load_completed(directory, cases)), 1)
            (directory / 'batch_two.jsonl').write_text(path.read_text())
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                grid.load_completed(directory, cases)
            (directory / 'batch_two.jsonl').unlink()
            changed = dict(cases[0], position_sigma_km=50)
            path.write_text(json.dumps(changed) + '\n')
            with self.assertRaisesRegex(ValueError, 'parameters'):
                grid.load_completed(directory, cases)

    @unittest.skipIf(EXE is None, 'Executable argument required for end-to-end test')
    def test_real_masked_case_resume_and_manifest_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            config = directory / 'base.cfg'
            config.write_text((REPO / 'examples/coverage_grid.cfg').read_text())
            output = directory / 'study'
            command = [sys.executable, str(REPO / 'tools/coverage_grid.py'), '--exe', str(EXE),
                       '--config', str(config), '--output', str(output), '--jobs', '1',
                       '--threads', '1', '--altitudes', '1000', '--inclinations', '0',
                       '--eccentricities', '0', '--positions', '100', '--velocities', '1']
            first = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            study = json.loads((output / 'study.json').read_text())
            self.assertTrue(study['metadata']['complete'])
            self.assertEqual(study['runs'][0]['samples'], 20000)
            self.assertEqual(study['runs'][0]['status'], 'earth_intersection')
            self.assertIsNone(study['runs'][0]['coverage_time_days'])
            raw = list((output / 'raw').glob('*.jsonl'))
            before = raw[0].read_bytes()
            with grid.output_lock(output):
                competing = subprocess.run(command, capture_output=True, text=True)
                self.assertNotEqual(competing.returncode, 0)
                self.assertIn('locked by another', competing.stderr)
            resumed = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            self.assertIn('1 already complete', resumed.stdout)
            self.assertEqual(len(list((output / 'raw').glob('*.jsonl'))), 1)
            self.assertEqual(raw[0].read_bytes(), before)
            config.write_text(config.read_text() + '\n# Configuration provenance changed\n')
            changed = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn('Resume rejected', changed.stderr)


if __name__ == '__main__':
    unittest.main()
