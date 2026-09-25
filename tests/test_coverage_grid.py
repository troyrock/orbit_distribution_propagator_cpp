"""Checkpoint integrity and real C++ runner workflow regression tests."""
import argparse
import copy
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
    def reuse_fixture(self, directory, altitudes):
        axes = [altitudes, [0], [0], [100], [1]]
        cases = grid.make_cases(axes)
        metadata = dict(zip(grid.AXES, axes), samples=20000, seed=20260919,
                        coverage_definition='test coverage contract')
        manifest = {'schema_version': 1, 'executable_sha256': 'exe-sha',
                    'config_sha256': 'config-sha', 'config_text': 'test config',
                    'metadata': metadata}
        records = [dict(case, status='earth_intersection', samples=20000, seed=20260919,
                        earth_intersections=1, particles_below_300km=1,
                        particles_below_500km=1, particles_below_1000km=1,
                        min_initial_perigee_km=-1, min_mean_perigee_km=None,
                        evaluated_epochs=0, cadence_days=0, evaluated_horizon_days=0,
                        coverage_time_days=None, coverage_lower_days=None,
                        coverage_upper_days=None, coverage_confirmation_days=None,
                        onset_max_gap_deg=None) for case in cases]
        directory.mkdir()
        grid.atomic_json(directory / 'manifest.json', manifest)
        grid.atomic_json(directory / 'study.json', {
            'metadata': dict(metadata, complete=True, expected_cases=len(cases),
                             completed_cases=len(cases)), 'runs': records})
        return manifest, cases

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

    def test_reuse_reindexes_parameters_and_resumes_without_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source'
            manifest, _ = self.reuse_fixture(source, [2000, 3000])
            original = (source / 'study.json').read_bytes()
            destination = root / 'extended'
            raw = destination / 'raw'
            raw.mkdir(parents=True)
            manifest['metadata']['perigee_altitudes_km'] = [1000, 2000, 3000]
            manifest['metadata']['altitude_interpolation'] = 'log_geocentric_perigee_radius'
            cases = grid.make_cases([manifest['metadata'][name] for name in grid.AXES])
            expected_manifest = copy.deepcopy(manifest)
            saved = grid.initialize_manifest(destination, raw, manifest, cases, source)
            completed = grid.load_completed(raw, cases)
            self.assertEqual(list(completed), ['case_000001', 'case_000002'])
            self.assertEqual(completed['case_000001']['perigee_altitude_km'], 2000)
            self.assertEqual(saved['reuse_study']['adopted_cases'], 2)
            self.assertEqual(saved['reuse_study']['study_sha256'], grid.digest(source / 'study.json'))
            self.assertEqual((source / 'study.json').read_bytes(), original)
            checkpoint = raw / 'batch_reused.jsonl'
            checkpoint_before = checkpoint.read_bytes()
            source.rename(root / 'source_moved')
            resumed = grid.initialize_manifest(destination, raw, expected_manifest, cases, None)
            self.assertEqual(resumed, saved)
            self.assertEqual(checkpoint.read_bytes(), checkpoint_before)
            self.assertEqual(len(list(raw.glob('batch_*.jsonl'))), 1)

    def test_reuse_rejects_hash_mismatch_and_inconsistent_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source'
            manifest, cases = self.reuse_fixture(source, [1000])
            for key in ('executable_sha256', 'config_sha256', 'config_text'):
                changed = copy.deepcopy(manifest)
                changed[key] = 'different'
                with self.subTest(key=key), self.assertRaisesRegex(ValueError, key):
                    grid.prepare_reuse(source, root / 'target', changed, cases)
            changed = copy.deepcopy(manifest)
            changed['metadata']['seed'] += 1
            with self.assertRaisesRegex(ValueError, 'metadata seed'):
                grid.prepare_reuse(source, root / 'target', changed, cases)
            with self.assertRaisesRegex(ValueError, 'source must differ'):
                grid.prepare_reuse(source, source, manifest, cases)
            study = json.loads((source / 'study.json').read_text())
            study['metadata']['coverage_definition'] = 'changed'
            grid.atomic_json(source / 'study.json', study)
            with self.assertRaisesRegex(ValueError, 'coverage_definition'):
                grid.prepare_reuse(source, root / 'target', manifest, cases)

    def test_reuse_rejects_duplicate_parameters_and_incomplete_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source'
            manifest, cases = self.reuse_fixture(source, [1000, 2000])
            study = json.loads((source / 'study.json').read_text())
            study['runs'][1] = dict(study['runs'][0], run_id='different_id')
            grid.atomic_json(source / 'study.json', study)
            with self.assertRaisesRegex(ValueError, 'duplicate parameter'):
                grid.prepare_reuse(source, root / 'target', manifest, cases)
            study['metadata']['complete'] = False
            grid.atomic_json(source / 'study.json', study)
            with self.assertRaisesRegex(ValueError, 'incomplete study'):
                grid.prepare_reuse(source, root / 'target', manifest, cases)

    def test_reuse_recovers_interrupted_adoption_and_detects_changed_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source'
            manifest, cases = self.reuse_fixture(source, [1000])
            destination = root / 'target'
            raw = destination / 'raw'
            raw.mkdir(parents=True)
            initial = copy.deepcopy(manifest)
            provenance, _ = grid.prepare_reuse(source, destination, manifest, cases)
            grid.atomic_json(destination / 'manifest.json', dict(manifest, reuse_study=provenance))
            grid.initialize_manifest(destination, raw, manifest, cases, None)
            self.assertEqual(len(grid.load_completed(raw, cases)), 1)
            study = json.loads((source / 'study.json').read_text())
            study['runs'][0]['min_initial_perigee_km'] = -2
            grid.atomic_json(source / 'study.json', study)
            with self.assertRaisesRegex(ValueError, 'source content differs'):
                grid.initialize_manifest(destination, raw, initial, cases, source)

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
            self.assertNotIn('altitude_interpolation', study['metadata'])
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
            extended = directory / 'extended'
            extension_command = command.copy()
            extension_command[extension_command.index('--output') + 1] = str(extended)
            extension_command[extension_command.index('--altitudes') + 1] = '0,1000'
            extension_command += ['--altitude-interpolation', 'log_geocentric_perigee_radius']
            reused = subprocess.run(extension_command + ['--reuse-study', str(output)],
                                    capture_output=True, text=True)
            self.assertEqual(reused.returncode, 0, reused.stderr)
            self.assertIn('2 cases; 1 already complete', reused.stdout)
            extension_study = json.loads((extended / 'study.json').read_text())
            self.assertTrue(extension_study['metadata']['complete'])
            self.assertEqual(extension_study['metadata']['altitude_interpolation'],
                             'log_geocentric_perigee_radius')
            self.assertEqual(extension_study['runs'][1], dict(study['runs'][0], run_id='case_000001'))
            self.assertEqual(len(list((extended / 'raw').glob('*.jsonl'))), 2)
            extension_resumed = subprocess.run(extension_command, capture_output=True, text=True)
            self.assertEqual(extension_resumed.returncode, 0, extension_resumed.stderr)
            self.assertIn('2 already complete', extension_resumed.stdout)
            self.assertEqual(raw[0].read_bytes(), before)
            config.write_text(config.read_text() + '\n# Configuration provenance changed\n')
            changed = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn('Resume rejected', changed.stderr)


if __name__ == '__main__':
    unittest.main()
