"""Scientific provenance and lossless merging checks for report extensions."""
import copy
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from package_coverage_extension import audit_regression_holdouts, merge_studies, original_reference
from test_coverage_validation import fixture


def manifest(document):
    return {'executable_sha256': 'same-native-executable', 'config_sha256': 'same-config',
            'config_text': 'same force model and sampling settings',
            'metadata': copy.deepcopy(document['metadata'])}


def reference_fixture():
    """Canonical reference metadata for early validation, without particle files."""
    positions = [10 ** (index / 5) for index in range(11)]
    velocities = [10 ** (-2 + index / 5) for index in range(11)]
    metadata = {'samples': 20000, 'seed': 20260919, 'force_model': 'j2_j2sq',
                'initial_type': 'osculating', 'output_type': 'mean',
                'orbit_keplerian_deg': [26560000, 0.02, 55, 20, 30, 10],
                'position_sigmas_km': positions, 'velocity_sigmas_m_s': velocities}
    rows = [{'run_id': f'reference_{index}', 'position_sigma_km': position,
             'velocity_sigma_m_s': velocity} for index, (position, velocity) in
            enumerate(itertools.product(positions, velocities))]
    return {'metadata': metadata, 'runs': rows}


class CoveragePackageTests(unittest.TestCase):
    def setUp(self):
        self.base = fixture([[1000, 3000], [0, 90], [0, 0.1]])
        self.extra = fixture([[7000, 30000], [0, 90], [0, 0.1]])
        for document in (self.base, self.extra):
            document['metadata']['earth_radius_km'] = 6378.137
        self.documents = [self.base, self.extra]
        self.manifests = list(map(manifest, self.documents))

    def test_reused_measurements_are_unchanged_and_have_unique_ids(self):
        before = copy.deepcopy(self.documents)
        result = merge_studies(self.documents, self.manifests)
        self.assertEqual(self.documents, before)
        self.assertEqual(result['metadata']['perigee_altitudes_km'], [1000, 3000, 7000, 30000])
        self.assertEqual(result['metadata']['altitude_interpolation'], 'log_geocentric_perigee_radius')
        self.assertEqual(len({row['run_id'] for row in result['runs']}), 64)
        for row in result['runs']:
            original = next(item for item in self.documents[row['source_study']]['runs']
                            if item['run_id'] == row['source_run_id'])
            for field, value in original.items():
                if field != 'run_id':
                    self.assertEqual(row[field], value)

    def test_changed_executable_or_force_settings_cannot_be_mixed(self):
        for key in ('executable_sha256', 'config_sha256', 'config_text'):
            with self.subTest(key=key):
                changed = copy.deepcopy(self.manifests)
                changed[1][key] = 'different'
                with self.assertRaisesRegex(ValueError, 'provenance'):
                    merge_studies(self.documents, changed)

    def test_changed_covariance_axes_or_seed_cannot_be_mixed(self):
        for key, value in [('uncertainty_convention', 'correlated'), ('seed', 42)]:
            with self.subTest(key=key):
                docs = copy.deepcopy(self.documents)
                docs[1]['metadata'][key] = value
                if key == 'seed':
                    for row in docs[1]['runs']:
                        row['seed'] = value
                with self.assertRaisesRegex(ValueError, 'metadata differs'):
                    merge_studies(docs, self.manifests)

    def test_overlapping_sources_cannot_duplicate_measurements(self):
        with self.assertRaisesRegex(ValueError, 'overlap'):
            merge_studies([self.base, self.base], [self.manifests[0], self.manifests[0]])

    def test_incomplete_source_cannot_be_packaged_as_complete(self):
        self.extra['runs'].pop()
        self.extra['metadata'].update(complete=False, completed_cases=31)
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            merge_studies(self.documents, self.manifests)

    def test_manifest_cannot_misrepresent_source_axes(self):
        self.manifests[1]['metadata']['perigee_altitudes_km'] = [7000, 29000]
        with self.assertRaisesRegex(ValueError, 'Manifest axes'):
            merge_studies(self.documents, self.manifests)

    def test_manifest_cannot_misrepresent_scientific_metadata(self):
        self.manifests[1]['metadata']['seed'] = 42
        with self.assertRaisesRegex(ValueError, 'Manifest metadata.*seed'):
            merge_studies(self.documents, self.manifests)

    def test_original_reference_rejects_a_different_experiment(self):
        mutations = {'seed': 42, 'samples': 1000, 'initial_type': 'mean',
                     'output_type': 'osculating', 'force_model': 'j2',
                     'orbit_keplerian_deg': [26560000, 0.03, 55, 20, 30, 10]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'reference.json'
            for key, value in mutations.items():
                with self.subTest(key=key):
                    document = reference_fixture()
                    document['metadata'][key] = value
                    path.write_text(json.dumps(document), encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'reference metadata differs'):
                        original_reference(path)

    def test_original_reference_requires_complete_unique_uncertainty_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'reference.json'
            for mutation in ('missing', 'duplicate_point', 'duplicate_id', 'bad_axis'):
                with self.subTest(mutation=mutation):
                    document = reference_fixture()
                    if mutation == 'missing':
                        document['runs'].pop()
                    elif mutation == 'duplicate_point':
                        document['runs'][-1].update(position_sigma_km=1, velocity_sigma_m_s=0.01)
                    elif mutation == 'duplicate_id':
                        document['runs'][-1]['run_id'] = document['runs'][0]['run_id']
                    else:
                        document['metadata']['position_sigmas_km'][1] = 0
                    path.write_text(json.dumps(document), encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'requires'):
                        original_reference(path)

    def test_legacy_holdout_audit_preserves_distinct_executable_provenance(self):
        grid = merge_studies(self.documents, self.manifests)
        held = fixture([[2000], [45], [0.05]])
        held['metadata']['earth_radius_km'] = 6378.137
        held_manifest = manifest(held)
        held_manifest['executable_sha256'] = 'previous-validated-implementation'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'study.json').write_text(json.dumps(held), encoding='utf-8')
            (root / 'manifest.json').write_text(json.dumps(held_manifest), encoding='utf-8')
            report = audit_regression_holdouts(grid, root)
            self.assertTrue(report['passed'])
            self.assertEqual(report['altitude_interpolation'], 'log_geocentric_perigee_radius')
            self.assertEqual(report['source']['executable_sha256'], 'previous-validated-implementation')
            self.assertEqual(len(report['source']['study_sha256']), 64)
            held_manifest['metadata']['seed'] = 42
            (root / 'manifest.json').write_text(json.dumps(held_manifest), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Regression manifest metadata'):
                audit_regression_holdouts(grid, root)


if __name__ == '__main__':
    unittest.main()
