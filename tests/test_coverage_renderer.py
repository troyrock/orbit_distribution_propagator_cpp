"""Portable coverage-report contracts; no plotting package or browser required.

Run: python tests/test_coverage_renderer.py -v
"""

from __future__ import annotations

import copy
import importlib.util
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("coverage_renderer", REPO / "tools/render_coverage_explorer.py")
renderer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(renderer)


def fixture():
    """Analytic test values only; intentionally not claimed as DSST output."""
    metadata = {
        "samples": 20000,
        "seed": 20260919,
        "perigee_altitudes_km": [1000, 3000],
        "inclinations_deg": [0, 90],
        "eccentricities": [0, 0.1],
        "position_sigmas_km": [1, 10, 100],
        "velocity_sigmas_m_s": [0.01, 0.1, 1],
        "coverage_definition": "Synthetic test criterion; not simulation results.",
        "model_scope": "Synthetic renderer test.",
        "synthetic_test_data": True,
    }
    runs = []
    axes = [metadata[name] for name, _ in renderer.AXES]
    for a, i, e, p, v in itertools.product(*axes):
        time = math.exp(a / 5000 + i / 180 + e * 2) * 70 / math.hypot(p, v * 3)
        runs.append({
            "run_id": str(len(runs)), "samples": 20000,
            "perigee_altitude_km": a, "inclination_deg": i, "eccentricity": e,
            "position_sigma_km": p, "velocity_sigma_m_s": v,
            "status": "ok", "coverage_time_days": time,
            "coverage_lower_days": 0.99 * time, "coverage_upper_days": time,
            "min_initial_perigee_km": a - p * 8,
            "particles_below_300km": 5 if a - p * 8 < 300 else 0,
            "particles_below_500km": 15 if a - p * 8 < 500 else 0,
        })
    return {"metadata": metadata, "runs": runs}


def reference_fixture():
    """Extended altitude grid with a separately preserved measured surface."""
    data = fixture()
    data["metadata"]["perigee_altitudes_km"][-1] = 30000
    data["metadata"]["altitude_interpolation"] = "log_geocentric_perigee_radius"
    data["metadata"]["earth_radius_km"] = 6378.137
    for run in data["runs"]:
        if run["perigee_altitude_km"] == 3000:
            run["perigee_altitude_km"] = 30000
    reference = {
        "label": "Original MEO / PDF", "description": "Preserved original measured times.",
        "metadata": copy.deepcopy(data["metadata"]), "runs": [],
    }
    reference["metadata"].update(perigee_altitudes_km=[19650.663], inclinations_deg=[55], eccentricities=[0.02])
    for run in data["runs"][:9]:
        reference_run = copy.deepcopy(run)
        reference_run.update(perigee_altitude_km=19650.663, inclination_deg=55, eccentricity=0.02)
        reference["runs"].append(reference_run)
    data["reference_surface"] = reference
    return data


class CoverageRendererTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="coverage-renderer-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.library = self.root / "plotly.min.js"
        # The unit tests verify report packaging; real rendering is browser-tested.
        self.library.write_text("window.Plotly = {}; // Local test stub\n", encoding="utf-8")
        self.data = fixture()

    def test_successful_fixture_and_all_unavailable_statuses(self):
        renderer.validate_data(self.data)
        for status in ("earth_intersection", "mean_earth_intersection", "no_phase_shear", "unobserved"):
            with self.subTest(status=status):
                candidate = copy.deepcopy(self.data)
                candidate["runs"][0].update(status=status, coverage_time_days=None,
                                            coverage_lower_days=None, coverage_upper_days=None)
                renderer.validate_data(candidate)
                self.assertIn(status, renderer.render(candidate, self.library))
        # Partial grids are valid: absent coordinates must remain browser holes.
        del self.data["runs"][0]
        renderer.validate_data(self.data)
        self.data["runs"] = []
        renderer.validate_data(self.data)

    def test_duplicate_ids_and_coordinates_are_rejected(self):
        for duplicate_id in (False, True):
            with self.subTest(duplicate_id=duplicate_id):
                candidate = copy.deepcopy(self.data)
                row = copy.deepcopy(candidate["runs"][0])
                if not duplicate_id:
                    row["run_id"] = "another-id"
                candidate["runs"].append(row)
                with self.assertRaisesRegex(ValueError, "unique|duplicate"):
                    renderer.validate_data(candidate)

    def test_nonfinite_values_and_nonpositive_log_axes_are_rejected(self):
        for value in (float("nan"), float("inf"), -float("inf"), True):
            for field in ("coverage_time_days", "coverage_lower_days", "perigee_altitude_km",
                          "min_initial_perigee_km"):
                with self.subTest(value=value, field=field):
                    candidate = copy.deepcopy(self.data)
                    candidate["runs"][0][field] = value
                    with self.assertRaises(ValueError):
                        renderer.validate_data(candidate)
        for name in ("position_sigmas_km", "velocity_sigmas_m_s"):
            candidate = copy.deepcopy(self.data)
            candidate["metadata"][name][0] = 0
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "positive"):
                renderer.validate_data(candidate)
        self.data["metadata"]["eccentricities"][0] = float("nan")
        with self.assertRaisesRegex(ValueError, "finite"):
            renderer.validate_data(self.data)

    def test_grid_coordinates_status_and_sample_count_contracts(self):
        cases = (("samples", 19999), ("status", "unknown"), ("perigee_altitude_km", 1777),
                 ("particles_below_300km", 20001), ("particles_below_500km", 2.5))
        for field, value in cases:
            with self.subTest(field=field):
                candidate = copy.deepcopy(self.data)
                candidate["runs"][0][field] = value
                with self.assertRaises(ValueError):
                    renderer.validate_data(candidate)
        self.data["metadata"]["inclinations_deg"] = [90, 0]
        with self.assertRaisesRegex(ValueError, "increasing"):
            renderer.validate_data(self.data)

    def test_brackets_require_order_and_contain_the_reported_time(self):
        time = self.data["runs"][0]["coverage_time_days"]
        for lower, upper in ((time + 1, time), (time + 1, time + 2), (0, time - 1), (-1, time)):
            with self.subTest(lower=lower, upper=upper):
                candidate = copy.deepcopy(self.data)
                candidate["runs"][0].update(coverage_lower_days=lower, coverage_upper_days=upper)
                with self.assertRaises(ValueError):
                    renderer.validate_data(candidate)
        self.data["runs"][0].update(coverage_lower_days=0, coverage_upper_days=None)
        renderer.validate_data(self.data)

    def test_embedded_json_round_trips_without_executable_html(self):
        malicious = "</script><script>window.injected=true</script>&<>\u2028\u2029"
        self.data["metadata"]["model_scope"] = malicious
        report = renderer.render(self.data, self.library)
        self.assertNotIn("<script>window.injected", report)
        self.assertNotIn("__COVERAGE_EXPLORER_DATA__", report)
        self.assertNotIn("__COVERAGE_EXPLORER_PLOTLY__", report)
        self.assertIn("\\u003c/script\\u003e", report)
        payload = re.search(r'<script id="coverage-data" type="application/json">(.*?)</script>',
                            report, re.DOTALL).group(1)
        self.assertEqual(json.loads(payload), self.data)
        self.assertIn(self.library.read_text(encoding="utf-8"), report)

    def test_altitude_interpolation_policy_is_explicit_and_validated(self):
        self.data["metadata"]["earth_radius_km"] = 6378.137
        for policy in ("linear_altitude", "log_geocentric_perigee_radius"):
            with self.subTest(policy=policy):
                self.data["metadata"]["altitude_interpolation"] = policy
                renderer.validate_data(self.data)
        self.data["metadata"]["altitude_interpolation"] = "unsupported"
        with self.assertRaisesRegex(ValueError, "altitude_interpolation"):
            renderer.validate_data(self.data)
        self.data["metadata"].update(altitude_interpolation="log_geocentric_perigee_radius",
                                     perigee_altitudes_km=[-7000, 3000])
        with self.assertRaisesRegex(ValueError, "positive geocentric"):
            renderer.validate_data(self.data)
        del self.data["metadata"]["earth_radius_km"]
        with self.assertRaisesRegex(ValueError, "explicit metadata.earth_radius_km"):
            renderer.validate_data(self.data)

    def test_reference_surface_preserves_the_measured_dataset(self):
        data = reference_fixture()
        data["metadata"]["initial_view"] = "reference"
        renderer.validate_data(data)
        report = renderer.render(data, self.library)
        payload = re.search(r'<script id="coverage-data" type="application/json">(.*?)</script>',
                            report, re.DOTALL).group(1)
        self.assertEqual(json.loads(payload)["reference_surface"], data["reference_surface"])

    def test_initial_view_requires_an_available_reference(self):
        self.data["metadata"]["initial_view"] = "reference"
        with self.assertRaisesRegex(ValueError, "requires a reference_surface"):
            renderer.validate_data(self.data)
        self.data["metadata"]["initial_view"] = "unknown"
        with self.assertRaisesRegex(ValueError, "initial_view"):
            renderer.validate_data(self.data)

    def test_reference_surface_rejects_incompatible_or_ambiguous_metadata(self):
        mutations = (
            lambda r: r.update(label=""),
            lambda r: r.update(description=123),
            lambda r: r.update(reference_surface={}),
            lambda r: r["metadata"].update(samples=19999),
            lambda r: r["metadata"].update(earth_radius_km=6000),
            lambda r: r["metadata"].update(perigee_altitudes_km=[40000]),
            lambda r: r["metadata"].update(inclinations_deg=[55, 60]),
            lambda r: r["metadata"].update(velocity_sigmas_m_s=[0.01, 0.1, 1, 2]),
        )
        for mutation in mutations:
            data = reference_fixture()
            mutation(data["reference_surface"])
            with self.subTest(reference=data["reference_surface"]["metadata"]), self.assertRaises(ValueError):
                renderer.validate_data(data)

    def test_template_and_bundle_errors_are_actionable(self):
        broken = self.root / "template.html"
        broken.write_text("<html>No insertion tokens</html>", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            renderer.render(self.data, self.library, broken)
        with self.assertRaisesRegex(FileNotFoundError, "Plotly JavaScript not found"):
            renderer.find_plotly(self.root / "missing.js")

    def test_cli_writes_nested_output_and_rejects_invalid_data(self):
        source, destination = self.root / "fixture.json", self.root / "nested/report.html"
        source.write_text(json.dumps(self.data), encoding="utf-8")
        command = [sys.executable, str(REPO / "tools/render_coverage_explorer.py"), str(source),
                   "--output", str(destination), "--plotly-js", str(self.library)]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(destination.is_file())
        self.assertIn("recorded grid nodes", result.stdout)
        self.data["runs"][0]["coverage_time_days"] = None
        source.write_text(json.dumps(self.data), encoding="utf-8")
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("successful coverage requires", result.stderr)


if __name__ == "__main__":
    unittest.main()
