"""Independent analytic fixtures for grid integrity and interpolation audits."""
import copy
import importlib.util
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "validate_coverage_grid.py"
SPEC = importlib.util.spec_from_file_location("coverage_validator", SCRIPT)
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def fixture(orbit_axes):
    axes = [*orbit_axes, [1.0, 100.0], [0.01, 1.0]]
    metadata = dict(zip(VALIDATOR.AXES, axes))
    metadata.update(samples=20000, seed=20260919, complete=True,
                    expected_cases=math.prod(map(len, axes)))
    rows = []
    for index, point in enumerate(itertools.product(*axes)):
        altitude, inclination, eccentricity, position, velocity = point
        time = math.exp(0.0001 * altitude + 0.004 * inclination + 2 * eccentricity
                        - 0.05 * math.log(position) - 0.1 * math.log(velocity))
        cadence = time / 100
        row = dict(zip(VALIDATOR.FIELDS, point))
        row.update(run_id=f"node_{index}", samples=20000, seed=20260919, status="ok",
                   coverage_time_days=time, coverage_lower_days=time - cadence,
                   coverage_upper_days=time, coverage_confirmation_days=time + 2 * cadence,
                   onset_max_gap_deg=4.0, cadence_days=cadence,
                   evaluated_horizon_days=time + 2 * cadence, evaluated_epochs=103,
                   min_initial_perigee_km=900.0, min_mean_perigee_km=901.0,
                   earth_intersections=0, particles_below_300km=0,
                   particles_below_500km=0, particles_below_1000km=1)
        rows.append(row)
    metadata["completed_cases"] = len(rows)
    return {"metadata": metadata, "runs": rows}


def mask(row):
    row.update(status="earth_intersection", earth_intersections=1, particles_below_300km=1,
               particles_below_500km=1, particles_below_1000km=1,
               min_initial_perigee_km=-1, min_mean_perigee_km=None, evaluated_epochs=0,
               cadence_days=0, evaluated_horizon_days=0)
    row.update({key: None for key in VALIDATOR.EVENT_FIELDS})


def power_law_fixture(orbit_axes):
    """Independent position-dominated scaling, with separable angular effects."""
    document = fixture(orbit_axes)
    radius = 6378.137
    document["metadata"]["earth_radius_km"] = radius
    for row in document["runs"]:
        time = ((radius + row["perigee_altitude_km"]) / radius) ** 2.5 * math.exp(
            0.004 * row["inclination_deg"] + 2 * row["eccentricity"]
            - 0.05 * math.log(row["position_sigma_km"])
            - 0.1 * math.log(row["velocity_sigma_m_s"]))
        scale = time / row["coverage_time_days"]
        for key in ("coverage_time_days", "coverage_lower_days", "coverage_upper_days",
                    "coverage_confirmation_days", "cadence_days", "evaluated_horizon_days"):
            row[key] *= scale
    return document


class CoverageValidationTests(unittest.TestCase):
    def setUp(self):
        self.study = fixture([[1000.0, 3000.0], [0.0, 90.0], [0.0, 0.1]])
        self.holdouts = fixture([[2000.0], [45.0], [0.05]])

    def test_log_linear_physics_is_reproduced_exactly(self):
        result = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.assertTrue(result["passed"])
        self.assertFalse(result["provisional"])
        self.assertEqual(result["statistics"]["n"], 4)
        for name in ("max_relative_error", "p95_relative_error", "rms_relative_error"):
            self.assertLess(result["statistics"][name], 2e-15)
        for row in result["comparisons"]:
            self.assertEqual(len(row["support"]), 8)
            self.assertAlmostEqual(sum(node["weight"] for node in row["support"]), 1)

    def test_explicit_linear_altitude_preserves_legacy_predictions(self):
        legacy = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.study["metadata"]["altitude_interpolation"] = "linear_altitude"
        explicit = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.assertEqual(legacy, explicit)

    def test_power_law_is_exact_in_log_radius_at_geometric_midpoint(self):
        radius = 6378.137
        midpoint = math.sqrt((radius + 3000) * (radius + 30000)) - radius
        study = power_law_fixture([[3000.0, 30000.0], [0.0, 90.0], [0.0, 0.1]])
        holdouts = power_law_fixture([[midpoint], [45.0], [0.05]])
        study["metadata"]["altitude_interpolation"] = "log_geocentric_perigee_radius"
        # Direct measurements need no interpolation policy; legacy metadata is valid.
        result = VALIDATOR.validate_grid(study, holdouts)
        self.assertTrue(result["passed"])
        self.assertEqual(result["altitude_interpolation"], "log_geocentric_perigee_radius")
        self.assertLess(result["statistics"]["max_relative_error"], 3e-15)
        for row in result["comparisons"]:
            self.assertEqual(len(row["support"]), 8)
            for node in row["support"]:
                self.assertAlmostEqual(node["weight"], 0.125, places=14)
        holdouts["metadata"]["altitude_interpolation"] = "linear_altitude"
        self.assertEqual(result, VALIDATOR.validate_grid(study, holdouts))
        del study["metadata"]["altitude_interpolation"]
        self.assertGreater(VALIDATOR.validate_grid(study, holdouts)["statistics"][
            "max_relative_error"], 0.05)

    def test_log_radius_endpoints_reproduce_measured_values(self):
        study = power_law_fixture([[1000.0, 30000.0], [0.0, 90.0], [0.0, 0.1]])
        holdouts = power_law_fixture([[1000.0, 30000.0], [90.0], [0.1]])
        study["metadata"]["altitude_interpolation"] = "log_geocentric_perigee_radius"
        result = VALIDATOR.validate_grid(study, holdouts)
        self.assertTrue(result["passed"])
        self.assertLess(result["statistics"]["max_relative_error"], 2e-15)
        self.assertTrue(all(len(row["support"]) == 1 for row in result["comparisons"]))

    def test_bad_interpolation_policy_and_log_radius_are_rejected(self):
        for policy in ("log_altitude", "unknown", None, 1, []):
            with self.subTest(policy=policy):
                bad = copy.deepcopy(self.study)
                bad["metadata"]["altitude_interpolation"] = policy
                with self.assertRaisesRegex(ValueError, "altitude_interpolation"):
                    VALIDATOR.validate_grid(bad, self.holdouts)
        for radius in (None, 0, -1, float("nan"), float("inf"), "6378.137", True):
            with self.subTest(radius=radius):
                bad = copy.deepcopy(self.study)
                bad["metadata"].update(altitude_interpolation="log_geocentric_perigee_radius",
                                       earth_radius_km=radius)
                with self.assertRaisesRegex(ValueError, "earth radius"):
                    VALIDATOR.validate_grid(bad, self.holdouts)

    def test_floating_axis_roundoff_and_node_endpoints(self):
        self.holdouts = fixture([[1000.0], [0.0], [0.0]])
        for row in self.holdouts["runs"]:
            row["velocity_sigma_m_s"] *= 1 + 5e-13
        result = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.assertTrue(result["passed"])
        self.assertTrue(all(len(row["support"]) == 1 for row in result["comparisons"]))

    def test_invalid_corner_keeps_a_hole(self):
        mask(self.study["runs"][0])
        result = VALIDATOR.validate_grid(self.study, self.holdouts)
        blocked = result["comparisons"][0]
        self.assertEqual(blocked["status"], "invalid_corners")
        self.assertEqual(blocked["invalid_corners"][0]["status"], "earth_intersection")
        self.assertIsNone(blocked["predicted_time_days"])
        self.assertEqual(result["statistics"]["n"], 3)
        self.assertTrue(result["passed"])

    def test_missing_corner_is_distinct_and_requires_partial_mode(self):
        self.study["runs"].pop(0)
        self.study["metadata"].update(complete=False, completed_cases=31)
        with self.assertRaisesRegex(ValueError, "allow-partial"):
            VALIDATOR.validate_grid(self.study, self.holdouts)
        result = VALIDATOR.validate_grid(self.study, self.holdouts, allow_partial=True)
        self.assertTrue(result["provisional"])
        self.assertEqual(result["measured"]["missing_cases"], 1)
        self.assertEqual(result["comparisons"][0]["status"], "missing_corners")
        self.assertEqual(result["comparisons"][0]["invalid_corners"], [])
        self.assertIsNone(result["comparisons"][0]["predicted_time_days"])

    def test_masked_holdout_is_not_a_prediction_or_error(self):
        mask(self.holdouts["runs"][0])
        result = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.assertEqual(result["comparisons"][0]["status"], "withheld_masked")
        self.assertIsNone(result["comparisons"][0]["predicted_time_days"])
        self.assertEqual(result["statistics"]["n"], 3)

    def test_no_supported_comparisons_are_not_a_false_validation_pass(self):
        for row in self.holdouts["runs"]:
            mask(row)
        result = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.assertTrue(result["data_integrity_passed"])
        self.assertFalse(result["passed"])
        self.assertEqual(result["interpolation_gate"], "not_evaluated")
        self.assertEqual(result["statistics"]["n"], 0)
        self.assertIsNone(result["statistics"]["max_relative_error"])

    def test_orbit_extrapolation_and_uncertainty_interpolation_are_not_invented(self):
        outside = fixture([[3500.0], [45.0], [0.05]])
        result = VALIDATOR.validate_grid(self.study, outside)
        self.assertEqual(result["comparison_status_counts"], {"outside_orbit_grid": 4})
        changed = copy.deepcopy(self.holdouts)
        changed["metadata"]["position_sigmas_km"] = [2.0, 100.0]
        for row in changed["runs"]:
            if row["position_sigma_km"] == 1:
                row["position_sigma_km"] = 2
        result = VALIDATOR.validate_grid(self.study, changed)
        self.assertEqual(result["comparison_status_counts"]["missing_uncertainty_axis"], 2)
        self.assertEqual(result["statistics"]["n"], 2)

    def test_complete_grid_must_have_unique_full_cartesian_product(self):
        self.study["runs"].pop()
        self.study["metadata"]["completed_cases"] = 31
        with self.assertRaisesRegex(ValueError, "full Cartesian"):
            VALIDATOR.validate_grid(self.study, self.holdouts, allow_partial=True)
        self.setUp()
        duplicate = copy.deepcopy(self.study["runs"][0])
        duplicate["run_id"] = "different_id_same_point"
        self.study["runs"][-1] = duplicate
        with self.assertRaisesRegex(ValueError, "duplicate parameter"):
            VALIDATOR.validate_grid(self.study, self.holdouts)

    def test_bad_event_arithmetic_is_rejected(self):
        mutations = {"coverage_upper_days": 100, "coverage_lower_days": 0,
                     "coverage_confirmation_days": 100, "evaluated_horizon_days": 100,
                     "onset_max_gap_deg": 5.01, "cadence_days": -1,
                     "evaluated_epochs": 2, "coverage_time_days": float("nan")}
        for key, value in mutations.items():
            with self.subTest(field=key):
                bad = copy.deepcopy(self.study)
                bad["runs"][0][key] = value
                with self.assertRaises(ValueError):
                    VALIDATOR.validate_grid(bad, self.holdouts)

    def test_particle_counts_and_seeds_are_checked(self):
        for key, value in (("samples", 19999), ("seed", 7)):
            with self.subTest(field=key):
                bad = copy.deepcopy(self.study)
                bad["runs"][0][key] = value
                with self.assertRaisesRegex(ValueError, "count or seed"):
                    VALIDATOR.validate_grid(bad, self.holdouts)
        self.holdouts["metadata"]["seed"] = 7
        for row in self.holdouts["runs"]:
            row["seed"] = 7
        with self.assertRaisesRegex(ValueError, "different random seeds"):
            VALIDATOR.validate_grid(self.study, self.holdouts)

    def test_five_percent_gate_reports_actual_error(self):
        row = self.holdouts["runs"][0]
        for key in ("coverage_time_days", "coverage_lower_days", "coverage_upper_days",
                    "coverage_confirmation_days", "cadence_days", "evaluated_horizon_days"):
            row[key] *= 1.1
        result = VALIDATOR.validate_grid(self.study, self.holdouts)
        self.assertFalse(result["passed"])
        self.assertEqual(result["interpolation_gate"], "failed")
        self.assertEqual(result["statistics"]["exceeding_five_percent"], 1)
        self.assertAlmostEqual(result["statistics"]["max_relative_error"], 1 - 1 / 1.1)

    def test_cli_writes_auditable_json_for_success_and_bad_data(self):
        with tempfile.TemporaryDirectory(prefix="coverage-validation-") as directory:
            paths = [Path(directory) / name for name in ("study.json", "held.json", "report.json")]
            paths[0].write_text(json.dumps(self.study), encoding="utf-8")
            paths[1].write_text(json.dumps(self.holdouts), encoding="utf-8")
            command = [sys.executable, "-B", str(SCRIPT), "--study", str(paths[0]),
                       "--holdouts", str(paths[1]), "--output", str(paths[2])]
            result = subprocess.run(command, capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(json.loads(paths[2].read_text())["passed"])
            paths[0].write_text("{}", encoding="utf-8")
            result = subprocess.run(command, capture_output=True, text=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            report = json.loads(paths[2].read_text())
            self.assertEqual(report["status"], "invalid_data")
            self.assertFalse(report["data_integrity_passed"])
            self.assertTrue(report["errors"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
