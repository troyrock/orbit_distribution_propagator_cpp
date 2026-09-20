"""Functional contract of the native coverage CSV/JSONL command.

Run with the executable path as the sole positional argument. All fixtures live
in an isolated temporary directory; no existing study outputs are changed.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


HEADER = ("run_id,perigee_altitude_km,inclination_deg,eccentricity,"
          "position_sigma_km,velocity_sigma_m_s")
EXECUTABLE = Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else None


class CoverageCliTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="distribution-coverage-cli-")
        self.root = Path(self.directory.name)
        self.config = self.root / "base.cfg"
        self.cases = self.root / "cases.csv"
        self.output = self.root / "result.jsonl"
        self.write_config()
        self.write_cases(["first,2000,45,0.05,10,0.1"])

    def tearDown(self):
        self.directory.cleanup()

    def write_config(self, **overrides):
        values = {
            "force_model": "j2_j2sq",
            "initial_type": "osculating",
            "output_type": "mean",
            "orbit_keplerian_deg": "26560000, 0.02, 55, 20, 30, 10",
            "uncertainty_coordinates": "rtn",
            "sigma": "1000, 1000, 1000, 0.01, 0.01, 0.01",
            "samples": 256,
            "seed": 20260919,
            "threads": 1,
            "duration_days": 1,
            "output_step_days": 0.1,
            "visual_samples": 0,
            "phase_bins": 8,
            "coverage_max_gap_deg": 45,
            "coverage_occupied_fraction": 1,
            "persistence": 3,
            "minimum_altitude_m": 1000000,
            "max_memory_mb": 128,
        }
        values.update(overrides)
        if 'empirical_samples_csv' in overrides:
            values.pop('sigma', None)
        self.config.write_text("".join(f"{key} = {value}\n" for key, value in values.items()),
                               encoding="ascii")

    def write_cases(self, rows, header=HEADER):
        self.cases.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="ascii")

    def invoke(self, *extra, full=True):
        arguments = [str(EXECUTABLE)]
        if full:
            arguments += ["--config", str(self.config), "--cases", str(self.cases),
                          "--output", str(self.output)]
        return subprocess.run([*arguments, *extra], capture_output=True, text=True, timeout=20)

    def rejected(self, result, message=None):
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("Error:", result.stderr)
        if message:
            self.assertIn(message, result.stderr)

    def test_help(self):
        result = self.invoke("--help", full=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--cases CASES.csv", result.stdout)
        self.assertIn(HEADER, result.stdout)
        self.assertFalse(self.output.exists())

    def test_argument_validation(self):
        for value in ("nan", "inf", "1.5", "0", "1025", "2x"):
            with self.subTest(value=value):
                self.rejected(self.invoke("--threads", value, full=False))
        for arguments in (("--threads",), ("--unknown", "value"), ()):
            with self.subTest(arguments=arguments):
                self.rejected(self.invoke(*arguments, full=False))
        self.assertFalse(self.output.exists())

    def test_never_overwrites_output(self):
        sentinel = b"existing completed study\n"
        self.output.write_bytes(sentinel)
        self.rejected(self.invoke(), "Output already exists")
        self.assertEqual(self.output.read_bytes(), sentinel)

    def test_exact_csv_header_and_column_count(self):
        for header in (HEADER + ",", HEADER + ",\r", HEADER.replace("run_id", "id"), ""):
            with self.subTest(header=header):
                self.write_cases(["first,2000,45,0.05,10,0.1"], header)
                self.rejected(self.invoke(), "Wrong cases CSV header")
                self.assertFalse(self.output.exists())
        for row in ("first,2000,45,0.05,10", "first,2000,45,0.05,10,0.1,",
                    "first,2000,45,0.05,10,0.1,\r"):
            with self.subTest(row=row):
                self.write_cases([row])
                self.rejected(self.invoke(), "Wrong case column count")
                self.assertFalse(self.output.exists())

    def test_invalid_cases_reject_before_output_creation(self):
        rows = [
            ["same,2000,45,0.05,10,0.1", "same,2000,45,0.05,10,0.1"],
            ['bad"id,2000,45,0.05,10,0.1'],
            [",2000,45,0.05,10,0.1"],
            ["bad,-1,45,0.05,10,0.1"],
            ["bad,2000,91,0.05,10,0.1"],
            ["bad,2000,45,1,10,0.1"],
            ["bad,2000,45,0.05,0,0.1"],
            ["bad,2000,45,0.05,10,-1"],
            ["bad,2000,45,nan,10,0.1"],
            ["bad,2000,45,0.05,inf,0.1"],
            ["bad,2000,45,0.05,10,0.1junk"],
            [],
        ]
        for case_rows in rows:
            with self.subTest(rows=case_rows):
                self.write_cases(case_rows)
                self.rejected(self.invoke())
                self.assertFalse(self.output.exists())

    def test_small_batch_records_every_case(self):
        self.write_cases(["first,2000,0,0,10,0.1", "second,3000,90,0.1,100,1"])
        result = self.invoke("--threads", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [json.loads(line) for line in self.output.read_text().splitlines()]
        self.assertEqual([row["run_id"] for row in rows], ["first", "second"])
        for row in rows:
            self.assertEqual(row["status"], "ok")
            self.assertEqual(row["samples"], 256)
            self.assertEqual(row["threads_used"], 2)
            self.assertEqual(row["earth_intersections"], 0)
            self.assertEqual(row["seed"], 20260919)
            self.assertGreater(row["min_initial_perigee_km"], 0)
            self.assertLessEqual(row["coverage_lower_days"], row["coverage_time_days"])
            self.assertEqual(row["coverage_time_days"], row["coverage_upper_days"])
            self.assertAlmostEqual(row["coverage_confirmation_days"] - row["coverage_time_days"],
                                   2 * row["cadence_days"], places=10)
            self.assertLessEqual(row["onset_max_gap_deg"], 45)

    def test_fixed_times_include_final_partial_interval(self):
        self.write_config(samples=8, duration_days=0.1, output_step_days=0.03)
        self.write_cases(["short,2000,45,0.05,1,0.01"])
        result = self.invoke("--fixed-times")
        self.assertEqual(result.returncode, 0, result.stderr)
        row = json.loads(self.output.read_text())
        self.assertEqual(row["status"], "unobserved")
        self.assertEqual(row["evaluated_epochs"], 5)
        self.assertEqual(row["samples"], 8)
        self.assertIsNone(row["coverage_time_days"])
        self.assertIsNone(row["coverage_confirmation_days"])
        self.assertAlmostEqual(row["cadence_days"], 0.03)
        self.assertAlmostEqual(row["evaluated_horizon_days"], 0.1)

    def test_earth_crossing_case_is_masked_without_redraw(self):
        self.write_config(samples=20000, phase_bins=72, coverage_max_gap_deg=5)
        self.write_cases(["tail,1000,55,0,100,0.01"])
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = json.loads(self.output.read_text())
        self.assertEqual(row["status"], "earth_intersection")
        self.assertEqual(row["samples"], 20000)
        self.assertEqual(row["earth_intersections"], 4)
        self.assertEqual(row["particles_below_300km"], 174)
        self.assertEqual(row["particles_below_500km"], 1037)
        self.assertEqual(row["particles_below_1000km"], 13247)
        self.assertAlmostEqual(row["min_initial_perigee_km"], -63.097, places=2)
        self.assertIsNone(row["coverage_time_days"])
        self.assertIsNone(row["min_mean_perigee_km"])
        self.assertEqual(row["evaluated_epochs"], 0)

    def test_empirical_input_is_not_silently_replaced(self):
        self.write_config(empirical_samples_csv="samples.csv")
        self.rejected(self.invoke(), "requires Gaussian covariance input")
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    if EXECUTABLE is None or not EXECUTABLE.is_file():
        raise SystemExit("Usage: test_coverage_cli.py /path/to/distribution_coverage")
    unittest.main(verbosity=2)
