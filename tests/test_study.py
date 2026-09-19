"""Stdlib regression tests for study orchestration and result interpretation.

Run: python tests/test_study.py EXE STUDY_SCRIPT [-v or unittest test name]
The suite runs six tiny successful child simulations and one invalid child.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest


if len(sys.argv) < 3:
    raise SystemExit(__doc__)
EXECUTABLE = Path(sys.argv.pop(1)).resolve()
STUDY_SCRIPT = Path(sys.argv.pop(1)).resolve()
if not EXECUTABLE.is_file() or not STUDY_SCRIPT.is_file():
    raise SystemExit("Executable and study script must both exist")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


class StudyRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="distribution-study-tests-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name).resolve()
        cls.config = cls.root / "tiny.cfg"
        cls.config.write_text(
            "force_model=kepler\ninitial_type=mean\noutput_type=mean\n"
            "orbit_equinoctial=26560000,0,0,0,0,0\n"
            "sigma=1000,0,0,0,0,0.001\n"
            "duration_days=0.5\noutput_step_days=0.25\n"
            "phase_bins=8\npersistence=2\n", encoding="utf-8")
        cls.output = cls.root / "complete"
        cls.base = [sys.executable, "-B", str(STUDY_SCRIPT), "--exe", str(EXECUTABLE),
                    "--config", str(cls.config)]
        command = cls.base + ["--output", str(cls.output), "--counts", "4,8", "--seeds", "11",
                              "--threads", "2", "--cadences", "0.25,0.5", "--benchmark-threads", "1,2"]
        completed = subprocess.run(command, capture_output=True, text=True, timeout=90, check=False)
        if completed.returncode:
            logs = "\n".join(path.read_text(encoding="utf-8", errors="replace")[-4000:]
                             for path in cls.output.glob("runs/*/process.log"))
            raise RuntimeError(f"Tiny study failed ({completed.returncode}):\n"
                               f"{completed.stdout}\n{completed.stderr}\n{logs}")
        cls.report = read_json(cls.output / "study.json")
        cls.module = runpy.run_path(str(STUDY_SCRIPT), run_name="study_regression_import")

    def test_completed_matrix_and_output_dimensions(self):
        report = self.report
        self.assertEqual(report["status"], "complete")
        self.assertIsNone(report["active_run"])
        self.assertEqual(report["planned_run_count"], 6)
        self.assertEqual(len(report["runs"]), 6)
        science = [row for row in report["runs"] if row["kind"] == "science"]
        self.assertEqual({(row["samples"], row["seed"], row["cadence_days"]) for row in science},
                         {(4, 11, .25), (4, 11, .5), (8, 11, .25), (8, 11, .5)})
        self.assertEqual(report["config_sha256"], hashlib.sha256(self.config.read_bytes()).hexdigest())
        self.assertEqual(report["executable_sha256"], hashlib.sha256(EXECUTABLE.read_bytes()).hexdigest())
        for row in report["runs"]:
            with self.subTest(run=row["run_id"]):
                child = Path(row["output_directory"])
                self.assertTrue(Path(row["log_path"]).is_file())
                self.assertFalse((child / "visualization.html").exists())
                run = read_json(child / "run.json")
                expected_frames = 3 if row["cadence_days"] == .25 else 2
                self.assertEqual(row["frame_count"], expected_frames)
                self.assertEqual(len(run["frames"]), expected_frames)
                self.assertEqual(run["metadata"]["samples"], row["samples"])
                self.assertEqual(run["metadata"]["visual_samples"], 0)
                for frame in run["frames"]:
                    self.assertEqual(frame["positions_m"], [])
                    self.assertEqual(sum(frame["metrics"]["histogram"]), row["samples"])
                for filename, field in (("metrics.csv", "metrics_sha256"),
                                        ("initial_samples.csv", "initial_samples_sha256")):
                    self.assertEqual(row[field], hashlib.sha256((child / filename).read_bytes()).hexdigest())

    def test_exact_thread_benchmark_comparison(self):
        benchmarks = [row for row in self.report["runs"] if row["kind"] == "benchmark"]
        self.assertEqual([row["requested_threads"] for row in benchmarks], [1, 2])
        first = benchmarks[0]
        for row in benchmarks:
            self.assertTrue(row["metrics_match_baseline"])
            self.assertTrue(row["initial_samples_match_baseline"])
            self.assertEqual(row["metrics_sha256"], first["metrics_sha256"])
            self.assertEqual(row["initial_samples_sha256"], first["initial_samples_sha256"])
            self.assertGreater(row["wall_seconds"], 0)
            self.assertGreater(row["simulation_elapsed_seconds"], 0)
            self.assertAlmostEqual(row["simulation_speedup_vs_first"],
                                   first["simulation_elapsed_seconds"] / row["simulation_elapsed_seconds"])

    def test_unobserved_events_and_unchecked_accuracy_stay_null(self):
        for row in self.report["runs"]:
            self.assertEqual(row["accuracy_checked_samples"], 0)
            self.assertIsNone(row["accuracy_max_position_m"])
            self.assertIsNone(row["accuracy_max_phase_rad"])
            for event in self.module["EVENTS"]:
                self.assertIsNone(row[f"{event}_time_s"])
                self.assertIsNone(row[f"{event}_lower_s"])
                self.assertIsNone(row[f"{event}_upper_s"])
        self.assertEqual(len(self.report["seed_scatter"]), 4)
        for group in self.report["seed_scatter"]:
            self.assertEqual(group["completed_seeds"], [11])
            for event in group["events"].values():
                self.assertEqual(event["observed_count"], 0)
                self.assertEqual(event["not_observed_count"], 1)
                self.assertIsNone(event["median_s"])
                self.assertIsNone(event["sample_stddev_s"])
        with (self.output / "study.csv").open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 6)
        for row in rows:
            self.assertEqual(row["mixing_time_s"], "")
            self.assertEqual(row["accuracy_max_position_m"], "")
            self.assertEqual(row["accuracy_checked_samples"], "0")

    def test_nonempty_output_is_rejected_without_mutation(self):
        before_json = (self.output / "study.json").read_bytes()
        before_csv = (self.output / "study.csv").read_bytes()
        completed = subprocess.run(self.base + ["--output", str(self.output), "--counts", "4", "--seeds", "11"],
                                   capture_output=True, text=True, timeout=30, check=False)
        self.assertEqual(completed.returncode, 2)
        self.assertIn("never reused", completed.stderr)
        self.assertEqual((self.output / "study.json").read_bytes(), before_json)
        self.assertEqual((self.output / "study.csv").read_bytes(), before_csv)

    def test_failed_child_retains_checkpoint_and_log(self):
        bad = self.root / "invalid.cfg"
        bad.write_text(self.config.read_text().replace("force_model=kepler", "force_model=invalid"), encoding="utf-8")
        output = self.root / "failed"
        command = [sys.executable, "-B", str(STUDY_SCRIPT), "--exe", str(EXECUTABLE),
                   "--config", str(bad), "--output", str(output), "--counts", "4", "--seeds", "11"]
        completed = subprocess.run(command, capture_output=True, text=True, timeout=30, check=False)
        self.assertEqual(completed.returncode, 1)
        report = read_json(output / "study.json")
        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["runs"], [])
        self.assertEqual(report["planned_run_count"], 1)
        self.assertIn("returned", report["error"])
        self.assertIsNotNone(report["active_run"])
        log = Path(report["active_run"]["log_path"])
        self.assertTrue(log.is_file())
        self.assertIn("force_model", log.read_text(encoding="utf-8"))
        with (output / "study.csv").open(newline="", encoding="utf-8") as source:
            self.assertEqual(list(csv.DictReader(source)), [])

    def test_zero_epoch_events_remain_numeric_zero(self):
        summary = {"sample_count": 4, "threads_used": 1, "frame_count": 2,
                   "elapsed_seconds": 1.0, "accuracy_checked_samples": 0,
                   "accuracy_max_position_m": 0, "accuracy_max_phase_rad": 0,
                   "dkw_cdf_error_95": .3, "analytic_mixing_time_s": 0}
        for event in self.module["EVENTS"]:
            summary[f"{event}_time_s"] = 0
            summary[f"{event}_bracket_s"] = [0, 0]
        row = self.module["extract_summary"](summary, 4)
        for event in self.module["EVENTS"]:
            self.assertIsNotNone(row[f"{event}_time_s"])
            self.assertEqual(row[f"{event}_time_s"], 0)
            self.assertEqual(row[f"{event}_lower_s"], 0)
            self.assertEqual(row[f"{event}_upper_s"], 0)
        self.assertEqual(row["analytic_mixing_time_s"], 0)
        self.assertIsNone(row["accuracy_max_position_m"])

    def test_seed_scatter_is_descriptive_and_preserves_missing_counts(self):
        common = {"kind": "science", "samples": 8, "cadence_days": .25}
        runs = [
            {**common, "seed": 11, "coverage_time_s": 0, "mixing_time_s": 0, "central95_wrap_time_s": None},
            {**common, "seed": 22, "coverage_time_s": None, "mixing_time_s": 10, "central95_wrap_time_s": None},
            {**common, "kind": "benchmark", "seed": 33, "coverage_time_s": 100,
             "mixing_time_s": 100, "central95_wrap_time_s": 100},
        ]
        groups = self.module["summarize_runs"](runs)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["completed_seeds"], [11, 22])
        events = groups[0]["events"]
        self.assertEqual(events["coverage"]["observed_count"], 1)
        self.assertEqual(events["coverage"]["not_observed_count"], 1)
        self.assertEqual(events["coverage"]["mean_s"], 0)
        self.assertIsNone(events["coverage"]["sample_stddev_s"])
        self.assertEqual(events["mixing"]["median_s"], 5)
        self.assertEqual(events["mixing"]["mean_s"], 5)
        self.assertAlmostEqual(events["mixing"]["sample_stddev_s"], math.sqrt(50))
        self.assertEqual(events["central95_wrap"]["not_observed_count"], 2)
        self.assertIsNone(events["central95_wrap"]["median_s"])


if __name__ == "__main__":
    unittest.main()
