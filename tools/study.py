#!/usr/bin/env python3
"""Run reproducible Gaussian ensemble, cadence, and thread convergence studies.

Uses only Python's standard library. Each child invocation gets a fresh output
directory. Completed outputs are retained if a later run fails; rerunning a
study requires a new or empty study directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone


EVENTS = ("coverage", "mixing", "central95_wrap")
CSV_FIELDS = [
    "run_id", "kind", "samples", "seed", "cadence_days", "requested_threads",
    "threads_used", "frame_count", "wall_seconds", "simulation_elapsed_seconds",
    "coverage_time_s", "coverage_lower_s", "coverage_upper_s",
    "mixing_time_s", "mixing_lower_s", "mixing_upper_s",
    "central95_wrap_time_s", "central95_wrap_lower_s", "central95_wrap_upper_s",
    "analytic_mixing_time_s", "dkw_cdf_error_95", "accuracy_checked_samples",
    "accuracy_max_position_m", "accuracy_max_phase_rad",
    "metrics_sha256", "initial_samples_sha256", "metrics_match_baseline",
    "initial_samples_match_baseline", "simulation_speedup_vs_first", "output_directory",
]


class StudyError(RuntimeError):
    """A completed or attempted invocation failed a study contract."""


def comma_integers(text: str, minimum: int, maximum: int) -> list[int]:
    fields = text.split(",")
    try:
        values = [int(field.strip()) for field in fields]
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected comma-separated integers") from error
    if any(value < minimum or value > maximum for value in values):
        raise argparse.ArgumentTypeError(f"Values must be between {minimum} and {maximum}")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("Duplicate values would repeat a study case")
    return values


def comma_cadences(text: str) -> list[float]:
    try:
        values = [float(field.strip()) for field in text.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected comma-separated cadences in days") from error
    if any(not math.isfinite(value) or value <= 0 for value in values):
        raise argparse.ArgumentTypeError("Cadences must be positive and finite")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("Duplicate cadence values")
    return values


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_json(path: Path) -> dict:
    def reject_constant(value: str) -> None:
        raise StudyError(f"Nonfinite JSON constant in {path}: {value}")

    with path.open(encoding="utf-8-sig") as source:
        value = json.load(source, parse_constant=reject_constant)
    if not isinstance(value, dict):
        raise StudyError(f"Expected a JSON object in {path}")
    return value


def numeric(value: object, name: str, nullable: bool = False) -> float | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise StudyError(f"Invalid numeric summary field {name}: {value!r}")
    if not math.isfinite(value) or value < 0:
        raise StudyError(f"Invalid nonnegative summary field {name}: {value!r}")
    return value


def summarize_runs(runs: list[dict]) -> list[dict]:
    groups: dict[tuple[int, float | None], list[dict]] = {}
    for row in runs:
        if row["kind"] == "science":
            groups.setdefault((row["samples"], row["cadence_days"]), []).append(row)
    summaries = []
    for (count, cadence), rows in groups.items():
        group = {"samples": count, "cadence_days": cadence,
                 "completed_seeds": [row["seed"] for row in rows], "events": {}}
        for event in EVENTS:
            values = [row[f"{event}_time_s"] for row in rows if row[f"{event}_time_s"] is not None]
            group["events"][event] = {
                "observed_count": len(values),
                "not_observed_count": len(rows) - len(values),
                "minimum_s": min(values) if values else None,
                "median_s": statistics.median(values) if values else None,
                "maximum_s": max(values) if values else None,
                "mean_s": statistics.mean(values) if values else None,
                "sample_stddev_s": statistics.stdev(values) if len(values) > 1 else None,
                "interpretation": "Descriptive scatter among observed seed runs only; not a confidence interval. Unobserved events are not zero.",
            }
        summaries.append(group)
    return summaries


def checkpoint(directory: Path, report: dict) -> None:
    report["updated_at_utc"] = utc_now()
    report["seed_scatter"] = summarize_runs(report["runs"])
    json_temp = directory / "study.json.tmp"
    with json_temp.open("w", encoding="utf-8", newline="\n") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")
    json_temp.replace(directory / "study.json")
    csv_temp = directory / "study.csv.tmp"
    with csv_temp.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["runs"])
    csv_temp.replace(directory / "study.csv")


def extract_summary(summary: dict, requested_count: int) -> dict:
    if summary.get("sample_count") != requested_count:
        raise StudyError("Child sample count differs from requested Gaussian ensemble size")
    record = {"threads_used": summary.get("threads_used"), "frame_count": summary.get("frame_count")}
    for field in ("threads_used", "frame_count"):
        if not isinstance(record[field], int) or isinstance(record[field], bool) or record[field] < 1:
            raise StudyError(f"Invalid child summary {field}")
    record["simulation_elapsed_seconds"] = numeric(summary.get("elapsed_seconds"), "elapsed_seconds")
    for event in EVENTS:
        event_time = numeric(summary.get(f"{event}_time_s"), f"{event}_time_s", nullable=True)
        bracket = summary.get(f"{event}_bracket_s")
        if event_time is None:
            if bracket is not None:
                raise StudyError(f"Null {event} time has a non-null bracket")
            lower = upper = None
        else:
            if not isinstance(bracket, list) or len(bracket) != 2:
                raise StudyError(f"Observed {event} needs a two-endpoint bracket")
            lower = numeric(bracket[0], f"{event}_lower_s")
            upper = numeric(bracket[1], f"{event}_upper_s")
            if lower > upper or event_time != upper:
                raise StudyError(f"Inconsistent {event} time/bracket")
        record.update({f"{event}_time_s": event_time, f"{event}_lower_s": lower,
                       f"{event}_upper_s": upper})
    for field in ("analytic_mixing_time_s", "dkw_cdf_error_95"):
        record[field] = numeric(summary.get(field), field, nullable=field.startswith("analytic"))
    checked = summary.get("accuracy_checked_samples")
    if not isinstance(checked, int) or isinstance(checked, bool) or checked < 0:
        raise StudyError("Invalid accuracy_checked_samples")
    record["accuracy_checked_samples"] = checked
    for field in ("accuracy_max_position_m", "accuracy_max_phase_rad"):
        # A disabled check does not establish zero numerical error.
        record[field] = numeric(summary.get(field), field) if checked else None
    return record


def execute_case(args: argparse.Namespace, report: dict, output: Path,
                 run_id: str, kind: str, count: int, seed: int,
                 cadence: float | None, threads: int) -> dict:
    case_directory = output / "runs" / run_id
    case_directory.mkdir(parents=True, exist_ok=False)
    result_directory = case_directory / "outputs"
    command = [str(args.exe), "--config", str(args.config), "--samples", str(count),
               "--seed", str(seed), "--threads", str(threads), "--output", str(result_directory),
               "--no-html", "--visual-samples", "0", "--accuracy-check", "0"]
    if cadence is not None:
        command.extend(["--output-step-days", format(cadence, ".17g")])
    invocation = {"run_id": run_id, "kind": kind, "samples": count, "seed": seed,
                  "cadence_days": cadence, "requested_threads": threads,
                  "command": command, "output_directory": str(result_directory),
                  "log_path": str(case_directory / "process.log")}
    report["active_run"] = invocation
    checkpoint(output, report)
    print(f"[{len(report['runs']) + 1}/{report['planned_run_count']}] {run_id}", flush=True)
    start = time.perf_counter()
    with (case_directory / "process.log").open("w", encoding="utf-8") as log:
        log.write("Command arguments: " + json.dumps(command) + "\n")
        log.flush()
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                   timeout=args.timeout_seconds, check=False)
    wall = time.perf_counter() - start
    if completed.returncode != 0:
        raise StudyError(f"{run_id}: executable returned {completed.returncode}; inspect {invocation['log_path']}")
    summary = strict_json(result_directory / "summary.json")
    record = {**invocation, "wall_seconds": wall, **extract_summary(summary, count),
              "metrics_sha256": sha256(result_directory / "metrics.csv"),
              "initial_samples_sha256": sha256(result_directory / "initial_samples.csv")}
    report["active_run"] = None
    report["runs"].append(record)
    print(f"  completed in {wall:.3f} s; simulation {record['simulation_elapsed_seconds']:.3f} s; "
          f"mixing onset {record['mixing_time_s']!r} s", flush=True)
    return record


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    result.add_argument("--exe", required=True, type=Path, help="Compiled distribution_propagator executable")
    result.add_argument("--config", required=True, type=Path, help="Gaussian input configuration")
    result.add_argument("--output", required=True, type=Path, help="New or empty study directory")
    result.add_argument("--counts", default=[2000, 5000, 20000],
                        type=lambda text: comma_integers(text, 2, 100000000))
    result.add_argument("--seeds", default=[11, 22, 33],
                        type=lambda text: comma_integers(text, 0, (1 << 64) - 1))
    result.add_argument("--threads", default=8, type=int, help="Science-run threads; 0 uses hardware concurrency")
    result.add_argument("--cadences", type=comma_cadences,
                        help="Comma-separated snapshot cadences in days; omission preserves config cadence")
    result.add_argument("--benchmark-threads", type=lambda text: comma_integers(text, 1, 1024),
                        help="Separate fixed-seed thread sweep; first entry is reference, normally 1")
    result.add_argument("--timeout-seconds", type=float, help="Optional wall-clock limit per child invocation")
    return result


def main(argv: list[str] | None = None) -> int:
    command_parser = parser()
    args = command_parser.parse_args(argv)
    if not 0 <= args.threads <= 1024:
        command_parser.error("--threads must be between 0 and 1024")
    if args.timeout_seconds is not None and (not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0):
        command_parser.error("--timeout-seconds must be positive and finite")
    args.exe, args.config, args.output = args.exe.resolve(), args.config.resolve(), args.output.resolve()
    for path in (args.exe, args.config):
        if not path.is_file():
            command_parser.error(f"Required file does not exist: {path}")
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        command_parser.error("--output must be new or empty; unfinished and completed runs are never reused")
    args.output.mkdir(parents=True, exist_ok=True)
    cadences = args.cadences if args.cadences is not None else [None]
    benchmarks = args.benchmark_threads or []
    report = {
        "schema_version": 1, "status": "running", "started_at_utc": utc_now(),
        "executable": str(args.exe), "executable_sha256": sha256(args.exe),
        "config": str(args.config), "config_sha256": sha256(args.config),
        "counts": args.counts, "seeds": args.seeds, "cadences_days": cadences,
        "science_threads": args.threads, "benchmark_threads": benchmarks,
        "planned_run_count": len(args.counts) * len(args.seeds) * len(cadences) + len(benchmarks),
        "notes": [
            "All requested particles contribute to scientific metrics; displayed particles and HTML are disabled.",
            "Accuracy repropagation is disabled for comparable timing; null accuracy fields mean not checked.",
            "Null event times mean not observed within the configured horizon/persistence; never replace them by zero.",
            "Seed scatter is descriptive, not a confidence interval; summaries use observed events and state missing counts.",
            "Cadence null means the input configuration cadence, not a zero time interval.",
            "Thread benchmarks use the first count, seed, and cadence and compare exact metric/sample file hashes.",
        ],
        "runs": [], "active_run": None,
    }
    checkpoint(args.output, report)
    try:
        for count in args.counts:
            for cadence_index, cadence in enumerate(cadences):
                for seed in args.seeds:
                    run_id = f"science_n{count}_seed{seed}_cadence{cadence_index:02d}"
                    execute_case(args, report, args.output, run_id, "science", count, seed, cadence, args.threads)
                    checkpoint(args.output, report)
        baseline = None
        for threads in benchmarks:
            run_id = f"benchmark_n{args.counts[0]}_seed{args.seeds[0]}_threads{threads}"
            row = execute_case(args, report, args.output, run_id, "benchmark", args.counts[0],
                               args.seeds[0], cadences[0], threads)
            if baseline is None:
                baseline = row
            row["metrics_match_baseline"] = row["metrics_sha256"] == baseline["metrics_sha256"]
            row["initial_samples_match_baseline"] = row["initial_samples_sha256"] == baseline["initial_samples_sha256"]
            elapsed = row["simulation_elapsed_seconds"]
            row["simulation_speedup_vs_first"] = baseline["simulation_elapsed_seconds"] / elapsed if elapsed else None
            checkpoint(args.output, report)
            if not row["metrics_match_baseline"] or not row["initial_samples_match_baseline"]:
                raise StudyError(f"Exact thread-determinism check failed for {run_id}")
        report["status"] = "complete"
        report["finished_at_utc"] = utc_now()
        checkpoint(args.output, report)
    except (OSError, ValueError, StudyError, subprocess.TimeoutExpired, KeyboardInterrupt) as error:
        report["status"] = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        report["error"] = str(error) or type(error).__name__
        checkpoint(args.output, report)
        print(f"Study {report['status']}: {report['error']}", file=sys.stderr)
        print(f"Completed artifacts and diagnostic logs remain in {args.output}", file=sys.stderr)
        return 130 if isinstance(error, KeyboardInterrupt) else 1
    print(f"Completed {len(report['runs'])} invocations. Results: {args.output / 'study.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
