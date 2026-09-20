#!/usr/bin/env python3
"""Independently audit coverage records and withheld log-time interpolation.

Only Python's standard library is used. Missing or masked support nodes never
produce a prediction. --allow-partial permits a provisional checkpoint audit;
it does not claim that unevaluated study cases have passed validation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
import math
from pathlib import Path


AXES = ("perigee_altitudes_km", "inclinations_deg", "eccentricities",
        "position_sigmas_km", "velocity_sigmas_m_s")
FIELDS = ("perigee_altitude_km", "inclination_deg", "eccentricity",
          "position_sigma_km", "velocity_sigma_m_s")
STATUSES = {"ok", "earth_intersection", "mean_earth_intersection",
            "unobserved", "no_phase_shear"}
MASKED = {"earth_intersection", "mean_earth_intersection"}
EVENT_FIELDS = ("coverage_time_days", "coverage_lower_days", "coverage_upper_days",
                "coverage_confirmation_days", "onset_max_gap_deg")
LIMIT = 0.05


def finite(value, context):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{context}: expected a finite number")
    return float(value)


def integer(value, context, lower=0, upper=None):
    if isinstance(value, bool) or not isinstance(value, int) or value < lower or (
            upper is not None and value > upper):
        raise ValueError(f"{context}: invalid integer {value!r}")
    return value


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12)


def axis_index(axis, value):
    matches = [index for index, node in enumerate(axis) if close(node, value)]
    if len(matches) != 1:
        raise ValueError(f"Value {value:.17g} does not uniquely match its declared grid axis")
    return matches[0]


def validate_record(row, samples, seed, context):
    if row.get("samples") != samples or row.get("seed") != seed:
        raise ValueError(f"{context}: particle count or seed differs from study metadata")
    integer(row["samples"], context + " samples", 20000, 20000)
    integer(row["seed"], context + " seed", 0, 2**64 - 1)
    status = row.get("status")
    if status not in STATUSES:
        raise ValueError(f"{context}: unknown status {status!r}")
    tail = [integer(row.get(key), context + " " + key, 0, samples) for key in
            ("earth_intersections", "particles_below_300km", "particles_below_500km",
             "particles_below_1000km")]
    if tail != sorted(tail):
        raise ValueError(f"{context}: perigee-tail counts are not nested")
    minimum = finite(row.get("min_initial_perigee_km"), context + " initial perigee")
    epochs = integer(row.get("evaluated_epochs"), context + " epoch count")
    cadence = finite(row.get("cadence_days"), context + " cadence")
    horizon = finite(row.get("evaluated_horizon_days"), context + " horizon")
    if cadence < 0 or horizon < 0:
        raise ValueError(f"{context}: negative cadence or horizon")
    if status == "earth_intersection":
        if not tail[0] or minimum >= 0 or epochs != 0:
            raise ValueError(f"{context}: initial Earth-intersection census is inconsistent")
        if row.get("min_mean_perigee_km") is not None:
            raise ValueError(f"{context}: masked initial case must not claim initialized mean perigee")
    else:
        mean_minimum = finite(row.get("min_mean_perigee_km"), context + " mean perigee")
        if tail[0] or minimum < 0:
            raise ValueError(f"{context}: initially Earth-crossing case was not masked")
        if (mean_minimum < 0) != (status == "mean_earth_intersection"):
            raise ValueError(f"{context}: mean Earth-intersection status is inconsistent")
    if status != "ok":
        if any(row.get(key) is not None for key in EVENT_FIELDS):
            raise ValueError(f"{context}: unavailable coverage must have null event fields")
        return
    time, lower, upper, confirmation, gap = [finite(row.get(key), context + " " + key)
                                            for key in EVENT_FIELDS]
    tolerance = 1e-9 + 1e-11 * max(1, time, confirmation, cadence)
    if epochs < 3 or cadence <= 0 or lower < 0 or upper < lower or time < 0:
        raise ValueError(f"{context}: invalid successful event bounds or persistence")
    if abs(time - upper) > tolerance:
        raise ValueError(f"{context}: coverage time differs from the upper onset bound")
    if upper - lower > cadence + tolerance:
        raise ValueError(f"{context}: onset bracket exceeds one cadence")
    if abs(confirmation - (upper + 2 * cadence)) > tolerance:
        raise ValueError(f"{context}: confirmation is not two cadences after onset")
    if abs(horizon - confirmation) > tolerance:
        raise ValueError(f"{context}: evaluated horizon differs from confirmation")
    if not 0 <= gap <= 5 + 1e-10:
        raise ValueError(f"{context}: reported onset gap violates five-degree coverage")


def validate_document(document, label, allow_partial=False):
    metadata, rows = document.get("metadata"), document.get("runs")
    if not isinstance(metadata, dict) or not isinstance(rows, list):
        raise ValueError(f"{label}: expected metadata object and runs list")
    samples = integer(metadata.get("samples"), label + " samples", 20000, 20000)
    seed = integer(metadata.get("seed"), label + " seed", 0, 2**64 - 1)
    axes = []
    for name in AXES:
        raw = metadata.get(name)
        if not isinstance(raw, list) or not raw:
            raise ValueError(f"{label}: missing or empty {name}")
        axis = [finite(value, label + " " + name) for value in raw]
        if any(axis[i] >= axis[i + 1] or close(axis[i], axis[i + 1])
               for i in range(len(axis) - 1)):
            raise ValueError(f"{label}: {name} must be distinct and strictly increasing")
        axes.append(axis)
    if (axes[0][0] < 0 or axes[1][0] < 0 or axes[1][-1] > 90 or axes[2][0] < 0 or
            axes[2][-1] >= 1 or axes[3][0] <= 0 or axes[4][0] <= 0):
        raise ValueError(f"{label}: axes outside the supported orbit/uncertainty domain")
    expected = math.prod(len(axis) for axis in axes)
    if metadata.get("expected_cases") != expected or metadata.get("completed_cases") != len(rows):
        raise ValueError(f"{label}: metadata case counts do not match axes and records")
    complete = metadata.get("complete")
    if not isinstance(complete, bool):
        raise ValueError(f"{label}: complete must be Boolean")
    if complete and len(rows) != expected:
        raise ValueError(f"{label}: complete study does not contain the full Cartesian grid")
    if not complete and not allow_partial:
        raise ValueError(f"{label}: incomplete study requires --allow-partial")
    identifiers, lookup = set(), {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("run_id"), str) or not row["run_id"]:
            raise ValueError(f"{label}: invalid run record or identifier")
        context = label + "/" + row["run_id"]
        if row["run_id"] in identifiers:
            raise ValueError(f"{context}: duplicate run identifier")
        identifiers.add(row["run_id"])
        key = tuple(axis_index(axis, finite(row.get(field), context + " " + field))
                    for axis, field in zip(axes, FIELDS))
        if key in lookup:
            raise ValueError(f"{context}: duplicate parameter point")
        validate_record(row, samples, seed, context)
        lookup[key] = row
    # Membership and uniqueness plus the Cartesian count establish completeness.
    summary = {"complete": complete, "expected_cases": expected, "validated_cases": len(rows),
               "missing_cases": expected - len(rows), "samples_per_case": samples, "seed": seed,
               "status_counts": dict(sorted(Counter(row["status"] for row in rows).items()))}
    return summary, axes, lookup


def brackets(axis, value):
    for index, node in enumerate(axis):
        if close(node, value):
            return [(index, 1.0)]
    for index in range(len(axis) - 1):
        if axis[index] < value < axis[index + 1]:
            weight = (value - axis[index]) / (axis[index + 1] - axis[index])
            return [(index, 1 - weight), (index + 1, weight)]
    raise ValueError("Withheld orbit lies outside the measured interpolation domain")


def predict_withheld(row, axes, lookup):
    result = {"run_id": row["run_id"], **{key: row[key] for key in FIELDS},
              "withheld_status": row["status"],
              "measured_time_days": row.get("coverage_time_days"),
              "predicted_time_days": None, "relative_error": None,
              "absolute_relative_error": None}
    if row["status"] != "ok":
        result["status"] = "withheld_masked" if row["status"] in MASKED else "withheld_unavailable"
        return result
    try:
        uncertainty = tuple(axis_index(axes[k], row[FIELDS[k]]) for k in (3, 4))
    except ValueError:
        result["status"] = "missing_uncertainty_axis"
        return result
    try:
        orbit = [brackets(axes[k], row[FIELDS[k]]) for k in range(3)]
    except ValueError:
        result["status"] = "outside_orbit_grid"
        return result
    missing, invalid, support = [], [], []
    for corner in itertools.product(*orbit):
        key = tuple(index for index, _ in corner) + uncertainty
        weight = math.prod(weight for _, weight in corner)
        node = lookup.get(key)
        if node is None:
            missing.append({field: axes[k][key[k]] for k, field in enumerate(FIELDS)})
        elif node["status"] != "ok" or node["coverage_time_days"] <= 0:
            invalid.append({"run_id": node["run_id"], "status": node["status"],
                            "coverage_time_days": node["coverage_time_days"]})
        else:
            support.append({"run_id": node["run_id"], "weight": weight,
                            "time_days": node["coverage_time_days"]})
    result["missing_corners"], result["invalid_corners"] = missing, invalid
    if missing or invalid:
        result["status"] = ("missing_and_invalid_corners" if missing and invalid else
                            "missing_corners" if missing else "invalid_corners")
        return result
    if row["coverage_time_days"] <= 0:
        result["status"] = "nonpositive_withheld_time"
        return result
    prediction = math.exp(math.fsum(node["weight"] * math.log(node["time_days"])
                                    for node in support))
    error = prediction / row["coverage_time_days"] - 1
    result.update(status="compared", predicted_time_days=prediction, relative_error=error,
                  absolute_relative_error=abs(error), support=support,
                  within_five_percent=abs(error) <= LIMIT)
    return result


def percentile(values, fraction):
    ordered = sorted(values)
    at = fraction * (len(ordered) - 1)
    lower = int(math.floor(at))
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (at - lower) * (ordered[upper] - ordered[lower])


def validate_grid(study, holdouts, allow_partial=False):
    measured, axes, lookup = validate_document(study, "measured", allow_partial)
    withheld, _, _ = validate_document(holdouts, "holdouts", allow_partial)
    if measured["seed"] != withheld["seed"]:
        raise ValueError("Measured and withheld studies use different random seeds")
    for name in ("model_scope", "uncertainty_convention", "altitude_definition", "earth_radius_km",
                 "raan_deg", "argument_of_perigee_deg", "mean_anomaly_deg", "coverage_definition"):
        if study["metadata"].get(name) != holdouts["metadata"].get(name):
            raise ValueError(f"Measured and withheld metadata differ for {name}")
    comparisons = [predict_withheld(row, axes, lookup) for row in holdouts["runs"]]
    errors = [row["absolute_relative_error"] for row in comparisons if row["status"] == "compared"]
    statistics = {"n": len(errors), "max_relative_error": max(errors) if errors else None,
                  "p95_relative_error": percentile(errors, 0.95) if errors else None,
                  "rms_relative_error": math.sqrt(math.fsum(e * e for e in errors) / len(errors))
                  if errors else None,
                  "exceeding_five_percent": sum(error > LIMIT for error in errors)}
    gate = "not_evaluated" if not errors else "failed" if max(errors) > LIMIT else "passed"
    partial = not (measured["complete"] and withheld["complete"])
    return {"schema_version": 1, "data_integrity_passed": True,
            "passed": gate == "passed", "provisional": partial,
            "status": "partial" if partial else "complete", "measured": measured,
            "holdouts": withheld, "interpolation_method": "Trilinear interpolation of log time",
            "relative_error_limit": LIMIT, "interpolation_gate": gate,
            "statistics": statistics,
            "comparison_status_counts": dict(sorted(Counter(row["status"] for row in comparisons).items())),
            "comparisons": comparisons,
            "scope_note": "Masked or unsupported holdouts are reported without predictions and excluded from error statistics. Passing applies only to the directly compared holdouts; partial results are provisional."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--holdouts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    try:
        report = validate_grid(json.loads(args.study.read_text(encoding="utf-8")),
                               json.loads(args.holdouts.read_text(encoding="utf-8")),
                               args.allow_partial)
    except (ValueError, KeyError, TypeError, OSError) as error:
        report = {"schema_version": 1, "passed": False, "data_integrity_passed": False,
                  "status": "invalid_data", "errors": [str(error)]}
    report["inputs"] = {"study": str(args.study.resolve()), "holdouts": str(args.holdouts.resolve())}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({key: value for key, value in report.items() if key != "comparisons"}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
