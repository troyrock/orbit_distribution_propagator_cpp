#!/usr/bin/env python3
"""Render native C++ coverage-grid results as a self-contained HTML explorer.

Only the Python standard library is required. Plotly's local JavaScript bundle
is embedded in the report, so viewing never requires a network connection.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
AXES = (
    ("perigee_altitudes_km", "perigee_altitude_km"),
    ("inclinations_deg", "inclination_deg"),
    ("eccentricities", "eccentricity"),
    ("position_sigmas_km", "position_sigma_km"),
    ("velocity_sigmas_m_s", "velocity_sigma_m_s"),
)


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def validate_data(data: Any) -> None:
    """Reject ambiguous grids; retain explicit invalid or missing-node gaps."""
    if not isinstance(data, dict) or not isinstance(data.get("metadata"), dict):
        raise ValueError("Input must contain a metadata object.")
    meta = data["metadata"]
    if not isinstance(data.get("runs"), list):
        raise ValueError("Input must contain a runs array.")
    if meta.get("initial_view", "grid") not in ("grid", "reference"):
        raise ValueError("metadata.initial_view must be grid or reference.")
    if meta.get("initial_view") == "reference" and "reference_surface" not in data:
        raise ValueError("The initial reference view requires a reference_surface.")
    if not isinstance(meta.get("samples"), int) or isinstance(meta["samples"], bool) or meta["samples"] < 1:
        raise ValueError("metadata.samples must be a positive integer.")
    for field, _ in AXES:
        axis = meta.get(field)
        if not isinstance(axis, list) or not axis or not all(finite_number(v) for v in axis):
            raise ValueError(f"metadata.{field} must be a nonempty array of finite numbers.")
        if any(b <= a for a, b in zip(axis, axis[1:])):
            raise ValueError(f"metadata.{field} must be strictly increasing.")
    if meta["position_sigmas_km"][0] <= 0 or meta["velocity_sigmas_m_s"][0] <= 0:
        raise ValueError("Uncertainty axes must be positive for logarithmic plotting.")
    if not 0 <= meta["eccentricities"][0] <= meta["eccentricities"][-1] < 1:
        raise ValueError("Eccentricity must lie in [0, 1).")
    if "earth_radius_km" in meta and (not finite_number(meta["earth_radius_km"]) or meta["earth_radius_km"] <= 0):
        raise ValueError("metadata.earth_radius_km must be positive and finite.")
    if meta.get("altitude_interpolation", "linear_altitude") not in ("linear_altitude", "log_geocentric_perigee_radius"):
        raise ValueError("Unsupported metadata.altitude_interpolation policy.")
    if meta.get("altitude_interpolation") == "log_geocentric_perigee_radius":
        if "earth_radius_km" not in meta:
            raise ValueError("Log-radius interpolation requires explicit metadata.earth_radius_km.")
        if any(not finite_number(value + meta["earth_radius_km"]) or value + meta["earth_radius_km"] <= 0
               for value in meta["perigee_altitudes_km"]):
            raise ValueError("Log-radius interpolation requires positive geocentric perigee radii that are finite.")
    seen, run_ids = set(), set()
    for row, run in enumerate(data["runs"]):
        if not isinstance(run, dict):
            raise ValueError(f"Run {row} must be an object.")
        run_id = run.get("run_id")
        if not isinstance(run_id, str) or not run_id or run_id in run_ids:
            raise ValueError(f"Run {row} must have a unique nonempty run_id.")
        run_ids.add(run_id)
        indices = []
        for field, value_field in AXES:
            value = run.get(value_field)
            if not finite_number(value):
                raise ValueError(f"Run {run_id}: {value_field} must be finite.")
            matches = [j for j, item in enumerate(meta[field]) if math.isclose(item, value, rel_tol=1e-9, abs_tol=1e-9)]
            if len(matches) != 1:
                raise ValueError(f"Run {run_id}: {value_field} is not an unambiguous grid coordinate.")
            indices.append(matches[0])
        location = tuple(indices)
        if location in seen:
            raise ValueError(f"Run {run_id}: duplicate grid coordinates.")
        seen.add(location)
        if run.get("samples", meta["samples"]) != meta["samples"]:
            raise ValueError(f"Run {run_id}: sample count differs from metadata.samples.")
        if run.get("status") not in ("ok", "earth_intersection", "mean_earth_intersection", "no_phase_shear", "unobserved"):
            raise ValueError(f"Run {run_id}: unsupported status.")
        t = run.get("coverage_time_days")
        if run["status"] == "ok" and (not finite_number(t) or t <= 0):
            raise ValueError(f"Run {run_id}: successful coverage requires a positive finite time.")
        for field in ("coverage_time_days", "coverage_lower_days", "coverage_upper_days"):
            value = run.get(field)
            if value is not None and (not finite_number(value) or value < 0):
                raise ValueError(f"Run {run_id}: {field} must be null or nonnegative and finite.")
        lo, hi = run.get("coverage_lower_days"), run.get("coverage_upper_days")
        if lo is not None and hi is not None and lo > hi:
            raise ValueError(f"Run {run_id}: onset bracket is reversed.")
        if run["status"] == "ok" and ((lo is not None and lo > t) or (hi is not None and hi < t)):
            raise ValueError(f"Run {run_id}: coverage time lies outside its onset bracket.")
        for field in ("particles_below_300km", "particles_below_500km"):
            count = run.get(field, 0)
            if not isinstance(count, int) or isinstance(count, bool) or not 0 <= count <= meta["samples"]:
                raise ValueError(f"Run {run_id}: {field} must be a count within the ensemble size.")
        perigee = run.get("min_initial_perigee_km")
        if perigee is not None and not finite_number(perigee):
            raise ValueError(f"Run {run_id}: min_initial_perigee_km must be finite or null.")
    if "reference_surface" in data:
        reference = data["reference_surface"]
        if not isinstance(reference, dict) or not isinstance(reference.get("label"), str) or not reference["label"].strip():
            raise ValueError("reference_surface must have a nonempty label.")
        if "reference_surface" in reference:
            raise ValueError("A reference surface cannot contain another reference surface.")
        if not isinstance(reference.get("description", ""), str):
            raise ValueError("reference_surface.description must be a string.")
        validate_data(reference)
        reference_meta = reference["metadata"]
        for field, _ in AXES[:3]:
            if len(reference_meta[field]) != 1:
                raise ValueError("The reference surface must describe exactly one nominal orbit.")
            value = reference_meta[field][0]
            if not meta[field][0] <= value <= meta[field][-1]:
                raise ValueError(f"Reference orbit {field} must lie within the explorer range.")
        for field in ("samples", "position_sigmas_km", "velocity_sigmas_m_s"):
            if reference_meta[field] != meta[field]:
                raise ValueError(f"Reference surface {field} must match the main grid.")
        if reference_meta.get("earth_radius_km", 6378.137) != meta.get("earth_radius_km", 6378.137):
            raise ValueError("Reference surface earth_radius_km must match the main grid.")


def find_plotly(explicit: Path | None) -> Path:
    if explicit is not None:
        if not explicit.is_file():
            raise FileNotFoundError(f"Plotly JavaScript not found: {explicit}")
        return explicit
    candidates = [
        REPO.parent / "outputs/coverage_surface_20000/.deps/plotly/package_data/plotly.min.js",
        REPO / "external/plotly.min.js",
    ]
    # Discovery only: importing Plotly is not required to render the report.
    candidates.extend(Path(p) / "plotly/package_data/plotly.min.js" for p in sys.path if p)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("No local Plotly JavaScript bundle found. Supply --plotly-js /path/to/plotly.min.js.")


def render(data: dict[str, Any], plotly_path: Path, template_path: Path | None = None) -> str:
    """Return an offline report; escape data so embedded HTML cannot execute."""
    validate_data(data)
    template = (template_path or REPO / "web/coverage_explorer.html").read_text(encoding="utf-8")
    data_token, plotly_token = "__COVERAGE_EXPLORER_DATA__", "__COVERAGE_EXPLORER_PLOTLY__"
    if template.count(data_token) != 1 or template.count(plotly_token) != 1:
        raise ValueError("Template must contain exactly one data token and one plotting-library token.")
    payload = json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    payload = payload.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    payload = payload.replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    library = plotly_path.read_text(encoding="utf-8").replace("</script", "<\\/script")
    return template.replace(plotly_token, library).replace(data_token, payload)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON result grid from the C++ coverage study")
    parser.add_argument("--output", type=Path, required=True, help="Standalone HTML output path")
    parser.add_argument("--plotly-js", type=Path, help="Local plotly.min.js bundle")
    parser.add_argument("--template", type=Path, help="Alternative explorer HTML template")
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text(encoding="utf-8-sig"))
        report = render(data, find_plotly(args.plotly_js), args.template)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8", newline="\n")
    except (OSError, ValueError) as error:
        parser.exit(2, f"Error: {error}\n")
    print(f"Wrote {args.output.resolve()} ({len(data['runs']):,} recorded grid nodes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
