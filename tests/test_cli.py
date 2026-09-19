"""End-to-end CLI contracts and independent numerical checks (Python stdlib only).

Run: python tests/test_cli.py path/to/distribution_propagator[.exe]
An optional unittest test name or -v may follow the executable path.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import unittest


MU = 3.986004418e14
A = 26_560_000.0
DAY = 86400.0
TAU = 2 * math.pi
if len(sys.argv) < 2:
    raise SystemExit(__doc__)
EXECUTABLE = Path(sys.argv.pop(1)).resolve()


def quantile(values, p):
    """Independent linear-interpolation sample quantile."""
    ordered = sorted(values)
    location = p * (len(ordered) - 1)
    lo = int(location)
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] + (location - lo) * (ordered[hi] - ordered[lo])


def cartesian_from_equinoctial(elements):
    """Independent Kepler-elements conversion, not the C++ longitude algorithm."""
    a, ex, ey, hx, hy, longitude = elements
    eccentricity = math.hypot(ex, ey)
    node = math.atan2(hy, hx)
    inclination = 2 * math.atan(math.hypot(hx, hy))
    periapsis_longitude = math.atan2(ey, ex)
    argument = periapsis_longitude - node
    mean_anomaly = math.remainder(longitude - periapsis_longitude, TAU)
    eccentric_anomaly = mean_anomaly
    for _ in range(30):
        delta = (eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly)
                 - mean_anomaly) / (1 - eccentricity * math.cos(eccentric_anomaly))
        eccentric_anomaly -= delta
        if abs(delta) < 1e-15:
            break
    n = math.sqrt(MU / a**3)
    ce, se = math.cos(eccentric_anomaly), math.sin(eccentric_anomaly)
    beta = math.sqrt(1 - eccentricity**2)
    radius_ratio = 1 - eccentricity * ce
    px, py = a * (ce - eccentricity), a * beta * se
    vx, vy = -a * n * se / radius_ratio, a * n * beta * ce / radius_ratio
    co, so = math.cos(node), math.sin(node)
    cw, sw = math.cos(argument), math.sin(argument)
    ci, si = math.cos(inclination), math.sin(inclination)
    u = (co*cw - so*sw*ci, so*cw + co*sw*ci, sw*si)
    v = (-co*sw - so*cw*ci, -so*sw + co*cw*ci, cw*si)
    return [px*u[k] + py*v[k] for k in range(3)] + [vx*u[k] + vy*v[k] for k in range(3)]


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def read_elements(path):
    """Read exported equinoctial states, permitting documented unit suffixes."""
    rows = read_csv(path)
    aliases = [("a_m", "a", "semimajor_axis_m"), ("ex",), ("ey",),
               ("hx",), ("hy",), ("lambda_rad", "lambda", "lm_rad", "mean_longitude_rad")]
    if not rows:
        raise AssertionError("No initial samples were exported")
    columns = []
    for names in aliases:
        match = next((name for name in names if name in rows[0]), None)
        if match is None:
            raise AssertionError(f"Missing equinoctial column {names}: {list(rows[0])}")
        columns.append(match)
    return [[float(row[column]) for column in columns] for row in rows]


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="distribution-cli-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.serial = 0

    def command(self, *args, success=True):
        result = subprocess.run([str(EXECUTABLE), *map(str, args)], capture_output=True,
                                text=True, timeout=120)
        message = f"Command: {args}\nstdout: {result.stdout}\nstderr: {result.stderr}"
        if success:
            self.assertEqual(result.returncode, 0, message)
        else:
            self.assertNotEqual(result.returncode, 0, message)
            self.assertTrue(result.stderr.strip() or result.stdout.strip(), message)
        return result

    def config(self, changes=None, raw=None):
        self.serial += 1
        path = self.root / f"case-{self.serial}.cfg"
        values = {"force_model": "kepler", "initial_type": "mean", "output_type": "mean",
                  "mu": repr(MU), "orbit_equinoctial": f"{A},0,0,0,0,0",
                  "uncertainty_coordinates": "equinoctial", "sigma": "1000,0,0,0,0,0.0001",
                  "samples": 32, "threads": 1, "visual_samples": 32, "seed": 481516,
                  "duration_days": 1.05, "output_step_days": 0.4, "phase_bins": 16,
                  "persistence": 1, "minimum_altitude_m": 1_000_000}
        for key, value in (changes or {}).items():
            if value is None:
                values.pop(key, None)
            else:
                values[key] = value
        path.write_text(raw if raw is not None else "\n".join(f"{k} = {v}" for k, v in values.items()) + "\n",
                        encoding="utf-8")
        return path

    def run_case(self, changes=None, extra=(), no_html=True):
        config = self.config(changes)
        output = self.root / f"output-{self.serial}"
        self.command("--config", config, "--output", output, *(('--no-html',) if no_html else ()), *extra)
        return output, json.loads((output / "run.json").read_text()), json.loads((output / "summary.json").read_text())

    def reject(self, config, *extra):
        output = self.root / f"invalid-{self.serial}"
        self.command("--config", config, "--output", output, *extra, success=False)
        self.assertTrue(not output.exists() or not any(output.iterdir()),
                        "Invalid runs must not leave plausible output artifacts")

    def test_help_example_and_cli_overrides(self):
        self.assertIn("--config", self.command("--help").stdout)
        example = self.root / "example.cfg"
        self.command("--write-example", example)
        self.assertIn("force_model", example.read_text())
        output, run, summary = self.run_case(extra=("--samples", 12, "--threads", 2, "--seed", 17,
            "--duration-days", 0.13, "--output-step-days", 0.05, "--accuracy-check", 2, "--export-states"))
        self.assertEqual(summary["sample_count"], 12)
        self.assertEqual(summary["frame_count"], 4)
        self.assertEqual(summary["threads_used"], 2)
        self.assertEqual(run["metadata"]["seed"], 17)
        self.assertAlmostEqual(run["frames"][-1]["time_s"], .13 * DAY, places=8)
        self.assertEqual(len(read_csv(output / "states.csv")), 12 * 4)
        self.assertLess(summary["accuracy_max_position_m"], .01)
        self.assertLess(summary["accuracy_max_phase_rad"], 1e-10)

    def test_output_contract_nonintegral_cadence_and_html(self):
        output, run, summary = self.run_case(no_html=False)
        for name in ("summary.json", "run.json", "metrics.csv", "initial_samples.csv", "visualization.html"):
            self.assertTrue((output / name).is_file(), name)
        self.assertFalse((output / "states.csv").exists())
        self.assertEqual(run["schema_version"], 1)
        self.assertEqual([frame["time_s"] for frame in run["frames"]], [0, .4*DAY, .8*DAY, 1.05*DAY])
        self.assertEqual(len(read_csv(output / "metrics.csv")), 4)
        self.assertEqual(len(read_elements(output / "initial_samples.csv")), 32)
        for key in ("coverage_time_s", "mixing_time_s", "analytic_mixing_time_s", "coverage_bracket_s",
                    "mixing_bracket_s", "sample_count", "frame_count", "threads_used", "elapsed_seconds",
                    "accuracy_max_position_m", "accuracy_max_phase_rad"):
            self.assertIn(key, summary)
        html = (output / "visualization.html").read_text(encoding="utf-8")
        self.assertNotIn("__DISTRIBUTION_DATA__", html)
        self.assertIn('id="simulation-data"', html)
        # Overwriting an existing result must fail without changing any bytes.
        before = {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()}
        config = self.config()
        self.command("--config", config, "--output", output, success=False)
        self.assertEqual(before, {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()})

    def test_deterministic_across_thread_counts(self):
        first, r1, s1 = self.run_case({"samples": 64, "threads": 1}, extra=("--export-states",))
        second, r4, s4 = self.run_case({"samples": 64, "threads": 4}, extra=("--export-states",))
        self.assertEqual(s1["threads_used"], 1)
        self.assertEqual(s4["threads_used"], 4)
        self.assertEqual(r1["frames"], r4["frames"])
        for name in ("initial_samples.csv", "metrics.csv", "states.csv"):
            self.assertEqual((first / name).read_bytes(), (second / name).read_bytes(), name)
        self.assertFalse((first / "visualization.html").exists())
        third, _, _ = self.run_case({"samples": 64, "seed": 481517})
        self.assertNotEqual((first / "initial_samples.csv").read_bytes(), (third / "initial_samples.csv").read_bytes())

    def test_empirical_ensemble_equal_weights_and_exact_kepler_shear(self):
        # Asymmetric, nonuniform phase offsets exercise equal weights and quantiles.
        elements = [[A + da, 0, 0, 0, 0, phase] for da, phase in
                    [(-300000, -.012), (-110000, .005), (-20000, -.004),
                     (30000, .019), (90000, -.006), (190000, .002), (310000, .025)]]
        csv_path = self.root / "posterior.csv"
        csv_path.write_text("# absolute equinoctial states, every row has equal weight\n" +
                            "\n".join(",".join(map(repr, row)) for row in elements) + "\n")
        output, run, summary = self.run_case({"sigma": None, "empirical_samples_csv": csv_path.name,
            "samples": 99, "duration_days": 18.5, "output_step_days": 6, "visual_samples": 20})
        self.assertEqual(summary["sample_count"], len(elements), "Empirical rows must not be resampled to configured samples")
        self.assertEqual(read_elements(output / "initial_samples.csv"), elements)
        nominal_rate = math.sqrt(MU / A**3)
        for frame in run["frames"]:
            time = frame["time_s"]
            phases = [e[5] + (math.sqrt(MU / e[0]**3) - nominal_rate) * time for e in elements]
            metrics = frame["metrics"]
            self.assertAlmostEqual(metrics["phase_sigma_rad"], statistics.stdev(phases), delta=2e-10)
            self.assertAlmostEqual(metrics["phase_q95_width_rad"], quantile(phases,.975)-quantile(phases,.025), delta=5e-10)
            resultants = [abs(sum(complex(math.cos(k*p),math.sin(k*p)) for p in phases))/len(phases) for k in range(1,5)]
            for actual, expected in zip(metrics["resultants"], resultants):
                self.assertAlmostEqual(actual, expected, delta=3e-10)
            histogram = [0] * 16
            for phase in phases:
                histogram[min(15, int((phase % TAU) / TAU * 16))] += 1
            self.assertEqual(metrics["histogram"], histogram)
            self.assertEqual(sum(metrics["histogram"]), len(elements))
            wrapped = sorted(p % TAU for p in phases)
            gap = max([wrapped[i+1]-wrapped[i] for i in range(len(wrapped)-1)] + [wrapped[0]+TAU-wrapped[-1]])
            self.assertAlmostEqual(metrics["max_gap_deg"], math.degrees(gap), delta=3e-8)
            self.assertEqual(len(frame["positions_m"]), len(elements))
            for point, initial in zip(frame["positions_m"], elements):
                propagated = initial.copy()
                propagated[5] += math.sqrt(MU / initial[0]**3) * time
                expected = cartesian_from_equinoctial(propagated)[:3]
                self.assertLess(math.dist(point, expected), .01)

    def test_mixed_units_correlated_rtn_covariance(self):
        sigmas = [400, 250, 80, .03, .02, .008]
        covariance = [[0.0] * 6 for _ in range(6)]
        for i in range(6):
            covariance[i][i] = sigmas[i]**2
        for i, j, rho in [(0,4,.7),(1,3,-.45)]:
            covariance[i][j] = covariance[j][i] = rho*sigmas[i]*sigmas[j]
        cov_path = self.root / "rtn_covariance.csv"
        cov_path.write_text("\n".join(",".join(map(repr,row)) for row in covariance))
        output, _, summary = self.run_case({"sigma": None, "covariance_csv": cov_path.name,
            "uncertainty_coordinates": "rtn", "samples": 4096, "visual_samples": 16,
            "duration_days": .000001, "output_step_days": .000001, "threads": 4})
        self.assertEqual(summary["sample_count"], 4096)
        nominal = [A, 0, 0, 0, math.sqrt(MU/A), 0]
        deviations = [[state[i]-nominal[i] for i in range(6)] for state in
                      map(cartesian_from_equinoctial,read_elements(output / "initial_samples.csv"))]
        columns = list(zip(*deviations))
        for values, sigma in zip(columns, sigmas):
            self.assertAlmostEqual(statistics.stdev(values)/sigma, 1, delta=.06)
            self.assertLess(abs(statistics.mean(values)), 5*sigma/math.sqrt(len(values)))
        for i, j, target in [(0,4,.7),(1,3,-.45)]:
            xi, xj = columns[i], columns[j]
            mi, mj = statistics.mean(xi), statistics.mean(xj)
            covariance_observed = sum((a-mi)*(b-mj) for a,b in zip(xi,xj))/(len(xi)-1)
            correlation = covariance_observed/(statistics.stdev(xi)*statistics.stdev(xj))
            self.assertAlmostEqual(correlation, target, delta=.04)

    def test_cartesian_empirical_and_rtn_offset_semantics(self):
        speed = math.sqrt(MU / A)
        absolute = [[A+50,0,0,0,speed+.002,0], [A-80,20,3,.001,speed-.004,.002]]
        offsets = [[s[i]-[A,0,0,0,speed,0][i] for i in range(6)] for s in absolute]
        for mode, rows in [("cartesian", absolute), ("rtn", offsets)]:
            with self.subTest(mode=mode):
                source = self.root / f"{mode}.csv"
                source.write_text("\n".join(",".join(map(repr,row)) for row in rows))
                output, _, _ = self.run_case({"sigma": None,"empirical_samples_csv": source.name,
                    "uncertainty_coordinates": mode,"duration_days": .001,"output_step_days": .001})
                recovered = list(map(cartesian_from_equinoctial,read_elements(output / "initial_samples.csv")))
                for actual, expected in zip(recovered, absolute):
                    self.assertLess(math.dist(actual[:3],expected[:3]),1e-5)
                    self.assertLess(math.dist(actual[3:],expected[3:]),1e-8)

    def test_zero_covariance_is_valid_and_does_not_mix(self):
        _, run, summary = self.run_case({"sigma": "0,0,0,0,0,0", "samples": 4})
        self.assertIsNone(summary["coverage_time_s"])
        self.assertIsNone(summary["mixing_time_s"])
        self.assertIsNone(summary["analytic_mixing_time_s"])
        for frame in run["frames"]:
            self.assertEqual(frame["metrics"]["phase_sigma_rad"],0)
            self.assertEqual(frame["metrics"]["phase_q95_width_rad"],0)
            self.assertEqual(frame["metrics"]["resultants"],[1,1,1,1])

    def test_initial_uniform_phase_meets_sustained_criteria_at_zero(self):
        source = self.root / "uniform.csv"
        source.write_text("\n".join(f"{A},0,0,0,0,{(i+.13)*TAU/64:.17g}" for i in range(64)))
        _, run, summary = self.run_case({"sigma": None,"empirical_samples_csv": source.name,
            "coverage_max_gap_deg": 6,"persistence": 2,"visual_samples": 8})
        self.assertEqual(summary["coverage_time_s"],0)
        self.assertEqual(summary["mixing_time_s"],0)
        self.assertEqual(summary["coverage_bracket_s"],[0,0])
        self.assertEqual(summary["mixing_bracket_s"],[0,0])
        for frame in run["frames"]:
            self.assertTrue(frame["metrics"]["coverage"])
            self.assertTrue(frame["metrics"]["mixed"])
            self.assertEqual(frame["metrics"]["histogram"],[4]*16)
            self.assertEqual(len(frame["positions_m"]),8)

    def test_invalid_random_draw_is_rejected_without_clipping_or_redraw(self):
        config = self.config({"sigma":"40000000,0,0,0,0,0","samples":32})
        output = self.root / "invalid-random"
        result = self.command("--config",config,"--output",output,success=False)
        self.assertIn("sampled orbit",result.stdout+result.stderr)
        self.assertTrue(not output.exists() or not any(output.iterdir()))

    def test_rejects_invalid_configuration_before_artifacts(self):
        changes = [{"samples": -1},{"samples": 1},{"duration_days": -2},{"output_step_days": 0},
                   {"sigma": "-1,0,0,0,0,0"},{"sigma": "1,2,3"},{"seed": "5junk"},
                   {"duration_days": "nan"},{"force_model": "made_up"},{"initial_type": "made_up"},
                   {"uncertainty_coordinates": "made_up"},{"unknown_key": 1},{"phase_bins": 2},
                   {"persistence": 0},{"sigma": None,"uncertainty_coordinates": "rtn"},
                   {"max_memory_mb": .000001}]
        for values in changes:
            with self.subTest(values=values):
                self.reject(self.config(values))
        for raw in ["samples 32\n", "samples=32\nsamples=64\n", "=32\n", "samples=\n"]:
            with self.subTest(raw=raw):
                self.reject(self.config(raw=raw))

    def test_rejects_bad_covariance_and_invalid_empirical_orbits(self):
        for name, content in [
            ("nonsymmetric", "1,1,0,0,0,0\n0,1,0,0,0,0\n0,0,1,0,0,0\n0,0,0,1,0,0\n0,0,0,0,1,0\n0,0,0,0,0,1"),
            ("indefinite", "1,2,0,0,0,0\n2,1,0,0,0,0\n0,0,1,0,0,0\n0,0,0,1,0,0\n0,0,0,0,1,0\n0,0,0,0,0,1"),
            ("wrong_shape", "1,0,0,0,0,0\n"),
        ]:
            with self.subTest(covariance=name):
                path = self.root / f"{name}.csv"
                path.write_text(content)
                self.reject(self.config({"sigma": None,"covariance_csv": path.name}))
        for name, content in [
            ("hyperbolic", f"{A},0,0,0,0,0\n{A},1.1,0,0,0,0\n"),
            ("low_perigee", f"{A},0,0,0,0,0\n6700000,0,0,0,0,0\n"),
            ("single_row", f"{A},0,0,0,0,0\n"),
            ("malformed", f"{A},0,0,0,0,0\ninvalid,row\n"),
        ]:
            with self.subTest(empirical=name):
                path = self.root / f"{name}.csv"
                path.write_text(content)
                self.reject(self.config({"sigma": None,"empirical_samples_csv": path.name}))


if __name__ == "__main__":
    unittest.main()
