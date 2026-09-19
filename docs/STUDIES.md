# Reproducible sample, cadence, and thread studies

`tools/study.py` is a Python 3.10+ standard-library driver for the compiled simulator.
It runs the Cartesian product of requested Gaussian sample counts, seeds, and
optional snapshot cadences, then collects the results. It does not reduce the
number of scientific particles: each run disables only displayed positions,
HTML generation, and the optional tighter-integration check.

For a normal convergence study in PowerShell, from the repository root:

```powershell
python tools/study.py `
  --exe build/distribution_propagator.exe `
  --config examples/meo_ball.cfg `
  --output outputs/meo-study `
  --counts 2000,5000,20000 `
  --seeds 11,22,33 `
  --threads 8 `
  --cadences 1,0.5 `
  --benchmark-threads 1,2,4,8
```

Use the executable/configuration paths appropriate to your build; a multi-config
build may place the executable under `build/Release/`. The simulator must support
the `--visual-samples` option. `--cadences` is in days; omitting it preserves the
configuration's cadence. The propagation horizon, covariance, force model,
integration tolerances, and detection thresholds remain those in the input
configuration. Input configurations for this driver must use Gaussian sampling;
the simulator rejects `--samples` for an empirical ensemble.

This command schedules 18 scientific runs and four thread benchmarks. Default
counts are `2000,5000,20000`, default seeds `11,22,33`, and default science worker
count eight. Set `--threads 0` to request hardware concurrency. The optional
`--timeout-seconds` limits each child process separately; no timeout is imposed
by default.

For a small workflow check before committing to a large study, use the controlled
Keplerian example with small counts:

```powershell
python tools/study.py --exe build/distribution_propagator.exe `
  --config examples/meo_energy_shear.cfg --output outputs/study-smoke `
  --counts 16,32 --seeds 11,22 --threads 2 --benchmark-threads 1,2
```

The second example checks orchestration and reproducibility, not adequate
statistical sampling. The first example is a study prescription, not a claim
that those runs have already been completed. Consult the root README for
delivered measurements and the actual example filenames.

## Files, failure handling, and reproducibility

The output directory must be new or empty. Each invocation receives its own
`runs/<case>/outputs/` directory and `runs/<case>/process.log`. The log records
the exact argument list and child-process output. The driver never resumes an
unfinished output directory, changes a previous scientific run, or suppresses
an executable failure.

`study.json` and `study.csv` are refreshed after each completed invocation.
Their temporary files are replaced atomically. JSON records the requested
matrix of cases, executable/configuration paths and SHA-256 hashes, actual
commands, timings, diagnostics, event intervals, and run outputs. Overall status
is `running`, `complete`, `failed`, or `interrupted`. If a child fails, the driver
stops, retains completed results and the active invocation/log, and exits with
a nonzero status. Launch a later retry into a different output directory.

Numeric fields use explicit units: seconds for event times and timings; days
for the configured cadence. The simulation time comes from `summary.json` and
excludes output writing. Driver wall time includes process startup and output.
They measure different work and should not be interchanged in benchmarks.

There are three important null conventions:

- An event time/interval is JSON `null` when it was not observed within the
  configured horizon and persistence requirement. It is not a zero-day event.
- Accuracy differences are `null` when `accuracy_checked_samples=0`. A disabled
  accuracy check is not evidence of zero numerical error.
- A `null` cadence means “use the cadence in the configuration.”

CSV represents these nulls as empty fields. An event observed at the initial
epoch remains numeric zero, with bracket `[0,0]`. Exact event brackets are kept
as separate lower/upper columns. No confidence intervals are manufactured from
those brackets: they bound saved-epoch detection, not Monte Carlo uncertainty.

## Comparing sample counts and seeds

For each sample count/cadence pair, `seed_scatter` in JSON reports the observed
event minimum, median, maximum, mean, and sample standard deviation across
seeds. The standard deviation is `null` with fewer than two observed events.
These are descriptive statistics only. The driver reports how many seeds did
and did not observe each event; a summary over observed events excludes the
unobserved ones and cannot stand in for an unconditional distribution of event
times when some runs are censored by the horizon.

Use the same seeds when increasing counts. The simulator's deterministic
Gaussian stream preserves sample prefixes, which helps distinguish effects of
larger sample count. Different seeds supply independent replicates; the small
default set of three provides an initial stability check, not an automatic
confidence statement. Compare the spread across seeds, changes from `N` to
`2N`, and the detection-time bracket width. Refine cadence around a potential
transition instead of interpreting a daily snapshot as a precise continuous
event time.

Cadence sweeps retain the selected integration tolerances. They examine
saved-epoch detection and diagnostic stability, not just plotting resolution.
They also change integration endpoints, so use the simulator's tighter
integration checks and reference tests separately to establish a numerical
error budget. The study driver sets `--accuracy-check 0` explicitly to avoid
mixing different validation costs into timing comparisons.

See [SCIENCE.md](SCIENCE.md) for statistical sample-count reasoning, circular
metrics, and the distinction between angular coverage, phase mixing, and a
geometrically thin ribbon.

## Thread benchmarks and equality checks

`--benchmark-threads` adds separate runs using the **first sample count, first
seed, and first cadence**. Only the worker count changes. The first entry is
the reference run; put one first to obtain conventional serial speedups.

The driver computes SHA-256 hashes of `metrics.csv` and `initial_samples.csv`
for each benchmark. It requires exact equality with the first benchmark's
hashes. A mismatch marks the study failed and preserves all evidence. This
gate checks initial sampling and full-ensemble numerical diagnostics; it does
not claim equality of every propagated Cartesian state because full state
export is disabled in these timing runs.

`simulation_speedup_vs_first` is the reference simulation time divided by the
current simulation time. Wall-clock times are also retained. Small tests can be
dominated by setup, scheduler, and fixed progress-polling costs; use a workload
large enough for propagation to dominate before drawing scaling conclusions.
Run on the same machine without competing heavy workloads, retain the exact
configuration/seed, and report the equality gate along with measured speedup.
