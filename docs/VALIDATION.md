# Validation and measured results — 2026-09-19

All physical experiments below use the declared Earth monopole/J2/J2-squared
model. They are conditional uncertainty experiments, not real-object forecasts.
The native dependency was clean at revision
`e653cd40e5e057cb89395afd8a9e3012ad59b12c`.

## Test evidence

The final CTest configuration contains **13 tests**; assertions remain active
in Release. Builds use GCC 15.2 and MSVC 19.51 through CMake/Ninja.
The final gates passed **13/13 with no skipped tests on both compilers**
(34.47 seconds GCC; 43.92 seconds MSVC on the shared machine).

- Math: 10 groups covering independent orbit geometry, energy/angular momentum,
  covariance factors, phase branches, lobes, uniform distributions and events.
- Configuration: 9 groups including covariance moments, RTN conventions,
  empirical branches, invalid inputs and memory preflight.
- Backend: 20-year analytic Kepler, initial mean/osculating inversion, native
  RK4 comparisons including eccentric/retrograde orbits, tolerance/cadence
  convergence and exact parallel agreement.
- Independent Java DSST: 32 comparisons at 0/1/30/365 days, J2 and J2-squared,
  mean/osculating inputs and outputs. Maximum position difference **8.781 mm**;
  declared new gate **20 mm**. [Fixture provenance](../tests/data/README.md).
- Seven existing upstream Java fixture groups pass with their original bounds.
  No upstream source, fixture or tolerance was modified.
- CLI: 11 executable-level tests; study driver: 7 tests covering dimensions,
  hashes, null/zero conventions and failed-run diagnostics.
- Browser QA: real 121-epoch data and synthetic data at desktop/mobile sizes;
  playback, camera controls, initial cloud, early banana, final ribbon and
  exact large seeds. No external requests or JavaScript errors.

## Broad initial ball demonstration

Run `build/distribution_propagator.exe --config examples/meo_ball.cfg --output
outputs/meo-ball`. The orbit has a=26,560 km, e=0.02 and i=55 degrees. The
osculating input has independent RTN position sigma=100 km per axis and velocity
sigma=1 m/s per axis. This broad input is illustrative, not a measured posterior.
Seed=20260919, N=5,000, eight workers, 121 epochs over 60 days.

| Quantity | Measured result |
|---|---:|
| Initial radial / along-track / cross-track sigma | 99.70 / 99.32 / 100.93 km |
| Same sigmas after 0.5 days | 141.26 / 1,995.45 / 100.60 km |
| First sustained coverage bracket | 7.0–7.5 days |
| Central-95% unwrapped width reaches one turn | 11.0–11.5 days |
| First sustained mixing bracket | 15.5–16.0 days |
| Initial Gaussian Kepler-shear estimate | 16.53 days |
| Tighter integration check on eight samples | max 2.99 micrometres position difference |
| Sampling/propagation/analysis/check time | 10.34 seconds, excluding writing |

Early growth is predominantly along-track. At day 60 nominal-frame radial and
along-track sigmas both approach 18,800 km because Cartesian covariance covers
the whole ring. Phase-matched radial residual sigma is only about **234 km**
and cross-track sigma **72.5 km**: global radial sigma is not tube thickness.
The screenshots show the initial ball, an arc at 2.5 days, and a full ribbon
at 60 days. Results: `outputs/meo-ball/`; screenshots: `build/viewer-qa-actual/`.
Generated artifacts are gitignored.

## Sample-count, seed and cadence convergence

The completed study used seeds 11,22,33, N=2,000/5,000/20,000 and cadences
0.5/0.25 days: **18 scientific runs, 162,000 particles**, plus four thread
benchmarks containing 8,000 additional particles.

```powershell
python tools/study.py --exe build/distribution_propagator.exe `
  --config examples/meo_ball.cfg --output outputs/convergence `
  --counts 2000,5000,20000 --seeds 11,22,33 --threads 8 `
  --cadences 0.5,0.25 --benchmark-threads 1,2,4,8
```

Ranges below describe **first sampled passing epochs across seeds**, not
confidence intervals. Subtract one cadence for each run's preceding endpoint.
Persistence stays fixed at three snapshots; its physical duration therefore
changes with cadence.

| N | Cadence (days) | Coverage epochs (days) | Mixing epochs (days) | Central-95% wrap epochs (days) |
|---:|---:|---:|---:|---:|
| 2,000 | 0.50 | 8.0–8.5 | 15.5–20.5 | 11.0–11.5 |
| 2,000 | 0.25 | 7.75–8.5 | 15.5–20.25 | 10.75–11.25 |
| 5,000 | 0.50 | 7.5–8.0 | 16.5–17.5 | 11.0 |
| 5,000 | 0.25 | 7.25–8.0 | 16.25–17.25 | 11.0 |
| 20,000 | 0.50 | 6.0–6.5 | 16.5 | 11.0 |
| 20,000 | 0.25 | 6.0–6.25 | 16.5 | 11.0 |

All three 20,000-sample quarter-day runs give a mixing bracket of
**16.25–16.5 days**. This supports at least 20,000 particles for this example
and these thresholds, not a universally sufficient count. The 2,000-particle
results fluctuate substantially. Coverage gets earlier with more particles
because more low-density tails are represented: all-bin occupancy is not a
sample-count-independent physical wrap time.

The complete JSON/CSV, exact child commands and executable hashes are under
`outputs/convergence/`. Completed study JSON SHA-256:
`b1ceea472fe47a0a2e1385437a6c799fcd1e49f95f0f89a1300ec638af7ee190`.
No confidence interval is manufactured from three seeds. See
[STUDIES.md](STUDIES.md) for the workflow and interpretation.

## Throughput

Identical N=2,000, seed=11, 121-epoch inputs:

| Workers | Simulation seconds | Speedup from one worker |
|---:|---:|---:|
| 1 | 18.477 | 1.00 |
| 2 | 11.909 | 1.55 |
| 4 | 6.250 | 2.96 |
| 8 | 4.006 | 4.61 |

Every thread count produced byte-identical initial samples and metrics. These
are single runs on a shared Windows machine, not controlled hardware benchmarks.
The 20,000-particle science runs took about 41–53 seconds at half-day cadence
and 81–86 seconds at quarter-day cadence.

Separately, retaining native Hansen polynomial tables reduced the one-year,
100-snapshot osculating workload from 4.154 to 0.382 seconds (about **10.9x**),
with unchanged derivative counts and byte-identical 32-case reference output.
[PERFORMANCE.md](PERFORMANCE.md) contains the focused benchmark and reproducer.

## Ten-year sensitivity experiment

```powershell
.\build\distribution_propagator.exe --config examples/meo_correlated.cfg `
  --samples 1000 --threads 4 --accuracy-check 16 --output outputs/meo-ten-year `
  --no-html --visual-samples 0
```

This covariance uses 1 km position sigmas, 0.1 m/s velocity sigmas and +0.5
correlation between radial position and along-track velocity. N=1,000 over
3,650 days at ten-day cadence completed in **57.99 s**. Coverage: **620–630
days**; central-95% wrapping: **760–770 days**; mixing: **1,030–1,040 days**.
Tighter integration of 16 samples changed positions by at most **0.190 mm**
and mean longitude by 7.28e-12 radians.

This small-sample long-horizon experiment has not received the convergence study
above. The independent Java comparison reaches one year, not ten years.
Numerical agreement does not remove physical errors from omitted perturbations.

## Windows runtime and compiler checks

The `nanosleep64` entry-point error came from incompatible GCC/pthread runtimes:
current UCRT libstdc++ imports the symbol, while the older Meld libwinpthread on
PATH does not export it. CMake now deploys the **matching compiler's** three
runtime DLLs beside each GCC executable. Launch tests passed with Meld first on
PATH and UCRT absent. Keep those DLLs when moving the executable.

MSVC was tested with Ninja in the Visual Studio developer environment. Its
MSBuild generator encountered a sandbox FileTracker access error during compiler
discovery; Ninja avoided that requirement. No privilege escalation was needed.
