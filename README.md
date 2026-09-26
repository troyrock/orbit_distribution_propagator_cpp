# Distribution propagator

C++ Monte Carlo propagation of initial Earth-orbit uncertainty using the native
Orekit DSST port. Follow nonlinear particle trajectories, measure orbital
coverage and phase mixing, and explore the cloud in a self-contained HTML viewer.

The equivalent native Python tool is
[orbit_distribution_propagator_python](https://github.com/troyrock/orbit_distribution_propagator_python).

The ball → banana → ribbon hypothesis is plausible when orbital-energy
uncertainty creates different orbital periods. Coverage and near-uniform mixing
are different events. Read [the scientific critique](docs/SCIENCE.md) and
[the algorithms and data structures](DESCRIPTION.md).

Measured example results, sample-count convergence and validation evidence are
in [docs/VALIDATION.md](docs/VALIDATION.md). For the broad demonstration, the
20,000-particle study found mixing at 16.25–16.5 days for all three tested seeds.

## Build

Requires CMake 3.20+, Git, a C++17 compiler and the external native DSST checkout.
Python 3 enables CLI tests and study scripts. Java is needed only to regenerate
reference fixtures; neither Java nor Python is required by the executable.

Clone this application, then obtain the public
[native DSST C++ dependency](https://github.com/troyrock/DSST-Cpp) at the validated
revision. These commands work in PowerShell and POSIX shells:

```text
git clone https://github.com/troyrock/orbit_distribution_propagator_cpp.git
cd orbit_distribution_propagator_cpp
git clone https://github.com/troyrock/DSST-Cpp.git external/DSST-Cpp
git -C external/DSST-Cpp checkout --detach e653cd40e5e057cb89395afd8a9e3012ad59b12c
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel 6
ctest --test-dir build -C Release --output-on-failure -j 1
```

CMake uses `external/DSST-Cpp` by default. To use an existing dependency checkout,
pass `-DDSST_SOURCE_DIR=/absolute/path/to/DSST-Cpp` when configuring. The build
does not download or modify the dependency. It records the actual dependency
Git revision in output metadata. The pinned revision is an Orekit 13.1.6 port.
If using a source archive without Git metadata, the recorded revision is `unknown`.

For Windows with MSYS2 UCRT GCC and Ninja, configure with the matching toolchain:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_CXX_COMPILER=C:/msys64/ucrt64/bin/g++.exe
cmake --build build --parallel 6
ctest --test-dir build --output-on-failure -j 1
```

For Windows GCC, the build copies its matching C++/GCC/pthread runtime DLLs
beside the executables, so running them does not depend on global `PATH` order.
Keep those three DLLs with the executable if you move it. This prevents the
`nanosleep64` entry-point error caused by older runtime DLLs from applications
such as Meld. For Visual Studio omit `-G Ninja` and the compiler argument, then
build with `--config Release` and test with `-C Release`. Use a separate build
directory when switching generators or compilers. Visual Studio places the
executable in `build/Release/`; Ninja and Makefiles place it in `build/`.

## Run and visualize

```powershell
.\build\distribution_propagator.exe --config examples/meo_ball.cfg --output outputs/meo-ball
Start-Process .\outputs\meo-ball\visualization.html
```

With Visual Studio use `.\build\Release\distribution_propagator.exe`. On Linux
or macOS use `./build/distribution_propagator`, then open the generated HTML in
your browser. No web server or network connection is needed for the viewer.

This deliberately broad demonstration uses isotropic 100 km position and 1 m/s
velocity standard deviations in RTN, around a 26,560 km semimajor-axis orbit
with eccentricity 0.02 and inclination 55°. It propagates 5,000 particles for
60 days with J2/J2-squared. This is **not a measured OD covariance**.

The offline viewer provides rotation, wheel zoom, shift-drag pan, orbital-plane
and edge-on views, cloud framing, exact-snapshot playback, a phase histogram and
time histories. Particle IDs remain stable. Displayed particles are a bounded
subset; metrics use all samples. Playback does not invent interpolated states.

The **Local cross-section** panel looks along the nominal direction of travel
to show the distribution's width in the perpendicular plane. The **Orbit angle**
slider places the section anywhere from 0° to 360° around the current reference
ellipse; an amber marker and highlighted arc locate it in the main orbit view.
Angles are true anomaly from perigee, or from the equinoctial x-axis for a
circular orbit. **Follow nominal position** is initially enabled; moving the
angle slider disables following. The **Full slice width** slider selects a
1°–40° angular neighborhood (10° initially).

**Remove reference-orbit curvature** is enabled by default. It aligns local
offsets after subtracting the reference ellipse at each selected particle's
angle, so a perfectly thin reference orbit has zero width. Turn it off to see
the raw projection onto one plane; a finite angular slice then includes
apparent width from the orbit's curvature. Both plot axes use the same distance
scale and automatically fit each slice. These are transverse sample offsets,
not a probability contour or a full-ensemble covariance.

The panel reports selected/displayed counts and flags sparse slices. For a
roughly uniform angular cloud, a 10° slice contains only about 42 of 1,500
displayed particles. Add `--visual-samples 5000` to a run with at least 5,000
particles to get about 139 in that slice, or increase the display count further
for denser shape inspection. Counts vary with angle and are not uniform in
true anomaly for eccentric orbits. Wider slices improve counts but average
over more orbital locations. Reports created before per-epoch reference
elements were added must be regenerated to enable this panel; their existing
3D view and playback remain supported.

| Example | Purpose |
|---|---|
| `examples/meo_ball.cfg` | Visible ball-to-ribbon demonstration |
| `examples/meo_energy_shear.cfg` | Controlled two-body phase-shear experiment |
| `examples/meo_correlated.cfg` | Smaller correlated covariance over ten years, conditional on J2/J2² |

```powershell
.\build\distribution_propagator.exe --write-example my-orbit.cfg
.\build\distribution_propagator.exe --config my-orbit.cfg --output outputs/my-orbit `
  --samples 20000 --threads 8 --accuracy-check 16
.\build\distribution_propagator.exe --help
```

Output directories must be absent or empty. Other overrides include `--seed`,
`--duration-days`, `--output-step-days`, `--no-html` and `--export-states`.
With no config the built-in broad equinoctial MEO demonstration is used.

## Inputs

Configuration uses one `key = value` per line, with `#` comments. Unknown and
duplicate keys are rejected. Paths are relative to the configuration file.
Units are SI, except explicit degree/day keys. All states share epoch zero
and a common inertial equatorial frame. No Earth-fixed transform or absolute
ephemeris date is implied. Choose exactly one nominal orbit representation:

```text
orbit_keplerian_deg = a_m, eccentricity, inclination_deg, RAAN_deg, argument_perigee_deg, mean_anomaly_deg
orbit_equinoctial = a_m, ex, ey, hx, hy, mean_longitude_rad
orbit_cartesian = x_m, y_m, z_m, vx_m_s, vy_m_s, vz_m_s
```

Here `ex=e*cos(ω+Ω)`, `ey=e*sin(ω+Ω)`, `hx=tan(i/2)*cos(Ω)`,
`hy=tan(i/2)*sin(Ω)`, and `lambda=M+ω+Ω`. Exact retrograde singularities,
nonfinite/unbound orbits and perigees below the configured floor are rejected.
The default 1,000 km altitude floor is a screen, not a universal no-drag guarantee.
Use `initial_type=osculating` for ordinary instantaneous OD states, or `mean`
for force-consistent DSST mean elements. `output_type=osculating` reconstructs
short-period geometry; `mean` is an explicitly labelled secular view.

Choose Gaussian `sigma` (six standard deviations) or `covariance_csv` (a full
six-by-six symmetric positive semidefinite matrix, headerless numeric rows).
Position/velocity correlation and zero-variance axes are retained without jitter.

| `uncertainty_coordinates` | Order |
|---|---|
| `equinoctial` | `[a,ex,ey,hx,hy,lambda]` in metres/dimensionless/radians |
| `cartesian` | Inertial `[x,y,z,vx,vy,vz]` in metres and metres/second |
| `rtn` | Nominal epoch `[dr,dt,dn,dvr,dvt,dvn]` in metres and metres/second |

RTN velocity entries are **inertial velocity errors resolved in the nominal RTN
basis**, not derivatives of rotating positions. Transform a covariance with a
different convention before use. Each sampled complete state is nonlinearly
converted to elements before any mean/osculating conversion. Invalid samples
stop with an index; no clipping or redrawing changes the specified posterior.

Alternatively set `empirical_samples_csv` instead of a covariance. Each
headerless six-value row has equal weight. Equinoctial/Cartesian rows are
absolute states; RTN rows are offsets. Every row is used; `--samples` is rejected.
Equinoctial longitudes preserve explicit winding numbers. Cartesian longitudes
are lifted near the nominal branch because Cartesian states lack winding history.

## Coverage, mixing and sampling

Metrics use continuous **mean longitude relative to the nominal orbit** and its
modulo-2π phase. For eccentric orbits, uniform mean phase is not uniform true
anomaly or spatial density. Default milestones are:

- Central 95% unwrapped width (97.5th minus 2.5th percentile) reaches 2π.
- Coverage: largest empty phase gap ≤5° and all 72 bins occupied.
- Mixing: coverage, all four resultants `R1..R4 ≤0.05`, and histogram total
  variation from uniform ≤0.15. Multiple harmonics detect symmetric lobes.

Crossings require three consecutive passing snapshots by default. Reports give
the first passing epoch and its preceding epoch as a **sampled onset bracket**,
not a continuous-time root or confidence interval. Refine cadence near the
crossing. Unobserved/unconfirmed events are `null`. These diagnostics depend
on count and thresholds, and a finite ensemble can later recur or fail a test.

An analytic Gaussian Kepler-shear comparison uses initial sampled mean elements,
including phase/rate covariance. It is a time-scale estimate, not the complete
numerical mixing condition. Mathematical Gaussian tails make “first nonzero
possibility anywhere” unsuitable as a finite-sample criterion.

Use **2,000–5,000 samples for exploration** and **20,000–50,000 for ordinary
studies**, then compare doubled counts and at least three independent seeds.
At 95% confidence the DKW one-epoch CDF bound is about 1.92% for 5,000 IID
samples and 0.96% for 20,000. It is not a wrap-time or multivariate shape bound.
Rare tails require more particles; a bin of probability `p` has relative
counting error approximately `sqrt((1-p)/(N*p))`.

Run repeatable sample/seed/cadence and thread studies with the standard-library
Python driver; [docs/STUDIES.md](docs/STUDIES.md) explains its outputs and checks:

```powershell
python tools/study.py --exe build/distribution_propagator.exe `
  --config examples/meo_ball.cfg --output outputs/my-study `
  --counts 2000,5000,20000 --seeds 11,22,33 --threads 8 `
  --cadences 0.5,0.25 --benchmark-threads 1,2,4,8
```

## Accuracy and performance

Models `kepler`, `j2`, and `j2_j2sq` use native DSST force formulas. An adaptive
Dormand–Prince 5(4) driver advances mean elements because the translated shell
only offers fixed-step RK4. Short-period terms are initialized before initial
mean conversion and reconstructed at exact output epochs. Default tolerances
are relative `1e-11`, absolute semimajor-axis `0.001 m`, and other-element
`1e-12`, with a one-day maximum step. Longitude error scaling remains bounded
as revolution count grows. No fast-math flags or simplified force equations are
used. Native Hansen polynomial objects are reused instead of rebuilt each rate
evaluation; every orbit-dependent coefficient is still recomputed.

`--accuracy-check N` reruns N particles with 10× tighter tolerances and a 4×
smaller maximum step, comparing position and phase at all output epochs.
This tests numerical convergence **within the selected model**. The independent
one-year Java Orekit fixture and its gates are documented in
[tests/data/README.md](tests/data/README.md).

No drag, higher gravity harmonics, tesseral resonances, Sun/Moon gravity, solar
radiation pressure, maneuvers, updates or process noise are included. Long-span
results are conditional model experiments, not certified real-object forecasts.
Add and validate forces needed by a real mission's error budget separately.

Each worker owns its mutable DSST objects; shared constructor caches are locked.
Samples are generated once from a specified `mt19937_64`/Box–Muller stream.
Worker count does not change results; larger counts preserve the same-seed
prefix. `max_memory_mb` preflights states, analysis, display and serialization
storage, but is not a hard OS cap. Reduce displayed particles or output cadence
before reducing the scientific count. HTML/JSON writing is serial.

## Outputs and tests

| File | Contents |
|---|---|
| `summary.json` | Events/brackets, count, runtime, integration work and numerical checks |
| `metrics.csv` | All-ensemble epoch statistics and RTN spreads |
| `initial_samples.csv` | Exact input equinoctial samples and stable IDs |
| `run.json` | Versioned settings/model provenance, metrics and display subset |
| `visualization.html` | Offline interactive report |
| `states.csv` | Optional complete Cartesian ensemble, position and velocity |

Global Cartesian/RTN covariance loses a local along-track interpretation after
wrapping. Metrics also include residual RTN spreads after matching the reference
ellipse at each particle's output mean longitude. This describes thickness
relative to that ellipse, not a nearest-point fit or proof of a thin tube.

Release-active tests cover orbit oracles, covariance moments, angle branch cuts,
equal-energy no-shear, non-Gaussian lobes, thread determinism, mean/osculating
round trips, tolerance refinement, Java parity and full CLI workflows.
This project's committed Java parity fixture runs without any Java installation
or external gravity data. The optional upstream Java fixture comparator is
also built when both its generated CSV and the adjacent Orekit Java source test
resource `potential/shm-format/eigen_cg03c_coef` are present. The build copies
those data into its test directory without changing their gates. A standalone
dependency checkout usually lacks the adjacent Java resource; CMake reports
the unavailable optional tests while retaining all application tests and this
project's independent Java parity fixture. Set
`-DDISTRIBUTION_OREKIT_FIXTURES=OFF` to disable the optional comparator explicitly.
`DISTRIBUTION_UPSTREAM_TESTS=ON` enables every upstream test; one current zonal
test passes a temporary to a mutable reference and GCC rejects it. Default
application/reference targets do not depend on that compiler-specific test.

For browser QA with Node, Playwright and Chrome/Chromium:

```powershell
node tools/test_viewer.cjs outputs/meo-ball/visualization.html --screenshot-dir build/viewer-qa
```

See [DESCRIPTION.md](DESCRIPTION.md) for extension points and
[docs/SCIENCE.md](docs/SCIENCE.md) for derivations and primary-source references.

## Coverage surface with orbital sliders

The additional `distribution_coverage` executable efficiently measures orbital
coverage for large parameter studies. The offline HTML report plots coverage
time against position and velocity uncertainty, with continuous sliders for
**perigee altitude**, inclination, and eccentricity. Values at simulated orbital
nodes are measured; intermediate slider positions are explicitly labeled as
interpolations of log time. Unsupported or Earth-intersecting cases leave gaps.

```powershell
cmake --build build --config Release
python tools/coverage_grid.py --exe build/distribution_coverage.exe `
  --config examples/coverage_grid.cfg --output outputs/coverage-orbits `
  --jobs 3 --threads 8
python tools/render_coverage_explorer.py outputs/coverage-orbits/study.json `
  --output outputs/coverage-orbits/explorer.html --plotly-js PATH/plotly.min.js
```

On multi-configuration generators, use `build/Release/distribution_coverage.exe`.
The driver and renderer use Python 3.10 or newer and its standard library; all
particle calculations run in C++. Supply a local Plotly JavaScript bundle (for example from an installed
Plotly package). It is embedded in the HTML, which works without network access.
The default study contains 27 orbital settings, each with an 11-by-11 uncertainty
surface and 20,000 particles per case. Repeating the driver command resumes
completed cases only when the executable, configuration, and axes match its
saved manifest. `--help` describes grid and concurrency options.
An operating-system lock prevents concurrent drivers from writing the same
study directory and is released automatically when the driver exits.

For a wider **1,000-30,000 km perigee-altitude** study, supply altitude nodes
explicitly. The optional reuse input imports compatible measured cases by their
parameters, avoiding duplicate particle calculations while retaining provenance:

```powershell
python tools/coverage_grid.py --exe build/distribution_coverage.exe `
  --config examples/coverage_grid.cfg --output outputs/coverage-orbits-30000 `
  --altitudes 1000,2000,3000,7000,12500,20000,24599.072067163703,30000 `
  --altitude-interpolation log_geocentric_perigee_radius `
  --reuse-study outputs/coverage-orbits --jobs 3 --threads 8
```

Reusing cases requires matching executable and configuration hashes. When
extending separately into disjoint altitude grids, `tools/package_coverage_extension.py`
can merge those grids after checking provenance and independent withheld
interpolation results. It can also embed the original measured MEO surface as
an **Original MEO / PDF** button: a = 26,560 km, perigee = 19,650.663 km,
i = 55 degrees, e = 0.02. That button retains the original measurements and
output epochs, providing a direct comparison with the PDF. Moving a slider
returns to the expanded grid. The optional common height/color scale supports
comparison across orbital settings; the fitted scale keeps each surface readable.

The fast phase flow is exact for this program's fixed degree-2 J2 plus J2-squared
mean equations after native initialization. It is restricted to those equations;
it does not add atmospheric or other omitted forces. See
[the coverage study guide](docs/COVERAGE_EXPLORER.md) for the coverage criterion,
domain masks, reproducibility, interpolation checks, and algorithm proof.
