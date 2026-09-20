# Coverage time across orbit and uncertainty parameters

`distribution_coverage` computes orbital phase coverage for Gaussian Monte Carlo
ensembles using the native C++ DSST force model. The associated grid runner and
HTML explorer compare position and velocity uncertainty across selected perigee
altitudes, inclinations, and eccentricities.

The intended study uses **20,000 particles per measured grid point**. Each
particle receives its own sampled state and its own native DSST initialization
and mean-phase rate. Interpolated slider positions are estimates between those
measured cases; they are not additional particle simulations.

## Build and reproduce

Configure the project and its DSST dependency as described in the main README.
Build the additional native target in the same build directory:

```powershell
cmake --build build --config Release --target distribution_coverage
python tools/coverage_grid.py --exe build/distribution_coverage.exe --config examples/coverage_grid.cfg --output outputs/coverage_orbit_explorer
```

The grid runner accepts `--jobs` for concurrent cases and `--threads` for the
native worker count, as well as `--altitudes`, `--inclinations`,
`--eccentricities`, `--positions`, and `--velocities` for explicit grid choices.
Use its `--help` output for argument syntax and defaults. Python coordinates the
study; all particle initialization and phase propagation use the C++ executable.

For direct native batch execution:

```powershell
build/distribution_coverage.exe --config BASE.cfg --cases CASES.csv --output NEW.jsonl --threads 4
```

The cases file has this exact six-column header, without a trailing comma:

```csv
run_id,perigee_altitude_km,inclination_deg,eccentricity,position_sigma_km,velocity_sigma_m_s
example,2000,45,0.05,10,0.1
```

Run IDs must be unique and contain only ASCII letters, digits, hyphens, and
underscores. Quoted CSV fields are unnecessary and unsupported. Values must be
finite; altitude and eccentricity must be nonnegative, inclination must be in
0–90 degrees, eccentricity must be less than one, and both sigmas must be
positive. The command refuses to overwrite an existing output file and validates
all CSV rows before starting. It writes one flushed JSON record per completed
case. A process failure can leave a partial batch; the runner must not classify
that batch as a completed study.

## Orbit and uncertainty definitions

Altitude means **nominal osculating perigee altitude above the configured Earth
radius**, not semimajor-axis altitude or the instantaneous altitude at the
initial epoch:

```text
a = (Earth radius + perigee altitude) / (1 - eccentricity).
```

The CSV command fixes right ascension of the ascending node at 20 degrees,
argument of perigee at 30 degrees, and initial mean anomaly at 10 degrees. These
replace the corresponding nominal elements from the base configuration. The
study base configuration must retain `initial_type = osculating` to give the
altitude its stated convention; the lower-level scanner also supports mean
initial elements when explicitly configured.

Each uncertainty value is a Gaussian **1σ standard deviation on each of three
axes**, not a three-dimensional RMS radius. Sampling uses the nominal epoch RTN
frame with offset order `[dr, dt, dn, dvr, dvt, dvn]`:

```text
P = diag(sigma_position², sigma_position², sigma_position²,
         sigma_velocity², sigma_velocity², sigma_velocity²).
```

Position values from the CSV are converted from kilometres to metres; velocity
values are already in metres per second. The errors in velocity are inertial
velocity-vector offsets resolved in RTN, not derivatives in a rotating frame.
The blocks are isotropic, so their distribution is also isotropic in inertial
Cartesian coordinates. No initial correlations are assumed.

The command replaces the base covariance with this diagonal RTN covariance,
sets the geometric altitude floor to zero, disables particle display output,
and requests mean output for phase evaluation. It rejects an empirical sample
file rather than silently replacing it. The force model, particle count, seed,
and coverage criteria come from the study configuration. The normal study uses
the same seed at each node, pairing underlying normal draws across comparisons.
This is one ensemble realization per node; it does not estimate seed-to-seed
statistical confidence intervals.

## Exact phase flow within the configured DSST model

For the currently supported fixed central gravity, degree-2 zonal J2, and
optional Zeis J2-squared terms, the mean equations preserve semimajor axis,
eccentricity magnitude, and inclination. The equinoctial eccentricity and
inclination vectors rotate; the native mean-longitude derivative depends on
their invariant magnitudes and therefore remains constant.

After the normal native osculating-to-mean conversion, the phase-only scanner
evaluates that native derivative once per particle and uses

```text
lambda_j(t) = lambda_j(0) + lambda_dot_j * t
relative_phase_j(t) = [lambda_j(0) - lambda_reference(0)]
                    + [lambda_dot_j - lambda_dot_reference] * t.
```

This is the analytic solution of the **same configured DSST mean equations**,
not a Kepler-only substitute or a fitted coverage formula. Evaluating relative
phase directly also avoids subtracting two large absolute accumulated
longitudes at every epoch. No Cartesian positions or short-period reconstruction
are needed for this diagnostic. Storage scales with particle count rather than
particle count multiplied by the number of saved epochs.

Each worker owns a reusable `MeanPhaseFactory`, avoiding repeated construction
of native force objects and invariant Hansen tables. Every particle still runs
the native osculating-to-mean iteration at its original tolerance (`1e-14`, up
to 200 iterations). Only the final preparation of unused short-period output
is omitted. Tests require exact equality with a freshly constructed backend,
including when inputs are processed in a different order or by separate workers.

The optimization is restricted to these force models. Adding higher zonals,
third bodies, radiation pressure, drag, time-dependent force parameters, or
maneuvers requires re-establishing the invariant-flow property or returning to
the general integrator. Exact integration of the selected mean equations does
not remove model truncation error or omitted environmental physics.

## Coverage and observation times

The default criterion retains the earlier uncertainty-surface definition:

1. Every one of 72 equal relative mean-longitude bins is occupied.
2. The maximum circular gap is at most 5 degrees.
3. Both conditions hold at three consecutive evaluated epochs.

The plotted time is the first qualifying epoch, not the later confirmation
epoch. The preceding epoch and first qualifying epoch form a sampled-onset
bracket. That bracket represents temporal sampling resolution, not a Monte
Carlo confidence interval or a proof that an earlier short event did not occur
between saved epochs. Coverage need not remain true forever and does not mean
the phase density is uniform or that a geometric tube has a specified width.

The scanner avoids sorting every ensemble at every epoch. With all bins
occupied and bin width equal to 5 degrees, any gap wholly inside a bin is at
most 5 degrees. Only the gap between each bin's maximum phase and the next
bin's minimum can violate the limit, including the last-to-first wrap. Keeping
bin counts and minima/maxima therefore gives the exact coverage Boolean in
O(particles + bins) operations. Other bin/fraction/gap settings and ambiguous
floating-point boundary cases use the full sorted-gap calculation. The reported
gap at the onset is always recomputed from all sorted particle phases.

Automatic observation cadence is `1 / (120 * sigma_phase_rate)`, where the
standard deviation is measured from the sampled native phase rates. The initial
horizon spans 192 cadences and is extended if necessary, subject to configured
limits. This choice sets observation times; it does not predict or replace the
measured coverage event. `--fixed-times` instead uses `duration_days` and
`output_step_days`, including a distinct final partial interval. The scanner
stops once persistence is established. Its reported evaluated horizon is the
last evaluated epoch, not an unexecuted prospective horizon.

## Geometric validity and low-altitude tails

A Gaussian distribution has unbounded tails. A nominal perigee of 1000 km does
not imply that every perturbed orbit remains above 1000 km or above Earth. The
domain audit retains every original draw and counts perigees below Earth and
below 300, 500, and 1000 km; it does not clip or redraw particles.

If any initially sampled bound ellipse intersects Earth, the entire case is
marked `earth_intersection` and its coverage time is null. A mean ellipse that
intersects Earth after native initialization produces
`mean_earth_intersection`. These masked cases are not simulated coverage
values. The ordinary distribution simulator retains its original stricter
perigee rejection contract. Hyperbolic or otherwise unsupported samples remain
errors rather than being silently repaired.

For example, with seed 20260919, 20,000 particles, nominal circular perigee at
1000 km, position σ of 100 km, and velocity σ of 0.01 m/s, four sampled ellipses
intersect Earth; their minimum perigee is approximately −63.1 km. The same
ensemble contains 174 samples below 300 km and 1037 below 500 km.

**A zero-altitude floor is a geometric check, not evidence that drag is
negligible.** Valid cases with low-perigee tails remain idealized gravity-only
calculations. The per-case counts expose where the omission may matter;
quantifying atmospheric drag would also require physical object and atmosphere
assumptions that this study does not supply.

## Sliders, grid resolution, and validation

The initial study grid has perigee altitudes `[1000, 2000, 3000]` km,
inclinations `[0, 45, 90]` degrees, and eccentricities `[0, 0.05, 0.1]`. Each
orbit node contains an 11-by-11 logarithmic uncertainty grid spanning position
σ of 1–100 km and velocity σ of 0.01–1 m/s. This totals 3267 measured cases and
65.34 million particle initializations before additional validation cases.

Continuous slider values between measured orbit nodes interpolate **logarithms
of coverage time**, with weights linear in altitude, inclination, and
eccentricity. The result is positive and matches every valid measured node.
The plot must distinguish those interpolated values from measured cases.
Neither coverage-time brackets nor low-altitude counts should be presented as
new measured diagnostics at interpolated slider positions.

Interpolation must not bridge a missing or invalid support node. An unavailable
interpolated value is distinct from a measured Earth-intersecting ensemble:
some interior orbits might be valid even when interpolation cannot be supported
by all surrounding nodes. Refinement near those boundaries can narrow that
unavailable region.

**Midpoint validation is pending until the generated study records its
results.** Withheld orbit midpoints should be directly simulated, compared with
the surrounding log-time interpolation, and included in the study metadata.
The intended refinement threshold is approximately 5% relative time error,
with observation cadence and finite-ensemble threshold jumps reported rather
than hidden. Physical phase-shear trends are smooth, but the last-gap and
persistence criteria can have small irregularities even with common draws;
denser nodes alone do not eliminate that sampling variability.

The independent standard-library validator checks the complete Cartesian case
count, each 20,000-particle record and seed, event arithmetic, and interpolation
against those withheld cases:

```powershell
python tools/validate_coverage_grid.py --study outputs/coverage_orbit_explorer/study.json --holdouts outputs/coverage_orbit_explorer/holdouts/study.json --output outputs/coverage_orbit_explorer/validation.json
```

Use `--allow-partial` only for a provisional audit while the main grid is still
running. The report lists every comparison and its supporting nodes or missing
and invalid corners, plus maximum, 95th-percentile, and RMS relative errors.
Masked cases receive no prediction and do not enter the error statistics. A
five-percent error-gate failure or malformed data produces a nonzero exit code;
an empty set of supported comparisons is reported as not evaluated.

Regression tests compare the optimized coverage Boolean against sorted phase
metrics, phase flow against the general native integration path, and complete
fixed-epoch coverage events against the established simulator across all three
force-model choices. They also exercise circular/equatorial and polar cases,
angle wrapping, bin boundaries, zero phase shear, initial coverage, final
partial epochs, input validation, and Earth-intersecting Gaussian tails.
