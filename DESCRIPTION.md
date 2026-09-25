# Implementation and numerical design

This document describes the implemented simulator, its data contracts, and its
extension points. [README.md](README.md) covers compilation and use.
[docs/SCIENCE.md](docs/SCIENCE.md) critiques the ball-to-banana-to-ribbon
hypothesis, derives phase-shear estimates, and discusses sample sizes and model
limitations. The simulator transports a specified initial distribution; it does
not perform orbit determination from observations.

## Modules and execution flow

| Module | Responsibility |
| --- | --- |
| `orbit.hpp`, `orbit.cpp` | State conversion, orbital validity, vector/RTN operations, covariance factorization. |
| `statistics.hpp`, `statistics.cpp` | Sample statistics, circular phase diagnostics, persistence intervals, DKW bounds. |
| `config.hpp`, `config.cpp` | Configuration parsing, ensemble sampling/import, output epochs, storage estimates. |
| `backend.hpp`, `backend.cpp` | Native DSST force configuration, adaptive mean-element integration, mean/osculating conversion. |
| `simulation.hpp`, `simulation.cpp` | Particle scheduling, retained states, nominal reference, diagnostics, accuracy checks. |
| `output.cpp`, `web/viewer.html` | JSON/CSV results and a self-contained HTML visualization. |
| `main.cpp` | Command-line options, invocation, and concise result reporting. |

CMake exposes `distribution_core`, `distribution_backend`, and
`distribution_simulation` libraries, plus the `distribution_propagator`
executable. The native DSST checkout is an external dependency linked through
`dsst_cpp`; normal execution requires no Java runtime. Java is used only to
generate independent reference fixtures.

The run reads and validates configuration, samples the initial distribution,
constructs output epochs, propagates the nominal reference, propagates every
particle, computes diagnostics over the full ensemble, optionally performs a
tighter integration comparison, and writes results. An invalid sample or a
failed propagation aborts the run with its particle/row context. The simulator
does not silently drop, clip, or redraw invalid members.

## Coordinates, units, and principal types

All particles share one epoch and an inertial geocentric equatorial frame. The
backend epoch is elapsed time zero, expressed in seconds. The selected force
models use a fixed Earth symmetry axis and do not require a calendar epoch.
Angles are radians internally, distances metres, velocities metres per second,
and `mu` cubic metres per square second. Configuration fields explicitly ending
in degrees or days are converted at the boundary.

`Elements` is `std::array<double,6>` with this exact order:

```text
[a, ex, ey, hx, hy, lambda]
ex = e*cos(omega+Omega)        ey = e*sin(omega+Omega)
hx = tan(i/2)*cos(Omega)       hy = tan(i/2)*sin(Omega)
lambda = M+omega+Omega        # continuous mean longitude
```

`lambda` is allowed to accumulate many revolutions. It is not true longitude
or true anomaly. The five other elements describe a mean orbit during
integration, or an osculating orbit when short-period terms have been applied.
This convention is regular at circular and equatorial prograde orbits; a
nearly 180-degree inclination reaches the unsupported equinoctial singularity.

`Cartesian` is `[x,y,z,vx,vy,vz]`; `Vector3` is three doubles; `Matrix6` is a
6-by-6 array. A covariance has the squared/cross-product units of its declared
coordinates. A Cartesian velocity covariance therefore cannot be substituted
for an equinoctial element covariance without a coordinate transformation.

| Type | Contents and interpretation |
| --- | --- |
| `BackendConfig` | Force selection/constants, initial/output state types, error tolerances, and step limits. |
| `BackendState` | The integrated `mean` elements and the selected output elements in `osculating`. With `output_type=mean`, **both fields contain mean elements**; retain the configured label. |
| `BackendStats` | Counts of accepted/rejected integration steps and native derivative evaluations. |
| `Config` | Backend settings, nominal orbit, covariance or empirical path, sample/worker/display counts, time grid, diagnostic thresholds, memory budget, and export options. |
| `SampleState` | Six output elements plus continuous mean longitude and mean semimajor axis. Keeping both mean quantities makes phase diagnostics independent of short-period output excursions. |
| `PhaseMetrics` | Unwrapped phase standard deviation and central 95% width, four circular harmonics, largest gap, occupancy, histogram total variation, and bin counts. |
| `Frame` | One output time, nominal reference state, complete-ensemble diagnostics, display particle positions, a sampled reference ellipse, and its six selected-output reference elements. |
| `EventInterval` | Whether a sustained event was found, its preceding/current output times, and the index at which persistence was confirmed. |
| `Simulation` | Effective configuration, initial members, all retained states/frames, event intervals, analytic estimate, timing, and validation statistics. |

The stable particle identifier is its index in `Simulation::initial`. Retained
states are indexed as `states[frame_index * sample_count + particle_id]`.

## Initial distribution and covariance sampling

Configuration accepts one nominal orbit representation: Keplerian elements
with explicitly degree-valued angles, equinoctial elements, or an absolute
Cartesian state. Uncertainty coordinates are separately declared as
`equinoctial`, `cartesian`, or `rtn`.

For Gaussian sampling, the user supplies six standard deviations or a complete
covariance. `covariance_factor` first extracts standard deviations `s_i` and
normalizes the covariance to correlation units:

```text
C_ij = P_ij/(s_i*s_j)
C = Lc*transpose(Lc)
L_ij = s_i*Lc_ij
delta = L*z, where the six components of z are independent N(0,1).
```

The lower-triangular factorization accepts positive semidefinite matrices,
including rank-deficient cases, up to scaled roundoff tolerances. Symmetry and
PSD checks operate in appropriate units instead of applying one absolute
threshold across metres and dimensionless elements. An exactly zero variance
requires an exactly zero covariance row. Deterministic axes remain
deterministic; no artificial diagonal jitter is added. Nonfinite, asymmetric,
or materially indefinite covariances are rejected.

The generator uses `std::mt19937_64`, an explicitly specified 53-bit uniform
mapping, and Box-Muller transformation with a cached second normal deviate.
Samples are generated serially before workers start. A seed therefore determines
the same particle sequence regardless of worker scheduling or requested thread
count. Increasing the Gaussian sample count preserves the previous prefix.
Cross-platform bitwise identity of transcendental functions is not promised.

Offsets are applied in the declared coordinates:

- **Equinoctial:** add all six offsets to the nominal elements, preserving
  continuous longitude.
- **Cartesian:** convert the nominal orbit to Cartesian coordinates, add the
  six offsets, and convert the resulting absolute state to elements.
- **RTN:** rotate position and velocity error vectors through the **same fixed
  nominal-epoch RTN basis**, add them to the nominal Cartesian state, then
  convert to elements.

For RTN, `R = r/|r|`, `N = (r cross v)/|r cross v|`, and `T = N cross R`.
The velocity entries are components of an inertial velocity error expressed in
that basis. They are not time derivatives of a displacement in a rotating frame;
the implementation does not add an `angular_velocity cross displacement` term.

Converting a Cartesian state yields a longitude reduced modulo one revolution.
Both Gaussian Cartesian/RTN samples and **empirical Cartesian samples** are
lifted to the branch nearest the nominal initial longitude:

```text
lambda = lambda_nominal + remainder(lambda_converted-lambda_nominal, 2*pi).
```

This prevents a narrow cloud straddling zero from appearing initially spread
over an entire revolution. A Cartesian snapshot cannot reveal an unknown number
of completed revolutions; this nearest-branch choice is the explicit convention.

Empirical input contains six numeric columns per non-comment row, without a
column header. Rows are equally weighted and all are used. Their meaning is:

| Coordinates | Empirical row meaning |
| --- | --- |
| `equinoctial` | Absolute elements; the user's explicitly supplied unwrapped longitude is preserved, including its revolution count. |
| `cartesian` | Absolute inertial position and velocity; the longitude branch is lifted near the nominal orbit after conversion. |
| `rtn` | Offsets from the nominal state using the fixed basis described above. |

At least two empirical members are required. This input can represent a
non-Gaussian or multimodal posterior, but arbitrary supplied rows do not
automatically satisfy independent-random-sample statistical assumptions.

## State conversion and validity

`to_cartesian` solves the equinoctial Kepler equation
`F-ex*sin(F)+ey*cos(F)=lambda` on a reduced angle. Its bracket has half-width
`e`; safeguarded Newton steps maintain the monotonic root bracket. Equinoctial
orbital-plane basis vectors then map position and velocity into inertial axes.
`from_cartesian` derives energy, angular momentum, the eccentricity vector, the
equinoctial plane basis, and eccentric/mean longitude.

Validation requires finite values, positive `mu` and semimajor axis, an elliptic
eccentricity, and a supported orbital-plane parameterization. The configured
perigee floor is `earth_radius_m + minimum_altitude_m`. It is checked on initial
members and on each saved particle output. This is an orbit-domain check at
saved epochs, not a continuous drag assessment or an atmospheric model.

## DSST adapter and adaptive integration

Each `Backend` owns its mutable native `dsst::DSSTPropagator` and native force
objects. The implemented selections are:

| `force_model` | Native force configuration |
| --- | --- |
| `kepler` | Newtonian central attraction. |
| `j2` | Central attraction and `DSSTZonal` with degree 2, order 0, unnormalized `C20=-J2`. |
| `j2_j2sq` | The above plus `DSSTJ2SquaredClosedForm` using the native Zeis model. |

No drag, third-body gravity, radiation pressure, higher-degree zonals, tesseral
gravity, maneuvers, measurements, or process noise are silently included.
Providing negligible-drag initial orbits does not make these omitted effects
universally negligible over arbitrary horizons.

The application uses native DSST mean-element derivatives and short-period
formulas. Its integration driver is an **adaptive Dormand-Prince embedded 5(4)**
method in `backend.cpp`. The upstream native shell's fixed-step RK4 remains an
independent implementation used in comparison tests; it is not the application's
production integration driver.

For the degree-2 zonal profile, the adapter retains native Hansen polynomial
tables per particle instead of rebuilding them at every derivative evaluation.
It calls the native `DSSTZonal::initializeStep`, `createUAnddU`, and
`computeMeanElementRates` methods; `createUAnddU` still recomputes all
orbit-dependent Hansen roots. Native J2-squared rates and Newtonian rates are
then added in the shell's original order. This cache changes allocation work,
not the DSST formulas or error controls. The 32-case before/after propagation
output was byte-identical. See [docs/PERFORMANCE.md](docs/PERFORMANCE.md) for
measured timings, reproduction instructions, and the cache/thread-safety scope.

At construction, the adapter initializes native short-period terms with
`beforeIntegration`. For an osculating initial state it then calls
`computeMeanState` with a `1e-14` conversion threshold and 200-iteration limit,
and initializes again on the resulting mean orbit. This explicit initialization
is necessary for the native shell's conversion lifecycle. A mean initial state
skips the osculating-to-mean conversion.

`advance(elapsed_seconds)` accepts finite, nondecreasing times and continues from
the previous **mean** state. It never converts a prior osculating output back to
mean at every snapshot. Adaptive integration uses these rules:

1. Assemble seven Dormand-Prince stages in fixed-size arrays using native DSST
   derivatives. Reuse the previous accepted final derivative as the next first
   derivative (FSAL); rejected steps also retain the valid current-state rate.
2. Compare embedded fourth- and fifth-order increments componentwise. Normalize
   by `absolute_tolerance + relative_tolerance*magnitude`, then take the maximum
   component error.
3. Use `magnitude=max(abs(current),abs(candidate))` for the first five elements.
   For longitude, use **magnitude 1**, so the permitted phase error does not grow
   merely because many orbital revolutions have accumulated.
4. Accept when the normalized error is at most one. Select the next step with
   the usual `0.9*error^(-1/5)` scale, bounded growth/shrink factors, and configured
   step limits. Clip to each exact output endpoint. Throw if the requested error
   cannot be achieved at the minimum step.
5. Accumulate accepted longitude increments with compensated summation to reduce
   long-horizon roundoff. Count accepted/rejected steps and derivative calls.

The semimajor-axis absolute tolerance is in metres; other element absolute
tolerances use their native units. Tolerances control local integration error
inside the selected force model, not the physical error of an orbit prediction.

With osculating output selected, the adapter reconstructs short-period terms at
the exact requested output date after reinitializing the native coefficient
slots. It does not interpolate those terms across widely separated saved epochs.
With mean output selected, it avoids that reconstruction and labels the result
accordingly. The instantaneous geometric cloud therefore depends on output type.

## Parallel execution, storage, and accuracy checks

Workers atomically claim particle identifiers. A worker constructs a backend for
one particle, advances it sequentially through all output times, stores its
results in disjoint frame-major array entries, then claims another particle.
No propagator instance is advanced concurrently by multiple threads. A mutex
protects construction of the upstream zonal coefficient cache; subsequent
propagation uses each particle's own mutable force objects.

The nominal orbit is propagated once. After all workers join, diagnostics are
computed serially in stable particle order, so their floating-point reduction
order is independent of the worker count. The first worker failure is retained
and rethrown with particle context after the workers finish. The reported native
step/evaluation totals are for ensemble propagation; elapsed run time also
includes setup, nominal propagation, diagnostics, and requested accuracy checks,
but excludes writing output files.

All particles at all output epochs are retained: the principal storage is
`O(sample_count * frame_count)`, with eight doubles per `SampleState`. Analysis
requires additional per-frame scratch vectors and a sorted phase copy. Histograms
use `O(frame_count * phase_bins)` storage. A conservative preflight estimate also
includes display positions, 181-point reference ellipses, output-string copies,
and worker overhead. Gaussian runs check the budget before allocating samples;
empirical runs check it as rows are admitted. A final check precedes retained
state allocation. This is an estimated memory budget, not an operating-system
resident-memory limit.

`visual_samples` controls only displayed positions: the first `min(visual_samples,N)`
stable particle identifiers are retained for viewing. Every diagnostic uses all
particles. For random Gaussian samples the displayed prefix is itself an IID
subset; for ordered empirical input, the prefix may not visually represent the
entire supplied distribution. Reorder empirical rows appropriately or display
the full ensemble when that distinction matters.

An optional accuracy check repropagates the first configured number of particles
with all error tolerances divided by ten and maximum step divided by four,
subject to the same minimum step. It reports the maximum Cartesian position
difference and continuous mean-longitude difference at saved epochs. This is a
numerical convergence diagnostic on selected members, not a global error bound,
an independent force-model validation, or an automatic acceptance gate.

## Phase diagnostics and event times

At each output time, phase is the particle's continuous **mean longitude minus
the nominal mean longitude**. Welford accumulation gives sample standard
deviations with denominator `N-1`. The central 95% width is `Q(0.975)-Q(0.025)`
of the unwrapped values, using linear interpolation at index `p*(N-1)`.
It may legitimately exceed `2*pi`; it is not a shortest circular 95% arc.

For circular statistics, phase is reduced to `[0,2*pi)`. Four harmonics are
computed with extended-precision sums:

```text
R_k = abs(sum(exp(i*k*phase))/N), k=1,2,3,4.
```

Sorting the reduced phases gives every adjacent gap, including the last-to-first
gap across zero. A fixed equal-angle histogram supplies the occupied-bin fraction
and total variation from uniform phase:

```text
TV = 0.5*sum_bins(abs(count_bin/N - 1/bin_count)).
```

The histogram branch is relative to the nominal phase. Binned diagnostics can
depend on bin edges, sample count, and cadence. Circular harmonic magnitudes
and largest gap are invariant under a common angular rotation.

At each frame, the implemented Boolean criteria are:

```text
coverage = maximum_gap <= configured_gap
           AND occupied_fraction >= configured_fraction

mixed    = coverage
           AND max(R1,R2,R3,R4) <= configured_resultant
           AND TV <= configured_TV

central95_wrap = unwrapped_central_95_percent_width >= 2*pi.
```

`first_sustained` requires each criterion at `persistence` successive saved
epochs. If the qualifying run begins at frame `j`, the reported interval is
`[time[j-1],time[j]]`, confirmed at `j+persistence-1`. An initially satisfied
event has interval `[0,0]`. The reported event time is the interval's upper
endpoint, not the later confirmation epoch. `null` in JSON means no sustained
event was observed within the simulated horizon.

These are sampled onset intervals. They do not exclude an earlier short event
between snapshots, prove permanent mixing afterward, or establish positive
probability at every mathematical angle. No threshold automatically classifies
the cloud as a geometrically thin ribbon.

Two complementary position summaries are retained:

- `rtn_sigma_m` projects particle-minus-nominal Cartesian position into the
  nominal instantaneous RTN basis. After wrapping, this can be large in more
  than one axis even when early growth was predominantly along-track.
- `tube_rtn_sigma_m` evaluates the nominal output ellipse at **each particle's
  output mean longitude**, subtracts that phase-matched reference position,
  and projects into its local RTN basis. It removes much of the phase separation
  before measuring element/plane differences. It is not a closest-point
  distance to the curve, a confidence contour, or a fitted tube surface.

Both are standard deviations about the residual mean, not RMS distances from
zero. `semimajor_sigma_m` always uses the mean semimajor axes.

## Analytic estimate and statistical interpretation

The analytic mixing estimate uses initial **mean** elements after any required
DSST input conversion. For each particle it computes
`n=sqrt(mu/a)/a`; stable one-pass sums estimate initial phase variance `Vp`,
mean-motion variance `Vn`, and their covariance `Cpn`. It solves

```text
Vp + 2*Cpn*t + Vn*t^2 = -2*log(mixing_max_resultant)
```

for the nonnegative future crossing, using a cancellation-resistant quadratic
root. It returns zero when the initial variance already exceeds that Gaussian
threshold, and an unavailable estimate when there is no mean-motion variance
and the threshold is not initially met. The stored negative unavailable sentinel
is serialized as JSON `null`.

This estimates only Gaussian Keplerian phase shear. It does not reproduce the
complete coverage/TV/four-harmonic criterion or establish validity for arbitrary
empirical multimodal input. `dkw_cdf_error_95` reports
`sqrt(log(40)/(2*N))`, a fixed-time one-dimensional IID CDF bound. Its applicability
message distinguishes Gaussian random draws from arbitrary empirical ensembles.
See [docs/SCIENCE.md](docs/SCIENCE.md) for derivations and convergence guidance.

## Output contract and visualization

For empirical ensembles, metadata `covariance` is `null`: the unused Gaussian
configuration matrix does not describe the supplied particles. Accuracy-error
fields are `null` when `accuracy_checked_samples` is zero, rather than reporting
an unmeasured zero error.

| Artifact | Contents |
| --- | --- |
| `run.json` | Schema version 1; metadata/configuration and native DSST revision; summary; per-frame display positions, reference ellipse/elements, and full-ensemble diagnostics. |
| `summary.json` | Sample/epoch/worker counts, elapsed time, event times and intervals, analytic estimate, DKW information, accuracy-check results, and native step counters. |
| `metrics.csv` | One row per output epoch, with time, unwrapped/circular phase statistics, RTN and phase-matched residual spreads, semimajor spread, and current-frame flags. |
| `initial_samples.csv` | Every absolute initial equinoctial member with a header and particle identifier. To reuse it as six-column empirical input, remove the header/identifier column. |
| `states.csv` | Optional full-ensemble Cartesian outputs, mean longitudes, and mean semimajor axes at every saved epoch. Interpret Cartesian coordinates using metadata's output state type. |
| `visualization.html` | Optional offline viewer with the same run JSON embedded; no external assets or network calls are required. |

Floating-point output uses 17 significant digits. JSON preserves the seed both
as a number and as `seed_string`, because JavaScript numbers cannot exactly
represent every unsigned 64-bit seed. String escaping includes `<` to keep
embedded data from terminating an HTML script element. Existing nonempty output
directories are rejected instead of overwritten.

CMake embeds `web/viewer.html` in a generated header; the C++ writer inserts the
run data into its single placeholder. The viewer renders saved snapshots in a
canvas, with playback, a time slider, camera rotation/pan/zoom, view presets,
whole-orbit/focused-cloud framing, the phase histogram, and metric history.
Playback steps through saved epochs; it does not invent interpolated particle
trajectories. Colors retain particle identities and do not encode probability
density. The reference curve is an instantaneous nominal ellipse, not a particle
trail or a force-integrated full revolution.

### Local cross-section geometry

Each frame also exports `reference_elements = [a,ex,ey,hx,hy,lambda]` for the
same instantaneous nominal orbit used by `reference_orbit_m`. Its mean or
osculating interpretation follows `metadata.output_type`. This additive schema
version 1 field lets the viewer follow plane/perigee precession and evaluate
positions and tangents analytically instead of differentiating the sampled
polyline. If these elements are missing or invalid, the new panel is disabled
with a regeneration message; the original scene, diagnostics, and playback
remain available. Initial `metadata.nominal_elements` are not substituted for
later reference geometry.

Let `P,Q` be orthonormal periapsis/transverse vectors in the nominal orbital
plane, `N=P cross Q`, `e` its eccentricity, and `p=a*(1-e^2)`. At selected true
anomaly `f`, the reference position, velocity direction, and section axes are:

```text
r_ref(f) = p/(1+e*cos(f)) * (cos(f)*P + sin(f)*Q)
v_ref(f) = sqrt(mu/p) * (-sin(f)*P + (e+cos(f))*Q)
T(f) = unit(v_ref(f))
X(f) = unit(T(f) cross N)       # horizontal, in-plane normal to travel
Y    = N                      # vertical, normal to the orbital plane
```

The two plot axes are perpendicular to travel. `X` differs from radial away
from the apsides of an eccentric orbit. For `e < 1e-12`, the angle origin uses
the equinoctial x-axis instead of an undefined periapsis. The 0°–360° slider
sets `f`; follow mode solves the nominal Kepler equation to locate the nominal
position at each saved epoch. Changing the angle manually disables following.

For each displayed Cartesian position `r_i`, projection into the current
reference plane gives its angular location `f_i`. A particle is selected when

```text
abs(atan2(sin(f_i-f), cos(f_i-f))) <= full_slice_width/2
```

with a small roundoff allowance at the boundary. The full width ranges from
1° to 40° and defaults to 10°. This wrapped angular selection excludes the
opposite orbital branch and is continuous across 0°/360°. Positions with an
undefined in-plane angle are not selected.

The default curvature correction subtracts the reference position at each
particle's projected angle and transports its local transverse components to
the chosen section:

```text
delta_i = r_i - r_ref(f_i)
x_i = dot(delta_i, X(f_i))
y_i = dot(delta_i, N)
```

Thus particles lying exactly on the reference ellipse have zero corrected
width. With correction off, the raw projection instead uses
`delta_i = r_i-r_ref(f)`, `x_i=dot(delta_i,X(f))`, and `y_i=dot(delta_i,N)`.
A circular reference segment with half-width `h` alone then produces an
inward offset `a*(1-cos(h))`: about 101 km for a 10° full slice at
`a=26,560 km`. The mode labels distinguish this geometric broadening from
distribution thickness. Neither mode propagates, interpolates, or changes
particle positions; the same saved epoch drives both views.

Cross-section points use only the exported display subset and preserve its
particle IDs/colors. Counts always state selected versus displayed samples;
fewer than 20 selected points trigger a sparse-slice message. Equal scales on
both axes preserve shape, while the automatic extent changes between slices.
Larger widths include more particles but combine more orbital locations. The
view is not a density estimate, confidence contour, or full-ensemble statistic.
Its true-angle matching and axes also differ from `tube_rtn_sigma_m`, which
uses output mean longitude and RTN axes over the complete ensemble. Existing
metrics, event detection, and saved-epoch playback semantics are unchanged.

## Coverage parameter studies

`MeanPhaseLaw` contains the initialized mean elements, epoch, and native constant
mean-longitude derivative. Each coverage worker owns a `MeanPhaseFactory` that
reuses native force objects and invariant Hansen tables, repeats the unchanged
native mean conversion for every sample, and skips unused short-period output
preparation. `coverage.cpp` stores relative initial phase and rate
in particle-major arrays, then evaluates the exact affine phase flow of the
currently supported fixed gravity models. `CoverageResult` records sampled onset
and confirmation, cadence, evaluated epochs, native rate spread, runtime, and
initial-perigee tail counts. Invalid whole ensembles have an explicit status and
no invented event time. The ordinary simulator's stricter sampler is unchanged;
the coverage audit uses a separate entry point that retains bound low-perigee
draws long enough to count and classify them.

The scanner keeps per-bin counts, minima, and maxima. For the default occupied
72-bin/5-degree condition, only inter-bin gaps can exceed the threshold, giving
an exact linear-time Boolean. It falls back to sorted phases when that proof
does not apply and always sorts once for the reported onset gap. It scans epochs
in order because the coverage criterion need not be monotonic. No frame-major
Cartesian archive or short-period output reconstruction is required.

`coverage_main.cpp` validates a six-column cases CSV before processing it and
flushes one JSONL record per ensemble. The standard-library Python grid driver
distributes native batches across processes, fingerprints the executable and
configuration, resumes matching per-case checkpoints, and exports ordered
JSON/CSV. The HTML renderer embeds validated data and a local Plotly bundle.
The browser interpolates log coverage time across orbital nodes while retaining
the original uncertainty axes. Metadata can explicitly select logarithmic
geocentric perigee radius as the altitude coordinate; inclination and eccentricity
remain linear. Legacy reports without that policy retain linear-altitude weights.
This coordinate follows the leading power-law dependence of phase shear on orbit
size and is independently checked against withheld C++ simulations.
Missing contributing nodes
produce holes; intermediate slider positions never acquire a simulated status.
An optional separate reference surface retains the original MEO/PDF measurements
and cadence. Selecting it shows those exact values; changing a slider returns to
the main grid. It never becomes an extra interpolation corner.

`package_coverage_extension.py` checks complete source grids, executable/config
fingerprints, shared covariance axes and scientific conventions before merging.
It assigns unique display IDs while retaining source IDs and all measured fields.
It requires a complete passing withheld interpolation audit before producing the
HTML, CSV, provenance summary and checksums. The original MEO reference includes
a census of every saved initial particle for truthful perigee-tail readouts.
Full formulas and data conventions are in
[docs/COVERAGE_EXPLORER.md](docs/COVERAGE_EXPLORER.md).

## Tests and extension points

Tests use throwing checks that remain active in Release builds. The suites cover:

- `test_math.cpp`: analytic circular and independent classical perifocal
  geometry, energy/angular-momentum/plane invariants, mixed-unit and singular
  covariance factors, angular branch cuts, opposing-lobe harmonics, exact
  uniform bins, multi-revolution widths, equal-energy phase behavior, persistence
  intervals, and DKW sample-count thresholds.
- `test_config.cpp`: statistically specified Gaussian/Cartesian moments,
  seed/prefix reproducibility, deterministic axes, independent RTN mapping,
  Cartesian empirical branch lifting, preserved empirical revolution counts,
  configuration rejection, endpoint schedules, and memory-budget preflight.
- `test_backend.cpp`: long-horizon Kepler invariants, mean/osculating conversion,
  native RK4 comparisons, tighter tolerance/output-cadence comparisons,
  independent concurrent backends, and invalid backend settings.
- `test_cli.py`: executable-level options, artifacts/schema, exact Kepler
  ensemble behavior, thread-count reproducibility, empirical/covariance input,
  no artificial mixing at zero covariance, and rejected invalid runs.
- `test_study.py`: actual study invocations, Cartesian products of run settings,
  byte-identical thread benchmarks, null/zero event conventions and failed-child
  checkpoint retention.
- Coverage grid, validation, renderer and packaging tests cover safe measurement
  reuse, explicit interpolation coordinates, power-law interpolation oracles,
  complete-grid/provenance checks, and exact reference-surface values. Browser
  tests exercise sliders, masks, CSV export, reference selection and plot scales.
- `backend_probe.cpp` and the versioned Orekit fixture: 32 comparisons spanning
  J2 and J2-plus-J2-squared, mean/osculating input/output, and 0/1/30/365-day
  epochs. Optional upstream fixture groups provide additional native-formula
  regression checks; skipped or unavailable fixtures are not passing evidence.

The README records validation results and their force-model/horizon boundaries.
Tests comparing the same equations establish numerical consistency; they do not
certify a real object's physical uncertainty over a year.

For new forces, extend `BackendConfig`, force construction, metadata, and an
independent oracle together. Time-dependent ephemerides, body rotation, and
parameter uncertainty require explicit new data contracts. Preserve independent
mutable force objects across particles. For new distributions, keep their
coordinate/weight semantics explicit; weighted or correlated samples require
corresponding changes to diagnostics and statistical bounds. For new event
criteria, add direct counterexamples such as opposing lobes, branch-cut clouds,
and transient threshold crossings. Streaming/chunked output can replace retained
states for very large studies, but must preserve whole-ensemble diagnostics,
stable identifiers, and the ability to validate numerical accuracy.
