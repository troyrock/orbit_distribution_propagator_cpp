# Orbital uncertainty: from a local cloud to a phase-mixed ribbon

## Assessment of the hypothesis

The proposed ball-to-banana-to-ribbon evolution is physically plausible for a
passive, bound Earth orbit with appreciable orbital-energy uncertainty. Different
energies imply different semimajor axes and orbital periods. Particles therefore
lose phase coherence even without atmospheric drag. Their position distribution
can bend around the Earth and eventually occupy every orbital phase.

This is a useful hypothesis to test, with several qualifications:

- Initial orbit determination usually supplies an anisotropic and correlated
  six-dimensional position/velocity distribution, not a literal ball. Position
  uncertainty alone is insufficient to predict phase spreading.
- In the two-body problem, identical semimajor axes imply identical mean
  motions. A distribution with exactly equal energy need not spread around the
  orbit at all. Uncertainty in eccentricity or plane orientation alone is not
  enough to guarantee eventual phase mixing.
- Initial phase/energy correlations can cause the cloud to contract before it
  expands. Multimodal input, resonances, and differential precession can produce
  several arcs, sheets, or a broad torus instead of a thin ribbon.
- A completed ribbon need not have constant spatial density. On an eccentric
  orbit, a phase-mixed object spends more time near apogee. Uniformity in mean
  anomaly is different from uniformity in true anomaly or distance along the
  orbit.
- In conservative dynamics the uncertainty is transported deterministically.
  A smooth six-dimensional Hamiltonian flow preserves phase-space volume while
  stretching and folding the distribution. Apparent spreading in three-dimensional
  position is not an assumption that the object receives random kicks.
- A Gaussian distribution has unbounded tails. A statement that an angle has
  mathematically nonzero probability is therefore a poor definition of the time
  a ribbon forms. A finite particle cloud also cannot establish zero probability
  between particles. Use explicit probability, angular-resolution, and persistence
  thresholds instead.

Non-Gaussian Cartesian orbital uncertainty from nonlinear dynamics is an
established phenomenon; a single propagated covariance ellipsoid can conceal
the curved shape. The NASA CARA study examines this issue with Monte Carlo
propagation, and Junkins, Akella, and Alfriend discuss the strong dependence on
coordinates. Neither reference establishes a universal time at which every
initial distribution becomes a ribbon. [NASA CARA study](https://ntrs.nasa.gov/api/citations/20120016932/downloads/20120016932.pdf?attachment=true),
[Junkins et al., author-institution publication record](https://scholars.library.tamu.edu/vivo/display/n162658SE).

## A two-body estimate before running a large ensemble

The formulas in this section are local derivations for interpreting the
experiment, not a substitute for its force model. Use metres, seconds, radians,
and an inertial geocentric frame.

Let `a` be semimajor axis, `mu` Earth's gravitational parameter, `n` mean motion,
and `T` the period:

```text
n = sqrt(mu/a^3),                 T = 2*pi/n
dn/da = -3*n/(2*a)
sigma_n ~= 3*n*sigma_a/(2*a).
```

For a sufficiently narrow ensemble, an unwrapped mean-phase deviation obeys

```text
delta_phi(t) ~= delta_phi(0) + t*delta_n
s_phi(t)^2 = s_phi(0)^2 + 2*t*Cov(phi(0),n) + t^2*sigma_n^2.
```

`phi` can be mean anomaly for nearly common orbital orientation, or mean
longitude for a nonsingular equinoctial description. In a perturbed model,
longitude also contains apsidal/nodal evolution, so the simple Keplerian rate
is only an estimate. Preserve continuous angles or count revolutions when
measuring the unwrapped spread; a standard deviation of angles reduced to
`[0,2*pi)` has a discontinuity at the branch cut.

For an initial inertial Cartesian covariance `P` in `[r_x,r_y,r_z,v_x,v_y,v_z]`,
the first-order semimajor-axis variance follows directly from orbital energy:

```text
energy = |v|^2/2 - mu/|r|
a = -mu/(2*energy)
g_r = 2*a^2*r/|r|^3
g_v = 2*a^2*v/mu
sigma_a^2 ~= [g_r,g_v] * P * transpose([g_r,g_v]).
```

All position/velocity cross-covariances matter. For example, near a circular
orbit, radial position and tangential velocity errors affect energy to first
order; radial velocity alone does not. For broad or non-Gaussian input, compute
each particle's `a` and `n` directly and use their sampled distribution instead
of linearizing this transformation.

If the unwrapped phase is Gaussian with variance `s_phi(t)^2`, reducing it modulo
`2*pi` gives a wrapped Gaussian. Its circular harmonic magnitudes are

```text
R_k(t) = |E[exp(i*k*phi(t))]| = exp(-k^2*s_phi(t)^2/2).
```

This formula is exact for the assumed Gaussian phase, but only approximate when
derived by linearizing mean motion about the nominal semimajor axis. It offers a
strong analytic check on a two-body Monte Carlo test. For a general distribution,
the expectation is its circular characteristic function and need not decay
monotonically.

Two useful, different reference times, neglecting initial phase spread and its
correlation with energy, are

```text
t_width = 2*pi/(6*sigma_n)       # unwrapped +/-3-sigma width equals one turn
t_R     = sqrt(-2*log(R_target))/sigma_n
```

The first is an intuitive tail-overlap scale. It is not a uniformity criterion:
at `t_width`, `R_1 = exp(-pi^2/18) ~= 0.578`. For `R_target=0.05`,
`t_R ~= 2.448/sigma_n`, about 2.34 times later. In orbital periods,

```text
t_R/T ~= 0.260*a/sigma_a       # R_target = 0.05
```

For illustration, a two-body orbit with `a=26,560 km` and
`mu=3.986004418e14 m^3/s^2` has an 11.966-hour period:

| Semimajor-axis standard deviation | +/-3-sigma width = one turn | Predicted R1 = 0.05 |
| --- | ---: | ---: |
| 1 km | 1,471 days | 3,439 days |
| 10 km | 147.1 days | 343.9 days |
| 100 km | 14.71 days | 34.39 days |

These are calculated scale estimates, not results for a measured orbit
determination solution. A tightly known orbit can require years or decades;
a deliberately broad demonstration ensemble can wrap in weeks. Long horizons
make force-model uncertainty and external perturbations increasingly important.

## What “all the way around” should mean

There is no unique physical transition time. Report the definition alongside
every number. The analysis should distinguish three questions:

1. **Angular coverage:** is there appreciable probability throughout a complete
   revolution at the chosen angular resolution?
2. **Phase mixing:** is the marginal distribution approximately uniform in a
   phase variable that advances uniformly in the two-body limit?
3. **Ribbon geometry:** is the radial/cross-plane thickness still small compared
   with orbital size, so that the ensemble actually resembles a thin ribbon?

Useful complementary diagnostics are:

| Diagnostic | Definition | What it can and cannot show |
| --- | --- | --- |
| Circular harmonics | `R_k = abs(sum(exp(i*k*phi))/N)`, for several `k` | Small values indicate weak low-order phase structure; `R1` alone misses two opposing lobes. |
| Largest angular gap | Sort wrapped phases; include the gap across zero | Measures unrepresented arcs at finite `N`; one remote tail particle can reduce a gap. |
| Occupied-bin fraction | Fraction of fixed equal-angle bins containing particles | Easy to inspect, but depends strongly on sample size and bin edges. |
| Histogram total variation | `0.5*sum(abs(count_b/N - 1/B))` in mean phase | Shows departure from a uniform phase marginal; has a positive sampling-noise floor. |
| Probability-content arc | Shortest circular arc containing a specified fraction, e.g. 99%, of particles | Describes bulk extent without using the most extreme particles. |
| Radial/normal thickness | Quantiles of signed deviations from a reference orbital curve | Distinguishes a thin ribbon from a thick torus or a dispersed plane family. |

The exact implemented diagnostics and configurable detection thresholds are
documented in the root README and DESCRIPTION. Diagnostic suggestions here are
not a claim that every listed extension is implemented.

A defensible default experiment can use 72 phase bins (5-degree resolution),
near-complete occupancy, a small maximum angular gap, and low circular
harmonics, and require the conditions at multiple successive output epochs.
Threshold choices are analysis conventions, not universal constants. A
“coverage” flag should not be renamed “uniform density” or “thin ribbon” without
the corresponding checks. Record the first observed crossing and the preceding
output epoch. That pair brackets a sampled threshold crossing; it is not proof
that no short earlier crossing occurred between outputs. Repeat with a finer
output grid near the event.

On a common eccentric Kepler ellipse, uniform mean anomaly gives

```text
p(true_anomaly=f) = (1-e^2)^(3/2) / (2*pi*(1+e*cos(f))^2).
```

This follows from Kepler's area law and `dM/df`. Consequently, density per unit
true anomaly peaks near apogee. Use mean longitude for phase-mixing tests and
physical position/true longitude for the visualization. If planes or apsides
have dispersed widely, different particles' individual true anomalies are not
angles along one common geometric orbit; inspect the three-dimensional cloud
and a stated reference-plane projection as well.

Once a cloud wraps, a Cartesian mean can be close to the Earth even though no
particle is near the Earth. A Cartesian covariance can also look broad in two
in-plane directions. Neither contradicts dominant early along-track shearing.
Continue to retain the particles and circular diagnostics rather than replacing
them with an ellipsoid or a single along-track standard deviation.

## How many particles are enough?

Use the required statistic to choose `N`, then check convergence. These are
starting recommendations, not accuracy guarantees:

- **2,000–5,000 particles:** pilot runs, visible shape, rough onset scale, and
  selection of the useful time interval.
- **20,000–50,000 particles:** normal studies of the phase marginal and ribbon
  onset after checking independent seeds and doubled sample count.
- **100,000 or more:** finer angular bins or moderate tail estimation. Rare-event
  probabilities may need substantially more particles or importance sampling.

For independent, equally weighted particles, the distribution-free
Dvoretzky–Kiefer–Wolfowitz–Massart bound for one scalar cumulative distribution
at a fixed time is

```text
P(sup_x |F_N(x)-F(x)| > epsilon) <= 2*exp(-2*N*epsilon^2)
N >= log(2/alpha)/(2*epsilon^2).
```

At 95% confidence, errors of 0.05, 0.02, and 0.01 require at least 738, 4,612,
and 18,445 samples respectively. This controls one-dimensional CDF error, not
three-dimensional density, a tail probability's relative error, or uncertainty
in a threshold time. A fixed angular branch cut allows a scalar phase CDF. For
simultaneous guarantees at `K` inspected times, a conservative union bound uses
`alpha/K`. Samples remain independent across particles under deterministic
propagation, although different times for the same particle are correlated.
[Massart's original paper](https://doi.org/10.1214/aop/1176990746),
[Wei and Dudley, primary research statement of the bound](https://arxiv.org/abs/1107.5356).

Useful self-derived sampling scales are:

```text
SE(probability estimate p_hat) ~= sqrt(p*(1-p)/N)
relative SE for a rare bin    ~= 1/sqrt(N*p)
E[R_k^2] under uniformity     = 1/N
```

For 72 uniform bins, 20,000 particles give about 278 particles per bin and
roughly 12% relative 95% sampling error in an individual bin. Therefore an
apparently uneven histogram is compatible with a uniform underlying phase
density. The large-`N` 95th percentile of a single uniform-sample harmonic
magnitude is approximately `sqrt(-log(0.05)/N)`: 0.0245 at `N=5,000`, or 0.0122
at `N=20,000`. A threshold of `R1=0.05` should be interpreted against this floor.

Under an exactly uniform distribution, a union bound gives
`P(any of B bins empty) <= B*exp(-N/B)`. Thus all-bin occupancy can occur with
only hundreds of particles; it does not establish uniformity. Under a strongly
nonuniform distribution, discovering low-probability bins can instead require
many more particles. For a tail probability of `1e-4`, 20,000 samples have only
two expected tail particles, far too few for a stable relative estimate.

For a final reported onset time, repeat at `N` and `2N`, use at least three
independent seeds, and refine the output cadence. Report the range of onset
times and the output-time bracket. That empirical range measures observed
Monte Carlo variability; it is not automatically a formal confidence interval.
Weighted posterior samples need weighted diagnostics and effective sample-size
reasoning. Deterministic or correlated sample designs do not inherit these
i.i.d. formulas unchanged.

## Force model and numerical accuracy boundaries

Negligible drag simplifies the question but does not remove Earth's nonspherical
gravity, lunisolar gravity, solar radiation pressure, resonances, or maneuvers.
J2 causes differential plane and apsidal precession; tesseral terms can matter
near resonances, and high area-to-mass objects may be sensitive to radiation
pressure. There is no single altitude or time span that makes every omitted
force negligible for every object. Results are conditional on the declared
force model and initial distribution. These are scientific sensitivity
considerations, not features silently added to a run.

The linked native dependency at `D:/orekit/DSST-cpp` is a C++17 port of Orekit
13.1.6 DSST, with its own lightweight orbit, date, frame, body, and force-model
adapters. Its README explicitly asks applications to supply external time,
frame, and celestial-ephemeris behavior and to validate mission-specific
propagation horizons. Its inspected `DSSTPropagator.hpp` integrates the mean
elements with fixed-step classical RK4. The presence of an integrator object
in its API must not be interpreted as proof that it implements Java Orekit's
adaptive integration behavior.

DSST separates slowly varying mean elements from short-period terms. Initial
orbit determination states are generally osculating: mark the state type
explicitly and use force-consistent conversion. Choose osculating output when
interpreting instantaneous physical geometry; mean output is useful for
secular analysis but does not include those short-period excursions. The
Java reference documentation describes these distinctions and warns against
sharing one propagator instance across threads. [Orekit DSSTPropagator documentation](https://www.orekit.org/site-orekit-latest/apidocs/org/orekit/propagation/semianalytical/dsst/DSSTPropagator.html).

Good numerical evidence includes: a Keplerian analytic check; separate testing
of sampling and circular statistics; comparison of successively smaller steps
for representative particles over the actual horizon; force-consistent
mean/osculating round trips; and agreement with an independently configured
reference for the model and horizon being claimed. Step convergence measures
integration error within the chosen model. It does not validate the physical
force model or prove complete C++/Java parity. A short-time upstream regression
fixture cannot by itself validate a multi-year study.

Performance improvements that preserve the defined calculation include
parallel independent particles, deterministic pre-generated samples, independent
mutable force models per worker, sequential output epochs for each particle,
immutable shared constants, bounded visualization subsampling, and avoiding
duplicate state conversion. Report benchmark configuration and accuracy checks
together. Increasing the integration step, dropping short-period terms, reducing
gravity degree, or simplifying an ephemeris changes the error/model budget and
requires explicit validation; it is not a free optimization.

## Interpretation of the deliverable

The simulator transports a specified initial distribution through a specified
dynamical model. It does not infer a real object's initial uncertainty from
measurements, add process noise automatically, or certify the truth of a
long-term prediction. Start with a controlled demonstration that shows phase
shearing, then insert the actual orbit and complete state covariance or posterior
samples. A statement such as “the distribution achieved the configured
5-degree coverage and mixing criteria between days X and Y” is reproducible.
“All uncertainty is uniformly spread around the orbit after day X” is generally
too strong.
