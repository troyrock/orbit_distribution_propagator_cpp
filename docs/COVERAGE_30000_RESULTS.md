# Coverage explorer through 30,000 km

Completed 2026-09-25 with the native C++ coverage executable and 20,000 particles
per case. The standalone HTML opens on the original MEO/PDF reference surface;
moving a slider selects the expanded orbital grid.

## Measured grid

- Perigee altitude: 1,000, 2,000, 3,000, 7,000, 12,500, 20,000,
  24,599.072067163703, and 30,000 km.
- Inclination: 0, 45, 90 degrees. Eccentricity: 0, 0.05, 0.1.
- Position sigma: 11 logarithmic nodes from 1 to 100 km.
- Velocity sigma: 11 logarithmic nodes from 0.01 to 1 m/s.
- Independent Gaussian per-axis RTN uncertainties, seed 20260919, initial
  osculating states, RAAN 20 degrees, argument of perigee 30 degrees, M 10 degrees.

The grid has **8,712 cases and 174,240,000 initial particle draws**. It reuses the
3,267 original low-altitude cases and adds 5,445 cases. There are 8,679 valid
coverage results and 33 whole-ensemble Earth-intersection masks. Masked particles
are not removed or redrawn to manufacture a valid result.

| Result | Days | Perigee km | Inclination | Eccentricity | Position sigma km | Velocity sigma m/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Minimum valid coverage | 0.252419 | 1,000 | 90 | 0.05 | 100 | 1 |
| Maximum valid coverage | 1553.451120 | 30,000 | 90 | 0.1 | 1 | 0.015848931924611134 |
| Original MEO/PDF maximum | 643.539443 | 19,650.663 | 55 | 0.02 | 1 | 0.025118864315095794 |

The original reference retains all 121 measured values and their original output
epochs. It is separate from the expanded grid's interpolation. It uses semimajor
axis 26,560 km, explaining the longer coverage times than the earlier
1,000-3,000 km explorer. A matched-input check with the explorer executable
reproduced its 643.539443-day maximum and onset bracket.

## Interpolation and refinement

The sliders interpolate log time in log geocentric perigee radius, linear
inclination, and linear eccentricity. Intermediate settings are labeled as
interpolated; unavailable supporting corners remain gaps.

The first extended grid produced approximately 5.75% errors at two withheld
cases near 24,599 km, e=0.025, position sigma 1 km, velocity sigma 1 m/s.
Targeted native runs isolated an altitude-dependent variation in the finite
sample coverage event. A complete altitude plane was added there. The original
holdouts were retained, and 72 new cases checked both smaller altitude intervals.
The 5% acceptance threshold was unchanged.

Of 288 extension/refinement holdouts, 282 have valid supporting corners:

- Maximum relative interpolation error: **4.504036%**.
- 95th percentile: **4.285983%**.
- RMS relative error: **1.740396%**.
- Comparisons above 5%: **0**.

The earlier 72-case holdout set was checked independently under the new
interpolation policy. Its 66 supported comparisons pass, with maximum error
4.454427%. Each set has six comparisons without valid supporting corners;
these remain unavailable and are excluded from error statistics. These checks
validate the sampled holdouts, not a global error bound over every slider value.

## Verification and model scope

All 45 focused Python tests pass, covering grid reuse/resume, interpolation,
rendering, and packaging provenance. Actual-report and synthetic browser checks
pass for slider endpoints, independent interpolation oracles, reference selection,
CSV export, fitted/common scales, masks, camera controls, responsive layouts,
and offline operation. The final screenshots were inspected visually.

The native propagation executable was unchanged. SHA-256:
`1fac5cd3856605ddc43c448cefd4d119c4925d10afe2f90039cb627799a7e361`.
Source datasets preserve executable/configuration hashes and per-case provenance.

Coverage requires all 72 relative mean-longitude bins occupied, largest gap at
most 5 degrees, and persistence for three output epochs. The reported time is
the first qualifying epoch, not the later confirmation. Sampling brackets are
temporal-resolution intervals, not Monte Carlo confidence intervals. Coverage
does not establish uniform phase probability or a particular tube thickness.

The model remains native DSST J2 plus J2-squared, with exact mean-phase flow after
native initialization. Drag, higher harmonics, third bodies, radiation pressure,
maneuvers, and process noise are omitted. Long-span results are conditional
experiments under that model.

See [the study guide](COVERAGE_EXPLORER.md) and
[the README](../README.md) for build and usage instructions. The report package
contains its complete dataset, measured CSV, validation JSON, provenance summary,
checksums, and a shareable HTML-only viewing experience.
