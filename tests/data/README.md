# Independent Orekit long-horizon fixture

`orekit_long_horizon.csv` was generated on 2026-09-19 by the original Java
Orekit 13.1.6 DSST implementation, independently of the C++ port and its
adaptive integrator. It contains 32 cases: both supported perturbing profiles
(`j2`, `j2_j2sq`), mean and osculating initial elements, 0/1/30/365-day
horizons, and mean and osculating output.

Input is equinoctial `[26560000 m, 0.01, -0.004, 0.4, 0.3, 0.8 rad]`, in
GCRF at J2000. Constants are `mu=3.986004418e14 m^3/s^2`, Earth equatorial
radius `6378137 m`, and `J2=1.08262668e-3`. Gravity is unnormalized degree 2,
order 0, with `C20=-J2`. Zonal short-period truncations are degree 2,
eccentricity power 1, and frequency 5. The optional second-order term uses
Java `DSSTJ2SquaredClosedForm(new ZeisModel(), provider)`.

Reference integration uses Hipparchus Dormand-Prince 8(5,3), maximum step
3600 seconds, minimum step 0.001 seconds, relative tolerances 1e-14, and
absolute tolerances `[1e-6 m,1e-14,1e-14,1e-14,1e-14,1e-14,1e-12 kg]`.
The osculating-to-mean fixed point uses epsilon 1e-14 and at most 200
iterations. Output short-period terms are evaluated at the exact output
epoch with `DSSTPropagator.computeOsculatingState`, rather than relying on
interpolation across a large integration step.

The tool reads the prebuilt local Orekit jar and cached Hipparchus jars.
Generation requires no downloads, UTC/EOP files, external gravity files,
or writes to Orekit or the Maven cache. Dependency jars are copied to
ignored workspace output before Java reads them to avoid Windows restricted
account real-path errors during Java ZIP close.

Provenance verified for this generation:

- Java: Eclipse Adoptium 17.0.20+8, Windows amd64.
- Hipparchus: 4.0.3.
- Orekit jar SHA-256:
  `b6ac8d760d112378c27484f261a19c94e491f790b8236086d968a77cabe49b3c`.
- Fixture SHA-256:
  `3654858704493f70875032ea381aace1d04ddf662f2482248d3f7da6d2c93015`.
- C++ DSST checkout: `e653cd40e5e057cb89395afd8a9e3012ad59b12c`, clean.

Default C++ backend comparisons observed:

| Horizon | Largest position difference | Largest velocity difference |
| --- | ---: | ---: |
| Initial epoch | 3.54e-8 m | 6.00e-12 m/s |
| 1 day | 2.08e-7 m | 3.05e-11 m/s |
| 30 days | 1.63e-4 m | 2.36e-8 m/s |
| 365 days | 8.79e-3 m | 1.28e-6 m/s |

The new fixture gate is 0.02 m position, 3e-6 m/s velocity, 1e-4 m
semimajor axis, and 8e-10 for the other equinoctial elements. The largest
observed one-year longitude difference is 3.32e-10 radians. These are
declared numerical regression bounds for this orbit and force configuration;
they are not real-orbit accuracy promises, coverage of all orbit shapes, or
evidence that omitted lunar/solar/tesseral/radiation forces are negligible.
No existing upstream fixtures or tolerances were changed.

Generation and comparison instructions are in
[`tools/java_reference/README.md`](../../tools/java_reference/README.md).
