# Development checkpoints

1. Audit the C++ DSST dependency; define coordinate, sampling and mixing conventions.
2. Implement and test nonlinear sampling, orbital conversions and circular statistics.
3. Integrate an accuracy-controlled DSST backend and deterministic ensemble execution.
4. Add an offline, interactive HTML viewer and functional CLI tests.
5. Run accuracy, sample-size, long-span and throughput experiments; document evidence.

Use separate informative commits for tested development checkpoints. Do not
change upstream numerical fixtures or relax tolerances to obtain a passing gate.

## Initial dependency findings

The translated propagator's stepping routine uses fixed-step RK4. Its public
DSST derivative and short-period APIs permit an adaptive integration driver
without replacing the force theory. Long-duration validation must therefore
include tolerance refinement and comparison against an independent reference,
and cannot be inferred from short Java fixture agreement alone.

## Completed validation checkpoints

- `b08e68c`: isolated repository and declared scientific/validation scope.
- `395d4a2`: orbit conversion, Gaussian/empirical sampling, PSD covariance,
  circular diagnostics, and direct mathematical/statistical tests.
- `ed74879`: adaptive native-force DSST backend and independent one-year Java
  reference with original source dependency unchanged.
- `9cf774c`: parallel simulation/CLI, offline viewer, native Hansen cache,
  reproducible outputs and compiler-matched Windows runtime DLL deployment.
- Subsequent checkpoint: sample/seed/cadence driver, dual-compiler gates,
  empirical/unchecked-result metadata, source formatting and recorded results.

Final measured evidence and limitations: [VALIDATION.md](VALIDATION.md).
