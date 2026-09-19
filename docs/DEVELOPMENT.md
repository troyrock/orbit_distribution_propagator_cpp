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
