# Distribution propagator

A C++ Monte Carlo experiment for the evolution of Earth-orbit uncertainty,
using the native Orekit DSST port in `D:/orekit/DSST-cpp`.

Development is in progress. The completed guide will include build commands,
configuration conventions, examples, accuracy checks and an offline HTML viewer.

## Scope

The simulator follows a fixed set of initial orbit samples without measurement
updates or atmospheric drag. It distinguishes orbital coverage from phase
mixing and retains nonlinear samples after a single Cartesian covariance loses
its geometric meaning. Initial Gaussian uncertainty and empirical ensembles
will be supported. Mean and osculating elements are explicitly distinguished.

The dependency is external and read-only. Its audited starting commit is
`e653cd40e5e057cb89395afd8a9e3012ad59b12c` (Orekit 13.1.6 C++ port).
