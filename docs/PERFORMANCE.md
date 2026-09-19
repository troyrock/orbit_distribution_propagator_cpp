# Measured native DSST cache optimization

Retaining each particle's native Hansen polynomial tables removed repeated
construction from every derivative evaluation. The numerical methods, force
configuration, adaptive tolerances, and output reconstruction were unchanged.
For the focused run below, 100 osculating snapshots over one year became about
10.9 times faster. All 32 independent Java-reference comparison outputs were
byte-identical before and after the change.

## Workload and measurements

These measurements were collected on 2026-09-19 using Windows amd64 and MSYS2
UCRT GCC 15.2.0, with `-O3`, on the local machine. Each timing advances 32
particles **sequentially**, over 365.25 days, with `force_model=j2_j2sq`, mean
initial elements, and otherwise default `BackendConfig` settings:

```text
particle i = [26560000 + 10*i, 0.01, -0.004, 0.4, 0.3, 0.8]
i = 0,...,31
output epochs j = (365.25 days)*j/frames, j = 1,...,frames
rtol = 1e-11; atol_a = 1e-3 m; atol_other = 1e-12
minimum step = 0.001 s; maximum step = 86400 s
```

The timer includes backend construction and propagation for all 32 particles.
It excludes compilation and output-file writing. The snapshots in this workload
exclude the initial epoch. These are backend timings, not full simulator runs:
they omit ensemble statistics, Cartesian visualization conversion, and HTML/JSON
export. They do not measure multithread speedup.

| Output | Frames per particle | Original seconds | Cached seconds | Speedup | Derivative evaluations, both versions |
| --- | ---: | ---: | ---: | ---: | ---: |
| Mean | 1 | 3.1682438 | 0.2275215 | 13.9x | 70,496 |
| Mean | 20 | 3.3490146 | 0.2269346 | 14.8x | 73,568 |
| Mean | 100 | 3.6142696 | 0.2418072 | 14.9x | 77,216 |
| Mean | 500 | 4.8150377 | 0.2922378 | 16.5x | 96,416 |
| Osculating | 1 | 3.5583923 | 0.2173580 | 16.4x | 70,496 |
| Osculating | 20 | 3.9621281 | 0.2589502 | 15.3x | 73,568 |
| Osculating | 100 | 4.1535569 | 0.3819979 | 10.9x | 77,216 |
| Osculating | 500 | 5.8446998 | 0.9985708 | 5.85x | 96,416 |

These were single, sequential measurements under concurrent machine activity,
not a controlled repeated benchmark with confidence intervals. The apparent
differences between one-frame mean and osculating timings illustrate that noise.
Treat the factors as workload-specific observations rather than portable speed
guarantees. After caching, exact-date short-period reconstruction accounts for
much more of the cost when many osculating snapshots are requested.

A separate diagnostic repeatedly evaluated the native zonal rate for one fixed
state 100,000 times: the original call took 4.49814 seconds and the cached native
method sequence took 0.12233 seconds, about 36.8x faster. That narrower diagnostic
excludes integrator, conversion, and snapshot costs; it is not the application
speedup. The end-to-end table above is the useful comparison.

## What is cached and what is recomputed

The original native `DSSTZonal::getMeanElementRate` calls
`createHansenObjects()` on every rate evaluation. The resulting object constructs
polynomial tables whose structure depends on the configured harmonic degree and
eccentricity truncation. For the supported degree-2 profile those sizes are fixed.

The adapter creates these tables once per particle. Each rate evaluation still
calls the original native methods:

1. `DSSTZonal::initializeStep` for the current auxiliary elements.
2. `DSSTZonal::createUAnddU` with the retained Hansen object. This recomputes the
   orbit-dependent Hansen roots and potential derivatives for the current state.
3. `DSSTZonal::computeMeanElementRates` for the zonal contribution.
4. `DSSTJ2SquaredClosedForm::getMeanElementRate` when selected, and the
   propagator's `elementRates` for Newtonian attraction.

Contributions retain the original summation order: zonal, J2-squared, Newtonian.
No tabulated trajectory, approximate secular replacement, frozen orbital rate,
larger error tolerance, reduced force set, or fewer integration stages is used.
The propagator still owns the native models and performs mean/osculating
conversion and short-period reconstruction.

Each particle owns its mutable force objects and Hansen roots. The only shared
mutable native cache reached by these supported profiles is populated/copied
during zonal construction; that operation is protected by the adapter's mutex.
Advancing separate instances therefore needs no shared rate-evaluation lock.
One instance must not be advanced concurrently. This audit does **not** establish
thread safety for arbitrary upstream tesseral/third-body configurations: their
other shared caches and truncation changes need review before adding support.

## Numerical preservation

The before/after comparison used `tools/backend_probe.cpp`: J2 and J2-plus-J2²,
mean and osculating input, 0/1/30/365-day epochs, and both output types. The full
32-row CSV, including orbital elements and Cartesian states, was byte-identical:

```text
SHA-256 e6bf45bb116d3df998207c9a306fb68c5b72b2139efc74f90a586c32f22ce2e9
```

The independent Java comparison remained at a maximum 8.781 mm position
difference, against the declared 20 mm regression gate. See
[the fixture provenance](../tests/data/README.md). Native RK4 comparisons for
moderate-eccentricity, high-eccentricity, and circular retrograde geometries,
conversion roundtrips, one-year tighter-tolerance comparisons, and concurrent
instance tests also passed after caching. These checks establish preservation
for the tested configurations, not a physical accuracy guarantee for all orbits.

## Reproduce the focused workload

Build the project in Release mode first as described in the README. The original
uncached backend is preserved in commit `ed74879`. The following instructions
compile that source and the current cached source into separate workspace output
files; they do not change the working source or the upstream DSST checkout.
The commands below use the GCC/MinGW library names from the measured build.

Save the following as `outputs/performance/benchmark.cpp`, creating its parent
directory first:

```cpp
#include "distribution/backend.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

int main() {
    std::cout << "output,frames,particles,seconds,rate_evaluations,snapshots\n";
    for (const char* output : {"mean", "osculating"}) {
        for (int frames : {1, 20, 100, 500}) {
            distribution::BackendConfig config;
            config.initial_type = "mean";
            config.output_type = output;
            std::size_t evaluations = 0;
            constexpr int particles = 32;
            const auto start = std::chrono::steady_clock::now();
            for (int i = 0; i < particles; ++i) {
                distribution::Backend backend(config,
                    {26560000.0 + i * 10.0, 0.01, -0.004, 0.4, 0.3, 0.8});
                for (int j = 1; j <= frames; ++j) {
                    backend.advance(365.25 * 86400.0 * j / frames);
                }
                evaluations += backend.stats().derivative_evaluations;
            }
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            std::cout << output << ',' << frames << ',' << particles << ','
                      << std::setprecision(8) << seconds << ',' << evaluations
                      << ',' << frames * particles << '\n';
        }
    }
}
```

From the repository root in PowerShell, adjust the compiler/DSST paths for your
machine. Use the same compiler and ABI as the Release build:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
git show ed74879:src/backend.cpp |
    Set-Content -Encoding utf8 outputs/performance/backend_uncached.cpp

foreach ($variant in @('uncached', 'cached')) {
    $source = if ($variant -eq 'uncached') {
        'outputs/performance/backend_uncached.cpp'
    } else { 'src/backend.cpp' }
    $object = "outputs/performance/backend_$variant.o"
    & g++ -std=c++17 -O3 -Iinclude -ID:/orekit/DSST-cpp/include `
        -c $source -o $object
    if ($LASTEXITCODE -ne 0) { throw "Backend compilation failed: $variant" }
    & g++ -std=c++17 -O3 -Iinclude outputs/performance/benchmark.cpp `
        $object build/dsst/libdsst_cpp.a `
        -o "outputs/performance/benchmark_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw "Benchmark compilation failed: $variant" }
    & ".\outputs\performance\benchmark_$variant.exe" |
        Tee-Object -FilePath "outputs/performance/benchmark_$variant.csv"
    & g++ -std=c++17 -O3 -Iinclude tools/backend_probe.cpp $object `
        build/libdistribution_core.a build/dsst/libdsst_cpp.a `
        -o "outputs/performance/probe_$variant.exe"
    if ($LASTEXITCODE -ne 0) { throw "Probe compilation failed: $variant" }
    & ".\outputs\performance\probe_$variant.exe" `
        "outputs/performance/probe_$variant.csv" tests/data/orekit_long_horizon.csv
    if ($LASTEXITCODE -ne 0) { throw "Parity comparison failed: $variant" }
}
Get-FileHash outputs/performance/probe_uncached.csv, `
    outputs/performance/probe_cached.csv -Algorithm SHA256
```

For a controlled performance study, run repeated alternated variants on an idle
machine and report median timings and spread. Keep the numerical parity check
and identical derivative counts alongside the timing result.
