#pragma once
#include <distribution/config.hpp>
#include <distribution/statistics.hpp>
#include <array>
#include <string>
#include <vector>

namespace distribution {
// Exact coverage Boolean when every equally spaced bin must be occupied and
// the permitted gap is at least a bin width. Falls back to sorted gaps for
// other thresholds or floating-point bin-edge ambiguities.
bool covers_orbit(const std::vector<double> &phase, std::size_t bins, double gap_radians,
                  double occupied_fraction = 1.0);

struct CoverageOptions {
    // Automatic grid uses the sampled native phase-rate sigma: cadence=1/(120*sigma),
    // initial horizon=192 cadences. Fixed mode uses Config's duration/cadence.
    bool automatic_times = true;
    std::size_t max_extensions = 3;
    double intervals_per_phase_sigma = 120;
};
struct CoverageResult {
    std::string status = "unobserved";
    std::size_t samples = 0, earth_intersections = 0, particles_below_300km = 0,
                particles_below_500km = 0, particles_below_1000km = 0,
                evaluated_epochs = 0, threads_used = 0;
    double min_initial_perigee_km = 0, min_mean_perigee_km = 0, cadence_s = 0,
           horizon_s = 0, phase_rate_sigma = 0, elapsed_seconds = 0;
    EventInterval coverage;
    double confirmation_s = 0;
    // Exact sorted maximum gap at the reported onset (not the bin bound).
    double onset_max_gap_deg = 0;
};
CoverageResult run_coverage(Config config, const CoverageOptions &options = {});
} // namespace distribution
