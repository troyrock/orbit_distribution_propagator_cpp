#pragma once
#include <distribution/config.hpp>
#include <distribution/statistics.hpp>
#include <filesystem>
namespace distribution {
struct SampleState {
    Elements output{};
    double mean_longitude = 0, mean_a = 0;
};
struct Frame {
    double time_s = 0;
    BackendState reference{};
    PhaseMetrics phase;
    std::array<double, 3> rtn_sigma_m{}, tube_rtn_sigma_m{};
    double semimajor_sigma_m = 0;
    bool coverage = false, mixed = false;
    std::vector<Vector3> positions, reference_orbit;
};
struct Simulation {
    Config config;
    std::vector<Elements> initial;
    std::vector<Frame> frames;
    std::vector<SampleState> states; // frame-major, all particles
    EventInterval coverage, mixing, central95_wrap;
    double analytic_mixing_time_s = -1, elapsed_seconds = 0, accuracy_max_position_m = 0,
           accuracy_max_phase_rad = 0;
    std::size_t threads_used = 0;
    BackendStats stats;
};
Simulation run_simulation(Config config);
void write_results(const Simulation &simulation, const std::filesystem::path &output);
} // namespace distribution
