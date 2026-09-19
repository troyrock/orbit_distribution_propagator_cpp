#pragma once
#include <distribution/backend.hpp>
#include <distribution/orbit.hpp>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>
namespace distribution {
struct Config {
    BackendConfig backend;
    Elements nominal = from_keplerian(26560000,0.02,55*pi/180,0.3,0.4,0.2);
    std::string uncertainty_coordinates = "equinoctial";
    Matrix6 covariance{};
    std::filesystem::path empirical_samples_csv;
    std::size_t samples=5000,threads=0,visual_samples=1500,phase_bins=72,persistence=3;
    std::uint64_t seed=20260919;
    double duration_days=120,output_step_days=1,minimum_altitude_m=1000000;
    double coverage_max_gap_deg=5,coverage_occupied_fraction=1,mixing_max_resultant=0.05;
    double mixing_max_tv=0.15,max_memory_mb=2048;
    std::size_t accuracy_check=0;
    bool write_html=true,export_states=false;
    Config();
};
Config read_config(const std::filesystem::path& path);
void validate_config(const Config& config);
std::vector<double> output_times(const Config& config);
long double estimated_memory_bytes(const Config& config,std::size_t sample_count);
std::vector<Elements> sample_initial(const Config& config);
std::string config_template();
}
