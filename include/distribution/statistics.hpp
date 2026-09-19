#pragma once
#include <distribution/orbit.hpp>
#include <cstddef>
#include <vector>
namespace distribution {
struct PhaseMetrics {
    double sigma=0, q95_width=0, max_gap=0, occupied_fraction=0, total_variation=0;
    std::array<double,4> resultants{};
    std::vector<std::size_t> histogram;
};
double quantile_sorted(const std::vector<double>& sorted,double fraction);
PhaseMetrics phase_metrics(const std::vector<double>& unwrapped_phase,std::size_t bins);
double sample_sigma(const std::vector<double>& values);
// Distribution-free simultaneous CDF error bound, P(sup |Fn-F| > eps) <= alpha.
double dkw_epsilon(std::size_t n,double alpha=0.05);
struct EventInterval {bool found=false;double lower_s=0,upper_s=0;std::size_t confirmed_at=0;};
EventInterval first_sustained(const std::vector<double>& times,const std::vector<bool>& flags,std::size_t persistence);
}
