#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>

namespace distribution {

// SI equinoctial elements: a, e*cos(omega+Omega), e*sin(omega+Omega),
// tan(i/2)*cos(Omega), tan(i/2)*sin(Omega), unwrapped mean longitude.
using Elements = std::array<double, 6>;

struct BackendConfig {
    double mu = 3.986004418e14;
    double earth_radius_m = 6378137.0;
    double j2 = 1.08262668e-3;
    std::string force_model = "j2_j2sq"; // kepler, j2, j2_j2sq
    std::string initial_type = "osculating";
    std::string output_type = "osculating";
    double relative_tolerance = 1e-11;
    double absolute_tolerance_m = 1e-3;
    double absolute_tolerance_elements = 1e-12;
    double min_step_s = 1e-3;
    double max_step_s = 86400.0;
};

struct BackendState {
    Elements mean{};
    // With output_type="mean", this contains mean elements as well, avoiding
    // expensive short-period reconstruction. Consumers must retain the label.
    Elements osculating{};
};

struct BackendStats {
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t derivative_evaluations = 0;
};

// Exact mean-longitude flow for the backend's fixed central + degree-2 zonal
// J2 + optional Zeis J2-squared models. The coefficients come from the native
// DSST rate evaluation, after the normal osculating-to-mean conversion. This
// is not Kepler-only propagation or a fit, and does not give osculating phase.
struct MeanPhaseLaw {
    Elements mean_at_epoch{};
    double epoch_s = 0.0;
    double longitude_rate_rad_s = 0.0;

    // Absolute elapsed time from Backend construction, at or after epoch_s.
    // Returns unwrapped mean longitude, in radians.
    double longitude(double elapsed_seconds) const;
};

// Each particle owns its mutable native DSST force objects. Construction
// serializes the upstream zonal coefficient cache; distinct instances may
// subsequently advance concurrently. Do not call one instance concurrently.
class Backend {
  public:
    Backend(const BackendConfig &config, const Elements &initial);
    ~Backend();
    Backend(Backend &&) noexcept;
    Backend &operator=(Backend &&) noexcept;
    Backend(const Backend &) = delete;
    Backend &operator=(const Backend &) = delete;

    BackendState advance(double elapsed_seconds);
    // Evaluates one native derivative at the current mean state; neither
    // advances time nor changes that state. The derivative statistic increases
    // by one. Calling immediately after construction retains the initialized
    // invariants without any numerical integration drift.
    MeanPhaseLaw mean_phase_law();
    const BackendStats &stats() const noexcept;
    double time() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace distribution
