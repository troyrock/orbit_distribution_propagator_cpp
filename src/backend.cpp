#include "distribution/backend.hpp"

#include <dsst/DSSTPropagator.hpp>
#include <dsst/forces/DSSTJ2SquaredClosedForm.hpp>
#include <dsst/forces/DSSTZonal.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace distribution {
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
using Orbit = dsst::utilities::EquinoctialOrbitData;

void require_positive(double value, const char *name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive and finite");
    }
}

void validate_elements(const Elements &values) {
    for (double value : values) {
        if (!std::isfinite(value))
            throw std::invalid_argument("nonfinite orbit element");
    }
    if (values[0] <= 0.0 || std::hypot(values[1], values[2]) >= 1.0) {
        throw std::invalid_argument("DSST backend requires a bound elliptic orbit");
    }
}

Elements elements(const Orbit &orbit) {
    return {orbit.a, orbit.equinoctialEx, orbit.equinoctialEy, orbit.hx, orbit.hy, orbit.lm};
}

Orbit make_orbit(const Elements &values, double time, double mu) {
    validate_elements(values);
    Orbit orbit;
    orbit.date = time;
    orbit.frame = std::string{"inertial_equatorial"};
    orbit.a = values[0];
    orbit.equinoctialEx = values[1];
    orbit.equinoctialEy = values[2];
    orbit.hx = values[3];
    orbit.hy = values[4];
    orbit.lm = values[5];
    orbit.mu = mu;
    orbit.e = std::hypot(values[1], values[2]);
    orbit.keplerianMeanMotion = std::sqrt(mu / orbit.a) / orbit.a;
    orbit.keplerianPeriod = 2.0 * pi / orbit.keplerianMeanMotion;
    // Solve the equinoctial Kepler equation on a reduced angle. Carry its
    // correction to the original unwrapped longitude, preserving revolutions.
    const double lm = std::remainder(orbit.lm, 2.0 * pi);
    double correction = 0.0;
    bool converged = false;
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double angle = lm + correction;
        const double f2 = values[1] * std::sin(angle) - values[2] * std::cos(angle);
        const double f1 = 1.0 - values[1] * std::cos(angle) - values[2] * std::sin(angle);
        const double f0 = correction - f2;
        const double shift = 2.0 * f0 * f1 / (2.0 * f1 * f1 - f0 * f2);
        correction -= shift;
        if (std::abs(shift) <= 2e-15) {
            converged = true;
            break;
        }
    }
    if (!converged)
        throw std::runtime_error("equinoctial Kepler equation did not converge");
    orbit.le = orbit.lm + correction;
    const double reduced_le = lm + correction;
    const double numerator = values[1] * std::sin(reduced_le) - values[2] * std::cos(reduced_le);
    const double denominator = 1.0 + std::sqrt(1.0 - orbit.e * orbit.e) -
                               values[1] * std::cos(reduced_le) - values[2] * std::sin(reduced_le);
    orbit.lv = orbit.le + 2.0 * std::atan2(numerator, denominator);
    return orbit;
}

std::mutex &construction_mutex() {
    static std::mutex value;
    return value;
}
} // namespace

struct Backend::Impl {
    BackendConfig config;
    dsst::DSSTPropagator propagator{std::any{}, "OSCULATING"};
    std::shared_ptr<dsst::forces::DSSTZonal> zonal;
    std::shared_ptr<dsst::forces::DSSTJ2SquaredClosedForm> j2_squared;
    std::unique_ptr<dsst::forces::DSSTZonal::HansenObjects> zonal_hansen;
    Elements current{};
    Elements first_rate{};
    BackendStats statistics;
    double elapsed = 0.0;
    double next_step = 0.0;
    double longitude_compensation = 0.0;
    bool have_first_rate = false;

    Impl(const BackendConfig &input_config, const Elements &initial) : config(input_config) {
        require_positive(config.mu, "mu");
        require_positive(config.earth_radius_m, "Earth radius");
        require_positive(config.relative_tolerance, "relative tolerance");
        require_positive(config.absolute_tolerance_m, "semimajor-axis absolute tolerance");
        require_positive(config.absolute_tolerance_elements, "element absolute tolerance");
        require_positive(config.min_step_s, "minimum step");
        require_positive(config.max_step_s, "maximum step");
        if (config.min_step_s > config.max_step_s)
            throw std::invalid_argument("minimum step exceeds maximum step");
        if (!std::isfinite(config.j2) || config.j2 < 0.0)
            throw std::invalid_argument("J2 must be finite and nonnegative");
        if (config.force_model != "kepler" && config.force_model != "j2" &&
            config.force_model != "j2_j2sq") {
            throw std::invalid_argument("force_model must be kepler, j2, or j2_j2sq");
        }
        if (config.initial_type != "mean" && config.initial_type != "osculating")
            throw std::invalid_argument("initial_type must be mean or osculating");
        if (config.output_type != "mean" && config.output_type != "osculating")
            throw std::invalid_argument("output_type must be mean or osculating");
        validate_elements(initial);
        propagator.setMu(config.mu);
        if (config.force_model != "kepler") {
            dsst::forces::SphericalHarmonicsProviderData gravity;
            gravity.ae = config.earth_radius_m;
            gravity.mu = config.mu;
            gravity.c20 = -config.j2;
            gravity.maxDegree = 2;
            gravity.maxOrder = 0;
            // Upstream DSSTZonal's constructor copies a process-global Vns
            // cache. Protect both population and copy, including constructors
            // that run after another worker has started propagating. All rates
            // and short-period calculations then use per-instance storage.
            {
                std::lock_guard<std::mutex> lock(construction_mutex());
                zonal = std::make_shared<dsst::forces::DSSTZonal>(gravity);
                propagator.addForceModel(zonal);
            }
            if (config.force_model == "j2_j2sq") {
                j2_squared = std::make_shared<dsst::forces::DSSTJ2SquaredClosedForm>(
                    dsst::forces::ZeisModel{}, gravity);
                propagator.addForceModel(j2_squared);
            }
        }
        auto orbit = make_orbit(initial, 0.0, config.mu);
        // The native shell does not initialize short-period terms before its
        // first computeMeanState call. Explicit initialization is essential.
        propagator.beforeIntegration(orbit);
        if (config.initial_type == "osculating") {
            orbit = propagator.computeMeanState(orbit, 1e-14, 200);
        }
        current = elements(orbit);
        validate_elements(current);
        propagator.beforeIntegration(make_orbit(current, 0.0, config.mu));
        if (zonal) {
            zonal_hansen = std::make_unique<dsst::forces::DSSTZonal::HansenObjects>(
                zonal->createHansenObjects());
        }
        next_step = std::min(config.max_step_s, std::max(config.min_step_s, 3600.0));
    }

    Elements rates(const Elements &state, double date) {
        const auto orbit = make_orbit(state, date, config.mu);
        const dsst::utilities::AuxiliaryElements auxiliary(orbit, dsst::DSSTPropagator::I);
        std::vector<double> native(6, 0.0);
        if (zonal) {
            // DSSTZonal::getMeanElementRate reconstructs invariant Hansen
            // polynomial tables on every call. Keep those tables per particle;
            // createUAnddU still recomputes every orbit-dependent Hansen root.
            // These are the same native method calls, in the same order, as the
            // uncached port. Degree-2 truncation sizes are fixed for all states.
            const auto context = zonal->initializeStep(auxiliary);
            const auto potential =
                zonal->createUAnddU(orbit.date, context, auxiliary, *zonal_hansen);
            native = zonal->computeMeanElementRates(context, potential);
        }
        if (j2_squared) {
            const auto second_order = j2_squared->getMeanElementRate(auxiliary);
            for (std::size_t i = 0; i < 6; ++i)
                native[i] += second_order[i];
        }
        // Keep the shell's summation order: zonal, J2-squared, Newtonian.
        const auto central =
            propagator.elementRates(propagator.getAllForceModels().back(), auxiliary);
        for (std::size_t i = 0; i < 6; ++i)
            native[i] += central[i];
        ++statistics.derivative_evaluations;
        Elements result{};
        for (std::size_t i = 0; i < result.size(); ++i) {
            result[i] = native.at(i);
            if (!std::isfinite(result[i]))
                throw std::runtime_error("nonfinite native DSST derivative");
        }
        return result;
    }

    BackendState advance(double target) {
        if (!std::isfinite(target) || target < elapsed)
            throw std::invalid_argument("output times must be finite and nondecreasing");
        // Dormand-Prince embedded 5(4), with first-same-as-last reuse. Array
        // storage keeps stage assembly allocation-free; forces remain the
        // exact native DSST models configured above.
        constexpr double a[7][7] = {
            {},
            {1.0 / 5.0},
            {3.0 / 40.0, 9.0 / 40.0},
            {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
            {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0},
            {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
            {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}};
        constexpr double c[7] = {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};
        constexpr double b4[7] = {
            5179.0 / 57600.0, 0.0,       7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0,
            187.0 / 2100.0,   1.0 / 40.0};
        std::size_t attempts = 0;
        while (elapsed < target) {
            if (++attempts > 10000000)
                throw std::runtime_error("integration step limit exceeded");
            const double remaining = target - elapsed;
            const double h = std::min({next_step, config.max_step_s, remaining});
            if (!(h > 0.0) || elapsed + h == elapsed)
                throw std::runtime_error("integration time resolution exhausted");
            std::array<Elements, 7> k{};
            if (!have_first_rate) {
                first_rate = rates(current, elapsed);
                have_first_rate = true;
            }
            k[0] = first_rate;
            Elements candidate{};
            Elements increment{};
            for (std::size_t stage = 1; stage < 7; ++stage) {
                Elements stage_state = current;
                for (std::size_t component = 0; component < 6; ++component) {
                    double sum = 0.0;
                    for (std::size_t j = 0; j < stage; ++j)
                        sum += a[stage][j] * k[j][component];
                    increment[component] = h * sum;
                    stage_state[component] += increment[component];
                }
                k[stage] = rates(stage_state, elapsed + c[stage] * h);
                if (stage == 6)
                    candidate = stage_state;
            }
            double error = 0.0;
            for (std::size_t component = 0; component < 6; ++component) {
                double weighted = 0.0;
                for (std::size_t stage = 0; stage < 7; ++stage) {
                    weighted += (a[6][stage] - b4[stage]) * k[stage][component];
                }
                // Do not relax the longitude tolerance as revolutions accumulate.
                const double magnitude = component == 5 ? 1.0
                                                        : std::max(std::abs(current[component]),
                                                                   std::abs(candidate[component]));
                const double absolute = component == 0 ? config.absolute_tolerance_m
                                                       : config.absolute_tolerance_elements;
                const double scale = absolute + config.relative_tolerance * magnitude;
                error = std::max(error, std::abs(h * weighted) / scale);
            }
            if (!std::isfinite(error))
                throw std::runtime_error("nonfinite integration error estimate");
            const double factor =
                error == 0.0 ? 5.0 : std::clamp(0.9 * std::pow(error, -0.2), 0.1, 5.0);
            if (error <= 1.0) {
                const double phase_delta = increment[5] - longitude_compensation;
                const double next_longitude = current[5] + phase_delta;
                longitude_compensation = (next_longitude - current[5]) - phase_delta;
                candidate[5] = next_longitude;
                current = candidate;
                elapsed = h == remaining ? target : elapsed + h;
                first_rate = k[6];
                ++statistics.accepted_steps;
                next_step = std::clamp(h * factor, config.min_step_s, config.max_step_s);
            } else {
                ++statistics.rejected_steps;
                if (h <= config.min_step_s || remaining < config.min_step_s) {
                    throw std::runtime_error("requested accuracy cannot be met at minimum step");
                }
                next_step = std::max(config.min_step_s, h * std::min(0.9, factor));
            }
        }
        BackendState result{current, current};
        if (config.output_type == "osculating") {
            const auto mean_orbit = make_orbit(current, elapsed, config.mu);
            // Evaluate short-period coefficients at the exact requested date.
            // No sparse-output interpolation and no repeated conversion back to
            // mean state occurs. Initialization clears old interpolation slots.
            propagator.beforeIntegration(mean_orbit);
            result.osculating = elements(propagator.computeOsculatingOrbit(mean_orbit));
            validate_elements(result.osculating);
        }
        return result;
    }
};

Backend::Backend(const BackendConfig &config, const Elements &initial)
    : impl_(std::make_unique<Impl>(config, initial)) {}
Backend::~Backend() = default;
Backend::Backend(Backend &&) noexcept = default;
Backend &Backend::operator=(Backend &&) noexcept = default;
BackendState Backend::advance(double elapsed_seconds) {
    return impl_->advance(elapsed_seconds);
}
const BackendStats &Backend::stats() const noexcept {
    return impl_->statistics;
}
double Backend::time() const noexcept {
    return impl_->elapsed;
}

} // namespace distribution
