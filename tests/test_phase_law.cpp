#include "distribution/backend.hpp"
#include "distribution/orbit.hpp"

#include <dsst/DSSTPropagator.hpp>
#include <dsst/forces/DSSTJ2SquaredClosedForm.hpp>
#include <dsst/forces/DSSTZonal.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using distribution::Backend;
using distribution::BackendConfig;
using distribution::Elements;
using Orbit = dsst::utilities::EquinoctialOrbitData;
constexpr double day = 86400.0;
constexpr double degree = distribution::pi / 180.0;

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char *message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr.precision(17);
        std::cerr << message << ": actual=" << actual << " expected=" << expected
                  << " tolerance=" << tolerance << '\n';
        throw std::runtime_error(message);
    }
}

// Independent native-shell route: deliberately do not use Backend::rates or
// the cached Hansen implementation that supplies the proposed phase law.
Orbit make_native_orbit(const Elements &y, double mu, double date) {
    Orbit o;
    o.date = date;
    o.frame = std::string{"inertial_equatorial"};
    o.mu = mu;
    o.a = y[0];
    o.equinoctialEx = y[1];
    o.equinoctialEy = y[2];
    o.hx = y[3];
    o.hy = y[4];
    o.lm = y[5];
    o.e = std::hypot(y[1], y[2]);
    o.keplerianMeanMotion = std::sqrt(mu / o.a) / o.a;
    o.keplerianPeriod = distribution::tau / o.keplerianMeanMotion;
    const double lm = std::remainder(o.lm, distribution::tau);
    double le = lm;
    for (int iteration = 0; iteration < 30; ++iteration) {
        const double shift = (le - y[1] * std::sin(le) + y[2] * std::cos(le) - lm) /
                             (1.0 - y[1] * std::cos(le) - y[2] * std::sin(le));
        le -= shift;
        if (std::abs(shift) < 2e-15)
            break;
    }
    o.le = o.lm + le - lm;
    o.lv = o.le + 2.0 * std::atan2(y[1] * std::sin(le) - y[2] * std::cos(le),
                                   1.0 + std::sqrt(1.0 - o.e * o.e) -
                                       y[1] * std::cos(le) - y[2] * std::sin(le));
    return o;
}

void configure_native(dsst::DSSTPropagator &native, const BackendConfig &config) {
    native.setMu(config.mu);
    if (config.force_model == "kepler")
        return;
    dsst::forces::SphericalHarmonicsProviderData gravity;
    gravity.mu = config.mu;
    gravity.ae = config.earth_radius_m;
    gravity.c20 = -config.j2;
    gravity.maxDegree = 2;
    gravity.maxOrder = 0;
    native.addForceModel(std::make_shared<dsst::forces::DSSTZonal>(gravity));
    if (config.force_model == "j2_j2sq")
        native.addForceModel(std::make_shared<dsst::forces::DSSTJ2SquaredClosedForm>(
            dsst::forces::ZeisModel{}, gravity));
}

Elements native_rates(dsst::DSSTPropagator &native, const Elements &mean,
                      double mu, double date) {
    const auto orbit = make_native_orbit(mean, mu, date);
    native.beforeIntegration(orbit);
    const dsst::utilities::AuxiliaryElements auxiliary(orbit, dsst::DSSTPropagator::I);
    Elements result{};
    for (const auto &force : native.getAllForceModels()) {
        const auto rate = native.elementRates(force, auxiliary);
        for (std::size_t component = 0; component < 6; ++component)
            result[component] += rate[component];
    }
    return result;
}

Elements sample(double perigee_km, double inclination_deg, double eccentricity) {
    return distribution::from_keplerian((6378137.0 + perigee_km * 1000.0) /
                                            (1.0 - eccentricity),
                                        eccentricity, inclination_deg * degree,
                                        20.0 * degree, 30.0 * degree, 10.0 * degree);
}

Elements rotated(const Elements &initial, double eccentricity_rotation,
                 double node_rotation, double phase) {
    Elements result = initial;
    result[1] = initial[1] * std::cos(eccentricity_rotation) -
                initial[2] * std::sin(eccentricity_rotation);
    result[2] = initial[1] * std::sin(eccentricity_rotation) +
                initial[2] * std::cos(eccentricity_rotation);
    result[3] = initial[3] * std::cos(node_rotation) - initial[4] * std::sin(node_rotation);
    result[4] = initial[3] * std::sin(node_rotation) + initial[4] * std::cos(node_rotation);
    result[5] = phase;
    return result;
}

void native_invariant_and_rotation_checks() {
    for (const std::string force : {"kepler", "j2", "j2_j2sq"}) {
        BackendConfig config;
        config.initial_type = "mean";
        config.output_type = "mean";
        config.force_model = force;
        dsst::DSSTPropagator native;
        configure_native(native, config);
        for (double altitude : {1000.0, 2000.0, 3000.0}) {
            for (double inclination : {0.0, 45.0, 90.0}) {
                for (double eccentricity : {0.0, 0.05, 0.1}) {
                    const auto initial = sample(altitude, inclination, eccentricity);
                    Backend backend(config, initial);
                    const auto law = backend.mean_phase_law();
                    for (double angle : {0.0, 0.7, 2.5, 5.0}) {
                        const auto mean = rotated(initial, angle, -0.6 * angle,
                                                  initial[5] + 100000.0 * angle);
                        const auto rate = native_rates(native, mean, config.mu,
                                                       100.0 * day * angle);
                        near(rate[0], 0.0, 0.0, "native semimajor-axis derivative");
                        near(mean[1] * rate[1] + mean[2] * rate[2], 0.0, 2e-22,
                             "native eccentricity magnitude derivative");
                        near(mean[3] * rate[3] + mean[4] * rate[4], 0.0, 2e-21,
                             "native inclination magnitude derivative");
                        near(rate[5], law.longitude_rate_rad_s, 3e-18,
                             "native longitude rate after arbitrary precession");
                    }
                }
            }
        }
    }
}

double compare_advance(const BackendConfig &config, const Elements &initial,
                       double phase_tolerance) {
    Backend backend(config, initial);
    const auto law = backend.mean_phase_law();
    const auto first = backend.advance(0.0);
    require(law.mean_at_epoch == first.mean, "phase law bypassed native mean conversion");
    double maximum = 0.0;
    for (double days : {1.0, 30.0, 365.0, 1000.0, 2000.0}) {
        const auto state = backend.advance(days * day);
        const double error = std::abs(state.mean[5] - law.longitude(days * day));
        maximum = std::max(maximum, error);
        near(state.mean[5], law.longitude(days * day), phase_tolerance,
             "adaptive long-horizon mean phase");
        near(state.mean[0], law.mean_at_epoch[0], 1e-8,
             "adaptive semimajor-axis invariant");
        near(std::hypot(state.mean[1], state.mean[2]),
             std::hypot(law.mean_at_epoch[1], law.mean_at_epoch[2]), 2e-8,
             "adaptive eccentricity magnitude invariant");
        near(std::hypot(state.mean[3], state.mean[4]),
             std::hypot(law.mean_at_epoch[3], law.mean_at_epoch[4]), 2e-8,
             "adaptive inclination magnitude invariant");
        const auto later_law = backend.mean_phase_law();
        near(later_law.longitude_rate_rad_s, law.longitude_rate_rad_s, 2e-15,
             "native longitude rate after adaptive precession");
        near(later_law.epoch_s, days * day, 0.0, "phase law rebasing epoch");
        near(later_law.longitude(days * day), state.mean[5], 0.0,
             "phase law rebasing position");
    }
    return maximum;
}

void adaptive_long_horizon_checks() {
    double maximum = 0.0;
    BackendConfig config;
    config.output_type = "mean";
    // Cover the complete requested parameter range, including both singular
    // classical limits; the actual propagation remains nonsingular.
    for (double altitude : {1000.0, 2000.0, 3000.0})
        for (double inclination : {0.0, 45.0, 90.0})
            for (double eccentricity : {0.0, 0.05, 0.1})
                maximum = std::max(maximum, compare_advance(
                    config, sample(altitude, inclination, eccentricity), 3e-8));

    // Isolate each force choice and exercise native mean initialization as
    // well as the osculating initialization used by the coverage experiment.
    const Elements meo{26560000.0, 0.01, -0.004, 0.4, 0.3, 0.8};
    for (const std::string force : {"kepler", "j2", "j2_j2sq"}) {
        config.force_model = force;
        for (const std::string initial_type : {"mean", "osculating"}) {
            config.initial_type = initial_type;
            maximum = std::max(maximum, compare_advance(config, meo, 3e-8));
        }
    }

    config.force_model = "j2_j2sq";
    config.initial_type = "osculating";
    config.relative_tolerance = 1e-13;
    config.absolute_tolerance_elements = 1e-14;
    config.absolute_tolerance_m = 1e-5;
    config.max_step_s = 3600.0;
    double fine_maximum = 0.0;
    for (const Elements initial : {sample(1000.0, 0.0, 0.0),
                                    sample(1000.0, 45.0, 0.1),
                                    sample(3000.0, 90.0, 0.05), meo})
        fine_maximum = std::max(fine_maximum, compare_advance(config, initial, 3e-9));
    std::cout.precision(12);
    std::cout << "2000-day phase residual: default=" << maximum
              << " rad; tightened / max 1-hour steps=" << fine_maximum << " rad\n";
}

void api_contract_checks() {
    BackendConfig config;
    config.output_type = "mean";
    Backend backend(config, sample(2000.0, 45.0, 0.05));
    Backend baseline(config, sample(2000.0, 45.0, 0.05));
    const auto initial = backend.advance(0.0);
    const auto evaluations = backend.stats().derivative_evaluations;
    const auto law = backend.mean_phase_law();
    require(backend.time() == 0.0, "phase-law query advanced time");
    require(backend.stats().derivative_evaluations == evaluations + 1,
            "phase-law native derivative not accounted for");
    require(law.mean_at_epoch == initial.mean, "phase-law query changed mean state");
    require(backend.advance(30.0 * day).mean == baseline.advance(30.0 * day).mean,
            "phase-law query changed adaptive propagation");
    const auto rebased = backend.mean_phase_law();
    for (double invalid : {-1.0, std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::quiet_NaN()}) {
        bool caught = false;
        try {
            static_cast<void>(rebased.longitude(invalid));
        } catch (const std::invalid_argument &) {
            caught = true;
        }
        require(caught, "invalid phase-law query accepted");
    }
    bool caught = false;
    try {
        static_cast<void>(rebased.longitude(day));
    } catch (const std::invalid_argument &) {
        caught = true;
    }
    require(caught, "query before rebased phase-law epoch accepted");
}
} // namespace

int main() {
    try {
        api_contract_checks();
        native_invariant_and_rotation_checks();
        adaptive_long_horizon_checks();
        std::cout << "Native DSST constant mean-phase law checks passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Phase-law test failed: " << error.what() << '\n';
        return 1;
    }
}
