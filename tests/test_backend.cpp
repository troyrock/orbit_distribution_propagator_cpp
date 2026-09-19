#include "distribution/backend.hpp"

#include <dsst/DSSTPropagator.hpp>
#include <dsst/forces/DSSTJ2SquaredClosedForm.hpp>
#include <dsst/forces/DSSTZonal.hpp>

#include <algorithm>
#include <cmath>
#include <future>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
constexpr double pi = 3.1415926535897932384626433832795;
using distribution::Backend;
using distribution::BackendConfig;
using distribution::Elements;
using Orbit = dsst::utilities::EquinoctialOrbitData;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr.precision(17);
        std::cerr << message << ": actual=" << actual << " expected=" << expected
                  << " tolerance=" << tolerance << '\n';
        throw std::runtime_error(message);
    }
}

Orbit orbit_from(const Elements& y, double mu, double time = 0.0) {
    Orbit o;
    o.date = time;
    o.frame = std::string{"inertial_equatorial"};
    o.mu = mu;
    o.a = y[0]; o.equinoctialEx = y[1]; o.equinoctialEy = y[2];
    o.hx = y[3]; o.hy = y[4]; o.lm = y[5];
    o.e = std::hypot(y[1], y[2]);
    o.keplerianMeanMotion = std::sqrt(mu / o.a) / o.a;
    o.keplerianPeriod = 2.0 * pi / o.keplerianMeanMotion;
    o.le = o.lm;
    for (int j = 0; j < 20; ++j) {
        const double shift = (o.le - y[1] * std::sin(o.le) + y[2] * std::cos(o.le) - o.lm)
                            / (1.0 - y[1] * std::cos(o.le) - y[2] * std::sin(o.le));
        o.le -= shift;
        if (std::abs(shift) < 1e-15) break;
    }
    o.lv = o.le + 2.0 * std::atan2(y[1] * std::sin(o.le) - y[2] * std::cos(o.le),
        1.0 + std::sqrt(1.0 - o.e * o.e) - y[1] * std::cos(o.le) - y[2] * std::sin(o.le));
    return o;
}

Elements values(const Orbit& o) { return {o.a, o.equinoctialEx, o.equinoctialEy, o.hx, o.hy, o.lm}; }

void configure_native(dsst::DSSTPropagator& native, const BackendConfig& config) {
    native.setMu(config.mu);
    dsst::forces::SphericalHarmonicsProviderData gravity;
    gravity.mu = config.mu; gravity.ae = config.earth_radius_m;
    gravity.c20 = -config.j2; gravity.maxDegree = 2; gravity.maxOrder = 0;
    native.addForceModel(std::make_shared<dsst::forces::DSSTZonal>(gravity));
    if (config.force_model == "j2_j2sq") {
        native.addForceModel(std::make_shared<dsst::forces::DSSTJ2SquaredClosedForm>(dsst::forces::ZeisModel{}, gravity));
    }
}

Elements sample() { return {26560000.0, 0.01, -0.004, 0.4, 0.3, 0.8}; }

void kepler_longitude_and_invariants() {
    BackendConfig config;
    config.force_model = "kepler";
    const auto initial = sample();
    Backend backend(config, initial);
    const double n = std::sqrt(config.mu / initial[0]) / initial[0];
    for (double date : {0.0, 86400.0, 365.25 * 86400.0, 20.0 * 365.25 * 86400.0}) {
        const auto state = backend.advance(date);
        for (std::size_t j = 0; j < 5; ++j) near(state.mean[j], initial[j], 1e-14, "Kepler invariant");
        near(state.mean[5], initial[5] + n * date, 3e-11, "continuous Kepler longitude");
        for (std::size_t j = 0; j < 6; ++j) near(state.mean[j], state.osculating[j], 1e-14, "Kepler mean equals osculating");
    }
    require(backend.stats().accepted_steps > 0, "missing integrator statistics");
    require(backend.stats().derivative_evaluations == 1 + 6 * (backend.stats().accepted_steps + backend.stats().rejected_steps), "FSAL evaluation count");
}

void conversion_roundtrip() {
    for (const std::string force : {"j2", "j2_j2sq"}) {
      for (const Elements initial : {sample(),
          Elements{26560000.0, 0.0, 0.0, 0.4, 0.3, 0.8},
          Elements{30000000.0, 0.6, 0.1, 0.99, 0.05, 2.8},
          Elements{30000000.0, 0.02, -0.04, 1.7, 0.1, -2.0}}) {
        BackendConfig config;
        config.force_model = force;
        config.initial_type = "mean";
        Backend from_mean(config, initial);
        const auto converted = from_mean.advance(0.0);
        require(std::abs(converted.osculating[0] - converted.mean[0]) > 1.0, "short-period corrections missing");
        config.initial_type = "osculating";
        Backend from_osculating(config, converted.osculating);
        const auto restored = from_osculating.advance(0.0);
        for (std::size_t j = 0; j < 6; ++j) {
            const double tolerance = j == 0 ? 5e-6 : 5e-13;
            near(restored.mean[j], initial[j], tolerance, "osculating-to-mean inverse");
            near(restored.osculating[j], converted.osculating[j], tolerance, "initial osculating state preserved");
        }
      }
    }
}

void native_rk4_comparison() {
    for (const std::string force : {"j2", "j2_j2sq"}) {
        BackendConfig config;
        config.force_model = force;
        config.initial_type = "mean";
        Backend adaptive(config, sample());
        dsst::DSSTPropagator native;
        configure_native(native, config);
        auto initial = orbit_from(sample(), config.mu);
        native.setInitialState(initial, "MEAN");
        const double target = 14.0 * 86400.0;
        const auto expected_orbit = native.propagate(target, 300.0);
        const auto expected = values(expected_orbit);
        const auto actual = adaptive.advance(target);
        for (std::size_t j = 0; j < 6; ++j) {
            near(actual.mean[j], expected[j], j == 0 ? 1e-5 : 5e-11, "native fixed-step RK4 comparison");
        }
        dsst::DSSTPropagator conversion(std::any{}, "OSCULATING");
        configure_native(conversion, config);
        conversion.beforeIntegration(expected_orbit);
        const auto expected_osculating = values(conversion.computeOsculatingOrbit(expected_orbit));
        for (std::size_t j = 0; j < 6; ++j) {
            near(actual.osculating[j], expected_osculating[j], j == 0 ? 5e-5 : 5e-11, "native osculating output comparison");
        }
    }
}

void tolerance_and_output_cadence() {
    BackendConfig normal;
    normal.initial_type = "mean";
    BackendConfig tight = normal;
    tight.relative_tolerance *= 0.01;
    tight.absolute_tolerance_elements *= 0.01;
    tight.absolute_tolerance_m *= 0.01;
    tight.max_step_s = 3600.0;
    Backend one_output(normal, sample());
    Backend many_outputs(normal, sample());
    Backend reference(tight, sample());
    const double target = 365.25 * 86400.0;
    const auto actual = one_output.advance(target);
    const auto truth = reference.advance(target);
    distribution::BackendState scheduled;
    for (int j = 1; j <= 100; ++j) scheduled = many_outputs.advance(target * static_cast<double>(j) / 100.0);
    for (std::size_t j = 0; j < 6; ++j) {
        const double tolerance = j == 0 ? 1e-3 : 3e-9;
        near(actual.mean[j], truth.mean[j], tolerance, "one-year tightened-tolerance convergence");
        near(scheduled.mean[j], truth.mean[j], tolerance, "output cadence does not reconvert mean state");
    }
    std::cout << "One-year convergence: phase error=" << std::abs(actual.mean[5] - truth.mean[5])
              << " rad, steps=" << one_output.stats().accepted_steps
              << ", reference steps=" << reference.stats().accepted_steps << '\n';
}

void concurrent_backends() {
    BackendConfig config;
    config.initial_type = "mean";
    Backend serial(config, sample());
    const auto expected = serial.advance(86400.0 * 30.0);
    std::vector<std::future<distribution::BackendState>> jobs;
    for (int j = 0; j < 8; ++j) {
        jobs.push_back(std::async(std::launch::async, [config]() {
            Backend backend(config, sample());
            return backend.advance(86400.0 * 30.0);
        }));
    }
    for (auto& job : jobs) {
        const auto actual = job.get();
        require(actual.mean == expected.mean && actual.osculating == expected.osculating,
                "parallel backend output differs from serial");
    }
}

void invalid_configuration() {
    BackendConfig config;
    config.max_step_s = 0.0;
    bool caught = false;
    try { Backend backend(config, sample()); } catch (const std::invalid_argument&) { caught = true; }
    require(caught, "zero maximum step accepted");
    config = BackendConfig{};
    Backend backend(config, sample());
    backend.advance(1.0);
    caught = false;
    try { backend.advance(0.0); } catch (const std::invalid_argument&) { caught = true; }
    require(caught, "backward output time accepted");
}
} // namespace

int main() {
    try {
        kepler_longitude_and_invariants();
        conversion_roundtrip();
        native_rk4_comparison();
        tolerance_and_output_cadence();
        concurrent_backends();
        invalid_configuration();
        std::cout << "DSST backend checks passed\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Backend test failed: " << exception.what() << '\n';
        return 1;
    }
}
