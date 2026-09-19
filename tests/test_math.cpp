#include <distribution/orbit.hpp>
#include <distribution/statistics.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace distribution;
constexpr double mu = 3.986004418e14;

// Throwing checks stay active when Release builds define NDEBUG.
void require(bool condition, const std::string& context) {
    if (!condition) throw std::runtime_error(context);
}

void near(double actual, double expected, double absolute_tolerance,
          double relative_tolerance, const std::string& context) {
    const double tolerance = absolute_tolerance + relative_tolerance * std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance) {
        std::ostringstream message;
        message << std::setprecision(17) << context << ": actual=" << actual
                << ", expected=" << expected << ", tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}

template<class Function>
void rejects(Function&& function, const std::string& context) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(context + ": expected std::invalid_argument");
}

Vector3 position(const Cartesian& state) { return {state[0], state[1], state[2]}; }
Vector3 velocity(const Cartesian& state) { return {state[3], state[4], state[5]}; }

void compare_states(const Cartesian& actual, const Cartesian& expected,
                    const std::string& context) {
    for (std::size_t k = 0; k < actual.size(); ++k) {
        near(actual[k], expected[k], k < 3 ? 3e-6 : 3e-9, 2e-12,
             context + " coordinate " + std::to_string(k));
    }
}

void test_circular_orbit() {
    constexpr double radius = 26560000.0;
    const double speed = std::sqrt(mu / radius);
    for (double phase : {0.0, pi / 2.0, pi, 3.0 * pi / 2.0}) {
        const Elements elements{radius, 0, 0, 0, 0, phase};
        const Cartesian expected{radius * std::cos(phase), radius * std::sin(phase), 0,
                                 -speed * std::sin(phase), speed * std::cos(phase), 0};
        compare_states(to_cartesian(elements, mu), expected, "analytic circular state");
        const auto reconstructed = from_cartesian(expected, mu);
        near(reconstructed[0], radius, 1e-6, 2e-15, "circular semimajor axis");
        near(std::hypot(reconstructed[1], reconstructed[2]), 0, 1e-14, 0,
             "circular eccentricity");
        near(std::remainder(reconstructed[5] - phase, tau), 0, 2e-14, 0,
             "circular phase");
    }
    const auto basis = rtn_basis(Cartesian{radius, 0, 0, 0, speed, 0});
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            near(basis[i][j], i == j ? 1.0 : 0.0, 1e-15, 0, "circular RTN basis");
        }
    }
}

// Independent classical perifocal construction avoids validating two mutually
// consistent but incorrect equinoctial conversion routines against each other.
struct ClassicalCase { double a, e, inclination, node, argument, true_anomaly; };

Cartesian perifocal_oracle(const ClassicalCase& c) {
    const double cn = std::cos(c.node), sn = std::sin(c.node);
    const double cw = std::cos(c.argument), sw = std::sin(c.argument);
    const double ci = std::cos(c.inclination), si = std::sin(c.inclination);
    const Vector3 p{cn * cw - sn * sw * ci, sn * cw + cn * sw * ci, sw * si};
    const Vector3 q{-cn * sw - sn * cw * ci, -sn * sw + cn * cw * ci, cw * si};
    const double parameter = c.a * (1 - c.e * c.e);
    const double cf = std::cos(c.true_anomaly), sf = std::sin(c.true_anomaly);
    const double radius = parameter / (1 + c.e * cf);
    const double scale = std::sqrt(mu / parameter);
    Cartesian state{};
    for (std::size_t k = 0; k < 3; ++k) {
        state[k] = radius * (cf * p[k] + sf * q[k]);
        state[k + 3] = scale * (-sf * p[k] + (c.e + cf) * q[k]);
    }
    return state;
}

void test_eccentric_inclined_orbits() {
    const std::array<ClassicalCase, 5> cases{{
        {26560000, .2, .7, 1.1, 2.3, -.9},
        {42164000, .01, .05, -.8, .3, 2.5},
        {42000000, .65, 1.2, 2.2, -1.7, pi},
        {70000000, .7, 2.6, -.3, 1.9, -2.2},
        {100000000, .9, 1.57, 3.1, .1, 1e-4}
    }};
    for (const auto& c : cases) {
        const double eccentric_anomaly = 2 * std::atan2(
            std::sqrt(1 - c.e) * std::sin(c.true_anomaly / 2),
            std::sqrt(1 + c.e) * std::cos(c.true_anomaly / 2));
        const double mean_anomaly = eccentric_anomaly - c.e * std::sin(eccentric_anomaly);
        const auto elements = from_keplerian(c.a, c.e, c.inclination, c.node,
                                             c.argument, mean_anomaly);
        const auto expected = perifocal_oracle(c);
        const auto actual = to_cartesian(elements, mu);
        compare_states(actual, expected, "independent perifocal oracle");
        const auto recovered = from_cartesian(expected, mu);
        near(recovered[0], c.a, 2e-5, 3e-13, "recovered semimajor axis");
        for (std::size_t k = 1; k < 5; ++k) {
            near(recovered[k], elements[k], 2e-13, 2e-13, "recovered equinoctial element");
        }
        near(std::remainder(recovered[5] - elements[5], tau), 0, 3e-13, 0,
             "recovered mean longitude");
        compare_states(to_cartesian(recovered, mu), expected, "Cartesian round trip");

        const auto r = position(actual), v = velocity(actual), h = cross(r, v);
        const double h_expected = std::sqrt(mu * c.a * (1 - c.e * c.e));
        near(dot(v, v) / 2 - mu / norm(r), -mu / (2 * c.a), 1e-7, 2e-12,
             "specific energy");
        near(norm(h), h_expected, 1e-4, 3e-13, "angular momentum magnitude");
        const Vector3 normal{std::sin(c.inclination) * std::sin(c.node),
                             -std::sin(c.inclination) * std::cos(c.node),
                             std::cos(c.inclination)};
        for (std::size_t k = 0; k < 3; ++k) {
            near(h[k] / norm(h), normal[k], 5e-14, 0, "orbital plane normal");
        }
        near(dot(normal, r), 0, 2e-7, 0, "position lies in orbital plane");
        near(dot(normal, v), 0, 2e-10, 0, "velocity lies in orbital plane");
        const auto basis = rtn_basis(actual);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                near(dot(basis[i], basis[j]), i == j ? 1 : 0, 2e-14, 0,
                     "RTN orthonormality");
            }
        }
    }
}

void test_invalid_orbits() {
    rejects([] { to_cartesian(Elements{7000000, 1, 0, 0, 0, 0}, mu); },
            "parabolic state");
    rejects([] { to_cartesian(Elements{7000000, 0, 0, 0, 0, 0}, -mu); },
            "negative gravitational parameter");
    rejects([] { from_cartesian(Cartesian{}, mu); }, "zero Cartesian state");
    rejects([] { from_cartesian(Cartesian{7000000, 0, 0, 1, 0, 0}, mu); },
            "zero angular momentum");
    rejects([] { from_keplerian(7000000, .01, pi, 0, 0, 0); },
            "retrograde singularity");
    rejects([] { validate_elements(Elements{7000000, .2, 0, 0, 0, 0}, mu, 6500000); },
            "perigee floor");
    rejects([] { to_cartesian(Elements{7000000, 0, 0, 0, 0,
                std::numeric_limits<double>::quiet_NaN()}, mu); }, "nonfinite longitude");
}

Matrix6 product_with_transpose(const Matrix6& a) {
    Matrix6 product{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            for (std::size_t k = 0; k < 6; ++k) product[i][j] += a[i][k] * a[j][k];
        }
    }
    return product;
}

void check_factor(const Matrix6& covariance) {
    const auto factor = covariance_factor(covariance);
    const auto reconstructed = product_with_transpose(factor);
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            const double natural_scale = std::sqrt(covariance[i][i]) *
                                         std::sqrt(covariance[j][j]);
            near(reconstructed[i][j], covariance[i][j],
                 2e-12 * natural_scale + 1e-300, 0, "PSD covariance reconstruction");
            if (j > i) near(factor[i][j], 0, 0, 0, "factor is lower triangular");
        }
    }
}

void test_covariance_factor() {
    // Twelve orders of magnitude in standard deviations exercise SI scaling;
    // the final row is exactly deterministic, with a zero variance.
    const std::array<double, 6> scales{1e6, 3e3, 7e2, 2e-3, 4e-6, 0};
    Matrix6 factor{{
        {{1, 0, 0, 0, 0, 0}},
        {{.3, std::sqrt(.91), 0, 0, 0, 0}},
        {{.2, -.15, std::sqrt(.9375), 0, 0, 0}},
        {{.4, -.3, 0, std::sqrt(.75), 0, 0}},
        {{-.2, .25, 0, .1, std::sqrt(.8875), 0}},
        {{0, 0, 0, 0, 0, 0}}
    }};
    for (std::size_t i = 0; i < 6; ++i)
        for (double& value : factor[i]) value *= scales[i];
    check_factor(product_with_transpose(factor));

    Matrix6 rank_one{};
    const std::array<double, 6> direction{1e6, -.001, .5, 0, 3, 1e-6};
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) rank_one[i][j] = direction[i] * direction[j];
    check_factor(rank_one);
    check_factor(Matrix6{});
}

void test_invalid_covariances() {
    Matrix6 covariance{};
    for (std::size_t i = 0; i < 6; ++i) covariance[i][i] = 1;
    auto invalid = covariance;
    invalid[0][1] = invalid[1][0] = 1.1;
    rejects([&] { covariance_factor(invalid); }, "indefinite covariance");
    invalid = covariance;
    invalid[0][1] = .25;
    invalid[1][0] = .2;
    rejects([&] { covariance_factor(invalid); }, "asymmetric covariance");
    invalid = covariance;
    invalid[2][2] = -1e-20;
    rejects([&] { covariance_factor(invalid); }, "negative variance");
    invalid = covariance;
    invalid[2][2] = 0;
    invalid[2][3] = invalid[3][2] = 1e-30;
    rejects([&] { covariance_factor(invalid); }, "nonzero covariance on deterministic axis");
    invalid = covariance;
    invalid[2][3] = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { covariance_factor(invalid); }, "nonfinite covariance");
}

void test_circular_statistics() {
    std::vector<double> two_lobes;
    for (std::size_t i = 0; i < 100; ++i) {
        two_lobes.push_back(0);
        two_lobes.push_back(pi);
    }
    const auto lobes = phase_metrics(two_lobes, 72);
    near(lobes.resultants[0], 0, 1e-14, 0, "opposing lobes R1");
    near(lobes.resultants[1], 1, 1e-14, 0, "opposing lobes R2 detects nonuniformity");
    near(lobes.resultants[2], 0, 1e-14, 0, "opposing lobes R3");
    near(lobes.resultants[3], 1, 1e-14, 0, "opposing lobes R4");
    near(lobes.max_gap, pi, 1e-14, 0, "opposing lobes angular gap");

    const auto across_zero = phase_metrics({-.1, .1}, 72);
    near(across_zero.max_gap, tau - .2, 2e-14, 0, "cluster straddling zero");
    const auto closure_gap = phase_metrics({.4, .6, 1.2, 2.5}, 72);
    near(closure_gap.max_gap, tau + .4 - 2.5, 2e-14, 0, "largest gap crosses branch cut");

    const auto singleton = phase_metrics({-7 * tau + .3}, 72);
    near(singleton.max_gap, tau, 1e-14, 0, "single-particle gap");
    near(singleton.sigma, 0, 0, 0, "single-particle spread");
    near(singleton.total_variation, 1 - 1.0 / 72, 2e-14, 0, "single-particle TV");
}

void test_uniform_phase_and_unwrapped_width() {
    std::vector<double> uniform;
    constexpr std::size_t bins = 72, points_per_bin = 5;
    constexpr std::size_t count = bins * points_per_bin;
    for (std::size_t i = 0; i < count; ++i)
        uniform.push_back(tau * (static_cast<double>(i) + .5) / count);
    const auto metrics = phase_metrics(uniform, bins);
    near(metrics.occupied_fraction, 1, 0, 0, "uniform occupancy");
    near(metrics.total_variation, 0, 1e-15, 0, "uniform total variation");
    near(metrics.max_gap, tau / count, 2e-14, 0, "uniform largest gap");
    for (double resultant : metrics.resultants)
        near(resultant, 0, 2e-14, 0, "uniform circular harmonic");
    for (auto occupancy : metrics.histogram)
        require(occupancy == points_per_bin, "uniform phase bin count");

    std::vector<double> several_turns;
    for (std::size_t i = 0; i <= 100; ++i)
        several_turns.push_back(-2 * tau + 4 * tau * static_cast<double>(i) / 100);
    const auto wrapped = phase_metrics(several_turns, bins);
    near(wrapped.q95_width, 3.8 * tau, 2e-13, 0, "unwrapped 95-percent width");
    require(wrapped.q95_width > tau, "unwrapped width retains multiple revolutions");
    const double expected_sigma = (4 * tau / 100) * std::sqrt(101.0 * 102.0 / 12.0);
    near(wrapped.sigma, expected_sigma, 2e-13, 0, "sample standard deviation convention");
}

void test_same_energy_has_no_keplerian_phase_shear() {
    constexpr double a = 26560000.0;
    const double n = std::sqrt(mu / (a * a * a));
    std::vector<double> initial_phase, final_phase;
    const double elapsed = 5 * 365.25 * day;
    for (std::size_t i = 0; i < 21; ++i) {
        const double offset = (static_cast<double>(i) - 10) * .001;
        auto elements = from_keplerian(a, .01 + .005 * i, .4 + .01 * i, 0, 0, 1 + offset);
        initial_phase.push_back(elements[5]);
        // Equal semimajor axes, not identical orbit shapes. The exact two-body
        // solution advances every member by the same mean phase n*t.
        elements[5] += n * elapsed;
        final_phase.push_back(elements[5]);
        const auto state = to_cartesian(elements, mu);
        const double energy = dot(velocity(state), velocity(state)) / 2 - mu / norm(position(state));
        near(energy, -mu / (2 * a), 1e-6, 3e-13, "equal-energy ensemble member");
    }
    const auto before = phase_metrics(initial_phase, 72);
    const auto after = phase_metrics(final_phase, 72);
    near(after.q95_width, before.q95_width, 2e-11, 0, "equal-energy width unchanged");
    near(after.sigma, before.sigma, 2e-11, 0, "equal-energy spread unchanged");
    near(after.max_gap, before.max_gap, 2e-11, 0, "equal-energy gap unchanged");
    for (std::size_t k = 0; k < 4; ++k)
        near(after.resultants[k], before.resultants[k], 2e-11, 0,
             "equal-energy harmonics unchanged");
}

void test_onset_intervals() {
    const std::vector<double> times{0, 10, 20, 30, 40, 50};
    const auto event = first_sustained(times, {false, true, false, true, true, true}, 2);
    require(event.found, "sustained event detected");
    near(event.lower_s, 20, 0, 0, "onset lower bracket precedes first true sample");
    near(event.upper_s, 30, 0, 0, "onset upper bracket is first true sample");
    require(event.confirmed_at == 4, "confirmation index follows persistence requirement");
    const auto initial = first_sustained(times, {true, true, true, false, false, false}, 3);
    require(initial.found && initial.confirmed_at == 2, "initial sustained event");
    near(initial.lower_s, 0, 0, 0, "initial lower bracket");
    near(initial.upper_s, 0, 0, 0, "initial upper bracket");
    require(!first_sustained(times, {false, true, false, true, false, true}, 2).found,
            "isolated crossings do not satisfy persistence");
    require(!first_sustained(times, {false, false, false, false, true, true}, 3).found,
            "incomplete final run does not satisfy persistence");
    const auto single = first_sustained(times, {false, true, false, false, false, false}, 1);
    require(single.found && single.confirmed_at == 1, "persistence one");
    near(single.lower_s, 0, 0, 0, "one-sample lower bracket");
    near(single.upper_s, 10, 0, 0, "one-sample upper bracket");
    rejects([&] { first_sustained(times, {true}, 1); }, "mismatched event lengths");
    rejects([&] { first_sustained(times, {true, true, true, true, true, true}, 0); },
            "zero persistence");
}

void test_dkw_and_invalid_statistics() {
    require(dkw_epsilon(738) < .05 && dkw_epsilon(737) > .05,
            "95-percent DKW count for five-percent CDF error");
    require(dkw_epsilon(4612) < .02 && dkw_epsilon(4611) > .02,
            "95-percent DKW count for two-percent CDF error");
    require(dkw_epsilon(18445) < .01 && dkw_epsilon(18444) > .01,
            "95-percent DKW count for one-percent CDF error");
    require(dkw_epsilon(20000, .01) > dkw_epsilon(20000, .05),
            "higher confidence enlarges DKW bound");
    near(sample_sigma({1, 2, 3}), 1, 1e-15, 0, "unbiased sample variance");
    rejects([] { dkw_epsilon(0); }, "zero DKW sample count");
    rejects([] { dkw_epsilon(10, 1); }, "invalid DKW confidence");
    rejects([] { phase_metrics({}, 72); }, "empty phase sample");
    rejects([] { phase_metrics({0, 1}, 3); }, "too few phase bins");
    rejects([] { phase_metrics({std::numeric_limits<double>::quiet_NaN()}, 72); },
            "nonfinite singleton phase");
    rejects([] { sample_sigma({std::numeric_limits<double>::infinity()}); },
            "nonfinite singleton statistic");
    rejects([] { quantile_sorted({}, .5); }, "empty quantile input");
    rejects([] { quantile_sorted({1, 2}, 1.1); }, "invalid quantile fraction");
}
}

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"analytic circular orbit", test_circular_orbit},
        {"independent eccentric/inclined orbit oracle", test_eccentric_inclined_orbits},
        {"invalid orbit inputs", test_invalid_orbits},
        {"scaled and semidefinite covariance", test_covariance_factor},
        {"invalid covariance rejection", test_invalid_covariances},
        {"circular harmonics and branch cut", test_circular_statistics},
        {"uniform phase and unwrapped width", test_uniform_phase_and_unwrapped_width},
        {"equal energy has no Keplerian phase shear", test_same_energy_has_no_keplerian_phase_shear},
        {"sustained onset intervals", test_onset_intervals},
        {"DKW bounds and invalid statistics", test_dkw_and_invalid_statistics}
    };
    std::size_t failed = 0;
    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "PASS: " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failed;
            std::cerr << "FAIL: " << test.first << ": " << error.what() << '\n';
        }
    }
    std::cout << tests.size() - failed << '/' << tests.size() << " math test groups passed\n";
    return failed == 0 ? 0 : 1;
}
