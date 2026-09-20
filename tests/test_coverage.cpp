#include <distribution/coverage.hpp>
#include <distribution/simulation.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace distribution;

// Release builds must execute all checks, independently of NDEBUG.
void require(bool condition, const std::string &context) {
    if (!condition)
        throw std::runtime_error(context);
}

void near(double actual, double expected, double tolerance, const std::string &context) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance) {
        std::ostringstream message;
        message << std::setprecision(17) << context << ": actual=" << actual
                << ", expected=" << expected << ", tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}

template <class Function> void rejects(Function &&function, const std::string &context) {
    try {
        function();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error(context + ": expected std::invalid_argument");
}

void compare_decisions(const std::vector<double> &phases, std::size_t bins,
                       const std::string &context) {
    // The established routine sorts every wrapped phase and evaluates all
    // adjacent gaps; it does not use the new bin-extrema optimization.
    const auto oracle = phase_metrics(phases, bins);
    const double width = tau / static_cast<double>(bins);
    std::vector<double> gaps{0.25 * width, width, 1.8 * width, tau,
                             oracle.max_gap,
                             std::nextafter(oracle.max_gap, 0.0),
                             std::nextafter(oracle.max_gap,
                                            std::numeric_limits<double>::infinity())};
    for (double gap : gaps) {
        if (!(gap > 0 && gap <= tau))
            continue;
        for (double fraction : {0.5, 0.9, 1.0}) {
            const bool expected = oracle.max_gap <= gap &&
                                  oracle.occupied_fraction >= fraction;
            const bool actual = covers_orbit(phases, bins, gap, fraction);
            require(actual == expected, context + ": decision differs from sorted oracle, bins=" +
                                            std::to_string(bins) + ", gap=" +
                                            std::to_string(gap) + ", fraction=" +
                                            std::to_string(fraction));
        }
    }
}

void test_random_clouds_against_sorted_oracle() {
    std::mt19937_64 engine(19840430);
    std::normal_distribution<double> normal(0, 1);
    std::uniform_real_distribution<double> uniform(-tau, tau);
    for (std::size_t bins : {4u, 7u, 36u, 72u, 360u}) {
        for (std::size_t count : {1u, 2u, 17u, 512u, 5000u}) {
            for (double scale : {0.01, 1.0, 8.0}) {
                std::vector<double> phases(count);
                for (auto &phase : phases)
                    phase = 1234.0 * tau + scale * normal(engine);
                compare_decisions(phases, bins, "Gaussian cloud");
            }
            std::vector<double> phases(count);
            for (auto &phase : phases)
                phase = uniform(engine) + 123456.0 * tau;
            compare_decisions(phases, bins, "Wrapped uniform cloud");
        }
    }
}

void test_bin_edges_and_known_gaps() {
    for (std::size_t bins : {4u, 7u, 72u, 360u}) {
        const double width = tau / static_cast<double>(bins);
        for (double shift : {0.0, -8.0 * tau, 1000000.0 * tau}) {
            std::vector<double> phases;
            for (std::size_t j = 0; j <= bins; ++j) {
                const double edge = static_cast<double>(j) * width + shift;
                phases.push_back(edge);
                phases.push_back(std::nextafter(edge,
                                                -std::numeric_limits<double>::infinity()));
                phases.push_back(std::nextafter(edge,
                                                std::numeric_limits<double>::infinity()));
            }
            compare_decisions(phases, bins, "Floating-point bin boundaries");
        }
        std::vector<double> dense, alternating;
        for (std::size_t j = 0; j < bins; ++j) {
            dense.push_back((static_cast<double>(j) + 0.25) * width);
            dense.push_back((static_cast<double>(j) + 0.75) * width);
            alternating.push_back((static_cast<double>(j) + (j % 2 ? 0.99 : 0.01)) * width);
        }
        require(covers_orbit(dense, bins, width), "Dense equally spaced cloud covers orbit");
        require(!covers_orbit(alternating, bins, width),
                "Occupying every bin alone does not guarantee the maximum-gap condition");
        compare_decisions(dense, bins, "Dense deterministic cloud");
        compare_decisions(alternating, bins, "Inter-bin gap failure");
    }
    compare_decisions({0, tau, -tau, 3 * tau}, 72, "Coincident phases across revolutions");
}

void test_coverage_input_validation() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    rejects([] { covers_orbit({}, 72, 0.1); }, "Empty ensemble");
    rejects([] { covers_orbit({0}, 3, 0.1); }, "Too few bins");
    for (double gap : {0.0, -1.0, tau + 0.01, nan, inf})
        rejects([&] { covers_orbit({0}, 72, gap); }, "Invalid gap");
    for (double fraction : {0.0, -1.0, 1.01, nan, inf})
        rejects([&] { covers_orbit({0}, 72, 0.1, fraction); }, "Invalid occupied fraction");
    for (double phase : {nan, inf, -inf})
        rejects([&] { covers_orbit({0, phase}, 72, 0.1); }, "Nonfinite phase");
}

Config ensemble_config(const std::string &force, double eccentricity, double inclination) {
    Config config;
    config.backend.force_model = force;
    config.backend.output_type = "mean";
    config.nominal = from_keplerian((config.backend.earth_radius_m + 3000000) /
                                       (1 - eccentricity),
                                   eccentricity, inclination * pi / 180, 20 * pi / 180,
                                   30 * pi / 180, 10 * pi / 180);
    config.uncertainty_coordinates = "rtn";
    config.covariance = {};
    for (std::size_t k = 0; k < 6; ++k)
        config.covariance[k][k] = k < 3 ? 1e10 : 1.0;
    config.samples = 512;
    config.threads = 2;
    config.visual_samples = 0;
    config.write_html = false;
    config.minimum_altitude_m = 0;
    config.duration_days = 4.375;
    config.output_step_days = 0.125;
    return config;
}

void compare_fixed_scanner(const Config &config, const std::string &context,
                           bool expected_event) {
    CoverageOptions options;
    options.automatic_times = false;
    const auto scan = run_coverage(config, options);
    const auto oracle = run_simulation(config);
    require(scan.samples == oracle.initial.size(), context + ": all particles retained");
    require(scan.earth_intersections == 0, context + ": valid Earth domain");
    require(scan.coverage.found == oracle.coverage.found, context + ": same coverage outcome");
    require(scan.coverage.found == expected_event, context + ": intended event exercised");
    near(scan.cadence_s, config.output_step_days * day, 1e-9, context + ": configured cadence");
    if (oracle.coverage.found) {
        require(scan.status == "ok", context + ": success status");
        require(scan.coverage.confirmed_at == oracle.coverage.confirmed_at,
                context + ": same persistence confirmation epoch");
        near(scan.coverage.lower_s, oracle.coverage.lower_s, 1e-8, context + ": onset lower bound");
        near(scan.coverage.upper_s, oracle.coverage.upper_s, 1e-8, context + ": onset upper bound");
        near(scan.confirmation_s, oracle.frames.at(oracle.coverage.confirmed_at).time_s,
             1e-8, context + ": confirmation time");
        const auto onset_index = oracle.coverage.confirmed_at + 1 - config.persistence;
        near(scan.onset_max_gap_deg, oracle.frames.at(onset_index).phase.max_gap * 180 / pi,
             2e-7, context + ": reported onset gap from every particle");
        require(scan.evaluated_epochs == oracle.coverage.confirmed_at + 1,
                context + ": stops after first confirmed event");
    } else {
        require(scan.status == "unobserved", context + ": unobserved status");
        require(scan.evaluated_epochs == oracle.frames.size(),
                context + ": evaluates every distinct saved epoch");
        near(scan.horizon_s, config.duration_days * day, 1e-8, context + ": final partial epoch");
    }
}

void test_fixed_scanner_matches_full_propagation() {
    for (const std::string force : {"kepler", "j2", "j2_j2sq"}) {
        for (const auto shape : {std::array<double, 2>{0.0, 0.0},
                                 std::array<double, 2>{0.05, 55.0},
                                 std::array<double, 2>{0.1, 90.0}}) {
            auto config = ensemble_config(force, shape[0], shape[1]);
            compare_fixed_scanner(config, force + ", e=" + std::to_string(shape[0]) +
                                              ", i=" + std::to_string(shape[1]), true);
        }
    }
    auto short_run = ensemble_config("j2_j2sq", 0.05, 55);
    short_run.samples = 64;
    short_run.duration_days = 0.01;
    short_run.output_step_days = 0.003;
    compare_fixed_scanner(short_run, "No event with a nonuniform final interval", false);

    auto already_covered = ensemble_config("kepler", 0.05, 55);
    already_covered.backend.initial_type = "mean";
    already_covered.uncertainty_coordinates = "equinoctial";
    already_covered.covariance = {};
    already_covered.covariance[5][5] = 100;
    already_covered.samples = 1024;
    already_covered.duration_days = 0.3;
    already_covered.output_step_days = 0.1;
    compare_fixed_scanner(already_covered, "Initially covered cloud without phase shear", true);
    const auto automatic_initial = run_coverage(already_covered);
    require(automatic_initial.status == "ok" && automatic_initial.coverage.found,
            "Automatic mode recognizes a previously covered zero-shear cloud");
    near(automatic_initial.coverage.lower_s, 0, 0,
         "Automatic already-covered cloud has zero lower onset");
    near(automatic_initial.coverage.upper_s, 0, 0,
         "Automatic already-covered cloud has zero upper onset");
    near(automatic_initial.confirmation_s, 2 * already_covered.output_step_days * day,
         1e-8, "Zero-shear confirmation uses the configured fixed cadence");
    require(automatic_initial.evaluated_epochs == 3,
            "Initially covered cloud still satisfies three-epoch persistence");
    already_covered.persistence = 5;
    compare_fixed_scanner(already_covered, "Insufficient epochs for persistence", false);
}

void test_no_shear_and_invalid_cadence() {
    auto config = ensemble_config("kepler", 0, 0);
    config.samples = 8;
    config.covariance = {};
    const auto result = run_coverage(config);
    require(result.status == "no_phase_shear" && !result.coverage.found,
            "Coincident particles with equal rates cannot spread around the orbit");
    for (double value : {0.0, -1.0, 1e-300, 9.0, 10001.0,
                         std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
        CoverageOptions options;
        options.intervals_per_phase_sigma = value;
        rejects([&] { run_coverage(config, options); }, "Invalid automatic cadence scaling");
    }
    for (std::size_t value : {11u, 64u}) {
        CoverageOptions options;
        options.max_extensions = value;
        rejects([&] { run_coverage(config, options); }, "Extension limit avoids integer overflow");
    }
    CoverageOptions excessive;
    excessive.intervals_per_phase_sigma = 10000;
    excessive.max_extensions = 10;
    rejects([&] { run_coverage(config, excessive); }, "Bound total epoch work before sampling");
}

void test_domain_audit_retains_gaussian_tail() {
    auto config = ensemble_config("j2_j2sq", 0, 55);
    config.nominal = from_keplerian(config.backend.earth_radius_m + 1000000, 0,
                                   55 * pi / 180, 20 * pi / 180, 30 * pi / 180, 10 * pi / 180);
    config.samples = 20000;
    for (std::size_t k = 3; k < 6; ++k)
        config.covariance[k][k] = 0.0001;
    const auto samples = sample_initial_for_domain_audit(config);
    require(samples.size() == 20000, "Domain audit retains every Gaussian draw");
    std::size_t earth = 0, below300 = 0, below500 = 0, below1000 = 0;
    double minimum = std::numeric_limits<double>::infinity();
    for (const auto &e : samples) {
        const double altitude = (e[0] * (1 - std::hypot(e[1], e[2])) -
                                 config.backend.earth_radius_m) / 1000;
        minimum = std::min(minimum, altitude);
        earth += altitude < 0;
        below300 += altitude < 300;
        below500 += altitude < 500;
        below1000 += altitude < 1000;
    }
    require(earth == 4 && below300 == 174 && below500 == 1037 && below1000 == 13247,
            "Known seeded domain census, independently verified with Cartesian energy/angular momentum");
    near(minimum, -63.0966, 0.001, "Known minimum initial perigee in km");
    const auto result = run_coverage(config);
    require(result.status == "earth_intersection" && result.samples == 20000,
            "Physically invalid ensemble is masked without replacement draws");
    require(result.earth_intersections == earth && result.particles_below_300km == below300 &&
                result.particles_below_500km == below500 &&
                result.particles_below_1000km == below1000,
            "Coverage result preserves the complete initial-domain census");
    near(result.min_initial_perigee_km, minimum, 1e-10, "Domain minimum propagated to result");
    require(result.evaluated_epochs == 0 && !result.coverage.found,
            "Earth-crossing ensemble produces no simulated coverage value");
    bool rejected = false;
    try {
        static_cast<void>(sample_initial(config));
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("samples are not clipped or redrawn") !=
                   std::string::npos;
    }
    require(rejected, "Ordinary sampler still rejects the original Earth-crossing draw");
}
} // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"random clouds match sorted coverage", test_random_clouds_against_sorted_oracle},
        {"bin edges and explicit gap failures", test_bin_edges_and_known_gaps},
        {"coverage input validation", test_coverage_input_validation},
        {"fixed scanner matches full propagation", test_fixed_scanner_matches_full_propagation},
        {"no phase shear and invalid cadence", test_no_shear_and_invalid_cadence},
        {"Gaussian domain census without redraw", test_domain_audit_retains_gaussian_tail}};
    try {
        for (const auto &test : tests) {
            test.second();
            std::cout << "PASS " << test.first << '\n';
        }
    } catch (const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "All " << tests.size() << " coverage test groups passed.\n";
}
