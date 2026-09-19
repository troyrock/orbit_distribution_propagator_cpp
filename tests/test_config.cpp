#include <distribution/config.hpp>
#include <distribution/statistics.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace distribution;

void require(bool condition, const std::string &description) {
    if (!condition)
        throw std::runtime_error(description);
}
void near(double actual, double expected, double tolerance, const std::string &description) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::ostringstream message;
        message << std::setprecision(17) << description << ": actual=" << actual
                << ", expected=" << expected << ", tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}
template <class Function> void rejects(Function &&function, const std::string &description) {
    try {
        function();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error(description + ": expected std::invalid_argument");
}

// Each executable invocation owns an isolated directory and removes only the
// exact files it created. No recursive cleanup or shared-build mutations.
class Fixtures {
    std::filesystem::path directory_;
    std::vector<std::filesystem::path> files_;

  public:
    Fixtures() {
        const auto unique = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        for (unsigned attempt = 0; attempt < 100; ++attempt) {
            directory_ = std::filesystem::temp_directory_path() /
                         ("distribution-config-test-" + std::to_string(unique) + "-" +
                          std::to_string(attempt));
            if (std::filesystem::create_directory(directory_))
                return;
        }
        throw std::runtime_error("Unable to create isolated test fixture directory");
    }
    ~Fixtures() {
        std::error_code ignored;
        for (const auto &file : files_)
            std::filesystem::remove(file, ignored);
        std::filesystem::remove(directory_, ignored);
    }
    std::filesystem::path write(const std::string &name, const std::string &text) {
        const auto path = directory_ / name;
        files_.push_back(path);
        std::ofstream output(path);
        if (!(output << text))
            throw std::runtime_error("Unable to write test fixture");
        return path;
    }
};

std::string csv(const std::vector<Elements> &rows) {
    std::ostringstream out;
    out << std::setprecision(17);
    for (const auto &row : rows) {
        for (std::size_t k = 0; k < 6; ++k)
            out << (k ? "," : "") << row[k];
        out << '\n';
    }
    return out.str();
}

Config simple_config() {
    Config c;
    c.backend.force_model = "kepler";
    c.backend.initial_type = "mean";
    c.backend.output_type = "mean";
    c.covariance = {};
    c.samples = 2000;
    return c;
}

// Statistical tests compare against independently specified first/second
// moments. Bounds are six standard errors, with a deterministic fixed seed.
void check_moments(const std::vector<Elements> &values, const Elements &center,
                   const Matrix6 &covariance, const std::string &description) {
    const std::size_t count = values.size();
    std::array<long double, 6> mean{};
    std::array<std::array<long double, 6>, 6> second{};
    std::array<double, 6> scale{};
    for (std::size_t k = 0; k < 6; ++k) {
        scale[k] = std::sqrt(covariance[k][k]);
        require(scale[k] > 0, "moment test expects six stochastic axes");
    }
    for (const auto &row : values) {
        std::array<long double, 6> z{};
        for (std::size_t k = 0; k < 6; ++k) {
            z[k] = (static_cast<long double>(row[k]) - center[k]) / scale[k];
            mean[k] += z[k];
        }
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                second[i][j] += z[i] * z[j];
    }
    for (auto &value : mean)
        value /= static_cast<long double>(count);
    const double mean_tolerance = 6 / std::sqrt(static_cast<double>(count));
    for (std::size_t i = 0; i < 6; ++i) {
        near(static_cast<double>(mean[i]), 0, mean_tolerance,
             description + " standardized mean " + std::to_string(i));
        for (std::size_t j = 0; j < 6; ++j) {
            const double expected = covariance[i][j] / scale[i] / scale[j];
            const long double empirical = (second[i][j] - count * mean[i] * mean[j]) / (count - 1);
            const double tolerance = 6 * std::sqrt((1 + expected * expected) / (count - 1));
            near(static_cast<double>(empirical), expected, tolerance,
                 description + " standardized covariance " + std::to_string(i) + "," +
                     std::to_string(j));
        }
    }
}

void test_gaussian_moments_and_reproducibility() {
    auto c = simple_config();
    c.samples = 50000;
    const Elements sigma{50000, .001, .002, .0002, .0003, .005};
    for (std::size_t k = 0; k < 6; ++k)
        c.covariance[k][k] = sigma[k] * sigma[k];
    c.covariance[0][5] = c.covariance[5][0] = -.65 * sigma[0] * sigma[5];
    c.covariance[1][2] = c.covariance[2][1] = .35 * sigma[1] * sigma[2];
    const auto sample = sample_initial(c);
    check_moments(sample, c.nominal, c.covariance, "correlated equinoctial Gaussian");

    c.samples = 40;
    const auto repeated = sample_initial(c);
    for (std::size_t i = 0; i < repeated.size(); ++i)
        require(repeated[i] == sample[i], "same seed preserves exact sample prefix");
    c.threads = 7;
    require(sample_initial(c) == repeated, "sampling does not depend on requested worker count");
    ++c.seed;
    require(sample_initial(c) != repeated, "different seed changes ensemble");
}

void test_cartesian_moments() {
    auto c = simple_config();
    c.samples = 30000;
    c.uncertainty_coordinates = "cartesian";
    const Elements sigma{50, 30, 20, .03, .02, .01};
    for (std::size_t k = 0; k < 6; ++k)
        c.covariance[k][k] = sigma[k] * sigma[k];
    c.covariance[0][4] = c.covariance[4][0] = .4 * sigma[0] * sigma[4];
    const auto nominal = to_cartesian(c.nominal, c.backend.mu);
    auto samples = sample_initial(c);
    for (auto &row : samples)
        row = to_cartesian(row, c.backend.mu);
    check_moments(samples, nominal, c.covariance, "Cartesian Gaussian after element conversion");
}

void test_rank_deficient_energy_distribution() {
    auto c = simple_config();
    c.covariance[1][1] = 1e-6;
    c.covariance[3][3] = 1e-8;
    c.covariance[5][5] = 1e-6;
    const auto sample = sample_initial(c);
    const double rate = std::sqrt(c.backend.mu / c.nominal[0]) / c.nominal[0];
    std::vector<double> before, after;
    for (const auto &row : sample) {
        require(row[0] == c.nominal[0], "zero semimajor variance remains exactly deterministic");
        require(row[2] == c.nominal[2] && row[4] == c.nominal[4],
                "zero-variance axes are not jittered");
        before.push_back(row[5]);
        after.push_back(row[5] + rate * 100 * day);
    }
    require(sample.front()[1] != sample.back()[1], "nonzero eccentricity uncertainty is sampled");
    near(sample_sigma(after), sample_sigma(before), 2e-12,
         "same-energy sampled distribution has no artificial Keplerian phase diffusion");
}

void test_empirical_rtn_mapping() {
    Fixtures files;
    auto c = simple_config();
    const double a = 26560000;
    c.nominal = from_keplerian(a, 0, pi / 2, 0, 0, 0);
    c.uncertainty_coordinates = "rtn";
    const std::vector<Elements> offsets{{10, 20, 30, .1, .2, .3}, {-40, 50, -60, -.4, .5, -.6}};
    c.empirical_samples_csv = files.write("rtn.csv", csv(offsets));
    const auto samples = sample_initial(c);
    require(samples.size() == 2, "empirical file sets the ensemble count");
    const double speed = std::sqrt(c.backend.mu / a);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        // This polar circular state has R=+x, T=+z, N=-y analytically.
        // The velocity offsets are vector components in that basis, not time
        // derivatives of a rotating displacement; no omega-cross-r term.
        const auto &d = offsets[i];
        const Cartesian expected{a + d[0], -d[2], d[1], d[3], -d[5], speed + d[4]};
        const auto actual = to_cartesian(samples[i], c.backend.mu);
        for (std::size_t k = 0; k < 6; ++k)
            near(actual[k], expected[k], k < 3 ? 3e-7 : 3e-10, "empirical fixed-basis RTN mapping");
    }
}

void test_empirical_cartesian_phase_branch() {
    Fixtures files;
    auto c = simple_config();
    const double a = 26560000;
    c.nominal = Elements{a, 0, 0, 0, 0, 0};
    c.uncertainty_coordinates = "cartesian";
    const double speed = std::sqrt(c.backend.mu / a);
    std::vector<Elements> states;
    for (double phase : {-.01, .01})
        states.push_back({a * std::cos(phase), a * std::sin(phase), 0, -speed * std::sin(phase),
                          speed * std::cos(phase), 0});
    c.empirical_samples_csv = files.write("branch.csv", csv(states));
    const auto samples = sample_initial(c);
    near(samples[0][5], -.01, 2e-14, "Cartesian empirical negative local phase");
    near(samples[1][5], .01, 2e-14, "Cartesian empirical positive local phase");
    near(phase_metrics({samples[0][5], samples[1][5]}, 72).q95_width, .019, 3e-14,
         "branch crossing must not masquerade as a wrapped initial distribution");
}

void test_empirical_absolute_elements() {
    Fixtures files;
    auto c = simple_config();
    const std::vector<Elements> absolute{{26560000, .01, .02, .1, .2, 7.1},
                                         {26561000, .011, .021, .11, .21, 7.2}};
    c.empirical_samples_csv =
        files.write("elements.csv", "# Absolute equinoctial states\n" + csv(absolute));
    require(sample_initial(c) == absolute,
            "equinoctial empirical rows are absolute and preserve turns");
}

void test_config_parser_and_covariance_path() {
    Fixtures files;
    const Matrix6 matrix{{{{100, 0, 0, 0, 0, .02}},
                          {{0, 1e-8, 0, 0, 0, 0}},
                          {{0, 0, 4e-8, 0, 0, 0}},
                          {{0, 0, 0, 0, 0, 0}},
                          {{0, 0, 0, 0, 0, 0}},
                          {{.02, 0, 0, 0, 0, 1e-4}}}};
    files.write("covariance.csv", csv(std::vector<Elements>(matrix.begin(), matrix.end())));
    const auto input =
        files.write("valid.cfg", "force_model = kepler\ninitial_type=mean\noutput_type=mean\n"
                                 "orbit_keplerian_deg=26560000,.02,55,20,30,10\n"
                                 "samples=20\nseed=18446744073709551615\n"
                                 "covariance_csv=covariance.csv\n");
    const auto c = read_config(input);
    require(c.covariance == matrix, "relative covariance path resolves beside config");
    require(c.samples == 20 && c.seed == 18446744073709551615ULL,
            "integer fields preserve full range");
    near(c.nominal[5], pi / 3, 1e-15, "degree-valued mean longitude conversion");
    near(std::hypot(c.nominal[3], c.nominal[4]), std::tan(55 * pi / 360), 1e-15,
         "degree-valued inclination conversion");
    require(sample_initial(c).size() == 20, "parsed covariance is usable");

    const std::vector<std::string> invalid{
        "unrecognized=1\n",
        "samples=10\nsamples=20\n",
        "samples=-1\n",
        "samples=2.5\n",
        "sigma=1,2,3,4,5\n",
        "sigma=1,2,-3,4,5,6\n",
        "uncertainty_coordinates=cartesian\n",
        "duration_days=1oops\n",
        "sigma=0,0,0,0,0,0\ncovariance_csv=covariance.csv\n",
        "orbit_equinoctial=26560000,0,0,0,0,0\norbit_keplerian_deg=26560000,0,0,0,0,0\n"};
    for (std::size_t i = 0; i < invalid.size(); ++i) {
        const auto path = files.write("invalid-" + std::to_string(i) + ".cfg", invalid[i]);
        rejects([&] { read_config(path); }, "invalid config case " + std::to_string(i));
    }
}

void test_output_schedule_and_validation() {
    auto c = simple_config();
    c.duration_days = 2.5;
    c.output_step_days = 1;
    require(output_times(c) == std::vector<double>({0, day, 2 * day, 2.5 * day}),
            "schedule includes exact end when cadence does not divide duration");
    c.duration_days = 2;
    require(output_times(c) == std::vector<double>({0, day, 2 * day}),
            "exact end is not duplicated");
    c.duration_days = .25;
    require(output_times(c) == std::vector<double>({0, .25 * day}),
            "long output cadence still includes both endpoints");
    c.output_step_days = 0;
    rejects([&] { validate_config(c); }, "zero cadence rejected before allocation");
    c.output_step_days = 1e-12;
    rejects([&] { validate_config(c); }, "excessive schedule rejected before allocation");
    c.output_step_days = 1;
    c.max_memory_mb = 0;
    rejects([&] { validate_config(c); }, "nonpositive memory budget");
}

template <class Function> void rejects_memory(Function &&function, const std::string &description) {
    try {
        function();
    } catch (const std::exception &error) {
        const std::string message = error.what();
        require(message.find("memory") != std::string::npos ||
                    message.find("Memory") != std::string::npos,
                description + ": expected an explicit memory-budget diagnostic");
        return;
    }
    throw std::runtime_error(description + ": expected budget rejection");
}

void test_memory_preflight() {
    auto c = simple_config();
    c.samples = 10000;
    c.max_memory_mb = .001;
    rejects_memory([&] { sample_initial(c); },
                   "Gaussian preflight checks budget before allocation");

    c.samples = 2;
    c.visual_samples = 0;
    c.write_html = false;
    c.phase_bins = 100000;
    c.duration_days = 1000;
    c.output_step_days = .1;
    c.max_memory_mb = 50;
    // Only ~1.3 MB of particle states, but about 8 GB of histogram bins.
    // A states-only estimate would incorrectly allow this configuration.
    rejects_memory([&] { sample_initial(c); }, "frame histogram storage is budgeted");

    Fixtures files;
    c = simple_config();
    c.empirical_samples_csv = files.write("budget.csv", csv({c.nominal, c.nominal}));
    c.max_memory_mb = .001;
    rejects_memory([&] { sample_initial(c); }, "empirical loading enforces memory budget");
}
} // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"Gaussian moments and deterministic sample prefix",
         test_gaussian_moments_and_reproducibility},
        {"Cartesian sampling moments", test_cartesian_moments},
        {"rank-deficient energy distribution", test_rank_deficient_energy_distribution},
        {"empirical RTN mapping", test_empirical_rtn_mapping},
        {"empirical Cartesian phase branch", test_empirical_cartesian_phase_branch},
        {"absolute empirical equinoctial samples", test_empirical_absolute_elements},
        {"config parser and covariance path", test_config_parser_and_covariance_path},
        {"output schedule and validation", test_output_schedule_and_validation},
        {"memory budget preflight", test_memory_preflight}};
    std::size_t failed = 0;
    for (const auto &test : tests) {
        try {
            test.second();
            std::cout << "PASS: " << test.first << '\n';
        } catch (const std::exception &error) {
            ++failed;
            std::cerr << "FAIL: " << test.first << ": " << error.what() << '\n';
        }
    }
    std::cout << tests.size() - failed << '/' << tests.size()
              << " configuration test groups passed\n";
    return failed ? 1 : 0;
}
