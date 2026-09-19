#include <cmath>
#include <distribution/simulation.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void help() {
    std::cout << R"(distribution_propagator --config FILE --output DIRECTORY [options]
  --write-example FILE      Write a commented demonstration configuration and exit
  --samples N              Override Gaussian ensemble size (empirical uses all rows)
  --seed N                 Override unsigned 64-bit seed
  --threads N              Worker count; 0 selects hardware concurrency
  --visual-samples N        Display subset limit; 0 is useful with --no-html
  --duration-days DAYS     Override propagation horizon
  --output-step-days DAYS  Override snapshot interval; final epoch always included
  --accuracy-check N       Repropagate first N particles with tighter tolerances
  --no-html                Skip the embedded offline HTML viewer
  --export-states          Write all propagated Cartesian states to CSV
  --help                   Show this help

All input states share an epoch and inertial equatorial frame. No drag or
measurement updates. Output directory must be absent or empty. Without --config,
the deliberately broad MEO demonstration defaults are used. Read docs/SCIENCE.md
before interpreting coverage as near-uniform mixing or physical predictability.
)";
}
std::uint64_t integer(const std::string &text) {
    if (text.empty() || text[0] == '-')
        throw std::invalid_argument("Expected nonnegative integer: " + text);
    std::size_t at = 0;
    const auto value = std::stoull(text, &at);
    if (at != text.size())
        throw std::invalid_argument("Invalid integer: " + text);
    return value;
}
double scalar(const std::string &text) {
    std::size_t at = 0;
    const double value = std::stod(text, &at);
    if (at != text.size() || !std::isfinite(value))
        throw std::invalid_argument("Invalid number: " + text);
    return value;
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path config_file, output = "outputs/run", example;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--help") {
                help();
                return 0;
            }
            if (a == "--config" || a == "--write-example") {
                if (++i >= argc)
                    throw std::invalid_argument("Missing option argument");
                if (a == "--config")
                    config_file = argv[i];
                else
                    example = argv[i];
            }
        }
        if (!example.empty()) {
            if (std::filesystem::exists(example))
                throw std::runtime_error("Example output already exists");
            std::ofstream out(example);
            out.exceptions(std::ios::failbit | std::ios::badbit);
            out << distribution::config_template();
            out.close();
            return 0;
        }
        auto config =
            config_file.empty() ? distribution::Config{} : distribution::read_config(config_file);
        bool overridden_samples = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--no-html") {
                config.write_html = false;
                continue;
            }
            if (arg == "--export-states") {
                config.export_states = true;
                continue;
            }
            if (++i >= argc)
                throw std::invalid_argument("Missing argument for " + arg);
            const std::string value = argv[i];
            if (arg == "--config")
                continue;
            if (arg == "--output")
                output = value;
            else if (arg == "--samples") {
                config.samples = static_cast<std::size_t>(integer(value));
                overridden_samples = true;
            } else if (arg == "--threads")
                config.threads = static_cast<std::size_t>(integer(value));
            else if (arg == "--visual-samples")
                config.visual_samples = static_cast<std::size_t>(integer(value));
            else if (arg == "--seed")
                config.seed = integer(value);
            else if (arg == "--duration-days")
                config.duration_days = scalar(value);
            else if (arg == "--output-step-days")
                config.output_step_days = scalar(value);
            else if (arg == "--accuracy-check")
                config.accuracy_check = static_cast<std::size_t>(integer(value));
            else
                throw std::invalid_argument("Unknown option: " + arg);
        }
        if (overridden_samples && !config.empirical_samples_csv.empty())
            throw std::invalid_argument("--samples cannot resample an empirical ensemble");
        if (std::filesystem::exists(output) && !std::filesystem::is_empty(output))
            throw std::runtime_error("Output directory must be empty: " + output.string());
        const auto result = distribution::run_simulation(config);
        distribution::write_results(result, output);
        std::cout << "Propagated " << result.initial.size() << " particles at "
                  << result.frames.size() << " epochs in " << result.elapsed_seconds << " s using "
                  << result.threads_used << " threads.\n";
        const auto event = [](const char *label, const distribution::EventInterval &value) {
            std::cout << label << ": ";
            if (value.found)
                std::cout << '[' << value.lower_s / distribution::day << ", "
                          << value.upper_s / distribution::day
                          << "] days (sampled onset bracket)\n";
            else
                std::cout << "not observed within horizon/persistence criterion\n";
        };
        event("Orbital phase coverage", result.coverage);
        event("Phase mixing", result.mixing);
        if (config.accuracy_check)
            std::cout << "Tighter integration check: maximum position difference "
                      << result.accuracy_max_position_m << " m; phase "
                      << result.accuracy_max_phase_rad << " rad.\n";
        std::cout << "Results: " << std::filesystem::absolute(output).string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
