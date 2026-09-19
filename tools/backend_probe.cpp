#include "distribution/backend.hpp"
#include "distribution/orbit.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) result.push_back(field);
    return result;
}

int compare(const std::string& actual_path, const std::string& reference_path) {
    std::ifstream actual(actual_path), reference(reference_path);
    if (!actual || !reference) throw std::runtime_error("Cannot read parity CSV");
    std::string a_line, r_line;
    std::getline(actual, a_line); std::getline(reference, r_line);
    std::size_t count = 0;
    double max_position = 0.0, max_velocity = 0.0;
    while (std::getline(reference, r_line)) {
        if (!std::getline(actual, a_line)) throw std::runtime_error("Missing actual parity row");
        const auto a = split(a_line), r = split(r_line);
        if (a.size() != 16 || r.size() != 16 || a[0] != r[0] || a[1] != r[1] || a[3] != r[3] ||
            std::stod(a[2]) != std::stod(r[2])) throw std::runtime_error("Mismatched parity row");
        double position_squared = 0.0, velocity_squared = 0.0;
        for (std::size_t i = 4; i < 16; ++i) {
            const double av = std::stod(a[i]), rv = std::stod(r[i]);
            if (!std::isfinite(av) || !std::isfinite(rv)) throw std::runtime_error("Nonfinite parity result");
            const double delta = av - rv;
            if (i >= 10 && i < 13) position_squared += delta * delta;
            if (i >= 13) velocity_squared += delta * delta;
            const double element_limit = i == 4 ? 1e-4 : 8e-10;
            if (i < 10 && std::abs(delta) > element_limit) throw std::runtime_error("Element parity limit exceeded at " + r_line);
        }
        max_position = std::max(max_position, std::sqrt(position_squared));
        max_velocity = std::max(max_velocity, std::sqrt(velocity_squared));
        ++count;
    }
    if (std::getline(actual, a_line) || count != 32) throw std::runtime_error("Parity fixture must contain 32 cases");
    // This is an independent one-year same-model numerical parity gate, not a
    // bound on real-orbit physical accuracy. Thresholds are fixed and versioned.
    if (max_position > 0.02 || max_velocity > 3e-6) throw std::runtime_error("Java Orekit Cartesian parity limit exceeded");
    std::cout << std::setprecision(10) << "Java Orekit parity: " << count << " cases, max position "
              << max_position << " m, max velocity " << max_velocity << " m/s\n";
    return 0;
}
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) { std::cerr << "Expected output CSV path and optional Java reference CSV\n"; return 2; }
  try {
    std::ofstream out(argv[1]);
    if (!out) return 2;
    out << std::setprecision(17);
    out << "model,initial_type,days,output_type,a,ex,ey,hx,hy,lm,x,y,z,vx,vy,vz\n";
    for (const std::string model : {"j2", "j2_j2sq"}) {
        for (const std::string initial_type : {"mean", "osculating"}) {
            distribution::BackendConfig config;
            config.force_model = model;
            config.initial_type = initial_type;
            distribution::Backend backend(config, {26560000.0, 0.01, -0.004, 0.4, 0.3, 0.8});
            for (double days : {0.0, 1.0, 30.0, 365.0}) {
                const auto result = backend.advance(days * 86400.0);
                for (const std::string output_type : {"mean", "osculating"}) {
                    const auto& y = output_type == "mean" ? result.mean : result.osculating;
                    out << model << ',' << initial_type << ',' << days << ',' << output_type;
                    for (double value : y) out << ',' << value;
                    for (double value : distribution::to_cartesian(y, config.mu)) out << ',' << value;
                    out << '\n';
                }
            }
        }
    }
    out.close();
    if (!out) return 2;
    return argc == 3 ? compare(argv[1], argv[2]) : 0;
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << '\n';
    return 1;
  }
}
