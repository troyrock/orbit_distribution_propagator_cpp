#include "viewer_template.hpp"
#include <cmath>
#include <distribution/simulation.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace distribution {
namespace {
std::string quote(const std::string &value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        switch (c) {
        case '"':
            out << "\\\"";
            break;
        case '\\':
            out << "\\\\";
            break;
        case '\n':
            out << "\\n";
            break;
        case '\r':
            out << "\\r";
            break;
        case '\t':
            out << "\\t";
            break;
        default:
            if (c < 32 || c == '<')
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<unsigned int>(c) << std::dec;
            else
                out << c;
        }
    }
    return out.str() + '"';
}
template <class T> void array(std::ostream &out, const T &values) {
    out << '[';
    bool first = true;
    for (const auto &x : values) {
        if (!first)
            out << ',';
        first = false;
        out << x;
    }
    out << ']';
}
template <class T> void nested_array(std::ostream &out, const T &values) {
    out << '[';
    bool first = true;
    for (const auto &row : values) {
        if (!first)
            out << ',';
        first = false;
        array(out, row);
    }
    out << ']';
}
void nullable(std::ostream &out, double x) {
    if (std::isfinite(x) && x >= 0)
        out << x;
    else
        out << "null";
}
void event_time(std::ostream &out, const EventInterval &e) {
    if (e.found)
        out << e.upper_s;
    else
        out << "null";
}
void event_bracket(std::ostream &out, const EventInterval &e) {
    if (e.found)
        out << '[' << e.lower_s << ',' << e.upper_s << ']';
    else
        out << "null";
}
void summary(std::ostream &out, const Simulation &s) {
    out << "{\"sample_count\":" << s.initial.size() << ",\"frame_count\":" << s.frames.size()
        << ",\"threads_used\":" << s.threads_used << ",\"elapsed_seconds\":" << s.elapsed_seconds
        << ",\"coverage_time_s\":";
    event_time(out, s.coverage);
    out << ",\"coverage_bracket_s\":";
    event_bracket(out, s.coverage);
    out << ",\"mixing_time_s\":";
    event_time(out, s.mixing);
    out << ",\"mixing_bracket_s\":";
    event_bracket(out, s.mixing);
    out << ",\"central95_wrap_time_s\":";
    event_time(out, s.central95_wrap);
    out << ",\"central95_wrap_bracket_s\":";
    event_bracket(out, s.central95_wrap);
    out << ",\"analytic_mixing_time_s\":";
    nullable(out, s.analytic_mixing_time_s);
    out << ",\"analytic_estimate_model\":\"Gaussian Kepler phase shear using initial sampled mean "
           "elements; not the full numerical mixing criterion\""
        << ",\"dkw_cdf_error_95\":" << dkw_epsilon(s.initial.size()) << ",\"dkw_applicability\":"
        << quote(s.config.empirical_samples_csv.empty()
                     ? "IID initial Gaussian sampling; pointwise-in-time CDF bound"
                     : "Only if empirical particles are IID draws; not guaranteed for arbitrary "
                       "supplied ensembles")
        << ",\"accuracy_checked_samples\":" << std::min(s.config.accuracy_check, s.initial.size())
        << ",\"accuracy_max_position_m\":";
    nullable(out, s.config.accuracy_check ? s.accuracy_max_position_m : -1.0);
    out << ",\"accuracy_max_phase_rad\":";
    nullable(out, s.config.accuracy_check ? s.accuracy_max_phase_rad : -1.0);
    out << ",\"accepted_steps\":" << s.stats.accepted_steps
        << ",\"rejected_steps\":" << s.stats.rejected_steps
        << ",\"derivative_evaluations\":" << s.stats.derivative_evaluations << '}';
}
void json(std::ostream &out, const Simulation &s) {
    const auto &c = s.config;
    const auto &b = c.backend;
    out << std::setprecision(17)
        << "{\"schema_version\":1,\"metadata\":{\"samples\":" << s.initial.size()
        << ",\"visual_samples\":" << std::min(c.visual_samples, s.initial.size())
        << ",\"force_model\":" << quote(b.force_model)
        << ",\"input_type\":" << quote(b.initial_type)
        << ",\"output_type\":" << quote(b.output_type)
        << ",\"uncertainty_coordinates\":" << quote(c.uncertainty_coordinates) << ",\"sampling\":"
        << quote(c.empirical_samples_csv.empty() ? "IID Gaussian"
                                                 : "equal-weight empirical ensemble")
        << ",\"seed\":" << c.seed << ",\"seed_string\":" << quote(std::to_string(c.seed))
        << ",\"mu\":" << b.mu << ",\"earth_radius_m\":" << b.earth_radius_m << ",\"j2\":" << b.j2
        << ",\"elapsed_seconds\":" << s.elapsed_seconds
        << ",\"dsst_revision\":" << quote(dsst_revision)
        << ",\"phase_definition\":\"Continuous mean longitude minus nominal mean longitude; "
           "histogram modulo 2pi. Not true anomaly.\""
        << ",\"phase_bins\":" << c.phase_bins << ",\"coverage_definition\":"
        << quote("Largest empty phase gap <= " + std::to_string(c.coverage_max_gap_deg) +
                 " degrees and occupied-bin fraction >= " +
                 std::to_string(c.coverage_occupied_fraction) + "; sustained for " +
                 std::to_string(c.persistence) + " snapshots")
        << ",\"mixing_definition\":"
        << quote("Coverage plus all R1..R4 <= " + std::to_string(c.mixing_max_resultant) +
                 " and histogram total variation <= " + std::to_string(c.mixing_max_tv) +
                 "; sustained for " + std::to_string(c.persistence) + " snapshots")
        << ",\"model_scope\":\"Earth monopole/J2/J2-squared as selected; no drag, Sun/Moon, SRP, "
           "higher zonals, tesseral gravity, maneuvers, measurement updates or process noise. "
           "Inertial equatorial axes, relative epoch.\""
        << ",\"nominal_elements\":";
    array(out, c.nominal);
    out << ",\"covariance\":";
    if (c.empirical_samples_csv.empty())
        nested_array(out, c.covariance);
    else
        out << "null"; // Supplied particles have no configured Gaussian covariance.
    out << ",\"settings\":{\"duration_days\":" << c.duration_days
        << ",\"output_step_days\":" << c.output_step_days
        << ",\"minimum_altitude_m\":" << c.minimum_altitude_m
        << ",\"relative_tolerance\":" << b.relative_tolerance
        << ",\"absolute_tolerance_m\":" << b.absolute_tolerance_m
        << ",\"absolute_tolerance_elements\":" << b.absolute_tolerance_elements
        << ",\"min_step_s\":" << b.min_step_s << ",\"max_step_s\":" << b.max_step_s
        << ",\"persistence\":" << c.persistence
        << ",\"coverage_max_gap_deg\":" << c.coverage_max_gap_deg
        << ",\"coverage_occupied_fraction\":" << c.coverage_occupied_fraction
        << ",\"mixing_max_resultant\":" << c.mixing_max_resultant
        << ",\"mixing_max_tv\":" << c.mixing_max_tv << ",\"max_memory_mb\":" << c.max_memory_mb
        << "}},\"summary\":";
    summary(out, s);
    out << ",\"frames\":[";
    bool first = true;
    for (const auto &f : s.frames) {
        if (!first)
            out << ',';
        first = false;
        out << "{\"time_s\":" << f.time_s << ",\"positions_m\":";
        nested_array(out, f.positions);
        out << ",\"reference_orbit_m\":";
        nested_array(out, f.reference_orbit);
        out << ",\"metrics\":{\"phase_sigma_rad\":" << f.phase.sigma
            << ",\"phase_q95_width_rad\":" << f.phase.q95_width << ",\"resultants\":";
        array(out, f.phase.resultants);
        out << ",\"max_gap_deg\":" << f.phase.max_gap * 180 / pi
            << ",\"occupied_fraction\":" << f.phase.occupied_fraction
            << ",\"total_variation\":" << f.phase.total_variation << ",\"histogram\":";
        array(out, f.phase.histogram);
        out << ",\"rtn_sigma_m\":";
        array(out, f.rtn_sigma_m);
        out << ",\"tube_rtn_sigma_m\":";
        array(out, f.tube_rtn_sigma_m);
        out << ",\"semimajor_sigma_m\":" << f.semimajor_sigma_m
            << ",\"coverage\":" << (f.coverage ? "true" : "false")
            << ",\"mixed\":" << (f.mixed ? "true" : "false") << "}}";
    }
    out << "]}\n";
}
void write_text(const std::filesystem::path &path, const std::string &contents) {
    std::ofstream out(path, std::ios::binary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << contents;
    out.close();
}
} // namespace
void write_results(const Simulation &s, const std::filesystem::path &output) {
    if (std::filesystem::exists(output) && !std::filesystem::is_empty(output))
        throw std::runtime_error("Output directory must be empty: " + output.string());
    std::filesystem::create_directories(output);
    std::ostringstream document;
    json(document, s);
    const auto data = document.str();
    write_text(output / "run.json", data);
    if (s.config.write_html) {
        std::string html = viewer_template;
        const std::string token = "__DISTRIBUTION_DATA__";
        const auto at = html.find(token);
        if (at == std::string::npos || html.find(token, at + token.size()) != std::string::npos)
            throw std::runtime_error("Invalid embedded viewer template");
        html.replace(at, token.size(), data);
        write_text(output / "visualization.html", html);
    }
    std::ostringstream initial;
    initial << std::setprecision(17) << "particle,a_m,ex,ey,hx,hy,mean_longitude_rad\n";
    for (std::size_t id = 0; id < s.initial.size(); ++id) {
        initial << id;
        for (double x : s.initial[id])
            initial << ',' << x;
        initial << '\n';
    }
    write_text(output / "initial_samples.csv", initial.str());
    std::ostringstream metrics;
    metrics << std::setprecision(17)
            << "time_s,time_days,phase_sigma_rad,phase_q95_width_rad,R1,R2,R3,R4,max_gap_deg,"
               "occupied_fraction,total_variation,r_sigma_m,t_sigma_m,n_sigma_m,tube_r_sigma_m,"
               "tube_t_sigma_m,tube_n_sigma_m,semimajor_sigma_m,coverage,mixed\n";
    for (const auto &f : s.frames) {
        metrics << f.time_s << ',' << f.time_s / day << ',' << f.phase.sigma << ','
                << f.phase.q95_width;
        for (double x : f.phase.resultants)
            metrics << ',' << x;
        metrics << ',' << f.phase.max_gap * 180 / pi << ',' << f.phase.occupied_fraction << ','
                << f.phase.total_variation;
        for (double x : f.rtn_sigma_m)
            metrics << ',' << x;
        for (double x : f.tube_rtn_sigma_m)
            metrics << ',' << x;
        metrics << ',' << f.semimajor_sigma_m << ',' << f.coverage << ',' << f.mixed << '\n';
    }
    write_text(output / "metrics.csv", metrics.str());
    if (s.config.export_states) {
        std::ofstream out(output / "states.csv");
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out << std::setprecision(17)
            << "time_s,particle,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,mean_longitude_rad,mean_a_m\n";
        for (std::size_t f = 0; f < s.frames.size(); ++f)
            for (std::size_t id = 0; id < s.initial.size(); ++id) {
                const auto &state = s.states[f * s.initial.size() + id];
                const auto cart = to_cartesian(state.output, s.config.backend.mu);
                out << s.frames[f].time_s << ',' << id;
                for (double x : cart)
                    out << ',' << x;
                out << ',' << state.mean_longitude << ',' << state.mean_a << '\n';
            }
        out.close();
    }
    std::ostringstream report;
    report << std::setprecision(17);
    summary(report, s);
    report << '\n';
    write_text(output / "summary.json", report.str());
}
} // namespace distribution
