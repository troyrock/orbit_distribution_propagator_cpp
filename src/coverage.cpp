#include <distribution/coverage.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace distribution {
namespace {
struct Bins {
    std::vector<double> low, high;
    std::vector<std::size_t> count;
    explicit Bins(std::size_t size) : low(size), high(size), count(size) {}
    void clear() {
        std::fill(low.begin(), low.end(), std::numeric_limits<double>::infinity());
        std::fill(high.begin(), high.end(), -std::numeric_limits<double>::infinity());
        std::fill(count.begin(), count.end(), 0);
    }
    void add(double value) {
        const double p = wrap_angle(value);
        const auto b = std::min(count.size()-1, static_cast<std::size_t>(p/tau*count.size()));
        low[b] = std::min(low[b], p); high[b] = std::max(high[b], p); ++count[b];
    }
};
double sorted_gap(std::vector<double> phase) {
    for (auto &p : phase) p = wrap_angle(p);
    std::sort(phase.begin(), phase.end());
    double gap = phase.front()+tau-phase.back();
    for (std::size_t i=1; i<phase.size(); ++i) gap=std::max(gap,phase[i]-phase[i-1]);
    return gap;
}
// -1 requires a full-gap fallback; 0/1 are exact decisions.
int bin_decision(const Bins &b, double gap, double fraction) {
    std::size_t occupied=0;
    for (auto count:b.count) occupied += count!=0;
    if (static_cast<double>(occupied)/b.count.size() < fraction) return 0;
    if (occupied != b.count.size() || gap < tau/b.count.size()) return -1;
    for (std::size_t j=0;j<b.count.size();++j) {
        const double boundary = j+1==b.count.size() ? b.low[0]+tau-b.high[j] : b.low[j+1]-b.high[j];
        if (boundary>gap) return 0;
        // Avoid assuming ideal bin boundaries when floating-point division
        // assigns endpoints a few ulps beyond a nominal bin width.
        if (b.high[j]-b.low[j]>gap) return -1;
    }
    return 1;
}
}
bool covers_orbit(const std::vector<double> &phase, std::size_t bins, double gap, double fraction) {
    if (phase.empty() || bins<4 || !(gap>0 && gap<=tau) || !(fraction>0 && fraction<=1))
        throw std::invalid_argument("Invalid coverage inputs");
    Bins b(bins); b.clear();
    for (double p:phase) { if(!std::isfinite(p)) throw std::invalid_argument("Nonfinite phase"); b.add(p); }
    const int decision=bin_decision(b,gap,fraction);
    return decision<0 ? sorted_gap(phase)<=gap : decision!=0;
}

CoverageResult run_coverage(Config config, const CoverageOptions &options) {
    validate_config(config);
    if (!(options.intervals_per_phase_sigma>=10 && options.intervals_per_phase_sigma<=10000) ||
        !std::isfinite(options.intervals_per_phase_sigma) || options.max_extensions>10 ||
        std::ceil(1.6*options.intervals_per_phase_sigma)*std::pow(2.0,options.max_extensions)>1000000)
        throw std::invalid_argument("Invalid coverage cadence scaling");
    const auto start=std::chrono::steady_clock::now();
    CoverageResult result;
    auto finish=[&] { result.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); return result; };
    // The phase scanner does not retain Cartesian states or frame-major output.
    // Use a two-frame memory preflight for the sampler, then retain O(N) arrays.
    Config sampling=config; sampling.duration_days=1; sampling.output_step_days=1; sampling.visual_samples=0;
    const auto initial=sample_initial_for_domain_audit(sampling);
    result.samples=initial.size();
    result.min_initial_perigee_km=std::numeric_limits<double>::infinity();
    for(const auto &e:initial) {
        const double h=(e[0]*(1-std::hypot(e[1],e[2]))-config.backend.earth_radius_m)/1000;
        result.min_initial_perigee_km=std::min(result.min_initial_perigee_km,h);
        result.earth_intersections += h<0;
        result.particles_below_300km += h<300;
        result.particles_below_500km += h<500;
        result.particles_below_1000km += h<1000;
    }
    if(result.earth_intersections) { result.status="earth_intersection"; return finish(); }
    config.backend.output_type="mean";
    Backend nominal(config.backend,config.nominal);
    const auto reference=nominal.mean_phase_law();
    std::vector<double> phase0(initial.size()), rate(initial.size()), mean_perigee(initial.size());
    result.threads_used=std::min(initial.size(), config.threads ? config.threads : std::size_t(std::max(1u,std::thread::hardware_concurrency())));
    std::atomic<std::size_t> next{0}; std::atomic<bool> failed{false};
    std::exception_ptr exception; std::mutex mutex;
    auto worker=[&] {
        try {
            MeanPhaseFactory factory(config.backend);
            while(!failed) {
                const auto id=next.fetch_add(1); if(id>=initial.size()) break;
                const auto model=factory.prepare(initial[id]);
                phase0[id]=model.mean_at_epoch[5]-reference.mean_at_epoch[5];
                rate[id]=model.longitude_rate_rad_s-reference.longitude_rate_rad_s;
                const auto &e=model.mean_at_epoch;
                mean_perigee[id]=(e[0]*(1-std::hypot(e[1],e[2]))-config.backend.earth_radius_m)/1000;
            }
        } catch(...) { failed=true; std::lock_guard<std::mutex> lock(mutex); if(!exception) exception=std::current_exception(); }
    };
    std::vector<std::thread> workers;
    try { for(std::size_t t=0;t<result.threads_used;++t) workers.emplace_back(worker); }
    catch(...) { failed=true; for(auto &t:workers)t.join(); throw; }
    for(auto &t:workers)t.join();
    if(exception)std::rethrow_exception(exception);
    result.min_mean_perigee_km=*std::min_element(mean_perigee.begin(),mean_perigee.end());
    if(result.min_mean_perigee_km<0) { result.status="mean_earth_intersection"; return finish(); }
    result.phase_rate_sigma=sample_sigma(rate);
    if(options.automatic_times && result.phase_rate_sigma<=std::numeric_limits<double>::min()) {
        result.status="no_phase_shear";
        if(covers_orbit(phase0,config.phase_bins,config.coverage_max_gap_deg*pi/180,config.coverage_occupied_fraction)) {
            result.status="ok";
            result.cadence_s=config.output_step_days*day;
            result.confirmation_s=(config.persistence-1)*result.cadence_s;
            result.horizon_s=result.confirmation_s;
            result.evaluated_epochs=config.persistence;
            result.coverage={true,0,0,config.persistence-1};
            result.onset_max_gap_deg=sorted_gap(phase0)*180/pi;
        }
        return finish();
    }
    result.cadence_s=options.automatic_times ? 1/(options.intervals_per_phase_sigma*result.phase_rate_sigma) : config.output_step_days*day;
    const auto initial_intervals=options.automatic_times ? static_cast<std::size_t>(std::ceil(1.6*options.intervals_per_phase_sigma))
        : static_cast<std::size_t>(std::ceil(config.duration_days*day/result.cadence_s));
    if(!std::isfinite(result.cadence_s) || result.cadence_s<=0 ||
        (options.automatic_times && result.cadence_s*initial_intervals*std::pow(2.0,options.max_extensions)>1e13))
        throw std::invalid_argument("Coverage time grid is nonfinite or exceeds the supported horizon");
    std::size_t intervals=initial_intervals, run=0, onset=0;
    Bins bins(config.phase_bins); std::vector<double> phases(initial.size());
    const double gap=config.coverage_max_gap_deg*pi/180;
    double previous=0, onset_lower=0, onset_time=0;
    for(std::size_t epoch=0;;++epoch) {
        if(epoch>intervals) {
            if(!options.automatic_times || intervals>=initial_intervals*(std::size_t(1)<<options.max_extensions)) break;
            intervals*=2;
        }
        const double t=options.automatic_times ? epoch*result.cadence_s : std::min(config.duration_days*day,epoch*result.cadence_s);
        if(epoch && t<=previous) break;
        bins.clear();
        for(std::size_t id=0;id<rate.size();++id) bins.add(phase0[id]+rate[id]*t);
        int decision=bin_decision(bins,gap,config.coverage_occupied_fraction);
        if(decision<0) { for(std::size_t id=0;id<rate.size();++id)phases[id]=phase0[id]+rate[id]*t; decision=sorted_gap(phases)<=gap; }
        ++result.evaluated_epochs; result.horizon_s=t;
        if(decision) {
            if(!run) { onset=epoch; onset_lower=previous; onset_time=t; }
            ++run;
            if(run>=config.persistence) {
                result.coverage={true,onset ? onset_lower : 0,onset_time,epoch};
                result.confirmation_s=t; result.status="ok";
                for(std::size_t id=0;id<rate.size();++id)phases[id]=phase0[id]+rate[id]*onset_time;
                result.onset_max_gap_deg=sorted_gap(phases)*180/pi;
                return finish();
            }
        } else run=0;
        previous=t;
    }
    return finish();
}
} // namespace distribution
