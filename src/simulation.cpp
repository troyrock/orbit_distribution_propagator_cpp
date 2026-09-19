#include <distribution/simulation.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <thread>

namespace distribution {
namespace {
std::array<double,3> sigmas(const std::array<std::vector<double>,3>& x) {
    return {sample_sigma(x[0]),sample_sigma(x[1]),sample_sigma(x[2])};
}
double estimate_mixing(const std::vector<SampleState>& states,std::size_t n,const Config& c) {
    // Gaussian Kepler shear approximation, using sampled INITIAL MEAN orbits.
    double mp=0,mn=0,vp=0,vn=0,cpn=0;
    const double origin=states[0].mean_longitude;
    for(std::size_t i=0;i<n;++i) {
        const double p=states[i].mean_longitude-origin,rate=std::sqrt(c.backend.mu/states[i].mean_a)/states[i].mean_a;
        const double dp=p-mp,dn=rate-mn,div=static_cast<double>(i+1);
        mp+=dp/div;mn+=dn/div;vp+=dp*(p-mp);vn+=dn*(rate-mn);cpn+=dp*(rate-mn);
    }
    vp/=static_cast<double>(n-1);vn/=static_cast<double>(n-1);cpn/=static_cast<double>(n-1);
    const double target=-2*std::log(c.mixing_max_resultant);
    if(vp>=target)return 0;
    if(vn<=std::numeric_limits<double>::min())return -1;
    const double root=std::sqrt(cpn*cpn+vn*(target-vp));
    // Stable positive quadratic root.
    return cpn>=0 ? (target-vp)/(root+cpn) : (root-cpn)/vn;
}
}
Simulation run_simulation(Config config) {
    validate_config(config);const auto start=std::chrono::steady_clock::now();
    Simulation s;s.config=config;s.initial=sample_initial(config);s.config.samples=s.initial.size();
    const auto n=s.initial.size();const auto times=output_times(config);
    const long double bytes=estimated_memory_bytes(config,n);
    if(bytes>config.max_memory_mb*1024*1024)throw std::runtime_error("Estimated state/display storage exceeds max_memory_mb; reduce samples/snapshots or increase the budget");
    s.states.resize(n*times.size());s.frames.resize(times.size());
    Backend nominal(config.backend,config.nominal);
    for(std::size_t f=0;f<times.size();++f) {s.frames[f].time_s=times[f];s.frames[f].reference=nominal.advance(times[f]);}
    s.threads_used=std::min(n,config.threads?config.threads:std::max(1u,std::thread::hardware_concurrency()));
    std::atomic<std::size_t> next{0},completed{0};std::atomic<bool> failed{false};
    std::mutex result_mutex;std::exception_ptr exception;
    auto worker=[&] {
        BackendStats totals;
        try {
            while(!failed) {
                const auto id=next.fetch_add(1);if(id>=n)break;
                try {
                    Backend backend(config.backend,s.initial[id]);
                    for(std::size_t f=0;f<times.size();++f) {
                        const auto state=backend.advance(times[f]);
                        validate_elements(state.osculating,config.backend.mu,config.backend.earth_radius_m+config.minimum_altitude_m);
                        s.states[f*n+id]={state.osculating,state.mean[5],state.mean[0]};
                    }
                    const auto& stats=backend.stats();totals.accepted_steps+=stats.accepted_steps;totals.rejected_steps+=stats.rejected_steps;totals.derivative_evaluations+=stats.derivative_evaluations;
                    completed.fetch_add(1);
                }catch(const std::exception& e) {throw std::runtime_error("Particle "+std::to_string(id)+": "+e.what());}
            }
        }catch(...) {failed=true;std::lock_guard<std::mutex> lock(result_mutex);if(!exception)exception=std::current_exception();}
        std::lock_guard<std::mutex> lock(result_mutex);s.stats.accepted_steps+=totals.accepted_steps;s.stats.rejected_steps+=totals.rejected_steps;s.stats.derivative_evaluations+=totals.derivative_evaluations;
    };
    std::vector<std::thread> workers;
    try {for(std::size_t i=0;i<s.threads_used;++i)workers.emplace_back(worker);}
    catch(...) {failed=true;for(auto& w:workers)w.join();throw;}
    double last_report=0;
    while(completed<n && !failed) {
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        // Progress is emitted from the main thread, never interleaved workers.
        if(elapsed-last_report>=5) {last_report=elapsed;std::cerr<<"Propagated "<<completed.load()<<" / "<<n<<" particles\r";}
    }
    for(auto& w:workers)w.join();
    if(exception)std::rethrow_exception(exception);
    std::cerr<<"Propagated "<<n<<" / "<<n<<" particles\n";
    std::vector<bool> coverage,mixing,central95;
    for(std::size_t f=0;f<times.size();++f) {
        auto& frame=s.frames[f];const auto& ref=frame.reference.osculating;
        const auto refcart=to_cartesian(ref,config.backend.mu);const auto basis=rtn_basis(refcart);
        std::vector<double> phase(n),axes(n);std::array<std::vector<double>,3> rtn,tube;
        for(int k=0;k<3;++k) {rtn[k].resize(n);tube[k].resize(n);}
        frame.positions.reserve(std::min(config.visual_samples,n));
        for(std::size_t id=0;id<n;++id) {
            const auto& state=s.states[f*n+id];phase[id]=state.mean_longitude-frame.reference.mean[5];axes[id]=state.mean_a;
            const auto cart=to_cartesian(state.output,config.backend.mu);
            if(id<config.visual_samples)frame.positions.push_back({cart[0],cart[1],cart[2]});
            Vector3 delta{cart[0]-refcart[0],cart[1]-refcart[1],cart[2]-refcart[2]};
            for(int k=0;k<3;++k)rtn[k][id]=dot(delta,basis[k]);
            // At the SAME output mean longitude, compare to the reference
            // ellipse. These residuals remove phase spread before describing
            // tube thickness; they are not a substitute for the full cloud.
            auto matched=ref;matched[5]=state.output[5];const auto match_cart=to_cartesian(matched,config.backend.mu);
            const auto match_basis=rtn_basis(match_cart);
            delta={cart[0]-match_cart[0],cart[1]-match_cart[1],cart[2]-match_cart[2]};
            for(int k=0;k<3;++k)tube[k][id]=dot(delta,match_basis[k]);
        }
        frame.phase=phase_metrics(phase,config.phase_bins);frame.rtn_sigma_m=sigmas(rtn);frame.tube_rtn_sigma_m=sigmas(tube);frame.semimajor_sigma_m=sample_sigma(axes);
        frame.coverage=frame.phase.max_gap<=config.coverage_max_gap_deg*pi/180 && frame.phase.occupied_fraction>=config.coverage_occupied_fraction;
        frame.mixed=frame.coverage && *std::max_element(frame.phase.resultants.begin(),frame.phase.resultants.end())<=config.mixing_max_resultant && frame.phase.total_variation<=config.mixing_max_tv;
        coverage.push_back(frame.coverage);mixing.push_back(frame.mixed);central95.push_back(frame.phase.q95_width>=tau);
        for(int k=0;k<=180;++k) {auto orbit=ref;orbit[5]=tau*k/180;const auto cart=to_cartesian(orbit,config.backend.mu);frame.reference_orbit.push_back({cart[0],cart[1],cart[2]});}
    }
    s.coverage=first_sustained(times,coverage,config.persistence);s.mixing=first_sustained(times,mixing,config.persistence);s.central95_wrap=first_sustained(times,central95,config.persistence);
    s.analytic_mixing_time_s=estimate_mixing(s.states,n,config);
    if(config.accuracy_check) {
        auto tighter=config.backend;tighter.relative_tolerance*=0.1;tighter.absolute_tolerance_m*=0.1;tighter.absolute_tolerance_elements*=0.1;tighter.max_step_s=std::max(tighter.min_step_s,tighter.max_step_s/4);
        for(std::size_t id=0;id<std::min(config.accuracy_check,n);++id) {
            Backend check(tighter,s.initial[id]);
            for(std::size_t f=0;f<times.size();++f) {
                const auto fine=check.advance(times[f]);const auto& coarse=s.states[f*n+id];
                const auto a=to_cartesian(fine.osculating,config.backend.mu),b=to_cartesian(coarse.output,config.backend.mu);
                s.accuracy_max_position_m=std::max(s.accuracy_max_position_m,norm({a[0]-b[0],a[1]-b[1],a[2]-b[2]}));
                s.accuracy_max_phase_rad=std::max(s.accuracy_max_phase_rad,std::abs(fine.mean[5]-coarse.mean_longitude));
            }
        }
    }
    s.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return s;
}
}
