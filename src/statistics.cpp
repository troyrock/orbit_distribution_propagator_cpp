#include <distribution/statistics.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace distribution {
double quantile_sorted(const std::vector<double>& x,double p) {
    if(x.empty() || !(p>=0 && p<=1)) throw std::invalid_argument("Invalid quantile input");
    const double at=p*static_cast<double>(x.size()-1);const auto i=static_cast<std::size_t>(at);
    return x[i]+(at-static_cast<double>(i))*(x[std::min(i+1,x.size()-1)]-x[i]);
}
double sample_sigma(const std::vector<double>& x) {
    for(double v:x) if(!std::isfinite(v))throw std::invalid_argument("Nonfinite statistic");
    if(x.size()<2) return 0;
    double mean=0,m2=0;std::size_t n=0;
    for(double v:x) {if(!std::isfinite(v))throw std::invalid_argument("Nonfinite statistic");const double d=v-mean;mean+=d/static_cast<double>(++n);m2+=d*(v-mean);}
    return std::sqrt(std::max(0.0,m2/static_cast<double>(x.size()-1)));
}
PhaseMetrics phase_metrics(const std::vector<double>& phase,std::size_t bins) {
    if(phase.empty() || bins<4) throw std::invalid_argument("Need samples and >=4 phase bins");
    PhaseMetrics m;m.sigma=sample_sigma(phase);m.histogram.resize(bins);
    auto sorted=phase;std::sort(sorted.begin(),sorted.end());
    m.q95_width=quantile_sorted(sorted,0.975)-quantile_sorted(sorted,0.025);
    std::array<long double,4> cs{},sn{};
    for(auto& x:sorted) {
        x=wrap_angle(x);
        ++m.histogram[std::min(bins-1,static_cast<std::size_t>(x/tau*static_cast<double>(bins)))];
        for(std::size_t k=0;k<4;++k) {cs[k]+=std::cos(static_cast<double>(k+1)*x);sn[k]+=std::sin(static_cast<double>(k+1)*x);}
    }
    for(std::size_t k=0;k<4;++k)m.resultants[k]=static_cast<double>(std::hypot(cs[k],sn[k])/phase.size());
    std::sort(sorted.begin(),sorted.end());
    m.max_gap=sorted.front()+tau-sorted.back();
    for(std::size_t i=1;i<sorted.size();++i)m.max_gap=std::max(m.max_gap,sorted[i]-sorted[i-1]);
    std::size_t occupied=0;for(auto count:m.histogram) {occupied+=count>0;m.total_variation+=0.5*std::abs(static_cast<double>(count)/static_cast<double>(phase.size())-1.0/static_cast<double>(bins));}
    m.occupied_fraction=static_cast<double>(occupied)/static_cast<double>(bins);
    return m;
}
double dkw_epsilon(std::size_t n,double alpha) {
    if(n==0 || !(alpha>0 && alpha<1))throw std::invalid_argument("Invalid DKW inputs");
    return std::sqrt(std::log(2/alpha)/(2*static_cast<double>(n)));
}
EventInterval first_sustained(const std::vector<double>& times,const std::vector<bool>& flags,std::size_t persistence) {
    if(times.size()!=flags.size() || persistence==0) throw std::invalid_argument("Invalid event series");
    std::size_t run=0;
    for(std::size_t i=0;i<flags.size();++i) {
        run=flags[i]?run+1:0;
        if(run>=persistence) {const auto first=i+1-persistence;return {true,times[first?first-1:0],times[first],i};}
    }
    return {};
}
}
