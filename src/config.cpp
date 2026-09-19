#include <distribution/config.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace distribution {
namespace {
std::string trim(std::string x) {
    const auto first=x.find_first_not_of(" \t\r\n");if(first==std::string::npos)return {};
    return x.substr(first,x.find_last_not_of(" \t\r\n")-first+1);
}
std::vector<double> numbers(std::string s) {
    std::replace(s.begin(),s.end(),',',' ');std::istringstream in(s);std::vector<double> out;double x;
    while(in>>x) {if(!std::isfinite(x))throw std::invalid_argument("Nonfinite configuration value");out.push_back(x);}
    if(!in.eof())throw std::invalid_argument("Invalid numeric list: "+s);
    return out;
}
Elements six(const std::string& text) {
    const auto v=numbers(text);if(v.size()!=6)throw std::invalid_argument("Expected exactly six numeric values");
    Elements e{};std::copy(v.begin(),v.end(),e.begin());return e;
}
double number(const std::string& text) {
    const auto v=numbers(text);if(v.size()!=1)throw std::invalid_argument("Expected one number");return v[0];
}
std::uint64_t integer(const std::string& text) {
    if(text.empty() || text[0]=='-')throw std::invalid_argument("Expected nonnegative integer");
    std::size_t consumed=0;const auto value=std::stoull(text,&consumed);
    if(consumed!=text.size())throw std::invalid_argument("Invalid integer: "+text);
    return value;
}
Matrix6 read_covariance(const std::filesystem::path& file) {
    std::ifstream in(file);if(!in)throw std::runtime_error("Cannot read covariance file: "+file.string());
    Matrix6 m{};std::string line;std::size_t row=0;
    while(std::getline(in,line)) {line=trim(line.substr(0,line.find('#')));if(line.empty())continue;
        if(row==6)throw std::invalid_argument("Covariance must have six rows");
        m[row++]=six(line);
    }
    if(row!=6)throw std::invalid_argument("Covariance must have six rows");
    return m;
}
// Box-Muller with explicitly specified uniform mapping: deterministic sample
// streams on the same IEEE-754 platform, independent of worker scheduling.
class Normal {
    std::mt19937_64 engine;bool spare_ready=false;double spare=0;
public:
    explicit Normal(std::uint64_t seed):engine(seed){}
    double get() {
        if(spare_ready) {spare_ready=false;return spare;}
        const double u=(static_cast<double>(engine()>>11)+0.5)/9007199254740992.0;
        const double v=(static_cast<double>(engine()>>11)+0.5)/9007199254740992.0;
        const double r=std::sqrt(-2*std::log(u)),angle=tau*v;
        spare=r*std::sin(angle);spare_ready=true;return r*std::cos(angle);
    }
};
Elements apply_offset(const Config& c,const Elements& offset) {
    if(c.uncertainty_coordinates=="equinoctial") {auto e=c.nominal;for(int i=0;i<6;++i)e[i]+=offset[i];return e;}
    auto cart=to_cartesian(c.nominal,c.backend.mu);
    if(c.uncertainty_coordinates=="cartesian")for(int i=0;i<6;++i)cart[i]+=offset[i];
    else {
        const auto basis=rtn_basis(cart);
        // Rotate inertial position/velocity error components by the SAME
        // nominal epoch basis. This is not the derivative of rotating position.
        for(int i=0;i<3;++i)for(int j=0;j<3;++j) {cart[i]+=basis[j][i]*offset[j];cart[i+3]+=basis[j][i]*offset[j+3];}
    }
    auto e=from_cartesian(cart,c.backend.mu);
    e[5]=c.nominal[5]+std::remainder(e[5]-c.nominal[5],tau);return e;
}
}
Config::Config() {
    const Elements sigma{100000,0.0001,0.0001,0.00005,0.00005,0.0001};
    for(int i=0;i<6;++i)covariance[i][i]=sigma[i]*sigma[i];
}
Config read_config(const std::filesystem::path& path) {
    std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read configuration: "+path.string());
    std::map<std::string,std::string> kv;std::string line;std::size_t ln=0;
    while(std::getline(in,line)) {++ln;line=trim(line.substr(0,line.find('#')));if(line.empty())continue;
        const auto at=line.find('=');if(at==std::string::npos)throw std::invalid_argument("Missing '=' at config line "+std::to_string(ln));
        const auto key=trim(line.substr(0,at)),value=trim(line.substr(at+1));
        if(key.empty() || value.empty() || !kv.emplace(key,value).second)throw std::invalid_argument("Empty or duplicate config key: "+key);
    }
    Config c;bool supplied_cov=false;
    const auto take=[&](const std::string& key)->std::string {auto it=kv.find(key);if(it==kv.end())return {};auto value=it->second;kv.erase(it);return value;};
    const auto scalar=[&](const char* key,double& dest){const auto v=take(key);if(!v.empty())dest=number(v);};
    const auto count=[&](const char* key,std::size_t& dest){const auto v=take(key);if(!v.empty()) {const auto n=integer(v);if(n>std::numeric_limits<std::size_t>::max())throw std::invalid_argument("Integer overflow");dest=static_cast<std::size_t>(n);}};
    const auto word=[&](const char* key,std::string& dest){const auto v=take(key);if(!v.empty())dest=v;};
    word("force_model",c.backend.force_model);word("initial_type",c.backend.initial_type);word("output_type",c.backend.output_type);word("uncertainty_coordinates",c.uncertainty_coordinates);
    scalar("mu",c.backend.mu);scalar("earth_radius_m",c.backend.earth_radius_m);scalar("j2",c.backend.j2);
    scalar("relative_tolerance",c.backend.relative_tolerance);scalar("absolute_tolerance_m",c.backend.absolute_tolerance_m);scalar("absolute_tolerance_elements",c.backend.absolute_tolerance_elements);
    scalar("min_step_s",c.backend.min_step_s);scalar("max_step_s",c.backend.max_step_s);
    scalar("duration_days",c.duration_days);scalar("output_step_days",c.output_step_days);scalar("minimum_altitude_m",c.minimum_altitude_m);
    scalar("coverage_max_gap_deg",c.coverage_max_gap_deg);scalar("coverage_occupied_fraction",c.coverage_occupied_fraction);
    scalar("mixing_max_resultant",c.mixing_max_resultant);scalar("mixing_max_tv",c.mixing_max_tv);scalar("max_memory_mb",c.max_memory_mb);
    count("samples",c.samples);count("threads",c.threads);count("visual_samples",c.visual_samples);count("phase_bins",c.phase_bins);count("persistence",c.persistence);count("accuracy_check",c.accuracy_check);
    const auto seed=take("seed");if(!seed.empty())c.seed=integer(seed);
    const auto kep=take("orbit_keplerian_deg"),eq=take("orbit_equinoctial"),cart=take("orbit_cartesian");
    if((!kep.empty())+(!eq.empty())+(!cart.empty())>1)throw std::invalid_argument("Specify only one nominal orbit format");
    if(!kep.empty()) {const auto v=six(kep);c.nominal=from_keplerian(v[0],v[1],v[2]*pi/180,v[3]*pi/180,v[4]*pi/180,v[5]*pi/180);}
    if(!eq.empty())c.nominal=six(eq);
    if(!cart.empty())c.nominal=from_cartesian(six(cart),c.backend.mu);
    const auto sig=take("sigma"),cov=take("covariance_csv"),emp=take("empirical_samples_csv");
    if((!sig.empty())+(!cov.empty())+(!emp.empty())>1)throw std::invalid_argument("Specify one of sigma, covariance_csv, empirical_samples_csv");
    if(!sig.empty()) {c.covariance={};const auto v=six(sig);for(int i=0;i<6;++i) {if(v[i]<0)throw std::invalid_argument("sigma cannot be negative");c.covariance[i][i]=v[i]*v[i];}supplied_cov=true;}
    if(!cov.empty()) {c.covariance=read_covariance(path.parent_path()/cov);supplied_cov=true;}
    if(!emp.empty())c.empirical_samples_csv=path.parent_path()/emp;
    if(c.uncertainty_coordinates!="equinoctial" && !supplied_cov && emp.empty())throw std::invalid_argument("Non-equinoctial uncertainty requires explicit sigma or covariance_csv");
    if(!kv.empty())throw std::invalid_argument("Unknown config key: "+kv.begin()->first);
    validate_config(c);return c;
}
void validate_config(const Config& c) {
    const auto positive=[](double x){return std::isfinite(x)&&x>0;};
    const auto& b=c.backend;
    if(!positive(b.mu) || !positive(b.earth_radius_m) || !positive(b.relative_tolerance) || !positive(b.absolute_tolerance_m) || !positive(b.absolute_tolerance_elements) || !positive(b.min_step_s) || !positive(b.max_step_s) || b.min_step_s>b.max_step_s || !std::isfinite(b.j2) || b.j2<0)
        throw std::invalid_argument("Invalid DSST constants, tolerance or step limits");
    if((b.force_model!="kepler" && b.force_model!="j2" && b.force_model!="j2_j2sq") || (b.initial_type!="mean" && b.initial_type!="osculating") || (b.output_type!="mean" && b.output_type!="osculating"))
        throw std::invalid_argument("Invalid force_model, initial_type or output_type");
    if(c.samples<2 || c.samples>100000000 || c.threads>1024 || c.phase_bins<4 || c.phase_bins>100000 || c.persistence==0)
        throw std::invalid_argument("Invalid sample/thread/bin/persistence count");
    if(!positive(c.duration_days) || !positive(c.output_step_days) || !positive(c.max_memory_mb) || c.duration_days*day>1e13 || c.duration_days/c.output_step_days>1000000)
        throw std::invalid_argument("Invalid/excessive duration, cadence or memory budget");
    if(!std::isfinite(c.minimum_altitude_m) || c.minimum_altitude_m<0)throw std::invalid_argument("minimum_altitude_m must be nonnegative");
    if(!(c.coverage_max_gap_deg>0 && c.coverage_max_gap_deg<=360 && c.coverage_occupied_fraction>0 && c.coverage_occupied_fraction<=1 && c.mixing_max_resultant>0 && c.mixing_max_resultant<1 && c.mixing_max_tv>0 && c.mixing_max_tv<1))
        throw std::invalid_argument("Invalid coverage/mixing thresholds");
    if(c.uncertainty_coordinates!="equinoctial" && c.uncertainty_coordinates!="cartesian" && c.uncertainty_coordinates!="rtn")throw std::invalid_argument("uncertainty_coordinates must be equinoctial, cartesian, or rtn");
    validate_elements(c.nominal,c.backend.mu,c.backend.earth_radius_m+c.minimum_altitude_m);
    covariance_factor(c.covariance);
}
std::vector<double> output_times(const Config& c) {
    const double end=c.duration_days*day,step=c.output_step_days*day;
    std::vector<double> times{0};const auto count=static_cast<std::size_t>(std::floor(end/step));
    for(std::size_t i=1;i<=count;++i) {const double t=static_cast<double>(i)*step;if(t<end)times.push_back(t);}
    times.push_back(end);return times;
}
long double estimated_memory_bytes(const Config& c,std::size_t n) {
    const long double frames=std::ceil(c.duration_days/c.output_step_days)+1;
    const long double displayed=std::min(c.visual_samples,n);
    // Retained states, initial elements, statistics scratch, histograms,
    // reference curves, JSON/HTML copies; reserve for force objects/overhead.
    const auto workers=c.threads?c.threads:std::max(1u,std::thread::hardware_concurrency());
    return frames*(static_cast<long double>(n)*64 + displayed*24 + 181*24 + static_cast<long double>(c.phase_bins)*sizeof(std::size_t)+1024)
           +static_cast<long double>(n)*128
           +frames*(displayed*90 + 181*90 + static_cast<long double>(c.phase_bins)*24+2048)*4
           +(static_cast<long double>(std::min(n,workers))*2+8)*1024*1024;
}
std::vector<Elements> sample_initial(const Config& c) {
    std::vector<Elements> out;
    if(!c.empirical_samples_csv.empty()) {
        std::ifstream in(c.empirical_samples_csv);if(!in)throw std::runtime_error("Cannot read empirical ensemble");
        std::string line;std::size_t row=0;
        while(std::getline(in,line)) {++row;line=trim(line.substr(0,line.find('#')));if(line.empty())continue;
            try {auto v=six(line);if(c.uncertainty_coordinates=="cartesian") {v=from_cartesian(v,c.backend.mu);v[5]=c.nominal[5]+std::remainder(v[5]-c.nominal[5],tau);}else if(c.uncertainty_coordinates=="rtn")v=apply_offset(c,v);
                if(estimated_memory_bytes(c,out.size()+1)>c.max_memory_mb*1024*1024)throw std::invalid_argument("Empirical ensemble exceeds max_memory_mb");
                validate_elements(v,c.backend.mu,c.backend.earth_radius_m+c.minimum_altitude_m);out.push_back(v);
            } catch(const std::exception& e) {throw std::runtime_error("Invalid empirical row "+std::to_string(row)+": "+e.what());}
        }
        if(out.size()<2)throw std::invalid_argument("Empirical ensemble requires >=2 equal-weight samples");
    } else {
        if(estimated_memory_bytes(c,c.samples)>c.max_memory_mb*1024*1024)throw std::invalid_argument("Estimated run storage exceeds max_memory_mb before sampling");
        out.reserve(c.samples);Normal rng(c.seed);const auto lower=covariance_factor(c.covariance);
        for(std::size_t sample=0;sample<c.samples;++sample) {
            Elements z{},offset{};for(auto& v:z)v=rng.get();
            for(int i=0;i<6;++i)for(int j=0;j<=i;++j)offset[i]+=lower[i][j]*z[j];
            try {auto e=apply_offset(c,offset);validate_elements(e,c.backend.mu,c.backend.earth_radius_m+c.minimum_altitude_m);out.push_back(e);}
            catch(const std::exception& e) {throw std::runtime_error("Invalid sampled orbit "+std::to_string(sample)+": "+e.what()+". Supply a physically valid posterior; samples are not clipped or redrawn.");}
        }
    }
    return out;
}
std::string config_template() {return R"(# SI units except degrees explicitly named and time in days.
force_model = j2_j2sq
initial_type = osculating
output_type = osculating
orbit_keplerian_deg = 26560000, 0.02, 55, 20, 30, 10
uncertainty_coordinates = equinoctial
# Deliberately broad demonstration; replace with your actual OD covariance.
sigma = 100000, 0.0001, 0.0001, 0.00005, 0.00005, 0.0001
samples = 5000
seed = 20260919
threads = 0
duration_days = 120
output_step_days = 1
visual_samples = 1500
phase_bins = 72
coverage_max_gap_deg = 5
coverage_occupied_fraction = 1
mixing_max_resultant = 0.05
mixing_max_tv = 0.15
persistence = 3
minimum_altitude_m = 1000000
relative_tolerance = 1e-11
absolute_tolerance_m = 0.001
absolute_tolerance_elements = 1e-12
min_step_s = 0.001
max_step_s = 86400
max_memory_mb = 2048
accuracy_check = 8
)";}
}
