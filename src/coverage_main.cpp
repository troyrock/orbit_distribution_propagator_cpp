#include <distribution/coverage.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<std::string> fields(const std::string &line) {
    std::vector<std::string> result; std::istringstream input(line); std::string value;
    while(std::getline(input,value,',')) { if(!value.empty() && value.back()=='\r')value.pop_back(); result.push_back(value); }
    // getline does not emit the final empty field after a trailing delimiter.
    // Preserve it so the exact header and six-column checks reject it.
    const auto end=line.size()-(!line.empty() && line.back()=='\r' ? 1u : 0u);
    if(end && line[end-1]==',')result.emplace_back();
    return result;
}
double number(const std::string &text) {
    std::size_t n=0; const double x=std::stod(text,&n);
    if(n!=text.size() || !std::isfinite(x))throw std::invalid_argument("Invalid finite scalar: "+text);
    return x;
}
void scalar(std::ostream &out,double value) { if(std::isfinite(value))out<<value; else out<<"null"; }
void write(std::ostream &out,const std::vector<std::string> &row,const distribution::CoverageResult &r,const distribution::Config &c) {
    out<<std::setprecision(17)<<"{\"run_id\":\""<<row[0]<<"\",\"perigee_altitude_km\":"<<number(row[1])
       <<",\"inclination_deg\":"<<number(row[2])<<",\"eccentricity\":"<<number(row[3])
       <<",\"position_sigma_km\":"<<number(row[4])<<",\"velocity_sigma_m_s\":"<<number(row[5])
       <<",\"samples\":"<<r.samples<<",\"seed\":"<<c.seed<<",\"status\":\""<<r.status<<"\""
       <<",\"coverage_time_days\":";
    if(r.coverage.found)out<<r.coverage.upper_s/distribution::day; else out<<"null";
    out<<",\"coverage_lower_days\":"; if(r.coverage.found)out<<r.coverage.lower_s/distribution::day;else out<<"null";
    out<<",\"coverage_upper_days\":"; if(r.coverage.found)out<<r.coverage.upper_s/distribution::day;else out<<"null";
    out<<",\"coverage_confirmation_days\":";if(r.coverage.found)out<<r.confirmation_s/distribution::day;else out<<"null";
    out<<",\"onset_max_gap_deg\":"; if(r.coverage.found)out<<r.onset_max_gap_deg;else out<<"null";
    out<<",\"min_initial_perigee_km\":"; scalar(out,r.min_initial_perigee_km);
    out<<",\"min_mean_perigee_km\":"; if(r.status=="earth_intersection")out<<"null";else scalar(out,r.min_mean_perigee_km);
    out<<",\"earth_intersections\":"<<r.earth_intersections<<",\"particles_below_300km\":"<<r.particles_below_300km
       <<",\"particles_below_500km\":"<<r.particles_below_500km<<",\"particles_below_1000km\":"<<r.particles_below_1000km
       <<",\"cadence_days\":"<<r.cadence_s/distribution::day<<",\"evaluated_horizon_days\":"<<r.horizon_s/distribution::day
       <<",\"phase_rate_sigma_rad_s\":"<<r.phase_rate_sigma<<",\"evaluated_epochs\":"<<r.evaluated_epochs
       <<",\"threads_used\":"<<r.threads_used<<",\"elapsed_seconds\":"<<r.elapsed_seconds<<"}\n";
    out.flush();
}
}
int main(int argc,char **argv) {
    try {
        std::filesystem::path config_path,cases_path,output_path;
        std::size_t threads=0; distribution::CoverageOptions options;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--help") {
                std::cout<<"distribution_coverage --config BASE.cfg --cases CASES.csv --output NEW.jsonl [--threads N] [--fixed-times]\n"
                    <<"CSV header: run_id,perigee_altitude_km,inclination_deg,eccentricity,position_sigma_km,velocity_sigma_m_s\n"
                    <<"Uses the exact native mean-phase flow for fixed J2/J2-squared gravity. No drag; invalid ensembles are explicitly masked.\n";
                return 0;
            }
            if(arg=="--fixed-times") { options.automatic_times=false; continue; }
            if(++i>=argc)throw std::invalid_argument("Missing argument");
            if(arg=="--config")config_path=argv[i];
            else if(arg=="--cases")cases_path=argv[i];
            else if(arg=="--output")output_path=argv[i];
            else if(arg=="--threads") { const double x=number(argv[i]); if(x<1||x>1024||std::floor(x)!=x)throw std::invalid_argument("Invalid threads"); threads=static_cast<std::size_t>(x); }
            else throw std::invalid_argument("Unknown option: "+arg);
        }
        if(config_path.empty()||cases_path.empty()||output_path.empty())throw std::invalid_argument("Require config, cases and output");
        if(std::filesystem::exists(output_path))throw std::invalid_argument("Output already exists");
        auto config=distribution::read_config(config_path);
        if(!config.empirical_samples_csv.empty())throw std::invalid_argument("Coverage grid requires Gaussian covariance input");
        config.uncertainty_coordinates="rtn";config.minimum_altitude_m=0;config.visual_samples=0;
        if(threads)config.threads=threads;
        std::ifstream input(cases_path);if(!input)throw std::runtime_error("Cannot open cases CSV");
        std::string line;std::getline(input,line);
        const auto header=fields(line);
        if(header!=std::vector<std::string>{"run_id","perigee_altitude_km","inclination_deg","eccentricity","position_sigma_km","velocity_sigma_m_s"})
            throw std::invalid_argument("Wrong cases CSV header");
        std::vector<std::vector<std::string>> cases;std::set<std::string> ids;
        while(std::getline(input,line)) {
            if(line.empty()||line=="\r")continue;
            auto row=fields(line);if(row.size()!=6)throw std::invalid_argument("Wrong case column count");
            if(row[0].empty()||row[0].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos||!ids.insert(row[0]).second)
                throw std::invalid_argument("Invalid or duplicate run ID");
            for(std::size_t j=1;j<6;++j)number(row[j]);
            if(number(row[1])<0||number(row[2])<0||number(row[2])>90||number(row[3])<0||number(row[3])>=1||number(row[4])<=0||number(row[5])<=0)
                throw std::invalid_argument("Case outside supported positive Gaussian / prograde grid domain");
            cases.push_back(std::move(row));
        }
        if(cases.empty())throw std::invalid_argument("No cases");
        if(!output_path.parent_path().empty())std::filesystem::create_directories(output_path.parent_path());
        std::ofstream out(output_path);out.exceptions(std::ios::failbit|std::ios::badbit);
        std::size_t complete=0;
        for(const auto &row:cases) {
            const double ecc=number(row[3]),a=(config.backend.earth_radius_m+number(row[1])*1000)/(1-ecc);
            config.nominal=distribution::from_keplerian(a,ecc,number(row[2])*distribution::pi/180,
                20*distribution::pi/180,30*distribution::pi/180,10*distribution::pi/180);
            config.covariance={};for(int j=0;j<6;++j) { const double sigma=j<3?number(row[4])*1000:number(row[5]); config.covariance[j][j]=sigma*sigma; }
            const auto result=distribution::run_coverage(config,options);
            write(out,row,result,config);
            std::cout<<++complete<<'/'<<cases.size()<<' '<<row[0]<<' '<<result.status<<" t_days="
                <<(result.coverage.found?result.coverage.upper_s/distribution::day:-1)<<" seconds="<<result.elapsed_seconds<<'\n'<<std::flush;
        }
        return 0;
    } catch(const std::exception &e) { std::cerr<<"Error: "<<e.what()<<'\n';return 1; }
}
