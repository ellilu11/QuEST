#include <Eigen/Dense>
//#include <boost/test/unit_test.hpp>
//#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <ctime>

#include "../src/configuration.h"
#include "../src/quantum_dot.h"
#include "../src/quantum_dot.cpp"
#include "../src/integrator/history.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/interactions/AIM/aim_interaction.h"
#include "../src/interactions/AIM/expansion.h"
#include "../src/interactions/AIM/grid.h"
#include "../src/interactions/AIM/normalization.h"
#include "../src/math_utils.h"

// namespace po = boost::program_options;
 
const double c0 = 299.792458, mu0 = 2.0133545e-04, hbar = 0.65821193;
const double omega = 2278.9013;
const double lambda = 2 * M_PI * c0 / omega;
const double prop_constant = // 1.00 / hbar;
                             mu0 / (4 * M_PI * hbar);
 
Eigen::Vector2cd source(double t, double mu, double sigsqr){
    return Eigen::Vector2cd(0, exp(-std::pow(t - mu, 2) / (2.0 * sigsqr)));
}

Eigen::Vector3d efld_d0_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(0, exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0); 
}

Eigen::Vector3d efld_d1_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(0, -(t - mu - delay) / sigsqr * exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0); 
}

Eigen::Vector3d efld_d2_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(0, (std::pow(t - mu - delay, 2) - sigsqr) / pow(sigsqr, 2) *
        exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0); 
}

Eigen::Vector3d analytic_EFIE_interaction(Eigen::Vector3d &efld_d0,
                                     Eigen::Vector3d &efld_d1,
                                     Eigen::Vector3d &efld_d2,
                                     Eigen::Vector3d &dr,
                                     double c0, double dist){
    Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
    Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
    Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3 * rr;

    return - pow(c0, 2) * prop_constant * hbar * 
       (i3rr * efld_d0 / std::pow(dist, 3) +
        i3rr * efld_d1 / (c0 * std::pow(dist, 2)) +
        irr * efld_d2 / (std::pow(c0, 2) * dist));
}

Eigen::Vector3d analytic_Laplace_interaction(Eigen::Vector3d &efld_d0, double dist){
    return prop_constant * hbar * Eigen::Matrix3d::Identity() * efld_d0 / dist;
}

Eigen::Vector3cd analytic_Helmholtz_interaction(Eigen::Vector3d &efld_d0, double c0, double dist){
    const std::complex<double> iu(0, 1);
    return prop_constant * hbar * Eigen::Matrix3cd::Identity() * 
        efld_d0 * std::exp( -iu * omega / c0 * dist) / dist;
}

int main(int argc, char *argv[]){

    /*po::options_description desc("Allowed options");
    desc.add_options()
        ("expansion", po::value<int>(), "expansion order")
        ("ds", po::value<double>(), "ds multiplier");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    */
    const int interp = 5;
    const int steps = atoi(argv[1]);
    const double dt = 1.0;
                     // 1e-2;   
    const double tmax = steps*dt;
    const double mu = tmax/2.0;
    const double sig = tmax/12.0;
                       // 0.01;
    const double sigsqr = sig * sig;
    // const double fmax = 6.0/(2.0*M_PI*sqrt(sigsqr));

    // set up dots and history of dots
    auto dots = std::make_shared<DotVector>(import_dots("dots.cfg"));
    const int ndots = (*dots).size();
    std::cout << "Running with " << ndots << " dots" << std::endl;

    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, 22, steps); 
    history->fill(Eigen::Vector2cd::Zero());

    const int nsrcs = ndots/2;
    const int ntrgs = ndots - nsrcs;
    for (int n = 0; n < nsrcs; ++n)
        for (int i = -22; i < steps; ++i) 
            // if ( n != OBS )
                history->array_[n][i][0] = source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm();

    // set up propagator
    std::shared_ptr<Propagation::Kernel<cmplx>> propagator = 
        // std::make_shared<Propagation::Laplace<cmplx>>();
        // std::make_shared<Propagation::Helmholtz>(prop_constant, omega / c0 );
        std::make_shared<Propagation::EFIE<cmplx>>(c0, prop_constant);

    auto direct_interaction = std::make_shared<DirectInteraction>(dots, history, propagator, interp, c0, dt);

    // set up AIM interaction
    const int expansion = atoi(argv[2]); //vm["expansion"].as<int>();
    const int border = 1;
    const double ds_base = 0.050 * lambda;
    const int ds_n = 10;

    const double ds = ds_base * pow( 2, atoi(argv[3])/(double)ds_n ); //vm["ds"].as<double>()*lambda;
    Eigen::Vector3d dsvec(ds, ds, ds); 
    AIM::Grid grid(dsvec, expansion, *dots);
    const int transit_steps = grid.max_transit_steps(c0, dt) + interp;
    AIM::Expansions::EFIE expansion_func(transit_steps, c0, dt);
    AIM::Normalization::Laplace normalization( prop_constant );
    // AIM::Normalization::Helmholtz normalization( omega / c0, prop_constant );
    auto aim_interaction = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt,
        expansion_func, normalization, omega);

    std::cout << "AIM expansion order: " << expansion << std::endl;
    std::cout << "ds/lambda: " << ds/lambda << std::endl;

    // compute fields
    std::vector<std::vector<std::complex<double>>> fld_obs_dir(
        steps, std::vector<std::complex<double>>(ntrgs));
    std::vector<std::vector<std::complex<double>>> fld_obs_aim(
        steps, std::vector<std::complex<double>>(ntrgs));
    std::vector<std::vector<std::complex<double>>> fld_anlytc(
        steps, std::vector<std::complex<double>>(ntrgs));
  
    clock_t start_time;
    start_time = std::clock();

    double dist, delay;
    for (int i = 0; i < steps; ++i) {

        const InteractionBase::ResultArray dir_array = direct_interaction->evaluate(i);
        // const InteractionBase::ResultArray aim_array = aim_interaction->evaluate(i);

        for (int itrg = 0; itrg < ntrgs; ++itrg){
        // for (int itrg = 0; itrg < ndots; ++itrg){
            for (int isrc = 0; isrc < nsrcs; ++isrc){
                Eigen::Vector3d dr(separation( (*dots)[nsrcs+itrg], (*dots)[isrc] ));
                dist = dr.norm();
                delay = dist / c0;

                if ( i == 0 ) std::cout << dist << std::endl;

                Eigen::Vector3d efld_d0 = efld_d0_source( i * dt, mu, sigsqr, delay );
                Eigen::Vector3d efld_d1 = efld_d1_source( i * dt, mu, sigsqr, delay );
                Eigen::Vector3d efld_d2 = efld_d2_source( i * dt, mu, sigsqr, delay );

                fld_anlytc[i][itrg] += analytic_EFIE_interaction(efld_d0, efld_d1, efld_d2, dr, c0, dist)[1];
                /*std::cout << itrg << " " << isrc << " " << 
                    analytic_EFIE_interaction(efld_d0, efld_d1, efld_d2, dr, c0, dist)[1] * 
                    (*dots)[nsrcs+itrg].dipole().norm() / hbar 
                << std::endl;*/
                // fld_anlytc[i][itrg] += analytic_Helmholtz_interaction(efld_d0, c0, dist)[1];
                
            }
            fld_obs_dir[i][itrg] = dir_array[nsrcs+itrg] * hbar / (*dots)[nsrcs+itrg].dipole().norm();
            // fld_obs_aim[i][itrg] = aim_array[nsrcs+itrg] * hbar / (*dots)[nsrcs+itrg].dipole().norm();
            // std::cout << fld_anlytc[i][itrg] << " " << fld_obs_dir[i][itrg] << " " << fld_obs_aim[i][itrg] << std::endl;
        
           // ** think about how to handle when \vec{d} not parallel to \vec{E} **
        }
    } 

    // output 
    std::cout << "Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s" << std::endl;
   
    std::ofstream outfile1, outfile2, outfile3;
    outfile1.open("outtest/direct_interaction_fld1.dat");
    outfile1 << std::scientific << std::setprecision(15);
    outfile2.open("outtest/aim_interaction_fld1.dat");
    outfile2 << std::scientific << std::setprecision(15);
    outfile3.open("outtest/error.dat", std::ios_base::app);
    outfile3 << std::scientific << std::setprecision(15);

    double normdiff_dir = 0;
    double normdiff_aim = 0;
    double normdiff_aimdir = 0;
    for (int itrg = 0; itrg < ntrgs; ++itrg){
        for (int i = 0; i < steps; ++i){
            outfile1 << real(fld_obs_dir[i][itrg]) << " " << real(fld_anlytc[i][itrg]) << " " << real(fld_obs_dir[i][itrg]) / real(fld_anlytc[i][itrg]) << std::endl; 
            normdiff_dir += pow( real(fld_obs_dir[i][itrg]) - real(fld_anlytc[i][itrg]), 2);

            outfile2 << real(fld_obs_aim[i][itrg]) << " " << real(fld_anlytc[i][itrg]) << " " << real(fld_obs_aim[i][itrg]) / real(fld_anlytc[i][itrg]) << std::endl; 
            
            normdiff_aim += pow( real(fld_obs_aim[i][itrg]) - real(fld_anlytc[i][itrg]), 2);
            normdiff_aimdir += pow( real(fld_obs_aim[i][itrg]) - real(fld_obs_dir[i][itrg]), 2);
        }
    }

    std::cout << "Direct interaction L2 norm diff: " << sqrt(normdiff_dir / (steps*ntrgs))<< std::endl;
    std::cout << "AIM interaction L2 norm diff: " << sqrt(normdiff_aim / (steps*ntrgs)) << std::endl;
    std::cout << "AIM-Direct interaction L2 norm diff: " << sqrt(normdiff_aimdir / (steps*ntrgs)) << std::endl;
    /*outfile3 << expansion << " " << ds/lambda 
        << " " << sqrt(normdiff_dir / (steps*ntrgs)) 
        << " " << sqrt(normdiff_aim / (steps*ntrgs)) 
        << " " << sqrt(normdiff_aimdir / (steps*ntrgs)) << std::endl;
*/
    outfile1.close();
    outfile2.close();    
    outfile3.close();

}
