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
#include "../src/math_utils.h"

// namespace po = boost::program_options;
 
const double c0 = 1.00;
const double //c0 = 299.792458, 
             mu0 = 2.0133545e-04, hbar = 0.65821193;
const double prop_constant = // 1.00 / hbar;
                             mu0 / (4 * M_PI * hbar);
const int interp = 5;
const int steps = 2048;
const double dt = // 0.1;
                  1;
const double tmax = steps*dt;
const double f_max = 1 / ( 20.0 * dt );
const double lambda = c0 / f_max;
const double omega = 2 * M_PI * f_max;

const double mu = tmax/2.0;
const double sig = // 6 / ( 2 * M_PI * f_max );
                   tmax/12.0;
const double sigsqr = sig * sig;
// const double fmax = 6.0/(2.0*M_PI*sqrt(sigsqr));

Eigen::Vector2cd source(double t, double mu, double sigsqr){
    return Eigen::Vector2cd(0, exp(-std::pow(t - mu, 2) / (2.0 * sigsqr)));
}

Eigen::Vector3d efld_d0_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0, 0); 
}

Eigen::Vector3d efld_d1_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(-(t - mu - delay) / sigsqr * exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0, 0); 
}

Eigen::Vector3d efld_d2_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d((std::pow(t - mu - delay, 2) - sigsqr) / pow(sigsqr, 2) *
        exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0, 0); 
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

Eigen::Vector3d analytic_Identity_interaction(Eigen::Vector3d &efld_d0){
    return hbar * Eigen::Matrix3d::Identity() * efld_d0;
}

Eigen::Vector3d analytic_Laplace_interaction(Eigen::Vector3d &efld_d0, double dist){
    return hbar * prop_constant * Eigen::Matrix3d::Identity() * efld_d0 / dist;
}

Eigen::Vector3cd analytic_Helmholtz_interaction(Eigen::Vector3d &efld_d0, double c0, double dist){
    const std::complex<double> iu(0, 1);
    return prop_constant * hbar * Eigen::Matrix3cd::Identity() * 
        efld_d0 * std::exp( -iu * omega / c0 * dist) / dist;
}

std::vector<std::complex<double>>
    analytic_evaluate( std::shared_ptr<DotVector> dots, int i ){

    int ndots = (*dots).size();
    std::vector<std::complex<double>> fld_anlytc(ndots);    
    double dist, delay;
 
    for (int itrg = 0; itrg < ndots; itrg++) {
        for (int isrc = 0; isrc < ndots; ++isrc){
            if ( itrg != isrc ) {
                Eigen::Vector3d dr(separation( (*dots)[itrg], (*dots)[isrc] ));
                dist = dr.norm();
                delay = dist / c0;

                Eigen::Vector3d efld_d0 = efld_d0_source( i * dt, mu, sigsqr, delay );
                Eigen::Vector3d efld_d1 = efld_d1_source( i * dt, mu, sigsqr, delay );
                Eigen::Vector3d efld_d2 = efld_d2_source( i * dt, mu, sigsqr, delay );

                // fld_anlytc[itrg] += analytic_Identity_interaction(efld_d0)
                //                       .dot((*dots)[itrg].dipole()) / hbar;
                //fld_anlytc[itrg] += analytic_Laplace_interaction(efld_d0, dist)
                //                         .dot((*dots)[itrg].dipole()) / hbar;
                fld_anlytc[itrg] += analytic_EFIE_interaction(efld_d0, efld_d1, efld_d2, dr, c0, dist)
                                       .dot((*dots)[itrg].dipole()) / hbar;
                // std::cout << fld_anlytc[itrg] << std::endl;
            } 
        }
    }
    return fld_anlytc;
}



int main(int argc, char *argv[]){

   // set up dots and history of dots
    auto dots = std::make_shared<DotVector>(import_dots("dots.cfg"));
    const int ndots = (*dots).size();
    std::cout << "Running with " << ndots << " dots" << std::endl;

    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, steps/2, steps); 
    history->fill(Eigen::Vector2cd::Zero());

    // for (int n = 0; n < ndots; ++n)
        for (int i = -steps/2; i < steps; ++i) 
            // if ( n != OBS )
            history->array_[0][i][0] = source(i * dt, mu, sigsqr) / (*dots)[0].dipole().norm();

    // set up propagator
    std::shared_ptr<Propagation::Kernel<cmplx>> propagator = 
        std::make_shared<Propagation::EFIE<cmplx>>(c0, prop_constant);
        // std::make_shared<Propagation::Identity<cmplx>>();
        // std::make_shared<Propagation::Laplace<cmplx>>(prop_constant);
        // std::make_shared<Propagation::Helmholtz>(prop_constant, omega / c0 );
    // Propagation::Laplace<cmplx> propagator_aim(prop_constant);
    // Propagation::EFIE<cmplx> propagator_aim(c0, prop_constant);

    // calculate analytic and direct solution
    const int expansion = 5; //vm["expansion"].as<int>();
    const int border = 1;
    const double ds = c0 * dt; //vm["ds"].as<double>()*lambda;
    Eigen::Vector3d dsvec = Eigen::Vector3d::Ones() * ds; 
    AIM::Grid grid(dsvec, expansion, *dots); // sort dots according to AIM scheme, even if not calculating via AIM!
 
    auto direct_interaction = std::make_shared<DirectInteraction>(dots, history, propagator, interp, c0, dt);

    std::vector<std::vector<std::complex<double>>> 
                                     fld_anlytc( steps, std::vector<std::complex<double>>(ndots) );
    std::vector<std::vector<double>> fld_dir( steps, std::vector<double>(ndots) );
 
    double dist, delay;
    for (int i = 0; i < steps; ++i) {
        
        fld_anlytc[i] = analytic_evaluate( dots, i );
        const InteractionBase::ResultArray array_dir = direct_interaction->evaluate(i);
        for (int itrg = 0; itrg < ndots; itrg++)
            fld_dir[i][itrg] = array_dir[itrg].real(); // * hbar / (*dots)[itrg].dipole().norm();
    }
 
    // set up AIM interaction
    const int transit_steps = grid.max_transit_steps(c0, dt) + interp;
    // AIM::Expansions::TimeDerivative expansion_func(transit_steps, dt);
   // AIM::Normalization::Helmholtz normalization( omega / c0, prop_constant );
    auto aim_interaction = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt,
        AIM::Expansions::Retardation(transit_steps),
        AIM::Normalization::unit,
        // AIM::Normalization::Laplace normalization( prop_constant );
        omega);

    std::cout << "AIM expansion order: " << expansion << std::endl;
    std::cout << "ds/lambda: " << ds/lambda << std::endl;

    // calculate AIM solution
    std::vector<std::vector<double>> fld_aim( steps, std::vector<double>(ndots) );
 
    clock_t start_time;
    start_time = std::clock();

    for (int i = 0; i < steps; ++i) {
       
        const InteractionBase::ResultArray array_aim = aim_interaction->evaluate(i); 
        for (int itrg = 0; itrg < ndots; itrg++)
            fld_aim[i][itrg] = array_aim[itrg].real(); // * hbar / (*dots)[itrg].dipole().norm();
    }

    // output
    clock_t elapsed_time = (std::clock() - start_time) / (double) CLOCKS_PER_SEC; 
    std::cout << "Elapsed time: " << elapsed_time << "s" << std::endl;
   
    std::ofstream outfile1, outfile2, outfile3;
    outfile1.open("outtest/direct_interaction_fld1.dat");
    outfile1 << std::scientific << std::setprecision(15);
    outfile2.open("outtest/aim_interaction_fld1.dat");
    outfile2 << std::scientific << std::setprecision(15);
    
    double fld_dir1, fld_aim1, fld_anlytc1;
    double normdiff_dir = 0;
    double normdiff_aim = 0;
    double normdiff_aimdir = 0;

    for (int itrg = 0; itrg < ndots; ++itrg){
        for (int i = 0; i < steps; ++i){
            fld_dir1 = fld_dir[i][itrg];
            fld_aim1 = fld_aim[i][itrg];
            fld_anlytc1 = real(fld_anlytc[i][itrg]);

            outfile1 << i*dt << " " << fld_dir1 << " " << fld_anlytc1 << " " << fld_dir1 / fld_anlytc1 << std::endl; 
            normdiff_dir += pow( fld_dir1 - fld_anlytc1, 2);

            outfile2 << i*dt << " " << fld_aim1 << " " << fld_anlytc1 << " " << fld_aim1 / fld_anlytc1 << std::endl; 
            normdiff_aim += pow( fld_aim1 - fld_anlytc1, 2);
            normdiff_aimdir += pow( fld_aim1 - fld_dir1, 2);

        }
    }

    std::cout << "Direct interaction L2 norm diff: " << sqrt(normdiff_dir / (steps*ndots))<< std::endl;
    std::cout << "AIM interaction L2 norm diff: " << sqrt(normdiff_aim / (steps*ndots)) << std::endl;
    std::cout << "AIM-Direct interaction L2 norm diff: " << sqrt(normdiff_aimdir / (steps*ndots)) << std::endl;
    
    outfile1.close();
    outfile2.close();    

}
