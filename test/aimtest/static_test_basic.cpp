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

const double c0 = 1.00;
const double //c0 = 299.7924
             mu0 = 2.0133545e-04, hbar = 0.65821193;
const double prop_constant = // 1.00 / hbar;
                             mu0 / (4 * M_PI * hbar);
const double dt = 1;
const int steps = 2048;
const int interp = 5;
const int ndots = 2;
 
const double tmax = steps*dt;
const double f_max = 1 / ( 20.0 * dt );
const double lambda = c0 / f_max;
const double omega = 2 * M_PI * f_max;

const double dip = // 0.00052917721; 
                  1;
int main(int argc, char *argv[]){

    std::ofstream outfile1, outfile2;
    outfile1.open("outtest/direct_interaction_fld1.dat");
    outfile1 << std::scientific << std::setprecision(15);
    outfile2.open("outtest/aim_interaction_fld1.dat");
    outfile2 << std::scientific << std::setprecision(15);
 
   // set up dots and history of dots
    std::shared_ptr<DotVector> dots = std::make_shared<DotVector>();
    std::cout << "Running with " << ndots << " dots" << std::endl;

    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, steps/2, steps); 

    dots->emplace_back(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(dip, 0, 0), Eigen::Vector3d(0, 0, 0));
    dots->emplace_back(Eigen::Vector3d(0.5, 0.6, 0.7), Eigen::Vector3d(dip, 0, 0), Eigen::Vector3d(0, 0, 0));

    //dots->emplace_back(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(dip, 0, 0));
    //dots->emplace_back(Eigen::Vector3d(0.5, 0.6, 0.7), Eigen::Vector3d(dip, 0, 0));

    const double dr = (dots->at(1).position() - dots->at(0).position()).norm();
    const double delay = dr / c0;
    
    const double mu = tmax/2.0;
    const double sigma = tmax/12.0;

    std::function<double(double)> src_fn, obs_fn;
    src_fn = [=](const double t){ return gaussian((t-mu) / sigma); };
    obs_fn = [=](const double t){ return src_fn(t - delay); };

    for (int i = -steps/2; i < steps; ++i) {
        history->array_[0][i][0] = Eigen::Vector2cd(0, src_fn(i * dt)) / (*dots)[0].dipole().norm();
        history->array_[1][i][0] = Eigen::Vector2cd(0, 0);
    }

    // set up propagator
    std::shared_ptr<Propagation::Kernel<cmplx>> propagator = 
        std::make_shared<Propagation::Identity<cmplx>>();

    // calculate analytic and direct solution
    const int expansion = 5;
    const int border = 1;
    const double ds = c0 * dt; //vm["ds"].as<double>()*lambda;
    Eigen::Vector3d dsvec = Eigen::Vector3d::Ones() * ds; 
    AIM::Grid grid(dsvec, expansion, *dots); // sort dots according to AIM scheme, even if not calculating via AIM!
 
    auto direct_interaction = std::make_shared<DirectInteraction>(dots, history, propagator, interp, c0, dt);

    double normdiff_dir = 0;
    for (int i = 0; i < steps; ++i) {
        
        const double fld_anlytc = obs_fn(i * dt) * (*dots)[1].dipole()[0];
        const double fld_dir = direct_interaction->evaluate(i)(1).real();
        normdiff_dir += pow( fld_dir - fld_anlytc, 2);
        outfile1 << fld_dir << " " << fld_anlytc << " " << fld_dir/fld_anlytc << std::endl;
    }

    // set up AIM interaction
    const int transit_steps = grid.max_transit_steps(c0, dt) + interp;
    auto aim_interaction = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt,
        AIM::Expansions::Retardation(transit_steps),
        AIM::Normalization::unit);

    std::cout << "AIM expansion order: " << expansion << std::endl;
    std::cout << "ds/lambda: " << ds/lambda << std::endl;

    double normdiff_aim = 0;
    for (int i = 0; i < steps; ++i) {
      
        const double fld_anlytc_aim = obs_fn(i * dt) * (*dots)[1].dipole()[0];
        const double fld_aim = aim_interaction->evaluate(i)(1).real();
        normdiff_aim += pow( fld_aim - fld_anlytc_aim, 2);
        outfile2 << fld_aim << " " << fld_anlytc_aim << " " << fld_aim/fld_anlytc_aim << std::endl;
    }
    std::cout << "Direct L2 norm diff: " << sqrt(normdiff_dir / (steps * ndots)) << std::endl;
    std::cout << "AIM L2 norm diff: " << sqrt(normdiff_aim / (steps*ndots)) << std::endl;
    
    outfile1.close();
    outfile2.close();    

}
