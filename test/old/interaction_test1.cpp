#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>

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

BOOST_AUTO_TEST_SUITE(interaction_test1a)

struct Universe {
    const double c0, mu0, hbar;
    const int interp;
    double eps;
    // bool rotating;
 
    std::shared_ptr<Propagation::Kernel<cmplx>> propagator;

    Universe()
        : c0(299.792458), mu0(2.0133545e-04), hbar(0.65821193), 
          interp(5),
          eps(1),
          propagator( std::make_shared
                <Propagation::EFIE<cmplx>>(c0, mu0 / ( 4 * M_PI * hbar))
            // <Propagation::Identity<cmplx>>()
            ){};

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

    Eigen::Vector3d analytic_interaction(Eigen::Vector3d &efld_d0,
                                         Eigen::Vector3d &efld_d1,
                                         Eigen::Vector3d &efld_d2,
                                         Eigen::Vector3d &dr,
                                         double c0, double dist){
        Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
        Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
        Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3 * rr;

    return - pow(c0, 2) * mu0 / ( 4 * M_PI ) *
           (i3rr * efld_d0 / std::pow(dist, 3) +
            i3rr * efld_d1 / (c0 * std::pow(dist, 2)) +
            irr * efld_d2 / (std::pow(c0, 2) * dist));
    }

};


BOOST_FIXTURE_TEST_CASE(interaction, Universe)
{
    const double tmax = 10.0;
    const double mu = tmax/2.0;
    const double sigsqr = 1.0;
    const double fmax = 6.0/(2.0*M_PI*sqrt(sigsqr));
    const double dt = 1/(20.0*fmax);
    const int steps = floor(tmax/dt);

    const double omega = 2278.9013;
    const double lambda = 2 * M_PI * c0 / omega;
 
    // set up dots and history of dots
    auto dots = std::make_shared<DotVector>(import_dots("dots.cfg"));
    const int ndots = (*dots).size();
    std::cout << "Running with " << ndots << " dots" << std::endl;

    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, 22, steps);
    history->fill(Eigen::Vector2cd::Zero());

    int OBS = ndots/2;

    for (int n = 0; n < ndots; ++n){
        if ( n != OBS ){
            for (int i = -22; i < steps; ++i) 
                history->array_[n][i][0] = source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm();
        }
    }    

    // set up direct interaction
    auto direct_interaction = std::make_shared<DirectInteraction>(dots, history, propagator, interp, c0, dt);

    // set up AIM grid and interaction
    const int expansion = 3;
    const int border = 1;
    const double ds = c0 * dt;
    Eigen::Vector3d dsvec(ds, ds, ds); 
    AIM::Grid grid(dsvec, expansion, *dots);
    const int transit_steps = grid.max_transit_steps(c0, dt) + interp;
    AIM::Expansions::EFIE expansion_func(transit_steps, c0, dt);
    AIM::Normalization::Laplace normalization( mu0 / ( 4 * M_PI * hbar));
    auto aim_interaction = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt,
        expansion_func, normalization, omega);

    // compute fields
    std::vector<std::complex<double>> fld_obs_dir(steps);
    std::vector<std::complex<double>> fld_obs_aim(steps);
    std::vector<std::complex<double>> fld_anlytc(steps);
     
    double dist, delay;
    for (int i = 1; i < steps; ++i) {
        for (int n = 0; n < ndots; ++n){
            if ( n != OBS ){
                Eigen::Vector3d dr(separation( (*dots)[n], (*dots)[OBS] ));
                dist = dr.norm();
                delay = dist / c0;

                Eigen::Vector3d efld_d0 = efld_d0_source( i * dt, mu, sigsqr, delay );
                Eigen::Vector3d efld_d1 = efld_d1_source( i * dt, mu, sigsqr, delay );
                Eigen::Vector3d efld_d2 = efld_d2_source( i * dt, mu, sigsqr, delay );

                fld_anlytc[i] += analytic_interaction(efld_d0, efld_d1, efld_d2, dr, c0, dist)[1];
            }
        }

        fld_obs_dir[i] = direct_interaction->evaluate(i)[OBS] * hbar / (*dots)[OBS].dipole().norm();
        fld_obs_aim[i] = aim_interaction->evaluate(i)[OBS] * hbar / (*dots)[OBS].dipole().norm();
    } 


    // output 
    std::ofstream outfile1, outfile2;
    outfile1.open("outtest/direct_interaction_fld1a.dat");
    outfile1 << std::scientific << std::setprecision(15);
    outfile2.open("outtest/aim_interaction_fld1a.dat");
    outfile2 << std::scientific << std::setprecision(15);

    double normdiff_dir = 0;
    double normdiff_aim = 0;
    for (int i = 0; i < steps; ++i){
        outfile1 << real(fld_obs_dir[i]) << " " << real(fld_anlytc[i]) << " " << real(fld_obs_dir[i]) / real(fld_anlytc[i]) << std::endl; 
        normdiff_dir += pow( real(fld_obs_dir[i]) - real(fld_anlytc[i]), 2);

        outfile2 << real(fld_obs_aim[i]) << " " << real(fld_anlytc[i]) << " " << real(fld_obs_aim[i]) / real(fld_anlytc[i]) << std::endl; 
        normdiff_aim += pow( real(fld_obs_aim[i]) - real(fld_anlytc[i]), 2);

    }

    outfile1.close();
    outfile2.close();    

    std::cout << "Direct interaction L2 norm diff: " << sqrt(normdiff_dir / steps) << std::endl;
    std::cout << "AIM interaction L2 norm diff: " << sqrt(normdiff_aim / steps) << std::endl;

}
BOOST_AUTO_TEST_SUITE_END()
