#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>

#include "../src/configuration.h"
#include "../src/pulse.h"
#include "../src/quantum_dot.h"
#include "../src/quantum_dot.cpp"
#include "../src/integrator/history.h"
#include "../src/integrator/integrator.h"
#include "../src/integrator/RHS/bloch_rhs.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/interacitons/pulse_interaction.h"
#include "../src/interactions/AIM/aim_interaction.h"
#include "../src/math_utils.h"

BOOST_AUTO_TEST_SUITE(interaction_test1a)

struct Universe {
    const double c0, mu0, hbar, propagation_constant;
    const int interp;
    double eps;
    bool rotating;
 
    Universe()
        : c0(299.792458), mu0(2.0133545e-04), hbar(0.65821193), 
          propagation_const(4 * M_PI * hbar),
          interp(5),
          eps(1),
          rotating(0);

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
    const double beta = 0.00;
 
    // initialize dots
    auto dots = std::make_shared<DotVector>(import_dots("dots.cfg"));
    const int ndots = (*dots).size();
    std::cout << "Running with " << ndots << " dots" << std::endl;

    auto rhs_funcs = rhs_functions(*dots, omega, beta, rotating);

    int OBS = ndots/2;

    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, 22, steps);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0), Eigen::Vector2cd(1,0) );

    auto history_efld = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, 22, steps);
    history_efld->fill(Eigen::Vector2cd::Zero());

    // set up propagator
    std::shared_ptr<Propagation::Kernel<cmplx>> propagator;
    propagator = rotating ? 
        make_shared<Propagation::RotatingEFIE>(c0, propagation_constant, omega) :
        make_shared<Propagation::EFIE<cmplx>>(c0, propagation_constant);
 
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

    // set up pulse interaction
    auto pulse = make_shared<Pulse>(read_pulse_config("pulse.cfg"));
    auto pulse_interaction = std::make_shared<PulseInteraction>(dots, pulse, hbar, dt, rotating);

    std::vector<std::shared_ptr<InteractionBase>> interactions{
        pulse_interaction, direct_interaction};

    // integrate to compute fields
    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs = 
        std::make_unique<Integrator::BlochRHS>(
            dt, history, history_efld, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver(
        config.dt, 18, 22, 3.15, history, std::move(bloch_rhs));
 
    // output 
    std::ofstream outfile;
    outfile.open("outtest/interaction_fld2.dat");
    outfile << std::scientific << std::setprecision(15);

    double normdiff = 0;
    for (int i = 0; i < steps; ++i){
        outfile << history_efld->array_[OBS][i][0][0] * hbar / (*dots)[OBS].dipole.norm() << std::endl;
        // normdiff_dir += pow( real(fld_obs_dir[i]) - real(fld_anlytc[i]), 2);

    }

    outfile.close();

    // std::cout << "L2 norm diff: " << sqrt(normdiff_dir) / steps << std::endl;

}
BOOST_AUTO_TEST_SUITE_END()
