#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>

#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
//#include "integrator/integrator_newton.h"
//#include "interactions/AIM/aim_interaction.h"
// #include "interactions/AIM/grid.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << setprecision(15) << scientific;

    // parameters
    const int num_src = atoi(argv[1]);
    const double tmax = 20000;
    const double dt = 5.0 / pow(10.0, atoi(argv[2]) ); 
                              // rotframe: sigma = 1.0ps -> dt <= 0.52e-1
                              // fixframe: omega = 2278.9013 mev/hbar -> dt <= 1.379e-4
    const int num_timesteps = tmax/dt;
    const int num_corrector_steps = 0;

    const int interpolation_order = 4;
    const bool interacting = atoi(argv[3]);
    const bool rotating = atoi(argv[4]);

    // constants
    const double c0 = 299.792458, hbar = 0.65821193, mu0 = 2.0133545e-04;
    const double omega = 2278.9013;
    double beta = // 0.0;
                  // 1.0e-1 / pow(omega,3);
                  1.79e-4 / pow( omega, 3 );

    const double k0 = omega/c0, lambda = 2.0*M_PI/k0;    

    // AIM
    const double ds = 0.050*lambda;
    Eigen::Vector3d grid_spacing(ds, ds, ds);
    const int expansion_order = 4;
    const int border = 1;

    cout << "Initializing..." << endl;
    std::cout << "  Interacting particles: " 
              << (interacting ? "TRUE" : "FALSE") << std::endl;
    std::cout << "  Frame: "
              << ((rotating) ? "Rotating" : "Fixed") << std::endl;
    std::cout << "  dt: " << dt << std::endl;
    std::cout << "  Num timesteps: " << num_timesteps << std::endl;
    std::cout << "  Num sources: " << num_src << std::endl;
    std::cout << "  Beta: " << beta * pow(omega,3) << std::endl;

    string idstr(argv[5]);
    auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
    qds->resize(num_src);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

    // == HISTORY ====================================================
    int task_idx = atoi(argv[5]);
    int min_time_to_keep =
        max_transit_steps_between_dots(qds, c0, dt) +
        interpolation_order;
    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_src, 22, num_timesteps, min_time_to_keep, 2, task_idx);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0) );
    // history->initialize_past( qd_path );

    // == INTERACTIONS ===============================================

    const double propagation_constant = mu0 / (4 * M_PI * hbar);

    auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse.cfg"));
 
    std::shared_ptr<InteractionBase> selfwise;
    std::shared_ptr<InteractionBase> pairwise;

    if (rotating) {
        Propagation::RotatingEFIE dyadic(c0, propagation_constant, omega, beta, 0.0);
        Propagation::SelfRotatingEFIE dyadic_self(c0, propagation_constant, omega, beta);

        selfwise = make_shared<DirectInteraction>(qds, history, dyadic_self,
                                                    interpolation_order, c0, dt, omega, rotating);
        pairwise = make_shared<DirectInteraction>(qds, history, dyadic,
                                                      interpolation_order, c0, dt, omega, rotating);
    
    } else {
        Propagation::EFIE<cmplx> dyadic(c0, propagation_constant, beta, 0.0);
        Propagation::SelfEFIE dyadic_self(c0, propagation_constant, beta);
       
        selfwise = make_shared<DirectInteraction>(qds, history, dyadic_self,
                                                    interpolation_order, c0, dt, omega, rotating);
        pairwise = make_shared<DirectInteraction>(qds, history, dyadic,
                                                      interpolation_order, c0, dt, omega, rotating); 
    }
 
    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
      make_shared<PulseInteraction>(qds, pulse1, interpolation_order, c0, dt, hbar, rotating),
      selfwise} ;

    if (interacting)
      interactions.push_back( pairwise );

    // == INTEGRATOR =================================================

    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            dt, history, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver_pc(
        dt, num_corrector_steps, 18, 22, 3.15, history, std::move(bloch_rhs));

    cout << "Solving P-C..." << endl;
    double start_time = omp_get_wtime();

    solver_pc.solve();

    double elapsed_time = omp_get_wtime() - start_time;

    cout << "Elapsed time: " << elapsed_time << "s" << std::endl;

    // == FIELD INTERACTIONS ===============================================

/*    *obs = *qds; // examine field at sources only
 
    std::vector<std::shared_ptr<InteractionBase>> interactions_fld{ 
        make_shared<PulseInteraction>(qds, obs, pulse1, interpolation_order, c0, dt, hbar, rotating),
        make_shared<DirectInteraction>(qds, obs, history, dyadic_self,
                                                    interpolation_order, c0, dt, omega, beta, hbar, rotating) };

    if(interacting) {
        std::shared_ptr<InteractionBase> pairwise_fld;

        pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                      interpolation_order, c0, dt, omega, beta, hbar, rotating);
        interactions.push_back( pairwise_fld );
    }
*/

  return 0;
}
