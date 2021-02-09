#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"
#include "interactions/self_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << setprecision(15) << scientific;

    // parameters
    const int num_src = atoi(argv[1]);
    const double tmax = 10;
    const double dt = 5.0 / pow(10.0, atoi(argv[2]) ); 
                              // rotframe: sigma = 1.0ps -> dt <= 0.52e-1
                              // fixframe: omega = 2278.9013 mev/hbar -> dt <= 1.379e-4
    const int num_timesteps = tmax/dt;
    const int num_corrector_steps = 0;

    const int interpolation_order = 4;
    const bool interacting = atoi(argv[3]);
    const bool rotating = atoi(argv[4]);
    const bool solve_type = atoi(argv[5]);

    // constants
    const double c0 = 299.792458, hbar = 0.65821193, mu0 = 2.0133545e-04;
    const double omega = 4823.67;
    double beta = // 0.0;
                  // 1.0e-1 / pow(omega,3);
                  1.79e-4 / pow( omega, 3 );

    const double k0 = omega/c0, lambda = 2.0*M_PI/k0;    

    // AIM
    const double ds = 5.0e-2*lambda;
    const double h = 0.5*ds; // FDTD spacing
    Eigen::Vector3d grid_spacing(ds, ds, ds);
    const int expansion_order = 3;
    const int border = 1;

    cout << "Initializing..." << endl;
    cout << "  Interacting particles: " 
              << (interacting ? "TRUE" : "FALSE") << std::endl;
    cout << "  Solve type: "
              << ((solve_type) ? "AIM" : "Direct") << std::endl;
	  cout << "  Frame: "
              << ((rotating) ? "Rotating" : "Fixed") << std::endl;
		
    cout << "  dt: " << dt << std::endl;
		cout << "  Simulation time: " << tmax << std::endl;
    cout << "  Num timesteps: " << num_timesteps << std::endl;
		cout << "  Interp order: " << interpolation_order << std::endl;
    cout << "  Num sources: " << num_src << std::endl;
		if ( solve_type ) {
			cout << "  AIM ds/lambda: " << ds/lambda << endl;
			cout << "  AIM expansion order: " << expansion_order << endl;
			cout << "  AIM border: " << border << endl;
		}
    std::cout << "  Beta: " << beta * pow(omega,3) << std::endl;

    string idstr(argv[6]);
    auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
//    cout << (*qds).size() << std::endl;
    qds->resize(num_src);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

    // == HISTORY ====================================================
    int task_idx = atoi(argv[6]);
    int min_time_to_keep =
        max_transit_steps_between_dots(qds, c0, dt) +
        interpolation_order;
    std::cout << "  Min time to keep: " << min_time_to_keep << std::endl;

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_src, 22, num_timesteps, min_time_to_keep, 2, task_idx);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0) );
    // history->initialize_past( qd_path );

    // == INTERACTIONS ===============================================

    const double propagation_constant = mu0 / (4 * M_PI * hbar);

    auto pulse0 = make_shared<Pulse>(read_pulse_config("pulse_miyajima.cfg"));
    // auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse1.cfg"));
 
    cout << "Setting up interactions..." << endl;
 
    std::clock_t start_time;
    start_time = std::clock();

   	std::shared_ptr<InteractionBase> selfwise;
    std::shared_ptr<InteractionBase> pairwise;

    if (solve_type) {
      AIM::Grid grid(grid_spacing, expansion_order, h, *qds); 
      const int transit_steps = grid.max_transit_steps(c0, dt) + 
                                  interpolation_order;

      if (rotating) {
        Propagation::RotatingEFIE dyadic(c0, propagation_constant, omega, beta, 0.0);
        Propagation::SelfRotatingEFIE dyadic_self(c0, propagation_constant, omega, beta);

        selfwise = make_shared<SelfInteraction>(qds, history, dyadic_self,
                                                    interpolation_order, c0, dt, omega);
        pairwise = make_shared<AIM::Interaction>(
            qds, history, dyadic, grid_spacing, 
            interpolation_order, expansion_order, border,
            c0, dt, 0, 
						AIM::Expansions::RotatingEFIE(transit_steps, c0, dt, omega),
						AIM::Expansions::Zero(transit_steps),
            AIM::Normalization::Helmholtz(omega/c0, propagation_constant),
            omega);
 
     } else {
        Propagation::EFIE<cmplx> dyadic(c0, propagation_constant, beta, 0.0);
        Propagation::SelfEFIE dyadic_self(c0, propagation_constant, beta);
       
        selfwise = make_shared<SelfInteraction>(qds, history, dyadic_self,
                                                    interpolation_order, c0, dt);
        pairwise = make_shared<AIM::Interaction>(
            qds, history, dyadic, grid_spacing, 
            interpolation_order, expansion_order, border,
            c0, dt, 0, 
						AIM::Expansions::EFIE(transit_steps, c0, dt),
						AIM::Expansions::Zero(transit_steps),
            // AIM::Expansions::EFIE_TimeDeriv2(transit_steps, c0, dt), // "analytic" expansion function
            // AIM::Expansions::EFIE_Retardation(transit_steps, c0), // fdtd expansion function
            AIM::Normalization::Laplace(propagation_constant)
           );
      }
    
    } else {
      if (rotating) {
          Propagation::RotatingEFIE dyadic(c0, propagation_constant, omega, beta, 0.0);
          Propagation::SelfRotatingEFIE dyadic_self(c0, propagation_constant, omega, beta);

          selfwise = make_shared<SelfInteraction>(qds, history, dyadic_self,
                                                      interpolation_order, c0, dt, omega);
          pairwise = make_shared<DirectInteraction>(qds, history, dyadic,
                                                        interpolation_order, c0, dt, omega);
      
      } else {
          Propagation::EFIE<cmplx> dyadic(c0, propagation_constant, beta, 0.0);
          Propagation::SelfEFIE dyadic_self(c0, propagation_constant, beta);
         
          selfwise = make_shared<SelfInteraction>(qds, history, dyadic_self,
                                                      interpolation_order, c0, dt);
          pairwise = make_shared<DirectInteraction>(qds, history, dyadic,
                                                        interpolation_order, c0, dt); 
      }

    }

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
      make_shared<PulseInteraction>(qds, pulse0, interpolation_order, c0, dt, hbar, rotating)
	    // make_shared<PulseInteraction>(qds, pulse1, interpolation_order, c0, dt, hbar, rotating),
		  , selfwise} ;

    if (interacting)
      interactions.push_back( pairwise );

    double elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

    cout << "Elapsed time: " << elapsed_time << "s" << std::endl;

    // == INTEGRATOR =================================================

    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            dt, history, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver_pc(
        dt, num_corrector_steps, 18, 22, 3.15, history, std::move(bloch_rhs));

    cout << "Solving P-C..." << endl;
    
    start_time = std::clock();

    solver_pc.solve();

    elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

    cout << "Elapsed time: " << elapsed_time << "s" << std::endl;

  return 0;
}
