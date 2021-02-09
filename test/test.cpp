#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include "../src/configuration.h"
#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"
#include "interactions/self_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << setprecision(15) << scientific;
    auto vm = parse_configs(argc, argv);

    /*cout << "Initializing..." << endl;
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
 	  if ( solve_type ) {
			cout << "  AIM ds/lambda: " << ds/lambda << endl;
			cout << "  AIM expansion order: " << expansion_order << endl;
			cout << "  AIM border: " << border << endl;
		}

    string idstr(argv[1]);
    auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
    qds->resize(num_src);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

    const bool getflux = 0;
   
    std::shared_ptr<DotVector> obs;
    obs = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
    // obs = make_shared<DotVector>(import_dots("./dots/obss_line.cfg"));

    cout << "  Num sources: " << num_src << std::endl;
    cout << "  Num observers: " << obs->size() << std::endl;	
    // Eigen::Vector3d pos = (*obs)[0].position();
    // cout << "  Observer radius: " << pos.norm() << std::endl;

    // == HISTORY ====================================================

    int task_idx = atoi(argv[6]);
    int window = 22;
    int min_time_to_keep =
        // num_timesteps + window;
        max_transit_steps_between_dots(qds, c0, dt) +
        interpolation_order;
    std::cout << "  Min time to keep: " << min_time_to_keep << std::endl;

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_src, window, num_timesteps, min_time_to_keep, 2, task_idx);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0) );
    // history->initialize_past( qd_path );

    // == INTERACTIONS ===============================================

    const double propagation_constant = mu0 / (4 * M_PI * hbar);

    auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse_miyajima.cfg"));
 
    cout << "Setting up interactions..." << endl;
 
    std::clock_t start_time;
    start_time = std::clock();

   	std::shared_ptr<InteractionBase> selfwise;
    std::shared_ptr<InteractionBase> pairwise;

*/
    /*if (solve_type) {
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

          selfwise = make_shared<SelfInteraction>(qds, nullptr, history, dyadic_self,
                                                     interpolation_order, c0, dt, omega);
          pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                                     interpolation_order, c0, dt, omega);
      
      } else {
          Propagation::EFIE<cmplx> dyadic(c0, propagation_constant, beta, 0.0);
          Propagation::SelfEFIE dyadic_self(c0, propagation_constant, beta);
         
          selfwise = make_shared<SelfInteraction>(qds, nullptr, history, dyadic_self,
                                                      interpolation_order, c0, dt);
          pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                                        interpolation_order, c0, dt); 
      }

    }

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
      make_shared<PulseInteraction>(qds, nullptr, pulse1, interpolation_order, dt, hbar, omega, rotating) };
		  // selfwise}; //no selfwise!

    if (interacting)
      interactions.push_back( pairwise );

    double elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

    cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;

    // == SRC-OBS INTERACTIONS ===============================================

    cout << "Setting up efld interactions..." << endl;
   	std::shared_ptr<InteractionBase> selfwise_fld;
    std::shared_ptr<InteractionBase> pairwise_fld;

    if (rotating) {
      Propagation::RotatingEFIE dyadic(c0, propagation_constant, omega, beta, 0.0);
      Propagation::SelfRotatingEFIE dyadic_self(c0, propagation_constant, omega, beta);

      selfwise_fld = make_shared<SelfInteraction>(qds, obs, history, dyadic_self,
                                                     interpolation_order, c0, dt, omega);
      pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                    interpolation_order, c0, dt, omega);
    } else {
      Propagation::EFIE<cmplx> dyadic(c0, propagation_constant, beta, 0.0);
      Propagation::SelfEFIE dyadic_self(c0, propagation_constant, beta);

      selfwise_fld = make_shared<SelfInteraction>(qds, obs, history, dyadic_self,
                                                      interpolation_order, c0, dt);
      pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                    interpolation_order, c0, dt); 
    }

    std::vector<std::shared_ptr<InteractionBase>> efld_interactions{ 
      make_shared<PulseInteraction>(qds, obs, pulse1, interpolation_order, dt, hbar, omega, rotating),
		  selfwise_fld } ; // no selfwise!

    if (interacting)
      efld_interactions.push_back( pairwise_fld );

    cout << "Setting up bfld interactions..." << endl;
    std::shared_ptr<InteractionBase> pairwise_bfld;

    if (rotating) {
      Propagation::RotatingMFIE dyadic_mfie(c0, propagation_constant, omega);

      pairwise_bfld = make_shared<DirectInteraction>(qds, obs, history, dyadic_mfie,
                                                     interpolation_order, c0, dt, omega);
    } else {
      Propagation::MFIE<cmplx> dyadic_mfie(c0, propagation_constant);

      pairwise_bfld = make_shared<DirectInteraction>(qds, obs, history, dyadic_mfie,
                                                     interpolation_order, c0, dt);
    }

    std::vector<std::shared_ptr<InteractionBase>> bfld_interactions{ 
      make_shared<PulseInteraction>(qds, obs, pulse1, interpolation_order, dt, hbar, omega, rotating),
		  } ; // no selfwise!

    if (interacting)
      bfld_interactions.push_back( pairwise_bfld );

    // == INTEGRATOR =================================================

    cout << "Setting up integrator..." << endl;
 
    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            hbar, mu0, c0, omega,
            dt, num_timesteps, 
            history, 
            std::move(interactions), std::move(efld_interactions), std::move(bfld_interactions), 
            std::move(rhs_funcs), obs, 
            task_idx, getflux);

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver_pc(
        dt, num_corrector_steps, 18, window, 3.15, history, std::move(bloch_rhs));

    start_time = std::clock();

    cout << "Solving P-C..." << endl;
 
    solver_pc.solve();

    elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

    cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;
*/
  return 0;
}
