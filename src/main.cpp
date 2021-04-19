#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
//#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
//#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"
#include "interactions/self_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << setprecision(15) << scientific;
    auto vm = parse_configs(argc, argv);

    string idstr(argv[1]);
    const bool rotating = // atoi(argv[2]);
      config.ref_frame == Configuration::REFERENCE_FRAME::ROT ? true : false;
    const bool rwa = 
      config.rwa == Configuration::RWA::TRUE ? true : false;

    const double dip = 0.002536;
    const double beta = // 1.00e-3;
                        config.mu0 * pow(dip,2) * pow(config.omega,3) 
                          / ( 6.0 * M_PI * config.hbar * config.c0 ); // 1.79e-04;

    cout << "Initializing..." << endl;
    cout << "  Omega: " << config.omega << std::endl;
    cout << "  Beta: " << beta << std::endl;
 
    cout << "  Solve type: "
              << ((config.sim_type == Configuration::SIMULATION_TYPE::FAST) 
                ? "FAST" 
                : "SLOW") << std::endl;
	  cout << "  Interacting particles: "
              << ((config.interacting == Configuration::INTERACTING::TRUE) 
                ? "TRUE" 
                : "FALSE") << std::endl;
	  cout << "  Self-interacting particles: "
              << ((config.self_interacting == Configuration::SELF_INTERACTING::TRUE) 
                ? "TRUE" 
                : "FALSE") << std::endl;
		cout << "  Reference frame: "
              << (rotating 
                ? "ROT" 
                : "FIX") << std::endl;
		cout << "  RWA: "
              << ((config.rwa == Configuration::RWA::TRUE) 
                ? "TRUE" 
                : "FALSE") << std::endl;
	
    cout << "  Timestep: " << config.dt << std::endl;
		cout << "  Simulation time: " << config.total_time << std::endl;
    cout << "  Num timesteps: " << config.num_timesteps << std::endl;
		cout << "  Interp order: " << config.interpolation_order << std::endl;
  	cout << "  Num corrector steps: " << config.num_corrector_steps << std::endl;
    if (config.sim_type == Configuration::SIMULATION_TYPE::FAST) {
			cout << "  AIM grid spacing: " << config.grid_spacing << endl;
			cout << "  AIM expansion order: " << config.expansion_order << endl;
			cout << "  AIM border: " << config.border << endl;
		}

    auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
    qds->resize(config.num_particles);
    auto rhs_funcs = rhs_functions(*qds, config.omega, beta, rotating);
  
    std::shared_ptr<DotVector> obs;
    obs = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
    // obs = make_shared<DotVector>(import_dots("./dots/obss_line.cfg"));
    cout << "  Num sources: " << config.num_particles << std::endl;
    cout << "  Num observers: " << obs->size() << std::endl;	
    
    const bool getflux = 0;
    // Eigen::Vector3d pos = (*obs)[0].position();
    // cout << "  Observer radius: " << pos.norm() << std::endl;

    // == HISTORY ====================================================

    int task_idx = stoi(idstr);
    int window = 22;
    int min_time_to_keep =
        // config.num_timesteps + window;
        max_transit_steps_between_dots(qds, config.c0, config.dt) +
        config.interpolation_order;
    std::cout << "  Min time to keep: " << min_time_to_keep << std::endl;

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        config.num_particles, window, config.num_timesteps, min_time_to_keep, 2, task_idx);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0) );
    // history->initialize_past( qd_path );

    // == INTERACTIONS ===============================================

    const double propagation_constant = config.mu0 / (4 * M_PI * config.hbar);

    auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse.cfg"));
 
    cout << "Setting up interactions..." << endl;
 
    std::clock_t start_time;
    start_time = std::clock();

   	std::shared_ptr<InteractionBase> selfwise;
    std::shared_ptr<InteractionBase> pairwise;

    if (config.sim_type == Configuration::SIMULATION_TYPE::FAST) {
      /*AIM::Grid grid(config.grid_spacing, config.expansion_order, h, *qds); 
      const int transit_steps = grid.max_transit_steps(config.c0, config.dt) + 
                                  config.interpolation_order;

      if (config.ref_frame == Configuration::REFERENCE_FRAME::ROT) {
        Propagation::RotatingEFIE dyadic(config.c0, propagation_constant, config.omega, beta, 0.0);
        Propagation::SelfRotatingEFIE dyadic_self(config.c0, propagation_constant, config.omega, beta);

        selfwise = make_shared<SelfInteraction>(qds, history, dyadic_self,
                                                    config.interpolation_order, config.c0, config.dt, config.omega);
        pairwise = make_shared<AIM::Interaction>(
            qds, history, dyadic, config.grid_spacing, 
            config.interpolation_order, config.expansion_order, config.border,
            config.c0, config.dt, 0, 
						AIM::Expansions::RotatingEFIE(transit_steps, config.c0, config.dt, config.omega),
						AIM::Expansions::Zero(transit_steps),
            AIM::Normalization::Helmholtz(config.omega/config.c0, propagation_constant),
            config.omega);
 
     } else {
        Propagation::EFIE<cmplx> dyadic(config.c0, propagation_constant, beta, 0.0);
        Propagation::SelfEFIE dyadic_self(config.c0, propagation_constant, beta);
       
        selfwise = make_shared<SelfInteraction>(qds, history, dyadic_self,
                                                    config.interpolation_order, config.c0, config.dt);
        pairwise = make_shared<AIM::Interaction>(
            qds, history, dyadic, grid_spacing, 
            config.interpolation_order, config.expansion_order, config.border,
            config.c0, config.dt, 0, 
						AIM::Expansions::EFIE(transit_steps, config.c0, config.dt),
						AIM::Expansions::Zero(transit_steps),
            // AIM::Expansions::EFIE_TimeDeriv2(transit_steps, config.c0, config.dt), // "analytic" expansion function
            // AIM::Expansions::EFIE_Retardation(transit_steps, config.c0), // fdtd expansion function
            AIM::Normalization::Laplace(propagation_constant)
           );
      }*/
    
    } else {
      if (config.ref_frame == Configuration::REFERENCE_FRAME::ROT) {
          Propagation::RotatingEFIE dyadic(config.c0, propagation_constant, config.omega, beta, 0.0);
          Propagation::SelfRotatingEFIE dyadic_self(config.c0, propagation_constant, config.omega, beta);

          selfwise = make_shared<SelfInteraction>(qds, nullptr, history, dyadic_self,
                                                     config.interpolation_order, config.c0, config.dt, config.omega);
          pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                                     config.interpolation_order, config.c0, config.dt, config.omega);
      
      } else {
          Propagation::EFIE<cmplx> dyadic(config.c0, propagation_constant, beta, 0.0);
          Propagation::SelfEFIE dyadic_self(config.c0, propagation_constant, beta);
         
          selfwise = make_shared<SelfInteraction>(qds, nullptr, history, dyadic_self,
                                                      config.interpolation_order, config.c0, config.dt);
          pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                                        config.interpolation_order, config.c0, config.dt); 
      }

    }

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
      make_shared<PulseInteraction>(qds, nullptr, pulse1, 
                                    config.interpolation_order, config.dt, config.hbar, config.omega, 
                                    rotating, rwa) };

    if (config.self_interacting == Configuration::SELF_INTERACTING::TRUE)
      interactions.push_back( selfwise );

    if (config.interacting == Configuration::INTERACTING::TRUE)
      interactions.push_back( pairwise );

    double elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

    cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;

    // == SRC-OBS INTERACTIONS ===============================================

    cout << "Setting up efld interactions..." << endl;
   	std::shared_ptr<InteractionBase> selfwise_fld;
    std::shared_ptr<InteractionBase> pairwise_fld;

    if (config.ref_frame == Configuration::REFERENCE_FRAME::ROT) {
      Propagation::RotatingEFIE dyadic(config.c0, propagation_constant, config.omega, beta, 0.0);
      Propagation::SelfRotatingEFIE dyadic_self(config.c0, propagation_constant, config.omega, beta);

      selfwise_fld = make_shared<SelfInteraction>(qds, obs, history, dyadic_self,
                                                     config.interpolation_order, config.c0, config.dt, config.omega);
      pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                    config.interpolation_order, config.c0, config.dt, config.omega);
    } else {
      Propagation::EFIE<cmplx> dyadic(config.c0, propagation_constant, beta, 0.0);
      Propagation::SelfEFIE dyadic_self(config.c0, propagation_constant, beta);

      selfwise_fld = make_shared<SelfInteraction>(qds, obs, history, dyadic_self,
                                                      config.interpolation_order, config.c0, config.dt);
      pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                      config.interpolation_order, config.c0, config.dt); 
    }

    std::vector<std::shared_ptr<InteractionBase>> efld_interactions{ 
      make_shared<PulseInteraction>(qds, obs, pulse1, 
                                    config.interpolation_order, config.dt, config.hbar, config.omega, 
                                    rotating, rwa) };

    if (config.self_interacting == Configuration::SELF_INTERACTING::TRUE)
      efld_interactions.push_back( selfwise_fld );

    if (config.interacting == Configuration::INTERACTING::TRUE)
      efld_interactions.push_back( pairwise_fld );

    cout << "Setting up bfld interactions..." << endl;
    std::shared_ptr<InteractionBase> pairwise_bfld;

    if (config.ref_frame == Configuration::REFERENCE_FRAME::ROT) {
      Propagation::RotatingMFIE dyadic_mfie(config.c0, propagation_constant, config.omega);

      pairwise_bfld = make_shared<DirectInteraction>(qds, obs, history, dyadic_mfie,
                                                     config.interpolation_order, config.c0, config.dt, config.omega);
    } else {
      Propagation::MFIE<cmplx> dyadic_mfie(config.c0, propagation_constant);

      pairwise_bfld = make_shared<DirectInteraction>(qds, obs, history, dyadic_mfie,
                                                     config.interpolation_order, config.c0, config.dt);
    }

    std::vector<std::shared_ptr<InteractionBase>> bfld_interactions{ 
      make_shared<PulseInteraction>(qds, obs, pulse1, 
                                    config.interpolation_order, config.dt, config.hbar, config.omega, 
                                    rotating, rwa) };

    // no self B-fld
    if (config.interacting == Configuration::INTERACTING::TRUE)
      bfld_interactions.push_back( pairwise_bfld );

    // == INTEGRATOR =================================================

    cout << "Setting up integrator..." << endl;
 
    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            config.hbar, config.mu0, config.c0, config.omega,
            config.dt, config.num_timesteps, 
            history, 
            std::move(interactions), std::move(efld_interactions), std::move(bfld_interactions), 
            std::move(rhs_funcs), obs, 
            task_idx, getflux);

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver_pc(
        config.dt, config.num_corrector_steps, 18, window, 3.15, history, std::move(bloch_rhs));

    start_time = std::clock();

    cout << "Solving P-C..." << endl;
 
    solver_pc.solve();

    elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

    cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;

  return 0;
}
