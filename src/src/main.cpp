#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>
#include <ctime>

#include "configuration.h"
#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    std::clock_t start_time;
    start_time = std::clock();

    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    std::unique_ptr<DurationLogger> durationLogger =
        config.report_time_data ? std::make_unique<DurationLogger>() : nullptr;

    cout << "Initializing..." << endl;
    std::cout << "  Running in "
              << ((config.sim_type == Configuration::SIMULATION_TYPE::FAST)
                      ? "FAST"
                      : "SLOW")
              << " mode" << std::endl;
    std::cout << "  Interacting particles: " 
              << config.interacting
              /*<< ((config.interacting == Configuration::INTERACTING::TRUE)
                      ? "TRUE"
                      : "FALSE") */
              << std::endl;
    std::cout << "  Frame: "
              << ((config.rotating) ? "Rotating" : "Fixed") << std::endl;

    auto qds = make_shared<DotVector>(import_dots(config.qd_path));
    qds->resize(config.num_particles);
    auto rhs_funcs = rhs_functions(*qds, config.omega, config.beta, config.rotating);

    if(durationLogger) {
      durationLogger->log_event("Import dots");
    }

    // == HISTORY ====================================================

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        config.num_particles, 22, config.num_timesteps);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0), 50 );
    // history->initialize_past( config.qd_path );

    auto history_efld = make_shared<Integrator::History<Eigen::Vector2cd>>(
        config.num_particles, 22, config.num_timesteps);
    history_efld->fill(Eigen::Vector2cd::Zero());

    if(durationLogger) {
      durationLogger->log_event("Initialize history");
    }

    // == INTERACTIONS ===============================================

    const double propagation_constant = config.mu0 / (4 * M_PI * config.hbar);

    std::shared_ptr<InteractionBase> pairwise;
    Propagation::RotatingEFIE dyadic(config.c0, propagation_constant, config.omega);

/*    dyadic = config.rotating ? 
        make_shared<Propagation::RotatingEFIE>(config.c0, propagation_constant,
                                              config.omega) :
        make_shared<Propagation::EFIE<cmplx>>(config.c0, propagation_constant);
*/
    if(config.sim_type == Configuration::SIMULATION_TYPE::FAST) {
      AIM::Grid grid(config.grid_spacing, config.expansion_order, *qds);
      const int transit_steps = grid.max_transit_steps(config.c0, config.dt) +
                                config.interpolation_order;

      pairwise = make_shared<AIM::Interaction>(
          qds, history, dyadic, config.grid_spacing,
          config.interpolation_order, config.expansion_order, config.border,
          config.c0, config.dt,
          AIM::Expansions::RotatingEFIE(transit_steps, config.c0, config.dt,
                                        config.omega),
          AIM::Normalization::Helmholtz(config.omega / config.c0,
                                        propagation_constant),
          config.omega);
    } else {
      pairwise = make_shared<DirectInteraction>(qds, history, dyadic,
                                                config.interpolation_order,
                                                config.c0, config.dt);
    }

    auto pulse1 = make_shared<Pulse>(read_pulse_config(config.pulse_path));

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
        make_shared<PulseInteraction>(qds, pulse1, config.hbar, config.dt, config.rotating)};

    if(config.interacting) interactions.push_back( pairwise );

    if(durationLogger) {
      durationLogger->log_event("Initialize interactions");
    }

    // == INTEGRATOR =================================================

    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            config.dt, history, history_efld, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver(
        config.dt, 18, 22, 3.15, history, std::move(bloch_rhs));

    cout << "Solving..." << endl;
    solver.solve(log_level_t::LOG_INFO);

    if(durationLogger) {
      durationLogger->log_event("Full solution");
    }

    // == OUTPUT =====================================================

    cout << "Writing output..." << endl;
   
    ofstream rhofile("rho.dat");
    rhofile << scientific << setprecision(15);

    ofstream fldfile("fld.dat");
    fldfile << scientific << setprecision(15);
   
     for(int t = 0; t < config.num_timesteps; ++t) {
      for(int n = 0; n < config.num_particles; ++n) {
        rhofile << history->array_[n][t][0].transpose() << " ";
        fldfile << history_efld->array_[n][t][0][0] * config.hbar / (*qds)[n].dipole().Eigen::Vector3d::norm() << " ";
      }
      rhofile << "\n";
      fldfile << "\n";
    }

    /*ofstream outfilebin("output.bin");
    for(int t = 0; t < config.num_timesteps; ++t) {
      for(int n = 0; n < config.num_particles; ++n) {
       outfilebin.write((char*)&history->array_[n][t][0].transpose(), 2*sizeof(double));
      }
    }*/

    cout << "Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s";

    if(durationLogger) {
      durationLogger->log_event("Write output");
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
