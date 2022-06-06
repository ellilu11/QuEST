#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
//#include <iterator>
#include <vector>

#include "common.h"
#include "configuration.h"
#include "math_utils.h"
#include "quantum_dot.h"
#include "integrator/RHS/bloch_rhs.h"
#include "integrator/RHS/bloch_rhs.cpp"
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
  cout << setprecision(18) << scientific;

  // ================ UNIT TESTING ================================

  /*const double TH0 = 0.5*M_PI;
  const double TH1 = 0.6*M_PI;
  const double PH0 = 0.1*M_PI;
  const double PH1 = 0;

  double rho00_0 = ( 1 - cos(TH0) ) / 2;
  double rho00_1 = ( 1 - cos(TH1) ) / 2;
  std::complex<double> rho01_0(cos(PH0)*sin(TH0)/2, sin(PH0)*sin(TH0)/2);
  std::complex<double> rho01_1(cos(PH1)*sin(TH1)/2, sin(PH1)*sin(TH1)/2);

  history->initialize_past( Eigen::Vector2cd(rho00_0, rho01_0),
                            Eigen::Vector2cd(rho00_1, rho01_1) );
 

  const int num_solutions = 3;
  const int num_interactions = num_solutions * (num_solutions - 1) / 2;

  Eigen::Matrix2cd sigp, sigm, Id;
  Eigen::Vector2cd sp, sm;

  Id = Eigen::Matrix2cd::Identity();
  sigp << 0, 0, 1, 0;
  sigm << 0, 1, 0, 0;
  sp << 1, 0;
  sm << 0, 1;

  std::cout << kroneckerProduct(sp, Id) << std::endl;

  const int power_num = pow(2, num_solutions);
  auto H0 = std::vector<matrix_cmplx>(num_solutions);
  auto Hint = std::vector<matrix_cmplx>(num_interactions);
  auto Hintp = std::vector<matrix_cmplx>(num_interactions);
  auto Hintm = std::vector<matrix_cmplx>(num_interactions);
  matrix_cmplx HintTot = Eigen::MatrixXcd::Zero(power_num,power_num);
 
  // construct pulse Hamiltonian
  for (int i = 0; i < num_solutions; ++i){
    H0[i] = 
      i == 0 ? sigp+sigm : Id;
    for (int j = 1; j < num_solutions; ++j){
      H0[i] = 
        i == j ?
          kroneckerProduct(H0[i], sigp+sigm).eval() :
          kroneckerProduct(H0[i], Id).eval();
    }
  }

  // construct interaction Hamiltonian
  for (int pair = 0; pair < num_interactions; ++pair){
    int src, obs;
    std::tie(src, obs) = idx2coord(pair);
    std::cout << src << " " << obs << std::endl;

    Hintp[pair] = 
      obs == 0 ? sigp : Id;   
    Hintm[pair] = 
      obs == 0 ? sigm : Id;   
 
    for (int i = 1; i < num_solutions; ++i){
      if ( i == obs ){
        Hintp[pair] = kroneckerProduct(Hintp[pair], sigp).eval();
        Hintm[pair] = kroneckerProduct(Hintm[pair], sigm).eval();
      } else if ( i == src ){
        Hintp[pair] = kroneckerProduct(Hintp[pair], sigm).eval();
        Hintm[pair] = kroneckerProduct(Hintm[pair], sigp).eval();
      } else {
        Hintp[pair] = kroneckerProduct(Hintp[pair], Id).eval();
        Hintm[pair] = kroneckerProduct(Hintm[pair], Id).eval();
      }
    }
    Hint[pair] = Hintp[pair] + Hintm[pair];
    HintTot += Hint[pair];
  }
 
  std::cout << HintTot << std::endl;

  matrix_cmplx rhodt = // -iu * ( H * rho );
    -iu * ( HintTot * init_tensor.adjoint() ); // plus damping terms

  std::cout << rhodt << std::endl;
  */
  // === PARAMETERS ======================================================= 
 
  auto vm = parse_configs(argc, argv);

  string idstr(argv[1]);
  // string dip_idstr(argv[2]);
  const bool rotating = 
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

 // === DOTS ==============================================================

  auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
  qds->resize(config.num_particles);
  
  std::shared_ptr<DotVector> obs;
  obs = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
  // obs = make_shared<DotVector>(import_dots("./dots/obssXZ.cfg"));
  cout << "  Num sources: " << config.num_particles << std::endl;
  cout << "  Num observers: " << obs->size() << std::endl;	
 
  int task_idx = stoi(idstr);
  int window = 22;
  int min_time_to_keep =
      // config.num_timesteps + window;
      max_transit_steps_between_dots(qds, config.c0, config.dt) +
      config.interpolation_order;
  std::cout << "  Min time to keep: " << min_time_to_keep << std::endl;
 
  // === HISTORY ==========================================================

  const int power_num = pow(2,config.num_particles); 
  auto history = make_shared<Integrator::History<matrix_cmplx>>( 
      config.num_particles, window, config.num_timesteps, min_time_to_keep, 2, task_idx);
  history->fill(Eigen::MatrixXcd::Zero(power_num,power_num));
 
  Eigen::Matrix2cd init;
  init << 1, 0, 0, 0;

  matrix_cmplx init_tensor = init;
  for (int i = 1; i < config.num_particles; ++i)
    init_tensor = kroneckerProduct(init_tensor, init).eval();

  history->initialize_past(init_tensor);

  // test trace functions
  // Eigen::MatrixXcd init2 = traceOverDims( init_tensor, config.num_particles, 0 );
  // std::cout << init2 << std::endl;

  // === INTERACTIONS =====================================================
 
  const double propagation_constant = config.mu0 / (4 * M_PI);

  auto pulse = make_shared<Pulse>(read_pulse_config("pulse.cfg"));

  cout << "Setting up interactions..." << endl;

  std::clock_t start_time;
  start_time = std::clock();

  std::shared_ptr<InteractionBase> pairwise;

  std::vector<std::shared_ptr<InteractionBase>> interactions{};
  
  Propagation::RotatingEFIE_nodelay dyadic(config.c0, propagation_constant, config.omega, beta, 0.0);
  pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                             config.interpolation_order, config.c0, config.dt, config.omega);

  if (config.interacting == Configuration::INTERACTING::TRUE)
    interactions.push_back( pairwise );

  double elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

  cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;

  // === SRC-OBS INTERACTIONS =============================================

  std::vector<std::shared_ptr<InteractionBase>> efld_interactions{};
  std::vector<std::shared_ptr<InteractionBase>> bfld_interactions{};
 
  // === INTEGRATOR =======================================================

  cout << "Setting up integrator..." << endl;

  std::unique_ptr<Integrator::RHS<matrix_cmplx>> bloch_rhs =
    std::make_unique<Integrator::BlochRHS>(
      config.hbar, config.mu0, config.c0, config.omega, config.dt, 
      config.num_timesteps, 
      history, 
      pulse,
      std::move(interactions), 
      std::move(efld_interactions), 
      std::move(bfld_interactions), 
      qds, 
      obs, 
      task_idx, 
      rotating);

  Integrator::PredictorCorrector<matrix_cmplx> solver_pc(
      config.dt, config.num_corrector_steps, 18, window, 3.15, history, std::move(bloch_rhs));

  start_time = std::clock();

  cout << "Solving P-C..." << endl;

  solver_pc.solve();

  elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

  cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;

  return 0;
}
