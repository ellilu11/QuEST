#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
//#include <iterator>
#include <vector>

#include "quantum_dot.h"

using namespace std;

typedef Eigen::MatrixXcd matrix_cmplx;

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
  

  const int num_particles = 10;
  const int power_num = pow(2,num_particles); 

  Eigen::Matrix2cd init;
  init << 1, 0, 0, 0;

  matrix_cmplx init_tensor = init;
  for (int i = 1; i < num_particles; ++i)
    init_tensor = kroneckerProduct(init_tensor, init).eval();

  Eigen::Matrix2cd sigp, sigm, Id;
  Id = Eigen::Matrix2cd::Identity();
  sigp << 0, 0, 1, 0;
  sigm << 0, 1, 0, 0;

  matrix_cmplx H = Eigen::MatrixXcd::Zero(power_num,power_num);

  const double pulse = 1;

  for (int i = 0; i < num_particles; ++i){
    matrix_cmplx H_i = 
      i == 0 ? sigp+sigm : Id;
    for (int j = 1; j < num_particles; ++j){
      if ( i == j )
        H_i = kroneckerProduct(H_i, sigp+sigm).eval();
      else
        H_i = kroneckerProduct(H_i, Id).eval();
    }
    H += H_i;
  }
  H *= pulse;

  std::cout << H << std::endl;

  matrix_cmplx rhodt = // -iu * ( H * rho );
    -iu * ( H * init_tensor.adjoint() ); // plus damping terms

  std::cout << rhodt << std::endl;
*/
  // === DOTS ==============================================================

  auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
  qds->resize(config.num_particles);
  
  int task_idx = stoi(idstr);
  int window = 22;
  int min_time_to_keep =
      // config.num_timesteps + window;
      max_transit_steps_between_dots(qds, config.c0, config.dt) +
      config.interpolation_order;
  std::cout << "  Min time to keep: " << min_time_to_keep << std::endl;
 
  Eigen::Matrix2cd init;
  init << 1, 0, 0, 0;

  matrix_cmplx init_tensor = init;
  for (int i = 1; i < num_particles; ++i)
    init_tensor = kroneckerProduct(init_tensor, init).eval();

  // === HISTORY ==========================================================

  auto history = make_shared<Integrator::History<matrix_cmplx>>( 
      config.num_particles, window, config.num_timesteps, min_time_to_keep, 2, task_idx);
  history->fill(Eigen::MaxtrixXcd::Zero());
  history->initialize_past(init_tensor);

  // === INTERACTIONS =====================================================
 
  const double propagation_constant = config.mu0 / (4 * M_PI * config.hbar);

  auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse.cfg"));

  cout << "Setting up interactions..." << endl;

  std::clock_t start_time;
  start_time = std::clock();

  std::shared_ptr<InteractionBase> selfwise;
  std::shared_ptr<InteractionBase> pairwise;

  if (config.sim_type == Configuration::SIMULATION_TYPE::FAST) {
   // AIM interaction 
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

  std::vector<std::shared_ptr<InteractionBase>> interactions{};
    /*make_shared<PulseInteraction>(qds, nullptr, pulse1, 
                                  config.interpolation_order, config.dt, config.hbar, config.omega, 
                                  rotating, rwa) };
*/
  if (config.self_interacting == Configuration::SELF_INTERACTING::TRUE)
    interactions.push_back( selfwise );

  if (config.interacting == Configuration::INTERACTING::TRUE)
    interactions.push_back( pairwise );

  double elapsed_time = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;

  cout << "  Elapsed time: " << elapsed_time << "s" << std::endl;

  // === INTEGRATOR =======================================================

  cout << "Setting up integrator..." << endl;

  std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
      std::make_unique<Integrator::BlochRHS>(
          config.hbar, config.mu0, config.c0, config.omega, config.dt, 
          config.num_timesteps, 
          history, 
          pulse,
          std::move(interactions), std::move(efld_interactions), std::move(bfld_interactions), 
          obs, 
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
