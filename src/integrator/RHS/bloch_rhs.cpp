#include "bloch_rhs.h"
#include "../../main.cpp"
#include "../../quantum_dot.h"

Integrator::BlochRHS::BlochRHS(
    const double hbar, 
    const double mu0,
    const double c0,
    const double omega,
    const double dt, 
    const int num_timesteps,
    const std::shared_ptr<Integrator::History<matrix_cmplx>> history,
    const std::shared_ptr<const Pulse> pulse,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<std::shared_ptr<InteractionBase>> efld_interactions,
    std::vector<std::shared_ptr<InteractionBase>> bfld_interactions,
    std::shared_ptr<DotVector> obss,
    const int task_idx,
    const bool getflux)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history),
      hbar(hbar), mu0(mu0), eps0(1.0/(mu0*c0*c0)), omega(omega), dt(dt),
      num_timesteps(num_timesteps),
      num_solutions(history->num_particles),
      pulse(std::move(pulse)),
      interactions(std::move(interactions)),
      efld_interactions(std::move(efld_interactions)),
      bfld_interactions(std::move(bfld_interactions)),
      rhs_functions(std::move(rhs_functions)),
      obss(std::move(obss)),
      num_obs(obss->size()),
      getflux(getflux)
{
  outfile.open("./out/fld" + std::to_string(task_idx) + ".dat");
  outfile << std::scientific << std::setprecision(15);
  build_Hamiltonian();

}

Integrator::BlochRHS::~BlochRHS()
{
  outfile.close();
}

void Integrator::BlochRHS::build_Hamiltonian()
{
  Eigen::Matrix2cd sigp, sigm, Id;

  Id = Eigen::Matrix2cd::Identity();
  sigp << 0, 0, 1, 0;
  sigm << 0, 1, 0, 0;

  H = Eigen::MatrixXcd::Zero(power_num,power_num);

  for (int i = 0; i < num_solutions; ++i){
    matrix_cmplx H_i = 
      i == 0 ? sigp+sigm : Id;
    for (int j = 1; j < num_solutions; ++j){
      if ( i == j )
        H_i = kroneckerProduct(H_i, sigp+sigm).eval();
      else
        H_i = kroneckerProduct(H_i, Id).eval();
    }
    H += H_i;
  }
}

matrix_cmplx Integrator::BlochRHS::bloch_rhs( const matrix_cmplx rho,
                                              const std::vector<cmplx> rabi,
                                              const double step) const
{

  // update the Hamiltonian using pulse and rabi here
  // assume for now that pulse for all dots is at origin, and that RWA applies
  H *= (*pulse)(Eigen::Vector3d::Zero(), step*dt, 0, 0) / 2;

  matrix_cmplx rhodt = -iu * ( H * rho.adjoint() - rho * H.adjoint() ); // plus damping terms

  return rhodt;
}

void Integrator::BlochRHS::evaluate(const int step) const
{
 // cmplx rabi = 0;
  // for (int i = 0; i < interactions.size(); ++i){ 

  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();
 
  auto projected_rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  history->set_value(step, 1) = bloch_rhs(
        history->get_value(step, 0),
	      projected_rabi, 
        step);
  }
}

void Integrator::BlochRHS::evaluate_present(const int step) const
{
  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_present(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  // for(int solution = 0; solution < num_solutions; ++solution) {
  history->set_value(step, 1) = bloch_rhs(
        history->get_value(step, 0),
	      projected_rabi,
        step);
  }
}

/*void Integrator::BlochRHS::evaluate_field(const int step)
{
  auto nil = InteractionBase::ResultArray::Zero(num_obs, 1).eval();

  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_field(step);
      };
  auto eval_and_sum_bfld =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_field(step, 1);
      };

  // int idx = 0; // efld_interactions.size() - 1;

  set_dipole_of_dots( obss, Eigen::Vector3d(hbar, 0, 0) );
  auto efldx = 
    // efld_interactions[idx]->evaluate_field(step); 
    std::accumulate(
      efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldx = 
    // bfld_interactions[idx]->evaluate_field(step, 1);
    std::accumulate(
      bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
 
  set_dipole_of_dots( obss, Eigen::Vector3d(0, hbar, 0) );
  auto efldy = 
    // efld_interactions[idx]->evaluate_field(step); 
    std::accumulate(
      efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldy = 
    // bfld_interactions[idx]->evaluate_field(step, 1);
    std::accumulate(
      bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
  
  set_dipole_of_dots( obss, Eigen::Vector3d(0, 0, hbar) );
  auto efldz = 
    // efld_interactions[idx]->evaluate_field(step); 
    std::accumulate(
      efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldz = 
    // bfld_interactions[idx]->evaluate_field(step, 1);
    std::accumulate(
      bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);

  for(int obs = 0; obs < num_obs; ++obs) {
    Eigen::Vector3cd efld(efldx[obs], efldy[obs], efldz[obs]);
    Eigen::Vector3cd bfld(bfldx[obs], bfldy[obs], bfldz[obs]);

    double efld_abs = 
      // std::sqrt( 
                 std::pow( efldx[obs].real(), 2) +
                 std::pow( efldy[obs].real(), 2) +
                 std::pow( efldz[obs].real(), 2); // );
                  
    outfile << efld_abs << " ";
            // << bfld.real().transpose() << " ";
  }
    
  outfile << std::endl;

}*/

