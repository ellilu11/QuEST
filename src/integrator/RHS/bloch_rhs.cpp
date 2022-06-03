#include "bloch_rhs.h"
#include "../../common.h"
#include "../../quantum_dot.h"
#include "../../math_utils.h"
#include "../../interactions/interaction.h"

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
    std::shared_ptr<DotVector> srcs,
    std::shared_ptr<DotVector> obss,
    const int task_idx,
    const bool rotating)
    : Integrator::RHS<matrix_cmplx>(dt, history),
      hbar(hbar), mu0(mu0), eps0(1.0/(mu0*c0*c0)), omega(omega), dt(dt),
      num_timesteps(num_timesteps),
      num_solutions(history->num_particles),
      pulse(std::move(pulse)),
      interactions(std::move(interactions)),
      num_interactions(interactions.size()),
      efld_interactions(std::move(efld_interactions)),
      bfld_interactions(std::move(bfld_interactions)),
      srcs(std::move(srcs)),
      obss(std::move(obss)),
      num_obs(obss->size()),
      rotating(rotating)
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

  const int power_num = pow(2, num_solutions);
  H0 = std::vector<matrix_cmplx>(num_solutions);
  Hint = std::vector<matrix_cmplx>(num_interactions);
  auto Hintp = std::vector<matrix_cmplx>(num_interactions);
  auto Hintm = std::vector<matrix_cmplx>(num_interactions);

  // construct pulse Hamiltonian
  for (int i = 0; i < num_solutions; ++i){
    H0[i] = 
      i == 0 ? sigp+sigm : Id;
    for (int j = 1; j < num_solutions; ++j){
      H0[i] = 
        i == j 
          ? kroneckerProduct(H0[i], sigp+sigm).eval()
          : kroneckerProduct(H0[i], Id).eval();
    }
  }

  // construct interaction Hamiltonian
  for (int pair = 0; pair < num_interactions; ++pair){
    int src, obs;
    std::tie(src, obs) = idx2coord(pair);
    // std::cout << src << " " << obs << std::endl;

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
  }
}

matrix_cmplx Integrator::BlochRHS::rhsfunc( const matrix_cmplx rho,
                                            const ResultArray Efld,
                                            const int step )
{

  const int power_num = std::pow(2,num_solutions);
  matrix_cmplx rhodt = Eigen::MatrixXcd::Zero(power_num,power_num);

  // pulse interactions
  for (int i = 0; i < num_solutions; ++i){
    cmplx rabi_pulse = (*srcs)[i].dipole().dot(
                          (*pulse)((*srcs)[i].position(), step*dt, rotating, 1) ) / hbar; 
    matrix_cmplx Hpulse = H0[i] * rabi_pulse;    
    rhodt += -iu * ( Hpulse * rho.adjoint() - rho * Hpulse.adjoint() ); 
  }

  // src-src interactions
  for (int pair = 0; pair < num_interactions; ++pair){
    int src, obs;
    std::tie(src, obs) = idx2coord(pair);
    Eigen::Vector3cd Efld_pair = Efld[pair];
    cmplx rabi = (*srcs)[obs].dipole().dot( Efld_pair ) / hbar;
    matrix_cmplx Hpair = Hint[pair] * rabi;
    rhodt += -iu * ( Hpair * rho.adjoint() - rho * Hpair.adjoint() ); 
  }

  // damping terms 
  //
  //

  return rhodt;
}

void Integrator::BlochRHS::evaluate(const int step)
{
  /*auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();
 */
  auto projected_efld = interactions[0]->evaluate(step);
  /*std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);
*/
  history->set_value(step, 1) = rhsfunc(
        history->get_value(step, 0),
	      projected_efld, 
        step);
}

void Integrator::BlochRHS::evaluate_present(const int step)
{
  /*auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_present(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();
*/
  auto projected_efld = interactions[0]->evaluate(step);
  /*std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);
  */
  history->set_value(step, 1) = rhsfunc(
        history->get_value(step, 0),
	      projected_efld,
        step);
}

void Integrator::BlochRHS::evaluate_field(const int step)
{
/*  auto nil = InteractionBase::ResultArray::Zero(num_obs, 1).eval();

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
*/
}

