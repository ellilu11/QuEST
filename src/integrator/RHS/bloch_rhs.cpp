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
  build_rhsmatrices();
}

Integrator::BlochRHS::~BlochRHS()
{
  outfile.close();
}

void Integrator::BlochRHS::build_rhsmatrices()
{
  Eigen::Matrix2cd Id, sigp, sigm, gd, ex;

  Id = Eigen::Matrix2cd::Identity();
  sigp << 0, 0, 1, 0;
  sigm << 0, 1, 0, 0;
  gd << 1, 0, 0, 0;
  ex << 0, 0, 0, 1;

  const int power_num = pow(2, num_solutions);
  const int numsq = pow(num_solutions, 2);
  H0 = std::vector<matrix_cmplx>(num_solutions);
  D0 = std::vector<matrix_cmplx>(num_solutions);
  Dp = std::vector<matrix_cmplx>(num_solutions);
  Dm = std::vector<matrix_cmplx>(num_solutions);
  Ex = std::vector<matrix_cmplx>(num_solutions);
  Hint = std::vector<matrix_cmplx>(numsq);

  mtrx1 = gd;
  for (int i = 1; i < num_solutions; ++i)
    mtrx1 = kroneckerProduct(mtrx1, gd).eval();  

  mtrx2 = ex;
  for (int i = 1; i < num_solutions; ++i)
    mtrx2 = kroneckerProduct(mtrx2, ex).eval();  

  // build pulse Hamiltonian and damping matrices
  for (int i = 0; i < num_solutions; ++i){
    H0[i] = 
      i == 0 ? sigp+sigm : Id;
    D0[i] = 
      i == 0 ? sigp*sigm : Id;
    Dp[i] = 
      i == 0 ? sigp : Id;
    Dm[i] = 
      i == 0 ? sigm : Id;
    Ex[i] = 
      i == 0 ? ex : Id;
 
    for (int j = 1; j < num_solutions; ++j){
      H0[i] = 
        i == j 
          ? kroneckerProduct(H0[i], sigp+sigm).eval()
          : kroneckerProduct(H0[i], Id).eval();
      D0[i] = 
        i == j 
          ? kroneckerProduct(D0[i], sigp*sigm).eval()
          : kroneckerProduct(D0[i], Id).eval();
      Dp[i] = 
        i == j 
          ? kroneckerProduct(Dp[i], sigp).eval()
          : kroneckerProduct(Dp[i], Id).eval();
      Dm[i] = 
        i == j 
          ? kroneckerProduct(Dm[i], sigm).eval()
          : kroneckerProduct(Dm[i], Id).eval();
      Ex[i] = 
        i == j 
          ? kroneckerProduct(Ex[i], ex).eval()
          : kroneckerProduct(Ex[i], Id).eval();
    }
  }

  // build interaction Hamiltonian matrices
  for (int obs = 0; obs < num_solutions; ++obs){
    for (int src = 0; src < num_solutions; ++src){
      if ( src !=obs ){
        int pair = coord2idxrect(src, obs, num_solutions);
       
        if ( obs == 0 )
          Hint[pair] = sigp;
        else if ( src == 0 )
          Hint[pair] = sigm;
        else 
          Hint[pair] = Id;    
 
        for (int i = 1; i < num_solutions; ++i){
          if ( i == obs )
            Hint[pair] = kroneckerProduct(Hint[pair], sigp).eval();
          else if ( i == src )
            Hint[pair] = kroneckerProduct(Hint[pair], sigm).eval();
          else 
            Hint[pair] = kroneckerProduct(Hint[pair], Id).eval();
        }
     }
    }
  }
}

matrix_cmplx Integrator::BlochRHS::rhsfunc( const matrix_cmplx rho,
                                            const ResultArray Efld,
                                            const int step)
{

  const int power_num = std::pow(2,num_solutions);
  matrix_cmplx rhodt = Eigen::MatrixXcd::Zero(power_num,power_num);

  // pulse interactions
  for (int i = 0; i < num_solutions; ++i){
    cmplx rabi_pulse = (*srcs)[i].dipole().dot(
                          (*pulse)((*srcs)[i].position(), step*dt, rotating, 1) ) / hbar; 
    matrix_cmplx Hpulse = H0[i] * rabi_pulse;    
    rhodt += 0.0 * -iu * ( Hpulse * rho - rho * Hpulse.adjoint() ); 
  }

  // src-src interactions
  matrix_cmplx Hpair = Eigen::MatrixXcd::Zero(power_num,power_num);
 
  for (int obs = 0; obs < num_solutions; ++obs){
    for (int src = 0; src < num_solutions; ++src){
      if ( src !=obs ){
        int pair = coord2idxrect(src, obs, num_solutions);
        if (num_interactions){
          cmplx rabi = (*srcs)[obs].dipole().dot( Efld[pair] ) / hbar;
          Hpair -= Hint[pair] * rabi;

          // trace-preserving term
          rhodt += 2.0 * rabi.imag() * Dm[src] * rho * Dp[obs];
          // if (step == 0)
          //  std::cout << src << " " << obs << " " << std::endl << Hint[pair] << std::endl << std::endl;
          //  std::cout << src << " " << obs << " " << std::endl << Dm[src] << std::endl << Dp[obs] << std::endl;
        }
      }
    }
  }
  rhodt += -iu * ( Hpair * rho - rho * Hpair.adjoint() ); 

  // damping 
  for (int i = 0; i < num_solutions; ++i)
    rhodt -= ( D0[i] * rho + rho * D0[i] - 2 * Dm[i] * rho * Dp[i] )
           / (*srcs)[i].damping().first;

  return rhodt;
}

matrix_cmplx Integrator::BlochRHS::gammafunc( const matrix_cmplx rho,
                                        const ResultArray Efld,
                                        const int step)
{

  double gamma1 = 0;
  double gamma2 = 0;

  /*for (int i = 0; i < num_solutions; ++i){
    cmplx rabi_pulse = (*srcs)[i].dipole().dot(
                          (*pulse)((*srcs)[i].position(), step*dt, rotating, 1) ) / hbar; 
    matrix_cmplx Hpulse = H0[i] * rabi_pulse;    
  }*/

  for (int obs = 0; obs < num_solutions; ++obs){
    for (int src = 0; src < num_solutions; ++src){
      if ( src !=obs ){
        int pair = coord2idxrect(src, obs, num_solutions);
        if (num_interactions){
          cmplx rabi = (*srcs)[obs].dipole().dot( Efld[pair] ) / hbar;
          gamma1 += ( (rho * Hint[pair]).trace() ).real() * 2.0 * rabi.imag();
        }
      }
    }
  }

  for (int i = 0; i < num_solutions; ++i)
    gamma2 += ( (rho * Ex[i]).trace() ).real() / (*srcs)[i].damping().first;

  return gamma1 * mtrx1 + gamma2 * mtrx2;
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

  ResultArray projected_efld;

  if (num_interactions)
    projected_efld = interactions[0]->evaluate(step);
  else
    projected_efld.setZero();
  /*  std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);
*/
  history->set_value(step, 1) = rhsfunc(
        history->get_value(step, 0),
	      projected_efld, 
        step);
  history->set_gamma(step) = gammafunc(
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
  ResultArray projected_efld;

  if (num_interactions)
    projected_efld = interactions[0]->evaluate(step);
  else
    projected_efld.setZero();
  /*  std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);
  */
  history->set_value(step, 1) = rhsfunc(
        history->get_value(step, 0),
	      projected_efld,
        step);
  history->set_gamma(step) = gammafunc(
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

