#include "bloch_rhs.h"

Integrator::BlochRHS::BlochRHS(
    const double hbar, 
    const double dt, 
    const int num_timesteps,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<std::shared_ptr<InteractionBase>> efld_interactions,
    std::vector<std::shared_ptr<InteractionBase>> bfld_interactions,
    std::vector<BlochFunctionType> rhs_functions,
    std::shared_ptr<DotVector> obss,
    const int task_idx)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history),
      hbar(hbar), 
      num_timesteps(num_timesteps),
      num_solutions(history->num_particles),
      interactions(std::move(interactions)),
      efld_interactions(std::move(efld_interactions)),
      bfld_interactions(std::move(bfld_interactions)),
      rhs_functions(std::move(rhs_functions)),
      obss(std::move(obss)),
      num_obs(obss->size())
{
  outfile.open("./out/fld" + std::to_string(task_idx) + ".dat");
  outfile << std::scientific << std::setprecision(15);
}

void Integrator::BlochRHS::evaluate(const int step) const
{
  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->set_value(solution, step, 1) = rhs_functions[solution](
        history->get_value(solution, step, 0),
	      projected_rabi[solution],
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

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->set_value(solution, step, 1) = rhs_functions[solution](
        history->get_value(solution, step, 0),
	      projected_rabi[solution],
        step);
  }
}

void Integrator::BlochRHS::evaluate_field(const int step)
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
  // const int num_ints = fld_interactions->size();

  set_dipole_of_dots( obss, Eigen::Vector3d(hbar, 0, 0) );
  auto efldx = // efld_interactions[1]->evaluate_field(step); 
      std::accumulate(
        efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldx = // bfld_interactions[1]->evaluate_field(step, 1);
      std::accumulate(
        bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
 
  set_dipole_of_dots( obss, Eigen::Vector3d(0, hbar, 0) );
  auto efldy = // efld_interactions[1]->evaluate_field(step); 
      std::accumulate(
        efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldy = // bfld_interactions[1]->evaluate_field(step, 1);
      std::accumulate(
        bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
  
  set_dipole_of_dots( obss, Eigen::Vector3d(0, 0, hbar) );
  auto efldz = // efld_interactions[1]->evaluate_field(step); 
      std::accumulate(
        efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldz = // bfld_interactions[1]->evaluate_field(step, 1);
      std::accumulate(
        bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
 
  for(int solution = 0; solution < num_obs; ++solution) {
    outfile << std::abs(efldx[solution]) << " " 
            << std::abs(efldy[solution]) << " " 
            << std::abs(efldz[solution]) << " "    
            << std::abs(bfldx[solution]) << " " 
            << std::abs(bfldy[solution]) << " " 
            << std::abs(bfldz[solution]) << " ";
 }

  outfile << "\n";
  if ( step > num_timesteps*0.90)
    outfile.flush();
}

