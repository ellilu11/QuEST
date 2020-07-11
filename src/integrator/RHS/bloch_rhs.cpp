#include "bloch_rhs.h"

Integrator::BlochRHS::BlochRHS(
    const double dt,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<BlochFunctionType> rhs_functions)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history),
      num_solutions(history->array_.shape()[0]),
      interactions(std::move(interactions)),
      rhs_functions(std::move(rhs_functions))
{
}

/*std::vector<cmplx> Integrator::BlochRHS::evaluate_rabi(const int step) const
{
// get rabi frequency
  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  std::vector<cmplx> projected_rabi_vec(num_solutions);
  for (int solution = 0; solution < num_solutions; ++solution)
    projected_rabi_vec[solution] = projected_rabi[solution];

  return projected_rabi_vec;
}*/

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

//  auto projected_rabi = evaluate_rabi(step);

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->array_[solution][step][1] = rhs_functions[solution](
        history->array_[solution][step][0],
	      projected_rabi[solution],
        step); // apply Liouville eq
/*    history->array_[solution][step][1] = rhs_functions[solution](
        history->array_[solution][step][0],
	    (step+1)*dt); // Bessel diff eq */
  }
}


