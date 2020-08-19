#include "bloch_rhs.h"

Integrator::BlochRHS::BlochRHS(
    const double dt,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<BlochFunctionType> rhs_functions)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history),
      num_solutions(history->num_particles),
      interactions(std::move(interactions)),
      rhs_functions(std::move(rhs_functions))
{
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
        step); // apply Liouville eq
  }

  if ( step <= 1 )
    std::cout << projected_rabi[0] << std::endl;

}

