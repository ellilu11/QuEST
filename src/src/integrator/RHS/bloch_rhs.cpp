#include "bloch_rhs.h"

Integrator::BlochRHS::BlochRHS(
    const double dt,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history_efld,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<BlochFunctionType> rhs_functions)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history, history_efld),
      num_solutions(history->array_.shape()[0]),
      num_fldsolutions(history_efld->array_.shape()[0]),
      interactions(std::move(interactions)),
      rhs_functions(std::move(rhs_functions))
{
}


void Integrator::BlochRHS::evaluate(const int step) const
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

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->array_[solution][step][1] = rhs_functions[solution](
        history->array_[solution][step][0], projected_rabi[solution]); // apply Liouville eq
  }
}

void Integrator::BlochRHS::evaluatefld(const int step) const
{
  // get (magnitude of) efld
/*  auto eval_and_sum_fld =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluatefld(step);
      };
  auto nilfld = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_flds = std::accumulate(
      interactions.begin(), interactions.end(), nilfld, eval_and_sum_fld);
*/
  auto projected_flds = interactions[interactions.size()-1]->evaluatefld(step); // only get the radiated field

  for(int solution = 0; solution < num_fldsolutions; ++solution)
    history_efld->array_[solution][step][0][0] = projected_flds[solution]; 

}


