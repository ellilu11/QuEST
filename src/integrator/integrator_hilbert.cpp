#include <string>
#include "integrator.h"

using namespace std;

constexpr int DIM = 3;
constexpr double T1 = 1000.0;
constexpr double T2 = 2000.0;

template <class soltype>
Integrator::HilbertIntegrator<soltype>::HilbertIntegrator(
    const double dt,
    const double beta,
    const double omega,
    const int interp_order,
    const bool rotating,
    const std::shared_ptr<Integrator::History<soltype>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions
    )
    : num_solutions(history->array_.shape()[0]),
      time_idx_ubound(history->array_.index_bases()[1] +
                      history->array_.shape()[1]),
      dt(dt), beta(beta), omega(omega), interp_order(interp_order), rotating(rotating),
      history(std::move(history)),
      interactions(std::move(interactions)),
      rabi(num_solutions),
      coeffs(interp_order)

{
   for(int solution = 0; solution < num_solutions; ++solution)
      y_prev.segment<DIM>(DIM*solution) = cmplx2real( history->array_[solution][0][0] ); // assign initial conditions

}

template <class soltype>
void Integrator::HilbertIntegrator<soltype>::solve(
    const log_level_t log_level)
{

 for(int step = 1; step < time_idx_ubound; ++step) {
    assert(0 <= step && step < time_idx_ubound);

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    for(int j = 0; j < interp_order; ++j) {
      history->array_[sol_idx][step][0][0] +=
        -iu * ( rabi[sol_idx][step-j] * std:conj( history->array_[sol_idx][step-j][0][1] ) * coeffs[j] - 
                std::conj( rabi[sol_idx][step-j] ) * history->array_[sol_idx][step-j][0][1] * coeffs[j] );
      history->array_[sol_idx][step][0][1] +=
        -iu * ( rabi[sol_idx][step-j] * ( 1.0 - 2.0 * history->array_[sol_idx][step-j][0][0] ) * coeffs[j];

    }
  }

}

