#include <string>
#include "integrator.h"

constexpr int NUM_CORRECTOR_STEPS = 20;
constexpr double EPS = 1e-12;

template <class soltype>
Integrator::PredictorCorrector<soltype>::PredictorCorrector(
    const double dt,
    const int n_lambda,
    const int n_time,
    const double radius,
    const std::shared_ptr<Integrator::History<soltype>> history,
    std::unique_ptr<Integrator::RHS<soltype>> rhs)
    : num_solutions(history->array_.shape()[0]),
      time_idx_ubound(history->array_.index_bases()[1] +
                      history->array_.shape()[1]),
      dt(dt),
      weights(n_lambda, n_time, radius),
      history(std::move(history)),
      rhs(std::move(rhs))
{
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve(
    const log_level_t log_level) const
{
  for(int step = 0; step < time_idx_ubound; ++step) {
    solve_step(step);
    if(log_level >= log_level_t::LOG_INFO) log_percentage_complete(step);
  }
}

double cmplx_norm(std::vector<cmplx> vec)
{
  double norm;
  for( int i = 0; i < vec.size(); ++i )
    norm += std::norm(vec[i]);
  return sqrt(norm);
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve_step(const int step) const
{
  assert(0 <= step && step < time_idx_ubound);

  predictor(step);
  rhs->evaluate(step);

  std::vector<cmplx> history_prev(num_solutions);
  std::vector<cmplx> history_diff(num_solutions);

  for( int sol_idx = 0; sol_idx < num_solutions; ++sol_idx )
    history_prev[sol_idx] = history->array_[sol_idx][step][0][1]; 
  // int m = 0;

  // for(int m = 0; m < NUM_CORRECTOR_STEPS; ++m) {
  do{   
    corrector(step);
    rhs->evaluate(step);
    for( int sol_idx = 0; sol_idx < num_solutions; ++sol_idx )
      history_diff[sol_idx] = history->array_[sol_idx][step][0][1] - history_prev[sol_idx]; 
    // if (step == 100000) std::cout << cmplx_norm(history_diff) << std::endl;  

   //m++;
  } while ( cmplx_norm(history_diff) > EPS );
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::predictor(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->array_[sol_idx][step][0] +=
          history->array_[sol_idx][start + h][0] * weights.ps(0, h) +
          history->array_[sol_idx][start + h][1] * weights.ps(1, h) * dt;
    }
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::corrector(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    history->array_[sol_idx][step][0] =
        weights.future_coef * history->array_[sol_idx][step][1] * dt;
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->array_[sol_idx][step][0] +=
          history->array_[sol_idx][start + h][0] * weights.cs(0, h) +
          history->array_[sol_idx][start + h][1] * weights.cs(1, h) * dt;
    }
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::log_percentage_complete(
    const int step) const
{
  if(step % (time_idx_ubound / 10) == 0) {
    std::cout << "\t" << static_cast<int>(10.0 * step / time_idx_ubound)
              << std::endl;
  }
}
