#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "RHS/rhs.h"
#include "history.h"
#include "logging.h"
#include "math_utils.h"
#include "weights.h"

namespace Integrator {
  template <class soltype>
  class PredictorCorrector;
}

template <class soltype>
class Integrator::PredictorCorrector {
 public:
  PredictorCorrector(const double,
                     const int,
                     const int,
                     const int,
                     const double,
                     const std::shared_ptr<Integrator::History<soltype>>,
                     std::unique_ptr<Integrator::RHS<soltype>>);
  void solve(const log_level_t = log_level_t::LOG_NOTHING) const;

  void output_writer(const int) const;

 private:
  int num_solutions, time_idx_ubound, num_corrector_steps;
  double dt;
  Weights weights;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::unique_ptr<Integrator::RHS<soltype>> rhs;

  void solve_step(const int) const;
  void predictor(const int) const;
  void corrector(const int) const;

  void log_percentage_complete(const int) const;
};

template <class soltype>
Integrator::PredictorCorrector<soltype>::PredictorCorrector(
    const double dt,
    const int num_corrector_steps,
    const int n_lambda,
    const int n_time,
    const double radius,
    const std::shared_ptr<Integrator::History<soltype>> history,
    std::unique_ptr<Integrator::RHS<soltype>> rhs)
    : num_solutions(history->num_particles),
      time_idx_ubound(history->num_timesteps),
      dt(dt),
      num_corrector_steps(num_corrector_steps),
      weights(n_lambda, n_time, radius),
      history(std::move(history)),
      rhs(std::move(rhs))
{
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve(
    const log_level_t log_level) const
{
  int num_logsteps = 100;
  int outstep = 5;

  for(int step = 0; step < time_idx_ubound; ++step) {
    solve_step(step);
    
    //if (step%outstep == 0)
    history->write_step_to_file(step);

    if (step%(time_idx_ubound/num_logsteps) == 0) 
      std::cout << step / (time_idx_ubound/num_logsteps) << std::endl;

    // if(log_level >= log_level_t::LOG_INFO) log_percentage_complete(step);
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve_step(const int step) const
{
  assert(0 <= step && step < time_idx_ubound);

  predictor(step);
  rhs->evaluate(step);

  for(int m = 0; m < num_corrector_steps; ++m) {
    corrector(step);
    rhs->evaluate(step);
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::predictor(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->set_value(sol_idx, step, 0) +=
          history->get_value(sol_idx, start + h, 0) * weights.ps(0, h) +
          history->get_value(sol_idx, start + h, 1) * weights.ps(1, h) * dt;
    }
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::corrector(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    history->set_value(sol_idx, step, 0) =
        weights.future_coef * history->get_value(sol_idx, step, 1) * dt;
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->set_value(sol_idx, step, 0) +=
          history->get_value(sol_idx, start + h, 0) * weights.cs(0, h) +
          history->get_value(sol_idx, start + h, 1) * weights.cs(1, h) * dt;
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

#endif
