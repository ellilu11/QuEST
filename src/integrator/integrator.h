#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <string>

#include "RHS/rhs.h"
#include "history.h"
#include "logging.h"
#include "math_utils.h"
#include "weights.h"
#include "../interactions/interaction.h"
#include "../lagrange_set.h"
#include "../pulse.h"

namespace Integrator {
  template <class soltype>
  class PredictorCorrector;

  template <class soltype>
  class NewtonJacobian;

  template <class soltype>
  class HilbertIntegrator;
}

using namespace std;

// PREDICTOR-CORRECTOR SOLVER

template <class soltype>
class Integrator::PredictorCorrector {
 public:
  PredictorCorrector(const double,
                     const int,
                     const int,
                     const double,
                     const std::shared_ptr<Integrator::History<soltype>>,
                     std::unique_ptr<Integrator::RHS<soltype>>);
  void solve(const log_level_t = log_level_t::LOG_NOTHING) const;
  void solve_step(const int) const;
 
 private:
  int num_solutions, time_idx_ubound;
  double dt;
  Weights weights;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::unique_ptr<Integrator::RHS<soltype>> rhs;

  void predictor(const int) const;
  void corrector(const int) const;

  void log_percentage_complete(const int) const;
};

constexpr int NUM_CORRECTOR_STEPS = 20;
constexpr double EPS = 1e-10;

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
    if ( step%10000 == 0 ) std::cout << step << std::endl;

//    if(log_level >= log_level_t::LOG_INFO) log_percentage_complete(step);
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
  int m = 0;

  for(int m = 0; m < NUM_CORRECTOR_STEPS; ++m) {
  // do{   
    corrector(step);
    rhs->evaluate(step);
    // for( int sol_idx = 0; sol_idx < num_solutions; ++sol_idx )
    //  history_diff[sol_idx] = history->array_[sol_idx][step][0][1] - history_prev[sol_idx]; 
      // if (step == 100000) std::cout << cmplx_norm(history_diff) << std::endl;  

  //  m++;
  } // while ( cmplx_norm(history_diff) > EPS && m < NUM_CORRECTOR_STEPS );
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

// HILBERT INTEGRATOR SOLVER
const double T1 = 10;
const double T2 = 20;

template <class soltype>
class Integrator::HilbertIntegrator {
 public:
  typedef boost::multi_array<cmplx, 2> ResultArray2D;

  HilbertIntegrator(const double, 
                  const double,
                  const double,
                  const int,
                  const std::shared_ptr<Integrator::History<soltype>>,
                  std::vector<std::shared_ptr<InteractionBase>>,
                  const std::shared_ptr<Pulse>
                  );
  void solve(const log_level_t = log_level_t::LOG_NOTHING);
  void solve_step(const int);
 
 private:
  int num_solutions, time_idx_ubound;
  const double dt;
  const double omega, omega0;
  const int interp_order;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::vector<std::shared_ptr<InteractionBase>> interactions;
  std::shared_ptr<Pulse> pulse;
  Interpolation::HilbertLagrangeSet interp;
  ResultArray2D rabis;

};

template <class soltype>
Integrator::HilbertIntegrator<soltype>::HilbertIntegrator(
    const double dt,
    const double omega,
    const double omega0,
    const int interp_order,
    const std::shared_ptr<Integrator::History<soltype>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    const std::shared_ptr<Pulse> pulse
    )
    : num_solutions(history->array_.shape()[0]),
      time_idx_ubound(history->array_.index_bases()[1] +
                      history->array_.shape()[1]),
      dt(dt), omega(omega), omega0(omega0), interp_order(interp_order),
      history(std::move(history)),
      interactions(std::move(interactions)),
      pulse(std::move(pulse)),
      interp(interp_order),
      rabis(boost::extents[interp_order+1][num_solutions]){}
 
template <class soltype>
void Integrator::HilbertIntegrator<soltype>::solve(
    const log_level_t log_level)
{

  for(int step = 0; step < time_idx_ubound - 1; ++step) {
   solve_step(step);
  }
}


template <class soltype>
void Integrator::HilbertIntegrator<soltype>::solve_step(const int step)
{
  assert(0 <= step && step < time_idx_ubound);

  const double t0 = 0.0; // -pulse->delay_();
  interp.evaluate_table_at_x((step+1)*dt, step*dt, t0, dt, omega);
    
  // get Rabi frequency
  auto eval_and_sum =
        [step](const InteractionBase::ResultArray &r,
               const std::shared_ptr<InteractionBase> &interaction) {
            return r + interaction->evaluate(step);
        };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  // store rabis for future use
  for(int j = interp_order; j > 0; --j) 
    rabis[j] = rabis[j-1];
  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx)
    rabis[0][sol_idx] = rabi[sol_idx];

  // Adam-Moulton integration
  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
 
   history->array_[sol_idx][step+1][0][0] = history->array_[sol_idx][step][0][0] 
      + 1.0/dt*interp.evaluations[0][0]
      - 1.0/dt*interp.evaluations[0][1];
   // history->array_[sol_idx][step+1][0][1] = history->array_[sol_idx][step][0][1]; 
     //+ dt * iu * history->array_[sol_idx][step][0][1] * omega0;

    std::cout << std::setprecision(15) << std::scientific << step*dt << ' '
      << std::real(history->array_[sol_idx][step][0][0]) << ' '
      << interp.evaluations[0][0] << ' ' 
      << interp.evaluations[0][1] << std::endl;

   /*for(int j = 0; j <= interp_order; ++j) {
     history->array_[sol_idx][step+1][0][0] +=
      // std::cos( omega / 1000 * (time+t0-j*dt) ) *
      interp.evaluations[0][j];

      cmplx RHO_00 = history->array_[sol_idx][step-j][0][0];
      cmplx RHO_01 = history->array_[sol_idx][step-j][0][1];

      history->array_[sol_idx][step+1][0][0] +=
        -iu * ( rabis[j][sol_idx] * std::conj( RHO_01 ) * interp.evaluations[0][j] - 
                std::conj( rabis[j][sol_idx] ) * RHO_01 * interp.evaluations[0][j] );
        // 0.0 * dt * (RHO_00 - 1.0) / T1;
          
      history->array_[sol_idx][step+1][0][1] +=
        -iu * ( rabis[j][sol_idx] * ( 1.0 - 2.0 * RHO_00 ) * interp.evaluations[0][j] );
        // 0.0 * dt * RHO_01 / T2;
      
    }*/
  } 
}

#endif
