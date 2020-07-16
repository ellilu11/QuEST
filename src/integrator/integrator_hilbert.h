#ifndef INTEGRATOR_NEWTON_H
#define INTEGRATOR_NEWTON_H

#include <string>

#include "RHS/rhs.h"
#include "history.h"
#include "logging.h"
#include "math_utils.h"
#include "weights.h"

template <class soltype>
class Integrator::HilbertIntegrator {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  NewtonJacobian(const double, 
                  const double,
                  const double,
                  const int,
                  const bool,
                  const std::shared_ptr<Integrator::History<soltype>>,
                  std::vector<std::shared_ptr<InteractionBase>>,
                  Interpolation::HilbertLagrangeSet coeffs,
                  );
  void solve(const log_level_t = log_level_t::LOG_NOTHING);
  void solve_step(const int);
 
 private:
  int num_solutions, time_idx_ubound;
  const double dt;
  const double omega;
  const int interp_order;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::vector<std::shared_ptr<InteractionBase>> interactions;
  Interpolation::HilbertLagrangeSet coeffs;
  ResultArray rabi;


};

#endif
