#ifndef INTEGRATOR_NEWTON_H
#define INTEGRATOR_NEWTON_H

#include <string>

#include "RHS/rhs.h"
#include "history.h"
#include "logging.h"
#include "math_utils.h"
#include "weights.h"

template <class soltype>
class Integrator::NewtonJacobian {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  NewtonJacobian(const double, 
                  const double,
                  const double,
                  const int,
                  const bool,
                  const std::shared_ptr<Integrator::History<soltype>>,
                  std::vector<std::shared_ptr<InteractionBase>>
                  );
  void solve(const log_level_t = log_level_t::LOG_NOTHING);
  void solve_step(const int);
 
 private:
  int num_solutions, time_idx_ubound;
  const double dt;
  const double beta;
  const double omega;
  const int interp_order;
  const bool rotating;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::vector<std::shared_ptr<InteractionBase>> interactions;

  Eigen::VectorXd y_vec, y_prev, b_vec, x_vec;
  Eigen::VectorXd rhs_vec, rhs_prev;
  Eigen::MatrixXd rhs_J;
  ResultArray rabi;
  boost::multi_array<cmplx, 2> coeffs; 

  void evaluate(const int);
  void update_rhs(const int);
  void update_J(const int);
  void init_J();

  Eigen::VectorXd cmplx2real(Eigen::Vector2cd);
  Eigen::Vector2cd real2cmplx(Eigen::VectorXd);
  int coord2idx(int, int);

  void log_percentage_complete(const int) const;
};

#endif
