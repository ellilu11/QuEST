#ifndef BLOCH_RHS_H
#define BLOCH_RHS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

#include "../../common.h"
#include "../../interactions/interaction.h"
#include "../../math_utils.h"
#include "../../quantum_dot.h"
#include "../../pulse.h"
#include "rhs.h"

namespace Integrator {
  class BlochRHS;
}

class Integrator::BlochRHS : public Integrator::RHS<matrix_cmplx> {
 public:
  BlochRHS(const double, const double, const double, const double, const double,
           const int,
           const std::shared_ptr<History<matrix_cmplx>>,
           const std::shared_ptr<const Pulse>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::shared_ptr<DotVector>,
           std::shared_ptr<DotVector>,
           const int, 
           const bool);
  ~BlochRHS();
 
  void build_rhsmatrices();
  matrix_cmplx rhsfunc( const matrix_cmplx,
                        const ResultArray,
                        const int );

  void evaluate(const int) override;
  void evaluate_present(const int) override;
  void evaluate_field(const int) override;
  // std::vector<BlochFunctionType> rhs_functions;

private:
  int num_solutions, num_obs;
  int num_timesteps;
  int num_interactions;
  const double hbar, mu0, eps0, omega, dt;
  std::vector<matrix_cmplx> H0, Hint, Hint2, D0, Dp, Dm;
  std::shared_ptr<const Pulse> pulse;

  std::vector<std::shared_ptr<InteractionBase>> interactions, efld_interactions, bfld_interactions;
  std::shared_ptr<DotVector> srcs, obss;
  const bool rotating;

  std::ofstream outfile, fluxfile;
};

#endif
