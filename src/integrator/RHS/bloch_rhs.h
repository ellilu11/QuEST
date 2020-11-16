#ifndef BLOCH_RHS_H
#define BLOCH_RHS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

#include "../../interactions/interaction.h"
#include "../../quantum_dot.h"
#include "rhs.h"
namespace Integrator {
  class BlochRHS;
}

class Integrator::BlochRHS : public Integrator::RHS<Eigen::Vector2cd> {
 public:
  BlochRHS(const double, 
           const double,
           const int,
           const std::shared_ptr<History<Eigen::Vector2cd>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<BlochFunctionType>,
           std::shared_ptr<DotVector>,
           const int);
  void evaluate(const int) const override;
  void evaluate_present(const int) const override;
  void evaluate_field(const int) override;

  std::vector<BlochFunctionType> rhs_functions;
private:
  int num_solutions, num_obs;
  int num_timesteps;
  double hbar;
  std::vector<std::shared_ptr<InteractionBase>> interactions, fld_interactions;
  std::shared_ptr<DotVector> obss;

  std::ofstream outfile;
};

#endif
