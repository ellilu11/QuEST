#ifndef BLOCH_RHS_H
#define BLOCH_RHS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

#include "../../interactions/interaction.h"
#include "../../math_utils.h"
#include "../../quantum_dot.h"
#include "rhs.h"
namespace Integrator {
  class BlochRHS;
}

class Integrator::BlochRHS : public Integrator::RHS<Eigen::Vector2cd> {
 public:
  BlochRHS(const double, const double, const double, const double,
           const int,
           const std::shared_ptr<History<Eigen::Vector2cd>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<std::shared_ptr<InteractionBase>>,
           std::vector<BlochFunctionType>,
           std::shared_ptr<DotVector>,
           const int, 
           const bool);
  ~BlochRHS();
  void evaluate(const int) const override;
  void evaluate_present(const int) const override;
  void evaluate_field(const int) override;
  std::vector<BlochFunctionType> rhs_functions;

private:
  int num_solutions, num_obs;
  int num_timesteps;
  const double hbar, mu0, eps0;
  std::vector<std::shared_ptr<InteractionBase>> interactions, efld_interactions, bfld_interactions;
  std::shared_ptr<DotVector> obss;
  // std::vector<double> area_elements;

  std::ofstream outfile, fluxfile;
  const bool getflux;
};

#endif
