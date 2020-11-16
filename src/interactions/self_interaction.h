#ifndef SELF_INTERACTION_H
#define SELF_INTERACTION_H

#include "history_interaction.h"
#include "../quantum_obs.h"

class SelfInteraction final : public HistoryInteraction {
 public:
  SelfInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::Kernel<cmplx> &,
      const int,
      const double,
      const double,
      const double = 0);

  const ResultArray &evaluate(const int) final;
  const ResultArray &evaluate_present(const int) final;
  const ResultArray &evaluate_field(const int) final;

 private:
  int num_src, num_obs, num_srcobs;

  boost::multi_array<cmplx, 2> coeffs. fldcoeffs;
  double omega;

  void build_coeff_table(Propagation::Kernel<cmplx> &);

};

#endif
