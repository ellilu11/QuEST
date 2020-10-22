#ifndef DIRECT_INTERACTION_H
#define DIRECT_INTERACTION_H

#include "history_interaction.h"
#include "../quantum_obs.h"

class DirectInteraction final : public HistoryInteraction {
 public:
  DirectInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::Kernel<cmplx> &,
      const int,
      const double,
      const double,
      const double = 0);

  const ResultArray &evaluate(const int) final;
  const ResultArray &evaluate_present_field(const int) final;
//  boost::multi_array<cmplx, 2> &coefficients() final { return coeffs; }

 private:
  int num_src, num_interactions;

  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coeffs;
  double omega;

  void build_coeff_table(Propagation::Kernel<cmplx> &);

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
  static int coord2idxsq(int, int, int); 
};

#endif
