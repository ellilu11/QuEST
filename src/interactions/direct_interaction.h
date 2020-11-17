#ifndef DIRECT_INTERACTION_H
#define DIRECT_INTERACTION_H

#include "history_interaction.h"
#include "../quantum_obs.h"

class DirectInteraction final : public HistoryInteraction {
 public:
  DirectInteraction(
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
  const ResultArray &evaluate_field(const int, const bool=0) final;
//  boost::multi_array<cmplx, 2> &coefficients() final { return coeffs; }

 private:
  int num_src, num_interactions;
  int num_obs, num_srcobs;

  std::vector<int> floor_delays, floor_delays_fld;
  boost::multi_array<cmplx, 2> coeffs;
  boost::multi_array<Eigen::Vector3cd, 2> fldcoeffs;
  double omega;

  void build_coeff_table(Propagation::Kernel<cmplx> &);
  void build_fldcoeff_table(Propagation::Kernel<cmplx> &);

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
  // static int coord2idxsq(int, int, int); 
};

#endif
