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
      const double,
      const double,
      const double,
      const bool);

  const ResultArray &evaluate(const int) final;
  const ResultArray &evaluatefld(const int, const int) final;
  boost::multi_array<cmplx, 2> &coefficients() final { return coeffs; }

 private:
  int num_src, num_obs, num_srcsrc, num_srcobs;

  // Eigen::Array<Eigen::Vector3cd, 2> 
  std::vector<int> floor_delays, floor_delays_srcobs;
  boost::multi_array<cmplx, 2> coeffs;
  boost::multi_array<Eigen::Vector3cd, 3> coeffs_srcobs;
  double omega, beta, hbar;
 
  void build_coeff_table(Propagation::Kernel<cmplx> &);
  void build_coeff_srcobs_table(Propagation::Kernel<cmplx> &);  

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
  static int coord2idxsq(int, int, int); 
};

#endif
