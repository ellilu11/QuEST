#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <fftw3.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

#include <iomanip>

#include "../../common.h"
#include "../history_interaction.h"
#include "grid.h"

namespace AIM {
  class AimInteraction;
}

class AIM::AimInteraction final : public HistoryInteraction {
 public:
  struct Expansion {
    size_t index;
    double weight;
  };

  AimInteraction(
      const std::shared_ptr<const DotVector> &,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &,
      const int,
      const double,
      const double,
      const Grid &,
      const int);

  ~AimInteraction();

  const ResultArray &evaluate(const int) final;
  void fill_gmatrix_table(SpacetimeVector<cmplx> &) const;
  Eigen::VectorXd q_vector(const Eigen::Vector3d &) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;
  Eigen::VectorXd solve_expansion_system(const Eigen::Vector3d &) const;

  Array<Expansion> expansions() const;

 //private:
  Grid grid;
  int box_order, max_transit_steps;

  Array<Expansion> expansion_table;
  SpacetimeVector<cmplx> fourier_table, source_table, obs_table;
  fftw_plan vector_forward_plan, vector_backward_plan;

  void fill_fourier_table();
  std::pair<fftw_plan, fftw_plan> vector_fft_plans();
  void fill_source_table(const int);
  Eigen::VectorXcd fast_multiply(const int, const int) const;
};

#endif
