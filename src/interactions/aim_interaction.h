#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <fftw3.h>
#include <algorithm>
#include <functional>
#include <numeric>

#include <iomanip>

#include "../common.h"
#include "history_interaction.h"

namespace AIM {
  class Grid;
  class AimInteraction;
}

class AIM::Grid {
 public:
  typedef Eigen::Array<int, 3, 2> BoundsArray;

  Grid();
  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector> &, const int);
  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector> &);
  std::vector<DotRange> box_contents_map(
      const std::shared_ptr<DotVector> &) const;

  // Geometry routines (grid <---> space)
  Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &) const;
  size_t coord_to_idx(const Eigen::Vector3i &) const;
  Eigen::Vector3i idx_to_coord(size_t) const;
  Eigen::Vector3d spatial_coord_of_box(const size_t) const;
  std::vector<size_t> expansion_box_indices(const Eigen::Vector3d &,
                                            const int) const;

  Eigen::Array3i dimensions;
  size_t num_boxes;
  double max_diagonal;

  int max_transit_steps(double c, double dt) const
  {
    return static_cast<int>(ceil(max_diagonal / (c * dt)));
  };

  BoundsArray calculate_bounds() const;

 private:
  Eigen::Array3d spacing;
  std::shared_ptr<DotVector> dots;
  int padding;
  BoundsArray bounds;

  void sort_points_on_boxidx() const;
};

class AIM::AimInteraction : public HistoryInteraction {
 public:
  AimInteraction(
      const std::shared_ptr<const DotVector> &,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &,
      const int,
      const double,
      const double,
      const Grid &,
      const int);

  const ResultArray &evaluate(const int);
  void fill_gmatrix_table(SpacetimeArray<double> &) const;
  Eigen::VectorXd q_vector(const Eigen::Vector3d &) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;
  Eigen::VectorXd solve_expansion_system(const Eigen::Vector3d &) const;
  std::vector<Eigen::VectorXd> expansion_table() const;

  // private:
  Grid grid;

  int box_order;

  SpacetimeArray<cmplx> fourier_table;
  std::vector<Eigen::VectorXd> expansions;

  void fill_fourier_table();
};
#endif
