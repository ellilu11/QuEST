#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <fftw3.h>
#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>

#include <iomanip>

#include "../common.h"
#include "../lagrange_set.h"
#include "../math_utils.h"
#include "../quantum_dot.h"
#include "interaction.h"

namespace AIM {
  class Grid;
  class AimInteraction;
}

class AIM::Grid {
 public:
  typedef Eigen::Array<int, 3, 2> BoundsArray;
  typedef std::pair<DotVector::iterator, DotVector::iterator> BoxRange;

  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector> &, const int);
  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector> &);
  Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &) const;
  size_t coord_to_idx(const Eigen::Vector3i &) const;
  Eigen::Vector3i idx_to_coord(size_t) const;
  Eigen::Vector3d spatial_coord_of_box(const size_t) const;
  std::vector<size_t> expansion_box_indices(const Eigen::Vector3d &,
                                            const int) const;

  Eigen::Array3i dimensions;
  size_t num_boxes;
  double max_diagonal;
  std::vector<BoxRange> boxes;

  int max_transit_steps(double c, double dt)
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
  void map_points_to_boxes();
};

class AIM::AimInteraction final : public Interaction {
 public:
  AimInteraction(const std::shared_ptr<DotVector> &,
                 const Eigen::Vector3d &,
                 const int,
                 const double,
                 const double);

  const ResultArray &evaluate(const int);

  // private:
  std::shared_ptr<const DotVector> dots;
  Grid grid;

  int interp_order;
  double c, dt;

  SpacetimeArray<cmplx> fourier_table;

  void fill_fourier_table();
};
#endif
