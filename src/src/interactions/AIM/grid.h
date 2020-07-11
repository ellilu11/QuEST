#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>
#include "math_utils.h"
#include "quantum_dot.h"

namespace AIM {
  class Grid;
}

class AIM::Grid {
 public:
  using BoundsArray = Eigen::Array<int, 3, 2>;
  using ipair_t = std::pair<int, int>;

  Grid(const Eigen::Array3d &,
       const int,
       const Eigen::Array3i &,
       const Eigen::Vector3i & = Eigen::Vector3i::Zero());
  Grid(const Eigen::Array3d &, const int, DotVector &);

  BoundsArray calculate_bounds(const DotVector &) const;
  std::array<int, 4> circulant_shape(const double,
                                     const double,
                                     const int = 0) const;

  std::vector<size_t> expansion_indices(const int) const;
  std::vector<size_t> expansion_indices(const Eigen::Vector3d &pos) const
  {
    return expansion_indices(associated_grid_index(pos));
  };

  int max_transit_steps(double c, double dt) const
  {
    double max_diagonal = (dimensions.cast<double>() * spacing).matrix().norm();
    return static_cast<int>(ceil(max_diagonal / (c * dt)));
  };

  std::vector<const_DotRange> box_contents_map(const DotVector &) const;
  std::vector<ipair_t> nearfield_box_pairs(const int, const DotVector &) const;
  std::vector<ipair_t> nearfield_point_pairs(const int,
                                             const DotVector &) const;

  inline auto size() const { return num_gridpoints; }
  inline const auto &shape() const { return dimensions; }
  // == Geometry routines (grid <---> space) ==================================

  inline Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &coord) const
  {
    return (coord.array() / spacing).cast<int>();
  }

  inline size_t associated_grid_index(const Eigen::Vector3d &coord) const
  {
    Eigen::Vector3i grid_coord = grid_coordinate(coord);
    return coord_to_idx(grid_coord - bounds.col(0).matrix());
  }

  inline size_t coord_to_idx(const Eigen::Vector3i &coord) const
  {
    return coord(2) + dimensions(2) * (coord(1) + dimensions(1) * coord(0));
  }

  inline Eigen::Vector3i idx_to_coord(size_t idx) const
  {
    const int nynz = dimensions(1) * dimensions(2);
    const int x = idx / nynz;
    idx -= x * nynz;
    const int y = idx / dimensions(2);
    const int z = idx % dimensions(2);

    return Eigen::Vector3i(x, y, z);
  }

  inline Eigen::Vector3d spatial_coord_of_box(const size_t box_id) const
  {
    Eigen::Vector3i dr = (idx_to_coord(box_id) + bounds.col(0).matrix());
    return dr.array().cast<double>() * spacing;
  }
 
  // size_t num_gridpoints;

 private:
  Eigen::Array3d spacing;
  int expansion_order;
  BoundsArray bounds;
  Eigen::Array3i dimensions;
  size_t num_gridpoints;
  void sort_points_on_boxidx(DotVector &) const;
};

#endif
