#include "zmatrix.h"

std::vector<double> deriv_0_lagrange_coefficients(
    const int interp_order,
    const double x)
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  std::vector<double> result(interp_order + 1, 1);

  for(int basis_id = -interp_order; basis_id <= 0; ++basis_id) {
    for(int m = -interp_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      result.at(basis_id + interp_order) *= (m - x)/(m - basis_id);
    }
  }

  return result;
}

std::vector<double> deriv_1_lagrange_coefficients(
    const int interp_order,
    const double x)
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  std::vector<double> d0(deriv_0_lagrange_coefficients(interp_order, x));
  std::vector<double> result(interp_order + 1, 0);

  for(int basis_id = -interp_order; basis_id <= 0; ++basis_id) {
    for(int m = -interp_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      result.at(basis_id + interp_order) +=
        d0.at(basis_id + interp_order)/(x - m);
    }
  }

  return result;
}


zmatrix::zmatrix(const std::vector<Eigen::Vector3d> &pts,
    const int history_length) :
  weights(boost::extents[pts.size()][pts.size()][history_length][3])
{}
