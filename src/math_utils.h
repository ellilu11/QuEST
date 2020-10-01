#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <vector>
#include "common.h"

std::vector<double> linspace(const double,
                             const double,
                             const size_t,
                             double* const = nullptr);
Eigen::Vector3d unit_normal(const double, const double);
double gaussian(const double);
//double gaussian(const double, const double, const double, const double, const int);
double skew_gaussian(const double, const double);

int grid_sequence(const size_t);
Eigen::Vector3i idx_to_coord(const size_t, int);
Eigen::Vector3i idx_to_delta(const size_t, int);

//size_t coord_to_idx(const Eigen::Vector3i &, int);
//size_t inv_grid_sequence(const int);
//size_t delta_to_idx(const Eigen::Vector3i &, int);

std::pair<int, double> split_double(const double);
double falling_factorial(const double, int);
double spherical_bessel(double, int);
#endif
