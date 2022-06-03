#ifndef COMMON_H
#define COMMON_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <complex>

typedef std::complex<double> cmplx;

constexpr cmplx iu(0, 1);

typedef Eigen::MatrixXcd matrix_cmplx;

typedef Eigen::Array<Eigen::Vector3cd, Eigen::Dynamic, 1> ResultArray;

template <class T>
using Array2 = boost::multi_array<T, 2>;

#endif
