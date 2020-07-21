#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

namespace Interpolation {
  class UniformLagrangeSet;
  class HilbertLagrangeSet;
  class DerivFive;
  typedef boost::multi_array<double, 2> InterpolationTable;
}

class Interpolation::UniformLagrangeSet {
 public:
  InterpolationTable evaluations;
  UniformLagrangeSet(const int);
  UniformLagrangeSet(const double, const int, const double = 1);

  void evaluate_derivative_table_at_x(const double, const double = 1);
  int order() const { return order_; }
 private:
 const int NUM_DERIVATIVES = 4;
 int order_;
};

class Interpolation::HilbertLagrangeSet {
 public:
  InterpolationTable evaluations;
  HilbertLagrangeSet(const int);

  void evaluate_table_at_x(const double, const double, const double, const double, const double);
  int order() const { return order_; }
 private:
  const int NUM_DERIVATIVES = 2;
  int order_;
};

class Interpolation::DerivFive {
 public:
  InterpolationTable evaluations;

  DerivFive(const double = 1);

  void evaluate_derivative_table_at_x(const double);

  double d0(const int, const double) const;
  double d1(const int, const double) const;
  double d2(const int, const double) const;

  constexpr static int order() { return 5; }
 private:
 const int NUM_DERIVATIVES = 2;
 double dt_;
};

#endif
