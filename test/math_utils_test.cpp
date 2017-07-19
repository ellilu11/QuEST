#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../src/math_utils.h"

BOOST_AUTO_TEST_SUITE(math_utils)

BOOST_AUTO_TEST_CASE(LinSpace)
{
  constexpr double TOLER = 1e-13;
  std::vector<double> ls = linspace(0, 1, 10);

  BOOST_CHECK_SMALL(ls.at(0), TOLER);
  for(int i = 1; i < 10; ++i) {
    BOOST_CHECK_CLOSE(ls.at(i), i/9.0, TOLER);
  }
}

BOOST_AUTO_TEST_CASE(UnitNormal)
{
  constexpr double TOLER = 1e-16;

  const Eigen::Vector3d n_xhat = unit_normal(M_PI_2, 0);
  const Eigen::Vector3d n_yhat = unit_normal(M_PI_2, M_PI_2);
  const Eigen::Vector3d n_zhat = unit_normal(0, 0);

  BOOST_CHECK_CLOSE(n_xhat(0), 1, TOLER);
  BOOST_CHECK_SMALL(n_xhat(1), TOLER);
  BOOST_CHECK_SMALL(n_xhat(2), TOLER);

  BOOST_CHECK_SMALL(n_yhat(0), TOLER);
  BOOST_CHECK_CLOSE(n_yhat(1), 1, TOLER);
  BOOST_CHECK_SMALL(n_yhat(2), TOLER);

  BOOST_CHECK_SMALL(n_zhat(0), TOLER);
  BOOST_CHECK_SMALL(n_zhat(1), TOLER);
  BOOST_CHECK_CLOSE(n_zhat(2), 1, TOLER);
}

BOOST_AUTO_TEST_CASE(GridSequence)
{
  BOOST_CHECK_EQUAL(grid_sequence(0), 0);
  BOOST_CHECK_EQUAL(grid_sequence(1), 1);
  BOOST_CHECK_EQUAL(grid_sequence(2), -1);
  BOOST_CHECK_EQUAL(grid_sequence(3), 2);
  BOOST_CHECK_EQUAL(grid_sequence(4), -2);
  BOOST_CHECK_EQUAL(grid_sequence(5), 3);
  BOOST_CHECK_EQUAL(grid_sequence(6), -3);
  BOOST_CHECK_EQUAL(grid_sequence(7), 4);
  BOOST_CHECK_EQUAL(grid_sequence(8), -4);
  BOOST_CHECK_EQUAL(grid_sequence(9), 5);
}

BOOST_AUTO_TEST_CASE(SplitDouble)
{
  auto x1 = split_double(10.5);
  BOOST_CHECK_EQUAL(x1.first, 10);
  BOOST_CHECK_EQUAL(x1.second, 0.5);

  auto x2 = split_double(3.0);
  BOOST_CHECK_EQUAL(x2.first, 3);
  BOOST_CHECK_EQUAL(x2.second, 0.0);

  auto x3 = split_double(-8.4);
  BOOST_CHECK_EQUAL(x3.first, -8);
  BOOST_CHECK_CLOSE(x3.second, -0.4, 1e-13);
}

BOOST_AUTO_TEST_SUITE_END()
