#include "math_utils.h"

int wrapmod(const int n, const int M) { return ((n % M) + M) % M; }
std::vector<double> linspace(const double low,
                             const double high,
                             const size_t n,
                             double* const step /*= nullptr*/)
{
  std::vector<double> xs(n);
  const double dx = (high - low) / (n - 1);
  if(step) *step = dx;

  for(size_t i = 0; i < n; ++i) {
    xs[i] = low + i * dx;
  }

  return xs;
}

Eigen::Vector3d unit_normal(double theta, double phi)
{
  Eigen::Vector3d rhat(std::sin(theta) * std::cos(phi),
                       std::sin(theta) * std::sin(phi), std::cos(theta));

  return rhat;
}

double gaussian(const double t) { return std::exp(-std::pow(t, 2) / 2); }
double skew_gaussian(const double alpha, const double t)
{
  return gaussian(t) * std::erfc(-alpha * t / std::sqrt(2));
}

int grid_sequence(const size_t n)
{
  return (1 - std::pow(-1, n) * (1 + 2 * n)) / 4;
}

std::pair<int, double> split_double(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

double falling_factorial(const double x, const int n)
{
  if(n == 0) return 1;
  double result = x;
  for(int i = 1; i < n; ++i) result *= x - i;

  return result;
}

double spherical_bessel(double arg, int order){
    switch (order){
        case 0:
            return sin( arg ) / arg;
            break;
        case 2:
            return ( 3.0 / pow( arg, 2 ) - 1.0 ) * sin( arg ) / arg - 3.0 * cos( arg ) / pow( arg, 2 );
            break;
    }
}
