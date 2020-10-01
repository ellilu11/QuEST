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


double gaussian(const double t) { return std::exp(-std::pow(t,2) / 2); }
double skew_gaussian(const double alpha, const double t)
{
  return gaussian(t) * std::erfc(-alpha * t / std::sqrt(2));
}

//double cmplx_norm(std::vector<cmplx> vec)
//{

//}

/*
double gaussian(const double t, const double a, const double b, const double c, const int deriv) { 
  
  double func = std::exp( -( a*std::pow(t, 2) + b*t + c ) );

  switch ( deriv ){
    case 0 :
      return func; 
    case 1 :
      return ( -2.0*a*t-b ) * func;
    case 2 :
      return ( pow(-2.0*a*t-b,2) - 2*a ) * func;
    case 3 :
      return ( pow(-2.0*a*t-b,3) - 6.0*a*(-2.0*a*t-b) ) * func;
  }

}

double skew_gaussian(const double alpha, const double t)
{
  return gaussian(t,1,0,0,0) * std::erfc(-alpha * t / std::sqrt(2));
}
*/
int grid_sequence(const size_t n)
{
  return (1 - std::pow(-1, n) * (1 + 2 * n)) / 4;
}

Eigen::Vector3i idx_to_coord(const size_t idx, int dim)
{
  int dimsq = pow(dim,2);
  const int x = idx / dimsq;
  idx -= x * dimsq;
  const int y = idx / dim;
  const int z = idx % dim;

  Eigen::Vector3i coord(x,y,z);

  return coord;
}

Eigen::Vector3i idx_to_delta(const size_t idx, int dim)
{
  Eigen::Vector3i coord = idx_to_coord(idx, dim);
  Eigen::Vector3i delta(grid_sequence(coord[0]), 
                        grid_sequence(coord[1]),
                        grid_sequence(coord[2]))

  return delta;
}

/*size_t inv_grid_sequence(const int delta)
{
  return delta > 0 ? 2*delta - 1 : -2*delta;
}
size_t coord_to_idx(const Eigen::Vector3i &coord, int dim)
{
  return pow(dim,2)*coord(0) + dim*coord(1) + coord(2);
}

size_t delta_to_idx(const Eigen::Vector3i &delta, int dim)
{
  Eigen::Vector3i n(inv_grid_sequence(delta[0]),
                    inv_grid_sequence(delta[1]),
                    inv_grid_sequence(delta[2]);

  return coord_to_idx( n, dim );
}*/

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
