#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

using std::cout;
using std::endl;

typedef std::array<double, 12> Qdot;

std::ostream &operator<<(std::ostream &out, const Qdot &qd)
{
  std::copy(qd.begin(), qd.end(), std::ostream_iterator<double>(out, " "));
  return out;
}

double dist(const Qdot &qd1, const Qdot &qd2)
{
  Eigen::Vector3d dr(qd2[0] - qd1[0], qd2[1] - qd1[1], qd2[2] - qd1[2]);
  return dr.norm();
}

double min_dist(const std::vector<Qdot> &dots)
{
  int idx1 = -1, idx2 = -1;
  double min = 1e100;
  for(size_t i = 0; i < dots.size() - 1; ++i) {
    for(size_t j = i + 1; j < dots.size(); ++j) {
      double dij = dist(dots[i], dots[j]);
      if(dij < min) {
        min = dij;
        idx1 = i;
        idx2 = j;
      }
    }
  }

  // cout << "Min dist: " << min << " (" << idx1 << "," << idx2 << ")" << endl;

  return min;
}

int main(int argc, char *argv[])
{
  const double c0 = 299.792458;
  const double omega = 2278.9013;
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 10000.0, T2 = 20000.0;
  const double dipx = 5.2917721e-4, dipy = 0.0, dipz = 0.0;

  const int num_dots = atoi(argv[1]);
 
  /*const double ds_base = 0.200 * lambda;
  const double ds_n = 10.0;
  const double ds = ds_base * pow( 2, atof(argv[2]) / ds_n );
 */
  // const int expansion = atoi(argv[1]);
  const double mult1 = 1.0, mult2 = 2.0;
  const double xlen = mult1*lambda, ylen = mult1*lambda, zlen = mult1*lambda;
  const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  const double x1 = mult2*lambda, y1 = mult2*lambda, z1 = mult2*lambda;

  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  // src dots
  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> xdist(-xlen / 2 + x0, xlen / 2 + x0),
      ydist(-ylen / 2 + y0, ylen / 2 + y0), zdist(-zlen / 2 + z0, zlen / 2 + z0);
  for(int i = 0; i < num_dots/2; ++i) {
    Qdot qd = {{xdist(generator), ydist(generator), zdist(generator), omega,
                T1, T2, dipx, dipy, dipz, 0.0, 0.0, 0.0}};
     dots.push_back(qd);
  }

  // trg dots
  const unsigned seed1 = seed + 1;
  std::default_random_engine generator1(seed1);
  std::uniform_real_distribution<double> xdist1(-xlen / 2 + x1, xlen / 2 + x1),
      ydist1(-ylen / 2 + y1, ylen / 2 + y1), zdist1(-zlen / 2 + z1, zlen / 2 + z1);
  for(int i = 0; i < num_dots/2; ++i) {
   Qdot qd1 = {{xdist1(generator1), ydist1(generator1), zdist1(generator1), omega,
                T1, T2, dipx, dipy, dipz, 0.0, 0.0, 0.0}};
    dots.push_back(qd1);
  } 

  /*std::sort(dots.begin(), dots.end(), [](const Qdot &a, const Qdot &b) {
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] <
           b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
  });*/

  min_dist(dots);

  std::ofstream fd("dots.cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
