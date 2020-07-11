#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
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

  cout << "Min dist: " << min << " (" << idx1 << "," << idx2 << ")" << endl;

  return min;
}

int main(int argc, char *argv[])
{
  const double c0 = 299.792458;
  const double omega = 2278.9013;
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 5.0, T2 = 10.0;
  const double dipx = 0.0, dipy = 5.2917721e-4, dipz = 0.0;

  const int num_dots = atoi(argv[1]);
  const double ds = 0.20*lambda;
  const double xlen = lambda, ylen = lambda, zlen = lambda;
  const int nx = xlen/ds, ny = ylen/ds, nz = zlen/ds;
  const int dots_per_box = num_dots / ( nx * ny * nz ); 

  const unsigned seed = // 1669889962;
    std::chrono::system_clock::now().time_since_epoch().count();
   // std::cout << seed << std::endl;
  std::default_random_engine generator(seed);
  
  std::vector<Eigen::Vector3d> coords(dots_per_box);
  for (int dot = 0; dot < dots_per_box; ++dot){
    std::uniform_real_distribution<double> xdist(-ds / 2, ds / 2),
        ydist(-ds / 2, ds / 2), zdist(-ds / 2, ds / 2);
    coords[dot] = {xdist(generator), ydist(generator), zdist(generator)};
   // std::cout << coords[dot] << std::endl;
 
  }

  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  for (int k = 0; k < nz; ++k){
    for (int j = 0; j < ny; ++j){
        for (int i = 0; i < nx; ++i){
           for (int dot = 0; dot < dots_per_box; ++dot){
                Qdot qd = {{coords[dot][0] + ds*i, coords[dot][1] + ds*j, coords[dot][2] + ds*k, omega,
                    T1, T2, dipx, dipy, dipz, 0.0, 0.0, 0.0}};
                dots.push_back(qd);
            }
        }
    }
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
