#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

using std::cout;
using std::endl;

typedef std::array<double, 9> Qdot;

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

bool lt_min(const Qdot &qd, const std::vector<Qdot> &dots, const double mindist){
    int num_dots = dots.size();
    bool flag = 0;

    for (int i = 0; i < num_dots; ++i){
        if ( dist(qd, dots[i]) < mindist ){
            flag = 1;
            break;
        }
    }

    return flag;
}

int main(int argc, char *argv[])
{
  const double c0 = 299.792458 / 1.54;
  const double omega = 4823.67;
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 8.31, T2 = 2*T1;
  const double dip = 0.002536; // 5.2917721e-4 * 1.0;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;
  const int nx = atoi(argv[1]);
  const int ny = atoi(argv[2]);
  const int num_dots = nx * ny;

  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
  
  const double dist0 = 0.2 * lambda; // 2.0 * lambda / ns; // 2.0e-1; // mean distance between dots;

  std::cout << "dist0: " << dist0 << std::endl;
 
  Qdot qd;
  std::vector<Qdot> dots;
  dots.reserve(num_dots);

	double x0 = -(nx-1)*dist0 / 2.0;
  double y0 = -(ny-1)*dist0 / 2.0;
  double z0 = 0;

  for(int i = 0; i < num_dots; ++i) {
    int idz = i;
    int iz = idz / ( nx*ny );
    idz -= iz*nx*ny;
    int iy = idz % ny;
    int ix = idz / ny;

    /*double theta = std::acos(2.0*u(generator)-1.0);
    double phi = 2.0*M_PI*v(generator);

    double dx0 = r0*std::sin(theta)*std::cos(phi);
    double dy0 = r0*std::sin(theta)*std::sin(phi);
    double dz0 = r0*std::cos(theta);
*/
  	qd = {{x0+ix*dist0, y0+iy*dist0, 0, omega,
                    T1, T2, dipx, dipy, dipz}};
    dots.push_back(qd);
  }

  min_dist(dots);

  std::string taskstr(argv[3]);

  // std::ofstream fd("./dots/obss.cfg");
  std::ofstream fd("./dots/dots"+taskstr+".cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
