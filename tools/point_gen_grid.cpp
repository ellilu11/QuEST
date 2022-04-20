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
  const double c0 = 299.792458;
  const double n = 1.54;
  const double c = c0 / n;
  const double omega = 4823.67;
  const double lambda = 2 * M_PI * c / omega;
  const double T1 = 0.0, T2 = 0.0;
  const double dip = 0.002536; // 5.2917721e-4 * 1.0;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;
  const int num_dots = atoi(argv[1]);
  const int ns = pow(num_dots, 1.0/2.0); // num dots per dimension
  const int nssq = ns*ns;

  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
  
  const double L = 10 * lambda;  
  const double dist0 = L / ns; // 2.0 * lambda / ns; // 2.0e-1; // mean distance between dots;
	const double r0 = 0.0*dist0; 

  std::cout << "dist0: " << dist0 << " r0: " << r0 << std::endl;
 
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> u(0, 1);
  std::uniform_real_distribution<double> v(0, 1);

  Qdot qd;
  std::vector<Qdot> dots;
  dots.reserve(num_dots);

	double x0 = -(ns-1)*dist0 / 2.0;
  double y0 = x0;
  double z0 = 0;

  for(int i = 0; i < num_dots; ++i) {
    int idz = i;
    int iz = idz / nssq;
    idz -= iz*nssq;
    int iy = idz / ns;
    int ix = idz % ns;

    double theta = std::acos(2.0*u(generator)-1.0);
    double phi = 2.0*M_PI*v(generator);

    double dx0 = r0*std::sin(theta)*std::cos(phi);
    double dy0 = r0*std::sin(theta)*std::sin(phi);
    double dz0 = r0*std::cos(theta);

  	qd = {{x0+ix*dist0+dx0, y0+iy*dist0+dy0, 0, omega,
                    T1, T2, dipx, dipy, dipz}};
    dots.push_back(qd);
  }

  min_dist(dots);

  std::string taskstr(argv[2]);

  std::ofstream fd("./dots/obss2D.cfg");
  // std::ofstream fd("./dots/dots"+taskstr+".cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
