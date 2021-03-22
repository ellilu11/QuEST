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
  const double omega = 2278.9013;
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 10000.0, T2 = 20000.0;
  const double dip = 5.2917721e-4 * 1.0;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;

  const int num_long = atoi(argv[1]);
  const int num_trans = atoi(argv[2]);
  const int num_dots = num_long * num_trans;

  std::cout << "Num long: " << num_long << " Num trans: " << num_trans << std::endl;
 
  const double num_src = 20;
  const double src_dz = 5.0e-3; // avg distance between srcs
  const double leng = num_src * src_dz;
  const double dz = leng / ( num_long - 1 );
  
  const double r0 = 1.0e-2;
  const double dph = 2*M_PI / num_trans;

  Qdot qd;
  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  const double z0 = -(num_long-1)*dz / 2.0;

  for(int iz = 0; iz < num_long; ++iz) {
    double z = z0 + iz*dz;  

    for(int iph = 0; iph < num_trans; ++iph){
      double phi = iph*dph;
      double x = r0*std::cos(phi);
      double y = r0*std::sin(phi);

      qd = {{x, y, z, omega,
                      T1, T2, dipx, dipy, dipz}};
      dots.push_back(qd);
    }
  }

  min_dist(dots);

  std::string taskstr(argv[3]);

  std::ofstream fd("./dots/obss"+taskstr+".cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
