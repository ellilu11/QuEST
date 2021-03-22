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

double spherical_area_element(
  const double r0, const double theta0, const double theta1, const double phi0, const double phi1)
{
  return -pow(r0, 2) * ( phi1 - phi0 ) * ( std::cos(theta1) - std::cos(theta0) );
}

int main(int argc, char *argv[])
{
  const double c0 = 299.792458;
  const double omega = 2278.9013;
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 10000.0, T2 = 20000.0;
  const double dip = 5.2917721e-4 * 1.0;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;

  const int nth = atoi(argv[1]);
  const int nph = atoi(argv[2]);
  const int num_dots = nph*(nth-1)+2;

  const double dph = 2.0*M_PI / nph;
  const double dth = M_PI / nth;
 
  std::cout << "Nth: " << nth << " Nph: " << nph << " Num dots: " << num_dots << std::endl;

  const double num_src = 5;
  const double src_dz = 5.0e-3; // avg distance between srcs
  const double r0 = 0.1; // ( num_src + 2 ) * src_dz / 2.0;

  std::cout << "Radius: " << r0 << std::endl;

  Qdot qd;
  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  double area_element = spherical_area_element(r0, 0, dth/2.0, 0, 2*M_PI );
  double total_area = area_element;
  dots.push_back({{r0, 0, 0, area_element, 0, 0, 0, 0, 0}});
  
  for(int iph = 0; iph < nph; ++iph){
    for(int ith = 1; ith < nth; ++ith){
      double phi = iph * dph;
      double theta = ith * dth;

      // cyclic permutation of (x,y,z)! (assuming dipole is along x)
      double x = r0 * cos( theta );
      double y = r0 * sin( theta ) * cos ( phi );
      double z = r0 * sin( theta ) * sin ( phi );
      
      area_element =
        spherical_area_element(r0, theta-dth/2.0, theta+dth/2.0, phi-dph/2.0, phi+dph/2.0 );
        qd = {{x, y, z, area_element, 0, 0, 0, 0, 0}};
      dots.push_back(qd);
 
      total_area += area_element;

    /*if ( step == 0 ) 
      std::cout << obs << " " << theta-dth/2.0 << " " << theta+dth/2.0 << " " 
                            << phi-dph/2.0 << " " << phi+dph/2.0 << std::endl;
 */
    }
  }
  area_element = spherical_area_element(r0, 0, dth/2.0, 0, 2*M_PI );
  dots.push_back({{-r0, 0, 0, area_element, 0, 0, 0, 0, 0}});
  total_area += area_element;

  std::cout << "Sum of area elements: " << total_area << std::endl;
  std::cout << "Actual surface area: " << 4.0*M_PI*r0*r0 << std::endl;

  min_dist(dots);

  std::ofstream fd("./dots/obss_sphsurf.cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
