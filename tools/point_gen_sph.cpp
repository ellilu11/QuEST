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

double spherical_volume_element(
  const double r0, const double r1, 
  const double theta0, const double theta1, const double phi0, const double phi1)
{
  return ( pow(r1,3) - pow(r0,3) ) / 3.0 * spherical_area_element( 1.0, theta0, theta1, phi0, phi1 );
}

int main(int argc, char *argv[])
{
  const double c0 = 299.792458;
  const double omega = 2278.9013;
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 10000.0, T2 = 20000.0;
  const double dip = 5.2917721e-4 * 1.0;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;

  const int nr = atoi(argv[1]);
  const int nth = atoi(argv[2]);
  const int nph = atoi(argv[3]);
 
  const double dph = 2.0*M_PI / nph;
  const double dth = M_PI / nth;

  const int num_dots = nr*(nph*(nth-1)+2);

  std::cout << "Nr: " << nr << " Nth: " << nth << " Nph: " << nph << " Num dots: " << num_dots << std::endl;

  const double num_src = 20;
  const double src_dz = 5.0e-3; // avg distance between srcs
  const double R = ( num_src + 2 ) * src_dz / 2.0;
  const double dr = R / nr;

  std::cout << "Radius: " << R << std::endl;

  Qdot qd;
  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  double total_volume = 0;

  for(int ir = 0; ir < nr; ++ir){
    double r = dr * ( ir + 0.5 );
 
    double volume_element =
      spherical_volume_element(r-dr/2.0, r+dr/2.0, 0, dth/2.0, 0, 2*M_PI );
    total_volume += 2.0*volume_element;   

    dots.push_back({{0, 0, r, volume_element, 0, 0, 0, 0, 0}});
    for(int iph = 0; iph < nph; ++iph){
      for(int ith = 1; ith < nth; ++ith){
        double phi = 2.0 * M_PI * iph / nph;
        double theta = M_PI * ith / nth;

        double x = r * sin( theta ) * cos ( phi );
        double y = r * sin( theta ) * sin ( phi );
        double z = r * cos( theta );

        volume_element = 
          spherical_volume_element(
            r-dr/2.0, r+dr/2.0, theta-dth/2.0, theta+dth/2.0, phi-dph/2.0, phi+dph/2.0 );
        total_volume += volume_element;   

        qd = {{x, y, z, volume_element, 0, 0, 0, 0, 0}};
        dots.push_back(qd);
      }
    }
    volume_element =
      spherical_volume_element(r-dr/2.0, r+dr/2.0, 0, dth/2.0, 0, 2*M_PI );
    dots.push_back({{0, 0, -r, volume_element, 0, 0, 0, 0, 0}});
  }

  std::cout << "Sum of volume elements: " << total_volume << std::endl;
  std::cout << "Actual volume: " << 4.0/3.0*M_PI*pow(R,3) << std::endl;

  min_dist(dots);

  std::ofstream fd("./dots/obss_sph.cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();


/*  std::ofstream fd2("./dots/dotsobs.cfg");
  fd2 << std::setprecision(12);

  for(const auto &d : obss) {
    fd2 << d << endl;
  }

  fd2.close();
*/
  return 0;
}
