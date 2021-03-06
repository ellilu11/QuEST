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
  const double omega = 4823.67; // 2278.9013;
  const double lambda = 0.26; // 2 * M_PI * c0 / omega;
  const double T1 = 20.0, T2 = 40.0;
  const double dip = 0.002536; // 5.2917721e-4;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;
  const int num_dots = atoi(argv[1]);

  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
    
  const double ds = 1.0 * lambda;// 0.10 * lambda;
  const double xlen = ds; // ds;
  const double ylen = xlen;
  const double zlen = xlen; // std::max( xlen, 0.015 * num_dots * ds );  
  const double MINDIST = 0.0030;

  std::cout << "xlen: " << xlen/lambda << " ylen: " << ylen/lambda << " zlen: " << zlen/lambda << std::endl;
  // std::cout << "rlen: " << R/ds << " zlen: " << zlen/ds << std::endl;
 
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> 
      xdist(-xlen/2, xlen/2), ydist(-ylen/2, ylen/2),
      zdist(-zlen/2, zlen/2);

  Qdot qd;
  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  for(int i = 0; i < num_dots; ++i) {

    bool lt_minflag = 1;
    while (lt_minflag){
 
        qd = {{xdist(generator), ydist(generator), zdist(generator), omega,
                    T1, T2, dipx, dipy, dipz}};
        lt_minflag = lt_min( qd, dots, MINDIST );
    }  
 
    dots.push_back(qd);
  }

  // additionally, generate evenly spaced observer "dots" on surface of outer sphere
  /*const int nph = atoi(argv[2]);
  const int nth = atoi(argv[3]);
  for(int iph = 0; iph < nph; ++iph){
    for(int ith = 1; ith < nth; ++ith){
        double phi = 2.0 * M_PI * iph / nph;
        double theta = M_PI * ith / nth;
        double x = R * sin( theta ) * cos ( phi );
        double y = R * sin( theta ) * sin ( phi );
        double z = R * cos( theta );
        qd = {{x, y, z, 0, 0, 0, 0, 0, 0}};
        obss.push_back(qd);
    }
  }
  // generate a dot at each pole of the sphere
  obss.push_back({{0, 0, R, 0, 0, 0, 0, 0, 0}});
  obss.push_back({{0, 0, -R, 0, 0, 0, 0, 0, 0}});
*/
  /*std::sort(dots.begin(), dots.end(), [](const Qdot &a, const Qdot &b) {
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] <
           b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
  });*/

  min_dist(dots);

  std::string taskstr(argv[2]);

  std::ofstream fd("./dots/dots"+taskstr+".cfg");
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
