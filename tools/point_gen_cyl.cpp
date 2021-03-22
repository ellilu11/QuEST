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
  const double lambda = 2 * M_PI * c0 / omega;
  const double T1 = 10.0, T2 = 20.0;
  const double dip = 0.002536; // 5.2917721e-4;
  const double dipx = dip, dipy = 0.0, dipz = 0.0;
  const int num_dots = atoi(argv[1]);

  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
    
  // const double ds = 0.0040 * lambda;
  // const double xlen = 0.10 * lambda;
  // const double ylen = xlen;
  const double R = 0.10 * lambda;
  const double zlen = 10 * lambda; // std::max( xlen, 0.015 * num_dots * ds );  
  const double MINDIST = 0.0030;

  // std::cout << "xlen: " << xlen/ds << " ylen: " << ylen/ds << " zlen: " << zlen/ds << std::endl;
  std::cout << "rlen: " << R/lambda << " zlen: " << zlen/lambda << std::endl;
 
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> 
      rsqdist(0, R*R), phdist(0, 2*M_PI), 
      // xdist(-xlen/2, xlen/2), ydist(-ylen/2, ylen/2),
      zdist(-zlen/2, zlen/2);
  // std::uniform_real_distribution<double> u(0, 1);
  // std::uniform_real_distribution<double> v(0, 1);

  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  for(int i = 0; i < num_dots; ++i) {
    double r = sqrt(rsqdist(generator));
    double ph = phdist(generator);

    /*double th = acos( 2*u(generator) - 1 );
    double ph = 2 * M_PI * v(generator);
    double dipx = dip * sin(th) * cos(ph);
    double dipy = dip * sin(th) * sin(ph);
    double dipz = dip * cos(th);*/

    bool lt_minflag = 1;
    Qdot qd;

    while (lt_minflag){
        qd = {{r * cos(ph), r * sin(ph), zdist(generator), omega,
                T1, T2, dipx, dipy, dipz}};
 
        // qd = {{xdist(generator), ydist(generator), zdist(generator), omega,
         //           T1, T2, dipx, dipy, dipz}};
        lt_minflag = lt_min( qd, dots, MINDIST );
    }  
 
    dots.push_back(qd);
  }

  min_dist(dots);

  std::string taskstr(argv[2]);

  std::ofstream fd("./dots/dots"+taskstr+".cfg");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
