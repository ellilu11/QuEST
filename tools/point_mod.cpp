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
#include "../src/quantum_dot.h"
#include "../src/quantum_dot.cpp"

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  std::string taskstr(argv[1]);

  auto dots = make_shared<DotVector>(import_dots("./dots/dots"+taskstr+".cfg"));
  // auto dots = make_shared<DotVector>(import_dots("../build/dots/dots"+taskstr+".cfg"));

  const double hbar = 0.65821193;
  const double omega0 = 4823.67;
  const double dip0 = 0.002536;
  const double c = 194.6704;
  const double mu0 = 2.0133545e-4;
  const double Gamma = pow(dip0,2) * pow(omega0,3) / 
                        ( 3.0 * M_PI * c * hbar / mu0 );
  // 0.120341;

  const double dip_factor = 2.50;
  const double dip = dip0 * dip_factor;
  // const vector<double> dipfactors{ 0.5, 0.7, 1.0, 1.5, 2.0 }; 

  const double T1 = 1 / ( Gamma * pow(dip_factor,2) );
  const double T2 = 2*T1;
 
  const double delta = 0.00 / hbar; // 0.759634; // 1.51927;
  std::cout << "Gamma: " << Gamma << std::endl;
  std::cout << "T1, T2: " << T1 << ", " << T2 << std::endl;
  std::cout << "Sigma omega: " << delta << std::endl;

  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> omega(omega0, delta);

  for(int i = 0; i < dots->size(); ++i){
    (*dots)[i].set_freq( omega(generator) );
    (*dots)[i].set_dipole( Eigen::Vector3d(dip,0,0) );
    (*dots)[i].set_damping( T1, T2 );
  }

  std::ofstream fd("./dots/dots"+taskstr+".cfg");
  // std::ofstream fd("./dots/dots"+taskstr+".cfg");
  fd << std::setprecision(12);

  std::ostream_iterator<QuantumDot> dot_iter(fd);
  std::copy ( dots->begin(), dots->end(), dot_iter );
  fd.close();

  return 0;
}
