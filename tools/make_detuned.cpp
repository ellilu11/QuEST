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
  auto dots = make_shared<DotVector>(import_dots("../build/dots/dots0.cfg"));

  const double omega0 = 2278.9013;
	const double sigma = 1.0;
  std::cout << "Sigma omega: " << sigma << std::endl;
 
  const unsigned seed =
      std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> omega(omega0, sigma);

  for(int i = 0; i < dots->size(); ++i)
    (*dots)[i].set_freq( omega(generator) );

  std::string taskstr(argv[1]);

  std::ofstream fd("./dots/dots"+taskstr+".cfg");
  // std::ofstream fd("./dots/dots"+taskstr+".cfg");
  fd << std::setprecision(12);

  std::ostream_iterator<QuantumDot> out_iter(fd);
  std::copy ( dots->begin(), dots->end(), out_iter );
  fd.close();

  return 0;
}
