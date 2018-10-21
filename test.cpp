#include <iostream>
#include <cmath>
#include <complex>

double funct(double _a, double _b) {
  return std::log(_a)/log(_b);
}

int main()
{
    double a = 20;
    double b = 10;
    double c = funct(a, b);
    std::cout << "log(a)/log(b) = " << c << '\n';
}
