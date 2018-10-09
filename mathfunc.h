#include <cmath>
#include <complex>

#ifndef FMM2D_MATHFUNC_H
#define FMM2D_MATHFUNC_H

const std::complex<double> I(0.0,1.0);

// int factorial(int n);

int binomCoeff(int n, int k);

double absCmplx(std::complex<double> z);

#endif
