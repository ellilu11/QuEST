#include <complex>

/*int factorial(int n){
    if (n == 1) return 1;
    else return n*factorial(n-1);
}*/


int binomCoeff(int n, int k){
    if (k == 0 || k == n)
        return 1;
    else
        return binomCoeff(n-1, k-1) + binomCoeff(n-1, k);
}

double absCmplx(std::complex<double> z){
    return sqrt(pow(real(z),2)+pow(imag(z),2));

}
