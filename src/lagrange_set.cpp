#include "lagrange_set.h"
#include <iostream>

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const int order)
    : evaluations(boost::extents[Interpolation::UniformLagrangeSet::NUM_DERIVATIVES][order + 1]),
      order_(order)
{
}

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const double x,
                                                      const int order,
                                                      const double dt /* = 1 */)
    : UniformLagrangeSet(order)
{
  assert(x > 0);  // Don't extrapolate!
  evaluate_derivative_table_at_x(x, dt);
}

void Interpolation::UniformLagrangeSet::evaluate_derivative_table_at_x(
    const double x, const double dt /* = 1 */)
{
  double y = 1.0 - x;
  double y2 = pow(y,2), y3 = pow(y,3), y4 = pow(y,4);

  switch(order_){
    case 4 : {

      evaluations[0][4] = 1.0/24.0 * (-2.0*y - y2 + 2.0*y3 + y4);
      evaluations[0][3] = -1.0/6.0 * (-3.0*y - y2 + 3.0*y3 + y4);
      evaluations[0][2] = 1.0/4.0 * (-6.0*y + y2 + 4.0*y3 + y4);
      evaluations[0][1] = -1.0/6.0 * (-6.0 - 5.0*y + 5.0*y2 + 5.0*y3 + y4);
      evaluations[0][0] = 1.0/24.0 * y * (6.0 + 11.0*y + 6.0*y2 + y3);

      evaluations[1][4] = 1.0/24.0 * (-2.0 - 2.0*y + 6.0*y2 + 4.0*y3);
      evaluations[1][3] = -1.0/6.0 * (-3.0 - 2.0*y + 9.0*y2 + 4.0*y3);
      evaluations[1][2] = 1.0/4.0 * (-6.0 + 2.0*y + 12.0*y2 + 4.0*y3);
      evaluations[1][1] = -1.0/6.0 * (-5.0 + 10.0*y + 15.0*y2 + 4.0*y3);
      evaluations[1][0] = 1.0/24.0 * (6.0 + 22.0*y + 18.0*y2 + 4.0*y3);
      
      evaluations[2][4] = 1.0/24.0 * (-2.0 + 12.0*y + 12.0*y2);
      evaluations[2][3] = -1.0/6.0 * (-2.0 + 18.0*y + 12.0*y2);
      evaluations[2][2] = 1.0/4.0 * (2.0 + 24.0*y + 12.0*y2);
      evaluations[2][1] = -1.0/6.0 * (10.0 + 30.0*y + 12.0*y2);
      evaluations[2][0] = 1.0/24.0 * (22.0 + 36.0*y + 12.0*y2);

      evaluations[3][4] = 1.0/24.0 * (12.0 + 24.0*y);
      evaluations[3][3] = -1.0/6.0 * (18.0 + 24.0*y);
      evaluations[3][2] = 1.0/4.0 * (24.0 + 24.0*y);
      evaluations[3][1] = -1.0/6.0 * (30.0 + 24.0*y);
      evaluations[3][0] = 1.0/24.0 * (36.0 + 24.0*y);

    }
  }
/*    for(int basis_id = 0; basis_id <= order_; ++basis_id) {
    double d0_product = 1, d1_product = 1, d2_product = 1, d3_product = 1;
    double d1_sum = 0, d2_sum = 0, d3_sum = 0;
    for(int m = 0; m <= order_; ++m) {
      if(m == basis_id) continue;

      double numer = (x - m);

      d0_product *= numer / (basis_id - m)  

      if(numer != 0) {
        d1_sum -= 1 / numer;  // Note the minus sign!
        d2_sum -= std::pow(numer, -2);
        d3_sum -= std::pow(numer, -3);
      }

      d0_product *= (x - m) / (basis_id - m);

      for(int l=0; l <=order_; ++l){
        if(l == m) continue;
        d1_product *= (x - l) / (basis_id - l);
      }
      d1_sum -= d1_product;


    }

//    std::cout << x << std::endl;

    evaluations[0][basis_id] = d0_product;

    evaluations[1][basis_id] = (evaluations[0][basis_id] * d1_sum);
    evaluations[2][basis_id] = evaluations[0][basis_id] *
        (d2_sum + std::pow(d1_sum, 2));
    evaluations[3][basis_id] = evaluations[0][basis_id] *
       (2.0 * d3_sum + 3.0 * d1_sum * d2_sum + std::pow(d1_sum, 3));

    evaluations[1][basis_id] = d1_product;

 }
*/

  for(int i = 0; i <= order_; ++i) {
    evaluations[1][i] *= std::pow(dt, -1);
    evaluations[2][i] *= std::pow(dt, -2);
    evaluations[3][i] *= std::pow(dt, -3);
    
/*    if ( x == 0 ) {
      for(int j = 0; j <= 3; ++j)
        std::cout << evaluations[j][i] << " ";
      std::cout << std::endl;
    }*/
  }

}

Interpolation::HilbertLagrangeSet::HilbertLagrangeSet(const int order)
    : evaluations(boost::extents[Interpolation::HilbertLagrangeSet::NUM_DERIVATIVES][order + 1]),
      order_(order)
{
}

void Interpolation::HilbertLagrangeSet::evaluate_table_at_x(
    const double t, const double t0, const double dt, const double omega)
{
  switch(order_){
    case 1 : {
      evaluations[0][1] = -( std::cos( omega*(t+t0) ) - std::cos( omega*(t-dt+t0) ) + dt*omega*std::sin( omega*(t-dt+t0) ) );
      evaluations[0][0] = std::cos( omega*(t+t0) ) - std::cos( omega*(t-dt+t0) ) + dt*omega*std::sin( omega*(t+t0) );
    }
  }

  for(int i = 0; i <= order_; ++i) {
    evaluations[0][i] *= std::pow(dt, -1) * std::pow(omega, -2);
    // evaluations[1][i] *= std::pow(dt, -1) * std::pow(omega, -1);
    
/*    if ( x == 0 ) {
      for(int j = 0; j <= 3; ++j)
        std::cout << evaluations[j][i] << " ";
      std::cout << std::endl;
    }*/
  }

}

Interpolation::DerivFive::DerivFive(const double dt)
    : evaluations(boost::extents[Interpolation::DerivFive::NUM_DERIVATIVES][5 + 1]),
      dt_{dt}
{
}

void Interpolation::DerivFive::evaluate_derivative_table_at_x(const double x)
{
  for(int basis_id = 0; basis_id <= 5; ++basis_id) {
    evaluations[0][basis_id] = d0(basis_id, x);
    evaluations[1][basis_id] = d1(basis_id, x) / std::pow(dt_, 1);
    evaluations[2][basis_id] = d2(basis_id, x) / std::pow(dt_, 2);
  }
}

// clang-format off
double Interpolation::DerivFive::d0(const int i, const double x) const
{
  switch(i) {
    case 0 : return 1+x*(-137.0/60+x*(15.0/8+x*(-17.0/24+x*(1.0/8-x/120))));
    case 1 : return x*(5+x*(-77.0/12+x*(71.0/24+x*(-7.0/12+x/24))));
    case 2 : return x*(-5+x*(107.0/12+x*(-59.0/12+x*(13.0/12-x/12))));
    case 3 : return x*(10.0/3+x*(-13.0/2+x*(49.0/12+x*(-1+x/12))));
    case 4 : return x*(-5.0/4+x*(61.0/24+x*(-41.0/24+x*(11.0/24-x/24))));
    case 5 : return x*(1.0/5+x*(-5.0/12+x*(7.0/24+x*(-1.0/12+x/120))));
    default: return 0;
  }
}

double Interpolation::DerivFive::d1(const int i, const double x) const
{
  switch(i) {
    // Again, note the negative sign
    case 0 : return 137.0/60+x*(-15.0/4+x*(17.0/8+x*(-1.0/2+x/24)));
    case 1 : return -5+x*(77.0/6+x*(-71.0/8+x*(7.0/3-(5*x)/24)));
    case 2 : return 5+x*(-107.0/6+x*(59.0/4+x*(-13.0/3+(5*x)/12)));
    case 3 : return -10.0/3+x*(13+x*(-49.0/4+x*(4-(5*x)/12)));
    case 4 : return 5.0/4+x*(-61.0/12+x*(41.0/8+x*(-11.0/6+(5*x)/24)));
    case 5 : return -1.0/5+x*(5.0/6+x*(-7.0/8+x*(1.0/3-x/24)));
    default : return 0;
  }
}

double Interpolation::DerivFive::d2(const int i, const double x) const
{
  switch(i) {
    case 0 : return 15.0/4+x*(-17.0/4+x*(3.0/2-x/6));
    case 1 : return -77.0/6+x*(71.0/4+x*(-7+(5*x)/6));
    case 2 : return 107.0/6+x*(-59.0/2+x*(13-(5*x)/3));
    case 3 : return -13+x*(49.0/2+x*(-12+(5*x)/3));
    case 4 : return 61.0/12+x*(-41.0/4+x*(11.0/2-(5*x)/6));
    case 5 : return -5.0/6+x*(7.0/4+x*(-1+x/6));
    default: return 0;
  }
}
// clang-format on
