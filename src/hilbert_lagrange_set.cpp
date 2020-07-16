#include <iostream>

Interpolation::HilbertLagrangeSet::HilbertLagrangeSet(const int order)
    : evaluations(boost::extents[Interpolation::NUM_DERIVATIVES][order + 1]),
      order_(order)
{
}

/*Interpolation::HilbertLagrangeSet::HilbertLagrangeSet(const double time,
                                                      const int order,
                                                      const double dt,
                                                      const double omega,
                                                      const double phi)
    : HilbertLagrangeSet(order)
{
  assert(x >= 0);  // Don't extrapolate!
}*/

void Interpolation::HilbertLagrangeSet::evaluate_table_at_x(
    const double x, const double dt, const double t, const double omega, const double phi /* = 1 */)
{
  switch(order_){
    case 1 : {
      evaluations[0][1] = -( std::cos( omega*(t+phi) ) - std::cos( omega*(t-dt+phi) ) + dt*omega*std::sin( omega*(t-t0+phi) ) );
      evaluations[0][0] = std::cos( omega*(t+phi) ) - std::cos( omega*(t-dt+phi) ) + dt*omega*std::sin( omega*(t+phi) );
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

