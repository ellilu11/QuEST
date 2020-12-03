#include "self_interaction.h"
#include "../math_utils.h"

SelfInteraction::SelfInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt,
    const double omega)
    : HistoryInteraction(
          std::move(dots), std::move(history), interp_order, c0, dt),
      num_src((this->dots)->size()),
      omega(omega), 
      coeffs(boost::extents[num_src][interp_order + 1])
{
   build_coeff_table(kernel);
}

void SelfInteraction::build_coeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int src = 0; src < num_src; ++src) {
    lagrange.evaluate_derivative_table_at_x(0.0, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads_self(
        kernel.coefficients(Eigen::Vector3d::Zero(), lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
 
    for(int i = 0; i <= interp_order; ++i){
       coeffs[src][i] = dip_src.dot( interp_dyads_self[i] * dip_src );
       // std::cout << coeffs[src][i] << std::endl;
    }

  }
}

const InteractionBase::ResultArray &SelfInteraction::evaluate(
    const int time_idx)
{
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho_src;

  results.setZero();
  past_terms_of_convolution.setZero();

  for(int src = 0; src < num_src; ++src) {
    for(int i = 1; i <= interp_order; ++i) {
      const int s =
          std::max(time_idx - i, -history->window);
      rho_src = (history->get_value(src, s, 0))[RHO_01];

      if ( !omega )
        past_terms_of_convolution[src] += 2.0 * std::real( rho_src * coeffs[src][i] );
      else { 
        past_terms_of_convolution[src] += rho_src * coeffs[src][i];
 
				/*const auto phi0 = std::exp( iu*omega*time0 );
        past_terms_of_convolution[src] += 2.0 * std::real( rho_src * coeffs[src][i] 
                        * phi0 ) * std::conj( phi0 );
        */
      }
    }

    const int s = std::max(time_idx, -history->window);
    rho_src = (history->get_value(src, s, 0))[RHO_01];

    if ( !omega )
      results[src] += 2.0 * std::real( rho_src * coeffs[src][0] );
    else { 
	    results[src] += rho_src * coeffs[src][0];

 	 	  /*const auto phi0 = std::exp( iu*omega*time0 );
      results[src] += 2.0 * std::real( rho_src * coeffs[src][0] 
                      * phi0 ) * std::conj( phi0 );
      */
		}
  } 
  
  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &SelfInteraction::evaluate_present_field(
    const int time_idx)
{
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho_src;

  results.setZero();
  
  // source self-interactions
  for(int src = 0; src < num_src; ++src) {
  
    const int s = std::max(time_idx, -history->window);
    rho_src = (history->get_value(src, s, 0))[RHO_01];

    if ( !omega )
      results[src] += 2.0 * std::real( rho_src * coeffs[src][0] );
    else {
	    results[src] += rho_src * coeffs[src][0];

  		/*const auto phi0 = std::exp( iu*omega*time0 );
      results[src] += 2.0 * std::real( rho_src * coeffs[src][0] 
                      * phi0 ) * std::conj( phi0 );
      */
		}
  } 

  results += past_terms_of_convolution;

  return results;
}

