#include "self_interaction.h"
#include "../math_utils.h"

SelfInteraction::SelfInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const DotVector> obss,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt,
    const double omega)
    : HistoryInteraction(
          std::move(dots), std::move(obss), std::move(history), interp_order, c0, dt),
      num_src((this->dots)->size()),
      num_obs((this->obss) ? (this->obss)->size() : 0),
      num_srcobs( num_src * num_obs ),
      omega(omega), 
      coeffs(boost::extents[num_src][interp_order + 1]),
      fldcoeffs(boost::extents[num_srcobs][interp_order + 1])
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
 
    for(int i = 0; i <= interp_order; ++i)
       coeffs[src][i] = dip_src.dot( interp_dyads_self[i] * dip_src );

  }
}

void SelfInteraction::build_fldcoeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int obs = 0; obs < num_obs; ++obs) {
    for(int src = 0; src < num_src; ++src) {
      int pair_idx = coord2idxsq( obs, src, num_src );

      Eigen::Vector3d dr(separation((*dots)[src], (*obss)[obs]));
      double dist = dr.norm();
      if ( dist == 0.0 ) {
        lagrange.evaluate_derivative_table_at_x(0.0, dt);

        std::vector<Eigen::Matrix3cd> interp_dyads_self(
            kernel.coefficients(Eigen::Vector3d::Zero(), lagrange));

        Eigen::Vector3d dip_src = (*dots)[src].dipole();

        for(int i = 0; i <= interp_order; ++i){
          fldcoeffs[pair_idx][i] = interp_dyads_self[i] * dip_src;
          std::cout << src << " " << obs << " " << fldcoeffs[pair_idx][i] << std::endl;
        }
      }
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
      else 
        past_terms_of_convolution[src] += rho_src * coeffs[src][i];
    }

    const int s = std::max(time_idx, -history->window);
    rho_src = (history->get_value(src, s, 0))[RHO_01];

    if ( !omega )
      results[src] += 2.0 * std::real( rho_src * coeffs[src][0] );
    else 
      results[src] += rho_src * coeffs[src][0];
  } 
  
  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &SelfInteraction::evaluate_present(
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
    else
      results[src] += rho_src * coeffs[src][0];
  } 

  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &SelfInteraction::evaluate_field(
    const int time_idx, const bool flag)
{
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho_src, rho_obs;

  results.setZero();
  past_terms_of_convolution.setZero();

  for(int obs = 0; obs < num_obs; ++obs) {
      
    Eigen::Vector3d dip_obs = (*obss)[obs].dipole();

    for(int src = 0; src < num_src; ++src) {
      int pair_idx = coord2idxsq( obs, src, num_src );
  
      for(int i = 0; i <= interp_order; ++i) {
        const int s =
            std::max(time_idx - i, -history->window);
        rho_src = (history->get_value(src, s, 0))[RHO_01];

        if ( !omega )
          results[obs] += 2.0 * std::real( dip_obs.dot ( 
                            rho_src * fldcoeffs[pair_idx][i] ) ) ;
        else
          results[obs] += dip_obs.dot ( 
                            rho_src * fldcoeffs[pair_idx][i] ) ;
      }
    }
  }

  return results;
}


