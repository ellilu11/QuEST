#include "direct_interaction.h"
#include "../math_utils.h"

DirectInteraction::DirectInteraction(
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
      num_obs((this->obss)->size()),
      num_interactions( num_src * (num_src - 1) / 2 ),
      num_srcobs( num_src * num_obs ),
      omega(omega), 
      floor_delays(num_interactions),
      coeffs(boost::extents[num_interactions][interp_order + 1]),
      floor_delays_fld(num_srcobs),
      fldcoeffs(boost::extents[num_srcobs][interp_order+1])
{
  build_coeff_table(kernel);
  build_fldcoeff_table(kernel);
}

void DirectInteraction::build_coeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    Eigen::Vector3d dr(separation((*dots)[src], (*dots)[obs]));
    double dist = dr.norm();

    std::pair<int, double> delay(split_double(dist / (c0 * dt)));

    floor_delays[pair_idx] = delay.first;

    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
    Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

    for(int i = 0; i <= interp_order; ++i)
      coeffs[pair_idx][i] = dip_obs.dot( interp_dyads[i] * dip_src );
  }
}

void DirectInteraction::build_fldcoeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int obs = 0; obs < num_obs; ++obs) {
    for(int src = 0; src < num_src; ++src) {
      int pair_idx = coord2idxsq( obs, src, num_src ); 

      Eigen::Vector3d dr(separation((*dots)[src], (*obss)[obs]));
      double dist = dr.norm();
      if ( dist == 0.0 ) continue;

      std::pair<int, double> delay(split_double(dist / (c0 * dt)));

      floor_delays_fld[pair_idx] = delay.first;

      lagrange.evaluate_derivative_table_at_x(delay.second, dt);

      std::vector<Eigen::Matrix3cd> interp_dyads(
          kernel.coefficients(dr, lagrange));

      Eigen::Vector3d dip_src = (*dots)[src].dipole();

      for(int i = 0; i <= interp_order; ++i){
        fldcoeffs[pair_idx][i] = interp_dyads[i] * dip_src;
      }
    }
  } 
}

const InteractionBase::ResultArray &DirectInteraction::evaluate(
    const int time_idx)
{
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho_src, rho_obs;

  results.setZero();
  past_terms_of_convolution.setZero();

  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    for(int i = 1; i <= interp_order; ++i) {
      const int s =
          std::max(time_idx - floor_delays[pair_idx] - i, -history->window);
      const double time = (time_idx - i) * dt;

      rho_obs = (history->get_value(obs, s, 0))[RHO_01];
      rho_src = (history->get_value(src, s, 0))[RHO_01];

      if ( !omega ){
        past_terms_of_convolution[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][i] );
        past_terms_of_convolution[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][i] );
			} else {
        past_terms_of_convolution[src] += rho_obs * coeffs[pair_idx][i] ;
	      past_terms_of_convolution[obs] += rho_src * coeffs[pair_idx][i] ;
	
    	}
    }
   
    const int s = std::max(time_idx - floor_delays[pair_idx], -history->window);
    rho_obs = (history->get_value(obs, s, 0))[RHO_01];
    rho_src = (history->get_value(src, s, 0))[RHO_01];

    if ( !omega ){
      results[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][0] );
      results[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][0] );
    } else {
      results[src] += rho_obs * coeffs[pair_idx][0] ;
	    results[obs] += rho_src * coeffs[pair_idx][0] ;
    } 
  }

  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &DirectInteraction::evaluate_present(
    const int time_idx)
{
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho_src, rho_obs;

  results.setZero();
  
  // source pairwise interactions
  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    const int s = std::max(time_idx - floor_delays[pair_idx], -history->window);
    rho_obs = (history->get_value(obs, s, 0))[RHO_01];
    rho_src = (history->get_value(src, s, 0))[RHO_01];

    if ( !omega ){
      results[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][0] );
      results[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][0] );
    } else {
	    results[src] += rho_obs * coeffs[pair_idx][0] ;
	    results[obs] += rho_src * coeffs[pair_idx][0] ;
			
    }    
  }
 
  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &DirectInteraction::evaluate_field(
    const int time_idx)
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
            std::max(time_idx - floor_delays_fld[pair_idx] - i, -history->window);
        rho_src = (history->get_value(src, s, 0))[RHO_01];

        if ( !omega )
          results[obs] += 2.0 * dip_obs.dot ( 
                            std::real( rho_src * fldcoeffs[pair_idx][i] ) ) ;
        else
          results[obs] += dip_obs.dot ( 
                            rho_src * fldcoeffs[pair_idx][i] ) ;
      }
    }
  }

  return results;
}


int DirectInteraction::coord2idx(int row, int col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

std::pair<int, int> DirectInteraction::idx2coord(const int idx)
{
  const int row = std::floor((std::sqrt(1 + 8 * idx) + 1) / 2.0);
  const int col = idx - row * (row - 1) / 2;

  return std::pair<int, int>(row, col);
}

/*int DirectInteraction::coord2idxsq(int row, int col, int rowlen)
{
  return row*rowlen + col;
}*/
