#include "direct_interaction.h"
#include "../math_utils.h"

DirectInteraction::DirectInteraction(
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
      num_interactions( num_src * (num_src - 1) / 2 ),
      omega(omega), 
      floor_delays(num_interactions),
      coeffs(boost::extents[num_interactions][interp_order + 1])
{
   build_coeff_table(kernel);
}

void DirectInteraction::build_coeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  // src-src pairwise coefficients
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

    // std::cout << "Pairwise coeffs: " << std::endl;
    for(int i = 0; i <= interp_order; ++i){
      coeffs[pair_idx][i] = dip_obs.dot( interp_dyads[i] * dip_src );
      // std::cout << i << " " << coeffs[pair_idx][i] << std::endl;
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

  // source pairwise interactions
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
				const auto phi = std::exp( iu*omega*time );
        past_terms_of_convolution[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][i] 
                        * phi ) * std::conj( phi );
                                                                                                      
        past_terms_of_convolution[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][i] 
                        * phi ) * std::conj( phi );
    	}  
    }
   
    const int s = std::max(time_idx - floor_delays[pair_idx], -history->window);
    rho_obs = (history->get_value(obs, s, 0))[RHO_01];
    rho_src = (history->get_value(src, s, 0))[RHO_01];

     if ( !omega ){
      results[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][0] );
      results[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][0] );
    } else {
  	  const auto phi0 = std::exp( iu*omega*time0 );
     	results[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][0] 
                      * phi0 ) * std::conj( phi0 );
                                                                                                    
      results[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][0] 
                      * phi0 ) * std::conj( phi0 );
      // if ( time_idx == 10 ) std::cout << phi0 << " " << std::conj( phi0 ) << std::endl;
    } 
  }
  
	if (time_idx == 50) std::cout << std::exp( iu*omega*time0 ) << std::endl;

  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &DirectInteraction::evaluate_present_field(
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
	  	const auto phi0 = std::exp( iu*omega*time0 );
      results[src] += 2.0 * std::real( rho_obs * coeffs[pair_idx][0] 
                      * phi0 ) * std::conj( phi0 );
                                                                                                    
      results[obs] += 2.0 * std::real( rho_src * coeffs[pair_idx][0] 
                      * phi0 ) * std::conj( phi0 );
    }    
  }
 
  results += past_terms_of_convolution;

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

int DirectInteraction::coord2idxsq(int row, int col, int rowlen)
{
  return row*rowlen + col;
}
