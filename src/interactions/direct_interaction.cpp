#include "direct_interaction.h"
#include "../math_utils.h"

DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const DotVector> obss,
    std::shared_ptr<const Integrator::History<matrix_cmplx>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt,
    const double omega)
    : HistoryInteraction(
          std::move(dots), std::move(obss), std::move(history), interp_order, c0, dt),
      // num_src((this->dots)->size()),
      // num_obs((this->obss) ? (this->obss)->size() : 0),
      // num_interactions( num_src * (num_src - 1) / 2 ),
      omega(omega), 
      floor_delays(num_interactions),
      coeffs(boost::extents[num_interactions][interp_order + 1]),
      floor_delays_srcobs(num_srcobs),
      fldcoeffs(boost::extents[num_srcobs][interp_order+1]),
      cross_coeffs(boost::extents[num_srcobs][interp_order+1])
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

    std::pair<int, double> delay(split_double(0.0));
    // std::pair<int, double> delay(split_double(dist / (c0 * dt)));

    floor_delays[pair_idx] = delay.first;

    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();

    for(int i = 0; i <= interp_order; ++i)
      coeffs[pair_idx][i] = interp_dyads[i] * dip_src;
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

      Eigen::Vector3d rhat = dr / dist;
 
      std::pair<int, double> delay(split_double(dist / (c0 * dt)));

      floor_delays_srcobs[pair_idx] = delay.first;

      lagrange.evaluate_derivative_table_at_x(delay.second, dt);

      std::vector<Eigen::Matrix3cd> interp_dyads(
          kernel.coefficients(dr, lagrange));

      Eigen::Vector3d dip_src = (*dots)[src].dipole();

      for(int i = 0; i <= interp_order; ++i){
        fldcoeffs[pair_idx][i] = interp_dyads[i] * dip_src;
        cross_coeffs[pair_idx][i] = 
            (rhat.cross(fldcoeffs[pair_idx][i])).conjugate(); // not sure why need conjugate
      
      }
    }
  } 
}

const InteractionBase::ResultArray &DirectInteraction::evaluate(
    const int time_idx)
{
  const double time0 = time_idx * dt;
  cmplx rho_src, rho_obs;
  matrix_cmplx rho_matrix;

  results.setZero();
  past_terms_of_convolution.setZero();

  for (int obs = 0; obs < num_src; ++obs){
    for (int src = 0; src < num_src; ++src){
      if ( src !=obs ){
 
        int pair = coord2idxrect(src, obs, num_src);
        int pairint = coord2idx(src, obs);

        for(int i = 1; i <= interp_order; ++i) {
          const int s =
              std::max(time_idx - floor_delays[pairint] - i, -history->window);
          const double time = (time_idx - i) * dt;

          // work out trace 
          rho_matrix = traceOverDims( history->get_value(s, 0), num_src, src );
          rho_src = rho_matrix(0,1);

          if ( !omega )
            past_terms_of_convolution[pair] = 2.0 * (rho_src * coeffs[pairint][i]).real() ;
          else
            past_terms_of_convolution[pair] = rho_src * coeffs[pairint][i] ;
        }
       
        const int s = std::max(time_idx - floor_delays[pairint], -history->window);
        rho_matrix = traceOverDims( history->get_value(s, 0), num_src, src );
        rho_src = rho_matrix(0,1);

        if ( !omega )
          results[pair] = 2.0 * ( rho_src * coeffs[pairint][0] ).real();
        else 
          results[pair] = rho_src * coeffs[pairint][0] ;
      }
    }
  }

  results += past_terms_of_convolution;
  return results;
}

const InteractionBase::ResultArray &DirectInteraction::evaluate_present(
    const int time_idx)
{
  const double time0 = time_idx * dt;
  cmplx rho_src, rho_obs;
  matrix_cmplx rho_matrix;

  results.setZero();
  
  // source pairwise interactions
  for (int obs = 0; obs < num_src; ++obs){
    for (int src = 0; src < num_src; ++src){
      if ( src !=obs ){
 
        int pair = coord2idxrect(src, obs, num_src);
        int pairint = coord2idx(src, obs);

        const int s = std::max(time_idx - floor_delays[pairint], -history->window);
        // rho_src = (history->get_value(src, s, 0))[RHO_01];
        rho_matrix = traceOverDims( history->get_value(s, 0), num_src, src );
        rho_src = rho_matrix(0,1);

        if ( !omega )
          results[pair] = 2.0 * ( rho_src * coeffs[pairint][0] ).real();
        else 
          results[pair] = rho_src * coeffs[pairint][0] ;
      }
    }
  }
     
  results += past_terms_of_convolution;
  return results;
}

const InteractionBase::ResultArray &DirectInteraction::evaluate_field(
    const int time_idx, const bool flag)
{
  /*constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho_src, rho_obs;

  results.setZero();
  past_terms_of_convolution.setZero();

  for(int obs = 0; obs < num_obs; ++obs) {
      
    Eigen::Vector3d dip_obs = (*obss)[obs].dipole();
    Eigen::Vector3d pos_obs = (*obss)[obs].position();
    Eigen::Vector3d rhat = pos_obs / pos_obs.norm();    

    for(int src = 0; src < num_src; ++src) {
      int pair_idx = coord2idxsq( obs, src, num_src ); 
      Eigen::Vector3d dr(separation((*dots)[src], (*obss)[obs]));

      Eigen::Vector3d pos_src = (*dots)[src].position();

      for(int i = 0; i <= interp_order; ++i) {
        const int s =
            // std::max(time_idx - i, -history->window);
            std::max(time_idx - floor_delays_srcobs[pair_idx] - i, -history->window);
        rho_src = (history->get_value(src, s, 0))[RHO_01];

        if ( flag ) {
          if ( !omega )
            results[obs] += 2.0 * std::real( dip_obs.dot ( 
                              rho_src * cross_coeffs[pair_idx][i] ) ) ;
          else
            results[obs] += dip_obs.dot ( 
                              rho_src * cross_coeffs[pair_idx][i] ); // * std::exp( iu*omega*rhat.dot(pos_src)/c0 );
        } else {

          if ( !omega )
            results[obs] += 2.0 * std::real( dip_obs.dot ( 
                              rho_src * fldcoeffs[pair_idx][i] ) ) ;
          else
            results[obs] += dip_obs.dot ( 
                              rho_src * fldcoeffs[pair_idx][i] ); // * std::exp( iu*omega*rhat.dot(pos_src)/c0 );

          }
      }
    }
  }
*/
  return results;
}


