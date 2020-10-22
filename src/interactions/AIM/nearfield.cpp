#include "nearfield.h"

AIM::Nearfield::Nearfield(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const double c0,
    const double dt,
    const double h,
    std::shared_ptr<const Grid> grid,
    std::shared_ptr<const Expansions::ExpansionTable> expansion_table,
    Expansions::ExpansionFunction expansion_function,
    Expansions::ExpansionFunction expansion_function_fdtd,
    Normalization::SpatialNorm normalization,
    std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs,
    const double omega)
    : AimBase(dots,
              history,
              interp_order,
              c0,
              dt,
							omega,
              h,
              grid,
              expansion_table,
              expansion_function,
              expansion_function_fdtd,
              normalization),
      interaction_pairs_{std::move(interaction_pairs)},
      shape_({{static_cast<int>(interaction_pairs_->size()),
               grid->max_transit_steps(c0, dt) + 12, 2}}),
      support_(shape_[0], {std::numeric_limits<int>::max(), 0}),
      coefficients_{coefficient_table()}
{
  if(interp_order != 5)
    throw std::runtime_error(
        "Lagrange interpolation of order != 5 is not supported in nearfield "
        "calculations");
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate(const int time_idx)
{
 	constexpr int RHO_01 = 1;
	const double time0 = time_idx * dt;
  cmplx rho0, rho1;

  results.setZero();
	past_terms_of_convolution.setZero(); 

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    for(int t = support_[pair_idx].begin + 1; t < support_[pair_idx].end; ++t) {
      const int s = std::max(
          time_idx - t, -history->window);
			const double time = time0 - t*dt;

      rho0 = (history->get_value(pair.first, s, 0))[RHO_01];
      rho1 = (history->get_value(pair.second, s, 0))[RHO_01];

      if ( !omega_ ){
      	past_terms_of_convolution[pair.first] += 2.0 * std::real( rho1 * coefficients_[pair_idx][t][0] );
      	past_terms_of_convolution[pair.second] += 2.0 * std::real( rho0 * coefficients_[pair_idx][t][1] );
      } else {
        past_terms_of_convolution[pair.first] += 2.0 * std::real( rho1 * coefficients_[pair_idx][t][0] 
                        * std::exp( iu*omega_*time) ) * std::exp( -iu*omega_*time );
                                                                                                      
        past_terms_of_convolution[pair.second] += 2.0 * std::real( rho0 * coefficients_[pair_idx][t][1] 
                        * std::exp( iu*omega_*time) ) * std::exp( -iu*omega_*time );
      }    
    }
   
		const int t = support_[pair_idx].begin;
    const int s = std::max(time_idx - t, -history->window);
    rho1 = (history->get_value(pair.second, s, 0))[RHO_01];
    rho0 = (history->get_value(pair.first, s, 0))[RHO_01];

    if ( !omega_ ){
      results[pair.first] += 2.0 * std::real( rho1 * coefficients_[pair_idx][t][0] );
      results[pair.second] += 2.0 * std::real( rho0 * coefficients_[pair_idx][t][1] );
    } else {
      results[pair.first] += 2.0 * std::real( rho1 * coefficients_[pair_idx][t][0] 
                      * std::exp( iu*omega_*time0) ) * std::exp( -iu*omega_*time0 );
                                                                                                    
      results[pair.second] += 2.0 * std::real( rho0 * coefficients_[pair_idx][t][1] 
                      * std::exp( iu*omega_*time0) ) * std::exp( -iu*omega_*time0 );
    }    
  }
  
  results += past_terms_of_convolution;

  return results;
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate_present_field(
    const int time_idx)
{
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  cmplx rho0, rho1;

  results.setZero();
  
  // source pairwise interactions
  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

		const int t = support_[pair_idx].begin;
    const int s = std::max(time_idx - t, -history->window);
		const double time = time0 - t*dt;

    rho1 = (history->get_value(pair.second, s, 0))[RHO_01];
    rho0 = (history->get_value(pair.first, s, 0))[RHO_01];

    if ( !omega_ ){
      results[pair.first] += 2.0 * std::real( rho1 * coefficients_[pair_idx][t][0] );
      results[pair.second] += 2.0 * std::real( rho0 * coefficients_[pair_idx][t][1] );
    } else {
      results[pair.first] += 2.0 * std::real( rho1 * coefficients_[pair_idx][t][0] 
                      * std::exp( iu*omega_*time) ) * std::exp( -iu*omega_*time );
                                                                                                    
      results[pair.second] += 2.0 * std::real( rho0 * coefficients_[pair_idx][t][1] 
                      * std::exp( iu*omega_*time) ) * std::exp( -iu*omega_*time );
    }    
  }
 
  results += past_terms_of_convolution;

  return results;
}

boost::multi_array<cmplx, 3> AIM::Nearfield::coefficient_table()
{
  boost::multi_array<cmplx, 3> coefficients(shape_);
  std::fill(coefficients.data(),
            coefficients.data() + coefficients.num_elements(), cmplx(0, 0));

  Interpolation::DerivFive lagrange(dt);
  std::cout << "    # Nearfield pairs: " << shape_[0] << std::endl;

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {

    const auto &pair = (*interaction_pairs_)[pair_idx];
    const auto &dot0 = (*dots)[pair.first], &dot1 = (*dots)[pair.second];

    boost::multi_array<Eigen::Vector3cd, 2> stencil_coeffs0(
      boost::extents[lagrange.order()+1][expansion_table->shape()[1]]);
    boost::multi_array<Eigen::Vector3cd, 2> stencil_coeffs1(
      boost::extents[lagrange.order()+1][expansion_table->shape()[1]]);

    for(size_t i = 0; i < expansion_table->shape()[2]; ++i) {
      const auto &e0 = (*expansion_table)[pair.first][0][i];

      for(size_t j = 0; j < expansion_table->shape()[2]; ++j) {
        const auto &e1 = (*expansion_table)[pair.second][0][j];

        if(e0.index == e1.index) continue;

        Eigen::Vector3d dr(grid->spatial_coord_of_box(e1.index) -
                           grid->spatial_coord_of_box(e0.index));

        // == Retardation quantities ==

        const double arg = dr.norm() / (c0 * dt);
        const auto delay = split_double(arg);
        const cmplx norm = normalization(dr);

        support_[pair_idx].begin =
            std::min(support_[pair_idx].begin, delay.first);
        support_[pair_idx].end = std::max(
            support_[pair_idx].end, delay.first + lagrange.order() + 1);

        lagrange.evaluate_derivative_table_at_x(delay.second);

        // == Expansion quantities ==

        const double innerprod =
            dot1.dipole().dot(e1.d0 * e0.d0 * dot0.dipole());

        const std::array<cmplx, 2> dyad{
            norm * std::pow(c0, 2) *
                dot0.dipole().dot(e0.del_sq * e1.d0 * dot1.dipole()),
            norm * std::pow(c0, 2) *
                dot1.dipole().dot(e1.del_sq * e0.d0 * dot0.dipole())};

        for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
          const int convolution_idx = delay.first + poly;
          const cmplx time =
              ( omega_ ?
                (lagrange.evaluations[2][poly] +
                2.0 * iu * omega_ * lagrange.evaluations[1][poly] -
                std::pow(omega_, 2) * lagrange.evaluations[0][poly]) :
                lagrange.evaluations[2][poly] ) *
              innerprod * norm;

            coefficients[pair_idx][convolution_idx][0] +=
							// lagrange.evaluations[0][poly] * innerprod * norm; // Laplace/Helmholtz kernel
              -time + (h_ ? 0.0 : dyad[0] * lagrange.evaluations[0][poly]);

            if(pair.first != pair.second)
              coefficients[pair_idx][convolution_idx][1] +=
								// lagrange.evaluations[0][poly] * innerprod * norm; // Laplace/Helmholtz kernel
                -time + (h_ ? 0.0 : dyad[1] * lagrange.evaluations[0][poly]);
        }

        /*for(int obs = 0; obs < expansion_table->shape()[1]; ++obs){ // obs pt
          if ( !h_ ) break;
          if ( obs == 13 || obs == 14 || obs == 16 || obs == 17 )
            continue;
          if ( obs == 22 || obs == 23 || obs == 25 || obs == 26 )
            continue;

          const auto &e0_obs = (*expansion_table)[pair.first][obs][i];
          const auto &e1_obs = (*expansion_table)[pair.second][obs][j];

          // std::cout << obs << " " << e0_obs.index << std::endl;
          // std::cout << obs << " " <<  << " " << e0_obs.d0 << std::endl;
					
          const Eigen::Vector3cd stencil0 = 
            e0_obs.d0 * e1.d0 * dot1.dipole() * std::pow(c0,2) * norm;
          const Eigen::Vector3cd stencil1 = 
            e1_obs.d0 * e0.d0 * dot0.dipole() * std::pow(c0,2) * norm;

          for(int poly = 0; poly < lagrange.order() + 1; ++poly){
            stencil_coeffs0[poly][obs] += lagrange.evaluations[0][poly] * stencil0;
            if ( pair.first != pair.second)
              stencil_coeffs1[poly][obs] += lagrange.evaluations[0][poly] * stencil1;
          }
        }*/

      } // j
    } // i
    
    /*for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
      if ( !h_ ) break;
        const int convolution_idx = 0 + poly; // assume dist between any two src/obs expansion points is << c0 * dt 
        coefficients[pair_idx][convolution_idx][0] += 
          dot0.dipole().dot( FDTD_Del_Del(stencil_coeffs0[poly]) );
        coefficients[pair_idx][convolution_idx][1] += 
          dot1.dipole().dot( FDTD_Del_Del(stencil_coeffs1[poly]) );
    }*/
    
/*    for(int t = support_[pair_idx].begin; t < support_[pair_idx].end; ++t)
      std::cout << pair_idx << " " << t << " " 
                << coefficients_[pair_idx][t][0] << " " 
                << coefficients_[pair_idx][t][1] << std::endl; 
*/ 
 } // pair_idx 

  return coefficients;
}
