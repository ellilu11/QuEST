#include "direct_interaction.h"
#include "../math_utils.h"

DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt,
    const double omega,
    const bool rotating)
    : HistoryInteraction(
          std::move(dots), std::move(history), interp_order, c0, dt),
      num_src((this->dots)->size()),
      num_interactions( num_src * (num_src - 1) / 2 ),
      omega(omega), 
      rotating(rotating),
      floor_delays(num_interactions),
      coeffs(boost::extents[num_interactions+num_src][interp_order + 1])
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

//    lagrange.evaluate_derivative_table_at_x(0.0, dt);
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
    Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

    // std::cout << "Pairwise coeffs: " << std::endl;
    for(int i = 0; i <= interp_order; ++i){
       coeffs[pair_idx][i] = dip_obs.dot( interp_dyads[i] * dip_src );
 //      std::cout << coeffs[pair_idx][i] << std::endl;
    }
  }
// src self coefficients
  for(int src = 0; src < num_src; ++src) {
    lagrange.evaluate_derivative_table_at_x(0.0, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads_self(
        kernel.coefficients(Eigen::Vector3d::Zero(), lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
 
    // std::cout << "Self coeffs: " << std::endl;
    for(int i = 0; i <= interp_order; ++i){
       coeffs[num_interactions+src][i] = dip_src.dot( interp_dyads_self[i] * dip_src );
       // std::cout << coeffs[num_interactions+src][i] << std::endl;
    }

  }
}

const InteractionBase::ResultArray &DirectInteraction::evaluate(
    const int time_idx)
{
  results.setZero();
  constexpr int RHO_01 = 1;
  const double time0 = time_idx * dt;
  
  if ( time_idx <= 1 )

  // source pairwise interactions
    for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
      int src, obs;
      std::tie(src, obs) = idx2coord(pair_idx);

      Eigen::Vector3d dip_src = (*dots)[src].dipole();
      Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

      for(int i = 0; i <= interp_order; ++i) {
        const int s =
            std::max(time_idx - floor_delays[pair_idx] - i, -history->window);
        const double time = time0 - i*dt;

        if ( !rotating ){
          results[src] += 2.0 * std::real( (history->get_value(obs, s, 0))[RHO_01] * coeffs[pair_idx][i] );
          results[obs] += 2.0 * std::real( (history->get_value(src, s, 0))[RHO_01] * coeffs[pair_idx][i] );
        } else {
          results[src] += 2.0 * std::real( (history->get_value(obs, s, 0))[RHO_01] * coeffs[pair_idx][i] 
                          * std::exp( iu*omega*time) ) * std::exp( -iu*omega*time );
                                                                                                        
          results[obs] += 2.0 * std::real( (history->get_value(src, s, 0))[RHO_01] * coeffs[pair_idx][i] 
                          * std::exp( iu*omega*time) ) * std::exp( -iu*omega*time );
          
        }    
        if (time_idx <= 1)
          std::cout << i << " " << (history->get_value(0, s, 0))[RHO_01]  << " " << (history->get_value(1, s, 0))[RHO_01] << std::endl; 
 
      }
    }

    // source self-interactions
    for(int src = 0; src < num_src; ++src) {
      for(int i = 0; i <= interp_order; ++i) {
        const int s =
            std::max(time_idx - i, -history->window);
        if ( !rotating ){
          results[src] += 2.0 * std::real( (history->get_value(src, s, 0))[RHO_01] * coeffs[num_interactions+src][i] );
        } else {
          results[src] += 2.0 * std::real( (history->get_value(src, s, 0))[RHO_01] * coeffs[num_interactions+src][i] 
                          * std::exp( iu*omega*time0) ) * std::exp( -iu*omega*time0 );
        }

      }
    } 

       // std::cout << results[0] << " " << results[1] << std::endl;

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
