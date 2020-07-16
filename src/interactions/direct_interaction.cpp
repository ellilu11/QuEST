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
    const double omega,
    const double beta,
    const double hbar,
    const bool rotating)
    : HistoryInteraction(
          std::move(dots), std::move(obss), std::move(history), interp_order, c0, dt),
      num_src((this->dots)->size()),
      num_obs((this->obss) ? (this->obss)->size() : 0),
      num_srcsrc( num_src * (num_src - 1) / 2 ),
      num_srcobs( num_src * num_obs ),
      omega(omega), beta(beta), hbar(hbar),
      floor_delays(num_srcsrc),
      floor_delays_srcobs(num_srcobs),
      coeffs(boost::extents[num_srcsrc+num_src][interp_order + 1]),
      coeffs_srcobs(boost::extents[num_srcobs][4][interp_order+1])
{
   // if (this->obss) std::cout << (this->obss)->size() << std::endl;
   build_coeff_table(kernel);
   build_coeff_srcobs_table(kernel);
}

void DirectInteraction::build_coeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);
  int num_nf = 0;
  double dist0 = 0.0 * c0 / omega;

  // src-src pairwise coefficients
  for(int pair_idx = 0; pair_idx < num_srcsrc; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    Eigen::Vector3d dr(separation((*dots)[src], (*dots)[obs]));
    double dist = dr.norm();
    if (dist < dist0) num_nf++;

    std::pair<int, double> delay(split_double(dist / (c0 * dt)));

    floor_delays[pair_idx] = delay.first;

    std::cout << "Dist: " << dist << std::endl;
    std::cout << "Delay: " << delay.second << std::endl; // gets fractional part!

//    lagrange.evaluate_derivative_table_at_x(0.0, dt);
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
    Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

    for(int i = 0; i <= interp_order; ++i){
       coeffs[pair_idx][i] = dip_obs.dot( interp_dyads[i] * dip_src );
   }
  }
/*    for(int j = 0; j <= 3; j++){
      std::cout << j << " deriv coeffs: " << std::endl;
      for(int i = 0; i <= interp_order; ++i)
        std::cout << lagrange.evaluations[j][i] << std::endl; 
    }
*/
     // src self coefficients
  for(int src = 0; src < num_src; ++src) {
    lagrange.evaluate_derivative_table_at_x(0.0, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(Eigen::Vector3d::Zero(), lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
 
    for(int i = 0; i <= interp_order; ++i){
       coeffs[num_srcsrc+src][i] = dip_src.dot( interp_dyads[i] * dip_src );
       // std::cout << src << " " << i << " " << coeffs[num_srcsrc+src][i] << std::endl;
    }
  }
  std::cout << "# Nearfield pairs: " << num_nf << "/" << num_srcsrc << std::endl;

}

void DirectInteraction::build_coeff_srcobs_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  // src-obs coefficients
  for(int obs = 0; obs < num_obs; ++obs){
    for(int src = 0; src < num_src; ++src){
      int pair_idx = coord2idxsq( obs, src, num_src );

      Eigen::Vector3d dr(separation((*dots)[src], (*obss)[obs]));
      double dist = dr.norm();

      std::pair<int, double> delay(split_double(dist / (c0 * dt)));
      floor_delays_srcobs[pair_idx] = delay.first;

      lagrange.evaluate_derivative_table_at_x(delay.second, dt);

      std::vector<Eigen::Matrix3cd> interp_dyads_deriv0(
          kernel.coeffs_deriv0(lagrange));
      std::vector<Eigen::Matrix3cd> interp_dyads_deriv1(
          kernel.coeffs_deriv1(lagrange));
      std::vector<Eigen::Matrix3cd> interp_dyads_deriv2(
          kernel.coeffs_deriv2(lagrange));
      std::vector<Eigen::Matrix3cd> interp_dyads_deriv3(
          kernel.coeffs_deriv3(lagrange));

      Eigen::Vector3d dip_src = (*dots)[src].dipole();

      for(int i = 0; i <= interp_order; ++i) {
        coeffs_srcobs[pair_idx][0][i] = interp_dyads_deriv0[i] * dip_src;
        coeffs_srcobs[pair_idx][1][i] = interp_dyads_deriv1[i] * dip_src;
        coeffs_srcobs[pair_idx][2][i] = interp_dyads_deriv2[i] * dip_src;
        coeffs_srcobs[pair_idx][3][i] = interp_dyads_deriv3[i] * dip_src;

      }  
    }
  }
/*  for(int i = 0; i <= interp_order; ++i)
    for(int j = 0; j <= 3; ++j<< )
      std::cout << "Deriv order: " << j << " " << coeffs_srcobs[2][j][i] << std::endl;
*/

}

const InteractionBase::ResultArray &DirectInteraction::evaluate(
    const int time_idx)
{
  results.setZero();
  constexpr int RHO_01 = 1;
  const double time = time_idx * dt;

  // source pairwise interactions
    for(int pair_idx = 0; pair_idx < num_srcsrc; ++pair_idx) {
      int src, obs;
      std::tie(src, obs) = idx2coord(pair_idx);

      Eigen::Vector3d dip_src = (*dots)[src].dipole();
      Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

      for(int i = 0; i <= interp_order; ++i) {
        const int s =
            std::max(time_idx - floor_delays[pair_idx] - i,
                     static_cast<int>(history->array_.index_bases()[1]));

        results[src] += (history->array_[obs][s][0])[RHO_01] * coeffs[pair_idx][i];
        results[obs] += (history->array_[src][s][0])[RHO_01] * coeffs[pair_idx][i];

/*        results[src] += std::real( (history->array_[obs][s][0])[RHO_01] ) * coeffs[pair_idx][i] * std::cos( omega*time ) -
                        std::imag( (history->array_[obs][s][0])[RHO_01] ) * coeffs[pair_idx][i] * std::sin( omega*time );

        results[obs] += std::real( (history->array_[src][s][0])[RHO_01] ) * coeffs[pair_idx][i] * std::cos( omega*time ) -
                        std::imag( (history->array_[src][s][0])[RHO_01] ) * coeffs[pair_idx][i] * std::sin( omega*time );
*/      
//        results[src] *= std::exp( -iu*omega*time );
//        results[obs] *= std::exp( -iu*omega*time );


      }
    }

    // source self-interactions
    for(int src = 0; src < num_src; ++src) {
      for(int i = 0; i <= interp_order; ++i) {
        const int s =
            std::max(time_idx - i,
                     static_cast<int>(history->array_.index_bases()[1]));
        results[src] += (history->array_[src][s][0])[RHO_01] * coeffs[num_srcsrc+src][i];

      }
    } 
  return results;
}


const InteractionBase::ResultArray &DirectInteraction::evaluatefld(
    const int time_idx, const int deriv_order)
{
  results.setZero();
  constexpr int RHO_01 = 1;
  const double time = time_idx * dt;

  // source-observer interactions
  for(int obs = 0; obs < num_obs; ++obs){
    for(int src = 0; src < num_src; ++src){

      int pair_idx = coord2idxsq( obs, src, num_src );

      Eigen::Vector3d dip_src = (*dots)[src].dipole();
      Eigen::Vector3d dip_obs = (*obss)[obs].dipole(); 

      for(int i = 0; i <= interp_order; ++i) {
        const int s =
          std::max(time_idx - floor_delays_srcobs[pair_idx] - i,
                    static_cast<int>(history->array_.index_bases()[1]));
        results[obs] += (history->array_[src][s][0])[RHO_01] * dip_obs.dot( coeffs_srcobs[pair_idx][deriv_order][i] );
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

int DirectInteraction::coord2idxsq(int row, int col, int rowlen)
{
  return row*rowlen + col;
}
