#include "direct_interaction.h"
#include "../math_utils.h"

DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const DotVector> obss,
    std::shared_ptr<const Integrator::History<Eigen::Vector4cd>> history,
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
      coeffs_srcobs(boost::extents[num_srcobs][interp_order+1])
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
  double dist0 = 0 * c0 / omega;

  // src-src pairwise coefficients
  for(int pair_idx = 0; pair_idx < num_srcsrc; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    Eigen::Vector3d dr(separation((*dots)[src], (*dots)[obs]));
    double dist = dr.norm();
    if (dist < dist0) num_nf++;

    std::pair<int, double> delay(split_double(dist / (c0 * dt)));

    floor_delays[pair_idx] = delay.first;

//    std::cout << "Delay: " << delay.second << std::endl;

//    lagrange.evaluate_derivative_table_at_x(delay.second, dt);
    lagrange.evaluate_derivative_table_at_x(0.0, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
    Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

    for(int i = 0; i <= interp_order; ++i){
       coeffs[pair_idx][i] = dip_obs.dot( interp_dyads[i] * dip_src );
    }
  }
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

      std::vector<Eigen::Matrix3cd> interp_dyads(
          kernel.coefficients(dr, lagrange));

      Eigen::Vector3d dip_src = (*dots)[src].dipole();

      for(int i = 0; i <= interp_order; ++i) 
        coeffs_srcobs[pair_idx][i] = interp_dyads[i] * dip_src;
    }
  }
}

const InteractionBase::ResultArray &DirectInteraction::evaluate(
    const int time_idx)
{
  results.setZero();
  constexpr int RHO_01 = 1;
  const double time = time_idx * dt;

  if (!obss){
    // source pairwise interactions
      for(int pair_idx = 0; pair_idx < num_srcsrc; ++pair_idx) {
        int src, obs;
        std::tie(src, obs) = idx2coord(pair_idx);

        Eigen::Vector3d dip_src = (*dots)[src].dipole();
        Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

        for(int i = 0; i <= 0; ++i) {
          const int s =
              std::max(time_idx - floor_delays[pair_idx] - i,
                       static_cast<int>(history->array_.index_bases()[1]));

          results[src] += (history->array_[obs][s][0])[RHO_01] * coeffs[pair_idx][i];
          results[obs] += (history->array_[src][s][0])[RHO_01] * coeffs[pair_idx][i];
         
          if (time_idx == 2000 && i == 0) 
            std::cout << (history->array_[obs][s][0])[RHO_01] << " "
                      << coeffs[pair_idx][i] << std::endl;


        }
      }

      // source self-interactions
      /*for(int src = 0; src < num_src; ++src) {
        for(int i = 0; i <= interp_order; ++i) {
          const int s =
              std::max(time_idx - i,
                       static_cast<int>(history->array_.index_bases()[1]));
          results[src] += (history->array_[src][s][0])[RHO_01] * coeffs[num_srcsrc+src][i];

        }
      } */ 
  } else {

      // source-observer interactions
      for(int obs = 0; obs < num_obs; ++obs){
        for(int src = 0; src < num_src; ++src){

          int pair_idx = coord2idxsq( obs, src, num_src );

          Eigen::Vector3d dip_src = (*dots)[src].dipole();
          Eigen::Vector3d dip_obs = (*obss)[obs].dipole(); 

          for(int i = 0; i <= interp_order; ++i) {
            if (obs != src) {
               const int s =
                  std::max(time_idx - floor_delays_srcobs[pair_idx] - i,
                           static_cast<int>(history->array_.index_bases()[1]));
                results[obs] += (history->array_[src][s][0])[RHO_01] * dip_obs.dot( coeffs_srcobs[pair_idx][i] );
            } else {
               const int s =
                  std::max(time_idx - i,
                           static_cast<int>(history->array_.index_bases()[1]));
                results[obs] += (history->array_[src][s][0])[RHO_01] * dip_obs.dot( coeffs_srcobs[pair_idx][i] );
            }
          }
        }
      } 
  }

  return results;
}


const InteractionBase::ResultArray &DirectInteraction::evaluatefld(
    const int time_idx)
{

  // constexpr int RHO_01 = 1;
  
  // std::vector<Eigen::Vector3cd> flds_vec(num_src+num_obs); 
  flds.setZero();
/*
  // source pairwise interactions
  for(int pair_idx = 0; pair_idx < num_srcsrc; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    for(int i = 0; i <= interp_order; ++i) {
      const int s =
          std::max(time_idx - floor_delays[pair_idx] - i,
                   static_cast<int>(history->array_.index_bases()[1]));

      Eigen::Vector3d dip_src = (*dots)[src].dipole();
      Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

      flds_vec[src] += hbar * (history->array_[obs][s][0])[RHO_01] * 
          ( fldcoeffs_obs[pair_idx][i] + iu * rrcoeffs[pair_idx] * dip_obs / pow( dip_obs.norm(), 2 ) );

      flds_vec[obs] += hbar * (history->array_[src][s][0])[RHO_01] * 
          ( fldcoeffs_src[pair_idx][i] + iu * rrcoeffs[pair_idx] * dip_src / pow( dip_src.norm(), 2 ) );

    }
  }

  // source self-interactions
  for(int dot = 0; dot < num_src; ++dot){
    Eigen::Vector3d dip_dot = (*dots)[dot].dipole();
      flds_vec[dot] += hbar * (history->array_[dot][time_idx][0])[RHO_01] *
        iu * beta * dip_dot / pow( dip_dot.norm(), 2 );
  }
  
  // source-observer interactions
  for(int obs = 0; obs < num_obs; ++obs){
    for(int src = 0; src < num_src; ++src){

      int pair_idx = coord2idxsq( obs, src, num_src );

      for(int i = 0; i <= interp_order; ++i) {
        const int s =
          std::max(time_idx - floor_delays_srcobs[pair_idx] - i,
                   static_cast<int>(history->array_.index_bases()[1]));

      Eigen::Vector3d dip_src = (*dots)[src].dipole();
 
      flds_vec[num_src+obs] += hbar * (history->array_[src][s][0])[RHO_01] * 
          ( fldcoeffs_srcobs[pair_idx][i] + iu * rrcoeffs_srcobs[pair_idx] * dip_src / pow( dip_src.norm(), 2 ) );
      }
    }
  } 

  // get (complex) magnitude of fields
  for(int dot = 0; dot < (num_src + num_obs); ++dot)
      flds[dot] = flds_vec[dot].norm();
      // flds[dot] = (*dots)[dot].dipole().dot(flds_vec[dot]) / hbar;*/

  return flds;
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
