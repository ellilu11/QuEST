#include "direct_interaction.h"
#include "../math_utils.h"

constexpr int DIM = 5;

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
      omega(omega), beta(beta), hbar(hbar), rotating(rotating),
      floor_delays(num_srcsrc),
      floor_delays_srcobs(num_srcobs),
      coeffs(boost::extents[num_srcsrc][3]),
      coeffs_srcobs(boost::extents[num_srcobs][3])
{
   build_coeff_table(kernel);
   build_coeff_srcobs_table(kernel);
}

void DirectInteraction::build_coeff_table(
    Propagation::Kernel<cmplx> &kernel)
{
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

    Eigen::Vector3d dip_src = (*dots)[src].dipole();
    Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

    std::vector<Eigen::Matrix3cd> propagator(kernel.coefficients(dr,dip_src));

    for(int i = 0; i < 1; ++i)
        coeffs[pair_idx][i] = dip_obs.dot( propagator[i] * dip_src );
  }
  std::cout << "# Nearfield pairs: " << num_nf << "/" << num_srcsrc << std::endl;

}

void DirectInteraction::build_coeff_srcobs_table(
    Propagation::Kernel<cmplx> &kernel)
{
  // src-obs coefficients
  for(int obs = 0; obs < num_obs; ++obs){
    for(int src = 0; src < num_src; ++src){
      int pair_idx = coord2idxsq( obs, src, num_src );

      Eigen::Vector3d dr(separation((*dots)[src], (*obss)[obs]));
      double dist = dr.norm();

      std::pair<int, double> delay(split_double(dist / (c0 * dt)));
      floor_delays_srcobs[pair_idx] = delay.first;

      Eigen::Vector3d dip_src = (*dots)[src].dipole();

      std::vector<Eigen::Matrix3cd> propagator(kernel.coefficients(dr,dip_src));

      for(int i = 0; i < 3; ++i) 
        coeffs_srcobs[pair_idx][i] = propagator[i] * dip_src;
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

        const int s = //time_idx - floor_delays[pair_idx];
                    std::max(time_idx - floor_delays[pair_idx] ,
                       static_cast<int>(history->array_.index_bases()[1]));

        Eigen::Vector3d dip_src = (*dots)[src].dipole();
        Eigen::Vector3d dip_obs = (*dots)[obs].dipole();

        for (int i = 0; i < 1; ++i){
            results[src] += history->array_[src][s][0][i+1] * coeffs[pair_idx][i];
            results[obs] += history->array_[obs][s][0][i+1] * coeffs[pair_idx][i];
        
        if (time_idx == 2000 && i == 0)
          std::cout << history->array_[src][s][0][i+1] << " "
                    << coeffs[pair_idx][i] << std::endl;

        }
        
      }

  } else {

      // source-observer interactions
      for(int obs = 0; obs < num_obs; ++obs){
        for(int src = 0; src < num_src; ++src){

          int pair_idx = coord2idxsq( obs, src, num_src );

          Eigen::Vector3d dip_src = (*dots)[src].dipole();
          Eigen::Vector3d dip_obs = (*obss)[obs].dipole(); 

          if (obs != src) {
            const int s =
              std::max(time_idx - floor_delays_srcobs[pair_idx],
                       static_cast<int>(history->array_.index_bases()[1]));

            for (int i = 0; i < 2; ++i)
                results[obs] += history->array_[src][s][0][i+1] * dip_obs.dot( coeffs_srcobs[pair_idx][i] );

          } else 
             results[obs] += 
                beta * dip_obs.dot( dip_src / dip_src.squaredNorm() ) * (
			        history->array_[src][time_idx][0][1] * iu * pow(omega,3) +
		            history->array_[src][time_idx][0][2] * 0.0 * pow(omega,2) -
                    history->array_[src][time_idx][0][3] * 0.0 * omega );

        }
      } 

      /*for(int dot = 0; dot < num_src; ++dot){
        Eigen::Vector3d dip_src = (*dots)[dot].dipole();
        Eigen::Vector3d dip_obs = (*obss)[dot].dipole(); 

        results[dot] += beta * dip_obs.dot( dip_src / pow( dip_src.norm(), 2 ) ) * (
			iu * (history->array_[dot][time_idx][0])[RHO_01] -
		       3.0 * (history->array_[dot][time_idx][1])[RHO_01] / omega );

	} */
      
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
