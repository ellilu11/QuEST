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
              h,
              grid,
              expansion_table,
              expansion_function,
              expansion_function_fdtd,
              normalization),
      omega_{omega},
      interaction_pairs_{std::move(interaction_pairs)},
      shape_({{static_cast<int>(interaction_pairs_->size()),
               grid->max_transit_steps(c0, dt) + 12, 2}}),
      support_(shape_[0], {std::numeric_limits<int>::max(), 0}),
      floor_delays(shape_[0]),
      coefficients_{coefficient_table()}
{
  if(interp_order != 5)
    throw std::runtime_error(
        "Lagrange interpolation of order != 5 is not supported in nearfield "
        "calculations");
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate(const int time_idx)
{
  results.setZero();
  constexpr int RHO_01 = 1;

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    for(int t = support_[pair_idx].begin; t < support_[pair_idx].end; ++t) {
      const int s = std::max(
          time_idx - t, -history->window);

      results[pair.first] += (history->get_value(pair.second, s, 0))[RHO_01] *
                             coefficients_[pair_idx][t][0];

      results[pair.second] += (history->get_value(pair.first, s, 0))[RHO_01] *
                              coefficients_[pair_idx][t][1];
    }
  }

  return results;
}

boost::multi_array<cmplx, 3> AIM::Nearfield::coefficient_table()
{
  boost::multi_array<cmplx, 3> coefficients(shape_);
  std::fill(coefficients.data(),
            coefficients.data() + coefficients.num_elements(), cmplx(0, 0));

  Interpolation::DerivFive lagrange(dt);
  std::cout << "# Nearfield pairs: " << shape_[0] << std::endl;
  
  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {

    //if (pair_idx % ( shape_[0] / 1000 ) == 0) 
    //    std::cout << pair_idx / ( shape_[0] / 1000 ) << std::endl; 

    const auto &pair = (*interaction_pairs_)[pair_idx];
    const auto &dot0 = (*dots)[pair.first], &dot1 = (*dots)[pair.second];

    for(size_t i = 0; i < expansion_table->shape()[2]; ++i) {
      const auto &e0 = (*expansion_table)[pair.first][0][i];

      for(size_t j = 0; j < expansion_table->shape()[2]; ++j) {
        const auto &e1 = (*expansion_table)[pair.second][0][j];

        if(e0.index == e1.index) continue;

        Eigen::Vector3d dr(grid->spatial_coord_of_box(e1.index) -
                           grid->spatial_coord_of_box(e0.index));

        // == Retardation quantities ==

        const double arg = dr.norm() / (c0 * dt);
        const auto split_arg = split_double(arg);
        floor_delays[pair_idx] = split_arg.first;

        support_[pair_idx].begin =
            std::min(support_[pair_idx].begin, split_arg.first);
        support_[pair_idx].end = std::max(
            support_[pair_idx].end, split_arg.first + lagrange.order() + 1);

        lagrange.evaluate_derivative_table_at_x(split_arg.second);

        // == Expansion quantities ==

        const double innerprod =
            dot1.dipole().dot(e1.d0 * e0.d0 * dot0.dipole());

        const std::array<cmplx, 2> dyad{
            normalization(dr) * std::pow(c0, 2) *
                dot0.dipole().dot(e0.del_sq * e1.d0 * dot1.dipole()),
            normalization(dr) * std::pow(c0, 2) *
                dot1.dipole().dot(e1.del_sq * e0.d0 * dot0.dipole())};

        for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
          const int convolution_idx = split_arg.first + poly;
          const cmplx time =
              ( omega_ ?
                (lagrange.evaluations[2][poly] +
                2.0 * iu * omega_ * lagrange.evaluations[1][poly] -
                std::pow(omega_, 2) * lagrange.evaluations[0][poly]) :
                lagrange.evaluations[2][poly] ) *
              innerprod * normalization(dr);

          if ( !expansion_function_fdtd ) {
            coefficients[pair_idx][convolution_idx][0] +=
                -time + dyad[0] * lagrange.evaluations[0][poly];

            if(pair.first == pair.second) continue;

            coefficients[pair_idx][convolution_idx][1] +=
                -time + dyad[1] * lagrange.evaluations[0][poly];
          }
        }
      } // j
    } // i

    // FDTD del_sq

    // First dot0 (src) -> dot1 (obs)
    for(size_t i = 0; i < expansion_table->shape()[2]; ++i) { // src expansion pt
      if ( !expansion_function_fdtd ) break;
      const auto &e0 = (*expansion_table)[pair.first][0][i];
      std::vector<Eigen::Vector3cd> fld_stencil(27);

      for(int obs = 0; obs < shape_[1]; ++obs){ // obs pt
       fld_stencil[obs] = Eigen::Vector3cd::Zero();

        for(size_t j = 0; j < expansion_table->shape()[2]; ++j) { // expansion pt of obs pt
          const auto &e1 = (*expansion_table)[pair.second][obs][j];

          if(e0.index == e1.index) continue;

          fld_stencil[obs] += e1.d0 * e0.d0 * dot0.dipole();
        }     
      }
      const cmplx dyad1 = 
            dot1.dipole().dot(FDTD_Del_Del(fld_stencil));
      for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
        const int convolution_idx = floor_delays[pair_idx] + poly;
        coefficients[pair_idx][convolution_idx][1] += 
          dyad1 * lagrange.evaluations[0][poly];
      }
    } 

    // Then dot1 (src) -> dot0 (obs)
    for(size_t i = 0; i < expansion_table->shape()[2]; ++i) { // src expansion pt
      if ( !expansion_function_fdtd ) break;
      const auto &e1 = (*expansion_table)[pair.second][0][i];
      std::vector<Eigen::Vector3cd> fld_stencil(27);

      for(int obs = 0; obs < shape_[1]; ++obs){ // obs pt
       fld_stencil[obs] = Eigen::Vector3cd::Zero();

        for(size_t j = 0; j < expansion_table->shape()[2]; ++j) { // expansion pt of obs pt
          const auto &e0 = (*expansion_table)[pair.first][obs][j];

          if(e0.index == e1.index) continue;

          fld_stencil[obs] += e0.d0 * e1.d0 * dot1.dipole();
        }     
      }
      const cmplx dyad0 = 
            dot0.dipole().dot(FDTD_Del_Del(fld_stencil));
      for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
        const int convolution_idx = floor_delays[pair_idx] + poly;
        coefficients[pair_idx][convolution_idx][0] += 
          dyad0 * lagrange.evaluations[0][poly];
      }
    } 

   } // pair_idx 

  return coefficients;
}
