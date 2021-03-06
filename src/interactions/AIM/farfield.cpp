#include "farfield.h"

AIM::Farfield::Farfield(
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
      table_dimensions_{grid->circulant_shape(c0, dt, interp_order)},
      propagation_table_{make_propagation_table()},
      source_table_{spacetime::make_vector3d<cmplx>(table_dimensions_)},
      obs_table_{spacetime::make_vector3d<cmplx>(table_dimensions_)},
      spatial_vector_transforms_{spatial_fft_plans()}
{
/*  std::cout << table_dimensions_[1] << " " 
            << table_dimensions_[2] << " "
            << table_dimensions_[3] << std::endl;
*/
  auto clear = [](auto &table) {
    std::fill(table.data(), table.data() + table.num_elements(), cmplx(0, 0));
  };

  clear(source_table_);
  clear(obs_table_);
}

void AIM::Farfield::fill_source_table(const int step)
{
  using namespace Expansions::enums;

  const int wrapped_step = step % table_dimensions_[0];
  const auto p = &source_table_[wrapped_step][0][0][0][0];

  std::fill(p, p + 3 * 8 * grid->size(), cmplx(0, 0));

  for(auto dot_idx = 0u; dot_idx < expansion_table->shape()[0]; ++dot_idx) {
    for(auto expansion_idx = 0u; expansion_idx < expansion_table->shape()[2];
        ++expansion_idx) {
      const Expansions::Expansion &e =
          (*expansion_table)[dot_idx][0][expansion_idx];
      Eigen::Vector3i coord = grid->idx_to_coord(e.index);

      Eigen::Map<Eigen::Vector3cd> grid_field(
          &source_table_[wrapped_step][coord(0)][coord(1)][coord(2)][0]);

      // This is the seam between what's stored in the History (density matrix
      // elements) and the electromagnetic source quantities. Ideally the AIM
      // code should not have knowledge of this to better encapsulate
      // "propagation," but this is good enough for now.
      Eigen::Vector3cd source_field = e.d0 * (*dots)[dot_idx].dipole() *
                                      (history->get_value(dot_idx, step, 0))[RHO_01];
      grid_field += source_field;

      /*if ( step == 50 )
        std::cout << dot_idx << " " << expansion_idx << " " <<
          coord.transpose() << " " << grid_field.transpose() << std::endl;*/
    }
  }
  //if ( step == 50 )
  //  std::cout << std::endl;
}

void AIM::Farfield::propagate(const int step)
{
  const auto wrapped_step = step % table_dimensions_[0];
  const auto nb = 8 * grid->size(); // = table_dimensions_[1] * table_dimensions_[2] * table_dimensions_[3]
  const std::array<int, 5> front = {{wrapped_step, 0, 0, 0, 0}};

  const auto s_ptr = &source_table_(front);
  fftw_execute_dft(spatial_vector_transforms_.forward,
                   reinterpret_cast<fftw_complex *>(s_ptr),
                   reinterpret_cast<fftw_complex *>(s_ptr));

  Eigen::Map<Eigen::Array3Xcd> observers(&obs_table_(front), 3, nb);
  observers = 0;

  for(int i = 0; i < table_dimensions_[0]; ++i) {
    // If (step - i) runs "off the end", just propagate src[0][...]
    auto wrap = std::max(step - i, 0) % table_dimensions_[0];

    Eigen::Map<Eigen::ArrayXcd> prop(&propagation_table_[i][0][0][0], nb);
    Eigen::Map<Eigen::Array3Xcd> src(&source_table_[wrap][0][0][0][0], 3, nb);

    // Use broadcasting to do the x, y, and z component propagation
    observers += src.rowwise() * prop.transpose();
  }

  const auto o_ptr = &obs_table_(front);
  fftw_execute_dft(spatial_vector_transforms_.backward,
                   reinterpret_cast<fftw_complex *>(o_ptr),
                   reinterpret_cast<fftw_complex *>(o_ptr));
}

void AIM::Farfield::fill_results_table(const int step)
{
  results = 0;
	const double time = step*dt;	

  for(auto dot_idx = 0u; dot_idx < expansion_table->shape()[0]; ++dot_idx) {
 
    // calculate non-FDTD (e.g. time derivative) field first
    Eigen::Vector3cd field = Eigen::Vector3cd::Zero();
    for(auto expansion_idx = 0u; expansion_idx < expansion_table->shape()[2];
        ++expansion_idx) {
      const Expansions::Expansion &e =
          (*expansion_table)[dot_idx][0][expansion_idx];
      Eigen::Vector3i coord = grid->idx_to_coord(e.index);
      field += 
        expansion_function(obs_table_, {{step, coord(0), coord(1), coord(2)}}, e);
      // Don't use a _wrapped_ step here; the expansion_function needs knowledge
      // of where it's being called in the complete timeline to accommodate
      // boundary conditions
      /*if ( step == 50 ) {
        int wrapped_step = step % table_dimensions_[0];
        Eigen::Map<Eigen::Vector3cd> obs_vec(&obs_table_[wrapped_step][coord(0)][coord(1)][coord(2)][0], 1);
      
        std::cout << dot_idx << " " << expansion_idx << " " <<
            coord.transpose() << " " << obs_vec.transpose() << std::endl;
      }*/
    }
    
    // then calculate FDTD (e.g. spatial derivative) field 
    //std::vector<Eigen::Vector3cd> fld_stencil(27); 
/*    boost::multi_array<Eigen::Vector3cd, 1> fld_stencil(boost::extents[27]);

    for(auto obs_idx = 0u; obs_idx < 27; ++obs_idx) { 
      if ( h_ == 0 ) break;
      if ( obs_idx == 13 || obs_idx == 14 || obs_idx == 16 || obs_idx == 17 )
        continue;
      if ( obs_idx == 22 || obs_idx == 23 || obs_idx == 25 || obs_idx == 26 )
        continue;

      fld_stencil[obs_idx] = Eigen::Vector3cd::Zero();
      for(auto expansion_idx = 0u; expansion_idx < expansion_table->shape()[2];
          ++expansion_idx) {
        const Expansions::Expansion &e =
            (*expansion_table)[dot_idx][obs_idx][expansion_idx];
        Eigen::Vector3i coord = grid->idx_to_coord(e.index);
        fld_stencil[obs_idx] += 
          expansion_function_fdtd(
            obs_table_, {{step, coord(0), coord(1), coord(2)}}, e);
      }

     // if ( step == 50 )
     //   std::cout << dot_idx << " " << obs_idx << " " << (fld_stencil[obs_idx]).transpose() << std::endl;
 
   }

    // use fields at stencil points to calculate FDTD
    Eigen::Vector3cd deldel_field = h_ ? FDTD_Del_Del( fld_stencil ) : Eigen::Vector3cd::Zero(); 
*/ 

    // finally sum fields and calculate Rabi freq
    // results(dot_idx) += 2.0 * std::real( (field+deldel_field).dot((*dots)[dot_idx].dipole()) );
		if ( omega_ )
			results(dot_idx) += field.dot((*dots)[dot_idx].dipole());
		else
			results(dot_idx) += 2.0 * std::real( field.dot((*dots)[dot_idx].dipole()) );
			//results(dot_idx) += 2.0 * std::real( field.dot((*dots)[dot_idx].dipole()) *
			//														std::exp( iu*omega_*time) ) * std::exp( -iu*omega_*time );		
  }
}

spacetime::vector<cmplx> AIM::Farfield::make_propagation_table() const
{

  spacetime::vector<cmplx> g_mat(table_dimensions_);

  const int num_gridpts =
      table_dimensions_[1] * table_dimensions_[2] * table_dimensions_[3];
  TransformPair circulant_plan = {
      fftw_plan_many_dft(3, &table_dimensions_[1], table_dimensions_[0],
                         reinterpret_cast<fftw_complex *>(g_mat.data()),
                         nullptr, 1, num_gridpts,
                         reinterpret_cast<fftw_complex *>(g_mat.data()),
                         nullptr, 1, num_gridpts, FFTW_FORWARD, FFTW_MEASURE),
      nullptr};

  fill_gmatrix_table(g_mat);

  // Transform the circulant vectors into their equivalently-diagonal
  // representation. Buckle up.
  
  fftw_execute(circulant_plan.forward);

  // This accounts for FFTW's *un*normalized transform -- it takes the least
  // amount of computational effort to put all of the normalizations here.

  Eigen::Map<Eigen::ArrayXcd> gs(g_mat.data(), g_mat.num_elements());
  gs /= num_gridpts;

  return g_mat;
}

void AIM::Farfield::fill_gmatrix_table(
    spacetime::vector<cmplx> &gmatrix_table) const
{  // Build the circulant vectors that define the G "matrices." Since the G
  // matrices are Toeplitz (and symmetric), they're uniquely determined by
  // their first row. The first row gets computed here then mirrored to make a
  // list of every circulant (and thus FFT-able) vector. This function needs to
  // accept a non-const reference to a spacetime::vector (instead of just
  // returning such an array) to play nice with FFTW and its workspaces.

  std::fill(gmatrix_table.data(),
            gmatrix_table.data() + gmatrix_table.num_elements(),
            cmplx(0.0, 0.0));

  Interpolation::UniformLagrangeSet interp(interp_order);

  for(int x = 0; x < grid->shape()[0]; ++x) {
    for(int y = 0; y < grid->shape()[1]; ++y) {
      for(int z = 0; z < grid->shape()[2]; ++z) {
        const size_t box_idx = grid->coord_to_idx({x, y, z});

        if(box_idx == 0) continue;

        const Eigen::Vector3d dr =
            grid->spatial_coord_of_box(box_idx) - grid->spatial_coord_of_box(0);

        const double arg = dr.norm() / (c0 * dt);
        const auto split_arg = split_double(arg);

        interp.evaluate_derivative_table_at_x(split_arg.second, dt);

        for(int p = 0; p < interp.order() + 1; ++p) {
          gmatrix_table[split_arg.first + p][x][y][z] =
              interp.evaluations[0][p] * normalization(dr);
        }
      }
    }
  }

  spacetime::fill_circulant_mirror(gmatrix_table);
}

TransformPair AIM::Farfield::spatial_fft_plans()
{
  // Set up FFTW plans to transform projected source distributions. Due to the
  // requirements of the circulant extension, these plans perform transforms of
  // length 2 n_{x,y,z} to accommodate the requisite zero padding. While they're
  // constructed to work on the head of `source_table` (that is, what would be
  // the I_0 source), the advanced FFTW interface allows them to stride forward
  // to equivalently transform the source currents at every timestep.

  constexpr int num_transforms = 3; 
  constexpr int transform_rank = 3; // # of dimensions to transform
  constexpr int dist_between_elements = 3;
  constexpr int dist_between_transforms = 1;

  auto make_plan = [&](const int sign) {
    return fftw_plan_many_dft(
        transform_rank, &table_dimensions_[1], num_transforms,
        reinterpret_cast<fftw_complex *>(source_table_.data()), nullptr,
        dist_between_elements, dist_between_transforms,
        reinterpret_cast<fftw_complex *>(source_table_.data()), nullptr,
        dist_between_elements, dist_between_transforms, sign, FFTW_MEASURE);
  };

  return {make_plan(FFTW_FORWARD), make_plan(FFTW_BACKWARD)};
}
