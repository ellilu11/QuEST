#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(const std::shared_ptr<const DotVector> dots,
                                   const std::shared_ptr<const DotVector> obss,
                                   const std::shared_ptr<const Pulse> pulse,
                                   const int interp_order,
                                   const double c0,
                                   const double dt,
                                   const double hbar,
                                   const bool rotating)
    : HistoryInteraction(
        std::move(dots), std::move(obss), std::move(history), interp_order, c0, dt), 
        pulse(std::move(pulse)), hbar(hbar), rotating(rotating),
        num_src((this->dots)->size()), 
        num_srcsrc( num_src * ( num_src - 1 ) / 2 ),
        coeffs(boost::extents[num_srcsrc+num_src][interp_order+1])
{
}

const InteractionBase::ResultArray &PulseInteraction::evaluate(
    const int time_idx)
{
  const double shift = 40.5;
  const double nPshift = 10/dt + ceil(shift);

  const double time = time_idx * dt;
  const double it = floor( time / dt - shift ) + 1;
  const double t = it - ( time / dt - shift );
  
  Interpolation::UniformLagrangeSet lagrange(t, interp_order, dt);
//  Interpolation::UniformLagrangeSet lagrange(time_idx, interp_order, dt);

  for(size_t dot = 0; dot < dots->size(); ++dot) {
    Eigen::Vector3d fld_pulse = Eigen::Vector3d::Zero();
    for(int i = 0; i <= interp_order; ++i){
      if ( (it - i) >= 0 && (it - i) <= nPshift ) 
        fld_pulse +=
          lagrange.evaluations[1][i] *
          (*pulse)( (*dots)[dot].position(), ( it - i ) * dt, 0, rotating );
        //time - i*dt, rotating);
    }
    results[dot] = fld_pulse.dot((*dots)[dot].dipole()) / hbar;
  }

  return results;
}

const InteractionBase::ResultArray &PulseInteraction::evaluatefld(
    const int time_idx, const int deriv_order)
{
  const double shift = 40.0;
  const double nPshift = 10/dt + ceil(shift);

  const double time = time_idx * dt;
  const double it = floor( time / dt - shift );
  const double t = ( time / dt - shift ) - it;
 
  Interpolation::UniformLagrangeSet lagrange(t, interp_order, dt);
//  Interpolation::UniformLagrangeSet lagrange(time_idx, interp_order, dt);

  for(size_t dot = 0; dot < dots->size(); ++dot) {
    Eigen::Vector3d fld_pulse = Eigen::Vector3d::Zero();
    for(int i = 0; i <= interp_order; ++i){
      if ( (it - i + 1) >= 0 && (it - i + 1) <= nPshift ) 
        fld_pulse +=
          lagrange.evaluations[deriv_order][i] * 
          (*pulse)( (*obss)[dot].position(), ( it - i + 1 ) * dt, 0, rotating );
      if (time_idx == 0)
        std::cout << "Deriv order: " << deriv_order << " " << i << " " << lagrange.evaluations[deriv_order][i] << std::endl;

    }
    results[dot] = fld_pulse.dot((*obss)[dot].dipole()) / hbar;
  }

//  if (deriv_order == 0)
//      std::cout << time_idx << " " << ((*pulse)( (*obss)[0].position(), time_idx * dt, rotating ))[0] << std::endl;
 

  return results;
}
