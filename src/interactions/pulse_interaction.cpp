#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(const std::shared_ptr<const DotVector> dots,
                                   const std::shared_ptr<const DotVector> obss,
                                   const std::shared_ptr<const Pulse> pulse,
                                   const int interp_order,
                                   const double c0,
                                   const double dt,
                                   const double hbar,
                                   const bool rotating)
    : InteractionBase(dots, obss, dt), pulse(std::move(pulse)), hbar(hbar), rotating(rotating),
        num_src((this->dots)->size()), 
        num_srcsrc( num_src * ( num_src - 1 ) / 2 ),
        coeffs(boost::extents[num_srcsrc+num_src][interp_order+1])
{
}

const InteractionBase::ResultArray &PulseInteraction::evaluate(
    const int time_idx)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < dots->size(); ++i) {
    results[i] =
        (*pulse)((*dots)[i].position(), time, 0, rotating).dot((*dots)[i].dipole()) / hbar;

//    if (time_idx == 200) std::cout << (*pulse)((*dots)[i].position(), time, rotating) << std::endl;

    //((*pulse)((*dots)[i].position(), time, rotating).dot((*dots)[i].dipole()) + iu *
    // (*pulse)((*dots)[i].position(), time, rotating).dot((*dots)[i].dipole_imag())) / hbar;
  }

  return results;
}

const InteractionBase::ResultArray &PulseInteraction::evaluatefld(
    const int time_idx, const int deriv_order)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < obss->size(); ++i) {
    results[i] =
        (*pulse)((*obss)[i].position(), time, 0, rotating).dot((*obss)[i].dipole()) / hbar;

  }

  return results;
}
