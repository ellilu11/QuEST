#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const DotVector> obss,
    const std::shared_ptr<const Pulse> pulse,
    const int interp_order,
    const double dt,
    const double hbar,
    const double omega,
    const bool rotating)
    : InteractionBase(dots, obss, dt), pulse(std::move(pulse)), hbar(hbar), omega(omega), rotating(rotating)
{
}

const InteractionBase::ResultArray &PulseInteraction::evaluate(
    const int time_idx)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < dots->size(); ++i)
    results[i] =
        (*pulse)((*dots)[i].position(), time, 0, rotating).dot((*dots)[i].dipole()) / hbar;

  return results;
}

const InteractionBase::ResultArray &PulseInteraction::evaluate_field(
    const int time_idx, const bool flag)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < obss->size(); ++i){
    auto pulse0 = (*pulse)((*obss)[i].position(), time, 0, rotating);
    if ( flag ) {
      Eigen::Vector3d wavevector = pulse->wavevec();
      results[i] =
        (wavevector.cross(pulse0)).dot((*obss)[i].dipole()) / ( omega * hbar );
    } else 
      results[i] =
        pulse0.dot((*obss)[i].dipole()) / hbar;
  } 
  return results;
}

