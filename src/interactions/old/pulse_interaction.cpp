#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(const std::shared_ptr<const DotVector> dots,
                                   const std::shared_ptr<const DotVector> obss,
                                   const std::shared_ptr<const Pulse> pulse,
                                   const double hbar,
                                   const double dt,
                                   const bool rotating)
    : InteractionBase(dots, obss, dt), pulse(std::move(pulse)), hbar(hbar), rotating(rotating)
{

}

const InteractionBase::ResultArray &PulseInteraction::evaluate(
    const int time_idx)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < dots->size(); ++i) {
    results[i] =
        (*pulse)((*dots)[i].position(), time, rotating).dot((*dots)[i].dipole()) / hbar;

//    if (time_idx == 200) std::cout << (*pulse)((*dots)[i].position(), time, rotating) << std::endl;

    //((*pulse)((*dots)[i].position(), time, rotating).dot((*dots)[i].dipole()) + iu *
    // (*pulse)((*dots)[i].position(), time, rotating).dot((*dots)[i].dipole_imag())) / hbar;
}

  return results;
}

/*
const InteractionBase::ResultArray &PulseInteraction::evaluatefld(
    const int time_idx)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < dots->size(); ++i)
    flds[i] = (*pulse)((*dots)[i].position(), time, rotating).norm();

  return flds;
}*/
