#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const DotVector> obss,
    const std::shared_ptr<const Pulse> pulse,
    const int interp_order,
    const double c0,
    const double dt,
    const double hbar,
    const bool rotating)
    : InteractionBase(dots, obss, dt), pulse(std::move(pulse)), hbar(hbar), rotating(rotating)
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
    const int time_idx)
{
  const double time = time_idx * dt;

  for(size_t i = 0; i < obss->size(); ++i)
    results[i] =
        (*pulse)((*obss)[i].position(), time, 0, rotating).dot((*obss)[i].dipole()) / hbar;

  return results;
}

