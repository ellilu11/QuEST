#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "../pulse.h"
#include "interaction.h"

class PulseInteraction : public InteractionBase {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector>,
                   const std::shared_ptr<const DotVector>,
                   const std::shared_ptr<const Pulse>, 
                   const int, 
                   const double, const double, const double, 
                   const bool, const bool);
  const ResultArray &evaluate(const int);
  const ResultArray &evaluate_present(const int time_idx)
  {
    return evaluate(time_idx);
  }
  const ResultArray &evaluate_field(const int, const bool=0);
 
 private:
  std::shared_ptr<const Pulse> pulse;
  const double hbar, omega;
  const bool rotating, rwa;
};

#endif
