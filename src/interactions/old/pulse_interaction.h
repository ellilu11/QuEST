#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "../pulse.h"
#include "interaction.h"

class PulseInteraction : public InteractionBase {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector>,
                   const std::shared_ptr<const DotVector>,
                   const std::shared_ptr<const Pulse>, const double, const double, const bool);
  virtual const ResultArray &evaluate(const int);

 private:
  std::shared_ptr<const Pulse> pulse;
  const double hbar;
  const bool rotating;
};

#endif
