#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "../pulse.h"
#include "interaction.h"

class PulseInteraction : public InteractionBase {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector>,
                   const std::shared_ptr<const Pulse>, 
                   const int, 
                   const double, const double, const double, 
                   const bool);
  const ResultArray &evaluate(const int);
  const ResultArray &evaluate_present_field(const int time_idx)
  {
    return evaluate(time_idx);
  }
  boost::multi_array<cmplx, 2> &coefficients() final { return coeffs; }

//  const Eigen::MatrixXd &evaluateJ(const int, const int);

 private:
  int num_src, num_interactions; 
  boost::multi_array<cmplx, 2> coeffs;
  std::shared_ptr<const Pulse> pulse;
  const double hbar;
  const bool rotating;
};

#endif
