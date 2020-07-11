#ifndef HISTORY_INTERACTION_H
#define HISTORY_INTERACTION_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include "../integrator/history.h"
#include "../lagrange_set.h"
#include "../math_utils.h"
#include "../quantum_dot.h"
#include "green_function.h"
#include "interaction.h"

class HistoryInteraction : public InteractionBase {
 public:
  HistoryInteraction(
      std::shared_ptr<const DotVector> dots,
      std::shared_ptr<const DotVector> obss,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history, // feed Gaussian pulse into here
      const int interp_order,
      const double c0,
      const double dt)
      : InteractionBase(std::move(dots), std::move(obss), dt),
        history(std::move(history)),
        interp_order(interp_order),
        c0(c0){};

//    const ResultArray &evaluate(const int);

 protected:
  std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history;
  int interp_order;
  const double c0;

//  int coord2idx(int, int);
//  std::pair<int, int> idx2coord(const int);
//  std::pair<int, double> split_Double(const double);
 
};

#endif
