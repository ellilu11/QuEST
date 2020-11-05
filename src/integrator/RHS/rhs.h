#ifndef RHS_H
#define RHS_H

#include <memory>
#include "../history.h"
#include "../../common.h"

namespace Integrator {
  template <class soltype>
  class RHS;
}

template <class soltype>
class Integrator::RHS {
 public:
	typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  RHS(const double dt, const std::shared_ptr<History<soltype>> history)
      : dt(dt), history(std::move(history)){};
  virtual const ResultArray evaluate(const int) const = 0;
  virtual void evaluate_present(const int) const = 0;

protected:
  double dt;
  std::shared_ptr<History<soltype>> history;
};

#endif
