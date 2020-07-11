#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <memory>

#include "../quantum_dot.h"

class InteractionBase {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;
//  typedef Eigen::Array<Eigen::Vector3cd, Eigen::Dynamic, 1> ResultArrayVec;

  InteractionBase(
    const std::shared_ptr<const DotVector> dots, 
    const std::shared_ptr<const DotVector> obss, const double dt)
      : dots(std::move(dots)), 
        obss(std::move(obss)),
        results(obss ? obss->size() : (dots ? dots->size() : 0) ), 
        dt(dt){};
  const cmplx &operator[](const int i) const { return results[i]; }
  virtual const ResultArray &evaluate(const int) = 0;
  virtual const ResultArray &evaluatefld(const int, const int) = 0;
  void modifyObs(std::shared_ptr<const DotVector> newobss){ obss = newobss; }
  virtual ~InteractionBase(){};

  std::shared_ptr<const DotVector> obss;

  virtual boost::multi_array<cmplx, 2> &coefficients() = 0;
//  boost::multi_array<cmplx, 2> coeffs;
 
 protected:
  std::shared_ptr<const DotVector> dots;
  ResultArray results;
  double dt;
};

#endif
