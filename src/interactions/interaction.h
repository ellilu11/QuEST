#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <memory>

#include "../quantum_dot.h"

class InteractionBase {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  InteractionBase(
    const std::shared_ptr<const DotVector> dots, 
    const std::shared_ptr<const DotVector> obss, 
    const double dt)
      : dots(std::move(dots)),
        obss(std::move(obss)),
        results(obss ? obss->size() : (dots ? dots->size() : 0) ),
        past_terms_of_convolution(dots ? dots->size() : 0), 
        dt(dt){};
  const cmplx &operator[](const int i) const { return results[i]; }
  virtual const ResultArray &evaluate(const int) = 0;
  virtual const ResultArray &evaluate_present(const int) = 0;
  virtual const ResultArray &evaluate_field(const int, const bool=0) = 0;
  virtual ~InteractionBase(){};

 protected:
  static int coord2idxsq(int row, int col, int rowlen)
  {
    return row*rowlen + col;
  } 

  std::shared_ptr<const DotVector> dots, obss;
  ResultArray results;
  ResultArray past_terms_of_convolution;
  double dt;
};

#endif
