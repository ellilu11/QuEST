#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <memory>

#include "../common.h"
#include "../quantum_dot.h"

class InteractionBase {
 public:
  typedef Eigen::Array<Eigen::Vector3cd, Eigen::Dynamic, 1> ResultArray;

  InteractionBase(
    const std::shared_ptr<const DotVector> dots, 
    const std::shared_ptr<const DotVector> obss, 
    const double dt)
      : dots(std::move(dots)),
        obss(std::move(obss)),
        num_src((this->dots)->size()),
        num_obs((this->obss) ? (this->obss)->size() : 0),
        num_interactions( num_src * (num_src-1) / 2 ),
        num_srcobs( num_src * num_obs ),
        results(obss ? obss->size() : (dots ? num_src*num_src : 0) ),
        past_terms_of_convolution(dots ? num_src*num_src : 0), 
        dt(dt){};
  const Eigen::Vector3cd &operator[](const int i) const { return results[i]; }
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
  const int num_src, num_obs, num_interactions, num_srcobs;
  ResultArray results;
  ResultArray past_terms_of_convolution;
  double dt;
};

#endif
