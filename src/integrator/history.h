#ifndef HISTORY_H
#define HISTORY_H

#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace Integrator {
  template <class soltype>
  class History;

  template <class soltype>
  using soltype_array = boost::multi_array<soltype, 3>;

  inline namespace history_enums {
    enum DIMENSION { PARTICLES, TIMES, DERIVATIVES };
    enum ORDER { DERIV_0, DERIV_1 };
  }
}

template <class soltype>
class Integrator::History {
 public:
  History(const int, const int, const int, const int = 2);
  soltype_array<soltype> array_;
  // boost::multi_array<soltype, 3> array_;

  void fill(const soltype &);
  void initialize_past(const soltype &, const int);

 private:
};

template <class soltype>
Integrator::History<soltype>::History(const int num_particles,
                                      const int window,
                                      const int num_timesteps,
                                      const int num_derivatives)
    : array_(boost::extents[num_particles][
          typename soltype_array<soltype>::extent_range(-window, num_timesteps)]
                          [num_derivatives])

{
}

template <class soltype>
void Integrator::History<soltype>::fill(const soltype &val)
{
  std::fill(array_.data(), array_.data() + array_.num_elements(), val);
}

/*
template <class soltype>
void integrator::history<soltype>::initialize_past(const std::string fname)
{
    int nparts = static_cast<int>(array_.shape()[PARTICLES]);
    
}*/

template <class soltype>
void Integrator::History<soltype>::initialize_past(const soltype &val, const int nparts)
{
  /*for(int n = 0; n < static_cast<int>(array_.shape()[PARTICLES]); ++n) {
    for(int t = array_.index_bases()[TIMES]; t <= 0; ++t) {
      array_[n][t][DERIV_0] = val;
    }
  }*/
  // int nparts = static_cast<int>(array_.shape()[PARTICLES]);
  
  for(int n = 0; n < nparts; ++n) {
    for(int t = array_.index_bases()[TIMES]; t <= 0; ++t) {
        array_[n][t][DERIV_0] = val;
      // if (n == nparts-1) array_[n][t][DERIV_0] = val2;
      // else array_[n][t][DERIV_0] = val1;
    }
  }
  

}

#endif
