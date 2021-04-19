#ifndef HISTORY_H
#define HISTORY_H

#include <assert.h>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include "../quantum_dot.h"

namespace Integrator {
  template <class soltype>
  class History;

  template <class soltype>
  using soltype_array = boost::multi_array<soltype, 3>;

  inline namespace history_enums {
    enum DIMENSION { PARTICLES, TIMES, DERIVATIVES };
    enum ORDER { DERIV_0, DERIV_1 };
  }  // namespace history_enums
}  // namespace Integrator

template <class soltype>
class Integrator::History {
 public:
  // History(const int, const int, const int, const int = 2);
  History(const int, const int, const int, const int = -1, const int = 2, const int = 0);
  ~History();

  void fill(const soltype &);
  void initialize_past(const soltype &);
  soltype &set_value(const int, const int, const int);
  soltype get_value(const int, const int, const int) const;

  void write_step_to_file(const int);

  int num_particles;
  int num_timesteps;
  int window;

 private:
  soltype_array<soltype> array_;
  std::ofstream outfile;
  int time_idx_in_array(const int) const;
  // friend int time_idx_in_array(const int);
  template <class B = soltype>
  typename std::enable_if<std::is_same<B, Eigen::Vector2cd>::value,
                          Eigen::RowVector2cd>::type
      prep_for_output(Eigen::Vector2cd);

  template <class B = soltype>
  typename std::enable_if<!std::is_same<B, Eigen::Vector2cd>::value,
                          soltype>::type prep_for_output(soltype);
};

template <class soltype>
Integrator::History<soltype>::History(const int num_particles,
                                      const int window,
                                      const int num_timesteps,
                                      const int min_time_idx_ubound,  // rename!
                                      const int num_derivatives,
                                      const int task_idx)
    : num_particles(num_particles), num_timesteps(num_timesteps), window(window)
{
  // TODO: need to make the minimum: max_transit_steps_between_dots +
  // interpolation_order
  int time_idx_ubound = std::max(min_time_idx_ubound, window) + 10;

  array_.resize(
      boost::extents[num_particles][time_idx_ubound][num_derivatives]);

  outfile.open("./out/out" + std::to_string(task_idx) + ".dat");
  outfile << std::scientific << std::setprecision(15);

}

template <class soltype>
Integrator::History<soltype>::~History()
{
  outfile.close();
}

template <class soltype>
void Integrator::History<soltype>::fill(const soltype &val)
{
  std::fill(array_.data(), array_.data() + array_.num_elements(), val);
}

template <class soltype>
void Integrator::History<soltype>::initialize_past(const soltype &val)
{
  for(int n = 0; n < static_cast<int>(array_.shape()[PARTICLES]); ++n) {
    for(int t = -window; t <= 0; ++t) { set_value(n, t, DERIV_0) = val; }
  }
}

template <class soltype>
int Integrator::History<soltype>::time_idx_in_array(const int time_idx) const
{
  int time_base = 0;  // static_cast<int>(array_.index_bases()[1]);
  int time_idx_ubound =
      static_cast<int>(array_.index_bases()[1] + array_.shape()[1]);

  return (time_idx >= time_base) ? 
    //(time_idx-1) % time_idx_ubound + 1;
    time_idx % time_idx_ubound
    : time_idx_ubound - abs(time_idx);
}

template <class soltype>
soltype Integrator::History<soltype>::get_value(const int particle_idx,
                                                const int time_idx,
                                                const int derivative_idx) const
{
  int t = time_idx_in_array(time_idx);
  return array_[particle_idx][t][derivative_idx];
}

template <class soltype>
soltype &Integrator::History<soltype>::set_value(const int particle_idx,
                                                 const int time_idx,
                                                 const int derivative_idx)
{
  int t = time_idx_in_array(time_idx);
  return array_[particle_idx][t][derivative_idx];
}

template <class soltype>
template <class B>
typename std::enable_if<std::is_same<B, Eigen::Vector2cd>::value,
                        Eigen::RowVector2cd>::type
Integrator::History<soltype>::prep_for_output(Eigen::Vector2cd element)
{
  return element.transpose();
}

template <class soltype>
template <class B>
typename std::enable_if<!std::is_same<B, Eigen::Vector2cd>::value,
                        soltype>::type
Integrator::History<soltype>::prep_for_output(soltype element)
{
  return element;
}

template <class soltype>
void Integrator::History<soltype>::write_step_to_file(const int timestep)
{
  for(int n = 0; n < num_particles; ++n){
    //outfile << prep_for_output(get_value(n, timestep, 0)) << " ";
    double rho00 = (get_value(n, timestep, 0))[0].real(); 
    double rho01_real = (get_value(n, timestep, 0))[1].real();
    double rho01_imag = (get_value(n, timestep, 0))[1].imag();
    // double rho01_abs = sqrt( pow(rho01_real,2) + pow(rho01_imag,2) );

    outfile // << timestep << " "
            << rho00 << " "
            << rho01_real << " "
            << rho01_imag << " "
            << (get_value(n, timestep, 1))[0].real() << " "
            << (get_value(n, timestep, 1))[1].real() << " "
            << (get_value(n, timestep, 1))[1].imag() << " ";

            // << rho01_abs << " ";
  }

  outfile << "\n";
  if(timestep > num_timesteps*0.90)
    outfile.flush();  // QUESTION: why is this needed?
}

#endif
