#ifndef HISTORY_H
#define HISTORY_H

#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include "../quantum_dot.h"
#include "../math_utils.h"
//#include "RHS/bloch_rhs.h"

namespace Integrator {
  template <class soltype>
  class History;

  template <class soltype>
  using soltype_array = boost::multi_array<soltype, 2>;

  inline namespace history_enums {
    enum DIMENSION { TIMES, DERIVATIVES };
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
  // void initialize_past(const soltype &, const soltype &);
  soltype &set_value(const int, const int);
  soltype get_value(const int, const int) const;

  void write_step_to_file(const int);

  int num_particles;
  int num_timesteps;
  int window;
  int time_idx_ubound;

 private:
  soltype_array<soltype> array_;
  std::ofstream outfile, outfile2, outfile3;
  int time_idx_in_array(const int) const;
  // friend int time_idx_in_array(const int);
  /*template <class B = soltype>
  typename std::enable_if<std::is_same<B, matrix_cmplx>::value,
                          Eigen::RowVector2cd>::type
      prep_for_output(Eigen::Vector2cd);

  template <class B = soltype>
  typename std::enable_if<!std::is_same<B, matrix_cmplx>::value,
                          soltype>::type prep_for_output(soltype);
  */
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
  time_idx_ubound = std::max(min_time_idx_ubound, window) + 10;

  array_.resize(
      boost::extents[time_idx_ubound][num_derivatives]);

  outfile.open("./out/out" + std::to_string(task_idx) + ".dat");
  outfile << std::scientific << std::setprecision(18);

  outfile2.open("./out/eig" + std::to_string(task_idx) + ".dat");
  outfile2 << std::scientific << std::setprecision(18);

  outfile3.open("./out/pur" + std::to_string(task_idx) + ".dat");
  outfile3 << std::scientific << std::setprecision(18);


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
  for(int t = -window; t <= 0; ++t) { 
    // std::cout << time_idx_in_array(t) << " ";
    set_value(t, DERIV_0) = val; 
  }
}

template <class soltype>
int Integrator::History<soltype>::time_idx_in_array(const int time_idx) const
{
  int time_base = 0;  // static_cast<int>(array_.index_bases()[1]);
  /*int time_idx_ubound =
      static_cast<int>(array_.index_bases()[1] + array_.shape()[1]);
  */
  return (time_idx >= time_base) ? 
    //(time_idx-1) % time_idx_ubound + 1;
    time_idx % time_idx_ubound
    : time_idx_ubound - abs(time_idx);
}

template <class soltype>
soltype Integrator::History<soltype>::get_value(const int time_idx,
                                                const int derivative_idx) const
{
  int t = time_idx_in_array(time_idx);
  return array_[t][derivative_idx];
}

template <class soltype>
soltype &Integrator::History<soltype>::set_value(const int time_idx,
                                                 const int derivative_idx)
{
  int t = time_idx_in_array(time_idx);
  return array_[t][derivative_idx];
}

/*
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
*/

template <class soltype>
void Integrator::History<soltype>::write_step_to_file(const int timestep)
{

  for (int i = 0; i < num_particles; ++i){
   Eigen::MatrixXcd out = traceOverDims(
                          get_value(timestep, 0), 
                          num_particles, i ); 

   outfile << (out(0,0)).real() << " "
            << (out(0,1)).real() << " "
            << (out(0,1)).imag() << " ";
  }
  outfile << "\n";
 
  if(timestep > num_timesteps*0.90)
    outfile.flush();  // QUESTION: why is this needed?
 

  // ==== Purity test ==========================
  cmplx trc = get_value(timestep, 0).trace();
  cmplx pur = (get_value(timestep, 0)*get_value(timestep, 0)).trace();
  // outfile3 << (rho.trace()).real() << " " 
  //         << ((rho*rho).trace()).real() << " ";
  outfile3 << trc.real() << " " 
           << pur.real() << " ";
  outfile3 << "\n";

  if(timestep > num_timesteps*0.90)
    outfile3.flush();  // QUESTION: why is this needed?
 
  // ==== Separability test ======================
  const int pow_num = pow(2, num_particles);
  const int pow_half = pow_num / 2;
 
  Eigen::Matrix4cd rho = 
    // traceOverDim(
      // traceOverDim(
        get_value(timestep, 0);
        // num_particles, 2);          
      // num_particles-1, 2);

  // partial matrix transposition
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      (rho.block<2,2>(2*i,2*j)).transposeInPlace();
  
  // solve for Eigenvalues
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  ces.compute( rho );
  // std::cout << rho << std::endl << std::endl;
  // std::cout << ces.eigenvalues() << std::endl << std::endl;

  for (int dim = 0; dim < pow_num; ++dim)
    outfile2 << (ces.eigenvalues()[dim]).real() << " "
             << (ces.eigenvalues()[dim]).imag() << " ";
  outfile2 << "\n";
 
  if(timestep > num_timesteps*0.90)
    outfile2.flush();  // QUESTION: why is this needed?
 
}

#endif
