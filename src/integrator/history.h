#ifndef HISTORY_H
#define HISTORY_H

#include <assert.h>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>
#include <chrono>
#include <random>
#include <string>
#include <limits>
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
  History(
    std::shared_ptr<const DotVector>,
    const int, const int, const int = -1, const int = 2, const int = 1, const int = 0,
    const double = 0, const double = 0);
 ~History();

  void fill(const soltype &);
  void initialize_past(const double);
  soltype &set_value(const int, const int, const int);
  soltype get_value(const int, const int, const int) const;
  soltype &set_gamma(const int, const int);
  soltype get_gamma(const int, const int) const;

  void write_step_to_file(const int);
 
  Integrator::soltype_array<soltype> zero_array();
  Integrator::soltype_array<soltype> read_history(const std::string);
  Integrator::soltype_array<soltype> mean_over_trials(const soltype_array<soltype>);
 
  int num_particles;
  int num_timesteps;
  int num_trials;
  int window;
  int task_idx;

  double omega_, c_;

 private:
  std::shared_ptr<const DotVector> dots;
  soltype_array<soltype> array_;
  soltype_array<soltype> gamma_;

  std::ofstream outfile, outfile2, outfile3, ofsdfile;
  std::ifstream ifphfile, ifrhfile;
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
Integrator::History<soltype>::History(std::shared_ptr<const DotVector> dots,
                                      const int window,
                                      const int num_timesteps,
                                      const int min_time_idx_ubound,  // rename!
                                      const int num_derivatives,
                                      const int num_trials,
                                      const int task_idx,
                                      const double omega, const double c)
    : dots(std::move(dots)), 
      num_particles(dots->size()),
      num_timesteps(num_timesteps), num_trials(num_trials),
      window(window), task_idx(task_idx),
      omega_(omega), c_(c)
{
  // TODO: need to make the minimum: max_transit_steps_between_dots +
  // interpolation_order
  int time_idx_ubound = std::max(min_time_idx_ubound, window) + 10;

  array_.resize(
      boost::extents[num_particles][time_idx_ubound][num_derivatives]);

  gamma_.resize(
      boost::extents[num_particles][time_idx_ubound][1]);

  outfile.open("./out/out" + std::to_string(task_idx) + ".dat");
  outfile << std::scientific << std::setprecision(18);

  outfile2.open("./out/gamma" + std::to_string(task_idx) + ".dat");
  outfile2 << std::scientific << std::setprecision(18);

  // outfile3.open("./rhoinit.out", std::ios_base::app);
  // outfile3 << std::scientific << std::setprecision(18);
  // ofsdfile.open("./seeds.out", std::ios_base::app);
  ifphfile.open("./phases.cfg");
  ifrhfile.open("./rhoinit.cfg");

}

template <class soltype>
Integrator::History<soltype>::~History()
{
  outfile.close();
  outfile2.close();
  outfile3.close();
  ifphfile.close();
  ifrhfile.close();
}

template <class soltype>
void Integrator::History<soltype>::fill(const soltype &val)
{
  std::fill(array_.data(), array_.data() + array_.num_elements(), val);
}

// FIXED INITIAL RHO
/*template <class soltype>
void Integrator::History<soltype>::initialize_past(const double sigma)
{
  Eigen::Vector2cd val(1,0);
  for(int n = 0; n < num_particles; ++n){
    if ( n < 50 ) val[0] = 0;

    for(int t = -window; t <= 0; ++t) 
      set_value(n, t, DERIV_0) = val; 
  }
}*/

// RANDOM INITIAL RHO USING TIME SEED
/*template <class soltype>
void Integrator::History<soltype>::initialize_past(const double sigma)
{
  const int seed = 
       std::chrono::system_clock::now().time_since_epoch().count();
  ofsdfile << seed << std::endl;
 */

  /*
  ifsdfile.seekg(std::ios::beg);
  for (int i = 0; i < task_idx; ++i)
    ifsdfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::string seedstr;
  ifsdfile >> seedstr;
  const int seed = std::stoi(seedstr);  
  ofsdfile << seed << std::endl;

  std::default_random_engine generator(seed);

  std::normal_distribution<double> rho_01r(0, sigma);
  std::normal_distribution<double> rho_01i(0, sigma);

   // cmplx rho_01(rho_01r(generator),rho_01i(generator));
}*/   

// RANDOM INITIAL RHO, NO ADDED PHASE
/*template <class soltype>
void Integrator::History<soltype>::initialize_past(const double sigma)
{
  cmplx rho_00(0,0);
  // cmplx rho_00(0.5,0);

  ifrhfile.seekg(std::ios::beg);
  for (int i = 0; i < task_idx*num_particles; ++i)
    ifrhfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  for(int n = 0; n < static_cast<int>(array_.shape()[PARTICLES]); ++n){
    std::string str, line;

    std::getline( ifrhfile, line );
    std::stringstream ss(line);    

    int i = 0;
    Eigen::Vector2d rho_01_vec;
    while( std::getline(ss, str, ' ') ){
      rho_01_vec[i] = std::stod(str);
      i++;
    }

    cmplx rho_01 = rho_01_vec[0]+iu*rho_01_vec[1];

    Eigen::Vector2cd val(rho_00, rho_01);
    // Eigen::Vector2cd val(rho_00, 0.5);

    for(int t = -window; t <= 0; ++t) 
      set_value(n, t, DERIV_0) = val; 
  }
}*/

// IDENTICAL RANDOM INITIAL RHO WITH ADDED PHASE
template <class soltype>
void Integrator::History<soltype>::initialize_past(const double sigma)
{

  // ===== GET INIT RHO =====
  ifrhfile.seekg(std::ios::beg);
  for (int i = 0; i < task_idx*num_particles; ++i)
    ifrhfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::string str, line;

  std::getline( ifrhfile, line );
  std::stringstream ss(line);    

  int i = 0;
  Eigen::Vector2d rho_01_vec;
  while( std::getline(ss, str, ' ') ){
    rho_01_vec[i] = std::stod(str);
    i++;
  }
  
  cmplx rho_00(0,0);
  cmplx rho_01(rho_01_vec[0],rho_01_vec[1]);

  // ===== GET PHASES =====
  // double phi = 2.0 * M_PI / (double)nsets;
  ifphfile.seekg(std::ios::beg);
  for (int i = 0; i < task_idx*num_particles; ++i)
    ifphfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  int nsets = 1;
  int npset = num_particles / nsets;
  std::vector<double> phis(nsets);
  for(int n = 0; n < nsets; ++n){
    std::string phstr;
    ifphfile >> phstr;
    double phi = std::stod(phstr);  
    std::cout << n << " " << phi << std::endl; 
    phis[n] = phi;
  }

  // ===== APPLY PHASES TO INIT RHO =====
  for(int n = 0; n < num_particles; ++n){
    Eigen::Vector2cd val(rho_00, rho_01);
    // int mod = n % nsets;
    // val[1] *= std::exp( iu * phi * (double)mod );
 
    // val[1] *= std::exp( iu * phis[n%nsets] );   
    val[1] *= std::exp( iu * phis[n/npset] );   
    
    for(int t = -window; t <= 0; ++t){ 
      if (n < 50)
        set_value(n, t, DERIV_0) = val; 
      else
        set_value(n, t, DERIV_0) = Eigen::Vector2cd(1,0);
    }
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
soltype Integrator::History<soltype>::get_gamma(const int particle_idx,
                                                const int time_idx) const
{
  int t = time_idx_in_array(time_idx);
  return gamma_[particle_idx][t][0];
}

template <class soltype>
soltype &Integrator::History<soltype>::set_gamma(const int particle_idx,
                                                 const int time_idx)
{
  int t = time_idx_in_array(time_idx);
  return gamma_[particle_idx][t][0];
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
  Eigen::Vector2cd gamma = Eigen::Vector2cd::Zero();
  for(int n = 0; n < num_particles; ++n){
    //outfile << prep_for_output(get_value(n, timestep, 0)) << " ";
    double rho01_real = (get_value(n, timestep, 0))[1].real();
    double rho01_imag = (get_value(n, timestep, 0))[1].imag();

    outfile // << timestep << " "
            << (get_value(n, timestep, 0))[0].real() << " "
            << rho01_real << " "
            << rho01_imag << " ";

    gamma += get_gamma(n, timestep);

  }

  outfile << "\n";

  outfile2 << gamma[0].real() << " "
           << gamma[1].real() << " " 
           << gamma[1].imag() << "\n";
 
  if(timestep > num_timesteps*0.90){
    outfile.flush();  // QUESTION: why is this needed?
    outfile2.flush();  // QUESTION: why is this needed?
  }

}

/*
template <class soltype>
void Integrator::History<soltype>::write_step_to_file(const int timestep)
{
  for(int n = 0; n < num_particles; ++n){
    //outfile << prep_for_output(get_value(n, timestep, 0)) << " ";
    double rho01_real = (get_value(n, timestep, 0))[1].real();
    double rho01_imag = (get_value(n, timestep, 0))[1].imag();

    outfile // << timestep << " "
            << (get_value(n, timestep, 0))[0].real() << " "
            << rho01_real << " "
            << rho01_imag << " ";

    Eigen::Vector2cd gamma = get_gamma(n, timestep);
    outfile2 << gamma[0].real() << " " 
             << gamma[1].real() << " " 
             << gamma[1].imag() << " ";

  }

  outfile << "\n";
  outfile2 << "\n";

  if(timestep > num_timesteps*0.90){
    outfile.flush();  // QUESTION: why is this needed?
    outfile2.flush();  // QUESTION: why is this needed?
  }

}
*/

template <class soltype>
Integrator::soltype_array<soltype> 
  Integrator::History<soltype>::zero_array()
{
  soltype_array<soltype> new_array_; 
  new_array_.resize(
      boost::extents[0][0][0]);
  // std::fill(new_array_.data(), new_array_.data() + new_array_.num_elements(), 
  //          Eigen::Vector2cd::Zero());

  return new_array_;
}

template <class soltype>
Integrator::soltype_array<soltype> 
  Integrator::History<soltype>::read_history( const std::string filestr )
{
  std::ifstream file;
  file.open(filestr);

  soltype_array<soltype> new_array_ = zero_array(); 

  int i = 0;
  std::string str, line;

  while( file.good() ){
    std::getline( file, line );
    std::stringstream ss(line);    

    while( std::getline(ss, str, ' ') ){
      int row = i / ( 3*num_particles );
      int col = i % ( 3*num_particles );

      int part = col / 3;
      int var = col % 3;

      int trial = row / num_timesteps;
      int step = row % num_timesteps;

      i++;

      double val = std::stod(str);
    
      switch(var) {
        case 0 :
          new_array_[part][step][trial] += Eigen::Vector2cd(val,0);
          break;
        case 1 :
          new_array_[part][step][trial] += Eigen::Vector2cd(0, val);
          break;
        case 2 :
          new_array_[part][step][trial] += Eigen::Vector2cd(0, iu*val);
          break;
      } 
    }
  }

  return new_array_;
}

template <class soltype>
Integrator::soltype_array<soltype> 
  Integrator::History<soltype>::mean_over_trials( const soltype_array<soltype> history_array_ )
{
  soltype_array<soltype> new_array_; 
  new_array_.resize(
      boost::extents[num_particles][num_timesteps][1]);

  for (int step = 0; step < num_timesteps; ++step){
    for (int part = 0; part < num_particles; ++part){
      for (int trial = 0; trial < num_trials; ++trial)
        new_array_[part][step][0] += history_array_[part][step][trial];
      new_array_[part][step][0] /= (double)num_trials;
      Eigen::Vector2cd rho_mean = new_array_[part][step][0];

    /*outfile3 << rho_mean[0].real() << " "
             << rho_mean[1].real() << " "
             << rho_mean[1].imag() << " ";
   */
    }
    // outfile3 << '\n';

    /*if(step > num_timesteps*0.90)
      outfile3.flush();  // QUESTION: why is this needed?
    */
  }

  return new_array_;
}


#endif
