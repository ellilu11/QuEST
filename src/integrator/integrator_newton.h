#ifndef INTEGRATOR_NEWTON_H
#define INTEGRATOR_NEWTON_H

#include <string>

#include "RHS/rhs.h"
#include "history.h"
#include "logging.h"
#include "math_utils.h"
#include "weights.h"

constexpr int DIM = 3;
// constexpr int NUM_ITERATIONS = 20;
constexpr double T1 = 100000.0;
constexpr double T2 = 200000.0;

template <class soltype>
Integrator::NewtonJacobian<soltype>::NewtonJacobian(
    const double dt,
    const double beta,
    const double omega,
    const int interp_order,
    const bool rotating,
    const std::shared_ptr<Integrator::History<soltype>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions
    )
    : num_solutions(history->array_.shape()[0]),
      time_idx_ubound(history->array_.index_bases()[1] +
                      history->array_.shape()[1]),
      dt(dt), beta(beta), omega(omega), interp_order(interp_order), rotating(rotating),
      history(std::move(history)),
      interactions(std::move(interactions)),
      y_vec(DIM*num_solutions),
      y_prev(DIM*num_solutions),
      b_vec(DIM*num_solutions),
      x_vec(DIM*num_solutions),
      rhs_vec(DIM*num_solutions),
      rhs_prev(DIM*num_solutions),
      rhs_J(DIM*num_solutions, DIM*num_solutions),
      coeffs(boost::extents[num_solutions*(num_solutions+1)/2][interp_order+1])

{
   for(int solution = 0; solution < num_solutions; ++solution)
      y_prev.segment<DIM>(DIM*solution) = cmplx2real( history->array_[solution][0][0] ); // assign initial conditions

}

template <class soltype>
void Integrator::NewtonJacobian<soltype>::solve(
    const log_level_t log_level)
{

 coeffs = (interactions[interactions.size()-1])->coefficients(); 

 for(int i = 0; i < num_solutions; ++i)
   cout << coeffs[i][0] << endl;

 ofstream Jfile("./outsr/J.dat");

 for(int step = 1; step < time_idx_ubound; ++step) {
    assert(0 <= step && step < time_idx_ubound);

    int niter = 0;
    y_vec = y_prev;
  
    if (step%10000 == 0) cout << step << " ";

    do {
      update_rhs(step);
      update_J(step);

      if (step%1000 == 1){
          Jfile << step << endl;
          Jfile << rhs_J << endl << endl;
            // (Eigen::MatrixXd::Identity(DIM*num_solutions,DIM*num_solutions) - dt*rhs_J) << endl;
      }

      evaluate(step);
      niter++;
 } while ( x_vec.norm() > EPS ); //&& niter < NUM_ITERATIONS);

  y_prev = y_vec;

   if (step%10000 == 0) cout << endl;
  }
}

/*
template <class soltype>
void Integrator::NewtonJacobian<soltype>::solve_step(const int step)
{
  assert(0 <= step && step < time_idx_ubound);

  y_vec = y_prev;
  
  if (step%1000 == 0) cout << step << " ";

  for(int m = 0; m < NUM_ITERATIONS; ++m) {
    update_rhs(step);
    update_J(step);
    evaluate(step);
  }

  y_prev = y_vec;

  if (step%1000 == 0) cout << endl;
}*/

template <class soltype>
void Integrator::NewtonJacobian<soltype>::update_rhs(const int step) 
{
    rhs_prev = rhs_vec;
    
    auto eval_and_sum =
        [step](const InteractionBase::ResultArray &r,
               const std::shared_ptr<InteractionBase> &interaction) {
            return r + interaction->evaluate(step);
        };
    auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

    rabi = std::accumulate(
        interactions.begin(), interactions.end(), nil, eval_and_sum);

//    std::cout << step << " " << rabi[0] << " " << rabi[1] << endl;

    // fill in func for all dots
    double f0, w, fr, fi, gr, gi, hr, hi;

    for(int solution = 0; solution < num_solutions; ++solution) {
        Eigen::VectorXd y_vec_seg = y_vec.segment<DIM>(DIM*solution);
        f0 = y_vec_seg[0]; 
        fr = y_vec_seg[1];
        fi = y_vec_seg[2];

        rhs_vec.segment<DIM>(DIM*solution) <<
          2.0 * ( fr * imag( rabi[solution] ) - fi * real( rabi[solution] ) ) + (1.0 - f0) / T1,
          (1.0 - 2.0 * f0) * imag( rabi[solution] ) - fr / T2 - (rotating ? 0.0 : (omega * fi)),
         -(1.0 - 2.0 * f0) * real( rabi[solution] ) - fi / T2 + (rotating ? 0.0 : (omega * fr)); 

    }
}

template <class soltype>
void Integrator::NewtonJacobian<soltype>::update_J(const int step) 
{
  double w, fr, fi;

  // ResultArray rabi_pulse = (interactions[0])->evaluate(step);
  // ResultArray rabi_direct = (interactions[1])->evaluate(step);

  ResultArray rabi_minus_self = rabi - (interactions[1])->evaluate(step);

  Eigen::MatrixXd rhs_J_pair(DIM*num_solutions,DIM*num_solutions);
  Eigen::MatrixXd rhs_J_self(DIM*num_solutions,DIM*num_solutions);

  int num_srcsrc = num_solutions * ( num_solutions - 1 ) / 2;

  for(int i = 0; i < num_solutions; ++i) {

    Eigen::VectorXd y_vec_seg = y_vec.segment<DIM>(DIM*i);
    w = 1.0 - 2.0 * y_vec_seg[0];
    fr = y_vec_seg[1];
    fi = y_vec_seg[2];

    // diagonal submatrices
    rhs_J_pair.block<DIM,DIM>(DIM*i,DIM*i) <<
     -1.0 / T1, 2.0 * imag( rabi_minus_self[i] ), -2.0 * real( rabi_minus_self[i] ),
     -2.0 * imag( rabi_minus_self[i] ), -1.0 / T2, -(rotating ? 0.0 : omega),
      2.0 * real( rabi_minus_self[i] ), rotating ? 0.0 : omega, -1.0 / T2; 

    double self_real = real( coeffs[num_srcsrc+i][0] );
    double self_imag = imag( coeffs[num_srcsrc+i][0] );

    rhs_J_self.block<DIM,DIM>(DIM*i,DIM*i) <<
      0.0, 4.0 * self_imag * fr, 4.0 * self_imag * fi,
     -2.0 * ( self_imag * fr + self_real * fi ),  self_imag * w, self_real * w,
     -2.0 * ( self_imag * fi - self_real * fr ), -self_real * w, self_imag * w;
    
    // off-diagonal submatrices
    for(int j = 0; j < i; ++j){
      int pair_idx = coord2idx(i,j); 

      double pair_real, pair_imag;
      
      for(int m = 0; m <= 0; ++m){
        pair_real += real( coeffs[pair_idx][m] );
        pair_imag += imag( coeffs[pair_idx][m] );
      }

      rhs_J_pair.block<DIM,DIM>(DIM*i,DIM*j) <<
        0.0, 2.0 * ( pair_imag * fr - pair_real * fi ), 2.0 * ( pair_real * fr + pair_imag * fi ),
        0.0,  pair_imag * w, pair_real * w,
        0.0, -pair_real * w, pair_imag * w;

      rhs_J_pair.block<DIM,DIM>(DIM*j,DIM*i) = 
        rhs_J_pair.block<DIM,DIM>(DIM*i,DIM*j);

      rhs_J_self.block<DIM,DIM>(DIM*i,DIM*j) = Eigen::MatrixXd::Zero(DIM,DIM);
      rhs_J_self.block<DIM,DIM>(DIM*j,DIM*i) = Eigen::MatrixXd::Zero(DIM,DIM);

     }
  }
  rhs_J = rhs_J_pair + rhs_J_self;
}

template <class soltype>
void Integrator::NewtonJacobian<soltype>::evaluate(const int step)
{
    b_vec = -y_vec + dt * rhs_vec + y_prev;
//    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> factor(rhs_J);
    Eigen::MatrixXd rhs_A = Eigen::MatrixXd::Identity(DIM*num_solutions,DIM*num_solutions) - dt * rhs_J;
    Eigen::FullPivLU<Eigen::MatrixXd> factor(rhs_A);
    x_vec = factor.solve( b_vec );

    // get condition number
    Eigen::BDCSVD<Eigen::MatrixXd> svd(rhs_A);
    Eigen::VectorXd svalues = svd.singularValues();
    if (step%10000 == 0)
        cout << svalues[0] / svalues[DIM*num_solutions-1] << " ";

    /*solver.analyzePattern(rhs_J);
    solver.factorize(rhs_J);
    x_vec = solver.solve( -y_vec + dt * rhs_vec + y_prev );
*/
    // x_vec = rhs_J * ( -y_vec + dt * rhs_vec + y_prev );
    y_vec += x_vec;

    for(int solution=0; solution < num_solutions; ++solution)
        history->array_[solution][step][0] = real2cmplx( y_vec.segment<DIM>(DIM*solution) );
  
} 

template <class soltype>
Eigen::VectorXd Integrator::NewtonJacobian<soltype>::cmplx2real(Eigen::Vector2cd cmplx_vec)
{
    Eigen::VectorXd vec(DIM);
    vec << real(cmplx_vec[0]), 
           real(cmplx_vec[1]), imag(cmplx_vec[1]);

    return vec;
}


template <class soltype>
Eigen::Vector2cd Integrator::NewtonJacobian<soltype>::real2cmplx(Eigen::VectorXd vec)
{
    Eigen::Vector2cd cmplx_vec;
    cmplx_vec[0] = vec[0];
    cmplx_vec[1] = vec[1] + iu * vec[2];
    
    return cmplx_vec;
}

template <class soltype>
int Integrator::NewtonJacobian<soltype>::coord2idx(int row, int col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

template <class soltype>
void Integrator::NewtonJacobian<soltype>::log_percentage_complete(
    const int step) const
{
  if(step % (time_idx_ubound / 10) == 0) {
    std::cout << "\t" << static_cast<int>(10.0 * step / time_idx_ubound)
              << std::endl;
  }
}

#endif
