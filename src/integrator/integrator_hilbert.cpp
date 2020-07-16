#include <string>
#include "integrator.h"

using namespace std;

constexpr int DIM = 3;
constexpr double T1 = 1000.0;
constexpr double T2 = 2000.0;

template <class soltype>
Integrator::HilbertIntegrator<soltype>::HilbertIntegrator(
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
      interactions(std::move(interactions))

{
   for(int solution = 0; solution < num_solutions; ++solution)
      y_prev.segment<DIM>(DIM*solution) = cmplx2real( history->array_[solution][0][0] ); // assign initial conditions

}

template <class soltype>
void Integrator::Hilbert<soltype>::solve(
    const log_level_t log_level)
{

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

