#include <Eigen/Dense>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>

#include "../src/configuration.h"
#include "../src/integrator/history.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/math_utils.h"
#include "../src/quantum_dot.cpp"
#include "../src/quantum_dot.h"

// namespace po = boost::program_options;

const double c0 = 299.792458, mu0 = 2.0133545e-04, hbar = 0.65821193;
const double omega = 2278.9013;
const double prop_constant =  // 1.00 / hbar;
    mu0 / (4 * M_PI * hbar);

Eigen::Vector2cd source(double t, double mu, double sigsqr)
{
  return Eigen::Vector2cd(0, exp(-std::pow(t - mu, 2) / (2.0 * sigsqr)));
}

Eigen::Vector3d efld_d0_source(double t, double mu, double sigsqr, double delay)
{
  return Eigen::Vector3d(exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0,
                         0);
}

Eigen::Vector3d efld_d1_source(double t, double mu, double sigsqr, double delay)
{
  return Eigen::Vector3d(-(t - mu - delay) / sigsqr *
                             exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)),
                         0, 0);
}

Eigen::Vector3d efld_d2_source(double t, double mu, double sigsqr, double delay)
{
  return Eigen::Vector3d((std::pow(t - mu - delay, 2) - sigsqr) /
                             pow(sigsqr, 2) *
                             exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)),
                         0, 0);
}

Eigen::Vector3d analytic_EFIE_interaction(Eigen::Vector3d &efld_d0,
                                          Eigen::Vector3d &efld_d1,
                                          Eigen::Vector3d &efld_d2,
                                          Eigen::Vector3d &dr,
                                          double c0,
                                          double dist)
{
  Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
  Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
  Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3 * rr;

  return -pow(c0, 2) * prop_constant * hbar *
         (i3rr * efld_d0 / std::pow(dist, 3) +
          i3rr * efld_d1 / (c0 * std::pow(dist, 2)) +
          irr * efld_d2 / (std::pow(c0, 2) * dist));
}

std::vector<std::complex<double>> analytic_evaluate(
    std::shared_ptr<DotVector> dots, int i, double dt, double mu, double sigsqr)
{
  int ndots = (*dots).size();
  std::vector<std::complex<double>> fld_anlytc(ndots);
  double dist, delay;

  for(int itrg = 0; itrg < ndots; itrg++) {
    for(int isrc = 0; isrc < ndots; ++isrc) {
      if(itrg != isrc) {
        Eigen::Vector3d dr(separation((*dots)[itrg], (*dots)[isrc]));
        dist = dr.norm();
        delay = dist / c0;

        Eigen::Vector3d efld_d0 = 2.0 * efld_d0_source(i * dt, mu, sigsqr, delay);
        Eigen::Vector3d efld_d1 = 2.0 * efld_d1_source(i * dt, mu, sigsqr, delay);
        Eigen::Vector3d efld_d2 = 2.0 * efld_d2_source(i * dt, mu, sigsqr, delay);

        fld_anlytc[itrg] +=
            analytic_EFIE_interaction(efld_d0, efld_d1, efld_d2, dr, c0, dist)
                .dot((*dots)[itrg].dipole()) /
            hbar;

      }
    }
  }
  return fld_anlytc;
}

int main(int argc, char *argv[])
{
  const int interp = 4;
  const int steps = 10;  // atoi(argv[1]);
  const double dt = 0.1;

  const double tmax = steps * dt;
  const double mu = tmax / 2.0;
  const double sig = tmax / 12.0;
  const double sigsqr = sig * sig;

  const double lambda = 2 * M_PI * c0 / omega;

  // == set up dots and history ==
  std::cout << "initializing dot vector..." << std::endl;
  auto dots = std::make_shared<DotVector>(import_dots("./dots0.cfg"));
  const int ndots = (*dots).size();
  std::cout << "Running with " << ndots << " dots" << std::endl;

  auto time_floor = max_transit_steps_between_dots(dots, c0, dt);
  std::cout << "max timesteps between dots: " << time_floor << std::endl;

  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      ndots, 22, steps, time_floor, 1);
  history->fill(Eigen::Vector2cd::Zero());

  // initializing history
  for(int n = 0; n < ndots; ++n)
    for(int i = -22; i < 0; ++i)
      // if ( n != OBS )
      history->set_value(n, i, 0) =
          source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm();

  // == prepare simulation ==
  // set up propagator
  auto propagator = Propagation::EFIE<cmplx>(c0, prop_constant, 0.0, 0.0);

  // calculate analytic and direct solution
  auto direct_interaction = std::make_shared<DirectInteraction>(
      dots, history, propagator, interp, c0, dt, omega, 0);

  std::vector<std::vector<std::complex<double>>> fld_anlytc(
      steps, std::vector<std::complex<double>>(ndots));
  std::vector<std::vector<std::complex<double>>> fld_dir(
      steps, std::vector<std::complex<double>>(ndots));
  std::vector<std::vector<std::complex<double>>> fld_diff(
      steps, std::vector<std::complex<double>>(ndots));

  for(int i = 0; i < steps; ++i) {
    
    // analytic electric field at each dot
    fld_anlytc[i] = analytic_evaluate(dots, i, dt, mu, sigsqr);
    // fix the value of the dots for this timestep
    for(int n = 0; n < ndots; n++)
      history->set_value(n, i, 0) =
          source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm();

    // calculate electric field at each dot
    const InteractionBase::ResultArray array_dir =
        direct_interaction->evaluate(i);

    for(int itrg = 0; itrg < ndots; itrg++) {
      fld_dir[i][itrg] = array_dir[itrg];
    }
     // std::cout << fld_dir[i][0] << " " << fld_dir[i][1] << std::endl;

 }
  std::ofstream analytic_out("analytic.dat");
  std::ofstream calculated_out("calculated.dat");

  analytic_out << std::scientific << std::setprecision(15);
  calculated_out << std::scientific << std::setprecision(15);
  double total_diff{0};
  double max_diff(0);
  for(int i = 0; i < steps; i++) {
    for(int dot = 0; dot < ndots; dot++) {
      analytic_out << fld_anlytc[i][dot].real() << " ";
      calculated_out << fld_dir[i][dot].real() << " ";

      max_diff =
          std::max(std::abs(fld_anlytc[i][dot].real() - fld_dir[i][dot].real()),
                   max_diff);
      total_diff +=
          std::abs(fld_anlytc[i][dot].real() - fld_dir[i][dot].real());
    }
    analytic_out << "\n";
    calculated_out << "\n";
  }

  std::cout << "max difference: " << max_diff << std::endl;
  std::cout << "total difference: " << total_diff << std::endl;

  return 0;
}
