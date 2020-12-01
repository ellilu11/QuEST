#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "../src/common.h"
#include "../src/quantum_dot.h"
#include "../src/integrator/history.h"
#include "../src/integrator/integrator.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/interactions/AIM/aim_interaction.h"

using namespace std;

const double c0 = 300; // 299.792458
const double mu0 = 2.0133545e-04, hbar = 0.65821193;
const double omega = 2278.9013;
const double lambda = 2 * M_PI * c0 / omega;
const double prop_constant = // 1.00 / hbar;
                             mu0 / (4 * M_PI * hbar);

double gauss(double t, double mu, double sigsqr) {
  double val = exp(-std::pow(t - mu, 2) / (2.0 * sigsqr));
  return val;
}

Eigen::Vector2cd source(double t, double mu, double sigsqr){
  Eigen::Vector2cd vec(0, gauss(t, mu, sigsqr));
  return vec;
}

Eigen::Vector3d d0_source(double t, double mu, double sigsqr){
  Eigen::Vector3d vec(gauss(t, mu, sigsqr), 0, 0); 
  return vec;
}

Eigen::Vector3d d1_source(double t, double mu, double sigsqr){
  Eigen::Vector3d vec(-(t-mu) / sigsqr * gauss(t, mu, sigsqr), 0, 0); 
  return vec;
}

Eigen::Vector3d d2_source(double t, double mu, double sigsqr){
  Eigen::Vector3d vec((std::pow(t-mu, 2) - sigsqr) / pow(sigsqr, 2) *
          gauss(t, mu, sigsqr), 0, 0); 
  return vec;
}

Eigen::Vector3d analytic_EFIE_evaluate(Eigen::Vector3d &fld_d0,
                                       Eigen::Vector3d &fld_d1,
                                       Eigen::Vector3d &fld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){
  Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
  Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
  Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3.0 * rr;

  return -pow(c0, 2) * prop_constant * hbar * 
     (i3rr * fld_d0 / std::pow(dist, 3) +
      i3rr * fld_d1 / (c0 * std::pow(dist, 2)) +
      irr * fld_d2 / (std::pow(c0, 2) * dist));
}

Eigen::Vector3d analytic_MFIE_evaluate(Eigen::Vector3d &fld_d1,
                                       Eigen::Vector3d &fld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){

  Eigen::Vector3d rhat = dr / dr.norm();
  return -pow(c0, 2) * prop_constant * hbar * rhat.cross( 
     fld_d1 / (std::pow(dist, 2)) + fld_d2 / (c0 * dist) );
}

Eigen::Vector3cd analytic_RotatingEFIE_evaluate
																			(Eigen::Vector3d &efld_d0,
                                       Eigen::Vector3d &efld_d1,
                                       Eigen::Vector3d &efld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){
  Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
  Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
  Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3.0 * rr;

  return -pow(c0, 2) * prop_constant * hbar * exp( -iu*omega*dist/c0 ) *
     (1.0 * i3rr * efld_d0 / std::pow(dist, 3) +
      1.0 * i3rr * ( efld_d1 + iu*omega*efld_d0 ) / (c0 * std::pow(dist, 2)) +
      1.0 * irr * ( efld_d2 + 2.0*iu*omega*efld_d1 - pow(omega,2)*efld_d0 ) / (std::pow(c0, 2) * dist));
}

int main(int argc, char *argv[]){

		const int rotating = atoi(argv[1]);
    const int interp = atoi(argv[2]);
    const double tmax = 10;
    const int steps = 2000;
    const int window = 22;
    const double dt = tmax/(double)steps;

		cout << "  Frame: " << (rotating ? "Rotating" : "Fixed") << endl;
		cout << "  Interp order: " << interp << endl;
    cout << "  Using timestep: " << dt << endl;

    const double mu = tmax/2.0;
    const double sig = tmax/10.0;
    const double sigsqr = sig * sig;

    auto dots = std::make_shared<DotVector>(import_dots("dots/dots0.cfg"));
    const int ndots = (*dots).size();
    cout << "  Setting up " << ndots << " dots" << endl;
    const int nsrcs = ndots;
 
    auto obss = std::make_shared<DotVector>(import_dots("dots/obss_sphsurf0.cfg"));
    const int nobss = (*obss).size();
    cout << "  Setting up " << nobss << " observers" << endl;
 
    int min_time_to_keep = steps + window - 10; // for this test, just make the history array store all timedata
    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
        ndots, window, steps, min_time_to_keep, 2, 0); 
    history->fill(Eigen::Vector2cd::Zero());
    for (int i = -window; i < steps; ++i)
      for (int n = 0; n < nsrcs; ++n)
        history->set_value(n, i, 0) = source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm() / 2.0;
 
    cout << "  Setting up direct interaction" << endl;
		std::shared_ptr<InteractionBase> pairwise;    

		if (rotating) {
			Propagation::RotatingEFIE propagator(c0, prop_constant, omega, 0.0, 0.0);
    	pairwise = std::make_shared<DirectInteraction>(
        dots, history, propagator, interp, c0, dt, omega);
	
		} else {
			Propagation::EFIE<cmplx> propagator(c0, prop_constant, 0.0, 0.0);
    	pairwise = std::make_shared<DirectInteraction>(
        dots, history, propagator, interp, c0, dt);
		}

    cout << "  Setting up src-obs interaction" << endl;
    std::shared_ptr<InteractionBase> pairwise_fld;

    if (rotating) {
      Propagation::RotatingEFIE dyadic(c0, propagation_constant, omega, beta, 0.0);
      pairwise_fld = make_shared<DirectInteraction>(dots, obss, history, dyadic,
                                                    interpolation_order, c0, dt, omega);
    } else {
      Propagation::EFIE<cmplx> dyadic(c0, propagation_constant, beta, 0.0);
      pairwise_fld = make_shared<DirectInteraction>(dots, obss, history, dyadic,
                                                    interpolation_order, c0, dt); 
    }

    std::shared_ptr<InteractionBase> pairwise_bfld;

    if (rotating) {
      Propagation::RotatingMFIE dyadic_mfie(c0, propagation_constant, omega);
      pairwise_bfld = make_shared<DirectInteraction>(dots, obss, history, dyadic_mfie,
                                                     interpolation_order, c0, dt, omega);
    } else {
      Propagation::MFIE<cmplx> dyadic_mfie(c0, propagation_constant);
      pairwise_bfld = make_shared<DirectInteraction>(dots, obss, history, dyadic_mfie,
                                                     interpolation_order, c0, dt);
    }

    cout << "  Calculating and writing solutions" << endl;
 
    cmplx fld_dir, fld_anl;

    ofstream outfile, errfile;
    fldfile.open("out/test/fld.dat");
    fldfile << scientific << setprecision(15);
    fluxfile.open("out/test/flux.dat"); //ios:app
    fluxfile << scientific << setprecision(15);

    double err_dir = 0;
		double anl_sum = 0;

		std::clock_t start_time;
    start_time = std::clock();

    for (int step = 0; step < steps; ++step) {

			const double time = step*dt;
       
      set_dipole_of_dots( obss, Eigen::Vector3d(hbar, 0, 0) );
      auto efldx = pairwise_efld->evaluate_field(step); 
      auto bfldx = pairwise_bfld->evaluate_field(step, 1);
     
      set_dipole_of_dots( obss, Eigen::Vector3d(0, hbar, 0) );
      auto efldy = pairwise_efld->evaluate_field(step); 
      auto bfldy = pairwise_bfld->evaluate_field(step, 1);
      
      set_dipole_of_dots( obss, Eigen::Vector3d(0, 0, hbar) );
      auto efldz = pairwise_efld->evaluate_field(step); 
      auto bfldz = pairwise_bfld->evaluate_field(step, 1);

      double flux_dir = 0;
      double flux_anl = 0;

      for (int iobs = 0; iobs < nobss; ++idot) {
        Eigen::Vector3cd efld(efldx[iobs], efldy[iobs], efldz[iobs]);
        Eigen::Vector3cd bfld(bfldx[iobs], bfldy[iobs], bfldz[iobs]);

        for (int isrc = 0; isrc < ndots; ++isrc) {
          Eigen::Vector3d dr(separation( (*dots)[isrc], (*obss)[iobs] ));
          double dist = dr.norm();
          double delay = dist / c0;

          Eigen::Vector3d fld_d0 = fld_d0_source( step*dt, mu+delay, sigsqr );
          Eigen::Vector3d fld_d1 = fld_d1_source( step*dt, mu+delay, sigsqr );
          Eigen::Vector3d fld_d2 = fld_d2_source( step*dt, mu+delay, sigsqr );

					if (rotating) {
 	    			efld_anl += analytic_RotatingEFIE_evaluate(efld_d0, efld_d1, efld_d2, dr, c0, dist);

					} else {
     				efld_anl += analytic_EFIE_evaluate(fld_d0, fld_d1, fld_d2, dr, c0, dist);
      			bfld_anl += analytic_MFIE_evaluate(fld_d1, fld_d2, dr, c0, dist);;
          }
        }
 
        Eigen::Vector3d poynting_dir = (efld.real()).cross(bfld.real()) / mu0;
        Eigen::Vector3d poynting_anl = efld_anl.cross(bfld_anl) / mu0;

        Eigen::Vector3d pos = (*obss)[iobs].position();
        Eigen::Vector3d nhat = pos / pos.norm();
        double area_element = (*obss)[obs].frequency();
        
        flux_anl += area_element * poynting_anl.dot(nhat);
 
        err_dir += pow( fld_dir_abs - fld_anl_abs, 2 );
				anl_sum += fld_anl_abs;
      }
	  
      fluxfile << fld_dir << " " << flux_anl << std::endl;

    }

		const int ndata = steps*ndots;
		const double relerr_dir = sqrt(err_dir) / anl_sum;

    std::cout << "  Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s" << std::endl;
   
    std::cout << "  Direct err: " << relerr_dir << std::endl;
    cout << endl;

    errfile << interp << " " << relerr_dir << endl;
    outfile.close();
    errfile.close();
}
