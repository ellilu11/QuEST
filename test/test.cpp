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

using namespace std;

const double c0 = 299.792458;
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
     (1.0 * i3rr * fld_d0 / std::pow(dist, 3) +
      1.0 * i3rr * fld_d1 / (c0 * std::pow(dist, 2)) +
      1.0 * irr * fld_d2 / (std::pow(c0, 2) * dist));
}

Eigen::Vector3d analytic_MFIE_evaluate(Eigen::Vector3d &fld_d1,
                                       Eigen::Vector3d &fld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){

  Eigen::Vector3d rhat = dr / dr.norm();
  return -prop_constant * hbar * rhat.cross( 
     1.0 * fld_d1 / (std::pow(dist, 2)) + 1.0 * fld_d2 / ( c0 * dist ) );
}

Eigen::Vector3cd analytic_RotatingEFIE_evaluate
																			(Eigen::Vector3d &fld_d0,
                                       Eigen::Vector3d &fld_d1,
                                       Eigen::Vector3d &fld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){
  Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
  Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
  Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3.0 * rr;

  return -pow(c0, 2) * prop_constant * hbar * exp( -iu*omega*dist/c0 ) *
     (1.0 * i3rr * fld_d0 / std::pow(dist, 3) +
      1.0 * i3rr * ( fld_d1 + iu*omega*fld_d0 ) / (c0 * std::pow(dist, 2)) +
      1.0 * irr * ( fld_d2 + 2.0*iu*omega*fld_d1 - pow(omega,2)*fld_d0 ) / (std::pow(c0, 2) * dist));
}

Eigen::Vector3cd analytic_RotatingMFIE_evaluate
  																		(Eigen::Vector3d &fld_d0,
                                       Eigen::Vector3d &fld_d1,
                                       Eigen::Vector3d &fld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){
  Eigen::Vector3d rhat = dr / dr.norm();
  return -prop_constant * hbar * exp( -iu*omega*dist/c0 ) * rhat.cross( 
      1.0 * ( fld_d1 + iu*omega*fld_d0 ) / (std::pow(dist, 2)) + 
      1.0 * ( fld_d2 + 2.0*iu*omega*fld_d1 - pow(omega,2)*fld_d0) / ( c0 * dist ) );
}

int main(int argc, char *argv[]){

		const int rotating = atoi(argv[1]);
    const int interp = atoi(argv[2]);
    const double tmax = 10;
    const int steps = 1000;
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
    cout << "  Setting up " << ndots << " sources" << endl;
    const int nsrcs = ndots;
 
    auto obss = std::make_shared<DotVector>(import_dots("dots/obss_sphsurf.cfg"));
    const int nobss = (*obss).size();
    cout << "  Setting up " << nobss << " observers" << endl;
    Eigen::Vector3d pos = (*obss)[0].position();
    const double r0 = pos.norm();
    cout << "  Observer radius: " << r0 << endl;
 
    int min_time_to_keep = steps + window - 10; // for this test, just make the history array store all timedata
    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
        ndots, window, steps, min_time_to_keep, 2, 0); 
    history->fill(Eigen::Vector2cd::Zero());
    for (int i = -window; i < steps; ++i)
      for (int n = 0; n < nsrcs; ++n)
        history->set_value(n, i, 0) = source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm() / 2.0;
 
/*    cout << "  Setting up direct interaction" << endl;
		std::shared_ptr<InteractionBase> pairwise;    

		if (rotating) {
			Propagation::RotatingEFIE propagator(c0, prop_constant, omega, 0.0, 0.0);
    	pairwise = std::make_shared<DirectInteraction>(
        dots, nullptr, history, propagator, interp, c0, dt, omega);
	
		} else {
			Propagation::EFIE<cmplx> propagator(c0, prop_constant, 0.0, 0.0);
    	pairwise = std::make_shared<DirectInteraction>(
        dots, nullptr, history, propagator, interp, c0, dt);
		}
*/
    cout << "  Setting up src-obs interaction" << endl;
    std::shared_ptr<InteractionBase> pairwise_fld;

    if (rotating) {
      Propagation::RotatingEFIE prop_efie(c0, prop_constant, omega, 0.0, 0.0);
      pairwise_fld = make_shared<DirectInteraction>(dots, obss, history, prop_efie,
                                                    interp, c0, dt, omega);
    } else {
      Propagation::EFIE<cmplx> prop_efie(c0, prop_constant, 0.0, 0.0);
      pairwise_fld = make_shared<DirectInteraction>(dots, obss, history, prop_efie,
                                                    interp, c0, dt); 
    }

    std::shared_ptr<InteractionBase> pairwise_bfld;

    if (rotating) {
      Propagation::RotatingMFIE prop_mfie(c0, prop_constant, omega);
      pairwise_bfld = make_shared<DirectInteraction>(dots, obss, history, prop_mfie,
                                                     interp, c0, dt, omega);
    } else {
      Propagation::MFIE<cmplx> prop_mfie(c0, prop_constant);
      pairwise_bfld = make_shared<DirectInteraction>(dots, obss, history, prop_mfie,
                                                     interp, c0, dt);
    }

    cout << "  Calculating and writing solutions" << endl;
 
    ofstream fldfile, fluxfile;
    fldfile.open("out/test/fld.dat");
    fldfile << scientific << setprecision(15);
    fluxfile.open("out/test/flux.dat"); //ios:app
    fluxfile << scientific << setprecision(15);

    Eigen::Vector3d poynting_dir, poynting_anl;

    double err_fld = 0;
    double err_flux = 0;
		
    std::clock_t start_time;
    start_time = std::clock();

    for (int step = 0; step < steps; ++step) {

			const double time = step*dt;
       
      set_dipole_of_dots( obss, Eigen::Vector3d(hbar, 0, 0) );
      auto efldx = pairwise_fld->evaluate_field(step); 
      auto bfldx = pairwise_bfld->evaluate_field(step, 1);
     
      set_dipole_of_dots( obss, Eigen::Vector3d(0, hbar, 0) );
      auto efldy = pairwise_fld->evaluate_field(step); 
      auto bfldy = pairwise_bfld->evaluate_field(step, 1);
      
      set_dipole_of_dots( obss, Eigen::Vector3d(0, 0, hbar) );
      auto efldz = pairwise_fld->evaluate_field(step); 
      auto bfldz = pairwise_bfld->evaluate_field(step, 1);

      double flux_dir = 0;
      double flux_anl = 0;

      for (int iobs = 0; iobs < nobss; ++iobs) {
        Eigen::Vector3d pos = (*obss)[iobs].position();
        Eigen::Vector3d nhat = pos / pos.norm();
        double area_element = (*obss)[iobs].frequency();

        Eigen::Vector3cd efld_dir(efldx[iobs], efldy[iobs], efldz[iobs]);
        Eigen::Vector3cd bfld_dir(bfldx[iobs], bfldy[iobs], bfldz[iobs]);
        Eigen::Vector3cd efld_anl = Eigen::Vector3cd::Zero();
        Eigen::Vector3cd bfld_anl = Eigen::Vector3cd::Zero();

        for (int isrc = 0; isrc < ndots; ++isrc) {
          Eigen::Vector3d dr(separation( (*dots)[isrc], (*obss)[iobs] ));
          double dist = dr.norm();
          double delay = dist / c0;

          Eigen::Vector3d fld_d0 = d0_source( step*dt, mu, sigsqr );
          Eigen::Vector3d fld_d1 = d1_source( step*dt, mu, sigsqr );
          Eigen::Vector3d fld_d2 = d2_source( step*dt, mu, sigsqr );

					if (rotating) {
 	    			efld_anl += analytic_RotatingEFIE_evaluate(fld_d0, fld_d1, fld_d2, dr, c0, dist);
      			bfld_anl += analytic_RotatingMFIE_evaluate(fld_d0, fld_d1, fld_d2, dr, c0, dist);
					} else {
     				efld_anl += analytic_EFIE_evaluate(fld_d0, fld_d1, fld_d2, dr, c0, dist);
      			bfld_anl += analytic_MFIE_evaluate(fld_d1, fld_d2, dr, c0, dist);
          }
        }

        if (rotating) { 
          poynting_dir = (efld_dir.real()).cross(bfld_dir.real()) / mu0;
          poynting_anl = (efld_anl.real()).cross(bfld_anl.real()) / mu0;
        } else {
          poynting_dir = ( efld_dir.cross(bfld_dir.conjugate()) + 
                           efld_dir.cross(bfld_dir) * exp(2.0*iu*omega*time) ).real() / ( 2.0*mu0 );
          poynting_anl = ( efld_anl.cross(bfld_anl.conjugate()) + 
                           efld_anl.cross(bfld_anl) * exp(2.0*iu*omega*time) ).real() / ( 2.0*mu0 );
        }

        double flux_per_area_dir = poynting_dir.dot(nhat);     
        double flux_per_area_anl = poynting_anl.dot(nhat);     

        if ( step == steps / 2 ) {
          fldfile << flux_per_area_dir << " " << flux_per_area_anl << endl;
          err_fld += pow( flux_per_area_dir - flux_per_area_anl, 2 );
        }

        flux_dir += area_element * flux_per_area_dir;
        flux_anl += area_element * flux_per_area_anl;
 
        // err_dir += pow( fld_dir_abs - fld_anl_abs, 2 );
      }

      Eigen::Vector3d fld_d2 = d2_source( step*dt, mu, sigsqr );
      double flux_anl_ff = mu0 * fld_d2.squaredNorm() / ( 6.0 * M_PI * c0 );
	
      fluxfile << flux_dir << " " << flux_anl << " " << flux_anl_ff << " " << flux_anl_ff / flux_anl << std::endl;
      err_flux += pow( flux_dir - flux_anl, 2 );

    }

		const int ndata = steps*nobss;
		err_fld = sqrt(err_fld) / nobss;
    err_flux = sqrt(err_flux) / steps;

    cout << "  Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s" << endl;
   
    cout << "  Flux per area err: " << err_fld << endl;
    cout << "  Total flux err: " << err_flux << endl;

    // errfile << interp << " " << relerr_dir << endl;
    fldfile.close();
    fluxfile.close();
}
