#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "../src/quantum_dot.h"
#include "../src/integrator/history.h"
#include "../src/integrator/integrator.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/interactions/AIM/aim_interaction.h"

using namespace std;

const double c0 = 299.792458, mu0 = 2.0133545e-04, hbar = 0.65821193;
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

Eigen::Vector3d efld_d0_source(double t, double mu, double sigsqr){
  Eigen::Vector3d vec(gauss(t, mu, sigsqr), 0, 0); 
  return vec;
}

Eigen::Vector3d efld_d1_source(double t, double mu, double sigsqr){
  Eigen::Vector3d vec(-(t-mu) / sigsqr * gauss(t, mu, sigsqr), 0, 0); 
  return vec;
}

Eigen::Vector3d efld_d2_source(double t, double mu, double sigsqr){
  Eigen::Vector3d vec((std::pow(t-mu, 2) - sigsqr) / pow(sigsqr, 2) *
          gauss(t, mu, sigsqr), 0, 0); 
  return vec;
}

Eigen::Vector3d analytic_Laplace_evaluate(Eigen::Vector3d &efld_d0,
                                       Eigen::Vector3d &dr,
                                       double dist){
  return prop_constant * hbar * efld_d0 / dist;
}

/*Eigen::Vector3d analytic_Helmholtz_evaluate(Eigen::Vector3d &efld_d0,
                                       Eigen::Vector3d &dr,
                                       double dist){
  return prop_constant * hbar * efld_d0 * exp( iu*omega/c0*dist) / dist;
}*/

Eigen::Vector3d analytic_EFIE_evaluate(Eigen::Vector3d &efld_d0,
                                       Eigen::Vector3d &efld_d1,
                                       Eigen::Vector3d &efld_d2,
                                       Eigen::Vector3d &dr,
                                       double c0, double dist){
    Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
    Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
    Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3 * rr;

    return -pow(c0, 2) * prop_constant * hbar * 
       (i3rr * efld_d0 / std::pow(dist, 3) +
        i3rr * efld_d1 / (c0 * std::pow(dist, 2)) +
        irr * efld_d2 / (std::pow(c0, 2) * dist));
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
       (i3rr * efld_d0 / std::pow(dist, 3) +
        i3rr * ( efld_d1 + iu*omega*efld_d0 ) / (c0 * std::pow(dist, 2)) +
        irr * ( efld_d2 + 2.0*iu*omega*efld_d1 - pow(omega,2)*efld_d0 ) / (std::pow(c0, 2) * dist));
}

int main(int argc, char *argv[]){

    const int interp = 5;
    const double tmax = 10;
    const int steps = 1000;
    const int window = 22;
    const double dt = tmax/(double)steps;
    cout << "  Interp order: " << interp << endl;
    const double mu = tmax/2.0;
    const double sig = tmax/10.0;
    const double sigsqr = sig * sig;
    // const double fmax = 6.0/(2.0*M_PI*sqrt(sigsqr));

    auto dots = std::make_shared<DotVector>(import_dots("dots_aimtest.cfg"));
    const int ndots = (*dots).size();
    cout << "  Setting up " << ndots << " dots" << endl;
    const int nsrcs = ndots / 2;
    const int nobss = ndots - nsrcs;
 
    int min_time_to_keep = steps + window - 10; // for this test, just make the history array store all timedata
    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
        ndots, window, steps, min_time_to_keep, 2, 0); 
    history->fill(Eigen::Vector2cd::Zero());
    for (int i = -window; i < steps; ++i) {
      for (int n = 0; n < nsrcs; ++n){
        history->set_value(n, i, 0) = source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm() / 2.0;
//        cout << i << " " << history->get_value(n,i,0)[1] << " ";
      }
//      cout << endl;
    }
 
    cout << "  Setting up AIM grid" << endl;
    const int border = 1;
    const int n_ds = 10;
    const double ds = 5.0e-2 * lambda * pow( 2, atoi(argv[1])/(double)n_ds );
    Eigen::Vector3d dsvec(ds, ds, ds); 
    std::cout << "    ds/lambda: " << ds/lambda << std::endl;
 
    // const double n_h = 10;
    const double h = 5.0e-1 * ds;
    // const double h = 5.0e-4 * ds * pow( 10, atoi(argv[1])/(double)n_h );
    int expansion = atoi(argv[2]);
    AIM::Grid grid(dsvec, expansion, h, *dots);
    int transit_steps = grid.max_transit_steps(c0, dt) + interp;
  
    cout << "  Setting up direct interaction" << endl;
    // Propagation::Laplace<cmplx> propagator(prop_constant);
    // Propagation::EFIE<cmplx> propagator(c0, prop_constant, 0.0, 0.0);
		Propagation::RotatingEFIE propagator(c0, prop_constant, omega, 0.0, 0.0);
    auto pairwise_dir = std::make_shared<DirectInteraction>(
        dots, history, propagator, interp, c0, dt, 
				omega);

    std::shared_ptr<InteractionBase> pairwise_aim, pairwise_fdtd;
    cout << "  Setting up AIM interaction" << endl;
    cout << "    AIM expansion order: " << expansion << std::endl;
    pairwise_aim = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt, 0,
				AIM::Expansions::RotatingEFIE(transit_steps, c0, dt, omega),
        // AIM::Expansions::EFIE(transit_steps, c0, dt), 
        AIM::Expansions::Zero(transit_steps),
        // AIM::Normalization::Laplace(prop_constant)
        AIM::Normalization::Helmholtz(omega/c0, prop_constant),
				omega
        );
 
    // expansion -= 2;
    // AIM::Grid grid2(dsvec, expansion, h, *dots);
    // transit_steps = grid2.max_transit_steps(c0, dt) + interp;
    cout << "  Setting up AIM-FDTD interaction" << endl;
    cout << "    AIM expansion order: " << expansion << std::endl;
    cout << "    h/ds: " << h/ds << std::endl;
    clock_t start_time;
		pairwise_fdtd = pairwise_aim;

    start_time = std::clock();
/*    pairwise_fdtd = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt, h,
        AIM::Expansions::EFIE_TimeDeriv2(transit_steps, c0, dt), 
        AIM::Expansions::EFIE_Retardation(transit_steps, c0),
        AIM::Normalization::Laplace(prop_constant)
        );*/
    std::cout << "  Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s" << std::endl;
 
    cout << "  Calculating and writing solutions" << endl;
 
    cmplx fld_dir, fld_aim, fld_fdtd, fld_anl;

    ofstream outfile, errfile;
    outfile.open("out/aimtest/out.dat");
    outfile << scientific << setprecision(15);
    errfile.open("out/aimtest/err.dat", ios::app);
    errfile << scientific << setprecision(15);

    double err_dir, err_aim, err_fdtd, err_fdtd_aim;
    double anl_sum, aim_sum;

    start_time = std::clock();

    for (int step = 0; step < steps; ++step) {

			const double time = step*dt;

      const InteractionBase::ResultArray array_dir = pairwise_dir->evaluate(step);
      const InteractionBase::ResultArray array_aim = pairwise_aim->evaluate(step); 
      const InteractionBase::ResultArray array_fdtd = pairwise_fdtd->evaluate(step); 
    
      for (int iobs = 0; iobs < nobss; ++iobs) {
        fld_dir = array_dir[nsrcs+iobs];
        fld_aim = array_aim[nsrcs+iobs];
        fld_fdtd = array_fdtd[nsrcs+iobs];
        fld_anl = 0;

        for (int isrc = 0; isrc < nsrcs; ++isrc) {
          // if ( isrc == iobs ) continue;
          Eigen::Vector3d dr(separation( (*dots)[nsrcs+iobs], (*dots)[isrc] ));
          double dist = dr.norm();
          double delay = dist / c0;

          Eigen::Vector3d efld_d0 = efld_d0_source( time, mu+delay, sigsqr );
          Eigen::Vector3d efld_d1 = efld_d1_source( time, mu+delay, sigsqr );
          Eigen::Vector3d efld_d2 = efld_d2_source( time, mu+delay, sigsqr );

          // fld_anl += analytic_EFIE_evaluate(efld_d0, efld_d1, efld_d2, dr, c0, dist).dot(
          //                     (*dots)[nsrcs+iobs].dipole() ) / hbar;
          fld_anl += real( analytic_RotatingEFIE_evaluate(efld_d0, efld_d1, efld_d2, dr, c0, dist).dot(
                               (*dots)[nsrcs+iobs].dipole() ) / hbar * exp( iu*omega*time ) ) 
											* exp( -iu*omega*time );
          // fld_anl += analytic_Laplace_evaluate(efld_d0, dr, dist).dot(
          //                      (*dots)[iobs].dipole() ) / hbar;
        }

        double fld_dir_abs = abs(fld_dir);
        double fld_aim_abs = abs(fld_aim);
        double fld_fdtd_abs = abs(fld_fdtd);
        double fld_anl_abs = abs(fld_anl);

        outfile << fld_anl_abs << " " << fld_dir_abs << " " << fld_aim_abs << " " << fld_fdtd_abs << " ";
        // errfile << fld_fdtd_abs / fld_aim_abs << " ";

        err_dir += pow( fld_dir_abs - fld_anl_abs, 2 );
        err_aim += pow( fld_aim_abs - fld_anl_abs, 2);
        err_fdtd += pow( fld_fdtd_abs - fld_anl_abs, 2);
        err_fdtd_aim += pow( fld_fdtd_abs - fld_aim_abs, 2 );
        anl_sum += fld_anl_abs;
        aim_sum += fld_aim_abs;
      }
      outfile << endl;
    }

    const int ndata = steps*nobss;
    const double relerr_dir = sqrt(err_dir) / anl_sum;
    const double relerr_aim = sqrt(err_aim) / anl_sum;
    const double relerr_fdtd = sqrt(err_fdtd) / anl_sum;
    const double relerr_fdtd_aim = sqrt(err_fdtd_aim) / aim_sum;

    std::cout << "  Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s" << std::endl;
   
    std::cout << "  Direct err: " << relerr_dir << std::endl;
    std::cout << "  AIM err: " << relerr_aim << std::endl;
    std::cout << "  FDAIM err: " << relerr_fdtd << std::endl;
    std::cout << "  FDAIM-AIM err: " << relerr_fdtd_aim << std::endl;
    cout << endl;

    errfile << expansion
				// << " " << h/ds 
				<< " " << ds/lambda
        << " " << relerr_dir
        << " " << relerr_aim
        << " " << relerr_fdtd
        << " " << relerr_fdtd_aim << std::endl;
    outfile.close();
    errfile.close();
}
