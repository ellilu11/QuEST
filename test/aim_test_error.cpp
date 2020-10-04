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
//#include "../src/interactions/AIM/expansion.h"
//#include "../src/interactions/AIM/grid.h"
//#include "../src/interactions/AIM/normalization.h"

using namespace std;

const double c0 = 299.792458, mu0 = 2.0133545e-04, hbar = 0.65821193;
const double omega = 2278.9013;
const double lambda = 2 * M_PI * c0 / omega;
const double prop_constant = // 1.00 / hbar;
                             mu0 / (4 * M_PI * hbar);
 
Eigen::Vector2cd source(double t, double mu, double sigsqr){
    return Eigen::Vector2cd(0, exp(-std::pow(t - mu, 2) / (2.0 * sigsqr)));
}

Eigen::Vector3d efld_d0_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(0, exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0); 
}

Eigen::Vector3d efld_d1_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(0, -(t - mu - delay) / sigsqr * exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0); 
}

Eigen::Vector3d efld_d2_source(double t, double mu, double sigsqr, double delay){
    return Eigen::Vector3d(0, (std::pow(t - mu - delay, 2) - sigsqr) / pow(sigsqr, 2) *
        exp(-std::pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0); 
}

Eigen::Vector3d analytic_EFIE_interaction(Eigen::Vector3d &efld_d0,
                                     Eigen::Vector3d &efld_d1,
                                     Eigen::Vector3d &efld_d2,
                                     Eigen::Vector3d &dr,
                                     double c0, double dist){
    Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
    Eigen::Matrix3d irr = Eigen::Matrix3d::Identity() - rr;
    Eigen::Matrix3d i3rr = Eigen::Matrix3d::Identity() - 3 * rr;

    return - pow(c0, 2) * prop_constant * hbar * 
       (i3rr * efld_d0 / std::pow(dist, 3) +
        i3rr * efld_d1 / (c0 * std::pow(dist, 2)) +
        irr * efld_d2 / (std::pow(c0, 2) * dist));
}

int main(int argc, char *argv[]){

    const bool fdtd = 0;

    const int interp = 5;
    const int tmax = 10;
    const int steps = 1000;
    const double dt = tmax/steps;
    const double mu = tmax/2.0;
    const double sig = tmax/10.0;
    const double sigsqr = sig * sig;
    // const double fmax = 6.0/(2.0*M_PI*sqrt(sigsqr));

    // set up dots and history of dots
    auto dots = std::make_shared<DotVector>(import_dots("./dots0.cfg"));
    const int ndots = (*dots).size();
    std::cout << "  Running with " << ndots << " dots" << std::endl;

    int min_time_to_keep = 
      max_transit_steps_between_dots(dots, c0, dt) + interp; 
    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
        ndots, 22, steps, min_time_to_keep, 2, 0); 
    history->fill(Eigen::Vector2cd::Zero());

    const int nsrcs = ndots/2;
    const int nobss = ndots - nsrcs;
    for (int n = 0; n < nsrcs; ++n)
        for (int i = -22; i < steps; ++i) 
            history->set_value(n, i, 0) = source(i * dt, mu, sigsqr) / (*dots)[n].dipole().norm();

    // set up direct interaction
    Propagation::EFIE<cmplx> propagator(c0, prop_constant, 0.0, 0.0);
    auto pairwise_dir = std::make_shared<DirectInteraction>(
        dots, history, propagator, interp, c0, dt, 0.0, 0);

    // set up AIM interaction
    const int expansion = atoi(argv[1]);
    const int border = 1;
    const double ds_base = 0.050 * lambda;
    const int ds_n = 10;

    const double ds = ds_base * pow( 2, atoi(argv[2])/(double)ds_n );
    Eigen::Vector3d dsvec(ds, ds, ds); 
    AIM::Grid grid(dsvec, expansion, *dots);
    const int transit_steps = grid.max_transit_steps(c0, dt) + interp;
    const double h = 0.5*ds;

    std::shared_ptr<InteractionBase> pairwise_aim;
    //if ( fdtd )
      pairwise_aim = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt, h,
        AIM::Expansions::EFIE_TimeDeriv2(transit_steps, c0, dt), 
        AIM::Expansions::EFIE_Retardation(transit_steps, c0),
        AIM::Normalization::Laplace(prop_constant)
        );
    /*else
      pairwise_aim = std::make_shared<AIM::Interaction>(
        dots, history, propagator, dsvec, interp, expansion, border, c0, dt, h,
        AIM::Expansions::EFIE(transit_steps, c0, dt), 
        NULL,
        AIM::Normalization::Laplace(prop_constant)
        );*/
 
    std::cout << "  AIM expansion order: " << expansion << std::endl;
    std::cout << "  ds/lambda: " << ds/lambda << std::endl;

    // calculate and output solutions
//    std::vector<Eigen::Vector3cd> fld_anl(nobss);
//    std::vector<Eigen::Vector3cd> fld_dir(nobss);
//    std::vector<Eigen::Vector3cd> fld_aim(nobss);
 
    cmplx fld_dir, fld_aim, fld_anl;

    std::ofstream outfile1, outfile2, outfile3;
    outfile1.open("out/aimtest/out_dir.dat");
    outfile1 << std::scientific << std::setprecision(15);
    outfile2.open("out/aimtest/out_aim.dat");
    outfile2 << std::scientific << std::setprecision(15);
    outfile3.open("outtest/error.dat", std::ios_base::app);
    outfile3 << std::scientific << std::setprecision(15);

    double err_dir, err_aim, err_aimdir;

    clock_t start_time;
    start_time = std::clock();

    for (int i = 0; i < steps; ++i) {
 
      const InteractionBase::ResultArray array_dir = pairwise_dir->evaluate(i);
      const InteractionBase::ResultArray array_aim = pairwise_aim->evaluate(i); 

      for (int iobs = 0; iobs < nobss; ++iobs) {
        fld_dir = array_dir[nsrcs+iobs];
        fld_aim = array_aim[nsrcs+iobs];

        for (int isrc = 0; isrc < nsrcs; ++isrc) {
          Eigen::Vector3d dr(separation( (*dots)[nsrcs+iobs], (*dots)[isrc] ));
          double dist = dr.norm();
          double delay = dist / c0;

          Eigen::Vector3d efld_d0 = efld_d0_source( i * dt, mu, sigsqr, delay );
          Eigen::Vector3d efld_d1 = efld_d1_source( i * dt, mu, sigsqr, delay );
          Eigen::Vector3d efld_d2 = efld_d2_source( i * dt, mu, sigsqr, delay );

          fld_anl += analytic_EFIE_interaction(efld_d0, efld_d1, efld_d2, dr, c0, dist).dot(
                                (*dots)[nsrcs+iobs].dipole() ) / hbar;
        }

        double fld_dir_abs = abs(fld_dir);
        double fld_aim_abs = abs(fld_aim);
        double fld_anl_abs = abs(fld_anl);

        outfile1 << real(fld_dir) << " " << imag(fld_dir) << fld_dir_abs << " ";
        outfile2 << real(fld_aim) << " " << imag(fld_aim) << fld_aim_abs << " ";

        err_dir += pow( fld_dir_abs - fld_anl_abs, 2 );
        err_aim += pow( fld_aim_abs - fld_anl_abs, 2 );
        err_aimdir += pow( fld_aim_abs - fld_dir_abs, 2 );

      }
    }

    std::cout << "  Elapsed time: " << (std::clock() - start_time) / (double) CLOCKS_PER_SEC << "s" << std::endl;
   
    std::cout << "  Direct interaction err: " << sqrt(err_dir / (steps*nobss)) << std::endl;
    std::cout << "  AIM interaction err: " << sqrt(err_aim / (steps*nobss)) << std::endl;
    std::cout << "  AIM-Direct interaction err: " << sqrt(err_aimdir / (steps*nobss)) << std::endl;
/*    outfile3 << expansion << " " << ds/lambda 
        << " " << sqrt(normdiff_dir / (steps*ntrgs)) 
        << " " << sqrt(normdiff_aim / (steps*ntrgs)) 
        << " " << sqrt(normdiff_aimdir / (steps*ntrgs)) << std::endl;
*/
    outfile1.close();
    outfile2.close();    
    outfile3.close();

}
