#include <complex>
#include <Eigen/Dense>
//#include <boost/test/unit_test.hpp>
//#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <ctime>

#include "../src/configuration.h"
#include "../src/quantum_dot.h"
#include "../src/quantum_dot.cpp"
#include "../src/integrator/history.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/math_utils.h"

using namespace std;

const double c0 = 299.7924,
             mu0 = 2.0133545e-04,
             hbar = 0.65821193;
const double prop_constant = mu0 / (4 * M_PI * hbar);
const double beta = 1.0;
const double omega = 2278.9013;
const double k0 = omega / c0;
const double lambda = 2 * M_PI / k0;
const double dist0 = 0.1 * omega / c0; 

const double dt = 1e-3;
const int steps = 10240;
const int interp = 3;
const double tmax = steps*dt;

const double mu = tmax/2.0;
const double sigma = tmax/12.0;
const double sigsqr = sigma * sigma;

// assume dipoles along x
Eigen::Vector3d P_d0( double t, double delay ){
    if ( 1 )//t <= (mu + delay) )
        return Eigen::Vector3d( exp(-pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0, 0 );
    else 
        return Eigen::Vector3d( 1.0, 0, 0 );
}

Eigen::Vector3d P_d1( double t, double delay ){
    if ( 1 )//t <= (mu + delay) )
        return Eigen::Vector3d( -(t - mu - delay) / sigsqr * exp(-pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0, 0 );
    else
        return Eigen::Vector3d::Zero();
}

Eigen::Vector3d P_d2( double t, double delay ){
    if ( 1 ) //t <= (mu + delay) )
        return Eigen::Vector3d( 
                ( pow( t - mu - delay, 2) - sigsqr ) / pow(sigsqr, 2) *
                exp(-pow(t - mu - delay, 2) / (2.0 * sigsqr)), 0, 0 );
    else
        return Eigen::Vector3d::Zero();

}

Eigen::Vector3cd pairwise_fld(Eigen::Vector3d &P_d0,
                              Eigen::Vector3d &P_d1,
                              Eigen::Vector3d &P_d2,
                              Eigen::Vector3d &dr){
    Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
    Eigen::Matrix3d imrr = Eigen::Matrix3d::Identity() - rr;
    Eigen::Matrix3d im3rr = Eigen::Matrix3d::Identity() - 3.0 * rr;

    double dist = dr.norm();

    return -pow(c0, 2) * prop_constant * hbar * exp( -iu * omega * dist / c0 ) * (
            im3rr * P_d0 / pow(dist, 3) +
            im3rr * ( P_d1 + iu * omega * P_d0 ) / ( c0 * pow(dist, 2) ) +
            imrr * ( P_d2 + 2.0 * iu * omega * P_d1 - pow(omega, 2) * P_d0 ) / ( pow(c0, 2) * dist ) );

}


Eigen::Vector3cd pairwise_fld_nf(Eigen::Vector3d &P_d0,
                                 Eigen::Vector3d &P_d1,
                                 Eigen::Vector3d &P_d2,
                                 Eigen::Vector3d &dr){
    Eigen::Matrix3d rr = dr * dr.transpose() / dr.squaredNorm();
    Eigen::Matrix3d iprr = Eigen::Matrix3d::Identity() + rr;
    Eigen::Matrix3d im3rr = Eigen::Matrix3d::Identity() - 3.0 * rr;

    double dist = dr.norm();
    // steady state assumption
    // return -pow(c0, 2) * prop_constant * hbar * (
     //   1.0 * ( im3rr / pow(dist, 3) - pow(omega, 2) * iprr / ( 2.0 * pow(c0, 2) * dist ) ) * P_d0 );

    // no steady state assumption
    return -pow(c0, 2) * prop_constant * hbar * (
            im3rr * P_d0 / pow(dist, 3) +
            iprr * ( P_d2 + 2.0 * iu * omega * P_d1 - pow(omega, 2) * P_d0 ) / ( 2.0 * pow(c0, 2) * dist ) );

}

Eigen::Vector3cd rr_fld(Eigen::Vector3d &P_d0,
                        Eigen::Vector3d &P_d1,
                        Eigen::Vector3d &P_d2){
    return 
        // hbar * beta /  pow( 5.2917721e-4, 2 ) 
        pow(c0, 2) * prop_constant * hbar * 2.0 / ( 3.0 * pow(c0, 3) ) 
        * ( iu * pow(omega, 3) * P_d0 + 3.0 * pow(omega, 2) * P_d1 - 3.0 * iu * omega * P_d2 );
}

int main(int argc, char *argv[]){

   // set up dots and history of dots
    auto dots = std::make_shared<DotVector>(import_dots("dots/dots0.cfg"));
    auto obss = std::make_shared<DotVector>(import_dots("dots/dotsobs.cfg"));
    const int ndots = (*dots).size();
    
    std::cout << "Running with " << ndots << " dots" << std::endl;
    // std::cout << "Sigma: " << sigma << std::endl;

    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(ndots, 22, steps); 
    history->fill(Eigen::Vector2cd::Zero());
    
    for (int n = 0; n < ndots; ++n)
        for (int i = -22; i < steps; ++i)
            // Again, assume dipoles along x
            history->array_[n][i][0] = Eigen::Vector2cd( 0, P_d0(i*dt, 0)[0] / (*dots)[n].dipole().norm() );

    // set up propagator
    Propagation::RotatingEFIE propagator(c0, prop_constant, omega, beta);

    // set up interactions
    /*const int expansion = 5;
    const int border = 1;
    const double ds = c0*dt; //vm["ds"].as<double>()*lambda;
    Eigen::Vector3d dsvec = Eigen::Vector3d::Ones() * ds; 
    // AIM::Grid grid(dsvec, expansion, *dots); // sort dots according to AIM scheme, even if not calculating via AIM!
 
    // auto direct_interaction = std::make_shared<DirectInteraction>(dots, dots, history, propagator, interp, c0, dt, omega, beta, hbar);

    const int transit_steps = grid.max_transit_steps(c0, dt) + interp;
    Propagation::EFIE<cmplx> propagator_aim(c0, prop_constant);

    auto aim_interaction = std::make_shared<AIM::Interaction>(
        dots, history, propagator_aim, dsvec, interp, expansion, border, c0, dt,
        AIM::Expansions::EFIE(transit_steps, c0, dt),
        AIM::Normalization::Laplace(prop_constant) );
    */

    // calculate field
    *obss = *dots;
    auto fld_interaction = std::make_shared<DirectInteraction>(dots, obss, history, propagator, interp, c0, dt, omega, beta, hbar);

    std::ofstream outfile1, outfile2;
    outfile1.open("outstatic/fld_dir.dat");
    outfile1 << std::scientific << std::setprecision(15);
    outfile2.open("outstatic/fld_anl.dat");
    outfile2 << std::scientific << std::setprecision(15);

    double normdiff = 0;
    double totalfld = 0;
    const bool rrflag = atoi(argv[1]);

    for (int i = 0; i < steps; ++i) {
        // direct interaction fld
        set_dipolevec(obss, Eigen::Vector3d(hbar,0,0));
        auto fldx = fld_interaction->evaluate(i);

        set_dipolevec(obss, Eigen::Vector3d(0,hbar,0));
        auto fldy = fld_interaction->evaluate(i);

        set_dipolevec(obss, Eigen::Vector3d(0,0,hbar));
        auto fldz = fld_interaction->evaluate(i);

        for (int nobs = 0; nobs < ndots; ++nobs){
            Eigen::Vector3cd fld_anl = Eigen::Vector3cd::Zero();

            // analytic fld
            for (int nsrc = 0; nsrc < ndots; ++nsrc){
                Eigen::Vector3d dr( separation( (*dots)[nsrc], (*dots)[nobs] ) );

                double delay = dr.norm() / c0;        

                Eigen::Vector3d obs_d0 = P_d0( i*dt, delay );
                Eigen::Vector3d obs_d1 = P_d1( i*dt, delay );
                Eigen::Vector3d obs_d2 = P_d2( i*dt, delay );

                if (nsrc != nobs) { // pairwise interactions
                    // cout << dr.norm() << endl;
                    fld_anl += pairwise_fld( obs_d0, obs_d1, obs_d2, dr );
                    //if (rrflag)
                    //    fld_anl += rr_fld( obs_d0, obs_d1, obs_d2 );
                } else
                    fld_anl += rr_fld( obs_d0, obs_d1, obs_d2 );

            }

            outfile1 << abs( fldx[nobs] ) << " " << abs( fldy[nobs] ) << " " << abs( fldz[nobs] ) << " ";
            // outfile1 << sqrt( pow( abs(fldx[n]), 2 ) + pow( abs(fldy[n]), 2 ) + pow( abs(fldz[n]), 2 ) ) << " ";
            outfile2 << abs( fld_anl[0] ) << " " << abs( fld_anl[1] ) << " " << abs( fld_anl[2] ) << " ";
            // outfile3 << abs( fldx[nobs] ) - abs( fld_anl[0] ) << " " << abs( fldy[nobs] ) - abs( fld_anl[1] ) << " " << abs( fldz[nobs] ) - abs( fld_anl[2] ) << " ";


            double fldAbs = sqrt( pow( abs(fldx[nobs]), 2 ) + pow( abs(fldy[nobs]), 2 ) + pow( abs(fldz[nobs]), 2 ) );
            double fldAbsAnl = sqrt( pow( abs(fld_anl[0]), 2 ) + pow( abs(fld_anl[1]), 2 ) + pow( abs(fld_anl[2]), 2 ) );
            totalfld += fldAbsAnl; // abs( fld_anl[0] );
            normdiff += pow( fldAbsAnl - fldAbs, 2 );
                //pow( abs( fldx[nobs] ) - abs( fld_anl[0] ), 2 );

        }

        outfile1 << "\n";
        outfile2 << "\n";
    }
    
    cout << "Rel normdiff (|E|): " << sqrt(normdiff) / totalfld << endl;

    outfile1.close();
    outfile2.close();

}
