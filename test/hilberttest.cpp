#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
//#include <omp.h>

#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
//#include "integrator/integrator_newton.h"
//#include "interactions/AIM/aim_interaction.h"
// #include "interactions/AIM/grid.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << setprecision(12) << scientific;

    // parameters
    const double omega = 2278.9013;

    const int num_src = atoi(argv[1]);
    const int num_obs = 0;
    const double dt = 1.0e-2; // 2.0*M_PI / ( omega / 1000 ) * 0.001; 
    const int tmax = 10;
    const double num_timesteps = tmax/dt;
    const int tmult = 1; // 50 * pow(10, atoi(argv[2]) );

    const int interpolation_order = 4;
    const bool solve_type = 0;  
    const bool interacting = atoi(argv[4]);
    const bool rotating = 0;

    // constants
    const double c0 = 299.792458, hbar = 0.65821193, mu0 = 2.0133545e-04, d0 = 5.2917721e-4 * 1.0;
    double beta = 0.0 / pow( omega, 3 );
	    // mu0 * pow( omega, 3 ) * pow( d0, 2 ) / ( 6.0 * M_PI * hbar * c0 );
        // ( atoi(argv[5]) ? 0.1 * pow(2, atoi(argv[5]) - 1) : 0.0 );

    const double k0 = omega/c0, lambda = 2.0*M_PI/k0;    

    // AIM
    const double ds = 0.050*lambda;
    Eigen::Vector3d grid_spacing(ds, ds, ds);
    const int expansion_order = 4;
    const int border = 1;

    cout << "Initializing..." << endl;
    std::cout << "  Timestep: " << dt << std::endl;
    std::cout << "  Num timesteps: " << num_timesteps << std::endl;
    std::cout << "  Num sources: " << num_src << std::endl;

    string idstr(argv[5]);
 
    auto qds = make_shared<DotVector>(import_dots("./dots/dots"+idstr+".cfg"));
    qds->resize(num_src);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

    auto obs = make_shared<DotVector>(import_dots("./dots/dotsobs.cfg"));
    obs->resize(num_obs);

    // std::cout << obs->size() << std::endl;

    // == HISTORY ====================================================

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_src, 22, num_timesteps);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0), num_src );
    // history->initialize_past( qd_path );

    // == INTERACTIONS ===============================================

    const double propagation_constant = mu0 / (4 * M_PI * hbar);
//    Propagation::SelfEFIE dyadic(c0, propagation_constant, beta);

/*    Propagation::EFIE<cmplx> dyadic;
    Propagation::EFIE<cmplx> dyadic_self;

    if (rotating) {
        dyadic = Propagation::RotatingEFIE(c0, propagation_constant, omega, beta);
        dyadic_self = Propagation::SelfRotatingEFIE(c0, propagation_constant, omega, beta);
    } else {
        dyadic = Propagation::EFIE<cmplx>(c0, propagation_constant, beta);
        dyadic_self = Propagation::SelfEFIE(c0, propagation_constant, beta);
    }
*/
    
/*    Propagation::SelfRotatingEFIE dyadic_self(c0, propagation_constant,
                                                omega, beta);
    Propagation::RotatingEFIE dyadic(c0, propagation_constant,
                                           omega, beta);
*/
    
      Propagation::SelfEFIE dyadic_self(c0, propagation_constant,
                                      beta);
      Propagation::EFIE<cmplx> dyadic(c0, propagation_constant,
                                            beta);

    auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse.cfg"));

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
        make_shared<PulseInteraction>(qds, nullptr, pulse1, interpolation_order, c0, dt, hbar, rotating),
        make_shared<DirectInteraction>(qds, nullptr, history, dyadic_self,
                                                    interpolation_order, c0, dt, omega, beta, hbar, rotating) };

    if (interacting) {
        std::shared_ptr<InteractionBase> pairwise;

        pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                                    interpolation_order, c0, dt, omega, beta, hbar, rotating);
        interactions.push_back( pairwise );
    }

    // == INTEGRATOR =================================================

    if (solve_type) {
      std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
          std::make_unique<Integrator::BlochRHS>(
              dt, history, std::move(interactions), std::move(rhs_funcs));

      Integrator::PredictorCorrector<Eigen::Vector2cd> solver_pc(
          dt, 18, 22, 3.15, history, std::move(bloch_rhs));

      cout << "Solving..." << endl;
      solver_pc.solve();

    } else {
      Integrator::HilbertIntegrator<Eigen::Vector2cd> solver_hilbert(
        dt, omega, omega, 1, history, std::move(interactions), pulse1);
      
      cout << "Solving..." << endl;
      solver_hilbert.solve();
      /*Integrator::NewtonJacobian<Eigen::Vector2cd> solver_nt(
        dt, beta, omega, interpolation_order, rotating, history, std::move(interactions));
      
      cout << "Solving..." << endl;
      solver_nt.solve();*/

    }
    // start_time = omp_get_wtime();

    // double elapsed_time = omp_get_wtime() - start_time;

    // cout << "Elapsed time: " << elapsed_time << "s" << std::endl;

    // == FIELD INTERACTIONS ===============================================

/*    *obs = *qds; // examine field at sources only
 
    std::vector<std::shared_ptr<InteractionBase>> interactions_fld{ 
        make_shared<PulseInteraction>(qds, obs, pulse1, interpolation_order, c0, dt, hbar, rotating),
        make_shared<DirectInteraction>(qds, obs, history, dyadic_self,
                                                    interpolation_order, c0, dt, omega, beta, hbar, rotating) };

    if(interacting) {
        std::shared_ptr<InteractionBase> pairwise_fld;

        pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                      interpolation_order, c0, dt, omega, beta, hbar, rotating);
        interactions.push_back( pairwise_fld );
    }
*/

    // == OUTPUT =====================================================

    cout << "Writing output..." << endl;

    string dotstr(argv[1]);
    string prfx = "./out/";
    string sffx = dotstr + "dots_" + idstr + ".dat";

    string rhostr = prfx + "rho_" + sffx; 
/*    string fldstr = prfx + "fld_" + sffx;
    string fldxstr = prfx + "fldx_" + sffx;
    string fldystr = prfx + "fldy_" + sffx;
    string fldzstr = prfx + "fldz_" + sffx;
*/
    ofstream rhofile(rhostr);
    rhofile << scientific << setprecision(15);

/*    ofstream fldfile(fldstr);
    fldfile << scientific << setprecision(15);
    ofstream fldxfile(fldxstr);
    fldxfile << scientific << setprecision(15);
    ofstream fldyfile(fldystr);
    fldyfile << scientific << setprecision(15);
    ofstream fldzfile(fldzstr);
    fldzfile << scientific << setprecision(15);
*/
//    const int NUM_DERIV = 4;

    for(int t = 0; t < num_timesteps; ++t) {

      if (t%tmult == 0){

        // print rho
        for(int n = 0; n < num_src; ++n)
          rhofile << history->array_[n][t][0][0].real() << " "
                  << history->array_[n][t][0][1].real() << " "
                  << history->array_[n][t][0][1].imag() << " ";
  //                << history->array_[n][t][0][2].real() << " ";
  //                << history->array_[n][t][0][2].imag() << " ";

        // print pol and derivatives
/*        for(int deriv = 0; deriv < NUM_DERIV; ++deriv){
          auto eval_and_sum = 
            [t,deriv](const InteractionBase::ResultArray &r,
            const std::shared_ptr<InteractionBase> &interaction){
            return r + interaction->evaluatefld(t,deriv);
          };
          auto nilfld = InteractionBase::ResultArray::Zero(obs->size(), 1).eval();

          set_dipolevec(obs, Eigen::Vector3d(hbar,0,0));
          auto fldx = 
            interactions_fld[interactions_fld.size()-1]->evaluatefld(t,deriv); 
            //std::accumulate(
            //  interactions_fld.begin(), interactions_fld.end(), nilfld, eval_and_sum);  
          
          set_dipolevec(obs, Eigen::Vector3d(0,hbar,0));
          auto fldy = 
            interactions_fld[interactions_fld.size()-1]->evaluatefld(t,deriv); 
            //std::accumulate(
            //  interactions_fld.begin(), interactions_fld.end(), nilfld, eval_and_sum);  

          set_dipolevec(obs, Eigen::Vector3d(0,0,hbar));
          auto fldz = 
            interactions_fld[interactions_fld.size()-1]->evaluatefld(t,deriv); 
            //std::accumulate(
            //  interactions_fld.begin(), interactions_fld.end(), nilfld, eval_and_sum);  

          for(int n = 0; n < obs->size(); ++n){
            //fldfile << fldx[n] << "," << fldy[n] << "," << fldz[n] << ",";
            fldxfile << real(fldx[n]) << " ";
            fldyfile << real(fldy[n]) << " ";
            fldzfile << real(fldz[n]) << " ";

            fldfile << sqrt( pow( abs(fldx[n]), 2 ) + pow( abs(fldy[n]), 2 ) + pow( abs(fldz[n]), 2 ) ) << " ";
          }
        }
*/      
      rhofile << "\n";
/*
      fldfile << "\n";
      fldxfile << "\n";
      fldyfile << "\n";
      fldzfile << "\n";
*/
      }
    }

  return 0;
}
