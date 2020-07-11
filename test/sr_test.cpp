#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>
// #include <omp.h>

#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
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
    const int num_src = atoi(argv[1]);
    const int num_obs = atoi(argv[2]);
    const double tmax = 10;
    const double dt = 1e-3; // rotframe: sigma = 1.0ps -> dt <= 0.52e-1
                              // fixframe: omega = 2278.9013 mev/hbar -> dt <= 1.379e-4
    const int num_timesteps = tmax/dt;

    const int interpolation_order = 5;
    const bool sim_type = atoi(argv[3]);   // change this!
    const bool interacting = atoi(argv[4]);
    const bool rotating = 1;

    // constants
    const double c0 = 299.792458, hbar = 0.65821193, mu0 = 2.0133545e-04;
    const double omega = 2278.9013, d0 = 5.2917721e-4 * 1.0;
    double beta = 1.0;
	    // mu0 * pow( omega, 3 ) * pow( d0, 2 ) / ( 6.0 * M_PI * hbar * c0 );
        // ( atoi(argv[5]) ? 0.1 * pow(2, atoi(argv[5]) - 1) : 0.0 );

    const double k0 = omega/c0, lambda = 2.0*M_PI/k0;    

    // AIM
    const double ds = 0.050*lambda;
    Eigen::Vector3d grid_spacing(ds, ds, ds);
    const int expansion_order = 4;
    const int border = 1;

    cout << "Initializing..." << endl;
    std::cout << "  Running in "
              << (sim_type
                      ? "FAST"
                      : "SLOW")
              << " mode" << std::endl;
    std::cout << "  Interacting particles: " 
              << (interacting
                      ? "TRUE"
                      : "FALSE")
              << std::endl;
    std::cout << "  Frame: "
              << ((rotating) ? "Rotating" : "Fixed") << std::endl;
    std::cout << "  Num timesteps: " << num_timesteps << std::endl;
    std::cout << "  Num sources: " << num_src << std::endl;
    std::cout << "  Num observers: " << num_obs << std::endl;
    std::cout << "  Beta: " << beta << std::endl;
    std::cout << "  ds/lambda: " << ds/lambda << std::endl;
    std::cout << "  AIM expansion order: " << expansion_order << std::endl;


    auto qds = make_shared<DotVector>(import_dots("dots.cfg"));
    qds->resize(num_src);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

    auto obs = make_shared<DotVector>(import_dots("dotsobs.cfg"));
    obs->resize(num_obs);

    // std::cout << obs->size() << std::endl;

    // == HISTORY ====================================================

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_src, 22, num_timesteps);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0), num_src );
    // history->initialize_past( qd_path );

    auto history_efld = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_src+num_obs, 22, num_timesteps);
    history_efld->fill(Eigen::Vector2cd::Zero());

    // == INTERACTIONS ===============================================

    const double propagation_constant = mu0 / (4 * M_PI * hbar);

    Propagation::RotatingEFIE dyadic(c0, propagation_constant,
                                             omega, beta);

    // Propagation::EFIE<cmplx> dyadic(c0, propagation_constant);

    /*dyadic = rotating ? 
        Propagation::RotatingEFIE(c0, propagation_constant,
                                              omega) :
        Propagation::EFIE<cmplx>(c0, propagation_constant);
*/
    std::shared_ptr<InteractionBase> pairwise;

    // AIM::Grid grid(grid_spacing, expansion_order, *qds);

    // double start_time = omp_get_wtime();

    if (sim_type) {
       /*const int transit_steps = grid.max_transit_steps(c0, dt) +
                                interpolation_order;


         AIM::Expansions::EFIE expansion_func(transit_steps, c0, dt);
  
       if (rotating)
            expansion_func = AIM::Expansions::RotatingEFIE(transit_steps, c0, dt, omega);
          else
            expansion_func = AIM::Expansions::EFIE(transit_steps, c0, dt);
     
        pairwise = make_shared<AIM::Interaction>(
              qds, history, dyadic, grid_spacing,
              interpolation_order, expansion_order, border,
              c0, dt,
              AIM::Expansions::RotatingEFIE(transit_steps, c0, dt, omega),
              // AIM::Expansions::EFIE(transit_steps, c0, dt),
              // AIM::Normalization::Laplace( propagation_constant ),
              AIM::Normalization::Helmholtz(omega / c0,
                                            propagation_constant),
              omega);*/

    } else {
        pairwise = make_shared<DirectInteraction>(qds, nullptr, history, dyadic,
                                                interpolation_order,
                                                c0, dt, omega, beta, hbar);
    }

    // cout << "Interaction setup elapsed time: " << omp_get_wtime() - start_time << "s" << std::endl;

    auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse.cfg"));

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
        make_shared<PulseInteraction>(qds, nullptr, pulse1, hbar, dt, rotating)};

    if(interacting) interactions.push_back( pairwise );

    // == INTEGRATOR =================================================

    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            dt, history, history_efld, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver(
        dt, 18, 22, 3.15, history, std::move(bloch_rhs));

    // start_time = omp_get_wtime();

    cout << "Solving..." << endl;
    solver.solve();

    // double elapsed_time = omp_get_wtime() - start_time;

    // cout << "Elapsed time: " << elapsed_time << "s" << std::endl;

    // == FIELD INTERACTIONS ===============================================

    *obs = *qds; // examine field at sources only
   
    std::shared_ptr<InteractionBase> pairwise_fld;
    pairwise_fld = make_shared<DirectInteraction>(qds, obs, history, dyadic,
                                                interpolation_order,
                                                c0, dt, omega, beta, hbar);
   
    // std::cout << pairwise_fld->obss->size() << std::endl;
 
    std::vector<std::shared_ptr<InteractionBase>> interactions_fld{ 
        make_shared<PulseInteraction>(obs, obs, pulse1, hbar, dt, rotating)};

    if(interacting) interactions_fld.push_back( pairwise_fld );

    // std::cout << (nullptr == 0) << std::endl;

    // == OUTPUT =====================================================

    cout << "Writing output..." << endl;

    // ofstream fluxfile("flux.dat", std::ios_base::app);
    // fluxfile << scientific << setprecision(15);

    ofstream rhofile("outsr/rho.dat");
    rhofile << scientific << setprecision(15);

    ofstream fldfile("outsr/fld.dat");
    fldfile << scientific << setprecision(15);

    // ofstream posfile("posdip.dat");
    // posfile << scientific << setprecision(15);

    for(int t = 0; t < num_timesteps; ++t) {

      if (t%1 == 0){

      auto eval_and_sum = 
        [t](const InteractionBase::ResultArray &r,
        const std::shared_ptr<InteractionBase> &interaction){
        return r + interaction->evaluate(t);
      };
      auto nilfld = InteractionBase::ResultArray::Zero(obs->size(), 1).eval();

      for(int n = 0; n < num_src; ++n)
        rhofile << history->array_[n][t][0][0].real() << " "
                << history->array_[n][t][0][1].real() << " "
                << history->array_[n][t][0][1].imag() << " ";

      set_dipolevec(obs, Eigen::Vector3d(hbar,0,0));
      auto fldx = interactions_fld[1]->evaluate(t); 
        // std::accumulate(
        //    interactions_fld.begin(), interactions_fld.end(), nilfld, eval_and_sum);  
      
      set_dipolevec(obs, Eigen::Vector3d(0,hbar,0));
      auto fldy = interactions_fld[1]->evaluate(t); 
        //std::accumulate(
        //    interactions_fld.begin(), interactions_fld.end(), nilfld, eval_and_sum);  

      set_dipolevec(obs, Eigen::Vector3d(0,0,hbar));
      auto fldz = interactions_fld[1]->evaluate(t); 
        //std::accumulate(
        //    interactions_fld.begin(), interactions_fld.end(), nilfld, eval_and_sum);  

      // double flux = 0;

      for(int n = 0; n < obs->size(); ++n){
        // fldfile << fldx[n] << "," << fldy[n] << "," << fldz[n] << ",";
        fldfile << sqrt( pow( abs(fldx[n]), 2 ) + pow( abs(fldy[n]), 2 ) + pow( abs(fldz[n]), 2 ) ) << " ";

        // double rad = (*obs)[n]->position().norm();
        // flux += ( norm(fldx[n]) + norm(fldy[n]) + norm(fldz[n]) ); // * (*obs)[n]->position[2] / rad;
      }
      fldfile << "\n";
 
      //  flux += std::norm( history_efld->array_[n][t][0][0] );
   
      rhofile << "\n";
      // fluxfile << flux / obs->size() << " ";
      }
    }
    // fluxfile << "\n";

    /*for(int n = 0; n < (num_src + num_obs); ++n){
        posfile << (*qds)[n].pos[0] << " " << (*qds)[n].pos[1] << " " << (*qds)[n].pos[2] << " ";
        posfile << (*qds)[n].dipr[0] << " " << (*qds)[n].dipr[1] << " " << (*qds)[n].dipr[2] << std::endl;
    }*/
 
  return 0;
}
