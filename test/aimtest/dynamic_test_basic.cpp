#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>
#include <omp.h>

#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << setprecision(12) << scientific;

    // parameters
    const int num_particles = atoi(argv[1]);
    const double tmax = 10;
    const double dt = 0.1e-1; // rotframe: sigma = 0.1ps -> dt <= 0.52e-2
                              // fixframe: omega = 2278.9013 mev/hbar -> dt <= 1.379e-4
    const int num_timesteps = tmax/dt;

    const int interpolation_order = 5;
    const bool sim_type = atoi(argv[2]);   // change this!
    const bool interacting = atoi(argv[3]);
    const bool rotating = 1;

    // constants
    const double c0 = 299.792458, hbar = 0.65821193, mu0 = 2.0133545e-04;
    const double omega = 2278.9013, beta = 0.00000;
    const double k0 = omega/c0, lambda = 2.0*M_PI/k0;    

    // AIM
    const double ds = 0.005*lambda;
    Eigen::Vector3d grid_spacing(ds, ds, ds);
    const int expansion_order = 5;
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
    std::cout << "  Num particles: " << num_particles << std::endl;
    std::cout << "  ds/lambda: " << ds/lambda << std::endl;
    std::cout << "  AIM expansion order: " << expansion_order << std::endl;

    auto qds = make_shared<DotVector>(import_dots("dots.cfg"));
    qds->resize(num_particles);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

    /*const double xlen = 16*ds, ylen = 16*ds, zlen = 0.020*num_particles*ds;
    Eigen::Vector3d rlen(xlen, ylen, zlen);
    const int nx = 3, ny = 3, nz = (0.001 * num_particles)-1;
    double dx = xlen/nx, dy = ylen/ny, dz = zlen/nz;

    std::vector<bool> NearestDot(num_particles);

    for (int i = 0; i < nx+1; ++i){
        for (int j = 0; j < ny+1; ++j){
            for (int k = 0; k < nz+1; ++k){

            Eigen::Vector3d r_grid(-xlen/2 + i*dx, -ylen/2 + j*dy, -zlen/2 + k*dz);

            double dist_prev = rlen.norm();
            int neardot_idx = 0;
            for (int n = 0; n < num_particles; ++n){
                Eigen::Vector3d dr = ((*qds)[n].pos - r_grid);
                double dist = dr.norm();
                if ( dist < dist_prev ){
                    dist_prev = dist;
                    neardot_idx = n;
                }
            }
            // assert( !NearestDot[neardot_idx] );
            NearestDot[neardot_idx] = 1;
            }
        }
    }

    std::cout << "Done finding nearest dots" << std::endl;
*/
    // == HISTORY ====================================================

    auto history = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_particles, 22, num_timesteps);
    history->fill(Eigen::Vector2cd::Zero());
    history->initialize_past( Eigen::Vector2cd(1,0), Eigen::Vector2cd(1,0) );
    // history->initialize_past( qd_path );

    auto history_efld = make_shared<Integrator::History<Eigen::Vector2cd>>(
        num_particles, 22, num_timesteps);
    history_efld->fill(Eigen::Vector2cd::Zero());

    // == INTERACTIONS ===============================================

    const double propagation_constant = mu0 / (4 * M_PI * hbar);

    Propagation::RotatingEFIE dyadic(c0, propagation_constant,
                                                omega);
    /*dyadic = rotating ? 
        Propagation::RotatingEFIE(c0, propagation_constant,
                                              omega) :
        Propagation::EFIE<cmplx>(c0, propagation_constant);
*/
    std::shared_ptr<InteractionBase> pairwise;

    AIM::Grid grid(grid_spacing, expansion_order, *qds);

    double start_time = omp_get_wtime();

    if (sim_type) {
       const int transit_steps = grid.max_transit_steps(c0, dt) +
                                interpolation_order;

        /*AIM::Expansions::RotatingEFIE expansion_func(transit_steps, c0, dt, omega);

         AIM::Expansions::EFIE expansion_func(transit_steps, c0, dt);
  
       if (rotating)
            expansion_func = AIM::Expansions::RotatingEFIE(transit_steps, c0, dt, omega);
          else
            expansion_func = AIM::Expansions::EFIE(transit_steps, c0, dt);
     */
        pairwise = make_shared<AIM::Interaction>(
              qds, history, dyadic, grid_spacing,
              interpolation_order, expansion_order, border,
              c0, dt,
              AIM::Expansions::RotatingEFIE(transit_steps, c0, dt, omega),
              AIM::Normalization::Helmholtz(omega / c0,
                                            propagation_constant),
              omega);

    } else {
        pairwise = make_shared<DirectInteraction>(qds, history, dyadic,
                                                interpolation_order,
                                                c0, dt);
    }

    cout << "Interaction setup elapsed time: " << omp_get_wtime() - start_time << "s" << std::endl;

    auto pulse1 = make_shared<Pulse>(read_pulse_config("pulse.cfg"));

    std::vector<std::shared_ptr<InteractionBase>> interactions{ 
        make_shared<PulseInteraction>(qds, pulse1, hbar, dt, rotating)};

    if(interacting) interactions.push_back( pairwise );

    // output files

    ofstream rhofile("rho.dat");
    rhofile << scientific << setprecision(15);

    ofstream fldfile("rabi.dat");
    fldfile << scientific << setprecision(15);

    ofstream posfile("posdip.dat");
    posfile << scientific << setprecision(15);

    /*ofstream rhofile2("outtest/rho_smpl.dat");
    rhofile << scientific << setprecision(15);

    ofstream fldfile2("outtest/rabi_smpl.dat");
    fldfile << scientific << setprecision(15);

    ofstream posfile2("outtest/posdip_smpl.dat");
    posfile << scientific << setprecision(15);*/

    ofstream iprfile("ipr.dat");
    posfile << scientific << setprecision(15);

    // write data only at dots closest to grid points
    for(int n = 0; n < num_particles; ++n){
        posfile << (*qds)[n].pos[0] << " " << (*qds)[n].pos[1] << " " << (*qds)[n].pos[2] << " ";
        posfile << (*qds)[n].dipr[0] << " " << (*qds)[n].dipr[1] << " " << (*qds)[n].dipr[2] << std::endl;

       /* if (NearestDot[n]){
            posfile2 << (*qds)[n].pos[0] << " " << (*qds)[n].pos[1] << " " << (*qds)[n].pos[2] << " ";
            posfile2 << (*qds)[n].dipr[0] << " " << (*qds)[n].dipr[1] << " " << (*qds)[n].dipr[2] << std::endl;

        }*/
    }
 
    // == INTEGRATOR & OUTPUT ========================================

    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            dt, history, history_efld, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver(
        dt, 18, 22, 3.15, history, std::move(bloch_rhs));

    start_time = omp_get_wtime();

    cout << "Solving..." << endl;

    std::vector<double> ipr(num_timesteps);
    std::vector<double> ipr2(num_timesteps);

     for(int t = 0; t < num_timesteps; ++t) {

      solver.solve_step(t);

      double denom = 0;

      for(int n = 0; n < num_particles; ++n) {
        rhofile << history->array_[n][t][0].transpose() << " ";
        fldfile << history_efld->array_[n][t][0][0] << " " ;;
       /* if (NearestDot[n]){
           rhofile2 << history->array_[n][t][0].transpose() << " ";
            fldfile2 << history_efld->array_[n][t][0][0] << " " ;;
        }*/
        double p = abs( history->array_[n][t][0][1] );
        ipr[t] += pow( p, 4 );
        denom += pow( p, 2 );

      }
      ipr2[t] = ipr[t] / pow( denom, 2 );
      iprfile << ipr2[t] << " " << ipr[t] << " " << pow( denom, 2 ) << "\n";      

      rhofile << "\n";
      fldfile << "\n";

      /*rhofile2 << "\n";
      fldfile2 << "\n";*/
    }

    double elapsed_time = omp_get_wtime() - start_time;

    cout << "Elapsed time: " << elapsed_time << "s" << std::endl;

  return 0;
}
