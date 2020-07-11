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
    const double tmax = 1;
    const double dt = 0.2e-2; // rotframe: sigma = 0.1ps -> dt <= 0.52e-2
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
    const double ds = 0.20*lambda;
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
    std::cout << "  Num particles: " << num_particles << std::endl;
    std::cout << "  ds/lambda: " << ds/lambda << std::endl;
    std::cout << "  AIM expansion order: " << expansion_order << std::endl;


    auto qds = make_shared<DotVector>(import_dots("dots.cfg"));
    qds->resize(num_particles);
    auto rhs_funcs = rhs_functions(*qds, omega, beta, rotating);

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
   
    ofstream posfile("dots100k_postsort.cfg");
    posfile << scientific << setprecision(15);

    for(int n = 0; n < num_particles; ++n){
        posfile << (*qds)[n].pos[0] << " " << (*qds)[n].pos[1] << " " << (*qds)[n].pos[2] << " ";
        posfile << 2278.9013 << " " << 10 << " " << 20 << " ";
        posfile << (*qds)[n].dipr[0] << " " << (*qds)[n].dipr[1] << " " << (*qds)[n].dipr[2] << std::endl;
    }
 
    /*ofstream outfilebin("output.bin");
    for(int t = 0; t < num_timesteps; ++t) {
      for(int n = 0; n < num_particles; ++n) {
       outfilebin.write((char*)&history->array_[n][t][0].transpose(), 2*sizeof(double));
      }
    }*/

  return 0;
}
