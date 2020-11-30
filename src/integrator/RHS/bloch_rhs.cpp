#include "bloch_rhs.h"

Integrator::BlochRHS::BlochRHS(
    const double hbar, 
    const double mu0,
    const double c0,
    const double omega,
    const double dt, 
    const int num_timesteps,
    const std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history,
    std::vector<std::shared_ptr<InteractionBase>> interactions,
    std::vector<std::shared_ptr<InteractionBase>> efld_interactions,
    std::vector<std::shared_ptr<InteractionBase>> bfld_interactions,
    std::vector<BlochFunctionType> rhs_functions,
    std::shared_ptr<DotVector> obss,
    const int task_idx,
    const bool getflux)
    : Integrator::RHS<Eigen::Vector2cd>(dt, history),
      hbar(hbar), mu0(mu0), eps0(1.0/(mu0*c0*c0)), omega(omega),
      num_timesteps(num_timesteps),
      num_solutions(history->num_particles),
      interactions(std::move(interactions)),
      efld_interactions(std::move(efld_interactions)),
      bfld_interactions(std::move(bfld_interactions)),
      rhs_functions(std::move(rhs_functions)),
      obss(std::move(obss)),
      num_obs(obss->size()),
      getflux(getflux)
{
  outfile.open("./out/fld" + std::to_string(task_idx) + ".dat");
  outfile << std::scientific << std::setprecision(15);

  // fluxfile.open("./out/flux" + std::to_string(task_idx) + ".dat");
  // fluxfile << std::scientific << std::setprecision(15);
}

Integrator::BlochRHS::~BlochRHS()
{
  outfile.close();
  fluxfile.close();
}

void Integrator::BlochRHS::evaluate(const int step) const
{
  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->set_value(solution, step, 1) = rhs_functions[solution](
        history->get_value(solution, step, 0),
	      projected_rabi[solution],
        step);
  }
}

void Integrator::BlochRHS::evaluate_present(const int step) const
{
  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_present(step);
      };
  auto nil = InteractionBase::ResultArray::Zero(num_solutions, 1).eval();

  auto projected_rabi = std::accumulate(
      interactions.begin(), interactions.end(), nil, eval_and_sum);

  for(int solution = 0; solution < num_solutions; ++solution) {
    history->set_value(solution, step, 1) = rhs_functions[solution](
        history->get_value(solution, step, 0),
	      projected_rabi[solution],
        step);
  }
}

void Integrator::BlochRHS::evaluate_field(const int step)
{
  auto nil = InteractionBase::ResultArray::Zero(num_obs, 1).eval();

  auto eval_and_sum =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_field(step);
      };
  auto eval_and_sum_bfld =
      [step](const InteractionBase::ResultArray &r,
             const std::shared_ptr<InteractionBase> &interaction) {
        return r + interaction->evaluate_field(step, 1);
      };

  int idx = 0;
  
  set_dipole_of_dots( obss, Eigen::Vector3d(hbar, 0, 0) );
  auto efldx = 
      // efld_interactions[idx]->evaluate_field(step); 
      std::accumulate(
        efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldx = 
      // bfld_interactions[idx]->evaluate_field(step, 1);
      std::accumulate(
        bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
 
  set_dipole_of_dots( obss, Eigen::Vector3d(0, hbar, 0) );
  auto efldy = 
      // efld_interactions[idx]->evaluate_field(step); 
      std::accumulate(
        efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldy = 
      // bfld_interactions[idx]->evaluate_field(step, 1);
      std::accumulate(
        bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);
  
  set_dipole_of_dots( obss, Eigen::Vector3d(0, 0, hbar) );
  auto efldz = 
      // efld_interactions[idx]->evaluate_field(step); 
      std::accumulate(
        efld_interactions.begin(), efld_interactions.end(), nil, eval_and_sum);
  auto bfldz = 
      // bfld_interactions[idx]->evaluate_field(step, 1);
      std::accumulate(
        bfld_interactions.begin(), bfld_interactions.end(), nil, eval_and_sum_bfld);

  double flux = 0.0;
  double energy = 0.0;
  double time = step*dt;

  for(int obs = 0; obs < num_obs; ++obs) {
    Eigen::Vector3cd efld(efldx[obs], efldy[obs], efldz[obs]);
    Eigen::Vector3cd bfld(bfldx[obs], bfldy[obs], bfldz[obs]);

    Eigen::Vector3d poynting = (efld.real()).cross(bfld.real()) / ( mu0 ); // Fix frame
    //Eigen::Vector3d poynting = ( efld.cross(bfld.conjugate()) + 
    //                             efld.cross(bfld) * std::exp(2.0*iu*omega*time) ).real() / ( 2.0*mu0 ); // Rot frame
 
    Eigen::Vector3d pos = (*obss)[obs].position();
    double r0 = pos.norm();
    Eigen::Vector3d rhat = pos / pos.norm();
    double measure = (*obss)[obs].frequency(); 

    //if ( step == 4800000 )
      /*outfile << std::abs(efldx[obs]) << " " 
              << std::abs(efldy[obs]) << " " 
              << std::abs(efldz[obs]) << " "    
              << std::abs(bfldx[obs]) << " " 
              << std::abs(bfldy[obs]) << " " 
              << std::abs(bfldz[obs]) << " "
              << poynting[0] << " " 
              << poynting[1] << " " 
              << poynting[2] << " "
              << rhat[0] << " "
              << rhat[1] << " "
              << rhat[2] << " ";*/

    if (getflux) {
      // flux_real = flux_real + measure * std::real( poynting.dot(rhat) );
      // flux_imag = flux_imag + measure * std::imag( poynting.dot(rhat) );
      double flux_per_area = poynting.dot(rhat) ;
      flux = flux + measure * flux_per_area;
      // if ( step == 4800000 )
      //  outfile << flux_per_area << std::endl;
    } else {
      double energy_density = ( efld.squaredNorm() * eps0 + bfld.squaredNorm() / mu0 ) / 2.0;
      energy = energy + measure * energy_density;
    }
  }
    
  if (getflux)
    std::cout << flux << std::endl;
  else 
    std::cout << energy << std::endl;

  outfile << "\n";
  if ( step > num_timesteps*0.90){
    // outfile.flush();
    // fluxfile.flush();
  }
}

