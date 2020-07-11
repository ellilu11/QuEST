#include "configuration.h"
Configuration config; //One configuration to rule them all...

using namespace std;
namespace po = boost::program_options;

po::variables_map parse_configs(int argc, char *argv[]) {
  string config_path;

  po::options_description cmd_line_description("Command line options");
  cmd_line_description.add_options()
    ("help", "print this help message")
    ("version,v", "print version string")
    ("config,c", po::value<string>(&config_path)->default_value("input.cfg"), "path to configuration file")
    ("fast,f",   po::bool_switch()->default_value(false), "employ fast methods to calculate potentials")
    ("time,t",   po::bool_switch(&config.report_time_data)->default_value(false), "report execution time data");
  ;

  po::options_description files_description("Configuration files");
  files_description.add_options()
    ("files.qd_path",    po::value<std::string>(&config.qd_path)->default_value("dots.cfg"), "path to dot configuration file")
    ("files.pulse_path", po::value<std::string>(&config.pulse_path)->default_value("pulse.cfg"), "path to pulse configuration file")
  ;

  po::options_description constants_description("Physical constants");
  constants_description.add_options()
    ("constants.c0",   po::value<double>(&config.c0)->required(), "speed of light in vacuum")
    ("constants.hbar", po::value<double>(&config.hbar)->required(), "reduced Planck constant")
    ("constants.mu0",  po::value<double>(&config.mu0)->required(), "vacuum permeability")
    ("constants.laser frequency",  po::value<double>(&config.omega)->required(), "(angular) frequency of incident pulse")
    ("constants.RR damping",  po::value<double>(&config.beta)->required(), "radiation reaction damping factor")    
  ;   

  po::options_description parameters_description("System parameters");
  parameters_description.add_options()
    ("parameters.num_particles", po::value<int>(&config.num_particles)->required(), "number of particles in the system")
    ("parameters.total_time",
     po::value<double>(&config.total_time)
         ->required()
         ->notifier([](const double total_time) {
           config.num_timesteps =
             static_cast<int>(std::ceil(total_time / config.dt));
         }),"total simulation duration")
    ("parameters.timestep",            po::value<double>(&config.dt)->required(), "timestep size")
    ("parameters.interpolation_order", po::value<int>(&config.interpolation_order)->required(), "order of the temporal Lagrange interpolants")
    ("parameters.fast",                po::bool_switch()->default_value(false), "in-file alias of --fast")
    ("parameters.interacting",         po::value<bool>(&config.interacting)->default_value(true), "interaction flag")
    ("parameters.rotating",            po::value<bool>(&config.rotating)->default_value(true), "rotating frame flag")
  ;

  po::options_description aim_description("AIM & Grid parameters");
  aim_description.add_options()
    ("AIM.grid spacing",    po::value<std::string>()->required()->multitoken(), "spacing between grid nodes (distance units)")
    ("AIM.expansion order", po::value<int>(&config.expansion_order)->required(), "spatial grid expansion order")
    ("AIM.border",          po::value<int>(&config.border)->required(), "infinity-norm distance within which to consider nearfield pairs")
  ;

  po::options_description cmdline_options, file_options;
  cmdline_options.add(cmd_line_description);
  file_options.add(files_description)
      .add(constants_description)
      .add(parameters_description)
      .add(aim_description);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    po::options_description visible("QuEST options");
    visible.add(cmd_line_description).add(file_options);
    cout << visible << "\n";

    throw CommandLineException();
  }

  if (vm.count("version")) {
    cout << "Quantum Electromagnetics Simulation Toolbox (QuEST), version 0.2"
         << endl;
    cout << "Compiled with " << __VERSION__ << " on " << __DATE__ << " ("
         << __TIME__ << ")" << endl;
    //cout << "Git revision: " << __GIT_HASH__ << " (branch: " << __GIT_BRANCH__
    //     << ")" << endl;
    throw CommandLineException();
  }

  ifstream ifs(config_path.c_str());
  if (!ifs) {
    cerr << "ERROR: " << config_path << " not found" << endl;
    throw CommandLineException();
  } else {
    po::store(po::parse_config_file(ifs, file_options), vm);
    po::notify(vm);

    std::istringstream iss(vm["AIM.grid spacing"].as<std::string>());
    std::vector<double> tokens{
      std::istream_iterator<double>(iss),
      std::istream_iterator<double>()
    };

    config.grid_spacing = Eigen::Vector3d(tokens.at(0), tokens.at(1), tokens.at(2));
    config.sim_type = static_cast<Configuration::SIMULATION_TYPE>(
        vm["fast"].as<bool>() || vm["parameters.fast"].as<bool>());
  }

  return vm;
}
