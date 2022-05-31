#include "worm.hpp"
#include <iostream>


int main(int argc, char** argv)
{
  std::cout << "# " << worm::code_name() << std::endl;

  // Creates the parameters for the simulation
  // If an hdf5 file is supplied, reads the parameters there
  std::cout << "# Initializing parameters..." << std::endl;
  alps::parameters_type<worm>::type parameters(argc, argv, "/parameters");
  worm::define_parameters(parameters);
  if (parameters.help_requested(std::cout)) {
    exit(0);
  }
  // XXX not possible to select help for the user specified model and lattice?

    
  std::cout << "# Constructing worm..." << std::endl;
  worm sim(parameters);
  std::cout << "# ...done" << std::endl;
    
  //return(0);

  // If needed, restore the last checkpoint
  std::string checkpoint_file = parameters["checkpoint"];
  if (parameters.is_restored()) {
    std::cout << "# Restoring checkpoint from " << checkpoint_file << std::endl;
    sim.load(checkpoint_file);
  }
  else {
    std::cout << "# Creating simulation..." << std::endl;
    sim.initialize();
    std::cout << "# ...done" << std::endl;
  }
  sim.print_params(std::cout);
  std::ofstream outfile("qc_params");
  sim.print_params(outfile);
  outfile.close();

  // return 0;

  // Run the simulation
  std::cout << "# Running simulation..." << std::endl;
  //sim.run(alps::stop_callback(size_t(parameters["timelimit"])));
  sim.run(alps::stop_callback(size_t(parameters["runtimelimit"])));

  // Checkpoint the simulation
  std::cout << "# Checkpointing simulation after testing..." << std::endl;
  try {
    sim.test_conf();
  }
  catch (const char* e) {
    cerr << e << endl;
    exit(1);
  }
  sim.save(checkpoint_file);

    
  
  alps::results_type<worm>::type results = alps::collect_results(sim);
  std::ofstream fstat("qc_statistics");
  sim.print_update_statistics(fstat);
  fstat.close();
  
  // Print results
  std::cout << "# Result:" << std::endl;
  std::cout << results << std::endl;

  // Saving to the output file
  std::string output_file = parameters["outputfile"];
  alps::hdf5::archive ar(output_file, "w");
  ar["/parameters"] << parameters;
  ar["/simulation/results"] << results;
  
  sim.output_final(results,ar);

  return 0;
}
