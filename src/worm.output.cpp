#include "worm.hpp"

void worm::print_params(std::ostream& os) const {
  os << "# run time limit                                    : " << runtimelimit << "\n";
  os << "# maximum number of sweeps                          : " << total_sweeps << "\n";
  os << "# maximum number of sweeps for thermalization       : " << thermalization_sweeps << "\n";
  //os << "# Nloop                                            : " << Nloop << "\n";
  os << "# inverse temperature                               : " << beta << "\n";
  //os << "# maximum occupation number                         : " << nmax << "\n";
  os << "# canonical measurement value                       : " << canonical << "\n";
  os << "# random number seed                                : " << seed << "\n";
  os << "# Ntest                                             : " << Ntest<< "\n";
  os << "# Nsave                                             : " << Nsave<< "\n";
  os << "# Nmeasure                                          : " << Nmeasure << "\n";
  os << "# energy shift                                      : " << E_off << "\n";
  os << "# weight worm configurations                        : " << C_worm << "\n";
  os << "# update probability insert worm                    : " << update_prob[update_tag::insertworm]    << "\t" << update_prob_cuml[update_tag::insertworm]    << "\n";
  os << "# update probability move worm                      : " << update_prob[update_tag::moveworm]      << "\t" << update_prob_cuml[update_tag::moveworm]    << "\n";
  os << "# update probability insert kink                    : " << update_prob[update_tag::insertkink]    << "\t" << update_prob_cuml[update_tag::insertkink]    << "\n";
  os << "# update probability delete kink                    : " << update_prob[update_tag::deletekink]    << "\t" << update_prob_cuml[update_tag::deletekink]    << "\n";
  os << "# update probability glue worm                      : " << update_prob[update_tag::glueworm]      << "\t" << update_prob_cuml[update_tag::glueworm]      << "\n";
  
  MyLatt->print_params(os);
  MyModel->print_params(os);
  
}



// Saves the state to the hdf5 file
void worm::save(alps::hdf5::archive & ar) const {
  
  /* Test configuration before saving:
   * if something went wrong since the last checkpoint, we don't
   * want to corrupt the last healthy checkpoint, too.
   */
  //test_conf();
  
  
  // Most of the save logic is already implemented in the base class
  alps::mcbase::save(ar);
    
  // We just need to add our own internal state
  ar["checkpoint/sweeps"] << sweeps;
}


void worm::load(alps::hdf5::archive & ar) {
  // Most of the load logic is already implemented in the base class
  alps::mcbase::load(ar);

  ar["checkpoint/sweeps"] >> sweeps;
  
  
  
  // Restore the internal state that came from parameters
  //Lx                    = size_t(parameters["Lx"]);
  //Ly                    = size_t(parameters["Ly"]);
  //Lz                    = size_t(parameters["Lz"]);
  //t_hop                 = double(parameters["t_hop"]);
  //U_on                  = double(parameters["U_on"]);
  //V_nn                  = double(parameters["V_nn"]);
  //mu                    = double(parameters["mu"]);
  beta                  = double(parameters["beta"]);
  E_off                 = double(parameters["E_off"]);
  C_worm                = double(parameters["C_worm"]);
  canonical             = int(parameters["canonical"]);
  seed                  = int(parameters["seed"]);
  //nmax                  = int(parameters["nmax"]);
  thermalization_sweeps = int(parameters["thermalization"]);
  Ntest = (unsigned long)(parameters["Ntest"]);
  Nsave = (unsigned long)(parameters["Nsave"]);
  Nmeasure = size_t(parameters["Nmeasure"]);
  
  
  // restore update probabilities????
  
}

void worm::print_conf(std::ostream& os) {
  os << "\n\n Printing operator string";
  for (SiteType i = 0; i < Nsites; i++) {
    os << "\nSite : " << i;
    int n = operator_string[i].size();
    int ii = 0;
    for (list<Element_t>::iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it, ++ii) {
      it->print();
      list<Element_t>::iterator assoc;
      for (SiteType j = zcoord[0][i]; j < zcoord[1][i]; j++)
      {
        assoc = it->get_assoc(j);
        std::cout << "\tassoc " << assoc->time() << "\t" << assoc->color();
      }
      if (ii > n) {
        os << "\n Error with list!\n";
        char ch; cin >> ch;
      }
    }
    os << "\n--------------------\n\n\n";
  }
  //os << "\n Worm head site : " << worm_head.site() << "\t time " << worm_head.time() << "\t dir "<< worm_dir << "\t rising " << worm_rising << "\n";
  //os << "\n Worm tail site : " << worm_tail.site() << "\t time " << worm_tail.time() ;
  os << "\n Worm head iterator time : "<< worm_head_it->time()  << "\n";
  os << "\n Worm tail iterator time : "<< worm_tail_it->time()  << "\n";
  os << "\n worm_at_stop = " << worm_at_stop << "\t worm_passes_nb_kink = " << worm_passes_nb_kink << "\t worm_measdensmat = " << worm_meas_densmat << "\n";
  //char ch; cin >> ch;
}
void worm::print_conf(){
  print_conf(std::cout);
}
void worm::print_update_statistics(std::ostream& os) const {
  std::vector<string> name;
  os << setprecision(10);
  name.push_back("INSERT WORM       ");
  name.push_back("MOVE WORM         ");
  name.push_back("INSERT KINK       ");
  name.push_back("DELETE KINK       ");
  name.push_back("GLUE WORM         ");
  os << "\n\n# UPDATE STATISTICS";
  os << "\n" << "# col 1 : all updates"
     << "\n" << "# col 2 : impossible updates"
     << "\n" << "# col 3 : rejected updates"
     << "\n" << "# col 4 : accepted updates"
     << "\n" << "# col 5 : acceptance factor with respect to Metropolis ratio only"
     << "\n" << "# col 6 : acceptance factor with respect to all attempts.\n";
  for (size_t i = 0; i < update_count; i++) {
    os << "\n" << name[i]
       << "\t" << update_statistics[total_attempted][i]
       << "\t" << update_statistics[impossible][i]
       << "\t" << update_statistics[rejected][i]
       << "\t" << update_statistics[accepted][i]
       << "\t" << update_statistics[accepted][i] / (update_statistics[total_attempted][i] - update_statistics[impossible][i])
       << "\t" << update_statistics[accepted][i] / update_statistics[total_attempted][i];
  }
  os << "\n\n";
  double s = 0;
  for (size_t i = 0; i < update_count; i++) {
    s += update_statistics[total_attempted][i];
  }
  os << "# Total number of steps : " << s << "\t or 10^" << log10(s) << "\n";
  
}



void worm::output_final(const alps::results_type<worm>::type & results,
                          alps::hdf5::archive& ar) {
  //std::string basename = alps::fs::remove_extensions(parameters["outputfile"]);

  
  // Jackknife analyses
  //using alps::accumulators::result_wrapper;
  //const result_wrapper& normed_A      = (with_reweight ? results["Obs_green_A"]    / results["Obs_green_A"].mean<vector<double> >()[0]    :  mZ_A * results["Obs_green_A"]     /results["Obs_sign_order0"]);
  //const result_wrapper& normed_B      = (with_reweight ? results["Obs_green_B"]    / results["Obs_green_B"].mean<vector<double> >()[0]    :  mZ_B * results["Obs_green_B"]     /results["Obs_sign_order0"]);
  
  //const result_wrapper& normed      = (with_reweight ? results["Obs_green"]    / results["Obs_green"].mean<vector<double> >()[0]    :  mZ * results["Obs_green"]     /results["Obs_sign_order0"]);
  //const result_wrapper& normed_k0   = (with_reweight ? results["Obs_green_k0"] / results["Obs_green_k0"].mean<vector<double> >()[0] :  mZ * results["Obs_green_k0"]  /results["Obs_sign_order0"]);
  //const result_wrapper& logged      = log(abs(normed_k0));
  //const result_wrapper& normed_k00   = (with_reweight ? results["Obs_green_k0_order0"] / results["Obs_green_k0_order0"].mean<vector<double> >()[0] :  mZ * results["Obs_green_k0"]  /results["Obs_sign_order0"]);
  //const result_wrapper& logged_order0      = log(abs(normed_k00));

  // save to archive
  //ar["/simulation/normed/Greenfun_A"]      << normed_A;
  //ar["/simulation/normed/Greenfun_B"]      << normed_B;
  //ar["/simulation/normed/log(Greenfun)"] << logged;

  // write ASCII files for convenience

}
