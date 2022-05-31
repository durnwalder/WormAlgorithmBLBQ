#include "worm.hpp"

#include <alps/params/convenience_params.hpp>

std::string worm::code_name() {
    return "Continuous-time quantum Monte Carlo simulation of the bilinear-biquadratic spin-1 model.";
}

// Defines the parameters for the worm simulation
void worm::define_parameters(parameters_type & parameters) {
    // If the parameters are restored, they are already defined
    if (parameters.is_restored()) {
        return;
    }
    
    alps::mcbase::define_parameters(parameters);
    alps::define_convenience_parameters(parameters)
        .description(worm::code_name())
        .define<size_t>("runtimelimit",         60,      "run time limit in seconds")
        .define<double>("beta",                 1.0,     "inverse temperature")
        .define<int>("sweeps",                  1000,    "maximum number of sweeps")
        .define<int>("thermalization",          200,     "number of sweeps for thermalization")
        .define<double>("E_off",                1.0,     "energy offset (technical parameter in move update")
        .define<double>("C_worm",               2.0,     "weight factor of worm configurations vs diagonal configurations (technical parameter in insertworm/glueworm updates")
        .define<int>("canonical",               -1,       "-1 for grand-canonical measurement; a positive number indicates the canonical measurement at this value")
        .define<unsigned long>("Ntest",         10000,    "test configuration after this number of updates "  )
        .define<unsigned long>("Nsave",         100000,   "save configuration after this number of updates"   )
        .define<size_t>("Nmeasure",             1,        "measure after this number of updates"   )
        .define<double>("p_insertworm",         1.0,      "update probability to insert worm in diagonal configuration")
        .define<double>("p_moveworm",           0.3,      "update probability to move worm head around")
        .define<double>("p_insertkink",         0.2,      "update probability to insert kink at position of worm head")
        .define<double>("p_deletekink",         0.2,      "update probability tp remove kink at position of worm head")
        .define<double>("p_glueworm",           0.3,      "update_probability to glue worm head and tail together and go back to diagonal configuration");

  
  
  Lattice_t::define_parameters(parameters);
  std::string LatticeClassifierName = parameters["lattice"].as<std::string>();
  if ( LatticeClassifierName == "chain_lattice" ) {
    Chain::define_custom_lattice_parameters(parameters);
  }
  else {
    throw std::runtime_error("Invalid Lattice \"" + LatticeClassifierName + "\"");
  }
  blbq_model::define_parameters(parameters);
  std::string ModelClassifierName = parameters["model"].as<std::string>();
  
}



// Creates a new simulation.
// We always need the parameters and the seed as we need to pass it to
// the alps::mcbase constructor. We also initialize our internal state,
// mainly using values from the parameters.
//worm::worm(parameters_type const & parms, std::size_t seed_offset) : alps::mcbase(parms, seed_offset)
worm::worm(parameters_type const & parameters, std::size_t seed_offset) : alps::mcbase(parameters, seed_offset)
    , runtimelimit(parameters["runtimelimit"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["thermalization"]))
    , total_sweeps(parameters["sweeps"])
    , beta(parameters["beta"])
    , E_off(parameters["E_off"])
    , C_worm(parameters["C_worm"])
    , canonical(parameters["canonical"])
    , seed(parameters["SEED"])
    , Ntest(parameters["Ntest"])
    , Nsave(parameters["Nsave"])
    , Nmeasure(parameters["Nmeasure"])
    , update_statistics(boost::extents[statistics_tag::statistics_count][update_tag::update_count])
    , MyGenerator(rd())
    , nbtime(ZC+1)
    , dtime(ZC+1)
    , nbval(ZC)
    , nbsite(ZC)
{
  MyGenerator.seed(seed);
  initialize_update_prob();
  
  for (size_t i=0; i < 10; i++) cout << rnd(MyGenerator) << "\n";
  
  std::string ClassifierName = parameters["lattice"].as<std::string>();
  if (ClassifierName == "chain_lattice" ) {
    MyLatt = unique_ptr<Lattice_t>(new Chain(parameters));
  }
  
  MyLatt->print();
  
  Nsites = MyLatt->get_Nsites();
  latt_dim = MyLatt->get_dim();
  nb.resize(boost::extents[Nsites][ ZC]);
  bond_index.resize(boost::extents[Nsites][ZC]);
  opposite_direction.resize(boost::extents[Nsites][ZC]);
  zcoord.resize(boost::extents[ZC][Nsites]);
  for (int n = 0; n < Nsites; n++ ) {
    if (n == 0 && !MyLatt->get_pbc(0)){
      zcoord[1][n] = 1;
    }
    else {
      zcoord[1][n] = 2;
    }
    if (n == Nsites-1 && !MyLatt->get_pbc(0)){
      zcoord[0][n] = 1;
    }
    else {
      zcoord[0][n] = 0;
    }
  }

  for (SiteIndex s = 0; s < Nsites; s++) {
    for (size_t j = 0 ; j < ZC;j++) {
      nb[s][j] = MyLatt->nb(s,j);
      opposite_direction[s][j] = MyLatt->get_opposite_direction(j, s);
      bond_index[s][j] = MyLatt->get_bond_index(s,j);
    }
  }
  
  std::string ModelClassifierName = parameters["model"].as<std::string>();
  if ( ModelClassifierName == "3color" ) {
    if(parameters.supplied("alpha"))
    {
      double a = parameters["alpha"].as<double>();
      if(!(a >= -1 && a <= 1))
      {
      throw std::runtime_error("Model "+ModelClassifierName+" is not negative-sign free for chosen parameter");
      }
    }
    else if (parameters.supplied("theta/pi"))
    {
      double tp = parameters["theta/pi"].as<double>();
      if(!(tp >= -0.75 && tp <= 0.25))
      {
      throw std::runtime_error("Model "+ModelClassifierName+" is not negative-sign free for chosen parameter");
      }
    }
    
    MyModel = unique_ptr<blbq_model>(new three_color(parameters, *MyLatt));
    measurements
      << alps::accumulators::LogBinningAccumulator<double >("Total_Energy")
      << alps::accumulators::LogBinningAccumulator<double >("Kinetic_Energy")
      << alps::accumulators::LogBinningAccumulator<double >("Potential_Energy")
      << alps::accumulators::LogBinningAccumulator<double >("Dimer_Order_Parameter")
      << alps::accumulators::LogBinningAccumulator<vector<double> >("String_Correlation")
      << alps::accumulators::LogBinningAccumulator<vector<double> >("SpinSpin_Correlation")
    ;
  }
  
#ifdef UNIFORM
  nmax = MyModel->get_nmax();
  nmin = MyModel->get_nmin();
  site_weight_diag.resize(nmax - nmin + 1);
  bond_weight_diag.resize(boost::extents[nmax - nmin + 1][nmax - nmin + 1]);
  site_weight_offdiag.resize(boost::extents[nmax - nmin + 1][nmax - nmin + 1]);
  bond_weight_offdiag.resize(boost::extents[nmax - nmin + 1][nmax - nmin + 1][nmax - nmin + 1][nmax - nmin + 1]);
  for (size_t n = 0; n < nmax - nmin + 1; n++ ) {
    site_weight_diag[n] = MyModel->site_weight_diag(0,n + nmin);
    for (size_t m = 0; m < nmax - nmin + 1; m++ ) {
      bond_weight_diag[n][m] = MyModel->bond_weight_diag( 0,n + nmin, m + nmin);
      site_weight_offdiag[n][m] = MyModel->site_weight_offdiag(0, n+nmin, m+nmin);
      for (size_t k = 0; k < nmax - nmin + 1; k++ ) {
        for (size_t l = 0; l < nmax - nmin + 1; l++ ) {
          bond_weight_offdiag[n][m][k][l] = MyModel->bond_weight_offdiag(0, n+nmin, m+nmin,k+nmin,l+nmin);
        }
      }
    }
  }
  std::cout << "Matrix Elements:" << std::endl;
  std::string element;
         for (size_t m1 = 0; m1 <= nmax; m1++)
        {
          for (size_t n1 = 0; n1 <= nmax; n1++)
          {
            for (size_t m2 = 0; m2 <= nmax; m2++)
            {
              for (size_t n2 = 0; n2 <= nmax; n2++)
              {
                if (n1 == n2 && m1 == m2)
                {
                  element = std::to_string(bond_weight_diag[n1][m1] + site_weight_diag[m1]*site_weight_diag[m2]);
                }
                else
                {
                  element = std::to_string(bond_weight_offdiag[n1][n2][m1][m2]);
                }
                std::string s(20-element.length(), ' ');
                std::cout <<  element.substr(0,5) << " ";
              }
            }
            std::cout << "\n";
          }
        }
#endif
  
  
  C_NBW = C_worm / (2*Nsites * beta);
  //C_NBW = C_worm / (4 * Nsites * beta);
  mWinding.resize(latt_dim);
  //trans.resize(zcoord + 1);
  //dtime.resize(zcoord + 1);
  //nbtime.resize(zcoord + 1);
    
  operator_string.resize(Nsites);
  //site_it.resize(Nsites);
  dummy_it.resize(Nsites);
  //for (SiteType i = 0 ; i < Nsites; i++) {
  //  site_it[i] = operator_string[i].begin();
  //}
 

  state.resize(Nsites);
  av_dns.resize(Nsites);
  av_state.resize(Nsites);
  av_state_sq.resize(Nsites);
  Potential_Energy_local.resize(Nsites-!MyLatt->get_pbc(0));
  if (MyLatt->get_pbc(0))
  {
    string_correlation.resize(ceil(Nsites/2));
    spinspin_correlation.resize(ceil(Nsites/2)); 
  }
  else 
  {
    string_correlation.resize(Nsites);
    spinspin_correlation.resize(Nsites); 
  }
  MCdiag = 0;
  
}







void worm::initialize_update_prob() {
  //initialize_update_prob(insertworm, update_count);
  initialize_update_prob(insertworm, moveworm);
  initialize_update_prob(moveworm, update_count);
}


void worm::initialize_update_prob(update_tag begin, update_tag end) {
  double norm = 0.;
  for (size_t upd = begin; upd < end; ++upd) {
    if (parameters.defined("p_" + update_names[upd])) {
      
      update_prob[upd] = double(parameters["p_" + update_names[upd]]);
      if (update_prob[upd] < 0) {
        throw std::runtime_error("Negative update probability: "
                                         + update_names[upd]);
      }
    }
    else {
      update_prob[upd] = 0.;
    }
    norm += update_prob[upd];
  }
  if (abs(norm - 1.) > 1e-10) {
    std::cerr << "Update probabilities do not add up to one.\n"
              << "Renormalizing..." << std::endl;
    for (size_t upd = begin; upd < end; ++upd) {
      update_prob[upd] /= norm;
      parameters["p_" + update_names[upd]] = update_prob[upd];
    }
  }
  norm = 0.;
  for (size_t upd = begin; upd < end; ++upd) {
    norm += update_prob[upd];
    update_prob_cuml[upd] = norm;
  }
}



void worm::initialize() {
  worm_diag = true;
  worm_at_stop = 0;
  worm_passes_nb_kink = 0;
  worm_dtime = 0.;
  new_measurement = true;
  for (int k = 0; k < latt_dim; k++) mWinding[k] = 0.;
  nrvertex.resize(Nsites-!MyLatt->get_pbc(0),0); 
  nrv_diff.resize(Nsites-!MyLatt->get_pbc(0)-1,0); 
  nrvertex_type1 = 0;
  nrvertex_type2 = 0;
  for (SiteType s = 0; s < Nsites; s++) {
    state[s] = 1;
    std::cout <<"# Initial state " << s << "\t" << state[s] << "\n";
  }
 
  Nprtcls = 0.;
  for (SiteIndex i = 0; i < Nsites; i++) {
    Element_t new_elem(state[i],state[i], i,  beta, 0);
    dummy_it[i] = operator_string[i].insert(operator_string[i].end(), new_elem);
    Nprtcls += dummy_it[i]->before() * beta;
  }
  std::cout << "# Setting initial associations ... ";
  for (SiteIndex i = 0; i < Nsites; i++) {
    for (std::size_t j = zcoord[0][i] ; j < zcoord[1][i]; j++) {
      dummy_it[i]->set_assoc(j, dummy_it[nb[i][j] ] );
    }
  }
  worm_head_it = dummy_it[0];
  worm_tail_it = dummy_it[0];
  
  std::cout << "...done. Computing potential energies...";
  Epot_tot = calc_potential_energy_nb() + calc_potential_energy_loc();
  std::cout << "...done.\n";
   cout << "# Potential Energy tot : " << Epot_tot << endl;
   std::cout << "\n# Finished intializing.\n";
}




void worm::measure() {
  sweeps++;
  if (sweeps > thermalization_sweeps) {
    double Ep = calc_potential_energy_measure();
    double dop = calc_dop();
            // this also sets state and number_of_particles observables
    double Ek = (nrvertex_type1+nrvertex_type2)/beta * (-1.);
    for (SiteType j = 0; j < Nsites-!MyLatt->get_pbc(0)-1; j++) {
      nrv_diff[j] = abs((nrvertex[j+1] - nrvertex[j]));
    }

    measurements["String_Correlation"] << string_correlation;
    measurements["SpinSpin_Correlation"] << spinspin_correlation;
    measurements["Kinetic_Energy"] << Ek;
    measurements["Potential_Energy"] << Ep;
    measurements["Total_Energy"] << Ep + Ek;
    measurements["Dimer_Order_Parameter"] << dop;
  }
  
}


// Returns a number between 0.0 and 1.0 with the completion percentage
double worm::fraction_completed() const {
  //std::cout << "# Fraction completed : " << (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps)) << "\n";
  if ( ( sweeps - thermalization_sweeps ) % 100000 == 0 ) {
    std::cout << "# Fraction completed : " << sweeps << "\t" << thermalization_sweeps << "\t" << total_sweeps << "\n";
    std::cout << "# Fraction completed : " << (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps)) << "\n";
  }

  return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}




void worm::test_conf() {
  for (SiteIndex i = 0; i < Nsites; i++) {
    for (list<Element_t>::iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it) {
      if (MyModel->range_fail(it->before() )) {
        cerr << "\n# TEST_CONF :  Error with density  (before) << " << i << "\t" << it->time() << "\t" << it->before() << "\t" << it->after() << "\n";
        throw exception();
      }
      if (MyModel->range_fail(it->after() )) {
        cerr << "\n# TEST_CONF : Error with density  (after) << " << i << "\t" << it->time() << "\t" << it->before() << "\t" << it->after() << "\n";
        throw exception();
      }
      list<Element_t>::iterator itn = it;
      ++itn; if (itn == operator_string[i].end()) itn = operator_string[i].begin();
      if ( (it->after() == it->before()) && (it->color() != 0)) {
        cerr << "\n# TEST_CONF : error in configuration (before==after) on site  " << i << endl;
        throw exception();
      }
      if ( (it->color() == 0) && ((it->after() != it->before()))) {
        cerr << "\n# TEST_CONF : error in configuration (dummy) on site  " << i << endl;
        throw exception();
      }
      if (it->after() != itn->before() ) {
        cerr << "\n# TEST_CONF : error in configuration (diag) on site  " << i << endl;
        throw exception();
      }
      if (it->time() > itn->time() && !(it == dummy_it[i] ) ) {
        cerr << "\n# TEST_CONF : chronology broken on site  " << i << endl;
        throw exception();
      }
      if (!MyLatt->get_pbc(0) && it->color()>0 && std::max(i, it->link())-std::min(i, it->link())>1)
      {
        cerr << "\n# TEST_CONF : Bond  over boundary for open boundary condition " << i << endl;
        throw exception();
      }
      
      if (it->color() == 1) {
        bool b = false;
        bool b2 = false;
        for (std::size_t j=zcoord[0][i] ; j < zcoord[1][i]; j++ ) {
          if ( b && abs(it->time() - it->get_assoc(j)->time() ) < 1e-16 && it->get_assoc(j)->color() == 1 ) b2 = true;
          if (abs(it->time() - it->get_assoc(j)->time() ) < 1e-16 && it->get_assoc(j)->color() == 1) b = true;
          if ( (it->link() == nb[i][j] ) && (it->get_assoc(j)->link() == i) && (it->get_assoc(j)->time() != it->time())) {
            cerr << "\n# TEST_CONF : site is linked and has different association time; site " << i <<  " link " << it->link() << " at time " << it->time() << " and " <<  it->get_assoc(j)->time() << endl;
            throw exception();
          }
        }
        if (b2) {
          cerr << "\n# TEST_CONF : error in configuration, found two links at same time on site " << i << " for time " << it->time() << "\t links : " << it->get_assoc(0)->time() << "\t" << it->get_assoc(1)->time();
          throw exception();
        }
        if (!b) {
          cerr << "\n# TEST_CONF : error in configuration, found no link at same time on site " << i << " for time " << it->time() << "\t links : " << it->get_assoc(0)->time() << "\t" << it->get_assoc(0)->color()  << "\t" <<  it->get_assoc(1)->time() << "\t" << it->get_assoc(1)->color() << endl;
          throw exception();
        }
      }
      if (it->color() != 0) {         // extensive check on associations
        for (std::size_t j=zcoord[0][i]; j <zcoord[1][i]; j++ ) {
          if (it->get_assoc(j)->time() < it->time() ) {
            cerr << "\n# TEST_CONF : assoc has earlier time than current link for current link on site " << i << " at time " << it->time() << " and nb " <<j << endl;
            throw exception();
          }
          list<Element_t>::iterator itp = it->get_assoc(j);
          SiteIndex adj_site = nb[i][j];
          if (itp == operator_string[adj_site].begin()) itp = operator_string[adj_site].end();
          --itp;
          if (!(itp == dummy_it[adj_site]) && (itp->time() > it->time())  ) {
            cerr << "\n# TEST_CONF : there is an earlier element which should have been associated to cursite " << i << " at time " << it->time() << " and nb dir  " <<j << " site " <<  adj_site << "\t times it->get_assoc(j) and its predecessor : " << it->get_assoc(j)->time() << "\t" << itp->time() << endl;
            throw exception();
          }
        }
      } // ... if it->color() != 0
    } // ... iterate over configuration list
  } // iterate over all sites
  if (!worm_diag) {
    if (worm_head_it->color() > -1 ) {
      cerr << "\n# TEST_CONF : worm_head_it points to wrong element : " << worm_head_it->color()   << endl;
      throw exception();
    }
    if (worm_tail_it->color() > -1 ) {
      cerr << "\n# TEST_CONF : worm_tail_it points to wrong element : " << worm_tail_it->color()   << endl;
      throw exception();
    }
  }
  if (worm_meas_densmat) {
    if (worm_tail_it->time() != worm_head_it->time() ) {
      cerr << "\n# TEST_CONF : worm_meas_densmat is set but worms are at unequal times : " << worm_tail_it->time()  << "\t" << worm_head_it->time() << endl;
      throw exception();
    }
  }
  if ( (worm_at_stop == +1)  && (!worm_passes_nb_kink) && (!worm_meas_densmat) ) {
    list<Element_t>::iterator itp = worm_head_it;
    if (itp == operator_string[worm_head_it->link()].begin()) itp = operator_string[worm_head_it->link()].end();
    --itp;
    if ( ! ( (itp->time() == worm_head_it->time()) || ( abs(itp->time()-beta + worm_head_it->time()) < 1e-10 ) ) )  {
      cerr << "\n# TEST_CONF : worm_at_stop == +1 but times are unequal times : " << worm_head_it->time()  << "\t" << itp->time() << "\t site : " << worm_head_it->link() << endl;
      throw exception();
    }
  }
  if ( (worm_at_stop == -1) && (!worm_passes_nb_kink) && (!worm_meas_densmat) ) {
    list<Element_t>::iterator itn = worm_head_it;
    ++itn;
    if (itn == operator_string[worm_head_it->link()].end()) itn = operator_string[worm_head_it->link()].begin();
    if (itn->time() != worm_head_it->time()) {
      cerr << "\n# TEST_CONF : worm_at_stop == -1 but times are unequal times : " << worm_head_it->time()  << "\t" << itn->time() << "\t site : " << worm_head_it->link() << endl;
      throw exception();
    }
  }
  if (worm_passes_nb_kink) {
    bool b = false;
    for (std::size_t j=zcoord[0][worm_head_it->link()]; j <zcoord[1][worm_head_it->link()]; j++ ) {
      SiteIndex adj_site = nb[worm_head_it->link()][ j];
      for (list<Element_t>::iterator it = operator_string[adj_site].begin(); it != operator_string[adj_site].end(); ++it) if (it->time() == worm_head_it->time()) b = true;
    }
    if (!b) {
      cerr << "\n# TEST_CONF : worm_passes_nb_kink is set but none is found; worm head time and site : " << worm_head_it->time() << "\t" << worm_head_it->link() << endl;
      throw exception();
    }
  }
  double Ep = calc_potential_energy_loc() + calc_potential_energy_nb();
  if (is_not_close(Ep, Epot_tot, 1e-8)) {
    cerr << "# Potential total energies do not match " << Ep << "\t" << Epot_tot << "\n";
    throw exception();
  }
  else {
#ifdef DEBUGMODE
    cout << "# Potential energies OK : " << Ep << "\t" << Epot_tot << "\t nrvertex : " << nrvertex << "\t Nprtcls : " << Nprtcls / Nsites / beta << "\n";
#endif
  }
}
  


double worm::calc_potential_energy_loc() {
  double En = 0.;
  double dt = 0.;
  for (int isite = 0; isite < Nsites; isite++) {
    double t_old = 0.;
    list<Element_t>::iterator it1 = operator_string[isite].begin();
    int n1 = it1->before();
    for (;;) {
      dt = it1->time() - t_old;
#ifdef UNIFORM
      En += site_weight_diag[n1 - nmin] * dt;
#else
      En += MyModel->site_weight_diag(isite, n1) * dt;
#endif
      //En += (0.5 * U_on * n1 * (n1-1) - mu * n1) * dt;
      
      n1 = it1->after();
      t_old = it1->time();
      ++it1;
      if (it1 == operator_string[isite].end()) break;
    }
  }
  return (En / beta );
}

double worm::calc_potential_energy_nb() {
  double En = 0.;
  double dt = 0.;
  for (int isite = 0; isite < Nsites-!MyLatt->get_pbc(0); isite++) {
    for (int iz = 0; iz < latt_dim; iz++) {
      list<Element_t>::iterator it1 = operator_string[isite].begin();
      int n1 = it1->before();
      //int jsite = nb[isite][ iz];
      SiteType jsite = nb[isite][ iz];
      list<Element_t>::iterator it2 = operator_string[jsite].begin();
      int n2 = it2->before();
      double t_old = 0.;
      for (;;) {
        for (;;) {
          if (it2 == operator_string[jsite].end()) break;
          if (it2->time() > it1->time()) break;
          dt = it2->time() - t_old;
          //En += dt * n1 * n2;
#ifdef UNIFORM
          En += bond_weight_diag[ n1][ n2] * dt;
#else
          En += MyModel->bond_weight_diag(bond_index[isite][ iz], n1, n2) * dt;
#endif
          //cout << "# calc_pot_nb (2) : isite, jsite, it1->time, it2->time, n1, n2, t_old, dt : " << isite << "\t" << jsite << "\t" << it1->time() << "\t" << it2->time() << "\t" << n1 << "\t" << n2 << "\t" << t_old << "\t" << dt << "\n";
          t_old = it2->time();
          n2 = it2->after();
          ++it2;
        }
        dt = it1->time() - t_old;
        //En += dt * n1 * n2;
#ifdef UNIFORM
        En += bond_weight_diag[n1][ n2] * dt;
#else
        En += MyModel->bond_weight_diag(bond_index[isite][ iz], n1, n2) * dt;
#endif
        //if (it2 != operator_string[jsite].end()) cout << "# calc_pot_nb (1) : isite, jsite, it1->time, it2->time, n1, n2, t_old, dt : " << isite << "\t" << jsite << "\t" << it1->time() << "\t" << it2->time() << "\t" << n1 << "\t" << n2 << "\t" << t_old << "\t" << dt << "\n";
        n1 = it1->after();
        t_old = it1->time();
        ++it1;
        if (it1 == operator_string[isite].end()) break;
      }
    }
  }
  //return (En * V_nn / beta);
  return (En  / beta);
}
double worm::calc_dop(){
  double dp = 0;
  double Eold = 0;
  double En = 0;
  for (SiteType j = 0; j < Nsites-!MyLatt->get_pbc(0); j++) {
    Eold = En;
    En = 0;
    int n = dummy_it[j]->before();
    En +=  MyModel->site_weight_diag(j, n,0);
    En += bond_weight_diag[ n][ dummy_it[nb[j][0]]->before()];
    if (j>0)
    {
      dp += abs(En - Eold + (nrvertex[j+1] - nrvertex[j])/beta);
    }
    Potential_Energy_local[j] = En;
  }
  //return (0.5*U_on*En1 + En2 + En3);
  return dp/(Nsites-1-!MyLatt->get_pbc(0));
}
double worm::calc_potential_energy_measure()
{
  //int32_t En1 = 0;
  //double En2 = 0.;
  //double En3 = 0.;
  double En = 0.;
  number_of_particles=0;
  for (SiteType j = 0; j < Nsites-!MyLatt->get_pbc(0); j++) {
    int n = dummy_it[j]->before();
    
    //En1 += n*(n-1);
    //En2 += (mu_eff[j]+mu)*n;
    //En2 += -mu*n;
#ifdef UNIFORM
    En +=  MyModel->site_weight_diag(j, n,0) ;
    int snb = 0;
    for (std::size_t iz = 0; iz < latt_dim; iz++) {
      En += bond_weight_diag[ n][ dummy_it[nb[j][iz]]->before()];
    }
#else
    En +=  MyModel->site_weight_diag(j, n,0) ;
    int snb = 0;
    for (std::size_t iz = 0; iz < ZC; iz++) {
      //snb += dummy_it[nb[j][i]]->before();
      En += MyModel->bond_weight_diag(bond_index[j][ iz], n, dummy_it[nb[j][iz]]->before());
    }
#endif
    //En3 += 0.5 * V_nn * n * snb;
    number_of_particles += n;
    state[j] = n;
  }
  //return (0.5*U_on*En1 + En2 + En3);
  return (En);
}



void worm::find_assoc_insert(const SiteIndex cursite, list<Element_t>::iterator  it, const int shift) {
  // it is the iterator to the newly insterted element
  // memory for its associations has been allocated already but the associations must be set correctly here
  // and we need to check on the neighbors if their associations need to be changed
  //std::vector<SiteIndex> nbs = MyModel->couplers(cursite);
  
  //MyModel->couplers(cursite, nbsite);
  
  // part 1 : iterators for the new element located on site cursite
  //std::cout<< "# Welcome to find_assoc_insert " << cursite << "\t" << *it << "\n";
  // we go one element up where we by assumption have a properly associated element
  list<Element_t>::iterator ito, itl, itp, itw;
  ito= it;
  ++ito;
  if (ito == operator_string[cursite].end()) ito = operator_string[cursite].begin();
  
  double t0 = it->time();
  it->time( t0 + shift * dtol*5);
 
  for (std::size_t j = zcoord[0][cursite] ; j < zcoord[1][cursite]; j++) {
    //SiteIndex s = nbsite[j];
    SiteIndex s = nb[cursite][j];
    itl = ito->get_assoc(j);
    itp = itl;
    if (itp == operator_string[s].begin()) itp = operator_string[s].end();
    --itp;
    it->set_assoc(j, itl);
    while (!t_between(it->time(), itp->time(), itl->time())) {
      itl = itp;
      if (itp == operator_string[s].begin()) itp = operator_string[s].end();
      --itp;
      it->set_assoc(j, itl);
    }
  }
  
  
  // part 2 : iterators on the neighbors might have to be set newly in the "opposite" direction
  for (std::size_t j = zcoord[0][cursite]  ; j < zcoord[1][cursite]; j++) {
    SiteIndex s = nb[cursite][j];
    //SiteIndex s = nbsite[j];
    std::size_t oppdir = opposite_direction[s][j];
    //list<Element_t>::iterator itl;
    itl = it->get_assoc(j);
    //list<Element_t>::iterator itp = itl;
    itp = itl;
    if (itp == operator_string[s].begin()) itp = operator_string[s].end();
    --itp;
    if (zcoord[0][s] <= oppdir && zcoord[1][s] > oppdir)
    {
    itw = itp->get_assoc(oppdir);
    if (itl->time() == it->time()) {
      itl->set_assoc(oppdir, it);
      //std::cout <<  setprecision(16) << "# find_assoc_insert setting itl asso to it in oppdir : " << cursite << "\t" << j << "\t" << s << "\t" << oppdir << "\t" << it->time() << "\t" << itl->time() << "\t shift : " << shift << "\t" << shift*dtol << "\t" << itl->time() - it->time() << "\n";
    }
    while ((itp->time() != itw->time())  && ( t_between(it->time(), itp->time(), itw->time()) )) {
      itp->set_assoc(oppdir, it);
      if (itp == operator_string[s].begin()) itp = operator_string[s].end();
      --itp;
      itw = itp->get_assoc(oppdir);
    }
    }
  }
  
  it->time(t0);
}




void worm::find_assoc_delete(const SiteIndex cursite, list<Element_t>::iterator  it) {
  if (operator_string[cursite].size() == 1) {
    std::cerr <<"# How to erase a list of length one when there must be a dummy?\n";
    return;
  }
 
  list<Element_t>::iterator itn=it;
  ++itn;
  if (itn == operator_string[cursite].end()) itn = operator_string[cursite].begin();
  list<Element_t>::iterator itp = it;
  if (itp == operator_string[cursite].begin()) itp=operator_string[cursite].end();
  --itp;
  // on nb sites, go down until you find an element that points at the current element
  //MyModel->couplers(cursite, nbsite);
  list<Element_t>::iterator it_end, it_begin, itt;
  for (std::size_t j = zcoord[0][cursite] ; j < zcoord[1][cursite]; j++) {
    //SiteIndex s = nbsite[j];
    SiteIndex s= nb[cursite][j];
    std::size_t oppdir = opposite_direction[s][j];
    it_end=it->get_assoc(j);
    it_begin=itp->get_assoc(j);
    itt=it_begin;
    if ((itt == it_end) && (itt == dummy_it[s])) {
      for (list<Element_t>::iterator itt=operator_string[s].begin(); itt != operator_string[s].end(); ++itt) {
        if (itt->get_assoc(oppdir) == it) itt->set_assoc(oppdir, itn);
      }
    }
    else {
      while (itt != it_end) {
        if (itt->get_assoc(oppdir) == it) itt->set_assoc(oppdir, itn);
        ++itt;
        if (itt==operator_string[s].end()) itt=operator_string[s].begin();
      }
      if (it_end->get_assoc(oppdir) == it) it_end->set_assoc(oppdir, itn);
    } // else
  } //for
}


void worm::update_spinspin_corr(){
  int head_site = worm_head_it->link();
  bool Ttype[3]  = { true, true , true };
  if (worm_head_it->color()>0){
    cout << "worm_head_it->color(): " <<worm_head_it->color()<< endl;
  }
  Ttype[worm_head_it->before()] = false;
  Ttype[worm_head_it->after()] = false;
  int tail_site = worm_tail_it->link();
  if (worm_tail_it->color()>0){
    cout << "worm_tail_it->color(): " <<worm_tail_it->color()<< endl;
  }
  Ttype[worm_tail_it->before()] = false;
  Ttype[worm_tail_it->after()] = false;
  int type;
  if (Ttype[0] && !Ttype[1] && !Ttype[2])
  {
    type = 0;
  }
  else if (Ttype[1] && !Ttype[0] && !Ttype[2])
  {
    type = 1;
  }
  else if (Ttype[2] && !Ttype[0] && !Ttype[1])
  {
    type = 2;
  }
  else
  {
    cout << "Error in update Spin Spin Corr!" << endl;
  }
  if (abs(head_site-tail_site)>1)
  {
    list<Element_t>::iterator head_it = operator_string[head_site].begin();
    list<Element_t>::iterator tail_it = operator_string[tail_site].begin();
    SiteType head_link;
    SiteType tail_link;
    if (head_site>tail_site)
    {
      head_link = head_site - 1;
      tail_link = tail_site + 1;
    }
    else
    {
      head_link = head_site + 1;
      tail_link = tail_site - 1;
    }
    spinspin_corr = 0;
    bool head_end = false;
    bool tail_end = false;
    double tau = 0;
    double delta;
    int sign = 1;
    //do not start at head but start at tau=0 with dummy
    for (SiteType j = std::min(head_site,tail_site)+1 ; j < std::max(head_site,tail_site); j++) 
    {
      if (dummy_it[j]->before() != type) sign = -sign;
    }
    while (!(head_end  && tail_end))
    {
      while (!(head_it->link() == head_link && (head_it->before() == type || head_it->after() == type)) && !head_end)
      {
        head_it++;
        if (head_it == operator_string[head_site].end()) head_end = true;
      }
      while (!(tail_it->link() == tail_link && (tail_it->before() == type || tail_it->after() == type)) && !tail_end)
      {
        tail_it++;
        if (tail_it == operator_string[tail_site].end()) tail_end = true;
      }
      if (!(head_end  && tail_end))
      {
        if ((tail_it->time()<head_it->time() || head_end) && !tail_end)
        {
          delta = tail_it->time() - tau;
          tau = tail_it->time();
          tail_it++;
          if (tail_it == operator_string[tail_site].end()) tail_end = true;
        }
        else{
          delta = head_it->time() - tau;
          tau = head_it->time();
          head_it++;
          if (head_it == operator_string[head_site].end()) head_end = true;
        }
        spinspin_corr += sign*delta;
        sign = -1*sign;
      }
      else {
        spinspin_corr = (beta-tau)*sign;
      }
    }
  }
  else
  {
    spinspin_corr = beta;
  }
  spinspin_corr = spinspin_corr/beta;
  if (type != 2)
  {
    spinspin_corr = 0;
  }
  /*
  if (abs(head_site-tail_site) == 2){
    cout << "distance = 2 spinspin_corr =" << spinspin_corr << endl;
    print_conf();
  }
  if (abs(head_site-tail_site) == 3){
    cout << "distance = 3 spinspin_corr =" << spinspin_corr << endl;
    print_conf();
  }
  */
}

