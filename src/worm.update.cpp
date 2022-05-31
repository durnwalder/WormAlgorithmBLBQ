#include "worm.hpp"
#include <stdlib.h> 

void worm::update() {
  
  stringstream lognr;
  lognr << "qc_log" << int(MCdiag);
  //ofstream outlog(lognr.str().c_str() );
  //for (size_t i =0; i < Nloop; i++) {
  size_t Nupd = 0;
  do {
    //double q = random();
    double q = rnd(MyGenerator);
    int a;
    if (worm_diag) {
      if (q < update_prob_cuml[insertworm]) {
        a = INSERTWORM();
        update_statistics[statistics_tag::total_attempted][update_tag::insertworm] += 1.;
        update_statistics[a][update_tag::insertworm] += 1.;
      }
    }
    else {
      if (q < update_prob_cuml[moveworm]) {
        a = MOVEWORM();
        update_statistics[statistics_tag::total_attempted][update_tag::moveworm] += 1.;
        update_statistics[a][update_tag::moveworm] += 1.;
      }
      else if (q < update_prob_cuml[insertkink]) {
        a = INSERTKINK();
        update_statistics[statistics_tag::total_attempted][update_tag::insertkink] += 1.;
        update_statistics[a][update_tag::insertkink] += 1.;
      }
      else if (q < update_prob_cuml[deletekink]) {
        a = DELETEKINK();
        update_statistics[statistics_tag::total_attempted][update_tag::deletekink] += 1.;
        update_statistics[a][update_tag::deletekink] += 1.;
      }
      else if (q < update_prob_cuml[glueworm]) {
        a = GLUEWORM();
        update_statistics[statistics_tag::total_attempted][update_tag::glueworm] += 1.;
        update_statistics[a][update_tag::glueworm] += 1.;
      }
    }
//#ifdef DEBUGMODE
    //print_conf(std::cout);
    /*
    try {
      test_conf();
    }
    catch (const exception& e) {
      cerr << e.what() << endl;
      print_conf(std::cout);
      exit(1);
    }
    */
    //char ch; std::cin >> ch;
//#endif
    
    //update_hist();
    //outlog << q << " a = " << a << " Epot = " << Epot_tot << "\t nrvertex = " << nrvertex << "\n";
    
    if (worm_diag) {
      mZ += 1.;
      MCdiag += 1;
      //std::cout << "# Diagonal configuration " << Nprtcls/beta/ Nsites << "\t" << nrvertex << "\n";
      //char chd; cin>>chd;
    }
    else {
      if (worm_meas_densmat && (sweeps >= thermalization_sweeps)) {
	    SiteType dd = MyLatt->dist(worm_head_it->link(), worm_tail_it->link());
      if (MyModel->get_name() == "original")
      {
        string_correlation[dd] += spinspin_corr;
        spinspin_correlation[dd] ++;
      }
      else
      {
        string_correlation[dd]++;
        spinspin_correlation[dd] += spinspin_corr;
      }
      }
      //measure_density_matrix();
    }
    Nupd++;
    if (Nupd == 10000) {
      //print_conf(std::cout);
      Nupd = 0;
    }
    
  } while (!worm_diag);
  //outlog.close();
  //if (MCdiag == 10) exit(0);
  //sweeps++;
}




// Changes compared to old code : insert both head and tail in the configuration list


int worm::INSERTWORM() {
#ifdef DEBUGMODE
  if (!worm_diag) {
    cerr << "\nNon diagonal configuration in INSERTWORM ? \n" ; 
  }
  std::cout << "# Welcome to INSERTWORM\n";
  print_conf(std::cout);
#endif
  worm_dtime = 0;
  mZ_dns += 1.;
  //SiteIndex isite=static_cast<SiteIndex>(random()*Nsites);
  SiteIndex isite=static_cast<SiteIndex>(rnd(MyGenerator)*Nsites);
  //double start_time=random()*beta;
  double start_time=rnd(MyGenerator)*beta;
#ifdef DEBUGMODE
  cout << "INSERTWORM : isite " << isite << "\t start_time "<< start_time << "\n";
#endif
  // find suitable place to make a worm pair
  list<Element_t>::iterator cycle_it = dummy_it[isite];
  list<Element_t>::iterator prev_it = cycle_it;
  if (operator_string[isite].size() != 1) // dummy element always there
  {
    
    if (prev_it ==  operator_string[isite].begin()) prev_it = operator_string[isite].end();
    --prev_it;

    // you cannot  create a worm at the same place where the previous one was removed...
    // you really need to go around in space-time. maybe for insulating phases 
    // the next couple of lines are time consuming
    while (!t_between(start_time, prev_it->time(), cycle_it->time()))
    {
      cycle_it = prev_it;
      //site_right(isite);
      if (prev_it == operator_string[isite].begin()) prev_it = operator_string[isite].end();
      --prev_it;
    }
  }
  
  //worm_it = site_it[isite];
  int worm_dir = ( (rnd(MyGenerator) < 0.5) ? 1 : -1 );
  int start_occ = cycle_it->before(); //test
  //StateType n_mid = n_out + ((random() < 0.5) ? 1 : -1);
  int start_mid = (start_occ + (rnd(MyGenerator) < 0.5 ? 1 : 2 ) ) % 3;
  int clr = abs(start_mid - start_occ);
#ifdef UNIFORM
  if (!range_check(start_mid)) return impossible;
#else
  if ( MyModel->range_fail(n_mid) ) return impossible;
#endif
 
  
  //int worm_dir = ( (random() < 0.5) ? 1 : -1 );
  
#ifdef UNIFORM
  double wgt = site_weight_offdiag[start_occ][ start_mid]; //test
  double wgt_sq = wgt * wgt;
#else
  double wgt = MyModel->site_weight_offdiag(isite, n_out, n_mid);
  double wgt_sq = wgt * wgt;
#endif
  
  //double ratio =  4* C_NBW * Nsites * beta * wgt_sq * update_prob[glueworm] / update_prob[insertworm];
  
  double ratio =  4* C_NBW * Nsites * beta * update_prob[glueworm] / update_prob[insertworm];
  
#ifdef DEBUGMODE
  cout << "INSERTWORM ratio " <<  ratio << "\t" <<  4* C_NBW * Nsites * beta << "\t" <<  wgt << "\t" <<  update_prob[glueworm] / update_prob[insertworm] << "\n";
  //char ch; cin >> ch;
#endif
  if (Metropolis(ratio)) {
#ifdef DEBUGMODE
      cout << "INSERTWORM accepted " << cycle_it->time() << "\t" << passed_dummy << "\t ratio : " << ratio << "\n";
#endif
    Element_t lower_kink(start_occ, start_mid, isite, start_time, -clr );

    Element_t upper_kink(start_mid, start_occ, isite, start_time, -clr );
    worm_head_it = operator_string[isite].insert(cycle_it, upper_kink);
    find_assoc_insert(isite, worm_head_it,0);
    worm_tail_it = operator_string[isite].insert(worm_head_it, lower_kink);
    find_assoc_insert(isite, worm_tail_it, -1);
    
    if (worm_dir == 1) {
      //worm_tail.init(n_out, n_mid, isite, start_time);
      //worm_head.init(n_mid, n_out, isite, start_time);
      worm_at_stop = 1;
      worm_passes_nb_kink = 0;
    }
    else {
      std::swap(worm_head_it, worm_tail_it);
      //worm_head.init(n_out, n_mid, isite, start_time);
      //worm_tail.init(n_mid, n_out, isite, start_time);
      worm_at_stop = -1;
      worm_passes_nb_kink = 0;
    }
#ifdef DEBUGMODE
    cout << "\nInitial head (it): " << *worm_head_it;
    cout << "\nInitial tail (it): " << *worm_tail_it;
    
        //char ch; cin >> ch;
#endif
    worm_diag = 0;
    new_measurement = true;
    worm_meas_densmat = true;
    //update_spinspin_corr();
    return accepted;
  }
  else {
    return rejected;
  }
}










int worm::MOVEWORM() {
#ifdef DEBUGMODE
  if (worm_diag) cerr << "\n# Configuration is diagonal in MOVEWORM ? " << "\n";
  std::cout << "# Welcome to MOVEWORM.\n";
#endif
  
  //int dir = (random() < 0.5 ? 1 : -1);
  int dir = (rnd(MyGenerator) < 0.5 ? 1 : -1);
  
  SiteIndex isite = worm_head_it->link();
  list<Element_t>::iterator itn = worm_head_it;
  ++itn; if (itn == operator_string[isite].end() ) itn = operator_string[isite].begin();
  list<Element_t>::iterator itp = worm_head_it;
  if (itp == operator_string[isite].begin() ) itp = operator_string[isite].end(); --itp;
  
      
   if (dir == +1) {
     if (worm_at_stop == -1) {
       // several cases : dummy, worm_tail, density matrix, interaction
       if (worm_meas_densmat) {
         if ( itn->color() < 0 ) return impossible; // cannot pass worm tail
         
         for (size_t j=zcoord[0][isite] ; j < zcoord[1][isite]; j++ ) {
            if (nb[isite][j] == worm_tail_it->link() ) {
              size_t oppdir = opposite_direction[nb[isite][j] ] [j]  ;
              list<Element_t>::iterator it_nb = worm_head_it->get_assoc(j);
              list<Element_t>::iterator it_up = it_nb;
              ++it_up; if (it_up == operator_string[worm_tail_it->link()].end() ) it_up = operator_string[worm_tail_it->link()].begin();
              worm_head_it->set_assoc(j, it_up);
              it_nb->set_assoc(oppdir, worm_head_it);
              break;
            }
         }
         
         //bool is_neighbor = false;
         //for (size_t j=0; j <ZC; j++ ) if (nb[worm_head_it->link()][j] == worm_tail_it->link()) is_neighbor = true;
         //if (is_neighbor) {
         // find_assoc_delete(isite, worm_head_it);
         // find_assoc_insert(isite, worm_head_it, +1);
         //}
         worm_at_stop = +1;
       }
       else if (worm_passes_nb_kink) {
         // do nothing ? just pass
        
         list<Element_t>::iterator it_nb, it_up;
         size_t j, oppdir;
         for (j=zcoord[0][worm_head_it->link()] ; j < zcoord[1][worm_head_it->link()]; j++ ) {
           SiteIndex adj_site = nb[worm_head_it->link()][j];
           it_nb = worm_head_it->get_assoc(j) ;
           it_up = it_nb;
           ++it_up; if (it_up == operator_string[adj_site].end()) it_up == operator_string[adj_site].begin();
           oppdir = opposite_direction[ nb[ adj_site][j] ] [j]  ;
           if (!is_not_close(it_nb->time(), worm_head_it->time(), 2*tol ) ) break;
         }
         worm_head_it->set_assoc(j, it_up);
         it_nb->set_assoc(oppdir, worm_head_it);
         //find_assoc_delete(isite, worm_head_it);
         //find_assoc_insert(isite, worm_head_it, +1);
         worm_at_stop = +1;
       }
       // when measuring density matrix on other site : no problem
       else if (itn->color() == 0) { // dummy
         PASS_DUMMY(dir, isite ); // no further problem
         return accepted;
       }
       else if (itn->color() == 1) {
         return impossible; // cannot pass interaction
       }
       else if (itn->color() == 2)
      {
        return (2); // cannot pass interaction
      }
       else {
         cerr << "\n# Error in MOVEWORM : " << dir << "\t" << itn->color() << "\t" << itn->time() << "\t" << worm_head_it->time() << "\t" << worm_tail_it->time() << "\n";
         exit(1);
       }
     }
     
    
     get_diaginfo_nb_pos(isite);
     
     double t_loc = ( t_between(worm_tail_it->time(), worm_head_it->time(), itn->time() ) ) ? worm_tail_it->time() : itn->time() ; // for density matrix
     if (worm_meas_densmat) t_loc = itn->time(); // exceptional case
     nbtime[ZC] = t_loc;
     double dt_loc = t_loc - worm_head_it->time();
     if (dt_loc < 0) dt_loc += beta;
     dtime[ZC] = dt_loc;
     
     size_t nb_next = ZC; double dt_min = dtime[ZC];
     double tmax = nbtime[ZC];
     if (t_between(worm_tail_it->time() , worm_head_it->time(), tmax) ) {
       tmax = worm_tail_it->time();
       dt_min = worm_tail_it->time() - worm_head_it->time();
       if (dt_min < 0) dt_min += beta;
     }
     if (worm_meas_densmat) {
       tmax = nbtime[ZC];
       dt_min = dtime[ZC];
     }
     for (size_t i=zcoord[0][worm_head_it->link()] ; i < zcoord[1][worm_head_it->link()]; i++){
      if (dtime[i] < dt_min) {
       dt_min = dtime[i];
       nb_next = i;
       tmax = nbtime[i];
     }
     }

#ifdef DEBUGMODE 
     cout << "# MOVEWORM dir == +1 tmax " << tmax << "\tworm_meas_densmat " << ( worm_meas_densmat ? "true" : "false") << "\n";
     for (int i=0 ; i < ZC+1; i++) cout << "# MOVEWORM dir == +1 nb, dtime, nbtime " << i << "\t" << dtime[i] << "\t" << nbtime[i] << "\n";
     //char ch; cin >> ch;
#endif

     double tau = ( (tmax > worm_head_it->time() ) ? tmax - worm_head_it->time() : beta + tmax - worm_head_it->time() );
     StateType n_R = itn->before();
     StateType n_L = itp->after();
     std::vector<int> nbval_temp = nbval;
     if (!MyLatt->get_pbc(0) && (isite == 0)) nbval_temp.erase(nbval_temp.begin()+1);
     else if (!MyLatt->get_pbc(0) && (isite == Nsites-1)) nbval_temp.erase(nbval_temp.begin());
     double E_L = Diag_energy(isite, n_L, nbval_temp );
     double E_R = Diag_energy(isite, n_R, nbval_temp );
     double ratio, pexp, en_nom, en_denom;
     if (E_R > E_L) {
       //pexp = -log(random())/E_off;
       pexp = -log(rnd(MyGenerator))/E_off;
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ?  E_off : 1.;
       en_nom = (worm_at_stop == 0) ? E_R - E_L + E_off : 1.;
     }
     else {
       //pexp = -log(random()) / (E_L - E_R + E_off);
       pexp = -log(rnd(MyGenerator)) / (E_L - E_R + E_off);
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ? (E_L - E_R + E_off) : 1.;
       en_nom = (worm_at_stop == 0) ? E_off : 1.;
     }
     if (pexp < dtol) return impossible;
     ratio = en_nom / en_denom;
#ifdef DEBUGMODE
     cout << "MOVEWORM : worm head site-time : " << worm_head_it->link() << "\t" << worm_head_it->time()   << "\n";
     cout << "MOVEWORM : dir " << dir << "\t tau " << tau << "\t pexp : " << pexp << "\t n_R " << n_R << "\t n_L " << n_L << "\t E_R : " << E_R << "\t E_L : " << E_L << "\t ratio : " << ratio << "\n";
#endif
     if (Metropolis(ratio)) {
       if (pexp > tau || abs(pexp-tau) < dtol) {
         double t_w = tmax; //worm_head.time() + tau;
         //if (t_w > beta) t_w -= beta;
        
         //Element_t new_elem = *worm_head_it;
         //find_assoc_delete(isite, worm_head_it);
         //operator_string[isite].erase(worm_head_it);
         //new_elem.time(t_w);
         //worm_head_it = operator_string[isite].insert(itn, new_elem);
         //find_assoc_insert(isite, worm_head_it, -1);
         //worm_head.time(t_w);
         worm_head_it->time(t_w);
         worm_dtime += tau;
         //Epot_nb += tau * (n_L - n_R) * snb * V_nn  / beta;
         Epot_tot += tau * (E_L - E_R) / beta;
         Nprtcls += tau * (n_L - n_R);
         worm_at_stop = -1;
         worm_passes_nb_kink = (nb_next == ZC ? 0 : 1);
         //worm_dir = +1;
         worm_meas_densmat =  (worm_head_it->time() == worm_tail_it->time() ? true : false);
         return accepted;
       }
       else {
         double t_w = worm_head_it->time() + pexp;
         if (t_w > beta) t_w -= beta;
         //worm_head.time(t_w);
         //Element_t new_elem = *worm_head_it;
         //find_assoc_delete(isite, worm_head_it);
         //operator_string[isite].erase(worm_head_it);
         //new_elem.time(t_w);
         //worm_head_it = operator_string[isite].insert(itn, new_elem);
         //find_assoc_insert(isite, worm_head_it, 0);
         worm_head_it->time(t_w);
         worm_dtime += pexp;
         //Epot_nb += pexp * (n_L - n_R) * snb * V_nn  / beta;
         Epot_tot += pexp * (E_L - E_R) / beta;
         Nprtcls += pexp * (n_L - n_R);
         worm_at_stop = 0;
         worm_passes_nb_kink = 0;
         //worm_dir = +1;
         worm_meas_densmat = false;
         return accepted;
       }
     }
     else {
       return rejected;
     }
   } // ... dir == +1
   else { // dir == -1
     if (worm_at_stop == +1) {
       // several cases : dummy, worm_tail, density matrix, interaction
       if ( worm_meas_densmat ) {
         if (itp->color() < 0) return impossible; // cannot pass worm tail
         
         list<Element_t>::iterator it_up = worm_head_it;
         ++it_up; if (it_up == operator_string[isite].end()) it_up = operator_string[isite].begin();
         for (size_t j=zcoord[0][isite] ; j < zcoord[1][isite]; j++ ) {
           if (nb[isite][j] == worm_tail_it->link()) {
             SiteIndex adj_site =nb[isite][j];
             size_t oppdir = opposite_direction[adj_site ] [j]  ;
             //list<Element_t>::iterator it_nb = worm_head_it->get_assoc(j);
             //if (it_nb == operator_string[adj_site].begin()) it_nb = operator_string[adj_site].end(); --it_nb;
             worm_head_it->set_assoc(j, worm_tail_it);
             worm_tail_it->set_assoc(oppdir, it_up);
             break;
          }
        }
         
         
         
         //bool is_neighbor = false;
         //for (size_t j=zcoord[0][] ; j < zcoord[1][]; j++ ) if (nb[worm_head_it->link()][j] == worm_tail_it->link()) is_neighbor = true;
         //if (is_neighbor) {
         // find_assoc_delete(isite, worm_head_it);
         //  find_assoc_insert(isite, worm_head_it, -1);
         //}
         worm_at_stop = -1;
         // else do nothing
       }
       else if (worm_passes_nb_kink) {
         // do nothing ?
         list<Element_t>::iterator it_up = worm_head_it;
         ++it_up; if (it_up == operator_string[isite].end()) it_up = operator_string[isite].begin();
         list<Element_t>::iterator it_nb;
         size_t j, oppdir;
         for (j=zcoord[0][isite] ; j < zcoord[1][isite]; j++ ) {
           SiteIndex adj_site =nb[isite][j];
           it_nb = worm_head_it->get_assoc(j);
           if (it_nb == operator_string[adj_site].begin()) it_nb = operator_string[adj_site].end(); --it_nb;
           oppdir = opposite_direction [nb[isite][j] ] [j];
           if (!is_not_close(it_nb->time(), worm_head_it->time(), 2*tol ) ) break;
         }
         worm_head_it->set_assoc(j, it_nb);
         it_nb->set_assoc(oppdir, it_up );
         //find_assoc_delete(isite, worm_head_it);
         //find_assoc_insert(isite, worm_head_it, -1);
         worm_at_stop = -1;
       }
       // when measuring density matrix on other side : no problem
       else if (itp->color() == 0) { // dummy
         //if (worm_head.time() != beta) cerr << "# MOVEWORM, dir == -1, wrong worm time to pass dummy " << worm_head.time() << "\n";
         PASS_DUMMY(dir, isite ); // no further problem
         return accepted;
       }
       else if (itp->color() == 1) {
         return impossible; // cannot pass interaction          // XXXX WHY NOT PASSINTERACTION???
       }
       else if (itp->color() == 2)
      {
        return (2); // cannot pass interaction
      }
       else  {
         cerr << "\n# Error in MOVEWORM : " << dir << "\t" << itp->color() << "\t" << itp->time() << "\t" << worm_head_it->time() << "\t" << worm_tail_it->time() << "\n";
         exit(1);
       }
     }
    
     get_diaginfo_nb_neg(isite);
     
     double t_loc = ( t_between(worm_tail_it->time(), itp->time(), worm_head_it->time() ) ) ? worm_tail_it->time() : itp->time() ; // for density matrix
     if (worm_meas_densmat) t_loc = itp->time(); // exceptional case
     nbtime[ZC] = t_loc;
     double dt_loc =  worm_head_it->time() - t_loc;
     if (dt_loc < 0) dt_loc += beta;
     dtime[ZC] = dt_loc;
    
     double tmax = itp->time();
    
     size_t nb_next = ZC; double dt_min = dtime[ZC];
     if (t_between(worm_tail_it->time(), tmax, worm_head_it->time()) ) { // for density matrix
       tmax = worm_tail_it->time();
       dt_min = worm_head_it->time() - worm_tail_it->time();
       if (dt_min < 0) dt_min += beta;
     }
     if (worm_meas_densmat) {
       tmax = nbtime[ZC];
       dt_min = dtime[ZC];
     }
     for (size_t i=zcoord[0][worm_head_it->link()] ; i < zcoord[1][worm_head_it->link()]; i++){
       if (dtime[i] < dt_min) { //is this correct?
       dt_min = dtime[i];
       nb_next = i;
       tmax = nbtime[i];
      }
     } 
#ifdef DEBUGMODE
     cout << "# MOVEWORM dir == -1, tmax " << tmax << "\tworm_meas_densmat " << ( worm_meas_densmat ? "true" : "false") << "\n";
     for (size_t i=0 ; i < ZC+1; i++) cout << "# MOVEWORM dir == -1, nb, dtime, nbtime " << i << "\t" << dtime[i] << "\t" << nbtime[i] << "\n";
     //char ch; cin >> ch;
#endif
     double tau = ( (tmax < worm_head_it->time() ) ? worm_head_it->time() - tmax : beta + worm_head_it->time() - tmax );
     StateType n_R = worm_head_it->after();
     StateType n_L = worm_head_it->before();
     std::vector<int> nbval_temp = nbval;
     if (!MyLatt->get_pbc(0) && (isite == 0))  nbval_temp.erase(nbval_temp.begin()+1);
     else if (!MyLatt->get_pbc(0) && (isite == Nsites-1)) nbval_temp.erase(nbval_temp.begin());
     double E_L = Diag_energy(isite, n_L, nbval_temp );
     double E_R = Diag_energy(isite, n_R, nbval_temp );
     double pexp, ratio, en_nom, en_denom;
     if (E_R > E_L) {
       //pexp = -log(random()) / (E_R - E_L + E_off);
       pexp = -log(rnd(MyGenerator)) / (E_R - E_L + E_off);
       //ratio = E_off/ (E_R - E_L + E_off);
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ? E_R - E_L + E_off : 1.;
       en_nom = (worm_at_stop == 0) ? E_off : 1.;
     }
     else {
       //pexp = -log(random()) / (E_off);
       pexp = -log(rnd(MyGenerator)) / (E_off);
       en_denom = (tau > pexp || abs(pexp-tau) < dtol) ? E_off : 1.;
       en_nom = (worm_at_stop == 0) ? E_L - E_R + E_off : 1.;
       //ratio = (E_L - E_R + E_off) / E_off;
     }
     if (pexp < dtol) return impossible;
     ratio = en_nom / en_denom;
#ifdef DEBUGMODE
     cout << "MOVEWORM : worm head site-time : " << worm_head_it->link() << "\t" << worm_head_it->time()  << "\n";
     cout << "MOVEWORM : dir " << dir << "\t tau " << tau << "\t pexp : " << pexp << "\tn_R " << n_R << "\tn_L " << n_L << "\tE_R : " << E_R << "\t E_L : " << E_L << "\t ratio : " << ratio << "\n";
#endif
     if (Metropolis(ratio)) {
       if (pexp > tau || abs(pexp-tau) < dtol) {
         //double t_w = worm_head.time() - tau;
         //if (t_w < 0) t_w += beta;
         double t_w = tmax;
         if (tmax == beta) t_w = 0;
         
         //Element_t new_elem = *worm_head_it;
         //find_assoc_delete(isite, worm_head_it);
         //perator_string[isite].erase(worm_head_it);
         //new_elem.time(t_w);
         //worm_head_it = operator_string[isite].insert(itn, new_elem);
         //find_assoc_insert(isite, worm_head_it, +1);
         worm_head_it->time(t_w);
         //worm_head.time(t_w);
         worm_dtime -= tau;
         worm_at_stop = +1;
         //Epot_nb += snb * tau*V_nn* (n_R - n_L)/beta;
         Epot_tot += (E_R - E_L) * tau / beta;
         Nprtcls += tau * (n_R - n_L);
         worm_passes_nb_kink = (nb_next == ZC ? 0 : 1);
         //worm_dir = -1;
#ifdef DEBUGMODE
         cout << "MOVEWORM : worm head site-time : " << worm_head_it->link() << "\t" << worm_head_it->time()  << "\n";
         cout << "MOVEWORM : dir " << dir << "\t tau " << tau << "\t pexp : " << pexp << "\tn_R " << n_R << "\tn_L " << n_L << "\tE_R : " << E_R << "\t E_L : " << E_L << "\t ratio : " << ratio << "\n";
         cout << "MOVEWORM : " << worm_head_it->time() << "\t" << worm_tail_it->time() << "\n";
#endif
         worm_meas_densmat =  (worm_head_it->time() == worm_tail_it->time() ? true : false);
         return accepted;
       }
       else {
         double t_w = worm_head_it->time() - pexp;
         if (t_w < 0) t_w += beta;
         //worm_head.time(t_w);
         //Element_t new_elem = *worm_head_it;
         //find_assoc_delete(isite, worm_head_it);
         //operator_string[isite].erase(worm_head_it);
         //new_elem.time(t_w);
         //worm_head_it = operator_string[isite].insert(itn, new_elem);
         //find_assoc_insert(isite, worm_head_it, 0);
         worm_head_it->time(t_w);
         worm_dtime -= pexp;
         worm_at_stop = 0;
         //Epot_nb +=  snb * pexp*V_nn* (n_R - n_L)/beta;
         Epot_tot += pexp * (E_R - E_L) / beta;
         Nprtcls  += pexp * (n_R - n_L);
         worm_passes_nb_kink = 0;
         //worm_dir = -1;
         worm_meas_densmat = false;
         return accepted;
       }
     }
     else {
       return rejected;
    }
  } // ... dir == -1
}



int worm::INSERTKINK() {
 
#ifdef DEBUGMODE
  if (worm_diag) cerr << "Diagonal configuration in INSERTKINK ? \n\n";
#endif
  if (worm_at_stop != 0) return impossible;
  if (worm_meas_densmat) return impossible;
  if (worm_passes_nb_kink == 1) return impossible;
  
  SiteIndex isite = worm_head_it->link();

  // choose a direction
  //int dir = ( (random() < 0.5) ? -1 : 1 );
  int dir = ( (rnd(MyGenerator) < 0.5) ? -1 : 1 );
  
  // choose a neighbor
  //size_t nbs = static_cast<size_t> (random()*zc);
  size_t nbs = static_cast<size_t> (rnd(MyGenerator)*ZC);
  SiteType adj_site = nb[isite][ nbs];
  if (!MyLatt->get_pbc(0) && ((nbs == 1 && isite == 0) || (nbs == 0 && isite == Nsites-1))){
    return impossible;
  }
#ifdef DEBUGMODE
  cout << "INSERTKINK dir : "<< dir << "\tnbs " << nbs <<"\tadj_site " << adj_site << "\n";
#endif
  
  StateType n_up = worm_head_it->after();
  StateType n_down = worm_head_it->before();
  
  //int n_up = (worm_it->before() );
  //site_left(worm_head.to() );
  //int n_down = (site_it[worm_head.to()]->after());
  //site_right(worm_head.to());
  
  StateType n_diff = n_up - n_down;
  
  // let's set the iterator on the adjacent site right
  list<Element_t>::iterator it_adj = worm_head_it->get_assoc(nbs);
  list<Element_t>::iterator itp_adj = it_adj;
  //site_it[adj_site] = site_it[worm_head.to()]->get_assoc(nbs);
  //list<Element_t>::iterator previt = site_it[adj_site];

  if (itp_adj == operator_string[adj_site].begin()) itp_adj = operator_string[adj_site].end();
  --itp_adj;
    while (!t_between(worm_head_it->time(), itp_adj->time(), it_adj->time() ) ) {
    it_adj = itp_adj;
    if (itp_adj == operator_string[adj_site].begin()) itp_adj = operator_string[adj_site].end();
    --itp_adj;
    }
  StateType nbs_out = it_adj->before();
  StateType nbs_mid = nbs_out + (rnd(MyGenerator) < 0.5 ? -1 : 1) * n_diff; // n_diff is the change in occupancy across the worm head

#ifdef DEBUGMODE
  cout << "INSERTKINK : n_up" << n_up << "\tn_down" << n_down << "\tnbs_out " << nbs_out << "\tnbs_mid " << nbs_mid << "\n";
#endif
  
#ifdef UNIFORM
  if (!range_check(nbs_mid)) return impossible;
#else
  if (MyModel->range_fail(nbs_mid)) return impossible;
#endif
    
    
  //int weight1 = max(nbs_mid, nbs_out) ;
  
  BondIndex b = bond_index[isite][ nbs];
    
#ifdef UNIFORM
  double weight = (dir == +1 ? bond_weight_offdiag[n_down][nbs_out][n_up][nbs_mid] : bond_weight_offdiag[n_down][nbs_mid][n_up][nbs_out]);
#else
  double weight_old = MyModel->site_weight_offdiag(isite, n_down, n_up);
  double weight_new = (dir == 1 ? MyModel->bond_weight_offdiag(b, n_down, n_up, nbs_out, nbs_mid) : MyModel->bond_weight_offdiag(b, n_down, n_up, nbs_mid, nbs_out) );
  weight_new *= MyModel->site_weight_offdiag(adj_site, nbs_mid, nbs_out);
#endif
int color = (dir == +1 ? (nbs_mid == n_up ? 2 : 1) : (nbs_mid == n_down ? 2 : 1)); 
  //double ratio = 2*zcoord * weight1 * t_hop * update_prob[deletekink] / (update_prob[insertkink]);

  double ratio = 4 * ZC * abs(weight) * update_prob[deletekink] / update_prob[insertkink]; //test
  /*
  if (is_not_close(weight_new/weight_old , max(nbs_mid, nbs_out) * 1., 1e-10 ) ) {
    cerr << "INSERTKINK wrong weights \n";
    exit(1);
  }
  */ //test
#ifdef DEBUGMODE
  cout << "INSERTKINK : ratio " << ratio << "\t" << weight_new << "\t" << weight_old << "\n";
  
  //char ch; cin >> ch;
#endif
  
  
  if (Metropolis(ratio)) {
#ifdef DEBUGMODE
    cout << "INSERTKINK accepted. worm_head_it " << worm_head_it->time() << "\tratio : " << ratio << "\n";
    cout << "INSERTKINK adj_site : " << adj_site << "\tit_adj " << it_adj->time() << "\n";
#endif
   
    // insert the interaction : iw on the worm site, change the properties of the worm_head_it and insert a new kink on the adj site + a new worm there
    Element_t new_elem1(nbs_mid, nbs_out, isite,     worm_head_it->time(),  color);
    Element_t new_elem2(nbs_out, nbs_mid, adj_site,  worm_head_it->time(), worm_tail_it->color());
    worm_head_it->color(color);
    worm_head_it->link(adj_site);
    
    if (dir ==1 ) {
      new_elem1.before(nbs_out);
      new_elem1.after(nbs_mid);
      new_elem2.before(nbs_mid);
      new_elem2.after(nbs_out);
      worm_head_it = operator_string[adj_site].insert(it_adj, new_elem2);
      find_assoc_insert(adj_site, worm_head_it, +1);
      it_adj = operator_string[adj_site].insert(worm_head_it, new_elem1);
      find_assoc_insert(adj_site, it_adj, 0);
      worm_at_stop = +1;
    }
    else {
      it_adj = operator_string[adj_site].insert(it_adj, new_elem1);
      find_assoc_insert(adj_site, it_adj, 0);
      list<Element_t>::iterator new_it = operator_string[adj_site].insert(it_adj, new_elem2);
      find_assoc_insert(adj_site, new_it, -1);
      worm_head_it->set_assoc(nbs, it_adj);
      worm_head_it = new_it;
      worm_at_stop = -1;
    }
        
    // worm jumps
    if (color == 1)
    {
      nrvertex_type1++;
    }
    else
    {
      nrvertex_type2++;
    }
        int mi = min(isite,adj_site);
    int ma = max(isite,adj_site);
    if (ma-mi>1)
    {
       nrvertex[ma]++;
    }
    else
    {
       nrvertex[mi]++;
    }
    worm_passes_nb_kink = 0;
    //worm_dir = dir;
    worm_meas_densmat = false;
    //update_spinspin_corr();
    return accepted;
  }
  else {
    return rejected;
  }
}


int worm::DELETEKINK() {
  if (worm_meas_densmat ) return impossible; // no density matrix measurement
  
#ifdef DEBUGMODE
  cout << "# Welcome to DELETEKINK : \n";
#endif

  list<Element_t>::iterator it = worm_head_it;                                    // it should become the iterator pointing to the kink to be deleted; needs to be changed depending on worm_at_stop == +/- 1
  SiteIndex isite = it->link();                                                   // worm is linked to itself
  StateType n_out, n_mid, n_new;                                                         // n_new is the occupation which will be at the position of the kink if it's deletion is possible 
  int color;
  if (worm_at_stop == +1 && worm_passes_nb_kink == 0) {                         // ie, the worm is just above the kink
    n_new = it->after();
    if (it == operator_string[isite].begin()){
      it = operator_string[isite].end();
    }
    --it;
    color = it->color();
    n_out = it->before();
    n_mid = it->after();
    if (color != +1 && color != +2) return impossible; // check if interaction, else return (worm tail)
    if (n_new != n_out) return (PASSINTERACTION(-1, it));      // this is a situation where the new occupation is different from the occupation after the kink -> kink can not be deleted but must be passed 
  } 
  else if (worm_at_stop == -1 && worm_passes_nb_kink == 0) { // dir==+1           // ie, the worm is just below the kink
    n_new = it->before();
    ++it; if (it == operator_string[isite].end()) it = operator_string[isite].begin();
    color = it->color();
    n_out = it->after();
    n_mid = it->before();
    if (color != +1 && color != +2) return impossible; // check if interaction, else return (worm tail)
    if (n_new != n_out) return (PASSINTERACTION(+1, it));       // this is a situation where the new occupation is different from the occupation after the kink -> kink can not be deleted but must be passed 
  }
  else {
    return impossible;
  }
  
  SiteIndex adj_site = it->link();                                    // site to wich the current site is linked
  std::size_t linkdir = 0;
    for (linkdir =0 ; linkdir < ZC; linkdir++) {
      if (nb[isite][linkdir] ==adj_site) break;
    }
  if (Nsites == 2)
  {
    linkdir = isite;
  }
  BondIndex b = bond_index[isite][ linkdir];               // bond index of the kink
  list<Element_t>::iterator itlink = it->get_assoc(linkdir);
  StateType n_A_link = itlink->after();                               // occupancy on kink after  the hopping event on the linked site
  StateType n_B_link = itlink->before();                              // occupancy on kink before the hopping event on the linked site
  
  // weight of tunneling element
#ifdef UNIFORM
  double weight = (worm_at_stop == -1 ? bond_weight_offdiag[n_B_link][n_mid][n_A_link][n_out] : bond_weight_offdiag[n_B_link][n_out][n_A_link][n_mid]);
#else
  double weight_old = MyModel->bond_weight_offdiag(b, it->before(), it->after(), n_B_link, n_A_link) * MyModel->site_weight_offdiag(isite, worm_head_it->before(), worm_head_it->after() );
  double weight_new = MyModel->site_weight_offdiag(adj_site, n_B_link, n_A_link);
#endif
  //double ratio = weight_new / weight_old * update_prob[insertkink] / (update_prob[deletekink] * 2* ZC );
  double ratio = 4 * ZC * abs(weight) * update_prob[deletekink] / update_prob[insertkink]; //todo: is this correct?
/*
  if (is_not_close(weight_old/weight_new , max(n_mid, n_out) * 1., 1e-10 ) ) {
    cerr << "DELETEKINK wrong weights " << weight_old << "\t" << weight_new << "\t" << n_out << "\t" << n_mid << "\n";
    exit(1);
  }
  
#ifdef DEBUGMODE
   cout << "DELETEKINK : n_out " << n_out << "\tn_mid " << n_mid  << "\n";
   cout << "DELETEKINK it->link() " << it->link() << "\tit->time()" << it->time() << "\n";
#endif
*/ 
#ifdef DEBUGMODE
  
   cout << "DELETEKINK : ratio " << ratio << "\t" << 1./ratio << "\n";
  cout << "DELETEKINK old weight, new weight : " << "\t" << weight_old << "\t" << weight_new << "\n";
  //if (is_not_close(weight_old/weight_new , max(nbs_mid, nbs_out) * 1., 1e-10 ) ) {
 //   cerr << "DELETEKINK wrong weights \n";
  //  exit(1);
  //}
   //char ch; cin >> ch;
#endif

  if (Metropolis(1/ratio)) {
    double t_w = it->time();
    int color = it->color();
    if (worm_at_stop == +1) {
      find_assoc_delete(isite, worm_head_it);
      operator_string[isite].erase(worm_head_it);
      find_assoc_delete(isite, it);
      operator_string[isite].erase(it);
    }
    else {
      find_assoc_delete(isite, it);
      operator_string[isite].erase(it);
      find_assoc_delete(isite, worm_head_it);
      operator_string[isite].erase(worm_head_it);
    }
    worm_head_it = itlink;
    worm_head_it->color(worm_tail_it->color());
    worm_head_it->link(adj_site);
    if (color == 1)
    {
      nrvertex_type1--;
    }
    else
    {
      nrvertex_type2--;
    }
       int mi = min(isite,adj_site);
    int ma = max(isite,adj_site);
    if (ma-mi>1)
    {
       nrvertex[ma]--;
    }
    else
    {
       nrvertex[mi]--;
    }    
    worm_at_stop = 0;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
    //update_spinspin_corr();
    return accepted;
  }
  else { // Metropolis rejection
    return rejected;
  }
  
}


void worm::PASS_DUMMY(const int dir, const SiteType isite) {
#ifdef DEBUGMODE
  cout << "\n# PASSING DUMMY... dir = " << dir << "\n";
#endif
  if (dir == +1) {
    //list<Element_t>::iterator it = dummy_it[isite];
    //if (it == operator_string[isite].begin()) it = operator_string[isite].end();
    //--it;
    StateType new_dens = worm_head_it->before();
    dummy_it[isite]->before(new_dens);
    dummy_it[isite]->after(new_dens);
    
    Element_t new_elem = *worm_head_it;
    find_assoc_delete(isite, worm_head_it);
    operator_string[isite].erase(worm_head_it);
    new_elem.time(0);
    worm_head_it = operator_string[isite].insert(operator_string[isite].begin(), new_elem);
    find_assoc_insert(isite, worm_head_it, 0);
    //site_right(isite);
    //worm_it = site_it[isite];
    worm_at_stop = +1;
    worm_passes_nb_kink = 0;
    //worm_dir = +1;
    worm_meas_densmat = false;
  }
  else {
    int new_dens = worm_head_it->after();
    dummy_it[isite]->before(new_dens);
    dummy_it[isite]->after(new_dens);
    //worm_it = dummy_it[isite];
    //site_it[isite] = worm_it;
    //site_left(isite);
    Element_t new_elem = *worm_head_it;
    find_assoc_delete(isite, worm_head_it);
    operator_string[isite].erase(worm_head_it);
    new_elem.time(beta);
    worm_head_it = operator_string[isite].insert(dummy_it[isite], new_elem);
    find_assoc_insert(isite, worm_head_it, -1);
    worm_at_stop = -1;
    worm_passes_nb_kink = 0;
    //worm_dir = -1;
    worm_meas_densmat = false;
  }
#ifdef DEBUGMODE
  print_conf(std::cout);
  std::cout << "# ... done with PASSDUMMY\n";
#endif
  
  return;
}

int worm::PASSINTERACTION(const int dir, const list<Element_t>::iterator it) {
  // here : it is the iterator to the kink that will be passed, not the worm_head_it
  SiteIndex isite = worm_head_it->link(); // current site
  StateType n_new = (dir == +1 ? worm_head_it->before() : worm_head_it->after()); // new occupation at kink
  StateType n_A = it->after();                                        // occupancy on kink after the hopping event on the current site  (before and after as always in the view of positive time [0, beta[)
  StateType n_B = it->before();                                       // occupancy on kink before the hopping event on the current site
  SiteIndex adj_site = it->link();                                    // site to wich the current site is linked
  std::size_t linkdir = 0;
    for (linkdir =0 ; linkdir < ZC; linkdir++) {
      if (nb[isite][linkdir] ==adj_site) break;
    }
  if (Nsites == 2)
  {
    linkdir = isite;
  }
  if (linkdir > 1)
  {
    cout << "linkdir: " << linkdir << endl; 
  }
  
  BondIndex b = bond_index[isite][ linkdir];               // bond index of the kink
  StateType n_A_link = it->get_assoc(linkdir)->after();               // occupancy on kink after  the hopping event on the linked site
  StateType n_B_link = it->get_assoc(linkdir)->before();              // occupancy on kink before the hopping event on the linked site
  
  Element_t new_elem = *worm_head_it;
  find_assoc_delete(isite, worm_head_it);
  operator_string[isite].erase(worm_head_it);
  new_elem.time(it->time());
  new_elem.link(adj_site);

  list<Element_t>::iterator itn = it->get_assoc(linkdir);
  
  

  if (n_A == n_A_link || n_B == n_B_link) // type 2 interaction
  {
    #ifdef DEBUGMODE
    if (it->color() != 2)
    {
      std::count << "# PASSINTERACTION color is not consistent with occupation numbers?\n";
    }
    #endif
    if (dir == +1)
    {
      it->before(n_new);
      it->get_assoc(linkdir)->before(n_new);
    }
    else
    {
      it->after(n_new);
      it->get_assoc(linkdir)->after(n_new);
      ++itn; if (itn == operator_string[isite].end() ) itn = operator_string[isite].begin();
    }
    StateType n_temp = new_elem.before(); // for type 2 interaction we have a change in direction -> occupatione before and after the worm head need to be switched
    new_elem.before(new_elem.after());
    new_elem.after(n_temp);
  }
  else if (n_A == n_B_link || n_B == n_A_link) // type 1 interaction
  {
    #ifdef DEBUGMODE
    if (it->color() != 1)
    {
      std::count << "# PASSINTERACTION color is not consistent with occupation numbers?\n";
    }
    #endif
    if (dir == +1)
    {
      it->before(n_new);
      it->get_assoc(linkdir)->after(n_new);
      ++itn; if (itn == operator_string[isite].end() ) itn = operator_string[isite].begin();
    }
    else
    {
      it->after(n_new);
      it->get_assoc(linkdir)->before(n_new);
    }
    worm_at_stop = -worm_at_stop;
  }
  else
  {
    std::cerr << "# Error in PASSINTERACTION!\n";
  }
  worm_head_it = operator_string[adj_site].insert(itn, new_elem);
  find_assoc_insert(adj_site, worm_head_it, worm_at_stop);
  worm_meas_densmat = false;
  return accepted;
}



int worm::GLUEWORM() {
#ifdef DEBUGMODE
  cout << "# Welcome to GLUEWORM  : worm_at_stop : " << worm_at_stop << "\n";
  if (worm_diag) cerr << "Diagonal configuration in GLUEWORM ? \n\n";
#endif
  size_t isite = worm_head_it->link();
  if (worm_at_stop == 0) return impossible;
  if (worm_passes_nb_kink) return impossible;
  if (worm_tail_it->time() != worm_head_it->time()) return impossible;
  list<Element_t>::iterator it_lower = worm_head_it;
  list<Element_t>::iterator it_upper = worm_head_it;
  if (worm_at_stop == -1) {
    ++it_upper; if (it_upper == operator_string[isite].end()) it_upper = operator_string[isite].begin();
    if (!(it_upper == worm_tail_it)) return impossible;
  }
  else {
    if (it_lower == operator_string[isite].begin()) it_lower = operator_string[isite].end();
    --it_lower;
    if (!(it_lower == worm_tail_it)) return impossible;
  }
  StateType n_out = it_upper->after();
  if (n_out != it_lower->before()) {
    std::cerr << "Really strange n_out in GLUEWORM : " << it_upper->after() << "\t" << it_lower->before()<< "\n";
    return impossible;
  }
  StateType n_mid = it_upper->before();
  if (n_mid != it_lower->after()) {
    std::cerr << "Really strange n_mid in GLUEWORM : " << it_upper->before() << "\t" << it_lower->after() << "\n";
    return impossible;
  }
#ifdef UNIFORM
  double wgt = site_weight_offdiag[n_out][ n_mid];
  double wgt_sq = wgt * wgt;
#else
  double wgt = MyModel->site_weight_offdiag(isite, n_out, n_mid);
  double wgt_sq = wgt * wgt;
#endif
  
  //double ratio = 4 * C_NBW * Nsites * beta * wgt_sq*  update_prob[glueworm] / update_prob[insertworm];
  double ratio = 4 * C_NBW * Nsites * beta *  update_prob[glueworm] / update_prob[insertworm];
  
  
#ifdef DEBUGMODE
  cout << "GLUEWORM  ratio " << ratio << "\t" << 1./ratio << "\n";
  //char ch; cin >> ch;
#endif
  if (Metropolis(1./ratio)) {
    //update_spinspin_corr();
#ifdef DEBUGMODE
    cout << "GLUEWORM accepted. ratio " << ratio << "\n";
#endif
    find_assoc_delete(isite, it_upper);
    operator_string[isite].erase(it_upper);
    find_assoc_delete(isite, it_lower);
    operator_string[isite].erase(it_lower);
    
    worm_head_it = dummy_it[isite];
    worm_tail_it = dummy_it[isite];
    worm_diag = 1;
    worm_at_stop = 0;
    worm_passes_nb_kink = 0;
    worm_meas_densmat = false;
    return accepted;
  } // Metropolis
  else { // Metropolis
#ifdef DEBUGMODE
    cout << "GLUEWORM rejected. ratio " << ratio << "\n";
#endif
    return rejected;
  } // else...Metropolis
}