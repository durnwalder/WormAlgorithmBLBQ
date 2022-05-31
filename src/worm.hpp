/* $Id: worm.hpp,v 1.1 2006/09/09 8:59:59 pollet Exp $ */

#ifndef worm_HPP
#define worm_HPP

#pragma once

#include "lattice.hpp"

#include "worm.Element.hpp"
#include "blbq_model.hpp"
#include <alps/mc/mcbase.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/multi_array.hpp>
#include <alps/accumulators.hpp>
#include <alps/mc/api.hpp>
#include <alps/mc/stop_callback.hpp>

#include <random>
#include <vector>
#include <valarray>
#include <list>
#include <numeric>
#include <time.h>
#include <iomanip>
#include <exception>
#include <cstdlib>
#include <fstream>
#include <functional> 


typedef boost::multi_array<double, 2> DMatrix;
typedef boost::multi_array<double, 2>::index DMatrix_index;

typedef boost::multi_array<double, 4> DTensor;
typedef boost::multi_array<double, 4>::index DTensor_index;


typedef boost::multi_array<int, 2> IMatrix;
typedef boost::multi_array<int, 2>::index IMatrix_index;



typedef Element<SiteIndex, StateType, ZC> Element_t;

enum statistics_tag {
    impossible,
    rejected,
    accepted,
    total_attempted,
    statistics_count
};
namespace {
    const std::string statistics_names[] = {
        "impossible",
        "rejected",
        "accepted",
        "total_attempted"
    };
    const std::string update_names[] = {
        "insertworm",
        "moveworm",
        "insertkink",
        "deletekink",
        "glueworm",
    };
}


using namespace std;


class worm : public alps::mcbase {
public :
  
  //worm(parameters_type const & parms, std::size_t seed_offset = 0);
  worm(parameters_type const & parms, std::size_t seed_offset = 0);
    
  static std::string code_name();
    
  // ALPS overloads:
  static void define_parameters(parameters_type & parameters);

  virtual void update();
  virtual void measure();
  virtual double fraction_completed() const;

  using alps::mcbase::save;
  virtual void save(alps::hdf5::archive & ar) const;

  using alps::mcbase::load;
  virtual void load(alps::hdf5::archive & ar);
  
  bool is_not_close(const double val1, const double val2, const double tol=1e-10) const {
    bool answ = false;
    if (fabs(val1) > 1.) {
      if (fabs(val2) < 0.001) return true;
      if (fabs(1. - val1/val2) > tol) return true;
    }
    else {
      if (fabs(val1 - val2) > tol) return true;
    }
    return answ;
  }
  
  
  bool Metropolis(const double& x) {
 //#ifdef DEBUGMODE
     if (x < 0 || std::isnan(x) ) {
       cerr << "x < 0 in Metropolis? " << x << "\n";
       char ch; cin >> ch;
     }
 //#endif
     if (x >= 1) return (true);
     if (rnd(MyGenerator) < x) return (true);
     //if (random() < x) return (true);
     return (false);
   }

  double Diag_energy(const SiteIndex s, const StateType n, const size_t dir, const size_t m) {
#ifdef UNIFORM
    double d = site_weight_diag[n-nmin];
    d  += bond_weight_diag[n-nmin][m-nmin];
    return d;
#else
    double d = MyModel->site_weight_diag(s,n);
    d += MyModel->bond_weight_diag( bond_index[s][dir],n, m);
    return d;
#endif
  }

  double Diag_energy(const SiteIndex s, const StateType n, const vector<StateType>& snb) {
#ifdef UNIFORM
    double d = site_weight_diag[n - nmin];
    for (size_t m=0 ; m < snb.size(); m++ ) d  += bond_weight_diag[n-nmin][snb[m]-nmin];
    return d;
#else
    double d = MyModel->site_weight_diag(s,n);
    for (size_t i=0 ; i < snb.size(); i++) d+= MyModel->bond_weight_diag( bond_index[s][i],n, snb[i]);
    return d;
#endif
  }
  
  
  /*
  double MCstep_total, MCstep_run;
  double MCdiag;
  long Nsave;
  long MCstep;
  int32_t MCmeas;
  double acc_insert, tot_insert;
  typedef Element element_type;
  worm() : number_of_updates(5), number_statistics(glueworm+1), MyGenerator(rd()) {};
  static void print_copyright(std::ostream&);
    
  
  int dostep();
  void read_params();
  
  
  //double dist_realsq(const SiteType s1);
  int get_dim() {return dim;}
  SiteType get_Nsites() {return Nsites;}
  double get_kinetic_energy() {if (new_measurement) update_en(); return (Ekin);} // -nrvertex/beta);}
  double get_potential_energy() {if (new_measurement) update_en(); return (Epot);} //{return (calc_potential_energy());}
  double get_energy() {return (Epot + Ekin);} //{return (calc_potential_energy() - nrvertex/beta);}
  long get_Ncan() {return (Ncan);}
  long get_Npart() { if (new_measurement) update_en();  return (number_of_particles); }// {double d; return(calc_number_of_particles(d));} 
  long get_Npart(double& d) { if (new_measurement) update_en(); return (number_of_particles); } //{return(calc_number_of_particles(d));}
  double get_worm_insert_ratio() {return (acc_insert / tot_insert);}
  
  
  void calc_winding();
  //long calc_number_of_particles(double& );
  double calc_local_energy(const SiteType) const;
  void geometry();
  void update_system(const SiteType, int&);
  void worm_cycle(double&);
  bool worm_diag, worm_meas_densmat;
  int worm_at_stop;
  int worm_passes_nb_kink;
  double worm_dtime;  

  
    
  
  void save(ostream&, vector <int32_t >& );
  void load(istream&, vector <int32_t >& );
 
    
  //random number generator    

  

  //uint32_t Ntherm;
  //uint32_t Nloop;
  uint32_t itherm, iloop;

  void update_av() {
    mZ_state += 1.;
    for (SiteType i = 0; i < Nsites; i++) av_state[i] += state[i];
    for (SiteType i = 0; i < Nsites; i++) av_state_sq[i] += state[i]*state[i];
    if (new_measurement) update_en();
  }

  void update_en() {
    Ekin = -nrvertex / beta;
    Epot = calc_potential_energy();
    new_measurement = false;
  }

  

  
  double time_to_next_interaction(const int isite, const double t0, const int dir, double& time_to_dummy) {
	if (dir == +1) {
	  time_to_dummy = ( (dummy_it[isite]->time() > t0) ? dummy_it[isite]->time() - t0 : beta - t0 + dummy_it[isite]->time() );
	  if (operator_string[isite].size() == 1) return (beta);
	  // we assume that the site_iterator already points to the next interaction
	  double tau = ( (site_it[isite]->time() > t0) ? site_it[isite]->time() - t0 : beta + site_it[isite]->time() - t0);
	  if (site_it[isite] == dummy_it[isite]) {
		double t1 = site_it[isite]->time();
		list<element_type>::iterator it = site_it[isite];
		++it;
		if (it == operator_string[isite].end()) it = operator_string[isite].begin();
		tau += ( (it->time() > t1) ? it->time() - t1 : beta - t1 + it->time() );
	  }
	  return (tau);
	} else {
	  time_to_dummy = ( (dummy_it[isite]->time() < t0) ? t0 - dummy_it[isite]->time() : beta + t0 - dummy_it[isite]->time() );
	  if (operator_string[isite].size() == 1) return (beta);
	  // we do not assume that the site iterator already points to the previous interaction
	  site_left(isite);
	  double tau = ( (site_it[isite]->time() < t0) ? t0 - site_it[isite]->time() : beta + t0 - site_it[isite]->time() );
	  // cout << "# time_to_next interaction (left)  " << t0 << "\t" << site_it[isite]->time() << "\t" << tau;
	  if (site_it[isite] == dummy_it[isite] ) {
		double t1 = site_it[isite]->time();
		list<element_type>::iterator it = site_it[isite];
		if (it == operator_string[isite].begin() ) it = operator_string[isite].end() ;
		--it;
		tau += ( (it->time() < t1) ? t1 - it->time() : beta + t1 -it->time() );		
	    // cout << "# time_to_next interaction (left) with dummy " << t1 << "\t" << it->time() << "\t" << tau;
	  }
	  // we put the iterator back
	  site_right(isite);
	  return (tau);
	}
  }
	
        

  void reset_av() {
    mZ = 0.;
    mG = 0.;
	mZ_dns = 0.;
    mZ_state = 0.;
    for (int32_t i = 0; i < Nsites; i++) {
      av_dns[i] = 0.;
      av_state[i] = 0.;
      av_state_sq[i] = 0.;
    }
	mZ_green = 0.;
	for (int32_t i = 0; i < Ntimes; i++) av_green[i] = 0.;
    tot_insert = 0.;
    acc_insert = 0.;
    MCdiag = 0;
	MCstep = 0;
    MCstep_run = 0.;
    MCmeas = 0;
    times_run = 0;
  }


  void print_av(std::ostream& os)   {
    double s = 0.;
	for (int i = 0; i < Nsites; i++) s += av_state[i]/mZ_state;
	
	double norm = Nsites * av_dns[0] / s ;

	if (dim == 3) {
      for (int k = 0; k < Ls[2]; k++) {
	    for (int j = 0; j < Ls[1]; j++) {
	      for (int i = 0; i < Ls[0]; i++) {
	        int32_t s = i + j*Ls[0] + k * Ls[0]*Ls[1];
	        os << i << "\t" << j << "\t" << k << "\t" << av_state[s]/mZ_state << "\t" << av_state_sq[s]/mZ_state << "\t" << av_dns[s]/mZ_dns << "\n" ;
	      }
	    }
      }
    } // dim == 3
    else if (dim == 1) {
	  for (int i = 0; i < Ls[0]; i++) {
	    os << i << "\t" << av_state[i]/mZ_state << "\t" << av_state_sq[i]/mZ_state << "\t" << av_dns[i]/mZ_dns << "\t" << av_dns[i]/norm << "\n";
	  }
	  os << "# norm : " << norm << "\t mZ_dns : " << mZ_dns << "\n";
	  os << "# mZ : " << mZ << "\t MG : " << mG << "\tratio : " << mZ/mG << "\n";
	} // dim == 1
	else {
	 os << "# print_av not implemented for this dimension.\n";
	}
  }

  void print_green(std::ostream& os) {
    double norm = av_green[Ntimes/2];
    for (int i = 0; i < Ntimes+1; i++) os << time_green[i] << "\t" << av_green[i]/norm << "\n";
	return;
  }

  void print_av_chem(std::ostream& os) {
    if (dim == 3) {
      //int32_t Npoints = min( min (Ls[0], Ls[1]), Ls[2]);
      //Npoints = (Npoints % 2 == 0 ? Npoints / 2 : (Npoints + 1) / 2);
      int32_t vol = 1;
      for (int k = 0; k < dim; k++) vol += Ls[k]*Ls[k]/4;
      int32_t Npoints = int(sqrt(vol*1.))+1;
      vector<int32_t> Nentry(Npoints);
      vector<double> prof(Npoints);
      for (int32_t i = 0; i < Npoints; i++) {
	Nentry[i] = 0;
	prof[i] = 0.;
      }
 
      double mz = Ls[2]/2.0 - 0.5;
      double my = Ls[1]/2.0 - 0.5;
      double mx = Ls[0]/2.0 - 0.5;
      for (SiteType k = 0; k < Ls[2]; k++) {
	SiteType z = k*Ls[0]*Ls[1];
	for (SiteType j =0; j < Ls[1]; j++) {
	  SiteType y = j*Ls[0];
	  for (SiteType i =0; i < Ls[0]; i++) {
	    SiteType p = i + y + z;
      	    double d = vc[2]*(k-mz)*(k-mz) + vc[1]*(j-my)*(j-my) + vc[0]*(i-mx)*(i-mx);
	    d = sqrt(d/vc[1]);
	    int32_t pos = int(d);
	    if (pos >= Npoints) {
	      os << "Increase Npoints " << pos << "\t" << Npoints << "\n";
	      pos = Npoints-1;
	    }
	    Nentry[pos]++;
	    prof[pos] += av_state[p];
	  }
	}
      }
	    
      for (int32_t i = 0; i < Npoints; i++) {
	if (Nentry[i] > 0) os << i <<  "\t" << prof[i]/mZ_state/Nentry[i] << "\t" << Nentry[i] << "\n";
      }

    }
    else {
      os << "print_av_chem only implemented for dim == 3.\n";
    }
  }


  void get_density(valarray<double>& m) {for (SiteType i = 0; i < Nsites; i++) m[i] = state[i];}
  void get_winding(valarray<double>& m)  {calc_winding(); for (int i =0; i < dim; i++) m[i] = mWinding[i];}
  double get_rho_sf() {
    // only call this function after having called get_winding or calc_winding;
    double r = 0.;
    for (int i =0; i < dim; i++) {
      r+= mWinding[i] * mWinding[i];
    }
    if (dim == 1) {
      r *= Ls[0] / beta;
    }
    else if (dim == 2) {
      r /= (2. * beta);
    }
    else if (dim == 3) {
      r /= (3. * beta * Ls[0]);
    }
    return (r);
  }
  
  time_t times_runtime, times_tot, times1, times2, dtimes, times_overall, times_run, times_log1, times_log2, timesm1, timesm2;
  int read_configuration, read_statistics;
  int32_t Ntest, Nmeasure, Nmeasure_green;

  void reset_update_statistics();
   
   
  */
  
  double calc_potential_energy_measure();
  double calc_dop();
  double calc_potential_energy_nb();
  double calc_potential_energy_loc();
  
  int INSERTWORM();
  int MOVEWORM();
  int INSERTKINK();
  int DELETEKINK();
  int GLUEWORM();
  int PASSINTERACTION(const int, list<Element_t>::iterator);
  void PASS_DUMMY(const int, const SiteType);
  int UPDATE_DENSITYMATRIX();
  
  void output_final(const alps::results_type<worm>::type& results, alps::hdf5::archive& ar);
  void initialize();
  void test_conf();
  void print_update_statistics(std::ostream&) const;
  void print_params(std::ostream& os) const;
  void print_conf(ostream&);
  void print_conf();
  
  void measure_density_matrix();
  void update_hist();
  void update_spinspin_corr();
  
protected:
  unsigned long sweeps;
  unsigned long thermalization_sweeps;
  unsigned long total_sweeps;
  size_t Nloop;
  size_t runtimelimit;
  
  int MCdiag;
  
  enum update_tag {insertworm, moveworm, insertkink, deletekink, glueworm, update_count};
  double update_prob_cuml[update_count];
  double update_prob[update_count];
private:
  void initialize_update_prob(update_tag begin, update_tag end);
  void initialize_update_prob();

    
private :
  static constexpr double dtol = 1e-15;
  // system parameters
  //size_t Lx, Ly, Lz;      // linear system sizes
  
  //double t_hop;           // hopping amplitude
  //double U_on;            // on-site dens-dens repulsion : U/2 n_i(n_i-1)
  //double V_nn;            // nearest neighbor dens-dens interaction V_nn n_i n_j
  //double mu;              // chemical potential
  
  std::unique_ptr<Lattice_t> MyLatt;
  std::unique_ptr<blbq_model> MyModel;
  
  double beta;            // inverse temperature
  //int nmax;               // maximum occupation number
  //int Ntimes;
  double E_off;           // energy offset, internal parameter of the code
  long int seed;          // seed of random number generator
  double C_worm, C_NBW;
  int canonical;          // canonical number of particles
  
  unsigned long Ntest;
  unsigned long Nsave;
  size_t Nmeasure;
  
  
  //vector<int> Ls;    // lattice dimensions
  //vector<double> mid;
  //SiteType Nsites; // total number of sites
  
  
  // internal variables
  std::vector<double> nrvertex;
  std::vector<double> nrv_diff;
  unsigned long long nrvertex_type1;
  unsigned long long nrvertex_type2;
  long long number_of_particles;
  double Ekin, Epot, Epot_tot, En, Nprtcls;
  bool passed_dummy;
  std::vector<int> state; // occupation per site
  std::vector<std::list<Element_t> > operator_string;
  
 
  
  //vector<double> dtime;
  //vector<double> nbtime;
  //std::vector<double> trans; // transition probabilies
  
  
  // system variables
  
  //boost::multi_array<size_t, 2> nb;
  
 
  // worm variables
  //list<Element>::iterator worm_it;  // iterator pointing at the next element the worm will encounter
  //vector<list<Element>::iterator > site_it; // iterator for site
  std::list<Element_t>::iterator worm_tail_it;  // iterator pointing to the tail
  std::list<Element_t>::iterator worm_head_it;  // iterator pointing to the tail
  std::vector<std::list<Element_t>::iterator > dummy_it; // iterator pointing at dummy elements useful for measuring diag prop



  //Worm_t worm_head;                     // in Russian called ira
  //Worm_t worm_tail;                     // in Russian called masha
  bool worm_diag;
  //int worm_dir;                               // direction of the worm, to higher of lower imag times
  //int worm_dir_init;
  //int worm_rising;                           // convention : density is increased or decreased between tail and head
  //int start_bos, start_occ;
  int worm_at_stop;
  int worm_passes_nb_kink;
  bool worm_meas_densmat;
  double worm_dtime;
  bool new_measurement;
  bool transformedHamiltonian;
  
  void find_assoc_insert(const SiteIndex, list<Element_t>::iterator, const int );
  void find_assoc_delete(const SiteIndex, list<Element_t>::iterator );
  bool shift_right_kink(double&, double, double&);
  bool shift_left_kink(double&, double, double&);
  //void site_right(const SiteType s);
  //void site_left(const SiteType s);
  bool t_between(const double, const double, const double) const;
  void worm_left();
  void worm_right();
  void worm_pass_next();
  void worm_pass_prev();
  void worm_changedir_toleft();
  void worm_changedir_toright();
  //int sum_state_pos(const SiteIndex, vector<double>&, vector<double>& );
  //int sum_state_neg(const SiteIndex, vector<double>&, vector<double>& );
  void get_diaginfo_nb_pos(const SiteIndex s);
  void get_diaginfo_nb_neg(const SiteIndex s);
  
  SiteType dist(const SiteType s1, const SiteType s2);
  
  
  // measurement variables
  double mZ_dns, mZ_state;
  double mZ, mG;
  double mZ_green;
  vector<double> av_dns;
  vector<double> av_state;
  vector<double> av_state_sq;
  vector<double> av_green;
  vector<double> time_green;
  vector<double> string_correlation;
  double spinspin_corr;
  vector<double> spinspin_correlation;
  vector<double> Potential_Energy_local;
  vector<double> mWinding;
  IMatrix winding_element;
  IMatrix zcoord;
  IMatrix zco;
 
  
  
 
  DMatrix update_statistics;
  vector<double> nbtime;
  vector<double> dtime;
  vector<StateType> nbval;
  vector<SiteType> nbsite;
  
  
  // lattice variables
  SiteIndex Nsites, latt_dim;
  IMatrix nb;
  IMatrix bond_index;
  IMatrix opposite_direction;
#ifdef UNIFORM
  StateType nmin,nmax;
  std::vector<double> site_weight_diag;
  DMatrix bond_weight_diag;
  DMatrix site_weight_offdiag;
  DTensor bond_weight_offdiag;
  bool range_check(const StateType n) {
    return ( (n >= nmin) && (n <= nmax));
  }
#endif

  
  random_device rd;
  mt19937 MyGenerator;
  uniform_real_distribution<double> rnd;
  
  
  
};


inline void worm::get_diaginfo_nb_pos(const SiteIndex s) {
  //list<Element_t>::iterator it = worm_head_it;
  //size_t zc = MyModel->get_number_couplers(s);
  //vector<SiteIndex> nbsite = MyModel->couplers(s);
  //MyModel->couplers(s, nbsite);
  list<Element_t>::iterator its, itp;
  for (std::size_t i = zcoord[0][s]; i < zcoord[1][s]; i++) {
    //SiteType snb = nbsite[i];
    SiteType snb = nb[s][i];
    its = worm_head_it->get_assoc(i);
    itp = its;
    if (itp == operator_string[snb].begin()) itp=operator_string[snb].end();
    --itp;
    while (! t_between(worm_head_it->time(), itp->time(), its->time()) ) {
      itp = its;
      ++its;
      if (its == operator_string[snb].end() ) its = operator_string[snb].begin();
    }
    if ( worm_head_it->time() == its->time() ) {
      //itp = its; // not needed?
      ++its;
      if (its == operator_string[snb].end() ) its = operator_string[snb].begin();
      dtime[i] = its->time() - worm_head_it->time();
      if (dtime[i] <= 0.) dtime[i] += beta;
      nbtime[i] = its->time();
      nbval[i] = its->before();
    }
    else {
      nbtime[i] = its->time();
      dtime[i] = its->time() - worm_head_it->time();
      if (dtime[i] <= 0.) dtime[i] += beta;
      nbval[i] = its->before();
    }
  }
}

inline void worm::get_diaginfo_nb_neg(const SiteIndex s) {
  //size_t zc = MyModel->get_number_couplers(s);
  //MyModel->couplers(s, nbsite);
  list<Element_t>::iterator its, itp;
  for (std::size_t i = zcoord[0][s]; i < zcoord[1][s]; i++) {
    //SiteType snb = nbsite[i];
    SiteType  snb =  nb[s][i];
    its = worm_head_it->get_assoc(i);
    itp = its;
    if (itp == operator_string[snb].begin()) itp=operator_string[snb].end();
    --itp;
    while (! t_between(worm_head_it->time(), itp->time(), its->time()) ) {
      its = itp;
      if (itp == operator_string[snb].begin() ) itp = operator_string[snb].end();
      --itp;
    }
#ifdef DEBUGMODE
    cout << "\n" << i << "\t" << snb << "\t" << worm_head_it->time() << "\t" << itp->time() << "\t" << its->time() << "\n";
#endif
    nbtime[i] = itp->time();
    dtime[i] = worm_head_it->time() - itp->time();
    if (dtime[i] <= 0.) dtime[i] += beta;
    nbval[i] = itp->after();
#ifdef DEBUGMODE
    cout << i << "\n-\t" << dtime[i] << "\t" << its->time() << "\t" << itp->time()   << "\toccup : " << itp->after() << "\n";
#endif
  }
}





inline bool worm::t_between(const double t0, const double t1, const double t2) const
// check if t0 is between t1 and t2, given that t1 happens before t2 (mod beta)
{
  if (t2 >t1)
  {
    if ((t2>=t0) && (t0 > t1)) return (true);
  }
  else
  {
    if ((t0 > t1) || (t0 <= t2)) return (true);
  }
  return (false);
}



#endif

