#pragma once

#include<vector>
#include<string>
#include<utility>
#include<memory>
#include <alps/mc/mcbase.hpp>

#include<iostream>

#include "defs.hpp"





class Lattice_t  {
public:
  Lattice_t() {}
  const std::size_t get_dim() const { return dim;}
  const std::size_t get_zcoord(const SiteIndex n) const { return zcoord[n];}
  const std::size_t get_Nsites() const { return Nsites;}
  const std::size_t get_Nbonds() const { return Nbonds;}
  std::size_t get_N_nb(const SiteIndex n) const { return neighbor[n].size();}
  std::size_t get_Ls(const size_t d) const { return Ls[d];}
  std::size_t get_pbc(const size_t d) const { return pbc[d];}
  
  SiteIndex operator[](const SiteIndex n) const { return site[n];}
  SiteIndex& operator[](const SiteIndex n) { return site[n];}


  SiteType  get_site(const SiteIndex n) const { return site[n];}
  SiteType& get_site(const SiteIndex n) { return site[n];}
  
  BondType  get_bond(const BondIndex n) const { return bond[n];}
  BondType& get_bond(const BondIndex n) { return bond[n];}
  
  SiteType  nb(const SiteIndex n, const std::size_t m) const { return neighbor[n][m];}
  SiteType& nb(const SiteIndex n, const std::size_t m) { return neighbor[n][m];}
  
  //void define_parameters(alps::params& params) {
  static void define_parameters(alps::params& params) {
    params.define<std::string>("lattice", "chain_lattice", "name of the lattice. Currently implemented are chain_lattice");
    
    /*
    if (params["lattice"].as<std::string>() == "chain_lattice" ) {
      Chain::define_custom_lattice_parameters(params);
    }
    else {
      throw std::runtime_error("Invalid Lattice \"" + ClassifierName + "\"");
    }
    */
  }
  
  
  static bool parameters_supplied(alps::params const& params) {
    return params.supplied("lattice");
  }
  
  
  //virtual static void define_custom_lattice_parameters(alps::params & ) = 0;
  virtual void initialize() = 0;
  
  virtual std::size_t get_direction(const SiteIndex, const SiteIndex) const = 0;
  virtual std::size_t get_opposite_direction(const std::size_t, const SiteIndex) const = 0;
  virtual BondIndex get_bond_index(const SiteIndex n, const size_t dir) const = 0;

  virtual const double get_winding(const size_t dir) const = 0;
  virtual const double get_winding(const BondType b) const = 0;
  
  virtual const SiteType dist(const SiteType s1, const SiteType s2) const = 0;
  virtual void print() = 0;
  
  virtual void print_params(std::ostream&) const = 0;
  
protected:
  std::string ClassifierName;
  std::size_t dim;
  std::vector<std::size_t> zcoord;
  std::vector<std::size_t> Ls;
  std::vector<bool> pbc;
  std::size_t Nsites;
  std::size_t Nbonds;
  std::vector<SiteType> site;
  std::vector<BondType> bond;
  std::vector<std::vector<SiteType> > neighbor;
  std::vector<double> winding_element;
};







class Chain : public Lattice_t {
public:
  Chain(alps::params const& params) {
    dim = 1;
    Ls.resize(dim);
    pbc.resize(dim);
    ClassifierName = params["lattice"].as<std::string>();
    Ls[0]  = params["Lx"];
    pbc[0] = params["pbc"];
    //if (!pbc[0]) throw "Open boundaries not implemented yet.";        // XXX Idea : always pbc but use zero weights on non-existing bonds ???
    Nsites = 1;
    for (std::size_t idim =0; idim < dim; idim++) Nsites *= Ls[idim];
    //Nbonds = (pbc[0] ? Nsites : Nsites-1);
    Nbonds = Nsites;
    site.resize(Nsites);
    bond.resize(Nbonds);
    zcoord.resize(Nsites);
    winding_element.resize(2*dim);
    initialize();
  }
  
  //void define_custom_lattice_parameters(alps::params& params) {
  static void define_custom_lattice_parameters(alps::params & params) {
      params.define<std::size_t>("Lx", 12, "length of the chain")
            .define<bool>("pbc", true, "periodic boundary conditions");
  }
  
  void initialize() override {
    for (SiteIndex i=0; i < Nsites; i++) zcoord[i] = 2*dim;
    neighbor.resize(Nsites);
    for (SiteIndex i=0; i < Nsites; i++) neighbor[i].resize(zcoord[i]);
    for (SiteIndex i=0; i < Ls[0]; i++) site[i] = i;
    for (std::size_t i=0; i < Ls[0]-1; i++) {
      neighbor[i][0] = i+1;
    }
    //if (pbc[0]) neighbor[Ls[0]-1][0] = 0;
    neighbor[Ls[0]-1][0] = 0;
    for (std::size_t i=0; i < Ls[0]-1; i++) neighbor[ neighbor[i][0]  ][1] = i;
    //if (pbc[0]) neighbor[0][1] = Ls[0]-1;
    neighbor[0][1] = Ls[0]-1;
    for (BondIndex i=0; i < Nbonds; i++) bond[i] = std::make_pair(site[i], neighbor[i][0]);
      
    winding_element[0] = 1.;
    winding_element[1] = -1.;
  }

  std::size_t get_direction(const SiteIndex n1, const SiteIndex n2) const override {
    for (std::size_t j=0; j < zcoord[n1] ; j++) {
      if (neighbor[n1][j] == n2) return j;
    }
    return (zcoord[n1]);
  }
  
  BondIndex get_bond_index(const SiteIndex n, const size_t dir) const override {
    return ( (dir == 0) ? n : neighbor[n][dir]);
  }

  const double get_winding(const size_t dir) const override { return winding_element[dir];}
  const double get_winding(const BondType b) const override { return ((b.second == neighbor[b.first][0]) ? winding_element[0] : winding_element[1]);}
  
  const SiteType dist(const SiteType s1, const SiteType s2) const override {
    SiteType s0;
    if (pbc[0])
    {
      s0 = (s2 < s1 ? s2 + Ls[0] - s1 : s2 - s1);			// careful with unsigned int !!
      if (s0 > Nsites/2)
      {
        s0 = Nsites - s0;
      }
    }
    else
    {
      s0 = (s2 < s1 ? s1 - s2 : s2 - s1);
    }
    if (s0 > 17){
      s0 = s0;
    }
    return s0;
  }

  void print() {
    std::cout << "# Nsites, dim, zcoord " << Nsites << " " << dim << " " << zcoord[0] << "\t" << ClassifierName << "\n";
    for (SiteIndex i=0; i < Ls[0]; i++) std::cout << site[i] << "\t" << neighbor[i][0] << "\t" << neighbor[i][1] << "\t" << get_bond_index(i, 0) << "\t" << get_bond_index(i,1) << "\n";
  }
  
  void print_params(std::ostream& os) const override {
    os << "# lattice: name                                     : " << ClassifierName << "\n";
    os << "# lattice: dim                                      : " << dim << "\n";
    os << "# lattice: pbc                                      : " << pbc[0] << "\n";
    os << "# lattice: length                                   : " << Ls[0] << "\t" << Nsites << "\t" << Nbonds << "\n";
  }
  
  std::size_t get_opposite_direction(const std::size_t j, const SiteIndex s) const override {
    return (1-j);
  }
  
};
