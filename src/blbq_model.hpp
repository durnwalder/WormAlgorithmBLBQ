#include <math.h>

class blbq_model {
public:
  blbq_model(alps::params const& params, Lattice_t& pl) : pLatt(pl),  nmin(0), nmax(2) { 
    model_name = params["model"].as<std::string>();
    double pi = 2*acos(0.0);
    if (params.supplied("theta/pi"))
    {
      alpha = sin(double(params["theta/pi"])*pi);
      gamma = cos(double(params["theta/pi"])*pi);
    }
    else
    {
      alpha = sin(-0.75*pi);
      gamma = cos(-0.75*pi);
    }
      D = params["D"];
  };

  virtual double site_weight_diag(const SiteIndex, const StateType, const bool b = true) const = 0;
  virtual double site_weight_offdiag(const SiteIndex, const StateType, const StateType) const = 0;

  virtual double bond_weight_diag(const BondIndex, const StateType, const StateType) const = 0;
  virtual double bond_weight_offdiag(const BondIndex, const StateType, const StateType, const StateType, const StateType) const = 0;

  void print_params(std::ostream& os) const {
        os << "# model : name                                      : " << model_name << "\n";
        os << "# bilinear coupling constant:                                       : " << gamma    << "\n";
        os << "# biquadratic coupling constant:                                       : " << alpha    << "\n";
        os << "# model : D                                      : " << D    << "\n";
        os << "# model : nmin, nmax                                : " << nmin << "\t" << nmax << "\n";
        os << "# Matrix Elements: " << "\n";
        std::string element;
         for (size_t m1 = 0; m1 <= 2; m1++)
        {
          for (size_t n1 = 0; n1 <= 2; n1++)
          {
            for (size_t m2 = 0; m2 <= 2; m2++)
            {
              for (size_t n2 = 0; n2 <= 2; n2++)
              {
                if (n1 == n2 && m1 == m2)
                {
                  element = std::to_string(bond_weight_diag(0,n1,m1) + site_weight_diag(0,m1)*site_weight_diag(0,m2));
                }
                else
                {
                  element = std::to_string(bond_weight_offdiag(0,n1,n2,m1,m2));
                }
                std::string s(20-element.length(), ' ');
                os <<  element.substr(0,5) << " ";
              }
            }
            os << "\n";
          }
        }
        
    }


  StateType get_nmax() const { return nmax;}
  StateType get_nmin() const { return nmin;}
  
  static void define_parameters(alps::params& params) {
    params
    .define<std::string>("model", "original", "Name of the Hamiltonian. Choose either original or 3color")
    .define<double>("theta/pi", -0.75, "The model original is negative sign free for -1.0 <= theta/pi <= -0.75. The model 3color is negative sign free for -0.75 <= theta/pi <= 0.25.")
    .define<double>("D", 0, "uniaxial field term");
  }

  static bool parameters_supplied(alps::params const& params) {
    return params.supplied("model");
  }

  bool range_fail(const StateType n) const {
    return n < nmin || n > nmax;
  }
  std::string get_name() const { return model_name;}
protected:
  std::string model_name;
  StateType nmax;
  StateType nmin;
  Lattice_t& pLatt;
  double alpha; //bilinear coupling constant
  double gamma; //linear coupling constant
  double D;  
};

class three_color : public blbq_model  {
  public:
    three_color(alps::params const& params, Lattice_t& pl)  : blbq_model(params, pl)
    {}
    double site_weight_offdiag(const SiteIndex s, const StateType n1, const StateType n2) const override {
        return 0;  
    }

    double site_weight_diag(const SiteIndex s, const StateType n, const bool b = true) const override {
          if (n == 2)
          {
            return 0;
          }
          else
          {
            return D;
          } 
    }

    double bond_weight_diag(const BondIndex b, const StateType n, const StateType m) const override {
          if (n == m)
          {
              return alpha;
          }
          else
          {
              return 0;
          }
    
    }
    double bond_weight_offdiag(const BondIndex b, const StateType n1, const StateType n2, const StateType m1, const StateType m2) const override {
        if (range_fail(n1)) return 0;
        if (range_fail(n2)) return 0;
        if (range_fail(m1)) return 0;
        if (range_fail(m2)) return 0;
          if (n1 != n2 && m1 != m2 && n1 == m2 && n2 == m1) //nonzero matrix elements of gammaC
            {
              if (gamma<0)
              {
                return gamma;
              }
              else
              {
                return -gamma;
              }
            }
            else if (n1 == n2 && m1 == m2 && n1 != m1) //nonzero matrix elements of gammaH
            {
              return -(gamma - alpha);
            }
            else
            {
              return 0;
            }
    }
};

class original : public blbq_model  {
  public:
  original(alps::params const& params, Lattice_t& pl)  : blbq_model(params, pl) {};
    double site_weight_offdiag(const SiteIndex s, const StateType n1, const StateType n2)  const override {
        if (range_fail(n1) || range_fail(n2)) return 0;
        return 1;
    }
    double site_weight_diag(const SiteIndex s, const StateType n, const bool b = true)  const override {
        if (n == 1) // n==2 for transformed Hamiltonian, n==1 for original Hamiltonian
        {
          return 0;
        }
        else
        {
          return D;
        }
        
    }

    double bond_weight_diag(const BondIndex b, const StateType n, const StateType m)  const override {
        double weight = 0;

        if(n == 1 && m == 1){
          return alpha;
        }

        if (n == (2-m)) // n==2 for transformed Hamiltonian, n==2-m for original Hamiltonian
        {
            weight += alpha;
        }
        if(n+m == 0 || n+m == 4)
        {
            weight += gamma;
        }
        if(n+m == 2)
        {
            weight += -gamma;
        }
        return weight;
    
    }
    double bond_weight_offdiag(const BondIndex b, const StateType n1, const StateType n2, const StateType m1, const StateType m2)  const override {
        if (range_fail(n1)) return 0;
        if (range_fail(n2)) return 0;
        if (range_fail(m1)) return 0;
        if (range_fail(m2)) return 0;

        if (m2 == 1 && n2 == 1 && m1 == 0 && n1 == 2){
            return gamma - alpha > 0 ? (alpha - gamma) : (gamma - alpha);
          }
        if(m2 == 1 && n2 == 1 && m1 == 2 && n1 == 0){
            return gamma - alpha > 0 ? (alpha - gamma) : (gamma - alpha);
          }
        


        if ((n1 == n2 + 1) && (m2 == m1 + 1)) {
          if (n1 == m1 || m2 == m1)
          {
            return gamma - alpha > 0 ? (alpha - gamma) : (gamma - alpha);
          }
          else
          {
            return gamma;
          }
        }
        else if  ((n2 == n1 + 1) && (m1 == m2 + 1)) {
          if (n1 == m1 || m2 == m1)
          {
            return gamma - alpha > 0 ? (alpha - gamma) : (gamma - alpha);
          }
          else
          {
            return gamma;
          }
        }
        else if ((n1 == n2 + 2) && (m2 == m1 + 2))
        {
          return alpha;
        }
        else if((n2 == n1 + 2) && (m1 == m2 + 2))
        {
          return alpha;
        }
        else{
          return 0;
        }
    }
    };
