#pragma once
#include "../linalglib/linalg.hpp"
#include <cstdlib>
#include <fstream>

constexpr int n_species = 10; 
constexpr int n_ele = 3;
constexpr double Ru = 8314;

struct ThermoEntry {
    double rho;
    double e;

    double p;
    double T;
    double R;
    double cv;
    double gamma;
    double dpdrho;
    double dpde;
    double a;
};


class Chemistry {

    private:

    double rho, e, gam, T, p, R_mix, cv_mix, tot_initial_moles;
    Vector CP_0, H_0, S_0, mu_k0, Xk, Yk, R_k, MW, theta_v, h_f, theta_f, q, initial_moles,
    T_coeff, Int_const, high_t_coeff, high_t_int_const, low_t_coeff, low_t_int_const, middle_t_coeff, middle_t_int_const, a;
    ThermoEntry thermo; 

    public:

    Chemistry();

    void compute_mass_fractions();
    void compute_h0();
    void compute_s0();
    void compute_mu0();
    double norm (double* v1, double* v2);
    bool safe_compute_molar_fractions();
    void compute_molar_fractions();
    void compute_equilibrium(); 


    ThermoEntry compute_equilibrium_thermodynamic_variables(double Rho, double E);
    pair<Vector, Vector> compute_equilibrium_concentrations(double Rho, double E); 

    // This functions are independant
    void write_thermochemical_table();
    void plot_concentrations_for_e_range();
    void print_molar_concentrations();
    void print_mass_concentrations(); 
    void print_thermodynamic_variables();
    
}; 