#include "gibbs.hpp"

/** This is the class constructor and set up all of the necessary coefficients and variables to use in the Gibbs
 *  Free-Energy minimization. 
 */
Chemistry::Chemistry() {

    q = Vector(n_species, 0.0);                     // Vector that holds charge of each species
    a = Vector(3 * n_species, 0.0);                 // Vector that holds stoichiometric coefficients of species
    CP_0 = Vector(n_species, 0.0);                  // Vector that holds CP_0 for each species
    H_0 = Vector(n_species, 0.0);                   // Vector that holds H_0 for each species
    S_0 = Vector(n_species, 0.0);                   // Vector that holds S_0 for each species
    mu_k0 = Vector(n_species, 0.0);                 // Vector that holds mu_k0 for each species
    Xk = Vector(n_species, 0.0);                    // Vector that holds molar fractions of each species
    Yk = Vector(n_species, 0.0);                    // Vector that holds mass fractions of each species
    R_k = Vector(n_species, 0.0);                   // Vector that holds gas constants for each species
    MW = Vector(n_species, 0.0);                    // Vector that holds molecular weight for each species
    theta_v = Vector(3, 0.0);                       // Vector that holds characteristic temperature of vibration for N2, O2, and NO
    theta_f = Vector(n_species, 0.0);               // Vector that holds enthalpies of formation for each species
    T_coeff = Vector(n_species * 7, 0.0);           // Vector that holds temperature coefficient depending on temperature range
    Int_const = Vector(n_species * 2, 0.0);         // Vector that holds integrations constant depending on temperature range
    high_t_coeff = Vector(n_species * 7, 0.0);      // Vector that holds high temp range coefficients
    high_t_int_const = Vector(n_species * 2, 0.0);  // Vector that holds high temp rang integration constants
    middle_t_coeff = Vector(n_species * 7, 0.0);    // Vector that holds middle temp range coefficients
    middle_t_int_const = Vector(n_species * 2, 0.0);// Vector that holds middle temp rang integration constants
    low_t_coeff = Vector(n_species * 7, 0.0);       // Vector that holds low temp range coefficients
    low_t_int_const = Vector(n_species * 2, 0.0);   // Vector that holds low temp rang integration constants

    initial_moles =  { 0.7808, 0.2095, 0.0, 0.0, 0.0, 0.0097, 0.0, 0.0, 0.0, 0.0 };
    tot_initial_moles = 0.0;
    for (int i = 0; i < n_species; ++i) {
        tot_initial_moles += initial_moles[i];
    }

    /**
     * Order of species is: N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-
     */

    MW = { 28.0134, 31.998, 30.008, 14.0067, 15.9994, 39.948, 39.9474514, 14.0061514, 15.9988514, 0.000548579903 };
    theta_v = { 3395.0, 2239.0, 2817.0 };
    theta_f = { 0.0, 0.0, 2.996120e+6, 3.362160e+7, 1.542000e+7, 0.0, 3.82155e7, 1.34337e8, 9.80594e7, 0.0000 };

    for (int i = 0; i < n_species; ++i) {
        R_k[i] = Ru / MW[i]; // Individual gas constants 
    }

    double R_sum = 0.0;

    for (int i = 0; i < n_species; ++i) {
        R_sum += initial_moles[i] * MW[i];
    }

    R_mix = Ru / R_sum;

    a = {2, 0, 1, 1, 0, 0, 0, 1, 0, 0,
         0, 2, 1, 0, 1, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 1, 1, 0, 0, 0}; 

    q = { 0, 0, 0, 0, 0, 0, 1, 1, 1, -1 }; 


    low_t_coeff = {2.210371497e+4, -3.818461820e+2, 6.082738360, -8.530914410e-3, 1.384646189E-05, -9.625793620e-9, 2.519705809e-12,
            -3.425563420e+4, 4.847000970e+2, 1.119010961, 4.293889240e-3, -6.836300520e-7, -2.023372700e-9, 1.039040018e-12,
            -1.143916503e+4, 1.536467592e+2, 3.431468730, -2.668592368e-3, 8.481399120e-6, -7.685111050e-9, 2.386797655e-12,
             0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0,
            -7.953611300e+3, 1.607177787e+2, 1.966226438, 1.013670310e-3, -1.110415423e-6, 6.517507500e-10, -1.584779251e-13,
             0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0 ,
            -57312.09170, 793.079147, -1.717121217, 0.01044184018, -1.180207501e-05, 6.52813478e-09, -1.44755813e-12,
            5237.07921, 2.299958315, 2.487488821, 2.737490756e-05, -3.134447576e-08, 1.850111332e-11, -4.447350984e-15,
            0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0 };
    low_t_int_const = {7.108460860e+2, -1.076003744e+1,
            -3.391454870e+3, 1.849699470e+1,
            9.098214410e+3, 6.728725490,
            5.610463780e+4, 4.193905036,
            2.840362437e+4, 8.404241820,
            -745.375, 4.37967491,
            179057.223, 29.4915095,
            225628.4738, 5.076830786,
            187935.2842, 4.39337676,
            -745.375, -11.72081224 };

    middle_t_coeff = { 5.877124060e+5, -2.239249073e+3, 6.066949220, -6.139685500e-4, 1.491806679e-7, -1.923105485e-11, 1.061954386e-15 ,
                 -1.037939022e+6, 2.344830282e+3, 1.819732036, 1.267847582e-3, -2.188067988e-7, 2.053719572e-11, -8.193467050e-16 ,
                 2.239018716e+5, -1.289651623e+3, 5.433936030, -3.656034900e-4, 9.880966450e-8, -1.416076856e-11, 9.380184620e-16 ,
                 8.876501380e+4, -1.071231500e+2, 2.362188287, 2.916720081e-4, -1.729515100e-7, 4.012657880e-11, -2.677227571e-15 ,
                 2.619020262e+5, -7.298722030e+2, 3.317177270, -4.281334360e-4, 1.036104594e-7, -9.438304330e-12, 2.725038297e-16 ,
                20.10538475, -0.0599266107, 2.500069401, -3.99214116e-08, 1.20527214e-11, -1.819015576e-15, 1.078576636e-19,
                -383596.54, 816.20197, 2.301342628, -4.95298377e-06, 1.205108477e-08, -2.185050286e-12, 1.265493898e-16,
                290497.0374, -855.790861, 3.47738929, -5.28826719e-04, 1.352350307e-07, -1.389834122e-11, 5.046166279e-16,
                -216651.3208, 666.545615, 1.702064364, 4.71499281e-04, -1.427131823e-07, 2.016595903e-11, -9.107157762e-16,
                0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0 };
    middle_t_int_const = {1.283210415e+4, -1.586640027e+1,
            -1.689010929e+4, 1.738716506e+1,
            1.750317656e+4, -8.501669090,
             5.697351330e+4, 4.865231506 ,
            3.392428060e+4, -6.679585350e-1,
            -744.993961, 4.37918011,
            177181.1455, 7.94750748,
            231080.9984, -1.994146545,
            183719.1966, 10.05690382,
            -745.375, -11.72081224 };

    high_t_coeff = { 8.310139160e+8, -6.420733540e+5, 2.020264635e+2, -3.065092046e-2, 2.486903333e-6, -9.705954110e-11, 1.437538881e-15 ,
             4.975294300e+8, -2.866106874e+5, 6.690352250e+1, -6.169959020e-3, 3.016396027e-7, -7.421416600e-12, 7.278175770e-17 ,
             -9.575303540e+8, 5.912434480e+5, -1.384566826e+2, 1.694339403e-2, -1.007351096e-6, 2.912584076e-11, -3.295109350e-16 ,
             5.475181050e+8, -3.107574980e+5, 6.916782740e+1, -6.847988130e-3, 3.827572400e-7, -1.098367709e-11, 1.277986024e-16 ,
             1.779004264e+8, -1.082328257e+5, 2.810778365e+1, -2.975232262e-3, 1.854997534e-7, -5.796231540e-12, 7.191720164e-17 ,
            -9.95126508e+08, 6.45888726e+05, -167.5894697, 0.02319933363, -1.721080911e-06, 6.53193846e-11, -9.740147729e-16,
            10068848.27, -6624.36128, 4.4469082, -3.017567664e-04, 2.612882069e-08, -1.201637769e-12, 2.299206903e-17,
            1.646092148e+07, -11131.65218, 4.97698664, -2.005393583e-04, 1.022481356e-08, -2.691430863e-13, 3.539931593e-18,
            -2.143835383e+08, 1.469518523e+05, -36.8086454, 5.03616454e-03, -3.087873854e-07, 9.18683487e-12, -1.074163268e-16,
            0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0 };
    high_t_int_const = {4.938707040e+6, -1.672099740e+3,
            2.293554027e+6, -5.530621610e+2,
            -4.677501240e+6, 1.242081216e+3,
             2.550585618e+6, -5.848769753e+2 ,
            8.890942630e+5, -2.181728151e+2,
            -5.07830034e+06, 1465.298484,
            234950.4137, -10.32262257,
            313628.4696, -17.06646380,
            -961420.896, 342.619308,
            -745.375, -11.72081224 };


}

/**
 *  ===================== Inner Functions =========================
 *  None of these functions need to be called outside of whatever program you need to use this class for. 
 *  They are only used inside the functions you would call. They are just the inner math.
 */

void Chemistry::compute_mass_fractions() {
        double MW_mix = 0.0;
        for (int i = 0; i < n_species; ++i) {
                MW_mix += Xk[i] * MW[i];
        }

        for (int i = 0; i < n_species; ++i) {
                Yk[i] = (Xk[i] * MW[i]) / MW_mix;
        }
}
void Chemistry::compute_h0() {

        for (int j = 0; j < n_species; ++j) {
                int idx = j * 7; 
                int jdx = j * 2;

                H_0[j] = Ru * T * (- T_coeff[idx] / (T * T) 
                                   + T_coeff[idx + 1] * log(T) / T                    
                                   + T_coeff[idx + 2] 
                                   + T_coeff[idx + 3] * T / 2
                                   + T_coeff[idx + 4] * T * T / 3
                                   + T_coeff[idx + 5] * T * T * T / 4
                                   + T_coeff[idx + 6] * T * T * T * T / 5
                                   + Int_const[jdx] / T);
        }
}
void Chemistry::compute_s0() {
        
        for (int j = 0; j < n_species; ++j) {
                int idx = j * 7; 
                int jdx = j * 2;

                S_0[j] = Ru * (- T_coeff[idx] / (2 * T * T) 
                                   - T_coeff[idx + 1] / T                    
                                   + T_coeff[idx + 2] * log(T)
                                   + T_coeff[idx + 3] * T
                                   + T_coeff[idx + 4] * T * T / 2
                                   + T_coeff[idx + 5] * T * T * T / 3
                                   + T_coeff[idx + 6] * T * T * T * T / 4
                                   + Int_const[jdx + 1] );
        }
}
void Chemistry::compute_mu0() {
        compute_h0();
        compute_s0();

        for (int j = 0; j < n_species; ++j) {
                mu_k0[j] = H_0[j] - T * S_0[j];
        }
}
double Chemistry::norm(double* v1, double* v2) {
        double result = 0.0;
        for (int i = 0; i < 14; ++i) {
                result += fabs(v1[i] - v2[i]) * fabs(v1[i] - v2[i]);
        }
        return sqrt(result); 
}
void Chemistry::compute_molar_fractions() {

        compute_mu0();
        int mat_size = 14;
        Vector X_new(mat_size, 0.0), X(mat_size, 0.0), dx(mat_size, 0.0); 
        double ph = tot_initial_moles / n_species; 

        X = {ph, ph, ph, ph, ph, ph, ph, ph, ph, ph, 1e-4, 1e-4, 1e-4, 1e-4};

        for (int i = 0; i < n_species; ++i) {
                X_new[i] = X[i] / 2; 
        }      
        
        double residual = norm(X_new.data(), X.data()); 

        while (residual > 1e-10) {
      
                X[10] = 0.0; X[11] = 0.0; X[12] = 0.0; X[13] = 0.0;

                Vector J(mat_size * mat_size, 0.0);  // 14x14 matrix, row-major layout 

                // Rows 0–9: species equations
                for (int i = 0; i < n_species; ++i) {

                        J[i * mat_size + i] = 1.0;  // Identity part

                        for (int k = 0; k < n_ele; ++k) {
                                J[i * mat_size + n_species + k] = -a[k * n_species + i];  // -a[k][i] 
                        }

                        J[i * mat_size + mat_size - 1] = -q[i]; 

                }       

                // Row 10–12: nitrogen, oxygen, argon atomic balance
                for (int k = 0; k < n_ele; ++k) {
                        for (int j = 0; j < n_species; ++j) {
                                J[(n_species + k) * mat_size + j] = a[k * n_species + j] * X[j];  // a[k][j] * X[j]
                        }
                // last 4 entries (cols 10–13) remain 0
                }

                // Row 13: charge neutrality
                for (int j = 0; j < n_species; ++j) {
                        J[(mat_size - 1) * mat_size + j] = q[j] * X[j];
                }
                // last 4 entries of row 13 are 0

                Vector F(mat_size, 0.0); 
                for (int i = 0; i < n_species; ++i) {
                        double Xi_safe = max(X[i], 1e-10); // Protect log from zero  
                        F[i] = -(mu_k0[i] + Ru * T * log(Xi_safe * p / 101325)) / (Ru * T);
                }

                F[10] = initial_moles[0] - (2 * X[0] + X[2] + X[3] + X[7]); // Nitrogen  
                F[11] = initial_moles[1] - (2 * X[1] + X[2] + X[4] + X[8]); // Oxygen
                F[12] = initial_moles[5] - (X[5] + X[6]); // Argon
                F[13] = -(X[6] + X[7] + X[8] - X[9]); // Electron 

                matrix_divide(J.data(), F.data(), dx.data(), mat_size, 1);

                
                // Update molar fractions safely
                for (int i = 0; i < mat_size; ++i) {

                        double dx_safe = min(max(dx[i], -50.0), 50.0); // Clamp dx 
                        X_new[i] = X[i] * exp(dx_safe);

                        if (isnan(X_new[i]) || isinf(X_new[i])) {
                                throw std::runtime_error("Nonphysical molar fraction computed (NaN or Inf)");
                        }

                        if (X_new[i] < 0) {
                                throw std::runtime_error("Negative molar fraction computed");
                        }
                }

                residual = norm(X_new.data(), X.data());
                X = X_new;
        }

        double sum = 0.0;

        for (int i = 0; i < n_species; ++i) {
                sum += X_new[i];
        }

        for (int i = 0; i < n_species; ++i) {
                Xk[i] = X_new[i] / sum; // Normalize molar fractions
        }
        
}
bool Chemistry::safe_compute_molar_fractions() {

        try {
                compute_molar_fractions();
                return true;
        }
        catch (const exception& ex) {
                cerr << "[Warning] Chemical equilibrium failed: " << ex.what() << std::endl;

                // Fallback: assume frozen air (N2 and O2 only)
                // N2 = 78%, O2 = 21%, Ar = 0.97%, everything else 0%
                Xk = { 0.7808, 0.2095, 0.0, 0.0, 0.0, 0.0097, 0.0, 0.0, 0.0, 0.0 };
                return false;
        }

}
void Chemistry::compute_equilibrium() {

        double e_new = 0, cv_new = 0;

        while (fabs(e_new - e) >= 1) {
   
                if (T > 200.0 && T < 1000.0) {
                        T_coeff = low_t_coeff;
                        Int_const = low_t_int_const;
                }
                else if (T >= 1000.0 && T < 6000.0) {
                        T_coeff = middle_t_coeff;
                        Int_const = middle_t_int_const;
                }
                else {
                        T_coeff = high_t_coeff;
                        Int_const = high_t_int_const;
                }

                bool equilibrium_success = safe_compute_molar_fractions();
                compute_mass_fractions();

                if (equilibrium_success == false) {
                        R_mix = 287.0;
                        cv_mix = 717;
                        T = e / cv_mix;
                        p = rho * R_mix * T;
                        break;
                }
                else {

                        e_new = 0.0;

                        for (int i = 0; i < 3; ++i) {
                                e_new += Yk[i] * (2.5 * R_k[i] * T + R_k[i] * theta_v[i] / (exp(theta_v[i] / T) - 1) + theta_f[i]); //diatomic species
                        }

                        for (int i = 3; i < n_species; ++i) {
                                e_new += Yk[i] * (1.5 * R_k[i] * T + theta_f[i]); // atomic species
                        }

                        cv_new = e_new / T;

                        T = T - 0.1 * (e_new - e) / cv_new;


                        double sum = 0.0;
                        for (int i = 0; i < n_species; ++i) {
                                sum += Xk[i] * MW[i];
                        }

                        R_mix = Ru / sum; // Mixture gas constant  
                        p = rho * R_mix * T;
                }
        }

        cv_mix = cv_new;
}
 

/** 
 *  ======================= Main Functions =========================
 *  All of these are the main functions that you would want to call.
 */

 

/** The following funcion computes the thermochemical variables that are inside the ThermoEntry struct. It minimizes the 
 *  Gibbs Free-Energy 5 separate times. The first time calculates the main therodynamic equation of state, while the remaining
 *  four calculations are executed from a perturbed state in order to calculate the derivatives of pressure WRT density and 
 *  internal energy using a basic centered finite-differencing method. 
 */
ThermoEntry Chemistry::compute_equilibrium_thermo_vars(double Rho, double E) {

        rho = Rho;
        e = E;
        T = e / 717;
        if (e > 1e7) T = 20000; 	
        p = rho * R_mix * T;  // Initial pressure estimate	  

        compute_equilibrium();
        gam = 1 + R_mix / cv_mix;

        //////  Compute derivatives ////

        // dp/de | rho
        double delta_e = 1e-5 * E;
        e = E + delta_e;
        T = e / 717;
        if (e > 1e7) T = 20000;
        p = rho * 287.0 * T;
        compute_equilibrium();
        double p1 = p;


        e = E - delta_e;
        T = e / 717;
        if (e > 1e7) T = 20000;
        p = rho * 287.0 * T;
        compute_equilibrium();
        double p2 = p;

        double dpde = (p1 - p2) / (2 * delta_e); // dp/de 


        // dp/drho| e
        double delta_rho = 1e-4 * rho;
        rho = Rho + delta_rho;
        T = E / 717;
        if (e > 1e7) T = 20000;
        p = rho * 287.0 * T;
        compute_equilibrium();
        p1 = p;

        rho = Rho - delta_rho;
        T = E / 717;
        if (e > 1e7) T = 20000;
        p = rho * 287.0 * T;
        compute_equilibrium();
        p2 = p;
        double dpdrho = (p1 - p2) / (2 * delta_rho); // dp/drho

        thermo.rho = rho; thermo.e = e;
        thermo.p = p; thermo.T = T;
        thermo.R = R_mix, thermo.cv = cv_mix;
        thermo.gamma = gam;
        thermo.dpdrho = dpdrho; thermo.dpde = dpde; 

        return thermo;
}



/** The following function takes in a density and internal energy and returns two vectors. The first contains the mass
 *  fractions and the second contains the molar fractions. The final entry of each vector also contains the computed
 *  equilibirum temperature. 
 */
pair<Vector, Vector> Chemistry::compute_equilibrium_fractions(double Rho, double E) {
        
        rho = Rho;
        e = E;
        T = e / 717;
        if (e > 1e7) T = 20000;
        p = rho * R_mix * T;

        compute_equilibrium();

        Vector mass_fractions(n_species + 1, 0.0), molar_fractions(n_species + 1, 0.0); 

        for (int i = 0; i < n_species; ++i) {
                mass_fractions[i] = Yk[i];
                molar_fractions[i] = Xk[i];
        }

        mass_fractions[n_species] = T;
        molar_fractions[n_species] = T; 

        return make_pair(mass_fractions, molar_fractions); 

}



/** 
*   The following function writes the thermochemical table to a .csv file for a range of rho and e values.  
*   The thermochemical table contains density, internal energy, pressure, temperature, R_mix, Cv_mix,
*   gamma_mix, dpdrho, and dpde. Only this function is necessary after calling the chemistry constructor 
*   "Chemistry chem;"
*/
void Chemistry::write_thermochemical_table() {
        int n = 500; 

        Vector e(n), rho(n); 
        for (int i = 0; i < n; ++i) {
                e[i] = 717 * 600 + (5e7 - 717 * 600) / (n - 1) * i;
                rho[i] = 1e-3 + (2 - 1e-3)/(n - 1) * i;
        }

        Chemistry chem; 
        ofstream file("thermochemical_table.csv");
        file << "rho, e, p, T, R, cv, gam, dpdrho, dpde" << endl;  
        ThermoEntry holder; 

        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                        holder = compute_equilibrium_thermo_vars(rho[i], e[j]);
                        file << holder.rho << ", " << holder.e << ", " << holder.p << ", " << holder.T << ", " << holder.R
                                << ", " << holder.cv << ", " << holder.gamma << ", " << holder.dpdrho << ", " << holder.dpde << endl; 
                }
                cout << fixed << setprecision(4) << 100.0 * i * n / (n * n) << "% complete" << endl;

        }

        file.close(); 
}

/** 
 *  The following function is a standalone function that only requires the the constructor "Chemistry chem;" be place in you main.cpp. 
 *  When you call this function is will create a .csv file for plotting in paraview that plots the mass and molar fractions of
 *  Each species against the equilibirum temperature until T_eq = 20,000 K. 
 */
void Chemistry::plot_fractions() {
        double rho = 0.1, e;
        double T_eq = 600;
        Vector maf(11), mof(11); 

        ofstream mass("../plotfiles/mass_fractions.csv");
	ofstream molar("../plotfiles/molar_fractions.csv");

        mass << "T, N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-" << endl;
        molar << "T, N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-" << endl;

        int counter = 0;
        while (T_eq < 20000) {
                e = 717 * 600 + 20000 * counter;
                auto result = compute_equilibrium_fractions(rho, e);
                maf = result.first;
                mof = result.second; 

                mass << maf[10]; // Temperature
		molar << mof[10]; // Temperature
                for (int i = 0; i < 10; ++i) mass << ", " << maf[i];  // Mass fractions 
                for (int i = 0; i < 10; ++i) molar << ", " << mof[i];  // Molar fractions 
                mass << "\n";
                molar << "\n";

                counter++; 
                T_eq = maf[10];  // Assuming maf[10] is temperature
                if (counter % 100 == 0) cout << T_eq << endl;
        }

        mass.close();
	molar.close();
}



/** 
 *  The following two functions are just to display molar and mass fractions if you only want to compute a state for a certain
 *  rho and e combination. You would call compute_equilibrium_fractions(rho, e) and then use these to display the results. Note
 *  that the molar and mass fractions are member functions, so they are stored in the Chemistry class. The 
 *  compute_equilibirum_fractions(rho, e) returns a vector for the plot_fractions() function because it needs to iterate over many
 *  densities and internal energies.
 */
void Chemistry::display_molar_fractions() {
        cout << "Molar Fractions: " << endl << endl;
        cout << "N2: " << Xk[0] << endl;
        cout << "O2: " << Xk[1] << endl;
        cout << "NO: " << Xk[2] << endl;
        cout << "N:  " << Xk[3] << endl;
        cout << "O:  " << Xk[4] << endl;
        cout << "Ar: " << Xk[5] << endl;
        cout << "Ar+: " << Xk[6] << endl;
        cout << "N+: " << Xk[7] << endl;
        cout << "O+:  " << Xk[8] << endl;
        cout << "e-:  " << Xk[9] << endl;
}
void Chemistry::display_mass_fractions() {
        cout << "Mass Fractions: " << endl << endl;
        cout << "N2:  " << Yk[0] << endl;
        cout << "O2:  " << Yk[1] << endl;
        cout << "NO:  " << Yk[2] << endl;
        cout << "N:   " << Yk[3] << endl;
        cout << "O:   " << Yk[4] << endl;
        cout << "Ar:  " << Yk[5] << endl;
        cout << "Ar+: " << Yk[6] << endl;
        cout << "N+:  " << Yk[7] << endl;
        cout << "O+:  " << Yk[8] << endl;
        cout << "e-:  " << Yk[9] << endl;
        cout << "Temperature: " << T << endl;
}



/** 
 *  The following functions displays the thermochemical variables that are found from a given density and internal energy.
 *  You would call 'compute_equilibirium_thermo_vars(rho, e)' and then use this function to display the results. The results 
 *  are stored in a member variables "thermo", but the function 'compute_equilibrium_thermo_vars(rho, e) returns a ThermoEntry 
 *  for the sake of creating the thermochemical table used in CFD simulations.
 */
void Chemistry::display_thermo_vars() {
        cout << endl << "Thermochemical Variables: " << endl << endl;
        cout << setw(30) << "Density: " << thermo.rho << " kg/m^3\n";
        cout << setw(30) << "Internal energy: " << thermo.e << " J/kg\n";
        cout << setw(30) << "Pressure: " << thermo.p << " Pa\n";
        cout << setw(30) << "Equilibrium Temperature: " << thermo.T << " K\n";
        cout << setw(30) << "Mixture Gas Constant R: " << thermo.R << " J/(kg K)\n";
        cout << setw(30) << "Mixture Specific Heat Cv: " << thermo.cv << " J/(kg K)\n";
        cout << setw(30) << "Mixture gamma: " << thermo.gamma << '\n';
        cout << setw(30) << "dp/drho: " << thermo.dpdrho << '\n';
        cout << setw(30) << "dp/de: " << thermo.dpde << '\n';
}