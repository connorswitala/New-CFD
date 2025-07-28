#include <fstream>
#include <iostream>
#include "../solverlib/solver.hpp" 
using namespace std;

int main(int argc, char* argv[]) {

    Chemistry chem; 
    chem.compute_equilibrium_fractions(0.1, 6e6);
    chem.display_mass_fractions();

    return 0;

}