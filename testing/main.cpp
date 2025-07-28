#include <fstream>
#include <iostream>
#include "../gibbslib/gibbs.hpp" 
using namespace std;

int main(int argc, char* argv[]) {

    Chemistry chem; 
    chem.write_thermochemical_table();

    return 0;

}