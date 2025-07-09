#include <fstream>
#include <iostream>
#include "../writefilelib/writefile.hpp"
using namespace std;

int main() {
    
    string filename = "../plotfiles/newgridtesting.dat";

    int Nx = 200, Ny = 100;

    RampGrid grid(Nx, Ny, 3, 0.75, 15); 

    // cout << "i-face x normal: " << grid.iface_xNorm(Nx / 2, Ny / 2) << "\t i-face y normal: " << grid.iface_yNorm(Nx / 2, Ny / 2);
    // cout << "j-face x normal: " << grid.jface_xNorm(Nx / 2, Ny / 2) << "\t j-face y normal: " << grid.jface_yNorm(Nx / 2, Ny / 2);
    // output_grid(grid, Nx, Ny, filename); 
    return 0;
}