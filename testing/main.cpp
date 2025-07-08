#include <fstream>
#include <iostream>
#include "../writefilelib/writefile.hpp"
using namespace std;

int main() {
    
    string filename = "newgridtesting.dat";

    int Nx = 100, Ny = 100;

    RampGrid grid(Nx, Ny, 3, 0.75, 15); 
    cout << 1 << endl;
    cout << grid.iArea(0,0) << endl;
    output_grid(grid, Nx, Ny, filename); 
    cout << 2 << endl;
    return 0;
}