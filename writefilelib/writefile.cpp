#include "writefile.hpp"


// Vector convert(const Vector& U) {
//     double gam = 1.4; 

// 	static Vector V(9);
// 	V[0] = U[0]; //density
// 	V[1] = U[1] / U[0]; // u-vel
// 	V[2] = U[2] / U[0]; // v-vel
// 	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (gam - 1); // pressure
//     V[4] = V[3] / (287.0 * V[0]); // temperature
//     V[5] = sqrt(1.4 * V[4] * 287); // speed of sound
//     V[6] = sqrt(V[1] * V[1] + V[2] * V[2]) / V[5]; // mach
//     V[7] = U[3]; // total energy
//     V[8] = V[3] / ((gam - 1) * V[0]); // internal energy

// 	return V;
// } 

// Vector convert_real(const Vector& U, ThermoEntry thermo) {

// 	static Vector V(9);
// 	V[0] = U[0]; // Density
// 	V[1] = U[1] / U[0]; // u-vel
// 	V[2] = U[2] / U[0]; // v-vel
// 	V[3] = (U[3] - 0.5 * V[0] * (V[1] * V[1] + V[2] * V[2])) * (thermo.gamma - 1); // pressure
//     V[4] = thermo.T; // temperature
//     V[5] = sqrt(thermo.gamma * V[4] * thermo.R); // speed of sound a
//     V[6] = sqrt(V[1] * V[1] + V[2] * V[2])/V[5]; // Mach
//     V[7] = U[3]; 
//     V[8] = thermo.e; 

// 	return V;
// } 

void output_grid(Grid& grid, const int& Nx, const int& Ny, string& filename) {
    
    ofstream file(filename);
    file << "Variables = \" x points\", \"y points\" \n";    
    file << "ZONE T=\"Grid\", I=" << Nx << ", J=" << Ny << ", F=POINT\n";

    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {  
            file << grid.xVertex(i, j) << " " << grid.yVertex(i, j) << endl;          
        }
    }

    file.close(); 
}

// void write_perf_data(const int Nx, const int Ny, Vector& U, Vector& U_inlet, Grid& grid, BCMap& BCs, string& gridtype,  string filename) {

//     ofstream file(filename);
//     file << "VARIABLES = \"x\", \"y\", \"density\", \"u-vel\", \"v-vel\", \"energy\", \"pressure\", \"temperature\", \"a\", \"mach\", \"e\" \n";
//     file << "ZONE T=\"Flow Field\", I=" << Nx+1 << ", J=" << Ny+1 << ", F=BLOCK\n";
//     file << "VARLOCATION=([3-11]=CELLCENTERED)\n";
//     file << "# GRIDTYPE = \"" << gridtype << "\"" << endl;  

//     // U_inlet values
//     file << "# U_inlet = ";
//     for (size_t i = 0; i < U_inlet.size(); ++i) {
//         file << U_inlet[i];
//         if (i != U_inlet.size() - 1) file << ", ";
//     }
//     file << "\n";

//     // Boundary condition types
//     file << "# BOUNDARY_CONDITIONS:\n";
//     file << "# Left   = " << BCTypeToString(BCs.left)   << "\n";
//     file << "# Right  = " << BCTypeToString(BCs.right)  << "\n";
//     file << "# Bottom = " << BCTypeToString(BCs.bottom) << "\n";
//     file << "# Top    = " << BCTypeToString(BCs.top)    << "\n";

//     // Write x values
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             file << grid.xVertex(i,j) << " "; 
//             file << "\n";
//         }
//     }


//     // Write y values
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             file << grid.yVertex(i,j) << " "; 
//             file << "\n";
//         }
//     }

//     // Density
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[0] << " ";
//             file << "\n";
//         }
//     }

//     // u-vel
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[1] << " ";
//             file << "\n";
//         }
//     }

//     // v-vel
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[2] << " ";
//             file << "\n";
//         }
//     }

//         // Energy
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[7] << " ";
//             file << "\n";
//         }
//     }

//     // Pressure
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[3] << " ";
//             file << "\n";
//         }
//     }

//     // Temperature
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[4] << " ";
//             file << "\n";
//         }
//     }

//     // Speed of sound
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[5] << " ";
//             file << "\n";
//         }
//     }

//     // Mach
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[6] << " ";
//             file << "\n";
//         }
//     }

//     // Internal Energy
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert(U[i][j])[8] << " ";
//             file << "\n";
//         }
//     }


//     file.close(); 
// }

// void write_real_data(const int Nx, const int Ny, Vector& U, Vector U_inlet, Grid& grid, BCMap& BCs, string& gridtype, vector<vector<ThermoEntry>>& cell_thermo, string filename) {

//     ofstream file(filename);
//     file << "VARIABLES = \"x\", \"y\", \"density\", \"u-vel\", \"v-vel\", \"energy\", \"pressure\", \"temperature\", \"a\", \"mach\", \"e\" \n";
//     file << "ZONE T=\"Flow Field\", I=" << Nx+1 << ", J=" << Ny+1 << ", F=BLOCK\n";
//     file << "VARLOCATION=([3-11]=CELLCENTERED)\n";
//     file << "# GRIDTYPE = \"" << gridtype << "\"" << endl;  

//     // U_inlet values
//     file << "# U_inlet = ";
//     for (size_t i = 0; i < U_inlet.size(); ++i) {
//         file << U_inlet[i];
//         if (i != U_inlet.size() - 1) file << ", ";
//     }
//     file << "\n";

//     // Boundary condition types
//     file << "# BOUNDARY_CONDITIONS:\n";
//     file << "# Left   = " << BCTypeToString(BCs.left)   << "\n";
//     file << "# Right  = " << BCTypeToString(BCs.right)  << "\n";
//     file << "# Bottom = " << BCTypeToString(BCs.bottom) << "\n";
//     file << "# Top    = " << BCTypeToString(BCs.top)    << "\n";
    
//      // Write x values
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             file << grid.Vertex(i,j).x << " "; 
//             file << "\n";
//         }
//     }

//     // Write y values
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             file << grid.Vertex(i,j).y << " "; 
//             file << "\n";
//         }
//     }

//     // Density
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[0] << " ";
//             file << "\n";
//         }
//     }

//     // u-vel
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[1] << " ";
//             file << "\n";
//         }
//     }

//     // v-vel
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[2] << " ";
//             file << "\n";
//         }
//     }


//     // Energy
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[7] << " ";
//             file << "\n";
//         }
//     }


//     // Pressure
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[3] << " ";
//             file << "\n";
//         }
//     }

//     // Temperature
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[4] << " ";
//             file << "\n";
//         }
//     }

//     // Speed of sound
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[5] << " ";
//             file << "\n";
//         }
//     }

//     // Mach
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[6] << " ";
//             file << "\n";
//         }
//     }


//     // Internal Energy
//     for (int j = 1; j <= Ny; ++j) {
//         for (int i = 1; i <= Nx; ++i) {
//             file << convert_real(U[i][j], cell_thermo[i][j])[8] << " ";
//             file << "\n";
//         }
//     }


//     file.close(); 
// }

// void cfd_centerline(const int Nx, const int Ny, Vector& U, Grid& grid, string filename) { 