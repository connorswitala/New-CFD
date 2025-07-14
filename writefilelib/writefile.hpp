#pragma once

#include <fstream> 
#include <iomanip>
#include <iostream>
#include "../gridlib/grid.hpp"
#include "../linalglib/linalg.hpp"

#define BLUE "\033[34m"
#define RED "\033[31m"
#define CYAN "\033[36m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define MAGENTA "\033[35m"
#define RESET "\033[0m"

// struct ThermoEntry {
// 	double rho, e, p, T, R, cv, gamma, dpdrho, dpde;
// };

// inline string BCTypeToString(BCType type) {
// 	switch (type) {
// 		case BCType::IsothermalWall: return "IsothermalWall";
// 		case BCType::AdiabaticWall: return "AdiabaticWall";
// 		case BCType::Inlet:         return "Inlet";
// 		case BCType::Outlet:        return "Outlet";
// 		case BCType::Symmetry:      return "Symmetry";
// 		case BCType::Undefined:     return "Undefined";
// 		default:                    return "Unknown";
// 	}
// }

// inline BCType stringToBCType(string& str) {
// 		if (str == "IsothermalWall") return BCType::IsothermalWall;
// 		if (str == "AdiabaticWall") return BCType::AdiabaticWall;
// 		if (str == "Inlet") return BCType::Inlet;
// 		if (str == "Outlet") return BCType::Outlet;
// 		if (str == "Symmetry") return BCType::Symmetry;
// 		return BCType::Undefined;
// };

// Vector convert(const Vector& U);
// Vector convert_real(const Vector& U, ThermoEntry thermo);

// void write_perf_data(const int Nx, const int Ny, Vector& U, Vector& U_inlet, Grid& grid, BCMap& BCs, string& gridtype, string filename);
// void write_real_data(const int Nx, const int Ny, Vector& U, Vector U_inlet, Grid& grid, BCMap& BCs, string& gridtype, vector<vector<ThermoEntry>>& cell_thermo, string filename); 

// void cfd_centerline(const int Nx, const int Ny, Vector& U, Grid& grid, string filename);
void output_grid(Grid& grid, const int& Nx, const int& Ny, string& filename);
