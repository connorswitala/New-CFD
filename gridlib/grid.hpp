#pragma once

#include "../linalglib/linalg.hpp" 
#include <fstream> 
#include <cmath> 
#include <memory>
#include <iomanip> 
 
using namespace std;

// This enum class is for setting boundary conditions types
enum class BCType {
	IsothermalWall,
	AdiabaticWall,
	Inlet,
	Outlet,
	Symmetry,
	Undefined
};

// This struct contains the boundary conditions types for each side of the grid (left, right, bottom, top) 
struct BCMap {

	BCType left;
	BCType right;
	BCType bottom;
	BCType top;

	BCMap(BCType left, BCType right, BCType bottom, BCType top) : left(left), right(right), bottom(bottom), top(top) {}

};

// This sets the inlet flow conditions from inputs in the UI
struct inlet_conditions {
	double rho, u, v, p, T, M, a;
};


struct Point {
    double x;
    double y;

    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {};
};

// double NewtonMethod(double max_dist, int n_points, double d_min);

class Grid {
public:
    virtual double Volume(int i, int j) const = 0;
    virtual double xCenter(int i, int j) const = 0;
    virtual double yCenter(int i, int j) const = 0;
    virtual double xVertex(int i, int j) const = 0;
    virtual double yVertex(int i, int j) const = 0;

    virtual double iface_xNorm(int i, int j) const = 0;
    virtual double iface_yNorm(int i, int j) const = 0;
    virtual double jface_xNorm(int i, int j) const = 0;
    virtual double jface_yNorm(int i, int j) const = 0;

    virtual double iArea(int i, int j) const = 0;
    virtual double jArea(int i, int j) const = 0;


    virtual ~Grid() = default;
};

class RampGrid : public Grid {
private:
    int Nx, Ny;
    double L, inlet_height, ramp_angle;


public:
    RampGrid(int Nx, int Ny, double L, double inlet_height, double ramp_angle);


    Vector x_vertices;
    Vector y_vertices;
    Vector x_cellCenters;
    Vector y_cellCenters;
    Vector iface_xNormals;
    Vector iface_yNormals;
    Vector jface_xNormals;
    Vector jface_yNormals;
    Vector iAreas, jAreas, cellVolumes;  


    double Volume(int i, int j) const override;
    double xCenter(int i, int j) const override;
    double yCenter(int i, int j) const override;
    double xVertex(int i, int j) const override;
    double yVertex(int i, int j) const override;

    double iface_xNorm(int i, int j) const override;
    double iface_yNorm(int i, int j) const override;
    double jface_xNorm(int i, int j) const override;
    double jface_yNorm(int i, int j) const override;

    double iArea(int i, int j) const override;
    double jArea(int i, int j) const override;

};


// class CylinderGrid : public Grid {
// private:
//     int Nx, Ny;
//     double Cylinder_Radius, R1, R2, d_min, theta1, theta2;

//     vector<vector<Point>> vertices;
//     vector<vector<Point>> cellCenters;
//     vector<vector<Point>> iNormals;
//     vector<vector<Point>> jNormals;
//     vector<vector<double>> iAreas, jAreas, cellVolumes;

// public:
//     CylinderGrid(int Nx, int Ny, double Cylinder_Radius, double R1, double R2, double d_min, double theta1, double theta2);

//     double Volume(int i, int j) const override;
//     Point Center(int i, int j) const override;
//     Point Vertex(int i, int j) const override;
//     double iArea(int i, int j) const override;
//     double jArea(int i, int j) const override;
//     Point iNorms(int i, int j) const override;
//     Point jNorms(int i, int j) const override;

// };

// class FlatPlateGrid : public Grid {
// private:
//     int Nx, Ny;
//     double Lx, Ly, dmin;

//     vector<vector<Point>> vertices;
//     vector<vector<Point>> cellCenters;
//     vector<vector<Point>> iNormals;
//     vector<vector<Point>> jNormals;
//     vector<vector<double>> iAreas, jAreas, cellVolumes;

// public:

//     FlatPlateGrid(int Nx, int Ny, double Lx, double Ly, double dmin);

//     double Volume(int i, int j) const override;
//     Point Center(int i, int j) const override;
//     Point Vertex(int i, int j) const override;
//     double iArea(int i, int j) const override;
//     double jArea(int i, int j) const override;
//     Point iNorms(int i, int j) const override;
//     Point jNorms(int i, int j) const override;
// };

// class DoubleConeGrid : public Grid {
// private:
//     int Nx, Ny;
//     double theta1, theta2, l1, l2, l3, l4, inlet_height;

//     vector<vector<Point>> vertices;
//     vector<vector<Point>> cellCenters;
//     vector<vector<Point>> iNormals;
//     vector<vector<Point>> jNormals;
//     vector<vector<double>> iAreas, jAreas, cellVolumes;

// public:

//     DoubleConeGrid(int Nx, int Ny, double l1, double l2, double l3, double l4, double theta1, double theta2, double inlet_height);

//     double Volume(int i, int j) const override;
//     Point Center(int i, int j) const override;
//     Point Vertex(int i, int j) const override;
//     double iArea(int i, int j) const override;
//     double jArea(int i, int j) const override;
//     Point iNorms(int i, int j) const override;
//     Point jNorms(int i, int j) const override;
// };

// class MirroredGrid : public Grid {
// private:
//     int Nx, Ny;
//     double theta1, theta2, l1, l2, l3, l4, inlet_height;

//     vector<vector<Point>> vertices;
//     vector<vector<Point>> cellCenters;
//     vector<vector<Point>> iNormals;
//     vector<vector<Point>> jNormals;
//     vector<vector<double>> iAreas, jAreas, cellVolumes;

// public:

//     MirroredGrid(int Nx, int Ny, double l1, double l2, double l3, double l4, double theta1, double theta2, double inlet_height);

//     double Volume(int i, int j) const override;
//     Point Center(int i, int j) const override;
//     Point Vertex(int i, int j) const override;
//     double iArea(int i, int j) const override;
//     double jArea(int i, int j) const override;
//     Point iNorms(int i, int j) const override;
//     Point jNorms(int i, int j) const override;
// };