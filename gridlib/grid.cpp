#include "grid.hpp" 



double NewtonMethod(double max_dist, int n_points, double d_min) {
	double k = 1, k_new = 1 / 2, ratio = fabs(k - k_new);
	double func, func_prime;

	while (ratio >= 0.00000000001) {
		func = d_min - max_dist * (exp(k / (n_points - 1)) - 1) / (exp(k) - 1);
		func_prime = -max_dist * (((1 / (n_points - 1) * exp(k / (n_points - 1))) * (exp(k) - 1) - (exp(k / (n_points - 1)) - 1) * exp(k)) / ((exp(k) - 1) * (exp(k) - 1)));
		k_new = k - func / func_prime;
		ratio = fabs(k - k_new);
		k = k_new;
	}

	return k;
}

/////////////////////////////////////////////////
///////////// Ramp Grid functions ///////////////
/////////////////////////////////////////////////

RampGrid::RampGrid(int Nx, int Ny, double L, double inlet_height, double ramp_angle)
	: Nx(Nx), Ny(Ny), L(L), inlet_height(inlet_height), ramp_angle(ramp_angle), 
		x_vertices( (Nx + 1) * (Ny + 1), 0.0),
		y_vertices( (Nx + 1) * (Ny + 1), 0.0),
		x_cellCenters(Nx * Ny, 0.0),
		y_cellCenters(Nx * Ny, 0.0),
		iface_xNormals((Nx + 1) * Ny, 0.0),
		iface_yNormals((Nx + 1) * Ny, 0.0),
		jface_xNormals(Nx * (Ny + 1), 0.0),
		jface_yNormals(Nx * (Ny + 1), 0.0),
		iAreas((Nx + 1) * Ny, 0.0),
		jAreas(Nx * (Ny + 1), 0.0),
		cellVolumes(Nx * Ny, 0.0)  {

		int i,j;
		double dx = L / Nx;
		double theta_rad = ramp_angle * 3.141592653 / 180.0;
		double slope = tan(theta_rad);
		double ramp_start_x = L / 3.0;
		double ramp_start_y = 0.0; 
		double ramp_end_x = 2.0 * L / 3.0;
		double ramp_length = ramp_end_x - ramp_start_x;
		double ramp_end_y = slope * ramp_length;  // height at end of ramp

		for (int i = 0; i <= Nx; ++i) {
			for (int j = 0; j <= Ny; ++j) {

				int x_index = j * (Nx + 1) + i;
				double x = i * dx;

				x_vertices[x_index] = x;  
				if ( (x < ramp_start_x) && (x + dx > ramp_start_x) ) ramp_start_x = x;
				if ( (x < ramp_end_x) && (x + dx > ramp_end_x) ) {
					ramp_end_x = x;
					ramp_length = ramp_end_x - ramp_start_x;
					ramp_end_y = slope * ramp_length;
				}
			}
		}

		for (int i = 0; i <= Nx; ++i) {

			double x = i * dx;
			double y_bot, y_top;

			y_top = inlet_height;

			if (x <= ramp_start_x) {
				y_bot = 0.0;
			} else if (x <= ramp_end_x) {
				double x_rel = x - ramp_start_x;
				y_bot = slope * x_rel;
			} else {
				y_bot = ramp_end_y;
			}

			for (int j = 0; j <= Ny; ++j) {

				int y_index = j * (Nx + 1) + i;
				double frac = static_cast<double>(j) / Ny;
				
				y_vertices[y_index] = y_bot + frac * (y_top - y_bot);
				
			}
		}
	

	// Edge vectors
	Point AB, BC, CD, DA;

	// Calculates cell centers and volumes
	for (i = 0; i < Nx; ++i) {
		for (j = 0; j < Ny; ++j) {

			int ij = j * (Nx + 1) + i;
			int iij = j * (Nx + 1) + i + 1;  
			int ijj = (j + 1) * (Nx + 1) + i;
			int iijj = (j + 1) * (Nx + 1) + i + 1;  

			DA = { x_vertices[ij] - x_vertices[ijj], y_vertices[ij] - y_vertices[ijj] }; 
			AB = { x_vertices[iij] - x_vertices[ij], y_vertices[iij] - y_vertices[ij] };
			BC = { x_vertices[iijj] - x_vertices[iij], y_vertices[iijj] - y_vertices[iij] };
			CD = { x_vertices[ijj] - x_vertices[iijj], y_vertices[ijj] - y_vertices[iijj] }; 

			int cell_ij = j * Nx + i; 

			x_cellCenters[cell_ij] = (x_vertices[ij] + x_vertices[iij] + x_vertices[iijj] + x_vertices[ijj]) / 4;
			y_cellCenters[cell_ij] =	(y_vertices[ij] + y_vertices[iij] + y_vertices[iijj] + y_vertices[ijj]) / 4;

			cellVolumes[cell_ij] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);

			// cout << "xcenter: " << fixed << setprecision(3) << x_cellCenters[cell_ij] 
			// << "\tycenter: " << fixed << setprecision(3) << y_cellCenters[cell_ij] 
			// << "\tcell volume: " << fixed << setprecision(3) << cellVolumes[cell_ij] << endl; 
		}
	}


	// Calculates geometries for i-faces
	for (i = 0; i <= Nx; ++i) {
		for (j = 0; j < Ny; ++j) {

			int ij = j * (Nx + 1) + i;
			int ijj = (j + 1) * (Nx + 1) + i; 

			int face_ij = j * (Nx + 1) + i;

			AB = { x_vertices[ijj] - x_vertices[ij], y_vertices[ijj] - y_vertices[ij] };

			iAreas[face_ij] = sqrt(AB.x * AB.x + AB.y * AB.y); 

			iface_xNormals[face_ij] = AB.y / fabs(iAreas[face_ij]);
			iface_yNormals[face_ij] = AB.x / fabs(iAreas[face_ij]);

			// cout << "i area: " << fixed << setprecision(3) << iAreas[face_ij] 
			// << "\ti-x norm: " << fixed << setprecision(3) << iface_xNormals[face_ij] 
			// << "\ti-y norm: " << fixed << setprecision(3) << iface_yNormals[face_ij] << endl;
		}
	}

	// Calculates geometries for j-faces
	for (i = 0; i < Nx; ++i) {
		for (j = 0; j <= Ny; ++j) {

			int ij = j * (Nx + 1) + i;
			int iij = j * (Nx + 1) + i + 1;

			CD = { x_vertices[iij] - x_vertices[ij], y_vertices[iij] - y_vertices[ij] };

			int face_ij = j * Nx + i;

			jAreas[face_ij] = sqrt(CD.x * CD.x + CD.y * CD.y);

			jface_xNormals[face_ij] = -CD.y / fabs(jAreas[face_ij]);
			jface_yNormals[face_ij] = CD.x / fabs(jAreas[face_ij]);

			// cout << "j area: " << fixed << setprecision(3) <<  jAreas[face_ij] 
			// << "\tj-x norm: " << fixed << setprecision(3) << jface_xNormals[face_ij] 
			// << "\tj-y norm: " << fixed << setprecision(3) << jface_yNormals[face_ij] << endl;
		}
	}

	cout << "Grid generation complete!" << endl;
}


// Cell centers
inline double RampGrid::xCenter(int i, int j) const {
	return x_cellCenters[j * Nx + i];
}

inline double RampGrid::yCenter(int i, int j) const {
	return y_cellCenters[j * Nx + i];
}

// Vertices
inline double RampGrid::xVertex(int i, int j) const {
	return x_vertices[j * (Nx + 1) + i];
}

inline double RampGrid::yVertex(int i, int j) const {
	return y_vertices[j * (Nx + 1) + i];
}

// Cell volumes
inline double RampGrid::Volume(int i, int j) const {
	return cellVolumes[j * Nx + i];
}

// i-face areas and normals
inline double RampGrid::iArea(int i, int j) const {
	return iAreas[j * (Nx + 1) + i];
}

inline double RampGrid::iface_xNorm(int i, int j) const {
	return iface_xNormals[j * (Nx + 1) + i];
}

inline double RampGrid::iface_yNorm(int i, int j) const {
	return iface_yNormals[j * (Nx + 1) + i];
}

// j-face areas and normals
inline double RampGrid::jArea(int i, int j) const {
	return jAreas[j * Nx + i];
}

inline double RampGrid::jface_xNorm(int i, int j) const {
	return jface_xNormals[j * Nx + i];
}

inline double RampGrid::jface_yNorm(int i, int j) const {
	return jface_yNormals[j * Nx + i];
}


/////////////////////////////////////////////////
/////////// Cylinder Grid functions /////////////
/////////////////////////////////////////////////


// CylinderGrid::CylinderGrid(int Nx, int Ny, double Cylinder_Radius, double R1, double R2, double d_min, double theta1, double theta2) :
// 	Nx(Nx), Ny(Ny), Cylinder_Radius(Cylinder_Radius), R1(R1), R2(R2), d_min(d_min), theta1(theta1), theta2(theta2),
// 	vertices(Nx + 1, vector<Point>(Ny + 1)), cellCenters(Nx, vector<Point>(Ny)), cellVolumes(Nx, vector<double>(Ny)),
// 	iAreas(Nx + 1, vector<double>(Ny)), jAreas(Nx, vector<double>(Ny + 1)), iNormals(Nx + 1, vector<Point>(Ny)), jNormals(Nx, vector<Point>(Ny + 1)) {


// 	const int Ntheta = Nx + 1, Nr = Ny + 1;
// 	double R_max, k1;

// 	vector<double> theta(Ntheta, 0.0);
// 	vector<vector<double>> r(Ntheta, vector<double>(Nr, 0.0));

// 	double dtheta = (theta2 - theta1) / (Ntheta - 1);

// 	for (int i = 0; i < Ntheta; ++i) {
// 		r[i][0] = Cylinder_Radius;
// 	}

// 	for (int i = 0; i < Ntheta; ++i) {
// 		theta[i] = theta2 - i * dtheta;
// 	}

// 	for (int i = 0; i < Ntheta; ++i) {

// 		R_max = R1 + (R2 - R1) * cos(theta[i]);
// 		k1 = NewtonMethod(R_max, Nr, d_min);

// 		for (int j = 1; j < Nr; ++j) {
// 			r[i][j] = r[i][0] + R_max * ((exp(k1 * j / (Nr - 1)) - 1) / (exp(k1) - 1));
// 		}

// 	}

// 	for (int i = 0; i <= Nx; ++i) {
// 		for (int j = 0; j <= Ny; ++j) {
// 			vertices[i][j].x = r[i][j] * cos(theta[i]);
// 			vertices[i][j].y = r[i][j] * sin(theta[i]);
// 		}
// 	}


// 	Point AB, BC, CD, DA;
// 	int i, j;

// 	// Calculates cell centers and volumes
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {
// 			DA = { vertices[i][j].x - vertices[i][j + 1].x, vertices[i][j].y - vertices[i][j + 1].y };
// 			AB = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };
// 			BC = { vertices[i + 1][j + 1].x - vertices[i + 1][j].x, vertices[i + 1][j + 1].y - vertices[i + 1][j].y };
// 			CD = { vertices[i][j + 1].x - vertices[i + 1][j + 1].x, vertices[i][j + 1].y - vertices[i + 1][j + 1].y };

// 			cellCenters[i][j] = { (vertices[i][j].x + vertices[i + 1][j].x + vertices[i + 1][j + 1].x + vertices[i][j + 1].x) / 4, (vertices[i][j].y + vertices[i + 1][j].y + vertices[i + 1][j + 1].y + vertices[i][j + 1].y) / 4 };
// 			cellVolumes[i][j] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);
// 		}
// 	}



// 	// Calculates geometries for i-faces
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {

// 			AB = { vertices[i][j + 1].x - vertices[i][j].x, vertices[i][j + 1].y - vertices[i][j].y };

// 			iAreas[i][j] = sqrt(AB.x * AB.x + AB.y * AB.y);

// 			iNormals[i][j].x = AB.y / fabs(iAreas[i][j]);
// 			iNormals[i][j].y = -AB.x / fabs(iAreas[i][j]);
// 		}
// 	}

// 	// Calculates geometries for j-faces
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {

// 			CD = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };

// 			jAreas[i][j] = sqrt(CD.x * CD.x + CD.y * CD.y);

// 			jNormals[i][j].x = -CD.y / fabs(jAreas[i][j]);
// 			jNormals[i][j].y = CD.x / fabs(jAreas[i][j]);
// 		}
// 	}
// }



// Point CylinderGrid::Center(int i, int j) const {
// 	return cellCenters[i][j];
// }

// Point CylinderGrid::Vertex(int i, int j) const {
// 	return vertices[i][j];
// }

// double CylinderGrid::Volume(int i, int j) const {
// 	return cellVolumes[i][j];
// }

// double CylinderGrid::iArea(int i, int j) const {
// 	return iAreas[i][j];
// }

// double CylinderGrid::jArea(int i, int j) const {
// 	return jAreas[i][j];
// }

// Point CylinderGrid::iNorms(int i, int j) const {
// 	return iNormals[i][j];
// }
// Point CylinderGrid::jNorms(int i, int j) const {
// 	return jNormals[i][j];
// }



// ///////////////////////////////////////////////////
// /////////////// Square Grid functions /////////////
// ///////////////////////////////////////////////////

// FlatPlateGrid::FlatPlateGrid(int Nx, int Ny, double Lx, double Ly, double dmin)
// 	: Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly), dmin(dmin),
// 	vertices(Nx + 1, vector<Point>(Ny + 1)), cellCenters(Nx, vector<Point>(Ny)),
// 	cellVolumes(Nx, vector<double>(Ny)), iAreas(Nx + 1, vector<double>(Ny, 0.0)),
// 	jAreas(Nx, vector<double>(Ny + 1, 0.0)), iNormals(Nx + 1, vector<Point>(Ny)), jNormals(Nx, vector<Point>(Ny + 1)) {


// 	int i, j;

// 	// Define lengths 
// 	double dx = Lx / (Nx + 1);
// 	double k = NewtonMethod(Ly, Ny, dmin);

// 	// Create x vertices
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {
// 			vertices[i][j].x = i * dx;
// 			vertices[i][j].y = 0 + Ly * ((exp(k * j / Ny) - 1) / (exp(k) - 1));
// 		}
// 	}

// 	// Edge vectors
// 	Point AB, BC, CD, DA;

// 	// Calculates cell centers and volumes
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {
// 			DA = { vertices[i][j].x - vertices[i][j + 1].x, vertices[i][j].y - vertices[i][j + 1].y };
// 			AB = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };
// 			BC = { vertices[i + 1][j + 1].x - vertices[i + 1][j].x, vertices[i + 1][j + 1].y - vertices[i + 1][j].y };
// 			CD = { vertices[i][j + 1].x - vertices[i + 1][j + 1].x, vertices[i][j + 1].y - vertices[i + 1][j + 1].y };

// 			cellCenters[i][j] = { (vertices[i][j].x + vertices[i + 1][j].x + vertices[i + 1][j + 1].x + vertices[i][j + 1].x) / 4,
// 									(vertices[i][j].y + vertices[i + 1][j].y + vertices[i + 1][j + 1].y + vertices[i][j + 1].y) / 4 };

// 			cellVolumes[i][j] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);
// 		}
// 	}



// 	// Calculates geometries for i-faces
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {

// 			AB = { vertices[i][j + 1].x - vertices[i][j].x, vertices[i][j + 1].y - vertices[i][j].y };

// 			iAreas[i][j] = sqrt(AB.x * AB.x + AB.y * AB.y);

// 			iNormals[i][j].x = AB.y / fabs(iAreas[i][j]);
// 			iNormals[i][j].y = AB.x / fabs(iAreas[i][j]);
// 		}
// 	}


// 	// Calculates geometries for j-faces
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {

// 			CD = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };

// 			jAreas[i][j] = sqrt(CD.x * CD.x + CD.y * CD.y);

// 			jNormals[i][j].x = CD.y / fabs(jAreas[i][j]);
// 			jNormals[i][j].y = CD.x / fabs(jAreas[i][j]);
// 		}
// 	}



// }



// Point FlatPlateGrid::FlatPlateGrid::Center(int i, int j) const {
// 	return cellCenters[i][j];
// }

// Point FlatPlateGrid::FlatPlateGrid::Vertex(int i, int j) const {
// 	return vertices[i][j];
// }

// double FlatPlateGrid::FlatPlateGrid::Volume(int i, int j) const {
// 	return cellVolumes[i][j];
// }

// double FlatPlateGrid::FlatPlateGrid::iArea(int i, int j) const {
// 	return iAreas[i][j];
// }

// double FlatPlateGrid::FlatPlateGrid::jArea(int i, int j) const {
// 	return jAreas[i][j];
// }


// Point FlatPlateGrid::FlatPlateGrid::iNorms(int i, int j) const {
// 	return iNormals[i][j];
// }
// Point FlatPlateGrid::FlatPlateGrid::jNorms(int i, int j) const {
// 	return jNormals[i][j];
// }



// ///////////////////////////////////////////////////
// /////////////// Double Cone Grid functions ////////
// ///////////////////////////////////////////////////

// DoubleConeGrid::DoubleConeGrid(int Nx, int Ny, double l1, double l2, double l3, double l4, double theta1, double theta2, double inlet_height) : Nx(Nx), Ny(Ny),
// l1(l1), l2(l2), l3(l3), l4(l4), theta1(theta1), theta2(theta2), inlet_height(inlet_height), vertices(Nx + 1, vector<Point>(Ny + 1)), cellCenters(Nx, vector<Point>(Ny)),
// cellVolumes(Nx, vector<double>(Ny)), iAreas(Nx + 1, vector<double>(Ny, 0.0)), jAreas(Nx, vector<double>(Ny + 1, 0.0)), iNormals(Nx + 1, vector<Point>(Ny)), jNormals(Nx, vector<Point>(Ny + 1)) {

// 	int i, j;

// 	// Important constants
// 	double deg_to_rads = 3.141592653 / 180.0;

// 	// Define lengths 
// 	double L_total = l1 + l2 + l3 + l4;
// 	double dx = L_total / (Nx + 1);
// 	double dy = inlet_height / (Ny + 1);
// 	double L, dy_ramp;

// 	// Create x vertices
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {
// 			vertices[i][j].x = i * dx;
// 		}
// 	}

// 	// Snap x-vertices to important boundary points.
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {
// 			if (vertices[i][j].x > l1 && vertices[i][j].x < l1 + dx) {
// 				vertices[i][j].x = l1;
// 			}
// 			else if (vertices[i][j].x > l1 + l2 && vertices[i][j].x < l1 + l2 + dx) {
// 				vertices[i][j].x = l1 + l2;
// 			}
// 			else if (i == Nx) {
// 				vertices[i][j].x = L_total;
// 			}
// 		}
// 	}

// 	// Create grid.
// 	for (i = 0; i <= Nx; ++i) {

// 		L = vertices[i][0].x;

// 		// y vertices for section 1
// 		if (L <= l1) {
// 			dy = inlet_height / (Ny + 1);
// 			for (j = 0; j <= Ny; ++j) {
// 				vertices[i][j].y = j * dy;
// 			}
// 		}

// 		// y vertices for section 2
// 		else if ((L > l1) && (L <= (l1 + l2))) {

// 			// Changes dy based on section
// 			if (L >= l1 && L < l1 + l2 + dx)  dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta1 * deg_to_rads);

// 			if (L == l1 + l2) dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta1 * deg_to_rads);

// 			vertices[i][0].y = vertices[i - 1][0].y + dy_ramp;
// 			dy = (inlet_height - vertices[i][0].y) / (Ny + 1);

// 			for (j = 1; j <= Ny; ++j) {
// 				vertices[i][j].y = vertices[i][0].y + j * dy;
// 			}
// 		}

// 		// y vertices for section 3
// 		else if ((L > l1 + l2) && (L <= L_total - l4)) {

// 			// Changes dy based on section
// 			if (L >= l2 && L < L_total)  dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta2 * deg_to_rads);

// 			if (L == L_total) dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta2 * deg_to_rads);

// 			vertices[i][0].y = vertices[i - 1][0].y + dy_ramp;
// 			dy = (inlet_height - vertices[i][0].y) / (Ny + 1);

// 			for (j = 1; j <= Ny; ++j) {
// 				vertices[i][j].y = vertices[i][0].y + j * dy;
// 			}
// 		}

// 		else {
// 			dy = (inlet_height - (l1 * tan(theta1 * deg_to_rads) + l2 * tan(theta2 * deg_to_rads))) / (Ny + 1);
// 			for (j = 0; j <= Ny; ++j) {
// 				vertices[i][j].y = l1 * tan(theta1 * deg_to_rads) + l2 * tan(theta2 * deg_to_rads) + j * dy;
// 			}
// 		}

// 	}

// 	// Edge vectors
// 	Point AB, BC, CD, DA;

// 	// Calculates cell centers and volumes
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {
// 			DA = { vertices[i][j].x - vertices[i][j + 1].x, vertices[i][j].y - vertices[i][j + 1].y };
// 			AB = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };
// 			BC = { vertices[i + 1][j + 1].x - vertices[i + 1][j].x, vertices[i + 1][j + 1].y - vertices[i + 1][j].y };
// 			CD = { vertices[i][j + 1].x - vertices[i + 1][j + 1].x, vertices[i][j + 1].y - vertices[i + 1][j + 1].y };

// 			cellCenters[i][j] = { (vertices[i][j].x + vertices[i + 1][j].x + vertices[i + 1][j + 1].x + vertices[i][j + 1].x) / 4,
// 									(vertices[i][j].y + vertices[i + 1][j].y + vertices[i + 1][j + 1].y + vertices[i][j + 1].y) / 4 };

// 			cellVolumes[i][j] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);
// 		}
// 	}

// 	// Calculates geometries for i-faces
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {

// 			AB = { vertices[i][j + 1].x - vertices[i][j].x, vertices[i][j + 1].y - vertices[i][j].y };

// 			iAreas[i][j] = sqrt(AB.x * AB.x + AB.y * AB.y);

// 			iNormals[i][j].x = AB.y / fabs(iAreas[i][j]);
// 			iNormals[i][j].y = AB.x / fabs(iAreas[i][j]);
// 		}
// 	}

// 	// Calculates geometries for j-faces
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {

// 			CD = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };

// 			jAreas[i][j] = sqrt(CD.x * CD.x + CD.y * CD.y);

// 			jNormals[i][j].x = -CD.y / fabs(jAreas[i][j]);
// 			jNormals[i][j].y = CD.x / fabs(jAreas[i][j]);
// 		}
// 	}

// }

// Point DoubleConeGrid::DoubleConeGrid::Center(int i, int j) const {
// 	return cellCenters[i][j];
// }

// Point DoubleConeGrid::DoubleConeGrid::Vertex(int i, int j) const {
// 	return vertices[i][j];
// }

// double DoubleConeGrid::DoubleConeGrid::Volume(int i, int j) const {
// 	return cellVolumes[i][j];
// }

// double DoubleConeGrid::DoubleConeGrid::iArea(int i, int j) const {
// 	return iAreas[i][j];
// }

// double DoubleConeGrid::DoubleConeGrid::jArea(int i, int j) const {
// 	return jAreas[i][j];
// }


// Point DoubleConeGrid::DoubleConeGrid::iNorms(int i, int j) const {
// 	return iNormals[i][j];
// }
// Point DoubleConeGrid::DoubleConeGrid::jNorms(int i, int j) const {
// 	return jNormals[i][j];
// }


// ///////////////////////////////////////////////////
// /////////////// Mirrored Grid functions ////////
// ///////////////////////////////////////////////////

// MirroredGrid::MirroredGrid(int Nx, int Ny, double l1, double l2, double l3, double l4, double theta1, double theta2, double inlet_height) : Nx(Nx), Ny(Ny),
// l1(l1), l2(l2), l3(l3), l4(l4), theta1(theta1), theta2(theta2), inlet_height(inlet_height), vertices(Nx + 1, vector<Point>(Ny + 1)), cellCenters(Nx, vector<Point>(Ny)),
// cellVolumes(Nx, vector<double>(Ny)), iAreas(Nx + 1, vector<double>(Ny, 0.0)), jAreas(Nx, vector<double>(Ny + 1, 0.0)), iNormals(Nx + 1, vector<Point>(Ny)), jNormals(Nx, vector<Point>(Ny + 1)) {

// 	int i, j;

// 	// Important constants
// 	double deg_to_rads = 3.141592653 / 180.0;

// 	// Define lengths 
// 	double L_total = l1 + l2 + l3 + l4;
// 	double dx = L_total / (Nx + 1);
// 	double dy = inlet_height / (Ny + 1);
// 	double L, dy_ramp;

// 	// Create x vertices
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {
// 			vertices[i][j].x = i * dx;
// 		}
// 	}

// 	// Snap x-vertices to important boundary points.
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {
// 			if (vertices[i][j].x > l1 && vertices[i][j].x < l1 + dx) {
// 				vertices[i][j].x = l1;
// 			}
// 			else if (vertices[i][j].x > l1 + l2 && vertices[i][j].x < l1 + l2 + dx) {
// 				vertices[i][j].x = l1 + l2;
// 			}
// 			else if (i == Nx) {
// 				vertices[i][j].x = L_total;
// 			}
// 		}
// 	}

// 	// Create grid.
// 	for (i = 0; i <= Nx; ++i) {

// 		L = vertices[i][0].x;

// 		// y vertices for section 1
// 		if (L <= l1) {
// 			dy = 2 * inlet_height / (Ny + 1);
// 			for (j = 0; j <= Ny; ++j) {
// 				vertices[i][j].y = j * dy;
// 			}
// 		}

// 		// y vertices for section 2
// 		else if ((L > l1) && (L <= (l1 + l2))) {

// 			// Changes dy based on section
// 			if (L >= l1 && L < l1 + l2 + dx)  dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta1 * deg_to_rads);

// 			if (L == l1 + l2) dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta1 * deg_to_rads);

// 			vertices[i][0].y = vertices[i - 1][0].y + dy_ramp;
// 			dy = (2 * (inlet_height - vertices[i][0].y)) / (Ny + 1);

// 			for (j = 1; j <= Ny; ++j) {
// 				vertices[i][j].y = vertices[i][0].y + j * dy;
// 			}
// 		}

// 		// y vertices for section 3
// 		else if ((L > l1 + l2) && (L <= L_total - l4)) {

// 			// Changes dy based on section
// 			if (L >= l2 && L < L_total)  dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta2 * deg_to_rads);

// 			if (L == L_total) dy_ramp = (vertices[i][0].x - vertices[i - 1][0].x) * tan(theta2 * deg_to_rads);

// 			vertices[i][0].y = vertices[i - 1][0].y + dy_ramp;
// 			dy = (2 * (inlet_height - vertices[i][0].y)) / (Ny + 1);

// 			for (j = 1; j <= Ny; ++j) {
// 				vertices[i][j].y = vertices[i][0].y + j * dy;
// 			}
// 		}

// 		else {
// 			dy = (2 * (inlet_height - (l1 * tan(theta1 * deg_to_rads) + l2 * tan(theta2 * deg_to_rads)))) / (Ny + 1);
// 			for (j = 0; j <= Ny; ++j) {
// 				vertices[i][j].y = l1 * tan(theta1 * deg_to_rads) + l2 * tan(theta2 * deg_to_rads) + j * dy;
// 			}
// 		}

// 	}

// 	// Edge vectors
// 	Point AB, BC, CD, DA;

// 	// Calculates cell centers and volumes
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {
// 			DA = { vertices[i][j].x - vertices[i][j + 1].x, vertices[i][j].y - vertices[i][j + 1].y };
// 			AB = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };
// 			BC = { vertices[i + 1][j + 1].x - vertices[i + 1][j].x, vertices[i + 1][j + 1].y - vertices[i + 1][j].y };
// 			CD = { vertices[i][j + 1].x - vertices[i + 1][j + 1].x, vertices[i][j + 1].y - vertices[i + 1][j + 1].y };

// 			cellCenters[i][j] = { (vertices[i][j].x + vertices[i + 1][j].x + vertices[i + 1][j + 1].x + vertices[i][j + 1].x) / 4,
// 									(vertices[i][j].y + vertices[i + 1][j].y + vertices[i + 1][j + 1].y + vertices[i][j + 1].y) / 4 };

// 			cellVolumes[i][j] = 0.5 * fabs(DA.x * AB.y - DA.y * AB.x) + 0.5 * fabs(BC.x * CD.y - BC.y * CD.x);
// 		}
// 	}

// 	// Calculates geometries for i-faces
// 	for (i = 0; i <= Nx; ++i) {
// 		for (j = 0; j < Ny; ++j) {

// 			AB = { vertices[i][j + 1].x - vertices[i][j].x, vertices[i][j + 1].y - vertices[i][j].y };

// 			iAreas[i][j] = sqrt(AB.x * AB.x + AB.y * AB.y);

// 			iNormals[i][j].x = AB.y / fabs(iAreas[i][j]);
// 			iNormals[i][j].y = AB.x / fabs(iAreas[i][j]);
// 		}
// 	}

// 	// Calculates geometries for j-faces
// 	for (i = 0; i < Nx; ++i) {
// 		for (j = 0; j <= Ny; ++j) {

// 			CD = { vertices[i + 1][j].x - vertices[i][j].x, vertices[i + 1][j].y - vertices[i][j].y };

// 			jAreas[i][j] = sqrt(CD.x * CD.x + CD.y * CD.y);

// 			jNormals[i][j].x = -CD.y / fabs(jAreas[i][j]);
// 			jNormals[i][j].y = CD.x / fabs(jAreas[i][j]);
// 		}
// 	}

// }

// Point MirroredGrid::MirroredGrid::Center(int i, int j) const {
// 	return cellCenters[i][j];
// }

// Point MirroredGrid::MirroredGrid::Vertex(int i, int j) const {
// 	return vertices[i][j];
// }

// double MirroredGrid::MirroredGrid::Volume(int i, int j) const {
// 	return cellVolumes[i][j];
// }

// double MirroredGrid::MirroredGrid::iArea(int i, int j) const {
// 	return iAreas[i][j];
// }

// double MirroredGrid::MirroredGrid::jArea(int i, int j) const {
// 	return jAreas[i][j];
// }


// Point MirroredGrid::MirroredGrid::iNorms(int i, int j) const {
// 	return iNormals[i][j];
// }
// Point MirroredGrid::MirroredGrid::jNorms(int i, int j) const {
// 	return jNormals[i][j];
// }