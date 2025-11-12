# MPI Structured Inviscid Implicit CFD Solver

## Important libraries

- [linalglib](./linalglib/): this library contains all linear algebra related tools like outer products, matrix multiplication, etc. The biggest thing in here is `typedef vector<double> Vector` which allows me to write `Vector A` instead of `vector<double> A` when defining a vector.

- [gibbslib](./gibbslib/): this is my old chemical equilibrium library that you can ignore

- [cfdsolver](./cfdsolver/): this contains the main file for cfd. Set initial conditions / one of two grids I have defined in my solver class, and then run the constructor / initialize MPI, and then solve.

- [solverlib](./solverlib/): this contains everything for the inviscid CFD solver. Also contains a 1D Sod solver.

- [explicit_sod](./explicit_sod/): This contains a serial version and an MPI version of the 1D Sod problem named "serial.cpp" and "parallel.cpp", respectively. No classes are used, its an all in one code that only relies on the linear algebra library. It solves the simulation explicitly. 

- [programs](./programs/): contains all of the built executables. To run the cfdsolver just type `mpirun -np <N> cfdsolver` from this directory. To run parallel Sod solver, `mpirun -np <N> sod`. To run serial Sod solver type `./ssod`.

- [plots](./plots/)], [plotfiles](./plotfiles/), [testing](./testing/), and [build](./build) are all for organizing non-code related things. You can ignore them.
