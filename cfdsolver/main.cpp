#include "../solverlib/solver.hpp"


int main(int argc, char* argv[]) {
    
    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime(); 

    int Nx = 50; 

    SodSolver1D solver(Nx);

    solver.solve();  

    double end_time = MPI_Wtime();    // End timer
    double elapsed = end_time - start_time;

    if (solver.rank == 0) {
        std::cout << "Elapsed time: " << elapsed << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}
