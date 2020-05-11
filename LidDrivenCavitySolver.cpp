#include <iostream>
#include <exception>
#include <iomanip>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono; 
  
#include "LidDrivenCavity.h"
#include "ProgramOptions.h"
#include "MpiConfig.h"
#include <mpi.h>

 // C++ Code (in parallel) to solve General 2D lid-driven cavity using
 // an Explicit Forward time-integration scheme. 
 // Arrays are stored in column-major format.
 // Written by Andrew Ng - 22/03/2020

int main(int argc, char **argv)
{
	// Get starting timepoint 
    auto begin = high_resolution_clock::now(); 
	
    int rank;
	int size;
    // Initialise MPI.
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get number of processors.
    //cout<<"RANK:" <<rank<<endl;
    //cout<<"SIZE:" <<size<<endl;
//Parameters from command input. Variables will be parsed through ProgramOptions.cpp file, called by FUNCTION ReadVals(...) below. 
	double Lx; 
	double Ly;
	int Nx;
	int Ny;
	int Px;
	int Py;
	double dt;
	double T;
	double Re;
    po::variables_map vm;
    bool status;	// to decide the input status, 1 for successful input
    status = OptionStatus(argc, argv, vm);
    /// Terminate the program if error occurs or help is called. 
    if (!status)
    {
        MPI_Finalize();
        return 0;
    }

    ReadVals(vm, Lx, Ly, Nx, Ny, Px, Py, dt, T, Re);// Reads in variables. 
    double dx = Lx/(Nx - 1);
    double dy = Ly/(Ny - 1);
    // Check if number of processors are compatible. KEY error check in order to avoid erronenous outputs. 
    if (!CheckWorkers(size, Px, Py))
    {
        if (rank == 0)
        {
            cout << "Number of processors does not match Px and Py." << endl;
            cout << "np = Px * Py must be fulfilled. Try again. " << endl;
        }
        MPI_Finalize();
        return 0;
    }
    		
    
    // Work before operating
    MPI_Comm My_grid;
    const int dims = 2;
    int sizes[dims] = {Py, Px};
    int periods[dims] = {0, 0};
    int reorder = 0;
    
    // New communicator based on Cartesian Topology
    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &My_grid);
    int coords[dims];
    /// Assign a rank to each coordinate
    MPI_Cart_coords(My_grid, rank, dims, coords);
    int Neighbour[4] = {0};
    /// Obtain neighborhood information for each rank
    FindNeighbour(My_grid, Neighbour);
    int nx; 
	int ny;	// number of grids for each rank
    double start[2] = {0.0, 0.0};
    // Distribute work to each process, calling function from MpiConfig.cpp
    DistributeWork(rank, Nx, Ny, Px, Py, dx, dy, coords, start, nx, ny);
    
    bool dt_stability;
    // Create a new instance of the LidDrivenCavity class for every rank. 
    LidDrivenCavity* solver = new LidDrivenCavity(My_grid, rank, coords, start, Neighbour,
            nx, ny, dt, T, Re, dx, dy, dt_stability);
    
    if (!dt_stability)//Condition to check for in order to enforce stability of solution. 
    {
        MPI_Finalize();
        return 0;
    }
    // Run the solver
    solver->Solve_full();

	cout<<"SOLVER CONCLUDED SUCCESSFULLY!"<<endl; //Statement to confirm successful operation. 
    // Output the result. Write to output.txt file.
    solver->Output_Results(Px, Py, Lx, Ly);
	cout<<"OUTPUT COMPLETE!"<<endl;
    // Exit
    MPI_Finalize();
	  // Get ending timepoint 
    auto stop = high_resolution_clock::now(); 
  
	  // Get duration. Substart timepoints to  
    // get durarion. To cast it to proper unit 
    // use duration cast method 
    auto duration = duration_cast<seconds>(stop - begin); 
  
    cout << "Time taken by code: "
         << duration.count() << " seconds" << endl; 
  
	return 0;
    
}




