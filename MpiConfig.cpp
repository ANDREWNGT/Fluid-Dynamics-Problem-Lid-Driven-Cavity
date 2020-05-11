#include <iostream>
#include <mpi.h>

using namespace std;
#include "MpiConfig.h"
//cpp file that contains functions to control MPI functionalities. Critical to ensuring processor interdependencies are catered for. 

///////////////////////////////////////////////////////////////////
//MEMBER FUNCTIONS
///////////////////////////////////////////////////////////////////

//Function to find neighbourhood in all 4 directions for current rank. 
void FindNeighbour(MPI_Comm mygrid, int *Neighbour)
{
	MPI_Cart_shift(mygrid, 0, -1, Neighbour ,    Neighbour + 2);	// top & bottom
	MPI_Cart_shift(mygrid, 1,  1, Neighbour + 1, Neighbour + 3);	// left & right
}



//Function to distribute work as evenly as possible to each processor. 
void DistributeWork(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py,
            const double &dx, const double &dy, int *coords, double *start, int &nx, int &ny)
{
    // ranks with smaller coordinates have 1 larger size in grids
    nx = (Nx-2) / Px;
    ny = (Ny-2) / Py;
    int x_remainder = (Nx-2) % Px;
    int y_remainder = (Ny-2) % Py;

    // Caution: start[i, j] and coordinates[j, i]. This is due to the axes defined in the problem.  
    // x direction 
    if (coordinates[1] < x_remainder){
        nx++;
        start[0] = (coordinates[1]*nx + 1)*dx;
    }
    else{
        start[0] = (Nx-1)*dx - (Px - coordinates[1])*dx*nx;
    }
    // y direction
    if (coordinates[0] < y_remainder){
        ny++;
        start[1] = (coordinates[0]*ny + 1)*dy;
    }
    else
    {
        start[1] = (Ny-1)*dy - (Py - coordinates[0])*dy*ny;
    }
    
    // Status report before commencing work. 
    cout<<"Rank: "<<rank<<", nx = "<< nx <<", ny = "<< ny << "; Commence at global coordinates(" << start[0] << ", " << start[1] <<"). Cartesian coordinates ["<<coordinates[0] <<", " << coordinates[1] << "]." << endl; 
}


 // Function to check that no. of processors is compatible with no. of partitions.
bool CheckWorkers(const int &np, const int &Px, const int &Py)
{
	if (np == Px * Py) {
        return true;
    }
	else{
        return false;
    }
}