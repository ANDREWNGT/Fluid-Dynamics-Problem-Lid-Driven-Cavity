#ifndef MPI_CONFIG
#define MPI_CONFIG

#include <mpi.h>
//Header file controlling all declarations of functions used to control MPI functionalities. 
void FindNeighbour(MPI_Comm mygrid, int *Neighbour);
void DistributeWork(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py,
	       	const double &dx, const double &dy, int *coordinates, double *start, int &nx, int &ny);
bool CheckWorkers(const int &np, const int &Px, const int &Py);

#endif