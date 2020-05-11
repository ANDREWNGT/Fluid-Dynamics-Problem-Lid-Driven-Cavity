#ifndef LID_DRIVEN_CAVITY
#define LID_DRIVEN_CAVITY

#include <string>
using namespace std;

#include "PoissonSolver.h"
#include <mpi.h>
class LidDrivenCavity
{
public:
    LidDrivenCavity(MPI_Comm My_grid, int rank, int *coordinates, double *start, int *Neighbour, int nx,
	       	int ny, double deltat, double finalt, double re, double dx, double dy, bool &dt_stability);
    ~LidDrivenCavity();
    
    //Setter functions
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    bool SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    
    //Methods employed by solver
    void Initialise();
    void PrintVorMatrix();
    void PrintStrMatrix();
    void BoundaryConditions();
    void InteriorVor_T();
    void InteriorVor_delT(); //After timestep
    void Boundary_vector(double *b, char matrix, char x);
    void Integrate();
    void Full_matrices();
    // Add any other public functions
    void Solve_full();
    //MPI Functions
    void MPI_Interaction(double *x, double *x_top, double *x_left,
            double *x_bottom, double *x_right);
    void Calculate_Velocity();
    //Methods for output 
    void Output_Results(int Px, int Py, double Lx, double Ly);

private:

    //MPI configuration 
    MPI_Comm My_grid;
    int rank;
    int coordinates[2];       // coordinate in Cartesian topology
    double start[2];     // start global coordinate
    int Neighbour[4];    // ranks of neighborhood

    double* v = nullptr;   //vorticity. array
    double* s = nullptr;   //streamfunction. array
    double *v_top = nullptr;
    double *v_left = nullptr;
    double *v_bottom = nullptr;
    double *v_right = nullptr;
    double *s_top = nullptr;
    double *s_left = nullptr;
    double *s_bottom = nullptr;
    double *s_right = nullptr;
    double *streamfunction_old = nullptr; // store previous streamfunction each step to obtain residual
    double *Vx =nullptr;	 // horizontal velocity 
    double *Vy =nullptr;	 // vertical velocity
    

    double dt;
    double T;
    unsigned int    Nx;
    unsigned int    Ny;
    double Lx;
    double Ly;
    double Re;
    double dy;
    double dx;
    const double U= 1.0; //Lid velocity is set as 1.0, given in question. 
    
    double *A = nullptr;	
    double *B = nullptr;    
    double *C = nullptr;    
    unsigned int size;
    unsigned int leading_a; 
	unsigned int leading_b; 
	unsigned int leading_c;
	
	PoissonSolver *Poisson_Solver = nullptr;	///< poisson solver object

};

#endif