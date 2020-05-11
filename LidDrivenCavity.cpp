#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include "cblas.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <mpi.h>

using namespace std;
//Used in the Output_Results Function
#define FORMAT(vOut, x, y, v, s, Vx, Vy) vOut << setw(10) << x << setw(10) << y << setw(12) << v << setw(12) << s << setw(12) << Vx << setw(12) << Vy << endl;

///////////////////////////////////////////////////////////////////
//CONSTRUCTOR
///////////////////////////////////////////////////////////////////
LidDrivenCavity::LidDrivenCavity(MPI_Comm My_grid, int rank, int *coordinates, double *start, int *Neighbour, int nx, int ny, double deltat, double finalt, double re, double dx, double dy, bool &dt_stability)
{   this-> My_grid = My_grid;
    this-> rank = rank;
    memcpy(this->coordinates, coordinates, 2*sizeof(int));
    memcpy(this->Neighbour, Neighbour, 4*sizeof(int));
    memcpy(this->start, start, 2*sizeof(double));
    Nx = nx;
    Ny = ny;
    size = Nx * Ny;
    leading_a = 1+Ny;
    leading_b = 3;
    leading_c = 1+2*Ny;
    T = finalt;
    Re = re;
    this->dx = dx;
    this->dy = dy;
    dt_stability = SetTimeStep(deltat);	// determine the stability restriction on dt
}

///////////////////////////////////////////////////////////////////
//DESTRUCTOR
///////////////////////////////////////////////////////////////////
LidDrivenCavity::~LidDrivenCavity()
{
	delete[] v;
	delete[] s;
	delete[] v_top;
	delete[] v_bottom;
	delete[] v_left;
	delete[] v_right;
	delete[] s_top;
	delete[] s_bottom;
	delete[] s_left;
	delete[] s_right;
	delete[] A;
	delete[] B;
	delete[] C;
}


///////////////////////////////////////////////////////////////////
//MEMBER FUNCTIONS
///////////////////////////////////////////////////////////////////
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;//
    Ly = ylen;//initialise array
}
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx=nx;//Given in section 2.
    Ny=ny;
}

bool LidDrivenCavity::SetTimeStep(double deltat)
{
    try
   {
	dt = deltat;
	/// check the restriction on time step dt < Re*dx*dy/4
	if (dt >= Re*dx*dy/4)
	{
        // error catching
		throw std::logic_error("Courant-Friedrichs-Lewy Condition has not been satisfied. USE different inputs. ");
	}
    }
   catch (const std::logic_error &e)
    {
	if (rank == 0){
    cout << "An error occured: " << e.what() << endl;

	return 0;
    }
    }
   return 1;
}
void LidDrivenCavity::SetFinalTime(double finalt)
{
     T= finalt;
}
void LidDrivenCavity::SetReynoldsNumber(double re)
{
     Re=re;
}
void LidDrivenCavity::Initialise()
{
    v=new double[size];
    s=new double[size];
    cblas_dscal(size, 0.0, v, 1);
    cblas_dscal(size, 0.0, s, 1);

    streamfunction_old = new double[size];
    Vx = new double[size];
    Vy = new double[size];
    // Allocate the memory for all 8 boundary vectors. Initialize all vectors to zero.
    v_top =  new double[Nx];
    v_bottom =  new double[Nx];
    v_left = new double[Ny];
    v_right =new double[Ny];

    s_top = new double[Nx];
    s_bottom = new double[Nx];
    s_left = new double[Ny];
    s_right = new double[Ny];
    cblas_dscal(Nx, 0.0, v_top, 1);
    cblas_dscal(Nx, 0.0, v_bottom, 1);
    cblas_dscal(Ny, 0.0, v_left, 1);
    cblas_dscal(Ny, 0.0, v_right, 1);
    cblas_dscal(Nx, 0.0, s_top, 1);
    cblas_dscal(Nx, 0.0, s_bottom, 1);
    cblas_dscal(Ny, 0.0, s_left, 1);
    cblas_dscal(Ny, 0.0, s_right, 1);
}
void LidDrivenCavity::Full_matrices()
{
	A = new double[size * leading_a];	                // store as symmetric banded matrix with bandwidth = Ny
	B = new double[size * leading_b];    	            // store as banded matrix with bandwidth = 1
	C = new double[size * leading_c];	                // store as banded matrix with bandwidth = Ny
	double A_MainDiag =2/(dx*dx) + 2/(dy*dy);	// diagonal entries
	double A_NearDiag =-1/(dy*dy);	            // entries above and below the diagonal
	double A_FarDiag = -1/(dx*dx);	            // entries in the off-diagonal blocks
	double B_NearDiag = 1/(2*dy);
	double C_FarDiag =  1/(2*dx);
	// Matrix A creation.
	for (unsigned int i = 0; i < size; i++)
	{
		A[leading_a-1 + i*leading_a] = A_MainDiag;
		if (i % Ny != 0)   A[leading_a-2 + i*leading_a] = A_NearDiag;
		if (i >= Ny)       A[i * leading_a] = A_FarDiag;
	}
	// Matrix B creation.
	for (unsigned int i = 0; i < size; i++)
	{
		if (i % Ny != 0)    B[i*leading_b] = B_NearDiag;
		if (i % Ny != Ny-1)  B[2 + i*leading_b] = -B_NearDiag;
	}
	/// Matrix C creation.
	for (unsigned int i = 0; i < size; i++)
	{
		if (i >= Ny)      C[i*leading_c] = C_FarDiag;
		if (i < size-Ny)  C[leading_c-1 + i*leading_c] = -C_FarDiag;
	}
}
// Function to output boundary vector b in the linear system y = Ax + b. Input matrix can be either A, B or C. 
// Arrays are stored in column major form. 
void LidDrivenCavity::Boundary_vector(double *b, char matrix, char x)
{
	double A_NearDiag = -1/(dy*dy);
	double A_FarDiag  = -1/(dx*dx);
	double B_NearDiag =  1/(2*dy);
	double C_FarDiag  =  1/(2*dx);
	if (matrix == 'A')
	{
		if (x == 's')//Streamfunction
		{
			for (unsigned int i = 0; i < Ny; i++)
			{
				b[i] += s_left[i] * A_FarDiag;
				b[i+(Nx-1)*Ny] += s_right[i] * A_FarDiag;
			}
			for (unsigned int i = 0; i < Nx; i++)
			{
				b[i*Ny] += s_bottom[i] * A_NearDiag;
				b[i*Ny + Ny-1] += s_top[i] * A_NearDiag;
			}
		}
		else if (x == 'v')//Vorticity
		{
			for (unsigned int i = 0; i < Ny; i++)
			{
				b[i] += v_left[i] * A_FarDiag;
				b[i+(Nx-1)*Ny] += v_right[i] * A_FarDiag;
			}
			for (unsigned int i = 0; i < Nx; i++)
			{
				b[i*Ny] += v_bottom[i] * A_NearDiag;
				b[i*Ny + Ny-1] += v_top[i] * A_NearDiag;
			}
		}
	}
	else if (matrix == 'B')
	{
		if (x == 's')//Streamfunction
		{
			for (unsigned int i = 0; i < Nx; i++)
			{
				b[i*Ny] = s_bottom[i] * (-B_NearDiag);
				b[i*Ny + Ny-1] = s_top[i] * B_NearDiag;
			}
		}
		else if (x == 'v')//Vorticity
		{
			for (unsigned int i = 0; i < Nx; i++)
			{
				b[i*Ny] = v_bottom[i] * (-B_NearDiag);
				b[i*Ny + Ny-1] = v_top[i] * B_NearDiag;
			}
		}
	}
	else if (matrix =='C')
	{
		if (x == 's')//Streamfunction
		{
			for (unsigned int i = 0; i < Ny; i++)
			{
				b[i] = s_left[i] * (-C_FarDiag);
				b[i+(Nx-1)*Ny] = s_right[i] * C_FarDiag;
			}
		}
		else if (x == 'v')//Vorticity
		{
			for (unsigned int i = 0; i < Ny; i++)
			{
				b[i] = v_left[i] * (-C_FarDiag);
				b[i+(Nx-1)*Ny] = v_right[i] * C_FarDiag;
			}
		} 
	}

}

//Solves for the vorticity boundary conditions, eqn 6 to 9.
//Avoid using cmath functions to reduce computation time
void LidDrivenCavity:: BoundaryConditions(){
        // Top boundary
    if (Neighbour[0] == -2)
    {
        for (unsigned int i = 0; i < Nx; i++)
        {
            v_top[i]=(s_top[i]-s[i*Ny+(Ny-1)])*2/(dy*dy)-2*U/dy;//Equation 6
        }
    }
    // Bottom wall
    if (Neighbour[2] == -2)
    {
        for (unsigned int i = 0; i < Nx; i++)
        {
            v_bottom[i]=(s_bottom[i]-s[i*Ny])*2/(dy*dy);//Equation 7
        }
    }
    // Left wall
    if (Neighbour[1] == -2)
    {
        for (unsigned int i = 0; i < Ny; i++)
        {
            v_left[i] =(s_left[i]-s[i])*2/(dx*dx); //Equation 8
        }
    }
    // Right wall
    if (Neighbour[3] == -2)
    {
        for (unsigned int i = 0; i < Ny; i++)
        {
            v_right[i] =(s_right[i]-s[(Nx-1)*Ny+i])*2/(dx*dx);//Equation 9
        }
    }

}
//Solves for the interior vorticity current time step, eqn 10.
void LidDrivenCavity:: InteriorVor_T(){
    for (unsigned int j = 0; j < Ny - 2; ++j) {            // index rows
		for (unsigned int i = 0; i < Nx - 2; ++i) {        // index cols
			v[Ny+1+(i*Ny)+j]=(-1)*((s[Ny+1+((i+1)*Ny)+j]-2*s[Ny+1+(i*Ny)+j]+s[Ny+1+((i-1)*Ny)+j])/(dx*dx)+(s[Ny+(i*Ny)+j]-2*s[Ny+1+(i*Ny)+j]+s[Ny+2+(i*Ny)+j])/(dy*dy));
		}
	}
}
//Solves for the interior vorticity at next time step, t+dt.  eqn 11. We will use cblas routines to minimise run time and produce a cleaner code. 
void LidDrivenCavity:: InteriorVor_delT()
{
	  	// Vector, derived from eqn 11.
	double *b1 = new double[size];	// viscosity term, Term on the RHS
	double *b2 = new double[size];	// 3rd product on the LHS of eqn 11
	double *b3 = new double[size];	// 2nd product on the LHS of eqn 11
	cblas_dscal(size, 0.0, b1, 1);
	cblas_dscal(size, 0.0, b2, 1);
	cblas_dscal(size, 0.0, b3, 1);

	// First calculate b1 = Av + b. A matrix is symmetrical
	Boundary_vector(b1, 'A', 'v');
	cblas_dsbmv (CblasColMajor, CblasUpper, size, Ny, 1.0, A, leading_a, v, 1, 1.0, b1, 1); // A is a symmetric banded matrix with bandwidth= Ny. Solves the equation: b = Av + b. Overwrites and outputs to b1.

	// Calculate b2
	double *temp1 = new double[size];	// temp1 = Cs + b, b refers to the boundary term
	double *temp2 = new double[size];	// temp2 = Bv + b
	cblas_dscal(size, 0.0,temp1,1);
	cblas_dscal(size, 0.0,temp2,1);

	// calculate temp1, C and B matrix are general banded.
	Boundary_vector(temp1, 'C', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, leading_c, s, 1, 1.0, temp1, 1);
	// calculate temp2
	Boundary_vector(temp2, 'B', 'v');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, leading_b, v, 1, 1.0, temp2, 1);
	//Multiply each vector element. Produces b3, the 2nd product on the LHS.
	for (unsigned int i = 0; i < size; i++)
	{
		b2[i] = temp1[i] * temp2[i];
	}
	// Frees up memory from temp1 and temp2
	delete[] temp1;
	delete[] temp2;

	// Calculate b3
	double *temp3 = new double[size];	///< temp3 = Cv + b, b is the boundary term
	double *temp4 = new double[size];	///< temp4 = Bs + b
	cblas_dscal(size,0.0,temp3,1);
	cblas_dscal(size,0.0,temp4,1);
	// calculate temp3
	Boundary_vector(temp3, 'C', 'v');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, leading_c, v, 1, 1.0, temp3, 1);
	// calculate temp4
	Boundary_vector(temp4, 'B', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, leading_b, s, 1, 1.0, temp4, 1);
	// Multiply each vector element. Produces b3, the 2nd product on the LHS.
	for (unsigned int i = 0; i < size; i++)
	{
		b3[i] = temp3[i] * temp4[i];
	}
	// release memory for temp3 and temp4
	delete[] temp3;
	delete[] temp4;

	// Update interior vorticity, i.e. vorticity = vorticity -(dt/Re)*b1 + dt*b2 - dt*b3. Rearrangement of eqn 11.
	cblas_daxpy (size, -dt/Re, b1, 1, v, 1);
	cblas_daxpy (size,  dt,    b2, 1, v, 1);
	cblas_daxpy (size, -dt,    b3, 1, v, 1);

	// release memory
	delete[] b1;
	delete[] b2;
	delete[] b3;
}
//Solves for the poisson problem, eqn 12.
void LidDrivenCavity::Integrate(){
    // Create a new PoissonSolver instance and solve Poisson problem
	if (Poisson_Solver == nullptr)//Makes sure this is the first time that the instance is created, helps with optimisation. 
	{
		Poisson_Solver = new PoissonSolver(Nx, Ny, dx, dy);
	}
	//Poisson_Solver->FeedData(Lx, Ly, Nx, Ny, dt);
	Poisson_Solver->SetBoundary(s_top, s_left, s_bottom, s_right);
	Poisson_Solver->Solve_Problem(s, v);
}
//Function to send and receive inputs and outputs between processors. Calls MPI functions. 
void LidDrivenCavity::MPI_Interaction(double *x, double *x_top, double *x_left, double *x_bottom, double *x_right)
{
    // Left to right
    MPI_Sendrecv(x + (Nx-1)*Ny, Ny, MPI_DOUBLE, Neighbour[3], 3, x_left, Ny, MPI_DOUBLE, Neighbour[1], 3, My_grid, MPI_STATUS_IGNORE);

    // Right to left
    MPI_Sendrecv(x, Ny, MPI_DOUBLE, Neighbour[1], 1, x_right, Ny, MPI_DOUBLE, Neighbour[3], 1, My_grid, MPI_STATUS_IGNORE);

    // Top to bottom
    // extract temporary bottom boundary vector
    double temp_bottom[Nx];
    cblas_dcopy(Nx, x, Ny, temp_bottom, 1);

    MPI_Sendrecv(temp_bottom, Nx, MPI_DOUBLE, Neighbour[2], 2, x_top, Nx, MPI_DOUBLE, Neighbour[0], 2, My_grid, MPI_STATUS_IGNORE);

    // Bottom to top
    // extract temporary top boundary vector
    double temp_top[Nx];
    cblas_dcopy(Nx, x + Ny-1, Ny, temp_top, 1);
    MPI_Sendrecv(temp_top, Nx, MPI_DOUBLE, Neighbour[0], 0, x_bottom, Nx, MPI_DOUBLE, Neighbour[2], 0, My_grid, MPI_STATUS_IGNORE);
}
 //Get velocity from stream function/vorticity. Relation is achieved through equation 3 of assignment.
void LidDrivenCavity::Calculate_Velocity()
{
	//Need to solve the equation ds/dy=u. To do this, express ds/dy as a finite difference using matrix B, and solve the corresponding matrix problem. Output to Vx!
	Boundary_vector(Vx, 'B', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, leading_b, s, 1, 1.0, Vx, 1);
	//Need to solve the equation ds/dx=-v. To do this, express ds/dx as a finite difference using matrix C, and solve the corresponding matrix problem. Output to Vy!
	Boundary_vector(Vy, 'C', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, leading_c, s, 1, 1.0, Vy, 1);
}
//Function to run through entire algorithm over the entire specified time period, t=0 to t=T.
//Measure of accuracy of our solution will be through the evolution of the relative residual error. I.e. A low residual error means that the accuracy of the solution at the time step is good. 
void LidDrivenCavity::Solve_full(){

    double t=0.0;
    double tolerance = 0.00001; //Terminate criteria. This function will run until either the entire time domain is completed
                          //or when a steady state solution is computed.
    double stream_func_2norm= 0.0;	// 2-norm of stream function
    double residual = 0.0;	// residual, take 2-norm of stream function
    double max_residual= 0.0;	// maximum residual along the whole domain

	SetDomainSize(Lx, Ly);
    //SetGridSize(Nx, Ny);
    SetTimeStep(dt);
    SetFinalTime(T);
    SetReynoldsNumber(Re);

	Full_matrices();
    Initialise();
    while (t<=T)
	{
        cblas_dcopy(size, s, 1, streamfunction_old, 1);	///< store streamfunction of the previous step
        BoundaryConditions();
        InteriorVor_T();
        MPI_Interaction(v, v_top, v_left, v_bottom, v_right);
        InteriorVor_delT();
        //cout<<"vorticity matrix at time "<< t <<endl;
        //PrintVorMatrix();
        for (unsigned int i=0; i<5;i++)
		{//Iteration for all 5 arrays.
            Integrate();
            MPI_Interaction(s, s_top, s_left, s_bottom, s_right);
        }
		t+=dt;
        stream_func_2norm = cblas_dnrm2(size, s, 1);
        cblas_daxpy(size, -1.0, s, 1, streamfunction_old, 1);
        residual = cblas_dnrm2(size, streamfunction_old, 1); //calculation of 2-norm with cblas routine
        residual /= stream_func_2norm;

        MPI_Reduce(&residual, &max_residual, 1, MPI_DOUBLE, MPI_MAX, 0, My_grid);//Outputs max value to max_residual.
        if (rank == 0)
        {
                cout << "t = " << setw(5) << t << setw(20) << "MAX_RESIDUAL: " << max_residual << endl;
        }
        MPI_Bcast(&max_residual, 1, MPI_DOUBLE, 0, My_grid);
        if(max_residual < tolerance)
        {
            if (rank == 0){
				cout << "TERMINATE CRITERIA ACHIEVED, EXIT SUCCESSFUL." << endl;

			}
			break;
		}
    }
    Calculate_Velocity();
}


//Outputs to Output.txt, all data to be used to plot on Matlab.
void LidDrivenCavity::Output_Results(int Px, int Py, double Lx, double Ly)
{
	for (int k = 0; k < Px *Py; k++)
	{
		if (k == rank)
		{
			// Write header and corner values. Corner values are known to be 0 all the time.
			if (k == 0)
			{
				ofstream vOut("Output.txt", ios::out | ios::trunc);
				FORMAT(vOut, "x", "y", "vort", "streamfunc", "Vx", "Vy");
				FORMAT(vOut, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
				FORMAT(vOut, 0.0, Ly, 0.0, 0.0, 0.0, 0.0);
				FORMAT(vOut, Lx, 0.0, 0.0, 0.0, 0.0, 0.0);
				FORMAT(vOut, Lx, Ly, 0.0, 0.0, 0.0, 0.0);
				vOut.close();
			}
			ofstream vOut("Output.txt", ios::out | ios::app);
			double x, y;
			// Write interior values of vorticity, streamfunction, Vx and Vy
			vOut.precision(5);
			for (unsigned int i = 0; i < Nx; i++)
			{
				for (unsigned int j = 0; j < Ny; j++)
				{
					x = start[0] + i*dx;
					y = start[1] + j*dy;
					FORMAT(vOut, x, y, v[i*Ny + j], s[i*Ny + j],
						       	Vx[i*Ny + j], -Vy[i*Ny + j]); //In order to be consistent with the given sign convention, we need to enforce -Vy.
				}
			}
			//__________________________________________________________________________________________________________//
			// Write boundary values
			// Top. Non zero Vx value, Vx=U!
			if (coordinates[0] == Py - 1)
			{
				y = Ly;
				for (unsigned int i = 0; i < Nx; i++)
				{
					x = start[0] + i*dx;
					FORMAT(vOut, x, y, v_top[i], 0.0, U, 0.0);
				}
			}
			// Bottom
			if (coordinates[0] == 0)
			{
				y = 0;
				for (unsigned int i = 0; i < Nx; i++)
				{
					x = start[0] + i*dx;
					FORMAT(vOut, x, y, v_bottom[i], 0.0, 0.0, 0.0);
				}
			}
			// left
			if (coordinates[1] == 0)
			{
				x = 0;
				for (unsigned int j = 0; j < Ny; j++)
				{
					y = start[1] + j*dy;
					FORMAT(vOut, x, y, v_left[j], 0.0, 0.0, 0.0);
				}
			}
			// right
			if (coordinates[1] == Px - 1)
			{
				x = Lx;
				for (unsigned int j = 0; j < Ny; j++)
				{
					y = start[1] + j*dy;
					FORMAT(vOut, x, y, v_right[j], 0.0, 0.0, 0.0);
				}
			}
			vOut.close();
		}
		MPI_Barrier(My_grid);
	}

	cout<<"OUTPUT COMPLETE!" <<endl;//STATEMENT TO SHOW COMPLETE OPERATION
}

//Used to print matrix, for checking and displaying solution.
void LidDrivenCavity::PrintVorMatrix() {
    cout.precision(4);
    for (unsigned int i = 0; i < Nx; ++i) {
        for (unsigned int j = 0; j < Ny; ++j) {
            cout << setw(6) << v[j*Nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
//Used to print matrix, for checking and displaying solution.
void LidDrivenCavity::PrintStrMatrix() {
    cout.precision(4);
    for (unsigned int i = 0; i < Nx; ++i) {
        for (unsigned int j = 0; j <Ny; ++j) {
            cout << setw(6) << s[j*Nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
