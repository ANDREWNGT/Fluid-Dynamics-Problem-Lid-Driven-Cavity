#include<cstring>
#include<iostream>
#include <iomanip>
#include<cmath>

using namespace std;

#include "PoissonSolver.h"
#include "cblas.h"

//Class that solves the Poisson Problem. Using cblas and lapack routines to achieve this. All arrays are arranged in column major format. 
#define F77NAME(x) x##_
extern "C"
{
	void F77NAME(dpbtrf) (const char &UPLO, const int &N, const int &KD, const double *AB,
		       const int &LDAB,	int &INFO);
	void F77NAME(dpbtrs) (const char &UPLO, const int &N, const int &KD, const int &NRHS,
		       	const double *AB, const int &LDAB, double *B, const int &LDB, int &INFO);
	void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, const int& nrhs, const double *AB, 
                        const int& ldab, int* ipiv, double * B,
                        const int& ldb, int& info);			
}

///////////////////////////////////////////////////////////////////
//CONSTRUCTOR
///////////////////////////////////////////////////////////////////
PoissonSolver::PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy)
{
	// accept the input parameter
	this->Nx = Nx;
	this->Ny = Ny;
	this->dx = dx;
	this->dy = dy;
	/// build matrix A of the system in banded and packed storage
	leading_a = Ny + 1;
	size = Nx * Ny;
	A = new double[size * leading_a];
	cblas_dscal(size*leading_a, 0.0, A, 1);	
	A_MainDiag = 2/(dx*dx) + 2/(dy*dy);	///< diagonal entries
	A_NearDiag = -1/(dy*dy);	///< entries above and below the diagonal 
	A_FarDiag = -1/(dx*dx);	///< entries in the off-diagonal block	
	for (unsigned int j = 0; j < size; j++)
	{
		A[leading_a-1 + j*leading_a] = A_MainDiag;
		if (j % Ny != 0)   {
			A[leading_a-2 + j*leading_a] = A_NearDiag;
		}
		if (j >= Ny){
			A[j * leading_a] = A_FarDiag;
		}
	}
	info = 1;
}
///////////////////////////////////////////////////////////////////
//DESTRUCTOR
///////////////////////////////////////////////////////////////////
PoissonSolver::~PoissonSolver()
{
	delete[] A;
	delete[] b;
}
///////////////////////////////////////////////////////////////////
//MEMBER FUNCTIONS
///////////////////////////////////////////////////////////////////
void PoissonSolver::FeedData(double xlen, double ylen, int nx, int ny, double deltat){
   double Lx = xlen;//
   double Ly = ylen;//initialise array
   Nx= nx;//Given in section 2. 
   Ny= ny;
   double dt=deltat;
   dx=Lx/(Nx-1);
   dy=Ly/(Ny-1);

   /// build matrix A of the system in banded and packed storage
	leading_a = Ny + 1;
	size = Nx * Ny;
	A = new double[size * leading_a];
	cblas_dscal(size*leading_a, 0.0, A, 1);	
	A_MainDiag = 2/(dx*dx) + 2/(dy*dy);	// diagonal entries
	A_NearDiag = -1/(dy*dy);	        // entries above and below the diagonal 
	A_FarDiag = -1/(dx*dx);	            // entries in the off-diagonal block	
	for (unsigned int j = 0; j < size; j++)
	{
		A[leading_a-1 + j*leading_a] = A_MainDiag;
		if (j % Ny != 0)   A[leading_a-2 + j*leading_a] = A_NearDiag;
		if (j >= Ny)       A[j * leading_a] = A_FarDiag;
	}
	info = 1;
}
//Accept boundary conditions from 4 walls and constructs boundary vector b.
void PoissonSolver::SetBoundary(const double *top, const double *left, const double *bottom, const double *right)
{	
	b = new double[size];
	cblas_dscal(size, 0.0, b, 1);
	for (unsigned int i = 0; i < Nx; i++)
	{
		b[i*Ny] += bottom[i] * A_NearDiag;        //b vector inputs for bottom wall 
		b[i*Ny + Ny-1] += top[i] * A_NearDiag; // b vector inputs for top wall
	}
	for (unsigned int i = 0; i < Ny; i++)
	{
		b[i] += left[i] * A_FarDiag;             //b vector inputs for left wall
		b[i+(Nx-1)*Ny] += right[i] * A_FarDiag;  //b vector inputs for right wall
	}
}
void PoissonSolver::Solve_Problem(double *x, const double *f)
{	//We will not use dgbmv as it is highly inefficient compared to using a combination of dpbtrf and dpbtrs.
	 // Check to see if Cholesky factorisation occurs only one time 
	if (info != 0)	
	{
		F77NAME(dpbtrf) ('U', size , Ny, A, leading_a, info); //A= U**T *U, output will overwrite A matrix. 
	}
	// Set RHS vector
	cblas_daxpy(size, -1.0, f, 1, b, 1);
	cblas_dscal(size, -1.0, b, 1);
	//// Solve Ax = RHS. Output is saved into b. 
	F77NAME(dpbtrs) ('U', size, Ny, 1, A, leading_a, b, size, info); //Solves A*X=B with a symmetric positive definite band matrix A. A has already undergone Cholesky factorization once before this line. 
	//Output overwrites b. 
	cblas_dcopy(size, b, 1, x, 1); //Outputs x. This function will be used to update streamfunction, hence parameter x represents array s in int main{} code. 
	
}
void PoissonSolver::PrintVector(int n, double* u) {
   for (int i = 0; i < n; ++i) {
       cout << u[i] << endl;
   }
   cout << endl;
}
void PoissonSolver::PrintMatrix(int n, int m, double* M) {
   for (int j = 0; j < m; ++j) {
       for (int i = 0; i < n; ++i) {
           cout << setw(6) << M[j+i*m] << " ";
       }
       cout << endl;
   }
}
