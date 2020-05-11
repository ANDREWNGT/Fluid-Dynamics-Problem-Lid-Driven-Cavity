#ifndef POISSON_SOLVER
#define POISSON_SOLVER

//Header file for Poisson Solver. Solving poisson problem with Cholesky factorisation.

class PoissonSolver
{
public:
	PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy);
	~PoissonSolver(); //Destructor. 
	
	void FeedData(double xlen, double ylen, int nx, int ny,double deltat);
	void PrintVector(int n, double* u);
	void PrintMatrix(int n, int m, double* M);
	void SetBoundary(const double *top, const double *left, const double *bottom, const double *right);
	void Solve_Problem(double *x, const double *f);
	
private:
	/// Size of the linear system
	unsigned int Nx;
	unsigned int Ny;
	double dx;
	double dy;
	int size;
	int info;     

	double *A = nullptr;	// store A matrix 
	double A_MainDiag;      //Values given in eqn 10 and 11. 
	double A_NearDiag;
	double A_FarDiag;	
	int leading_a;	        // leading dimension of A	
	double *b = nullptr;	// store boundary vector for cblas routines
};

#endif































