#ifndef SOLVE_BVP_H
#define SOLVE_BVP_H

// Declaration of the BVP Solver class
// This class can be used to compute the solution of linear and non-linear
// second order boundary value problems with Dirichlet boundary conditions
// R. Sheehan 5 - 4 - 2013

// This can be used to compute the solutions of linear BVP of the form
// y''(x) = p(x) y'(x) + q(x) y(x) + r(x), a < x < b, y(a) = \alpha, y(b) = \beta
// and non-linear BVP of the form
// y''(x) = f(x, y, y'), a < x < b, y(a) = \alpha, y(b) = \beta

class BVP_Solver{
public:
	// Constructors
	BVP_Solver(); 
	BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int )); // linear bvp constructor
	BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int )); // non-linear bvp constructor

	// methods
	void solve_linear_bvp();
	void solve_nonlinear_bvp_newton_raphson(double toler = 1.0e-6);
	void solve_nonlinear_bvp_shooting(double toler = 1.0e-6);

	void output_solution(std::string filename); 

	void clear(); // delete the arrays associated with the object

private:
	int N;
	int Ntotal; 
	int order; 
	int ndim; 

	double xl; // lower endpoint of domain
	double xu; // upper endpoint of domain
	
	double *alpha; // boundary conditions

	double **sol; // array to hold the solution being computed
	
	// the right-hand side functions are stored as elements of the vector func_vec
	coordfunc *func_vec; 
};

#endif