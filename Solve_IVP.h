#ifndef SOLVE_IVP_H
#define SOLVE_IVP_H

// Declaration of the IVP Solver class
// This class can be used to solve IVP of the first, second and third orders
// The methods implemented here include Runge-Kutta order four and a predictor-corrector scheme based on 
// a combination of Adams-Bashforth (Predictor) and Adams-Moulton (Corrector). 
// The ABM-PC method uses RK4 to compute the first four steps of its solution. 
// R. Sheehan 26 - 3 - 2013

class IVP_Solver{

public:
	// Constructors
	IVP_Solver(); // Default
	IVP_Solver(int Nsteps, double a, double b, double ic, double (*func1)(double *, int )); // First order IVP
	IVP_Solver(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int )); // Second order IVP and Pairs of IVP
	IVP_Solver(int Nsteps, double a, double b, double ic1, double ic2, double ic3, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int )); // Third order IVP and Triples of IVP

	// Setters
	void set_ivp2(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int )); // set the parameters required to solve a second order IVP

	// Getters
	inline double **get_sol(){return sol;}

	// Methods
	void RK4_Solve(); // Solver based on Runge-Kutta Method Order Four
	void ABM_PC_Solve(); // Adams-Bashforth-Moulton Predictor Corrector

	void output_solution(std::string filename); 

	void clear(); // delete the arrays associated with the object

	// End-user does not need access to these methods
private:
	double *RK4_Step(double *x); // Step Calculator based on Runge-Kutta Method Order Four
	double *AB_Step(double **x); // Step Calculator based on four step Adams-Bashforth method
	double *AM_Step(double **x); // Step Calculator based on three step Adams-Moulton method

private:
	int N; // number of steps at which the solution is to be computed
	int order; // order of the equations being solved 1^{st}, 2^{nd}, 3^{rd} etc...
	int ndim; // num. dimensions over which solution is sought ndim = order + 1 always

	double h; // node spacing
	double xl; // lower endpoint of domain
	double xu; // upper endpoint of domain
	
	double *alpha; // initial conditions

	double **sol; // array to hold the solution being computed
	
	// the right-hand side functions are stored as elements of the vector func_vec
	coordfunc *func_vec; 
};

#endif