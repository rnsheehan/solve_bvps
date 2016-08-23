#ifndef EXAMPLES_H
#define EXAMPLES_H

namespace examples{
	// Linear Boundary Value Problem 1
	double p(double *x, int n);
	double q(double *x, int n);
	double r(double *x, int n);

	double g11(double *x, int n); 
	double g12(double *x, int n);

	// Linear Boundary Value Problem 2
	double p1(double *x, int n);
	double q1(double *x, int n);
	double r1(double *x, int n);

	double g22(double *x, int n);
	double g32(double *x, int n);

	void tri_diag_test(); 

	void Q1(); 
	void Q1A(); 
	void Q2(); 
	void Q2A(); 
	void Q3(); 

}

#endif