#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions of the test functions

// Linear Boundary Value Problem 1
double examples::p(double *x, int n)
{
	return (-2.0/x[1]); 
}

double examples::q(double *x, int n)
{
	return 0.0; 
}

double examples::r(double *x, int n)
{
	return 0.0; 
}


double examples::g11(double *x, int n)
{
	// this is equivalent to u_{2}(x)

	return (x[3]); 
}

double examples::g12(double *x, int n)
{
	// this is equivalent to f(x, y, y')

	return ((-2.0*x[3])/x[1]);
}

// Linear Boundary Value Problem 2
double examples::p1(double *x, int n)
{
	return 0.0; 
}

double examples::q1(double *x, int n)
{
	double S = 1000.0;
	double E = 3.0e7;
	double I = 625.0; 

	return ( S / ( E * I ) ); 
}

double examples::r1(double *x, int n)
{
	double Q = 25.0; // changed from 100
	double E = 3.0e7;
	double I = 625.0; 
	double l = 120.0; 
	double k = ( Q / ( 2.0 * E * I ) );

	return ( k * x[1] * ( x[1] - l ) ); // k x ( x - l )
}

double examples::g22(double *x, int n)
{
	
	return ( ( q1(x,n) * x[2] ) + ( r1(x,n) ) ); // k_{1} w + k_{2} x ( x - l )
}

double examples::g32(double *x, int n)
{
	double t3 = ( 1.0 + template_funcs::DSQR(x[3]) ); // ( 1 + (w')^{2} )
	double t4 = sqrt( t3 * template_funcs::DSQR(t3) ); // t3^{3/2}

	return ( g22(x,n) * t4 ); 
}

void examples::tri_diag_test()
{
	int n=7; 

	double *a = new(double [n+1]); 
	double *b = new(double [n+1]); 
	double *c = new(double [n+1]); 
	double *d = new(double [n+1]); 

	a[1] = 0.0; a[2] = -0.9; a[3] = -0.909091; a[4] = -0.916667;
	a[5] = -0.923077; a[6] = -0.928571; a[7] = -0.933333; 

	b[1] = 2.0; b[2] = 2.0; b[3] = 2.0; b[4] = 2.0;
	b[5] = 2.0; b[6] = 2.0; b[7] = 2.0;

	c[1] = -1.11111; c[2] = -1.1; c[3] = -1.09091; c[4] = -1.08333;
	c[5] = -1.07692; c[6] = -1.07143; c[7] = 0.0;

	d[1] = 97.7778; d[2] = 0.0; d[3] = 0.0; d[4] = 0.0;
	d[5] = 0.0; d[6] = 0.0; d[7] = 0.0;

	double *sol=lin_slv::solver_tri_diag_system(a,b,c,d,n);

	for(int i=1; i<=n; i++){
		std::cout<<sol[i]<<"\n";
	}
	std::cout<<"\n";

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
}

void examples::Q1()
{
	// Solve the following linear BVP
	// \begin{equation*}
	//	\begin{cases}
	//	\displaystyle \frac{d^{2}u}{d r^{2}}+\frac{2}{r}\frac{d u}{d r}\,=\,0\quad R_{1}\,\leq\,r\,\leq\,R_{2}\\
	//	u(R_{1})\,=\,V_{1}\quad u(R_{2})\,=\,0
	//	\end{cases}
	// \end{equation*}

	int npts; 
	double a, b, alpha, beta; 

	npts = 77; 
	a = 2.0; b=4.0; alpha = 110.0; beta = 0.0; 

	BVP_Solver test(npts,a,b,alpha,beta,p,q,r);

	test.solve_linear_bvp(); 

	test.output_solution("Solution_Q1.txt"); 
}

void examples::Q1A()
{
	// Solve the following linear BVP using the methods for non-linear equations
	// \begin{equation*}
	//	\begin{cases}
	//	\displaystyle \frac{d^{2}u}{d r^{2}}+\frac{2}{r}\frac{d u}{d r}\,=\,0\quad R_{1}\,\leq\,r\,\leq\,R_{2}\\
	//	u(R_{1})\,=\,V_{1}\quad u(R_{2})\,=\,0
	//	\end{cases}
	// \end{equation*}

	int npts; 
	double a, b, alpha, beta; 

	npts = 77; 
	a = 2.0; b=4.0; alpha = 110.0; beta = 0.0; 

	BVP_Solver test(npts, a, b, alpha , beta, g11, g12);

	//test.solve_nonlinear_bvp_shooting(); 
	test.solve_nonlinear_bvp_newton_raphson(); 

	test.output_solution("Solution_Q1A.txt"); 
}

void examples::Q2()
{
	// Solve the following linear BVP from structural engineering
	/*
	
	\begin{equation*}
        \begin{cases}
        \displaystyle w''(x)\,=\,\frac{S}{E I}w+\frac{q\,x}{2\,E\,I}(x-l),\,0\,\leq\,x\,\leq\,l\\
        w(0)\,=\,0,\,w\left(l\right)\,=\,0
        \end{cases}
    \end{equation*}
	
	*/

	int npts; 
	double a, b, alpha, beta; 

	npts = 107; 
	a = 0.0; b=120.0; alpha = 0.0; beta = 0.0; 

	BVP_Solver test(npts,a,b,alpha,beta,p1,q1,r1);
	
	test.solve_linear_bvp(); 	

	test.output_solution("Solution_Q2.txt");
}

void examples::Q2A()
{
	// Solve the following linear BVP from structural engineering using the non-linear methods
	/*
	
	\begin{equation*}
        \begin{cases}
        \displaystyle w''(x)\,=\,\frac{S}{E I}w+\frac{q\,x}{2\,E\,I}(x-l),\,0\,\leq\,x\,\leq\,l\\
        w(0)\,=\,0,\,w\left(l\right)\,=\,0
        \end{cases}
    \end{equation*}
	
	*/

	int npts; 
	double a, b, alpha, beta; 

	npts = 107; 
	a = 0.0; b=120.0; alpha = 0.0; beta = 0.0; 

	BVP_Solver test(npts,a,b,alpha,beta,g11,g22);

	test.solve_nonlinear_bvp_shooting(); 
	//test.solve_nonlinear_bvp_newton_raphson(); 

	test.output_solution("Solution_Q2A.txt");
}

void examples::Q3()
{
	// Solve the following non-linear BVP from structural engineering using the non-linear methods
	/*
	
	\begin{equation*}
        \begin{cases}
        \displaystyle (1+(w'(x))^{2})^{-3/2}w''(x)\,=\,\frac{S}{E I}w+\frac{q\,x}{2\,E\,I}(x-l),\,0\,\leq\,x\,\leq\,l\\
        w(0)\,=\,0,\,w\left(l\right)\,=\,0
        \end{cases}
    \end{equation*}
	
	*/

	int npts; 
	double a, b, alpha, beta; 

	npts = 107; 
	a = 0.0; b = 120.0; alpha = 0.0; beta = 0.0; 

	BVP_Solver test(npts,a,b,alpha,beta,g11,g32); 

	test.solve_nonlinear_bvp_shooting();
	//test.solve_nonlinear_bvp_newton_raphson(1.0e-6);

	test.output_solution("Solution_Q3.txt");	
}

