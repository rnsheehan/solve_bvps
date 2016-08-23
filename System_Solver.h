#ifndef SYSTEM_SOLVER_H
#define SYSTEM_SOLVER_H

// Solve the system of linear equations using LU Decomposition
// R. Sheehan 25 - 3 - 2013

namespace lin_slv{

	void ludcmp(double **a, int *indx, int size, double *d); 
	void lubksb(double **a, int *indx, double *b, int size);
	void luinv(double **a, double **y, double *col, int *indx, int size);
	void ludet(double **a, double *d, int size);

	double *solve_system(double **a, double *b, int size); 

	double **find_inverse(double **a, int size); 

	//double **sherman_morrison_inverse(double **a, double *u, double *v, int size); 

	double *solver_tri_diag_system(double *a, double *b, double *c, double *d, int size); 
}

#endif