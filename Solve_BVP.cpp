#ifndef ATTACH_H
#include "Attach.h"
#endif

// Source Definitions for the BVP Solver class

// Constructors
BVP_Solver::BVP_Solver()
{
	// Default constructor for the class
	order = ndim = N = Ntotal = 0; 

	xl = xu = 0.0; 

	alpha = nullptr; 
	sol = nullptr; 
	func_vec = nullptr; 
}

BVP_Solver::BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))
{
	// Constructor for a linear BVP
	// func1 = p(x), func2 = q(x), func3 = r(x)

	try{

		bool c1 = Npts > 3 ? true : false; 
		bool c2 = fabs(a-b) > 1.0e-9 ? true : false; 

		if(c1 && c2){

			order = 2; 
			ndim = order; // for linear bvp only computing {x, y(x)}
			N = Npts; // number of nodes internal nodes of the mesh
			Ntotal = N+2; 

			xl = std::min(a, b); xu = std::max(a, b); 

			// this will hold the boundary conditions
			alpha = new (double [order+1]); 

			alpha[1] = lbc; alpha[2] = rbc; 

			// This will hold {x, y}
			sol = new (double *[Ntotal+1]); 

			for(int i=1; i<=Ntotal; i++){
				sol[i] = new (double [ndim+1]); 
			}

			sol[1][1] = xl; sol[1][2] = alpha[1]; // {x_{0}, \alpha}
			sol[Ntotal][1] = xu; sol[Ntotal][2] = alpha[2]; // {x_{N+1}, \beta}

			// this will hold p, q and r
			func_vec = new (coordfunc [3+1]); 
			func_vec[1] = func1; // p(x)
			func_vec[2] = func2; // q(x)
			func_vec[3] = func3; // r(x)		
		}
		else{
			std::string reason = "Error: BVP_Solver::BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))\n"; 
			if(!c1) reason += "Number of sub-divisions is too small Npts = " + template_funcs::toString(Npts) + "\n"; 
			if(!c2) reason += "Endpoints of solution domain are not correctly defined a = " + template_funcs::toString(a, 4) + ", b = " + template_funcs::toString(a, 4) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: BVP_Solver::BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

BVP_Solver::BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ))
{
	// non-linear bvp constructor

	try{

		bool c1 = Npts > 3 ? true : false; 
		bool c2 = fabs(a-b) > 1.0e-9 ? true : false; 

		if(c1 && c2){

			order = 2; 
			ndim = order+1; // for non-linear bvp computing {x, y(x), y'(x)} when the shooting method is used
			N = Npts; 
			Ntotal = N+2; 

			xl = std::min(a, b); xu = std::max(a, b); 

			// this will hold the boundary conditions
			alpha = new (double [order+1]); 

			alpha[1] = lbc; alpha[2] = rbc; 

			// This will hold {x, y, y'}
			sol = new (double *[Ntotal+1]); 

			for(int i=1; i<=Ntotal; i++){
				sol[i] = new (double [ndim+1]); // extra column is only written to when using the shooting method
			}

			sol[1][1] = xl; sol[1][2] = alpha[1]; // {x_{0}, \alpha}
			sol[Ntotal][1] = xu; sol[Ntotal][2] = alpha[2]; // {x_{N+1}, \beta}

			// this will hold p, q and r
			func_vec = new (coordfunc [order+1]); 
			func_vec[1] = func1; // u_{2}(x) = y'(x)
			func_vec[2] = func2; // f(x, y, y')
		}
		else{
			std::string reason = "Error: BVP_Solver::BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ))\n"; 
			if(!c1) reason += "Number of sub-divisions is too small Npts = " + template_funcs::toString(Npts) + "\n"; 
			if(!c2) reason += "Endpoints of solution domain are not correctly defined a = " + template_funcs::toString(a, 4) + ", b = " + template_funcs::toString(a, 4) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: BVP_Solver::BVP_Solver(int Npts, double a, double b, double lbc, double rbc, double (*func1)(double *, int ), double (*func2)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

// Methods
void BVP_Solver::solve_linear_bvp()
{
	// Solve the linear bvp by the finite difference method
	// This method returns {x, y(x) } as columns of sol
	
	try{

		bool c1 = order > 0 ? true : false; 
		bool c2 = ndim >= order ? true : false; 
		bool c3 = N > 3 ? true : false; 

		if(c1 && c2 && c3){
			
			double *a = new(double [N+1]); 
			double *b = new(double [N+1]); 
			double *c = new(double [N+1]); 
			double *d = new(double [N+1]); 
			double *pos = new(double [N+1]); 

			double *x = new(double [2]);

			double h = ((xu-xl)/(N+1)); // node spacing
			double h_2 = 0.5*h; // node spacing / 2
			double h_sqr = template_funcs::DSQR(h); // node spacing^{2}

			double hpval;

			// define the array to hold the internal node positions	
			for(int i=1; i<=N; i++){
				pos[i] = xl + static_cast<double>(i*h); 
			}

			// Define the vectors a, b, c, d that form the system of linear equations
			// whose solution approximates the solution of the linear bvp
			for(int i=1; i<=N; i++){
		
				x[1] = pos[i]; // point to the current position

				hpval = h_2*func_vec[1](x,N); // (h/2)*p(x) is used twice for each i

				// Define a
				if(i>1){
					a[i] = -(1.0+hpval); 
				}

				// Define b
				b[i] = 2.0+h_sqr*func_vec[2](x,N); 

				// Define c
				if(i<N){
					c[i] = -(1.0-hpval); 
				}

				// Define d
				if(i==1){
					d[i] = -(h_sqr*func_vec[3](x,N))+(1.0+hpval)*alpha[1]; 
				}
				else if(i==N){
					d[i] = -(h_sqr*func_vec[3](x,N))+(1.0-hpval)*alpha[2]; 
				}
				else{
					d[i] = -(h_sqr*func_vec[3](x,N)); 
				}
			}

			double *y = lin_slv::solver_tri_diag_system(a, b, c, d, N); // Compute the solution at each of the nodes

			// Store the solution for output
			for(int i=1; i<=N; i++){
				sol[i+1][1] = pos[i]; 
				sol[i+1][2] = y[i]; 
			}

			// Print the solution
			/*for(int i=1; i<=(N+2);i++){
				std::cout<<sol[i][1]<<" , "<<sol[i][2]<<"\n";
			}
			std::cout<<"\n";*/

			delete[] a;
			delete[] b; 
			delete[] c; 
			delete[] d; 
			delete[] x; 
			delete[] y; 
			delete[] pos; 
		}
		else{
			std::string reason = "Error: void BVP_Solver::solve_linear_bvp()\n"; 
			if(!c1) reason += "Parameter order = " + template_funcs::toString(order) + " has not been assigned correctly\n"; 
			if(!c2) reason += "Parameter ndim = " + template_funcs::toString(ndim) + " has not been assigned correctly\n"; 
			if(!c3) reason += "Parameter N = " + template_funcs::toString(N) + " has not been assigned correctly\n"; 
			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void BVP_Solver::solve_linear_bvp()\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void BVP_Solver::solve_nonlinear_bvp_newton_raphson(double toler)
{
	// Solve a non-linear BVP using the newton-raphson method
	// This method returns {x, y(x), y'(x)} as columns of sol

	try{

		bool c1 = order > 0 ? true : false; 
		bool c2 = ndim > order ? true : false; 
		bool c3 = N > 3 ? true : false; 
		bool c4 = fabs(toler) > 1.0e-12 ? true : false; 

		if(c1 && c2 && c3 && c4){

			int j, iter, max_iter = 100;

			bool cgt = false; 

			std::string file; 

			// arrays required to compute the solution
			double *a = new(double [N+1]); 
			double *b = new(double [N+1]); 
			double *c = new(double [N+1]); 
			double *d = new(double [N+1]); 
			double *pos = new(double [N+1]); 

			double *y = new(double [N+1]); 
			double *v = new(double [N+1]); 
			double *x = new(double [ndim+1]);

			// parameters used repeatedly through the calculation
			double lr = xu-xl;
			double h = (lr/(N+1)); // node spacing
			double h_2 = 0.5*h; // node spacing / 2
			double h_sqr = template_funcs::DSQR(h); // node spacing^{2}
			double two_h = 2.0*h; 
			double two_h_sqr = 2.0*h_sqr; 
			double h_prod = (h_2)*(h_sqr); 

			double error, t2, m; 
			double dy, dy2, dy3, dy4, dy5, df1, df2, v1, v2; // variables used inside the main loop
	
			// Can't have m = 0 for initial approximation
			if(alpha[2] == alpha[1]){
				m = 2.0/lr; 
			}
			else{
				m = ((alpha[2]-alpha[1])/lr); //slope of initial approximation
			}

			// define the array to hold the internal node positions	
			// also define the initial approximation to the solution
			for(int i=1; i<=N; i++){
				pos[i] = xl + static_cast<double>(i*h); 
				y[i] = alpha[1] + m*( pos[i] - xl ); 
			}

			// Compute the solution iteratively
			iter = 1;
			while(iter < max_iter){
		
				// Define the vectors a, b, c, d that form the linear system of equations whose solution improves the estimate of y
				for(j=1; j<=N; j++){

					if( j == 1 ){
						// Deal with approximation to solution outside mesh
						// Include effect of lhbc

						dy = y[j+1] - alpha[1]; // y_{j+1} - y_{j-1}
						dy2 = dy/two_h; // y_{j+1} - y_{j-1} / 2 h

						// define df
						dy4 = y[j+2] + 2.0*( alpha[1] - y[j+1] ) - ( alpha[1] - h*( ( y[j] - alpha[1] ) / ( pos[j] - xl ) ) ); 

						// define the r.h.s. vector
						x[1] = pos[j]; x[2] = y[j]; x[3] = dy2; // updated at each step

						d[j] = -1.0*( 2.0*y[j] - alpha[1] - y[j+1] + h_sqr*(func_vec[2](x,3)) );
					}
					else if( j == 2 ){
						// Include effect of lhbc

						dy = y[j+1] - y[j-1]; // y_{j+1} - y_{j-1}
						dy2 = dy/two_h; // y_{j+1} - y_{j-1} / 2 h

						// define df
						dy4 = y[j+2] + 2.0*( y[j-1] - y[j+1] ) - alpha[1]; 

						// define the r.h.s. vector
						x[1] = pos[j]; x[2] = y[j]; x[3] = dy2; // updated at each step

						d[j] = -1.0*( 2.0*y[j] - y[j-1] - y[j+1] + h_sqr*(func_vec[2](x,3)) );
					}
					else if( j == N-1 ){
						// Include effect of rhbc

						dy = y[j+1] - y[j-1]; // y_{j+1} - y_{j-1}
						dy2 = dy/two_h; // y_{j+1} - y_{j-1} / 2 h

						// define df
						dy4 = alpha[2] + 2.0*( y[j-1] - y[j+1] ) - y[j-2]; 

						// define the r.h.s. vector
						x[1] = pos[j]; x[2] = y[j]; x[3] = dy2; // updated at each step

						d[j] = -1.0*( 2.0*y[j] - y[j-1] - y[j+1] + h_sqr*(func_vec[2](x,3)) );
					}
					else if( j == N ){
						// Deal with approximation to solution outside mesh
						// Include effect of rhbc

						dy = alpha[2] - y[j-1]; // y_{j+1} - y_{j-1}
						dy2 = dy/two_h; // y_{j+1} - y_{j-1} / 2 h

						// define df
						dy4 = ( alpha[2] + h*( (alpha[2]-y[j]) / (xu-pos[j]) ) ) + 2.0*( y[j-1] - alpha[2] ) - y[j-2]; 

						// define the r.h.s. vector
						x[1] = pos[j]; x[2] = y[j]; x[3] = dy2; // updated at each step

						d[j] = -1.0*( 2.0*y[j] - y[j-1] - alpha[2] + h_sqr*func_vec[2](x,3) );
					}
					else{
						// All points in between

						dy = y[j+1] - y[j-1]; // y_{j+1} - y_{j-1}
						dy2 = dy/two_h; // y_{j+1} - y_{j-1} / 2 h

						// define df
						dy4 = y[j+2] + 2.0*( y[j-1] - y[j+1] ) - y[j-2]; 

						// define the r.h.s. vector
						x[1] = pos[j]; x[2] = y[j]; x[3] = dy2; // updated at each step

						d[j] = -1.0*( 2.0*y[j] - y[j-1] - y[j+1] + h_sqr*func_vec[2](x,3) );
					}

					// define df / dy and df / dy'
					dy3 = 0.5*dy; 
					dy5 = dy4 / two_h_sqr;

					x[3] = dy2 + dy5; v1 = func_vec[2](x,3); 
					x[3] = dy2 - dy5; v2 = func_vec[2](x,3); 
					df1 = (v1 - v2); 

					if( fabs(dy4) < 1.0e-9 ){
						df2 = 1.0; 
					}
					else{
						df2 = (h_prod*df1)/dy4; 
					}

					// Define the sub and super diagonals
					a[j] = -1.0 - df2; c[j] = -1.0 + df2; 
			
					// Define the main diagonal
					x[3] = dy2;
					x[2] = y[j] + dy3; v1 = func_vec[2](x,3); 
					x[2] = y[j] - dy3; v2 = func_vec[2](x,3);			
					b[j] = 2.0 + (h_sqr/dy)*(v1 - v2); 	
				}

				// Solve the system of equations to find the update vector
				v = lin_slv::solver_tri_diag_system(a, b, c, d, N);

				// Add the update vector to the old solution 
				for(j=1; j<=N; j++){
					y[j] += v[j]; 
				}

				// Check for convergence relative to the infinity norm
				error=0.0; 
				for(j=1; j<=N; j++){
					t2 = fabs(v[j]); 
					if( t2 > error ){
						error = t2; 
					}
				}

				if(error < toler){
					cgt = true; 
					break; 
				}

				iter++; 
			}

			if(cgt){
				std::cout<<"The newton-raphson method has converged within "<<iter<<" iterations\n"; 
			}
			else{
				std::cout<<"The newton-raphson method has failed converged within "<<max_iter<<" iterations\n"; 
			}

			// Store the solution for output
			for(int i=1; i<=N; i++){
				sol[i+1][1] = pos[i]; 
				sol[i+1][2] = y[i]; 
				sol[i+1][3] = 0.0; 
			}

			file = "Computed_Solution_Iteration_"+template_funcs::toString(iter)+".txt"; 

			output_solution(file);

			delete[] a;
			delete[] b; 
			delete[] c; 
			delete[] d; 
			delete[] x; 
			delete[] y; 
			delete[] v; 
			delete[] pos;
		}
		else{
			std::string reason = "Error: void BVP_Solver::solve_nonlinear_bvp_newton_raphson(double toler)\n"; 
			if(!c1) reason += "Parameter order = " + template_funcs::toString(order) + " has not been assigned correctly\n"; 
			if(!c2) reason += "Parameter ndim = " + template_funcs::toString(ndim) + " has not been assigned correctly\n"; 
			if(!c3) reason += "Parameter N = " + template_funcs::toString(N) + " has not been assigned correctly\n"; 
			if(!c4) reason += "Parameter toler = " + template_funcs::toString(toler,6) + " is too small\n"; 
			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void BVP_Solver::solve_nonlinear_bvp_newton_raphson(double toler)\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void BVP_Solver::solve_nonlinear_bvp_shooting(double toler)
{
	// Solve a non-linear BVP using the shooting method
	// This method returns {x, y(x)} as columns of sol

	try{

		bool c1 = order > 0 ? true : false; 
		bool c2 = ndim > order ? true : false; 
		bool c3 = N > 3 ? true : false; 
		bool c4 = fabs(toler) > 1.0e-12 ? true : false; 

		if(c1 && c2 && c3 && c4){
		
			// This method works by iteratively computing the solution to a second order IVP
			// This will require an IVP_Solver object

			int iter = 1, max_iter = 500; // iteration control

			bool cgt = false; 

			std::string file; 

			double **sol_ptr; // pointer to the solution from the other object
			double lr = xu-xl; // domain length
			double t0, t1, t2, dt; // shooting parameters
			double yold, yolder, dy; // solution values to be stored	

			IVP_Solver eqn; // IVP Solver object

			while(iter < max_iter){

				// Compute the shooting parameter
				if(iter == 1){
					// Define shooting parameter for first iteration
					t0 = (alpha[2] - alpha[1]) / lr; 
					t2 = t0; 
				}
				else if(iter == 2){
					// Define shooting parameter for second iteration
					t1 = t0 + (alpha[2]-sol[Ntotal][2])/lr; 
					t2 = t1; 
				}
				else{
					// Continue computing shooting parameters, with secant method, until convergence
					dy = yold - yolder; 
			
					if(fabs(dy) > 1.0e-12){
						// Compute the required value of dt by the secant method
						dt = ((yold - alpha[2])*(t1-t0))/dy; 
					}
					else{
						// solution is not likely to improve
						dt = 1.0; 
						break; 
					}

					t2 = t1 - dt; // update the value of the shooting parameter
					t0 = t1; 
					t1 = t2; 

				}

				// Compute the solution of the IVP with the given shooting parameter
				eqn.set_ivp2( Ntotal, xl, xu, alpha[1], t2, func_vec[1], func_vec[2] ); 
				eqn.RK4_Solve(); 
		
				//file = "Computed_Solution_Iteration_"+toString(iter)+".txt"; 
				//eqn.output_solution(file);
		
				// Store the computed solution
				sol_ptr = eqn.get_sol();

				for(int i=1; i<=Ntotal; i++){
					for(int d=1; d<=ndim; d++){
						sol[i][d] = sol_ptr[i][d]; 
					}
				}
		
				eqn.clear(); // clear the memory used by the IVP solver

				// Store the required value from the computed solution
				// this value is y(b, t_{k})
				if(iter == 1){
					yolder = sol[Ntotal][2]; 
				}
				else if(iter == 2){
					yold = sol[Ntotal][2]; 
				}
				else{
					yolder = yold; 
					yold = sol[Ntotal][2]; 
				}

				// Test for convergence
				if(iter > 2 && fabs(sol[Ntotal][2]-alpha[2]) < toler){
					cgt = true; 
					break; 
				}
		
				iter++; 
			}

			if(cgt){
				std::cout<<"The shooting method has converged within "<<iter<<" iterations\n"; 
				std::cout<<"t = "<<t2<<", y(b) = "<<sol[Ntotal][2]<<", beta = "<<alpha[2]<<"\n";
			}
			else{
				std::cout<<"The shooting method has failed converged within "<<max_iter<<" iterations\n"; 
				std::cout<<"t = "<<t2<<", y(b) = "<<sol[Ntotal][2]<<", beta = "<<alpha[2]<<"\n";
			}

			file = "Computed_Solution_Iteration_"+template_funcs::toString(iter)+".txt"; 

			output_solution(file);
		}
		else{
			std::string reason = "Error: void BVP_Solver::solve_nonlinear_bvp_shooting(double toler)\n"; 
			if(!c1) reason += "Parameter order = " + template_funcs::toString(order) + " has not been assigned correctly\n"; 
			if(!c2) reason += "Parameter ndim = " + template_funcs::toString(ndim) + " has not been assigned correctly\n"; 
			if(!c3) reason += "Parameter N = " + template_funcs::toString(N) + " has not been assigned correctly\n"; 
			if(!c4) reason += "Parameter toler = " + template_funcs::toString(toler,6) + " is too small\n"; 
			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void BVP_Solver::solve_nonlinear_bvp_shooting(double toler)\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void BVP_Solver::output_solution(std::string filename)
{
	// output the computed solution to a particular file

	// ensure that filename has a non null-string value
	if(filename == ""){
		filename = "BVP_Solver_Output.txt"; 
	}

	std::ofstream write; 
	write.open(filename.c_str(),std::ios_base::out|std::ios_base::trunc);

	if(write.is_open()){
		
		for(int i=1; i<=(N+2); i++){
			for(int d=1; d<=ndim; d++)
				if(d == ndim){
					write<<std::setprecision(12)<<sol[i][d];
				}
				else{
					write<<std::setprecision(12)<<sol[i][d]<<",";
				}
			write<<"\n";
		}

		write.close();

	}
	else{
		std::cerr<<"Could not open "<<filename<<"\n";
	}
}

void BVP_Solver::clear()
{
	// delete the arrays associated with the object

	if(alpha != nullptr) delete[] alpha; 
	if(sol != nullptr) delete[] sol; 
	if(func_vec != nullptr) delete[] func_vec; 
}