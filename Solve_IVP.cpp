#ifndef ATTACH_H
#include "Attach.h"
#endif

// Constructors

IVP_Solver::IVP_Solver()
{
	// Default constructor

	order = ndim = N = 0; 

	h = xl = xu = 0.0;

	alpha = nullptr; 
	sol = nullptr; 
	func_vec = nullptr; 
}

IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic, double (*func1)(double *, int ))
{
	// Constructor for first order IVP

	try{

		bool c1 = Nsteps > 3 ? true : false; 
		bool c2 = fabs(a-b) > 1.0e-9 ? true : false; 

		if(c1 && c2){

			order = 1; 
			ndim = order+1; 
			N = Nsteps; 

			xl = std::min(a, b); xu = std::max(a, b); h = ((xu-xl)/(N-1)); 

			alpha = new (double [order+1]); 

			alpha[1] = ic; 

			sol = new (double *[N+1]); 

			for(int i=1; i<=N; i++){
				sol[i] = new (double [ndim+1]); 
			}

			func_vec = new (coordfunc [order+1]); 
			func_vec[1] = func1; 
		
		}
		else{
			std::string reason = "Error: IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic, double (*func1)(double *, int ))\n"; 
			if(!c1) reason += "Number of steps is too small Nsteps = " + template_funcs::toString(Nsteps) + "\n"; 
			if(!c2) reason += "Endpoints of solution domain are not correctly defined a = " + template_funcs::toString(a, 4) + ", b = " + template_funcs::toString(a, 4) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic, double (*func1)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int ))
{
	// Constructor for second order IVP

	try{

		bool c1 = Nsteps > 3 ? true : false; 
		bool c2 = fabs(a-b) > 1.0e-9 ? true : false; 

		if(c1 && c2){

			order = 2; 
			ndim = order+1; 
			N = Nsteps; 

			xl = std::min(a, b); xu = std::max(a, b); h = ((xu-xl)/(N-1)); 

			alpha = new (double [order+1]); 

			alpha[1] = ic1; alpha[2] = ic2; 

			sol = new (double *[N+1]); 

			for(int i=1; i<=N; i++){
				sol[i] = new (double [ndim+1]); 
			}

			func_vec = new (coordfunc [order+1]); 
			func_vec[1] = func1; 
			func_vec[2] = func2;		
		}
		else{
			std::string reason = "Error: IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic, double (*func1)(double *, int ))\n"; 
			if(!c1) reason += "Number of steps is too small Nsteps = " + template_funcs::toString(Nsteps) + "\n"; 
			if(!c2) reason += "Endpoints of solution domain are not correctly defined a = " + template_funcs::toString(a, 4) + ", b = " + template_funcs::toString(a, 4) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic1, double ic2, double ic3, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))
{
	// Constructor for third order IVP

	try{

		bool c1 = Nsteps > 3 ? true : false; 
		bool c2 = fabs(a-b) > 1.0e-9 ? true : false; 

		if(c1 && c2){

			order = 3; 
			ndim = order+1; 
			N = Nsteps; 

			xl = std::min(a, b); xu = std::max(a, b); h = ((xu-xl)/(N-1)); 

			alpha = new (double [order+1]); 

			alpha[1] = ic1; alpha[2] = ic2; alpha[3] = ic3; 

			sol = new (double *[N+1]); 

			for(int i=1; i<=N; i++){
				sol[i] = new (double [ndim+1]); 
			}

			func_vec = new (coordfunc [order+1]); 
			func_vec[1] = func1; 
			func_vec[2] = func2; 
			func_vec[3] = func3; 
		
		}
		else{
			std::string reason = "Error: IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic, double (*func1)(double *, int ))\n"; 
			if(!c1) reason += "Number of steps is too small Nsteps = " + template_funcs::toString(Nsteps) + "\n"; 
			if(!c2) reason += "Endpoints of solution domain are not correctly defined a = " + template_funcs::toString(a, 4) + ", b = " + template_funcs::toString(a, 4) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: IVP_Solver::IVP_Solver(int Nsteps, double a, double b, double ic1, double ic2, double ic3, double (*func1)(double *, int ), double (*func2)(double *, int ), double (*func3)(double *, int ))\n"; 
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

void IVP_Solver::set_ivp2(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int ))
{
	// Set the parameters required to compute the solution of a second order IVP

	try{

		bool c1 = Nsteps > 3 ? true : false; 
		bool c2 = fabs(a-b) > 1.0e-9 ? true : false; 

		if(c1 && c2){

			order = 2; 
			ndim = order+1; 
			N = Nsteps; 

			xl = std::min(a, b); xu =std::max(a, b); h = ((xu-xl)/(N-1)); 

			alpha = new (double [order+1]); 

			alpha[1] = ic1; alpha[2] = ic2; 

			sol = new (double *[N+1]); 

			for(int i=1; i<=N; i++){
				sol[i] = new (double [ndim+1]); 
			}

			func_vec = new (coordfunc [order+1]); 
			func_vec[1] = func1; 
			func_vec[2] = func2;	
		}
		else{
			std::string reason = "Error: void IVP_Solver::set_ivp2(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int ))\n"; 
			if(!c1) reason += "Number of steps is too small Nsteps = " + template_funcs::toString(Nsteps) + "\n"; 
			if(!c2) reason += "Endpoints of solution domain are not correctly defined a = " + template_funcs::toString(a, 4) + ", b = " + template_funcs::toString(a, 4) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void IVP_Solver::set_ivp2(int Nsteps, double a, double b, double ic1, double ic2, double (*func1)(double *, int ), double (*func2)(double *, int ))\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	} 
}

double* IVP_Solver::RK4_Step(double *x)
{
	// Single Step Calculator for Runge-Kutta Method Order Four
	// x contains the solution from the previous step

	try{		

		bool c1 = x != nullptr ? true : false; 
		
		if(c1){		

			int i,j,d,m; 

			static const int nconsts = 4; // I think this was meant for something else, probably no need for it to be static const

			double *k = new (double [order+1]); // store the stepper constants for each step
			double *xhere = new (double [ndim+1]); // store the inputs for the function evaluations
			double *deltay = new (double [order+1]); // store the amounts by which the solution is updated
			double *xnew = new (double [ndim+1]); // store the value of the solution at the next step

			for(i=1; i<=order; i++){
				deltay[i] = 0.0; 
			}

			// Compute the four stepper constants	
			for(j=1; j<=nconsts; j++){

				// Define the inputs to the function evaluations for each of the stepper constants
				if(j==1){
					for(int m=1; m<=ndim; m++){
						xhere[m] = x[m]; 
					}
				}
				else if(j==2 || j==3){
					for(m=1; m<=ndim; m++){
						if(m==1){
							xhere[m] = x[m]+0.5*h; 
						}
						else{
							xhere[m] = x[m]+0.5*k[m-1]; 	
						}
					}
				}
				else if(j==4){
					for(m=1; m<=ndim; m++){
						if(m==1){
							xhere[m] = x[m] + h; 
						}
						else{
							xhere[m] = x[m]+k[m-1]; 	
						}
					}
				}

				// Compute the values of the stepper constants and the values by which the solution is to be updated
				for(d=1; d<=order; d++){

					k[d] = h*func_vec[d](xhere,ndim); 

					//std::cout<<k[d]<<" "; 
			
					if(j==1){
						deltay[d] = k[d]; 
					}
					else if(j==2 || j==3){
						deltay[d] += 2.0*k[d]; 
					}
					else if(j==4){
						deltay[d] += k[d]; 
					}
				}
				//std::cout<<endl;
			}
			//std::cout<<endl;
	
			// Update the value of the solution from the previous step
			for(d=1; d<=order; d++){
				//std::cout<<x[d+1]+(deltay[d]/6.0)<<" ";
				deltay[d] /= 6.0; 
				deltay[d] += x[d+1]; 
				/*std::cout<<deltay[d]<<endl;*/
			}

			// Store the solution at the next step
			for(d=1; d<=ndim; d++){
				if(d==1){
					xnew[d] = x[d] + h; 
				}
				else{
					xnew[d] = deltay[d-1]; 
				}
				/*std::cout<<xnew[d]<<" , ";*/
			}
			/*std::cout<<endl;*/

			delete[] k; 
			delete[] xhere;
			delete[] deltay; 

			return xnew; 
		}
		else{
			std::string reason = "Error: double* IVP_Solver::RK4_Step(double *x)\n"; 
			reason += "Input vector x was not correctly allocated\n";  
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: double* IVP_Solver::RK4_Step(double *x)\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double* IVP_Solver::AB_Step(double **x)
{
	// Single Step Calculator for four step Adams-Bashforth Method Order Four
	// x contains the solutions from the previous four steps
	// x[1] = {t_{i}, x_{i}, y_{i}, z_{i}, ....}
	// x[2] = {t_{i-1}, x_{i-1}, y_{i-1}, z_{i-1}, ....}
	// x[3] = {t_{i-2}, x_{i-2}, y_{i-2}, z_{i-2}, ....}
	// x[4] = {t_{i-3}, x_{i-3}, y_{i-3}, z_{i-3}, ....}

	// order equals the number of solutions being sought
	// ndim = order + 1 

	try{		

		bool c1 = x != nullptr ? true : false; 
		
		if(c1){		

			int j,d; 

			static const int nconsts = 4; // nconsts = nsteps for the PC stepper

			double k, k1 = (h / 24.0); 
			double mult[nconsts+1] = {0.0, 55.0, -59.0, 37.0, -9.0}; // store the multipliers for the function evaluations needed by the stepper

			double *deltay = new (double [order+1]); // store the amounts by which the solution is updated
			double *xnew = new (double [ndim+1]); // store the value of the solution at the next step

			for(j=1; j<=order; j++){
				deltay[j] = 0.0; 
			}

			// Compute the four stepper constants	
			for(d=1; d<=order; d++){
		
				for(j=1; j<=nconsts; j++){

					k = mult[j]*func_vec[d](x[j],ndim); 

					//std::cout<<k<<" ";

					if(j==1){
						deltay[d] = k; 
					}
					else{
						deltay[d] += k; 
					}
				}

				// Update the value of the solution from the previous step
				deltay[d] *= k1; // multiply by h/24
				deltay[d] += x[1][d+1]; // Add the updated value to the value from y_{i}

				/*std::cout<<deltay[d]<<endl;*/
			}
			//std::cout<<endl;

			// Store the solution at the next step
			for(d=1; d<=ndim; d++){
				if(d==1){
					xnew[d] = x[1][d] + h; 
				}
				else{
					xnew[d] = deltay[d-1]; 
				}
				/*std::cout<<xnew[d]<<" , ";*/
			}
			/*std::cout<<endl;*/

			delete[] deltay; 

			return xnew; 
		}
		else{
			std::string reason = "Error: double* IVP_Solver::AB_Step(double **x)\n"; 
			reason += "Input vector x was not correctly allocated\n";  
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: double* IVP_Solver::AB_Step(double **x)\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double* IVP_Solver::AM_Step(double **x)
{
	// Single Step Calculator for three step Adams-Moulton Method Order Four
	// x contains the solutions from the previous four steps
	// x[1] = {t_{i+1}, x_{i+1}, y_{i+1}, z_{i+1}, ....}
	// x[2] = {t_{i}, x_{i}, y_{i-1}, z_{i}, ....}
	// x[3] = {t_{i-1}, x_{i-1}, y_{i-1}, z_{i-1}, ....}
	// x[4] = {t_{i-2}, x_{i-2}, y_{i-2}, z_{i-2}, ....}

	// order equals the number of solutions being sought
	// ndim = order + 1 

	try{		

		bool c1 = x != nullptr ? true : false; 
		
		if(c1){		

			int j,d; 

			static const int nconsts = 4; // nconsts = nsteps for the PC stepper

			double k, k1 = (h / 24.0); 
			double mult[nconsts+1] = {0.0, 9.0, 19.0, -5.0, 1.0}; // store the multipliers for the function evaluations needed by the stepper

			double *deltay = new (double [order+1]); // store the amounts by which the solution is updated
			double *xnew = new (double [ndim+1]); // store the value of the solution at the next step

			for(j=1; j<=order; j++){
				deltay[j] = 0.0; 
			}

			// Compute the four stepper constants	
			for(d=1; d<=order; d++){
		
				for(j=1; j<=nconsts; j++){

					k = mult[j]*func_vec[d](x[j],ndim); 

					/*std::cout<<k<<" , "<<fval<<endl;*/

					if(j==1){
						deltay[d] = k; 
					}
					else{
						deltay[d] += k; 
					}
				}

				// Update the value of the solution from the previous step
				deltay[d] *= k1; // multiply by h/24
				deltay[d] += x[2][d+1]; // Add the update y value to the value from y_{i}, not y_{i+1}

				/*std::cout<<deltay[d]<<endl;*/
			}
			/*std::cout<<endl;*/

			// Store the solution at the next step
			for(d=1; d<=ndim; d++){
				if(d==1){
					xnew[d] = x[1][d]; 
				}
				else{
					xnew[d] = deltay[d-1]; // d-1 because xnew is 1 larger than deltay
				}
				/*std::cout<<xnew[d]<<" , ";*/
			}
			/*std::cout<<endl<<endl;*/

			delete[] deltay; 

			return xnew; 
		}
		else{
			std::string reason = "Error: double* IVP_Solver::AM_Step(double **x)\n"; 
			reason += "Input vector x was not correctly allocated\n";  
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: double* IVP_Solver::AM_Step(double **x)\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

void IVP_Solver::RK4_Solve()
{
	// Compute the solution of an IVP using Runge-Kutta order four

	try{		
		bool c1 = order > 0 ? true : false; 
		bool c2 = ndim > order ? true : false; 
		bool c3 = N > 3 ? true : false; 

		if(c1 && c2 && c3){

			int i,j,d; 

			double *ynew;
			double *yold = new (double [ndim+1]);

			// Initialise the solution
			yold[1] = xl; 
			for(j=1; j<=order;j++){
				yold[j+1] = alpha[j]; 
			}

			// Store the solution
			for(d=1; d<=ndim; d++){
				sol[1][d] = yold[d]; 
			}

			// Step the solution forward
			for(i=2; i<=N; i++){
		
				ynew = RK4_Step(yold); // Compute the solution at the next using RK4

				for(d=1; d<=ndim; d++){
					sol[i][d] = ynew[d]; 
				}
		
				yold = ynew; 
			}

			delete[] yold; 		
		}
		else{
			std::string reason = "Error: void IVP_Solver::RK4_Solve()\n"; 
			if(!c1) reason += "Parameter order = " + template_funcs::toString(order) + " has not been assigned correctly\n"; 
			if(!c2) reason += "Parameter ndim = " + template_funcs::toString(ndim) + " has not been assigned correctly\n"; 
			if(!c3) reason += "Parameter N = " + template_funcs::toString(N) + " has not been assigned correctly\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void IVP_Solver::RK4_Solve()\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

void IVP_Solver::ABM_PC_Solve()
{
	// Compute the solution of an IVP using Adams-Bashforth-Moulton predictor corrector

	try{		
		bool c1 = order > 0 ? true : false; 
		bool c2 = ndim > order ? true : false; 
		bool c3 = N > 3 ? true : false; 

		if(c1 && c2 && c3){

			int i,j,d,nsteps = 4; 

			double *ynew;
			//double *yold = new (double [ndim+1]);

			double **ABMinput = new (double *[nsteps+1]); // use this to store the solution from previous four steps

			for(i=1; i<=nsteps; i++){
				ABMinput[i] = new (double [ndim+1]); 
			}

			// Initialise the solution
			ABMinput[1][1] = xl; 
			for(j=1; j<=order;j++){
				ABMinput[1][j+1] = alpha[j]; 
			}

			// Store the solution
			for(d=1; d<=ndim; d++){
				sol[1][d] = ABMinput[1][d]; 
			}

			// Step the solution forward
			for(i=2; i<=N; i++){

				if(i<5){
					ynew = RK4_Step(ABMinput[1]); // Compute the solution at the next step using RK4

					ABMinput[1] = ynew;
				}
				else{
					// Store the previously computed solutions
					for(d=1; d<=ndim; d++){
						ABMinput[1][d] = sol[i-1][d]; // x[1] = {t_{i}, x_{i}, y_{i}, z_{i}, ....}
						ABMinput[2][d] = sol[i-2][d]; // x[2] = {t_{i-1}, x_{i-1}, y_{i-1}, z_{i-1}, ....}
						ABMinput[3][d] = sol[i-3][d]; // x[3] = {t_{i-2}, x_{i-2}, y_{i-2}, z_{i-2}, ....}
						ABMinput[4][d] = sol[i-4][d]; // x[4] = {t_{i-3}, x_{i-3}, y_{i-3}, z_{i-3}, ....}
					}	

					ynew = AB_Step(ABMinput); // Compute the solution using the AB predictor

					// Update ABMinput
					for(d=1; d<=ndim; d++){
						ABMinput[1][d] = ynew[d]; // x[1] = {t_{i+1}, x_{i+1}, y_{i+1}, z_{i+1}, ....}
						ABMinput[2][d] = sol[i-1][d]; // x[2] = {t_{i}, x_{i}, y_{i}, z_{i}, ....}
						ABMinput[3][d] = sol[i-2][d]; // x[3] = {t_{i-1}, x_{i-1}, y_{i-1}, z_{i-1}, ....}
						ABMinput[4][d] = sol[i-3][d]; // x[4] = {t_{i-2}, x_{i-2}, y_{i-2}, z_{i-2}, ....}
					}

					ynew = AM_Step(ABMinput); // Update the solution using the AM corrector
				}

				for(d=1; d<=ndim; d++){
					sol[i][d] = ynew[d]; 
				}
			}

			delete[] ABMinput;		
		}
		else{
			std::string reason = "Error: void IVP_Solver::ABM_PC_Solve()\n"; 
			if(!c1) reason += "Parameter order = " + template_funcs::toString(order) + " has not been assigned correctly\n"; 
			if(!c2) reason += "Parameter ndim = " + template_funcs::toString(ndim) + " has not been assigned correctly\n"; 
			if(!c3) reason += "Parameter N = " + template_funcs::toString(N) + " has not been assigned correctly\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void IVP_Solver::ABM_PC_Solve()\n"; 
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

void IVP_Solver::output_solution(std::string filename)
{
	// output the computed solution to a particular file

	// ensure that filename has a non null-string value
	if(filename == ""){
		filename = "IVP_Solver_Output.txt"; 
	}

	std::ofstream write; 
	write.open(filename.c_str(), std::ios_base::out|std::ios_base::trunc);

	if(write.is_open()){
		
		for(int i=1; i<=N; i++){
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
		// Note on the difference between cerr, cout and clog
		// http://www.cplusplus.com/forum/beginner/72020/	
		std::cerr<<"Could not open "<<filename<<"\n";
	}
}

void IVP_Solver::clear()
{
	// delete the arrays associated with the object

	if(alpha != nullptr) delete[] alpha; 
	if(sol != nullptr) delete[] sol; 
	if(func_vec != nullptr) delete[] func_vec; 
}