#ifndef ATTACH_H
#define ATTACH_H

// Put all #include statements here
// R. Sheehan 8 - 3 - 2013

// Basic libraries required for interaction
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// definitions of constant pi
static const double p = (atan(1.0));
static const double Two_PI = (8.0*p);
static const double PI = (4.0*p);
static const double PI_2 = (2.0*p);
static const double PI_3 = ((4.0/3.0)*p);
static const double PI_4 = (p);
static const double PI_5 = ((4.0/5.0)*p);
static const double PI_6 = ((2.0/3.0)*p);

typedef double (*coordfunc)(double *x, int n); // this maps a vector x into a real number 

// Files needed to perform calculations
// Order in which they are included matters
#include "Templates.h" // access template functions through the namespace access specifier
#include "Useful.h"
#include "PD_Func.h"
#include "System_Solver.h"
#include "Solve_IVP.h"
#include "Solve_BVP.h"
#include "Examples.h"

#endif