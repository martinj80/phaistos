// integrate.cpp --- Wrapper for fortran integration code
// Copyright (C) 2006-2008 Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.


#include "integrate.h"

// Computes a definite integral
// Integrate func from a to b (possibly infinite interval) using a technique
// from the Fortran library QUADPACK
double integrate_quad(double(*func)(double *), double a, double b, int *evaluations,
                      double epsabs, double epsrel, int limit) {

     double res;
     double abserr;
     int ier;
  
     int *iwork = new int[limit];
     int lenw = limit*4;
     int last;
     double *work = new double[lenw];
     int neval;
  
     dqags_(func, &a, &b, &epsabs, &epsrel, &res, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  
     delete[] work;
     delete[] iwork;

     if (evaluations) {
          *evaluations = neval;
     }
     return res;
}

// Global variables necessary to allow use to pass additional arguments with function
double (*quad_func)(double *, double *);
void *quad_extra_args = NULL;

// Single-argument funciton wrapper
// Calls the original function including additional arguments
double _quad_single_argument_function(double *x) {
     return quad_func(x, (double *)quad_extra_args);
}

// Version of quad taking additional arguments of type double
double integrate_quad(double(*func)(double *, double *), double extra_arguments[], double a, double b, int *evaluations,
            double epsabs, double epsrel, int limit) {
     quad_func = func;                           // Original function pointer is saved
     quad_extra_args = (void *)extra_arguments;            // Additional argument values are saved

     // quad is called on single-argument-wrapper
     return integrate_quad(&_quad_single_argument_function, a, b, evaluations, epsabs, epsrel, limit);
}



extern "C" int MAIN__() { return 1; }
