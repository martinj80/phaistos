// test_optimization.cpp --- Code for testing optimization functions
// Copyright (C) 2008 Wouter Boomsma
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


#include "utils/optimize.h"

// Different ways to specify a function

// Normal function
double rosen(Vector_nD x) {
     return(100*((x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])) + ((1-x[0])*(1-x[0])));
}


// function object
class Rosen {
public:
     double operator()(Vector_nD x) {
	  return rosen(x);
     }
};


// function object that allows optimization directly on already existing variables.
class Rosen_noarg {
     
     // pointer to variables x
     double *x;
public:

     // Constructor
     Rosen_noarg(double *x) {
	  this->x = x;
     }

     // Function evaluation (evaluates variable x)
     double operator()() {
	  return(100*((x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])) + ((1-x[0])*(1-x[0])));
     }
};



int main() {
     Powell powell;

     std::cout << "Optimizing Rosen function...\n";
     powell.optimize(&rosen, Vector_nD(2, 0.0, 0.0));
     std::cout << powell << "\n\n";

     std::cout << "Optimizing Rosen function object...\n";
     Rosen r;
     powell.optimize(&r, Vector_nD(2, 0.0, 0.0));
     std::cout << powell << "\n\n";

     std::cout << "Optimizing Rosen function in place...\n";
     std::cout << "Initial values of x:\n";
     double x[2] = {0.0, 0.0};	    // variables to be modified
     std::cout << "x[0]=" << x[0] << " " << "x[1]=" << x[1] << "\n";
     std::vector<double*> start;    // start values 
     start.push_back(&x[0]);
     start.push_back(&x[1]);
     Rosen_noarg rosen_noarg(x);    // Initialize special function object with variable x 
     powell.optimize(&rosen_noarg, start);
     std::cout << powell << "\n";
     std::cout << "Final values of x:\n";
     std::cout << "x[0]=" << x[0] << " " << "x[1]=" << x[1] << "\n\n";
}
