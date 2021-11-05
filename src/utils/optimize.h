// optimize.h --- Brent and Powell function minimizers
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

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "vector_nd.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <valarray>
#include <cstdarg>
#include <limits>

namespace phaistos {

//! Function base class.
//! Supports various ways of function calls from optimize functions
template <typename TYPE, typename FUNCTIONTYPE>
class Function {
public:
     //! Overload () operator
     //!
     //! \param x function argument
     //! \return function evaluation
     virtual double operator()(TYPE x) =0;

     //! Desctructor
     virtual ~Function(){};
};

//! Value Function
//!
//! \tparam TYPE Function argument type
//! \tparam FUNCTIONTYPE Function object
template <typename TYPE, typename FUNCTIONTYPE>
class ValueFunction: public Function<TYPE,FUNCTIONTYPE> {
private:

     //! Function object
     FUNCTIONTYPE func;
public:

     //! Constructor
     //!
     //! \param func Function object
     ValueFunction(FUNCTIONTYPE func) {
          this->func = func;          
     }
     
     //! Overload () operator
     //!
     //! \param x function argument
     //! \return function evaluation
     double operator()(TYPE x) {
          return (*func)(x);
     }
};

//! Pointer Function. Used with function that works on an existing
//! array of input values without taking arguments.
//!
//! \tparam TYPE Function argument type
//! \tparam FUNCTIONTYPE Function object
template <typename TYPE, typename FUNCTIONTYPE>
class PointerFunction: public Function<TYPE,FUNCTIONTYPE> {
private:

     //! Function object
     FUNCTIONTYPE func;

     //! Pointer to input value representation
     double **x_pointer;

     //! Size of input value vector
     int size;
public:

     //! Constructor
     //!
     //! \param func Function object
     //! \param size Input array size
     //! \param x_pointer Input array
     PointerFunction(FUNCTIONTYPE func, int size, double **x_pointer) {
          this->func = func;
          this->size = size;
          this->x_pointer = x_pointer;
     }

     //! Overload () operator - single value
     //!
     //! \param x function argument
     //! \return function evaluation     
     double operator()(double x) {
          **x_pointer = x;
          return (*func)();
     }

     //! Overload () operator - vector of input values
     //!
     //! \param x function argument
     //! \return function evaluation     
     double operator()(Vector_nD x) {
          for (int i=0; i<size; i++) {
               *x_pointer[i] = x[i]; 
          }
          return (*func)();
     }
};

     
//! Bracket.
//! Find points a, b, c that bracket the minimum of a given function
//!
//! \tparam TYPE Function argument type
template <typename TYPE>
class Bracket {
private:
     //@{
     //! Constants
     double gold_step;
     double grow_limit;
     double tiny_number;
     int max_iterations;
     //@}
     
     //! If b>=0, return fabs(a), otherwise return -fabs(a)
     inline double sign(double a, double b) {
          return ((b>=0.0) ? std::fabs(a) : -std::fabs(a));
     }

     //! Length function on value 
     static double len(double v) {return v;}

     //! Length function on vector
     static double len(Vector_nD v) {return v.norm();}

public:
     
     //@{
     //! Three bracketing x values a<b<c or c<b<a
     TYPE a, b, c;
     //@}

     //@{
     //! Corresponding function values
     double fa, fb, fc;
     //@}

     //! Error status
     bool status;

     //! Number of function evaluations     
     int evaluations;

     //! Constructor
     //!
     //! \param f Function pointer
     //! \param a bracketing x value
     //! \param b bracketing x value
     Bracket(double(*f)(double), double a=0.0, double b=1.0) {
          init(a, b);
          create_bracket(f);
     }

     //! Constructor
     //!
     //! \param f Function object
     //! \param a bracketing x value
     //! \param b bracketing x value
     //! \param c bracketing x value
     template <typename FUNCTIONTYPE>
     Bracket(Function<TYPE, FUNCTIONTYPE> f, double a, double b, double c) {
          init(a, b);
          this->a = a; this->fa = f(a);
          this->b = b; this->fb = f(b);
          this->c = c; this->fc = f(c);
          this->status = true;
          this->evaluations = 3;
     }

     //! Constructor
     //!
     //! \param func General functor object
     //! \param a bracketing x value
     //! \param b bracketing x value
     //! \param c bracketing x value
     template <typename FUNCTIONTYPE>
     Bracket(FUNCTIONTYPE func, double a, double b, double c) {
          ValueFunction<double, FUNCTIONTYPE> f(func);          
          init(a, b);
          this->a = a; this->fa = f(a);
          this->b = b; this->fb = f(b);
          this->c = c; this->fc = f(c);
          this->status = true;
          this->evaluations = 3;
     }
     
     //! Constructor
     //!
     //! \param f Function object specialized for double argument value
     //! \param a bracketing x value
     //! \param b bracketing x value
     template <typename FUNCTIONTYPE>
     Bracket(Function<double, FUNCTIONTYPE> &f, TYPE a=-1.0, TYPE b=1.0) {
          init(a, b);
          create_bracket(f);
     }

     //! Constructor
     //!
     //! \param f Function object specialized for Vector_nD argument value
     //! \param a bracketing x value
     //! \param b bracketing x value
     template <typename FUNCTIONTYPE>
     Bracket(Function<Vector_nD, FUNCTIONTYPE> &f, Vector_nD a=Vector_nD(), Vector_nD b=Vector_nD()) {
          init(a, b);
          create_bracket(f);
     }

     //! Initializer
     //!
     //! \param a bracketing value
     //! \param b bracketing value
     void init(TYPE &a=TYPE(), TYPE &b=TYPE()) {
          // Constants
          this->gold_step = 1.618034;
          this->grow_limit = 110.0;
          this->tiny_number = 1e-21;
          this->max_iterations = 50;
          
          this->a = a;
          this->b = b;
     
          this->status = true;
          this->evaluations = 0;
     }

     //! Construct a bracket
     //!
     //! \param func Function object
     template <typename FUNCTIONTYPE>     
     void create_bracket(FUNCTIONTYPE &func) {
          // std::cout << "********************************Bracketing\n";
          TYPE unit_vector = ((b-a)/(len(b-a)));
          TYPE origin = a;
     
          fa = func(a); evaluations++;
          fb = func(b); evaluations++;
     
          if (fb > fa) {
               // Switch roles so f(a)>f(b) - downhill from a to b
               TYPE tmp = a;
               a = b;
               b = tmp;
               double ftmp = fa;
               fa = fb;
               fb = ftmp;
          }

          double ax = len(a-origin);
          double bx = len(b-origin);
          double cx;

          double xmin = ax;
          double fmin = fa;
          
          // Guess a value
          cx = bx + (bx-ax)*gold_step;
          fc = func(origin + unit_vector * cx); evaluations++;
          
          int iterations=0;
          while (fb >= fc) {
               /* std::cout << ax << " " << bx << " " << cx << "\n"; */
               /* std::cout << fa << " " << fb << " " << fc << "\n"; */
               if (iterations > max_iterations) {
                    status = false;
                    a = origin + unit_vector * xmin;
                    fa = fmin;
                    return;
               }
               iterations++;

               double r = (bx-ax)*(fb-fc);
               double q = (bx-cx)*(fb-fa);
               double u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(fmax(std::fabs(q-r), tiny_number), q-r));
               double ulim = bx+(cx-bx)*grow_limit;
               double fu;
          
               if ((bx-u)*(u-cx) > 0.0) {
                    // Parabolic u is between b and c.
                    // std::cout << "Parabolic u is between b and c.\n";
                    fu=func(origin + unit_vector * u); evaluations++;
                    if (fu < fc) {
                         // std::cout << "Minimum is between b and c.\n";
                         // Minimum between b and c
                         ax=bx;
                         bx=u;
                         fa=fb;
                         fb=fu;
                         break;
                    } else if (fu > fb) {
                         // std::cout << "Minimum is between a and u.\n";
                         // Minimum between a and u
                         cx=u;
                         fc=fu;
                         break;
                    } 
                    // std::cout << "Parabolic fit of no use - use default magnification\n";
                    // Parabolic fit of no use - use default magnification
                    u=cx+(cx-bx)*gold_step;
                    fu=func(origin + unit_vector * u); evaluations++;
               } else if ((cx-u)*(u-ulim) > 0.0) {
                    // Parabolic fit between c and allowed limit
                    // std::cout << "Parabolic fit between c and allowed limit\n";
                    fu=func(origin + unit_vector * u); evaluations++;
                    if (fu < fc) { 
                         bx = cx;
                         cx = u;
                         u = cx + (cx-bx)*gold_step;
                         fb = fc;
                         fc = fu;
                         fu = func(origin + unit_vector * u); evaluations++;

                    }
               } else if ((u-ulim)*(ulim-cx) >= 0.0) {
                    // set u to limit
                    // std::cout << "set u to limit\n";
                    u=ulim;
                    fu=func(origin + unit_vector * u); evaluations++;
               } else {
                    // Reject parabolic fit - use default magnification
                    // std::cout << "Reject parabolic fit - use default magnification\n";
                    u=cx+(cx-bx)*gold_step;
                    fu=func(origin + unit_vector * u); evaluations++;
               }

               if (fu < fmin) {
                    fmin = fu;
                    xmin = u;
               }
               
               if (fu != fc) {   // Only update a if f(u) is different thatn f(c)
                    ax = bx;
                    fa = fb;
               }
               
               bx = cx;
               fb = fc;

               cx = u;
               fc = fu;
          }

          // Swap so a<b<c
          if (ax > cx) {
               c = origin + unit_vector * ax;
               b = origin + unit_vector * bx;
               a = origin + unit_vector * cx;
               double tmp = fa;
               fa = fc;
               fc = tmp;
          } else {
               a = origin + unit_vector * ax;
               b = origin + unit_vector * bx;
               c = origin + unit_vector * cx;
          }

          if (std::fabs(fa - fb) < tiny_number) { // If fa==fb call create bracket to get other boundary right
               TYPE tmpx = a;
               double tmpf = fa;
               a = c;
               fa = fc;
               c = tmpx;
               fc = tmpf;
               create_bracket(func);
          }

          return;
     }

     //! Overload output operator     
     friend std::ostream& operator<<(std::ostream& o, const Bracket<TYPE> &b) {
          o << "bracket = (" << b.a << ", " << b.b << ", " << b.c << ")" << "\t(";
          o << "evaluations=" << b.evaluations << ", ";
          o << "Status=";
          if (b.status) {
               o << "Completed";
          } else {
               o << "Terminated abnormally";
          }
          o << ")";
          return o;
     }
};



//! Brent optimization.
//! Optimize function of one variable.
//! Finds minimum of function given a bracketing triplet (a, b, c) (see numerical recipes)
//! Returns minimum x value and optionally the corresponding function value
//! NOTE: The tolerance should be no smaller than the square root of the machines precision of double
template <typename TYPE>
class Brent {
private:
     //@{
     //! Constants
     double min_tol;
     double golden_ratio;
     double tiny_number;
     //@}

     //! Maximum number of iterations allowed
     int max_iterations;

     //! Tolerance in x
     double xtol;

     //! Fractional tolerance in the function value (used for termination)
     double ftol;

     //! Length function on value 
     static double len(double v) {return v;}

     //! Length function on vector     
     static double len(Vector_nD v) {return v.norm();}

public:

     //! minimal x-value      
     TYPE xmin;

     //! Minimal function value
     double fmin;

     //! Number of evaluations done
     int evaluations;

     //! error status
     bool status;

     //! Default constructor
     Brent() {};

     //! Initialization
     //!
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     //! \param max_iterations Maximum number of allowed iterations
     void init(double xtol, double ftol, int max_iterations) {
          // Constants
          this->min_tol = 1.0e-11;
          this->golden_ratio = 0.3819660;
          this->tiny_number = 1e-25;
          

          evaluations = 0;
          this->max_iterations = max_iterations;
          this->xtol = xtol;
          this->ftol = ftol;
          
          this->status = true;
          
          std::numeric_limits<double> l;
          if (this->xtol<0) {
               this->xtol = sqrt(l.epsilon());
          }
          if (this->ftol<0) {
               this->ftol = sqrt(l.epsilon());
          }
     }

     //! Optimize - given point and direction - for functions already wrapped
     //!
     //! \param func Function object
     //! \param point starting point
     //! \param direction Direction in which to move
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     //! \param max_iterations Maximum number of allowed iterations
     template <typename FUNCTIONTYPE>
     void optimize(Function<TYPE, FUNCTIONTYPE> &func, TYPE point, TYPE direction, double xtol=-1, double ftol=-1, int max_iterations = 500) {
          compute(func, Bracket<TYPE>(func, point, point+(direction/len(direction))), xtol, max_iterations);
     }


     //! Optimize - given bracket
     //!
     //! \param func Function object
     //! \param bracket Bracket object
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     //! \param max_iterations Maximum number of allowed iterations
     template <typename FUNCTIONTYPE>
     void optimize(FUNCTIONTYPE func, Bracket<TYPE> bracket, double xtol=-1, double ftol=-1, int max_iterations = 500) {
          ValueFunction<TYPE, FUNCTIONTYPE> f(func);
          compute(f, bracket, xtol, ftol, max_iterations);
     }

     //! Optimize - default bracketing
     //!
     //! \param func Function object
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     //! \param max_iterations Maximum number of allowed iterations
     template <typename FUNCTIONTYPE>
     void optimize(FUNCTIONTYPE func, double xtol=-1, double ftol=-1, int max_iterations = 500) {
          ValueFunction<TYPE, FUNCTIONTYPE> f(func);
          compute(f, Bracket<TYPE>(f), xtol, ftol, max_iterations);
     }

     //! Compute - main method
     //!
     //! \param func Function object
     //! \param bracket Bracket object
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     //! \param max_iterations Maximum number of allowed iterations
     template <typename FUNCTIONTYPE>
     void compute(Function<TYPE, FUNCTIONTYPE> &func, Bracket<TYPE> bracket, double xtol=-1, double ftol=-1, int max_iterations = 500) {
          // std::cout << "********************************Brent\n";
          double a, b;                // Minimum is bracketed between a and b
          double x;                   // x-value corresponding to least function value so far
          double w;                   // x-value corresponding to next-to-least function value so far
          double v;                   // Previous value of w
          double u;                   // x-value at which the function has been evaluated most recently
          double xm;                  // midpoint between a and b (function is not evaluated here)
          double fx, fw, fv, fu;      // function values
          double xtol1, xtol2;        // Tolerance values
          double d = 0.0;             // Step-size
          double e = 0.0;             // Step-size of step i-2

          init(xtol, ftol, max_iterations);

          if (bracket.status == false) {
               status=false;
               std::cerr << "Unconstrained bracket - no optimization done\n";
               return;
          }

          TYPE unit_vector = ((bracket.b-bracket.a)/(len(bracket.b-bracket.a)));
          TYPE origin = bracket.a;
          
          // Ensure that a and b are in ascending order
          if ((bracket.a-origin) < (bracket.c-origin)) {
               a = len(bracket.a-origin);
               b = len(bracket.c-origin);
          } else {
               a = len(bracket.c-origin);
               b = len(bracket.a-origin);
          }

//           std::cout << bracket;
          
          // Initialize
          x=w=v=len(bracket.b-origin);
          fw=fv=fx=func(origin + unit_vector*x);
          evaluations++;


          int i;
          for (i=0; i<max_iterations; i++) {

               xtol1 = this->xtol*std::fabs(x) + min_tol;
               xtol2 = 2*xtol1;

               xm = 0.5 * (a+b);

               // Termination test (convergence)
               if (std::fabs(x-xm) <= (xtol2 - 0.5*(b-a))) {
                    xmin = origin + unit_vector * x;
                    fmin = fx;
                    return;
               }

               
               if (std::fabs(e) <= xtol1) {
                    // Golden step
                    if (x >= xm) {
                         e = a - x;
                    } else {
                         e = b - x;
                    }
                    d = golden_ratio * e;
               } else {
                    // Parabolic fit
                    double r = (x-w)*(fx-fv);
                    double q = (x-v)*(fx-fw);
                    double p = (x-v)*q - (x-w)*r;
                    q = 2.0 * (q-r);
                    if (q > 0.0) {
                         p = -p;
                    }
                    q = std::fabs(q);
                    double etemp = e;
                    e = d;
                    if ((std::fabs(p) < std::fabs(0.5*q*etemp)) && (p > q*(a-x)) && (p < q*(b-x))) {
                         // Use parabolic step
                         d = p/q;
                         u = x + d;
                         if ((u-a) < xtol2 or (b-u) < xtol2) {
                              if ((xm-x) >= 0) {
                                   d = xtol1;
                              } else {
                                   d = -xtol1;
                              }
                         }
                    } else {
                         // Use golden step 
                         if (x >= xm) {
                              e = (a-x);
                         } else {
                              e = (b-x);
                         }
                         d = golden_ratio * e;
                    }
               }

               // Use xtol1 if d is too small
               if (std::fabs(d) >= this->xtol) {
                    u = x + d;
               } else {
                    if (d >= 0) {
                         u = x + xtol1;
                    } else {
                         u = x - xtol1;
                    }
               }
          
               fu = func(origin + unit_vector * u); evaluations++;         // The function evaluation

               if (fu <= fx) {
                    if (u >= x) {
                         a = x;
                    } else {
                         b = x;
                    }
                    v = w;
                    w = x;
                    x = u;
                    fv = fw;
                    fw = fx;
                    fx = fu;
               } else {
                    if (u < x) {
                         a = u;
                    } else {
                         b = u;
                    }
                    if ((fu <= fw) || (w==x)) {
                         v = w;
                         w = u;
                         fv = fw;
                         fw = fu;
                    } else if (fu <= fv || v==x || v==w ) {
                         v = u;
                         fv = fu;
                    }
               }
          }

          xmin = origin + unit_vector * x;
          fmin = fx;
          return;
     }
          
     
     //! Overload output operator     
     friend std::ostream& operator<<(std::ostream& o, const Brent &b) {
       o << "(";
          o << "x=" << b.xmin << ", ";
          o << "f(x)=" << b.fmin << ", ";
          o << "evaluations=" << b.evaluations << ", ";
          o << "Status=";
          if (b.status) {
               o << "Completed";
          } else {
               o << "Terminated abnormally";
          }
          o << ")";
          return o;
     }
};



//! Direction set used by powell
class DirectionSet {
public:
     //! Internal data
     Vector_nD *data;

     //! Size of set
     int size;

     //! Constructor - Create direction set from unit vectors
     //!
     //! \param size Size of set
     DirectionSet(int size) {
          init(size);
          init_as_unit_vectors(size);
     }

     // Constructor - Create direction set from vector list
     //!
     //! \param size Size of set
     //! \param vector_list List of vectors
     DirectionSet(int size, Vector_nD *vector_list) {
          init(size);
          for (int i=0; i<size; i++) {
               data[i] = vector_list[i];
          }
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     DirectionSet(const DirectionSet &other) {
          init(other.size);
          for (int i=0; i<this->size; i++) {
               this->data[i] = other.data[i];
          }
     }

     //! Destructor
     ~DirectionSet() {
          delete[] this->data;
     }

     //! Initialization
     //!
     //! \param size Size of set
     void init(int size) {
          this->size = size;
          data = new Vector_nD[size];
     }
          
     //! Create unit vectors
     //!
     //! \param size Size of set
     void init_as_unit_vectors(int size) {
          for (int i=0; i<size; i++) {
               Vector_nD v(size);
               for (int j=0; j<size; j++) {
                    if (i==j) {
                         v[j] = 1;
                    } else {
                         v[j] = 0;
                    }
               }
               data[i] = v;
          }
     }
     
     //! Overload [] operator
     //!
     //! \param index Index into set
     //! \return direction vector
     Vector_nD operator[](const int index) const {
          return data[index];
     }

     //! Overload [] operator - non const
     //!
     //! \param index Index into set
     //! \return direction vector
     Vector_nD& operator[](const int index) {
          return data[index];
     }
     
     //! Overload output operator     
     friend std::ostream& operator<<(std::ostream& o, const DirectionSet &d) {
          for (int i=0; i<d.size; i++){
               std::cout << d.data[i] << "\n";
          }
          return o;
     }
};


//! Modified Powell's method
//! Minimizes function of n variables given a starting point and an initial direction matrix (defaults to unit vectors)
class Powell {
private:

     // Brent optimizer object
     Brent<Vector_nD> brent_optimizer;

     //! Tolerance constant 
     double tiny_number;

public:

     //! Minimum x-value      
     Vector_nD xmin;

     //! Minimum f(x) value     
     double fmin;

     //! Change in function value due to optimization
     double deltaf;

     //! Number of executed evaluations
     int evaluations;

     //! Number of executed iterations (an iteration can consist of several evaluations)
     int iterations;

     //! x-value tolerance
     double xtol;

     //! Fractional tolerance in the function value (used for termination)
     double ftol;
     
     //! error status
     bool status;

     //! Default constructor
     Powell() {
          this->init();
     };

     //! Initializer
     void init() {
          this->tiny_number = 1e-25;
     }     

     //! Optimize - Main method
     //!
     //! \param func Function object
     //! \param starting_point Initial x value
     //! \param directions Set of directions
     //! \param max_iterations Maximum number of allowed iterations
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     template <typename FUNCTIONTYPE>
     void compute(FUNCTIONTYPE &func, Vector_nD starting_point, DirectionSet directions,
                  int max_iterations=1000, double xtol=-1, double ftol=-1) {
          // std::cout << "********************************Powell\n";
          // std::cout.flush();
          evaluations = 0;
          this->status = true;

          xmin.clear();
          
          // Initialize tolerances
          this->xtol = xtol;
          this->ftol = ftol;
          std::numeric_limits<double> l;
          if (this->xtol<0) {
               this->xtol = sqrt(l.epsilon());
          }
          if (this->ftol<0) {
               this->ftol = sqrt(l.epsilon());
          }

          assert(directions.size == starting_point.size);

          xmin = starting_point;
          fmin = func(xmin); evaluations++;
          deltaf = fmin;

          for (iterations=1; ; iterations++) {
               Vector_nD x_old = xmin;                  // Keep track of previous values 
               double fx_old = fmin;

               int i_largest = -1;                // index corresponding to direction with largest decrease
               double delta_f_largest = 0.0;       // function decrease corresponding to iLargest

               // Call brent for all directions
               for (int i=0; i<directions.size; i++) {
                    Vector_nD direction = directions[i];
                    brent_optimizer.optimize(func, xmin, direction, xtol);
                    if (brent_optimizer.status == true) {
                         if ((fmin - brent_optimizer.fmin) > delta_f_largest) {
                              i_largest = i;
                              delta_f_largest = fmin - brent_optimizer.fmin;
                         }
                         xmin = brent_optimizer.xmin;     // Update x, fx and evaluations
                         fmin = brent_optimizer.fmin;
                    }
                    evaluations += brent_optimizer.evaluations;
               }

               // Termination criterium
               if (2.0 * (fx_old-fmin) < this->ftol*(std::fabs(fx_old)+std::fabs(fmin))+tiny_number) {
                    break;
               }

               // Maximum number of iterations exceeded
               if (iterations >= max_iterations) {
                    status = false;
                    break;
               }

               // Find extrapolated point based on Pn-P0 direction
               Vector_nD extrapolated = xmin*2.0 - x_old;
               double f_extrapolated = func(extrapolated); evaluations++;
               Vector_nD direction = xmin-x_old;

               // Check whether it is better to keep old direction set
               if (f_extrapolated < fx_old) {
                    if ((2.0*(fx_old-2.0*fmin+f_extrapolated) *
                         (fx_old - fmin - delta_f_largest) * (fx_old - fmin - delta_f_largest)) < 
                        (delta_f_largest * (fx_old - f_extrapolated)*(fx_old - f_extrapolated))) {
                         
                         brent_optimizer.optimize(func, xmin, direction, xtol); // Call brent on new direction
                         if (brent_optimizer.status == true) {
                              xmin = brent_optimizer.xmin; // This is different than in NR - They save the old value of x before doing this step (why?)
                              fmin = brent_optimizer.fmin;
                              directions[i_largest] = directions[directions.size-1]; // Add new direction to direction set
                              directions[directions.size-1] = direction;
                         }
                         evaluations += brent_optimizer.evaluations;
                    }
               }
          }
          deltaf = deltaf - fmin;
     }

     
     //! Optimize - Main method - automatically find direction set
     //!
     //! \param func Function object
     //! \param starting_point Initial x value
     //! \param max_iterations Maximum number of allowed iterations
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     template <typename FUNCTIONTYPE>
     void optimize(FUNCTIONTYPE func, Vector_nD starting_point,
                   int max_iterations=1000, double xtol=-1, double ftol=-1) {
          ValueFunction<Vector_nD, FUNCTIONTYPE> f(func);
          compute(f, starting_point, DirectionSet(starting_point.size), max_iterations, xtol, ftol);
     }
     

     //! Optimize - Main method
     //!
     //! \param func Function object
     //! \param starting_point Initial x value
     //! \param directions Set of directions
     //! \param max_iterations Maximum number of allowed iterations
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     template <typename FUNCTIONTYPE>
     void optimize(FUNCTIONTYPE func, std::vector<double *> starting_point, DirectionSet directions,
                   int max_iterations=1000, double xtol=-1, double ftol=-1) {
          std::cout << "Starting Powell optimization\n";
          std::cout.flush();
          double **x_pointers = new double*[starting_point.size()];
          for (unsigned int i=0; i<starting_point.size(); i++) {
               x_pointers[i] = starting_point[i];
          }
          PointerFunction<Vector_nD, FUNCTIONTYPE> f(func, starting_point.size(), x_pointers);
          compute(f, Vector_nD(starting_point), directions, max_iterations, xtol, ftol);
          for (unsigned int i=0; i<starting_point.size(); i++) {
               *x_pointers[i] = xmin[i];
          }
          delete[] x_pointers;
     }

     //! Optimize - Main method - automatically find direction set
     //!
     //! \param func Function object
     //! \param starting_point Initial x value
     //! \param max_iterations Maximum number of allowed iterations
     //! \param xtol Tolerance in x
     //! \param ftol Fractional tolerance in the function value (used for termination)
     template <typename FUNCTIONTYPE>
     void optimize(FUNCTIONTYPE func, std::vector<double *> starting_point,
                   int max_iterations=1000, double xtol=-1, double ftol=-1) {
          optimize(func, starting_point, DirectionSet(2), max_iterations, xtol, ftol);
     }
     
     
     //! Overload output operator     
     friend std::ostream& operator<<(std::ostream& o, const Powell &p) {
          o << "(";
          o << "xmin=" << p.xmin << ", ";
          o << "f(xmin)=" << p.fmin << ", ";
          o << "delta_f=" << p.deltaf << ", ";
          o << "iterations=" << p.iterations << ", ";
          o << "evaluations=" << p.evaluations << ", ";
          o << "Status=";
          if (p.status) {
               o << "Completed";
          } else {
               o << "Terminated abnormally";
          }
          o << ")";
          return o;
     }
};

}

#endif
