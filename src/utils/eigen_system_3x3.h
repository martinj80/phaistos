// eigen_system_3x3.h --- Efficient calculation of eigen vectors using analytical solution
//                 Based on: David Eberly
//                           Eigensystems for 3x3 symmetric matrices (revisited)
// Copyright (C) 2011 Wouter Boomsma
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


#include "vector_matrix_3d.h"
#include "math.h"

#ifndef EIGEN_SYSTEM_3X3_H
#define EIGEN_SYSTEM_3X3_H

namespace phaistos {

//! Calculation of eigen vectors and values for 3x3 matrices
//! Based on: David Eberly
//!           Eigensystems for 3x3 symmetric matrices (revisited)
class EigenSystem3x3 {

     //! Saved for efficiency reasons: inv_3 = 1.0/3.0
     //! NOTE: used by compute_roots, which requires double precision 
     static const double inv_3;

     //! Saved for efficiency reasons: sqrt_3 = sqrt(3.0)
     //! NOTE: used by compute_roots, which requires double precision 
     static const double sqrt_3;

public:

     //! Eigen values
     Vector_3D eigen_values;

     //! Eigen vectors
     Vector_3D eigen_vectors[3];

     //! From the original description:
     //! The maximum-magnitude entries for all three of the M_i = A - lambda_i *I are computed. If all three maxima are
     //! smaller than Mathf::ZERO TOLERANCE, the three eigenvalues are assumed to be the same and the eigenvectors
     //! are (1; 0; 0), (0; 1; 0), and (0; 0; 1). If at least one maximum-magnitude entry is larger than the tolerance, the
     //! matrix producing the maximum of these is processed further to construct the eigenvectors of A. The function
     //! PositiveRank returns true when the matrix has positive rank and also returns the maximum-magnitude
     //! entry and the row of the matrix in which this entry occurs.
     //!
     //! \param A Input matrix
     EigenSystem3x3(const Matrix_3D &A) {

          eigen_values.initialized = true;

          Matrix_3D A_scaled(A);

          // Find maximum value (A is symmetric so only 6 values are checked)
          double max_value = std::fabs(A(0,0));          
          double value;
          value = std::fabs(A(0,1));
          if (value > max_value)
               max_value = value;
          value = std::fabs(A(0,2));
          if (value > max_value)
               max_value = value;
          value = std::fabs(A(1,1));
          if (value > max_value)
               max_value = value;
          value = std::fabs(A(1,2));
          if (value > max_value)
               max_value = value;
          value = std::fabs(A(2,2));
          if (value > max_value)
               max_value = value;

          if (max_value > 1.0) {
               double max_value_inverse = 1.0/max_value;
               A_scaled *= max_value_inverse;
          }

          // Calculate roots
          double roots[3];
          compute_roots(A_scaled, roots);

          eigen_values[0] = roots[0];
          eigen_values[1] = roots[1];
          eigen_values[2] = roots[2];

          double max_entry[3];
          Vector_3D max_row[3];
          for (unsigned int i=0; i<3; ++i) {
               // M = A - eigen_values[i]*I
               Matrix_3D M = A_scaled;
               M(0,0) -= eigen_values[i];
               M(1,1) -= eigen_values[i];
               M(2,2) -= eigen_values[i];

               // Check whether matrix has zero rank. If so, all eigenvalues are assumed
               // to have the same value, and standard basis vectors can be used
               // as eigen vectors (spherical symmetry)
               if (!positive_rank(M, max_entry[i], max_row[i])) {

                    // If scaling was done, rescale eigenvalues back to original size
                    if (max_value > 1.0) {
                         eigen_values *= max_value;
                    }

                    eigen_vectors[0] = Vector_3D(1,0,0);
                    eigen_vectors[1] = Vector_3D(0,1,0);
                    eigen_vectors[2] = Vector_3D(0,0,1);

                    return;
               }
          }

          double total_max = max_entry[0];
          int i=0;
          if (max_entry[1] > total_max) {
               total_max = max_entry[1]; i=1;
          }
          if (max_entry[2] > total_max) {
               total_max = max_entry[2]; i=2;
          }

          if (i==0) {
               compute_vectors(A_scaled, max_row[0], 1, 2, 0);
          } else if (i==1) {
               compute_vectors(A_scaled, max_row[1], 2, 0, 1);
          } else {
               compute_vectors(A_scaled, max_row[2], 0, 1, 2);
          }

          // If scaling was done, rescale eigenvalues back to original size
          if (max_value > 1.0) {
               eigen_values *= max_value;
          }

          // double error = 0.0;
          // for (unsigned int i=0; i<3; ++i) {
          //      Vector_3D result = A*eigen_vectors[i] - eigen_values[i]*eigen_vectors[i];
          //      double length = result.norm();
          //      if (length > error) {
          //           error = length;
          //      }
          // }
          // std::cout << "Error: " << error << "\n";
     }

     //! Calculate the roots of a the cubic characteristic polynomium
     //! Roots are returned in increasing order
     //!
     //! \param A input matrix
     //! \param roots Output: roots in increasing order
     inline void compute_roots(const Matrix_3D &A, double roots[3]) {
          
          // The characteristic equation is 0 = lambda^3 - c_2*lambda^2 + c_1*lambda - c0, with
          double c0 = (A(0,0)*A(1,1)*A(2,2) + 2.0*A(0,1)*A(0,2)*A(1,2) - 
                       A(0,0)*A(1,2)*A(1,2) - A(1,1)*A(0,2)*A(0,2) - A(2,2)*A(0,1)*A(0,1));
          double c1 = (A(0,0)*A(1,1) - A(0,1)*A(0,1) + A(0,0)*A(2,2) - 
                       A(0,2)*A(0,2) + A(1,1)*A(2,2) - A(1,2)*A(1,2));
          double c2 = A(0,0) + A(1,1) + A(2,2);

          double c2_div_3 = c2*inv_3;
          double a_div_3 = (c1 - (c2*c2_div_3))*inv_3;

          // a is either 0 or negative
          if (a_div_3 > 0.0)
               a_div_3 = 0.0;

          double half_neg_b = 0.5*(c0 + c2_div_3*(2.0*c2_div_3*c2_div_3 - c1));

          double q = half_neg_b*half_neg_b + a_div_3*a_div_3*a_div_3;

          // q is eigher 0 or negative
          if (q > 0.0)
               q = 0.0;

          // Find roots
          double angle = atan2(sqrt(-q), half_neg_b)*inv_3;
          double cos_angle = cos(angle);
          double sin_angle = sin(angle);
          double rho = sqrt(-a_div_3);
          double root0 = c2_div_3 + 2.0*rho*cos_angle;
          double root1 = c2_div_3 - rho*(cos_angle + sqrt_3*sin_angle);
          double root2 = c2_div_3 - rho*(cos_angle - sqrt_3*sin_angle);

          // Sort roots
          if (root1 >= root0) {
               roots[0] = root0;
               roots[1] = root1;
          } else {
               roots[0] = root1;
               roots[1] = root0;
          }
          if (root2 >= roots[1]) {
               roots[2] = root2;
          } else {
               roots[2] = roots[1];
               if (root2 >= roots[0]) {
                    roots[1] = root2;
               } else {
                    roots[1] = roots[0];
                    roots[0] = root2;
               }
          }
     }

     //! Calculate the three eigen vectors based on a primary direction u2
     //!
     //! \param A Input matrix
     //! \param u2 Primary direction
     //! \param i0 Index of vector 0 in eigen_vectors matrix
     //! \param i1 Index of vector 1 in eigen_vectors matrix
     //! \param i2 Index of vector 2 in eigen_vectors matrix
     inline void compute_vectors(const Matrix_3D &A, Vector_3D &u2, int i0, int i1, int i2) {

          // Construct an orthonormal basis {u0, u1, u2}
          u2 /= u2.norm();
          Vector_3D u0, u1;
          // For numerical stability, we use the largest coordinates
          if (std::fabs(u2[0]) >= std::fabs(u2[1])) {

               // x-coordinate or z-coordinate are largest - project to that plane
               double inv_length = 1.0/sqrt(u2[0]*u2[0]+u2[2]*u2[2]);
               u0[0] = -u2[2]*inv_length;
               u0[1] = 0.0;
               u0[2] = u2[0]*inv_length;

               // Cross product - saves a few calculations to the the 0 component of u
               u1[0] = u2[1]*u0[2];
               u1[1] = u2[2]*u0[0] - u2[0]*u0[2];
               u1[2] = -u2[1]*u0[0];

          } else {

               // y-coordinate or z-coordinate are largest - project to that plane
               double inv_length = 1.0/sqrt(u2[1]*u2[1]+u2[2]*u2[2]);
               u0[0] = 0.0;
               u0[1] = u2[2]*inv_length;
               u0[2] = -u2[1]*inv_length;

               // Cross product - saves a few calculations due to the the 0 component of u
               u1[0] = u2[1]*u0[2] - u2[2]*u0[1];
               u1[1] = -u2[0]*u0[2];
               u1[2] = u2[0]*u0[1];
          }

          
          // The eigenvector v2 is perpendicular to u2 (M2*V2=0), therefore
          // v2 = c0*u0+c1*u1 [8], where c0=v2*u0 and c1=v2*u1. Since u0, u1 and u2 
          // are of unit length, c0^2+c1^2=1.
          // we multiply [8] with A
          //   A*v2 = lambda2*V2 = c0*A*u0 + c1*A*u1
          // Dotting this with u0 and u1, we get two equations
          //   u0^T*A*u0*c0 + u0^T*A*u1*c1 = lambda2*c0
          //   u0^T*A*u1*c0 + u1^T*A*u1*c1 = lambda2*c1
          // Subtracting the right hand side from the left
          //   (u0^T*A*u0 - lambda2)*c0 + u0^T*A*u1*c1 = 0
          //   u0^T*A*u1*c0 + (u1^T*A*u1 - lambda2)*c1 = 0

          Vector_3D Au0 = A*u0;

          // p00 = (lambda2 - u0^T*A*u0)
          double p00 = eigen_values[i2] - u0*Au0;
          // p01 = u1^T*Au0 = u0^T*Au1 (A is symmetric)
          double p01 = u1*Au0;
          // p11 = (lambda2 - u1^T*A*u1)
          double p11 = eigen_values[i2] - u1*(A*u1);

          // The system now becomes:
          //     p00*c0+p01*c1 = 0
          //     p01*c0+p11*c1 = 0
          // In addition, we have the relationship c0^2+c1^2=1. This means
          // that it is sufficient to solve one of the equations. We choose
          // the one which has the largest coefficient.
          double max_value = std::fabs(p00);
          int max_value_row = 0;

          double value = std::fabs(p01);
          if (value > max_value) {
               max_value = value;
          }
          value = std::fabs(p11);
          if (value > max_value) {
               max_value = value;
               max_value_row = 1;
          }

          if (max_value >= Math<double>::zero_tolerance) {

               if (max_value_row == 0) {

                    double inv_length = 1.0/sqrt(p00*p00+p01*p01);
                    eigen_vectors[i2] = (p01*inv_length)*u0 + (p00*inv_length)*u1;
                    // p00*=inv_length;
                    // p01*=inv_length;
                    // eigen_vectors[i2] = p01*u0 + p00*u1;
               } else {

                    double inv_length = 1.0/sqrt(p11*p11+p01*p01);
                    eigen_vectors[i2] = (p11*inv_length)*u0 + (p01*inv_length)*u1;
                    // p01*=inv_length;
                    // p11*=inv_length;
                    // eigen_vectors[i2] = p11*u0 + p01*u1;
               }

          } else  {

               // If the maximum entry is not above zero, any vectors in the plane
               // are valid. We choose u1 if the maximum element is in the first row
               // and u0 if the maximum element is in the second row
               if (max_value_row == 0) {
                    eigen_vectors[i2] = u1;
               } else {
                    eigen_vectors[i2] = u0;
               }
          }

          // We have found the eigenvector v2. The two other eigenvectors v0 and v1
          // lie in the plane perpendicular to v2. We know that u2 lies in the plane,
          // and construct a vector s as the final vector (S=u2xv2)
          Vector_3D s = cross_product(u2, eigen_vectors[i2]);

          // Using the same strategy as before, we have:
          //     v0 = c0*u2 + c1*S
          //     c0^2+c1^2 = 1
          // We multiply with A:
          //     A*v0 = v0*lambda0 = c0*A*u2 + c1*A*s
          // Dotting with u2 and s, we obtain two equations
          //     u2^T*A*u2*c0 + u2^T*A*s*c1 = lambda0*c0
          //     s^T*A*u2*c0 + s^T*A*s*c1 = lambda0*c1
          // For convenience:
          //     p00*c0+p01*c1 = 0
          //     p01*c0+p11*c1 = 0
          // with:
          //     p00 = u2^T*A*u2 - lambda0
          //     p01 = u2^T*A*s
          //     p11 = s^T*A*s - lambda0

          Vector_3D Au2 = A*u2;
          p00 = eigen_values[i0] - u2*Au2;
          p01 = s*Au2;
          p11 = eigen_values[i0] - s*(A*s);

          max_value = std::fabs(p00);
          max_value_row = 0;

          value = std::fabs(p01);
          if (value > max_value) {
               max_value = value;
          }
          value = std::fabs(p11);
          if (value > max_value) {
               max_value = value;
               max_value_row = 1;
          }

          if (max_value >= Math<double>::zero_tolerance) {

               if (max_value_row == 0) {
                    double inv_length = 1.0/sqrt(p00*p00+p01*p01);
                    // eigen_vectors[i0] = (p01*inv_length)*u2 + (p00*inv_length)*s;
                    p00 *= inv_length;
                    p01 *= inv_length;
                    eigen_vectors[i0] = p01*u2 + p00*s;
               } else {
                    double inv_length = 1.0/sqrt(p11*p11+p01*p01);
                    // eigen_vectors[i0] = (p11*inv_length)*u2 + (p01*inv_length)*s;
                    p01 *= inv_length;
                    p11 *= inv_length;
                    eigen_vectors[i0] = p11*u2 + p01*s;
               }

          } else  {
               // If the maximum entry is not above zero, any vectors in the plane
               // are valid. We choose u1 if the maximum element is in the first row
               // and u0 if the maximum element is in the second row
               if (max_value_row == 0) {
                    eigen_vectors[i0] = s;
               } else {
                    eigen_vectors[i0] = u2;
               }
          }

          // The last eigenvector is the cross product of the two others
          eigen_vectors[i1] = cross_product(eigen_vectors[i2], eigen_vectors[i0]);
     }


     //! Find the entry with maximum magnitude in the matrix
     //! Computes the maximum value, and the row in which this occurs
     //! Returns whether the rank is positive
     //!
     //! \param M Input matrix
     //! \param max_value Value with maximum magnitude in matrix (output)
     //! \param max_row Row containing element with maximum magnitude in matrix (output)
     //! \return True if max_value is above zero tolerance
     inline bool positive_rank(const Matrix_3D &M,
                               double &max_value,
                               Vector_3D &max_row) {

          // Rows
          max_value = -1.0;
          int max_row_index=-1;
          for (unsigned int i=0; i<3; ++i) {
               // Columns
               for (unsigned int j=0; j<3; ++j) {
                    double abs_value = std::fabs(M(i,j));
                    if (abs_value > max_value) {
                         max_value = abs_value;
                         max_row_index = i;
                    }
               }
          }

          max_row = M.row_vector(max_row_index);
          
          return (max_value > Math<double>::zero_tolerance);
     }
};

}

#endif
