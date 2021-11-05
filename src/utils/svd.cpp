// svd.cpp --- Wrapper for the LAPACK dgesvd fortran function
// Copyright (C) 2005-2008 Wouter Boomsma
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


#include "svd.h"

namespace phaistos {

/* Wrapper for the LAPACK dgesvd fortran function */
int dgesvd(double sigma[][3], double u[][3], double s[3], double vt[][3]) {
     int i,j;
     int size = 3;

     char JOBU = 'A';
     char JOBVT = 'A';
     __CLPK_integer M = size;
     __CLPK_integer N = size;
     __CLPK_integer LDA = size;
     __CLPK_integer LDU = size;
     __CLPK_integer LDVT = size;

     /* Unfolded matrices for fortran call */
     double *A = new double[size*size];
     double *S = new double[size*size];
     double *U = new double[size*size];
     double *VT = new double[size*size];
     double *WORK = new double[(N+M)*32];
     __CLPK_integer LWORK = (N+M)*32;
     __CLPK_integer INFO;

     for (i=0; i<size; i++) {
          for (j=0; j<size; j++) {
               A[i*size+j] = sigma[j][i];
          }
     }

     // Fortran function call
     dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);

     // Transform data into C-arrays
     if (INFO == 0) {
          for (i=0; i<size; i++) {
               s[i] = S[i];
               for (j=0; j<size; j++) {
                    sigma[j][i] = A[i*size+j];
                    u[j][i] = U[i*size+j];
                    vt[j][i] = VT[i*size+j];
               }
          }
     }

     delete[] A;
     delete[] S;
     delete[] U;
     delete[] VT;
     delete[] WORK;
     
     return INFO;
}


/* Rotation matrix based on singular value decomposition */
Matrix_3D& svd_rotation_matrix(Matrix_3D& Sigma, Matrix_3D& Gamma, Vector_3D& singular_values) {

     double sigma[3][3];
     double u[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
     double vt[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
     double s[3] = {0,0,0};

     Sigma.as_array(sigma);

     if (dgesvd(sigma, u, s, vt)==1) {
          fprintf(stderr, "SVD call failed\n");
          exit(1);
     } 

     Matrix_3D U = Matrix_3D(u);
     Matrix_3D VT = Matrix_3D(vt);

     /* Test for reflection (determinant == -1) */
     if ((U.determinant()*VT.determinant()) + 1.0 < 0.0001) {
          Matrix_3D S = Matrix_3D(Vector_3D(1,0,0),
                                  Vector_3D(0,1,0),
                                  Vector_3D(0,0,-1));
          VT = S*VT;
          Gamma = U*VT;

          // Also change sign of smallest singular value
          s[2] *= -1;
     } else {
          Gamma = U*VT;
     }

     singular_values = s;
     
     return Gamma;
}

}
