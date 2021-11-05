// svd.h --- Wrapper for the LAPACK dgesvd fortran function
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

#ifndef SVD_H
#define SVD_H

#include "utils/vector_matrix_3d.h"

// Platform dependent includes/definitions
// OS X
#ifdef __APPLE__
#include <Accelerate/Accelerate.h> 

// Linux
#else
#define __CLPK_integer int
extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n,
               double *a, int *lda, double *s, double *u,
               int *ldu, double *vt, int *ldvt, double *work,
               int *lwork, int *info);
#endif

namespace phaistos {

int dgesvd(double sigma[][3], double u[][3], double vt[][3]);

//! Single value decomposition (Sigma = U*S*VT)
//!
//! \param Sigma Input matrix
//! \param Gamma Output matrix (U*VT)
//! \param singular_values diagonal of S
//! \return Gamma result matrix
Matrix_3D& svd_rotation_matrix(Matrix_3D& Sigma, Matrix_3D& Gamma, Vector_3D& singular_values);

}

#endif
