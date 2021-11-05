// matrix.h --- matrix class. Includes interfaces to some BLAS and LAPACK routines
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


#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "utils/vector_nd.h"

// Platform dependent includes/definitions
// OS X
#ifdef __APPLE__
#include <Accelerate/Accelerate.h> 
extern "C" void dgemm_(const char*, const char*, const int*, const int*, const int*,
                       const double*, const double*, const int*, const double*,
                       const int*, const double*, double*, const int*);

// Linux
#else
#define __CLPK_integer int
#define __CLPK_doublereal double
extern "C" void dgemm_(const char*, const char*, const int*, const int*, const int*,
                       const double*, const double*, const int*, const double*,
                       const int*, const double*, double*, const int*);

extern "C" void dpotrf_(const char*, const int*, double*, const int*, int*);

extern "C" void dgetrf_(const int*, const int*, double*, const int*, int*, int*);
extern "C" void dgetri_(const int*, double*, const int*, const int*, double*, const int*, int*);
extern "C" void dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);

extern "C" void dposv_(const char*, const int*, const int*, double*,
                       const int*, double*, const int*, int*);
extern "C" void dtrtrs_(const char*, const char*, const char*, const int*, const int*, const double*,
                        const int*, double*, const int*, int*);
extern "C" void dsyev_(const char*, const char*, const int*, double*,
                       const int*, double*, double*, const int*, int*);
#endif


namespace phaistos {

class Matrix;


Matrix matrix_multiplication(double *A, double *B, int rows_A, int cols_A, int rows_B, int cols_B);


//! Matrix class.
//! The internal representation is M[i][j] = data[column_j*nRows+row_i] for compatibility with BLAS and LAPACK
class Matrix {
protected:

     //! Internal data
     double *data;

public:

     //! Number of rows
     __CLPK_integer n_row;

     //! Number of columns
     __CLPK_integer n_col;

     //! Construct matrix
     //!
     //! \param n_row Number of rows
     //! \param n_col Number of columns
     void construct(int n_row,int n_col) {
       this->n_row = n_row;
       this->n_col = n_col;
       if (data!=NULL)
            delete[] data;
       data = new double[n_col*n_row];
     }

     //! Constructor (default)
     Matrix(): data(NULL),n_row(0),n_col(0) {};

     
     //! Constructor
     //!
     //! \param n_row Number of rows
     //! \param n_col Number of columns
     //! \param val Value to initialize with
     Matrix(int n_row, int n_col, double val=UNINITIALIZED): data(NULL) {
          construct(n_row, n_col);

          if (is_initialized(val)) {
               for (int i=0; i<n_row*n_col; i++) {
                    this->data[i] = val;
               }
          }
     }

     //! Constructor - from 3D matrix
     //!
     //! \param matrix Input 3D matrix
     Matrix(const Matrix_3D &matrix): data(NULL) {
          construct(3, 3);
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = matrix[i][j];
               }
          }       
     }

     //! Constructor - from vector of vectors
     //!
     //! \param vec_matrix Vector of vector input matrix
     Matrix(std::vector<std::vector<double> > vec_matrix): data(NULL) {
          this->n_row = vec_matrix.size();
          this->n_col = 0;
          if (n_row>0)
               n_col = vec_matrix[0].size();
          data = new double[n_row*n_col];

          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = vec_matrix[i][j];
               }
          }
     }     

     //! Constructor - from 2D array
     //!
     //! \param n_row Number of rows
     //! \param n_col Number of columns
     //! \param array_matrix 2D input array
     Matrix(int n_row, int n_col, double **array_matrix): data(NULL) {
          construct(n_row, n_col);
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = array_matrix[i][j];
               }
          }
     }          
     
     //! Constructor - from 1D array
     //!
     //! \param n_row Number of rows
     //! \param n_col Number of columns
     //! \param array_matrix 1D input array
     Matrix(int n_row, int n_col, double *array_matrix): data(NULL) {
          this->n_row = n_row;
          this->n_col = n_col;
          data = new double[n_row*n_col];

          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = array_matrix[i*n_col+j];
               }
          }
     }          

     //! Constructor - vector of Vector_nD instances
     //!
     //! \param array_matrix std::vector of Vector_nD
     Matrix(std::vector<Vector_nD> *array_matrix): data(NULL) {
          this->n_row = array_matrix->size();
          this->n_col = (*array_matrix)[0].size;
          data = new double[n_row*n_col];

          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = (*array_matrix)[i][j];
               }
          }
     }


     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     Matrix(const Matrix &other): data(NULL) {
          construct(other.n_row, other.n_col);

          double *other_data = other.get_array();
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = other_data[j*n_row+i];
               }
          }
     }
     
     //! Destructor
     ~Matrix() {
          delete[] data;
     }


     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current matrix (*this)
     Matrix &operator=(const Matrix &other) {
          construct(other.n_row, other.n_col);
          
          int s=other.n_row*other.n_col;
          for (int i=0; i<s; i++)
               data[i] = other.data[i];
          return *this;
     };


     //! Overload in-place += operator (with matrix)
     //!
     //! \param other Source object to add to current matrix
     //! \return Current iterator (this)
     const Matrix& operator+=(const Matrix& other) {
          for (int i=0; i<n_row*n_col; i++) {
               this->data[i] += other.data[i];
          }
          return *this;
     }

     //! Overload in-place += operator (with value)
     //!
     //! \param value Numeric value to add to all entries
     //! \return Current iterator (this)
     const Matrix& operator+=(const double value) {
          for (int i=0; i<n_row*n_col; i++) {
               this->data[i] += value;
          }
          return *this;
     }

     //! Overload in-place -= operator (with matrix)
     //!
     //! \param other Source object to subtract from current matrix
     //! \return Current iterator (this)
     const Matrix& operator-=(const Matrix& other) {
          for (int i=0; i<n_row*n_col; i++) {
               this->data[i] -= other.data[i];
          }
          return *this;
     }

     //! Overload in-place -= operator (with value)
     //!
     //! \param value Numeric value to subtract from all entries
     //! \return Current iterator (this)
     const Matrix& operator-=(const double value) {
          for (int i=0; i<n_row*n_col; i++) {
               this->data[i] -= value;
          }
          return *this;
     }
     
     //! Overload () operator (row,column)
     //!
     //! \param row Row index
     //! \param column Column index
     //! \return vector element
     double& operator()(int row, int column) const {
          return data[column*n_row+row];
     }

     //! Initialize to zero
     void init_to_zero() {
          for (int i=0; i<n_row*n_col; i++) {
               this->data[i] = 0.0;
          }
     }

     //! Initialize to identity matrix
     void init_to_identity() {
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    if (i==j) {
                         this->data[i*n_col+j] = 1; 
                    } else {
                         this->data[i*n_col+j] = 0; 
                    }
               }
          }
     }

     //! Get internal array
     //!
     //! \return Internal 1D array
     double *get_array() const {
          return data;
     }

     //! Calculate transposed matrix
     //!
     //! \return Transposed matrix
     Matrix transpose() {
          Matrix new_matrix(n_col, n_row);
          double *new_data = new_matrix.get_array();
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    new_data[i*n_col+j] = data[j*n_row+i]; 
               }
          }
          return new_matrix;
     }

     
     //! Extract row from matrix
     //!
     //! \param i Row index
     //! \return Row vector
     Vector_nD get_row(int i) {
          Vector_nD row(n_col);
          int start_index = i;
          int step_size = n_row;
          for (int i=0; i<n_col; i++) {
               int index = start_index + i*step_size;
               row[i] = data[index];
          }
          return row;
     }

     //! Set a row
     //!
     //! \param i Row index
     //! \param v Input vector
     void set_row(int i, const Vector_nD &v) {
          int start_index = i;
          int step_size = n_row;
          for (int i=0; i<n_col; i++) {
               int index = start_index + i*step_size;
               data[index] = v[i];
          }       
     }
     
     //! Extract column from matrix
     //!
     //! \param i Column index
     //! \return Column vector
     Vector_nD get_column(int i) {
          Vector_nD column(n_row);
          int start_index = i*n_row;
          int step_size = 1;
          for (int i=0; i<n_row; i++) {
               int index = start_index + i*step_size;
               column[i] = data[index];
          }

          return column;
     }

     //! Set a column
     //!
     //! \param i Column index
     //! \param v Input vector
     void set_column(int i, const Vector_nD &v) {
          int start_index = i*n_row;
          int step_size = 1;
          for (int i=0; i<n_row; i++) {
               int index = start_index + i*step_size;
               data[index] = v[i];
          }       
     }

     //! Overload * operation (matrix multiplication)
     //!
     //! \param other Matrix to multiply with
     //! \return New matrix
     Matrix operator*(Matrix other) {
          double *A = data;
          double *B = other.get_array();
          return matrix_multiplication(A,B,
                                       n_row, n_col, other.n_row, other.n_col);
     }

     //! Overload * operation (matrix*vector)
     //!
     //! \param v Vector to multiply with
     //! \return New vector
     Vector_nD operator*(Vector_nD v) {
          double *A = data;
          double *B = v.get_array();
          return matrix_multiplication(A,B,
                                       n_row, n_col, v.size, 1).get_column(0);
     }
     //! Overload * operation (matrix*number)
     //!
     //! \param val Value to multiply with
     //! \return New matrix
     Matrix operator*(double val) {
          Matrix res(*this);
          for (int i=0; i<n_row*n_col; i++) {
               res.data[i] *= val;
          }
          return res;
     }
     
     //! Overload + operator
     //!
     //! \param other Matrix operant
     //! \return New matrix
     Matrix operator+(const Matrix& other) const {
          Matrix res(*this);
          for (int i=0; i<n_row*n_col; i++) {
               res.data[i] += other.data[i];
          }
          return res;
     }
     
     //! Overload - operator
     //!
     //! \param other Matrix operant
     //! \return New matrix
     Matrix operator-(const Matrix& other) const {
          Matrix res(*this);
          for (int i=0; i<n_row*n_col; i++) {
               res.data[i] -= other.data[i];
          }
          return res;
     }
     
     //! Cholesky decomposition
     //!
     //! \param type Choose between 'U' for upper and 'L' for lower triangular matrix
     //! \return New matrix
     Matrix cholesky(char type='L') {
          if (n_row != n_col)
               throw Matrix::MatrixError("Error (cholesky): Cannot compute cholesky decomposition for non-square matrix");

          Matrix new_matrix(*this);
          double *A = new_matrix.get_array();
          __CLPK_integer info;

          // Call LAPACK function
          dpotrf_(&type, &n_row, A, &n_row, &info);

          // Check error status
          if (info<0) {
               std::ostringstream s;           
               s << "Error (cholesky): Argument " << info << " has illegal value";
               throw Matrix::MatrixError(s.str());
          } else if (info>0) {
               std::ostringstream s;           
               s << "Error (cholesky): The leading minor of order " << info << " is not positive definite";
               throw Matrix::MatrixError(s.str());
          }

          // Set remaing entries to zero
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    if (((type=='L') && (i<j)) ||
                        ((type=='U') && (i>j)))
                         A[j*n_row+i] = 0.0;
               }
          }       
          
          return new_matrix;
     }


     //! LU decomposition of a general matrix
     //!
     //! \param pivots Optional output of number of pivots
     //! \return New matrix
     Matrix LU(__CLPK_integer *pivots=NULL) {

          if (n_row != n_col)
               throw Matrix::MatrixError("Error (LU): Cannot compute LU decomposition for non-square matrix");

          __CLPK_integer M = n_row;
          __CLPK_integer N = n_col;
          __CLPK_integer LDA = M;
          Matrix new_matrix(*this);
          double *A = new_matrix.get_array();
          __CLPK_integer info;
          
          int min_MN = M;
          if (N<M)
               min_MN = N;

          bool new_pivots = false;
          if (!pivots) {
               pivots = new __CLPK_integer[min_MN];
               new_pivots = true;
          }

          // Call LAPACK function
          dgetrf_(&M, &N, A, &LDA, pivots, &info);

          // Check error status
          if (info<0) {
               std::ostringstream s;
               s << "Error (LU): Argument " << info << " has illegal value";
               throw Matrix::MatrixError(s.str());
          } else if (info>0) {
               std::ostringstream s;
               //std::cout << *this << "\n";
               s << "Error (LU): U(" << info << "," << info << ") is exactly zero. Completed decomposition, but the division by zero will occur if used to solve a linear system.";
               throw Matrix::MatrixError(s.str());
          }

          if (new_pivots)
               delete[] pivots;
          
          return new_matrix;
     }
     

     //! Compute determinant using LU decomposition
     //!
     //! \return determinant
     double determinant() {

          __CLPK_integer *pivots = new __CLPK_integer[n_row];
          
          Matrix new_matrix(this->LU(pivots));
          
          double det=1;
          for (int i=0;i<n_row;i++) {
               det *= (pivots[i] != i+1)?-new_matrix(i,i):new_matrix(i,i);
          }
          delete[] pivots;
          return det;
     }

     
     //! Solves the linear system AX=B for positive definite square A.
     //!
     //! \param A A matrix (in AX=B system)
     //! \param B B matrix (in AX=B system)
     //! \param type Choose between 'U' for upper and 'L' for lower triangular matrix
     //! \return New matrix
     friend Matrix solve_ax_b_pos_def(const Matrix &A, const Matrix &B, char type) ;


     //! Solves the linear system AX=B where A is a triangular matrix. 
     //!
     //! \param A A matrix (in AX=B system)
     //! \param B B matrix (in AX=B system)
     //! \param type 'L': A is lower triangular, 'U': A is upper triangular 
     //! \param form 'N': no transpose, 'T': Transpose, 'C': Conjugate transpose
     //! \param diag 'U': diagonal contains only ones (unit triangular), 'N' non-unit triangular
     //! \return New matrix
     friend Matrix solve_ax_b_triangular(const Matrix &A, const Matrix &B, char type, char form, char diag);

     
     //! Solves the linear system AX=B for general square matrix A 
     //!
     //! \param A A matrix (in AX=B system)
     //! \param B B matrix (in AX=B system)
     //! \return New matrix
     friend Matrix solve_ax_b(const Matrix &A, const Matrix &B);


     //! Calculate the inverse of a general matrix using LU decomposition
     //!
     //! \return New matrix
     Matrix inverse() {

          __CLPK_integer N = n_col;
          __CLPK_integer LDA = N;
          __CLPK_integer *pivots = new __CLPK_integer[N];
          __CLPK_integer LWORK = N;
          __CLPK_doublereal *WORK = new __CLPK_doublereal[N];
          __CLPK_integer info=0;
          
          Matrix new_matrix = LU(pivots);
          
          dgetri_(&N, new_matrix.get_array(), &LDA, pivots, WORK, &LWORK, &info);

          delete[] pivots;
          delete[] WORK;

          return new_matrix;
     }

     //! Calculate eigen values and vectors of a symmetric matrix
     //!
     //! \param eigen_values Eigen values (output)
     //! \param eigen_vectors Eigen vectors (output)
     void eigen_sym(Vector_nD *eigen_values=NULL, Matrix *eigen_vectors=NULL) {

          char JOBZ = 'V';
          char UPLO = 'U';
          __CLPK_integer N = n_col;

          bool delete_eigenVectors = false;
          if (!eigen_vectors) {
               eigen_vectors = new Matrix(*this);
               delete_eigenVectors = true;
          } else {
               *eigen_vectors = *this;
          }
          __CLPK_doublereal *A = eigen_vectors->get_array();
          
          __CLPK_integer LDA = N;

          bool delete_eigenValues = false;
          if (!eigen_values) {
               eigen_values = new Vector_nD(N);
               delete_eigenValues = true;
          } 
          __CLPK_doublereal *W = eigen_values->get_array();
          
          __CLPK_integer INFO;
          __CLPK_integer LWORK = -1;

          __CLPK_doublereal WORK_tmp;
          dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, &WORK_tmp, &LWORK, &INFO);
          LWORK = (int)WORK_tmp;

          __CLPK_doublereal *WORK = new __CLPK_doublereal[LWORK];
          dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);

          if (delete_eigenVectors) {
               delete eigen_vectors;
          }
          if (delete_eigenValues) {
               delete eigen_values;
          }

          delete[] WORK;
     }

     
     //! Create a full symmetric matrix from a lower triangular one
     void symmetrize_from_L() {
          assert(n_row==n_col);
          for (int i=0;i<n_row;i++)
               for (int j=i+1;j<n_col;j++)
                    data[j*n_row+i] = data[i*n_row+j];
     }

     //! Create a full symmetric matrix from an upper triangular one
     void symmetrize_from_U() {
          assert(n_row==n_col);
          for (int i=0;i<n_row;i++)
               for (int j=0;j<i;j++)
                    data[j*n_row+i] = data[i*n_row+j];
     }
     
     //! Overload output operator
     friend std::ostream& operator<<(std::ostream& Out, const Matrix &m) {
          Out << "(";
          for (int i=0; i<m.n_row; i++) {
               if (i>0) {
                    Out << ",\n ";
               }
               Out << "(";
               for (int j=0; j<m.n_col; j++) {
                    if (j>0) {
                         Out << ", ";
                    }
                    Out << m.data[j*m.n_row+i];
               }
               Out << ")";

          }
          Out << ")";
          return Out;
     }

     friend Matrix operator-(const Matrix &m);
     friend Matrix log(const Matrix &m);
     friend Matrix exp(const Matrix &m);
     friend Matrix sq(const Matrix &m);
     friend double sum(const Matrix &m);
     friend double trace(const Matrix &m);
     friend Matrix reciproc(const Matrix &m);
     friend Matrix entrywise_multiplication(const Matrix &m1, const Matrix &m2);
     friend inline void save_matrix(std::string filename, const Matrix &A, bool append);
     
     //! Local exception class
     class MatrixError : public std::exception {
     private:
          //! Message
          std::string msg;
     public:
          //! Override what function
          virtual const char* what() const throw() {
               return(msg.c_str());
          }
          
          //! Constructor
          //!
          //! \param s Input message
          MatrixError(const std::string &s = "") : msg(s) {}
          
          //! Destructor
          ~MatrixError() throw() {}
     };
};

//! Overload negation operator
//!
//! \param m Matrix to operate on
//! \return New matrix
inline Matrix operator-(const Matrix &m) {
     Matrix res(m);
     for (int i=0; i<m.n_row*m.n_col; i++) {
          res.data[i] = -m.data[i];
     }
     return res;     
}

//! Logarithm of all elements in matrix
//!
//! \param m Matrix to operate on
//! \return New matrix
inline Matrix log(const Matrix &m) {
     Matrix res(m);
     for (int i=0; i<m.n_row*m.n_col; i++) {
          res.data[i] = std::log(m.data[i]);
     }
     return res;
}

//! exponent of all elements in matrix
//!
//! \param m Matrix to operate on
//! \return new matrix


inline Matrix exp(const Matrix &m) {
     Matrix res(m);
     for (int i=0;i<m.n_row*m.n_col;i++) {
          res.data[i] = std::exp(m.data[i]);
     }
     return res;
}


//! Sum of all elements in matrix
//!
//! \param m Matrix to operate on
//! \return Sum of all elements
inline double sum(const Matrix &m) {
     double res=0.0;
     for (int i=0; i<m.n_row*m.n_col; i++) {
          res+= m.data[i];
     }
     return res;
}


//! Trace of matrix m
//!
//! \param m Matrix to operate on
//! \return trace og matrix m
inline double trace(const Matrix &m) {
     double res=0.0;
     assert(m.n_row==m.n_col);
     for (int i=0; i<m.n_row; i++) {
          res+= m.data[i*m.n_row+i];
     }
     return res;
}


//! Square all elements in matrix
//!
//! \param m Matrix to operate on
//! \return New matrix
inline Matrix sq(const Matrix &m) {
     Matrix res(m);
     for (int i=0; i<m.n_row*m.n_col; i++) {
          res.data[i] = m.data[i]*m.data[i];
     }
     return res;
}

//! 1/x for all entries x in matrix
//!
//! \param m Matrix to operate on
//! \return New matrix
inline Matrix reciproc(const Matrix &m) {
    Matrix res(m);
    for (int i=0; i<m.n_row*m.n_col; i++) {
          res.data[i] = 1.0/m.data[i];
    }
    return res;
}

//! Multiply two matrices element by element
//!
//! \param m1 First matrix to operate on
//! \param m2 Second matrix to operate on
//! \return New matrix
inline Matrix entrywise_multiplication(const Matrix &m1, const Matrix &m2) {
    Matrix res(m1);
    for (int i=0; i<m1.n_row*m1.n_col; i++) {
          res.data[i] = m2.data[i]*m1.data[i];
    }
    return res;
}


//! Matrix multiplication (internal)
//!
//! \param A 1D array representation of Matrix
//! \param B 1D array representation of Matrix
//! \param rows_A Number of rows in A
//! \param cols_A Number of columns in A
//! \param rows_B Number of rows in B
//! \param cols_B Number of columns in B
//! \return New matrix
inline Matrix matrix_multiplication(double *A, double *B,
                                    int rows_A, int cols_A, int rows_B, int cols_B) {
     const char TRANSA = 'N';
     const char TRANSB = 'N';
     const int LDA = rows_A;
     const int LDB = rows_B;
          
     int M = rows_A;
     int N = cols_B;
     int K = cols_A;
     double ALPHA = 1.0;
     double BETA = 0.0;
          
     Matrix new_matrix(M,N);
     double *C = new_matrix.get_array();

     // Call BLAS function
     dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA,
            B, &LDB, &BETA, C, &M);

     return new_matrix;
}


//! Solves the linear system AX=B for positive definite square A.
//!
//! \param A A matrix (in AX=B system)
//! \param B B matrix (in AX=B system)
//! \param type Choose between 'U' for upper and 'L' for lower triangular matrix
//! \return New matrix
inline Matrix solve_ax_b_pos_def(const Matrix &A, const Matrix &B, char type='L') {
     if (A.n_row != A.n_col)
          throw Matrix::MatrixError("Error (solve_ax_b): Cannot solve for non-square matrix");
     if (A.n_row != B.n_row)
          throw Matrix::MatrixError("Error (solve_ax_b): Matrices are not conformable");
          
     Matrix new_matrix(B);
     Matrix A_copy(A);
          
     __CLPK_integer N = A.n_row;
     __CLPK_integer NRHS = B.n_col; 
     __CLPK_integer LDA = N;
     __CLPK_integer LDB = N;
     __CLPK_integer info=0;

     dposv_(&type, &N, &NRHS, A_copy.get_array(), &LDA, new_matrix.get_array(), &LDB, &info);

     if (info<0) {
          std::ostringstream s;        
          s << "Error (solve_ax_b_pos_def): Argument " << info << " has illegal value";
          throw Matrix::MatrixError(s.str());
     }

     return new_matrix;
}


//! Solves the linear system AX=B where A is a triangular matrix. 
//!
//! \param A A matrix (in AX=B system)
//! \param B B matrix (in AX=B system)
//! \param type 'L': A is lower triangular, 'U': A is upper triangular 
//! \param form 'N': no transpose, 'T': Transpose, 'C': Conjugate transpose
//! \param diag 'U': diagonal contains only ones (unit triangular), 'N' non-unit triangular
//! \return New matrix
inline Matrix solve_ax_b_triangular(const Matrix &A, const Matrix &B, char type='L', char form='N', char diag='N') {
     if (A.n_row != A.n_col)
          throw Matrix::MatrixError("Error (solve_ax_b): Cannot solve for non-square matrix");
     if (A.n_row != B.n_row)
          throw Matrix::MatrixError("Error (solve_ax_b): Matrices are not conformable");
          
     Matrix new_matrix(B);
     Matrix A_copy(A);
          
     __CLPK_integer N = A.n_row;
     __CLPK_integer NRHS = B.n_col; 
     __CLPK_integer LDA = N;
     __CLPK_integer LDB = N;
     __CLPK_integer info=0;

     dtrtrs_(&type, &form, &diag, &N, &NRHS, A_copy.get_array(), &LDA, new_matrix.get_array(), &LDB, &info);
          
     if (info<0) {
          std::ostringstream s;        
          s << "Error (solve_ax_b_triangular): Argument " << -info << " has illegal value";
          throw Matrix::MatrixError(s.str());
     } else if (info>0) {
          std::ostringstream s;        
          s << "Error (solve_ax_b_triangular): Diagonal element " << info << " of A is zero. A is singular and the solution for x cannot be computed.";
          throw Matrix::MatrixError(s.str());
     }

     return new_matrix;
}


//! Solves the linear system AX=B for general square matrix A 
//!
//! \param A A matrix (in AX=B system)
//! \param B B matrix (in AX=B system)
//! \return New matrix
inline Matrix solve_ax_b(const Matrix &A, const Matrix &B) {
     if (A.n_row != A.n_col)
          throw Matrix::MatrixError("Error (solve_ax_b): Cannot solve for non-square matrix");
     if (A.n_row != B.n_row)
          throw Matrix::MatrixError("Error (solve_ax_b): Matrices are not conformable");
          
     Matrix new_matrix(B);
     Matrix A_copy(A);
          
     __CLPK_integer N = A.n_row;
     __CLPK_integer NRHS = B.n_col; 
     __CLPK_integer LDA = N;
     __CLPK_integer LDB = N;
     __CLPK_integer info=0;

     __CLPK_integer maxMN = N;
          
     __CLPK_integer *pivots = new __CLPK_integer[maxMN];
          
     dgesv_(&N, &NRHS, A_copy.get_array(), &LDA, pivots, new_matrix.get_array(), &LDB, &info);

     if (info<0) {
          std::ostringstream s;        
          s << "Error (solve_ax_b): Argument " << info << " has illegal value";
          throw Matrix::MatrixError(s.str());
     } else if (info>0) {
          std::ostringstream s;        
          s << "Error (solve_ax_b): U(" << info << "," << info << ") is exactly zero. Completed decomposition, but the solution could not be computed.";
          throw Matrix::MatrixError(s.str());
     }

     delete[] pivots;
     return new_matrix;
}


//! Loads matrix (A) from file (filename) 
//!
//! \param filename filename of matrix
//! \return A matrix

inline Matrix load_matrix(std::string filename) {
     std::ifstream file_handle;
     std::vector < std::vector < double > > loaded_matrix;
     std::string read_line;
     std::vector < double > read_row;
     file_handle.open( filename.c_str() );
     if(!file_handle) {
          return Matrix(0,0);
     }
     do {
          getline(file_handle, read_line);
          std::istringstream row_elements(read_line);
          double index;
          do {
               row_elements >> index;
               read_row.push_back(index);
          } while(!row_elements.eof());
          if(read_row.size()>0)
               loaded_matrix.push_back(read_row);
          read_row.clear();
     } while(!file_handle.eof());

     return Matrix(loaded_matrix);
}

//! Saves matrix (A) to file (filename) 
//!
//! \param filename output filename of matrix
//! \param A matrix;
//! \param append Whether to append to file instead of overwriting it.
inline void save_matrix(std::string filename, const Matrix &A, bool append=true) {
     int i=0,j=0;
     std::ofstream file_stream;
     if(append) {
          file_stream.open(filename.c_str(), std::ios::app);
     } else {
          file_stream.open(filename.c_str());
     }
     assert(file_stream); 

     for(;i<A.n_row;i++) {
          for(j=0;j<A.n_col;j++) {
               file_stream << A.data[j*A.n_row+i]  << " ";
          }
          file_stream << std::endl;
     }
     file_stream.close();
}



//! Diagonal matrix class
class DiagonalMatrix: public Matrix {

     //! Constant zero
     const double zero;

public:
     
     //! Constructor
     //!
     //! \param size Size of square matrix
     //! \param val Value to initialize with
     DiagonalMatrix(int size, double val=UNINITIALIZED)
          : zero(0.0) {
          construct(size, 1);

          if (is_initialized(val)) {
               for (int i=0; i<size; i++) {
                    this->data[i] = val;
               }
          }

          // Articially override the n_col value to be the same as n_row
          this->n_col = size;
     }

     //! Overload () operator (row,column)
     //!
     //! \param row Row index
     //! \param column Column index
     //! \return vector element
     const double& operator()(int row, int column) const {
          if (row != column)
               return zero;
          else 
               return data[row];
     }

     //! Overload () operator (row,column)
     //!
     //! \param row Row index
     //! \param column Column index
     //! \return vector element
     double& operator()(int row, int column) {
          assert(row == column);
          return data[row];
     }

     //! Calculated the inverse of a general matrix
     //!
     //! \return New DiagonalMatrix
     DiagonalMatrix inverse() {
          DiagonalMatrix new_matrix(n_row);

          for (int i=0; i<n_row; i++) {
               new_matrix.data[i] = 1.0/this->data[i];;
          }          
          return new_matrix;
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current matrix (*this)
     DiagonalMatrix &operator=(const DiagonalMatrix &other) {
          construct(other.n_row, 1);
          
          for (int i=0; i<n_row; i++)
               data[i] = other.data[i];
          return *this;
     };

};

}

#endif
