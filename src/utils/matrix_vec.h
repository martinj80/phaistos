// matrixvec.h --- matrix of vector_3D objects
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


#ifndef MATRIXVEC_H
#define MATRIXVEC_H

#include <vector>
#include <iostream>
#include <sstream>
#include <utils/vector_matrix_3d.h>
#include <utils/matrix.h>

namespace phaistos {

//! Matrix of vector_3D objects
class MatrixVec {

     //! Internal data (array of Vector_3D objects)
     Vector_3D *data;

public:

     //! Number of rows
     int n_row;

     //! Number of columns     
     int n_col;

     //! Construct matrix
     //!
     //! \param n_row Number of rows
     //! \param n_col Number of columns
     void construct(int n_row,int n_col) {
          this->n_row = n_row;
          this->n_col = n_col;
          if (data!=NULL)
               delete[] data;
          data = new Vector_3D[n_col*n_row];
          for (int i=0; i<n_row*n_col; i++) {
               data[i] = Vector_3D(0.0, 0.0, 0.0);
          }
     }
     
     //! Constructor (default)
     MatrixVec(): data(NULL),n_row(0),n_col(0) {};

     //! Constructor
     //!
     //! \param n_row Number of rows
     //! \param n_col Number of columns
     MatrixVec(int n_row, int n_col): data(NULL) {
          construct(n_row, n_col);
     }

     //! Destructor
     ~MatrixVec() {
          delete[] data;
     }

     //! Get internal array
     //!
     //! \return Internal 1D array
     Vector_3D *get_array() const {
          return data;
     }

     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     MatrixVec(const MatrixVec &other): data(NULL) {
          construct(other.n_row, other.n_col);

          Vector_3D *other_data = other.get_array();
          for (int i=0; i<n_row; i++) {
               for (int j=0; j<n_col; j++) {
                    data[j*n_row+i] = other_data[j*n_row+i];
               }
          }
     }
     
     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current matrix (*this)
     MatrixVec &operator=(const MatrixVec &other) {
          construct(other.n_row, other.n_col);
          
          int s=other.n_row*other.n_col;
          for (int i=0; i<s; i++)
               data[i] = other.data[i];
          return *this;
     };

     //! Overload () operator (row,column)
     //!
     //! \param row Row index
     //! \param column Column index
     //! \return vector element
     Vector_3D& operator()(int row, int column) {
          return data[column*n_row+row];
     }

     //! Overload * operation (matrix multiplication)
     //!
     //! \param other Matrix to multiply with
     //! \return New matrix
     Matrix operator*(MatrixVec &other) {
          Matrix res(n_row, other.n_col);
          assert(n_col == other.n_row);
          for (int i=0; i<n_row; i++) {
               for (int j2=0; j2<other.n_col; j2++) {
                    double sum = 0.0;
                    for (int j=0; j<n_col; j++) {
                         sum += ((*this)(i,j))*other(j,j2);
                    }
                    res(i,j2) = sum;
               }
          }
          return res;
     }

     //! Overload output operator
     friend std::ostream& operator<<(std::ostream& o, const MatrixVec &m) {
          o << "(";
          for (int i=0; i<m.n_row; i++) {
               if (i>0) {
                    o << ",\n ";
               }
               o << "(";
               for (int j=0; j<m.n_col; j++) {
                    if (j>0) {
                         o << ", ";
                    }
                    o << m.data[j*m.n_row+i];
               }
               o << ")";

          }
          o << ")";
          return o;
     }     
};

}

#endif
