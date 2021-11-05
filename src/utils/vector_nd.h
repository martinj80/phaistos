// vector_nd.h --- n-dimensional vector class
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

#ifndef VECTOR_ND_H
#define VECTOR_ND_H

#include "utils.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdarg>

namespace phaistos {

//! n-dimensional vector class
class Vector_nD {
private:

     //! Internal data
     double *data;
public:

     //! Number of elements in vector
     int size;

     //! Constructor
     Vector_nD(): data(NULL){
          size = -1;
     }

     //! Constructor: Only allocate memory
     Vector_nD(int size): data(NULL) {
          init(size);
     }

     //! Constructor: Initialize based on variable argument list.
     //! THIS METHOD HAS SOME ISSUES WHEN USED WITH NON DIGIT FLOAT CONSTANTS. e.g. Vector_nD v(2,10,10) - always call using Vector_nD v(2, 10.0, 10.0)
     //! First element is explicit to avoid ambiguity with previous constructor
     Vector_nD(int size, double first_element, ...): data(NULL) {
          init(size);
          
          va_list ap;
          this->data[0] = first_element;
          va_start(ap,first_element);
          for (int i=1; i<size; i++) {
               this->data[i] = va_arg(ap, double);
          }
          va_end(ap);
     }

     
     //! Constructor: From array
     Vector_nD(int size, double *data): data(NULL) {
          init(size);
          for (int i=0; i<this->size; i++) {
               this->data[i] = data[i];
          }
     }

     //! Constructor: From vector
     Vector_nD(std::vector<double> data): data(NULL) {
          init(data.size());
          for (int i=0; i<this->size; i++) {
               this->data[i] = data[i];
          }
     }

     //! Constructor: From pointer vector
     Vector_nD(std::vector<double *> data): data(NULL) {
          init(data.size());
          for (int i=0; i<this->size; i++) {
               this->data[i] = *data[i];
          }
     }

     //! Copy constructor
     Vector_nD(const Vector_nD &other): data(NULL) {
          init(other.size);
          for (int i=0; i<this->size; i++) {
               this->data[i] = other.data[i];
          }
     }

     //! Initialize
     void init(int size) {
          this->size = size;
	  if (data!=NULL)
	       delete[] data;
          this->data = new double[this->size];
     }


     //! Clear contents
     void clear() {
          this->size = -1;
          if (this->data)
               delete[] this->data;
	  this->data = NULL;
     }

     //! Destructor
     ~Vector_nD() {
          if (data != NULL) {
               delete[] data;
          }
	  this->data = NULL;
     }

     //! Check whether vector is initialized
     bool initialized() {
          return (size!=-1);
     }
     
     //! Overload [] indexing operator (const)
     double operator[](const int index) const {
          return data[index];
     }

     //! Overload [] indexing operator (non-const)
     double& operator[](const int index) {
          return data[index];
     }

     //! Overload + operator (Vector_nD + Vector_nD)
     Vector_nD operator+(const Vector_nD& v2) const {
          Vector_nD res(std::min(this->size, v2.size));
          for (int i=0; i<res.size; i++)  {
               res.data[i] = data[i] + v2.data[i];
          }
          return res;
     }

     //! Overload + operator (Vector_nD + scalar value)
     Vector_nD operator+(const double value) const {
          Vector_nD res(this->size);
          for (int i=0; i<res.size; i++)  {
               res.data[i] = data[i] + value;
          }
          return res;
     }

     //! Overload - operator (Vector_nD - Vector_nD)
     Vector_nD operator-(const Vector_nD& v2) const {
          Vector_nD res(std::min(this->size, v2.size));
          for (int i=0; i<res.size; i++)  {
               res.data[i] = this->data[i] - v2.data[i];
          }
          return res;
     }

     //! Overload - operator (Vector_nD - scalar value)
     Vector_nD operator-(const double value) const {
          Vector_nD res(this->size);
          for (int i=0; i<res.size; i++)  {
               res.data[i] = this->data[i] - value;
          }
          return res;
     }

     //! Overload - operator (unary negate)
     Vector_nD operator-() const {
          Vector_nD res(size);
          for (int i=0; i<res.size; i++)  {
               res.data[i] = -data[i];
          }
          return res;
     }

     //! Overload * operator (dot product)
     double operator*(const Vector_nD& v2) const {
          double res = 0.0;
          assert(this->size == v2.size);
          for (int i=0; i<size; i++)  {
               res += data[i] * v2.data[i];
          }
          return res;
     }

     //! Overload * operator (Vector_nD * scalar value))
     Vector_nD operator*(const double value) const {
          Vector_nD res(this->size);
          for (int i=0; i<res.size; i++)  {
               res[i] = data[i] * value;
          }
          return res;
     }

     //! Overload / operator (Vector_nD / scalar value))
     Vector_nD operator/(const double value) const {
          Vector_nD res(this->size);
          for (int i=0; i<res.size; i++)  {
               res[i] = data[i] / value;
          }
          return res;
     }

     //! Overload / operator (entrywise division)
     Vector_nD operator/(const Vector_nD& v2) const {
          assert(this->size == v2.size);
          Vector_nD res(this->size);
          for (int i=0; i<res.size; i++)  {
               res[i] = data[i] / v2.data[i];
          }
          return res;
     }


     //! Overload assignment operator
     const Vector_nD& operator=(const Vector_nD& v2) {
          if (this->size != v2.size) {
               if (this->size != -1) {
                    if (data) {
                         delete[] data;
                    }
		    data = NULL;
               }
               init(v2.size);
          }
          for (int i=0; i<size; i++)  {
               data[i] = v2.data[i];
          }
          return *this;
     }

     //! Overload assignment operator (array)
     const Vector_nD& operator=(const double *array) {
          assert(size != -1);
          for (int i=0; i<size; i++)  {
               data[i] = array[i];
          }
          return *this;
     }

     //! Overload assignment operator (scalar value)
     const Vector_nD& operator=(const double value) {
          assert(size != -1);
          for (int i=0; i<size; i++)  {
               data[i] = value;
          }
          return *this;
     }
     
     //! Overload += operator
     const Vector_nD& operator+=(const Vector_nD& v2) {
          for (int i=0; i<v2.size; i++)  {
               data[i] += v2.data[i];
          }
          return *this;
     }

     //! Overload += operator (with scalar value)
     const Vector_nD& operator+=(const double value) {
          for (int i=0; i<size; i++)  {
               data[i] += value;
          }
          return *this;
     }

     //! Overload -= operator
     const Vector_nD& operator-=(const Vector_nD& v2) {
          for (int i=0; i<v2.size; i++)  {
               data[i] -= v2.data[i];
          }
          return *this;
     }

     //! Overload -= operator (with scalar value)
     const Vector_nD& operator-=(const double value) {
          for (int i=0; i<size; i++)  {
               data[i] -= value;
          }
          return *this;
     }

     //! Overload *= operator (with scalar value)
     const Vector_nD& operator*=(const double value) {
          for (int i=0; i<size; i++)  {
               data[i] *= value;
          }
          return *this;
     }

     //! Overload /= operator (with scalar value)
     const Vector_nD& operator/=(const double value) {
          for (int i=0; i<size; i++)  {
               data[i] /= value;
          }
          return *this;
     }

     //! Overload > operator
     bool operator>(const Vector_nD& v2) {
          return (this->norm() > v2.norm());
     }

     //! Overload >= operator
     bool operator>=(const Vector_nD& v2) {
          return (this->norm() >= v2.norm());
     }

     //! Overload < operator
     bool operator<(const Vector_nD& v2) {
          return (this->norm() < v2.norm());
     }

     //! Overload <= operator
     bool operator<=(const Vector_nD& v2) {
          return (this->norm() <= v2.norm());
     }

     //! Calculate length of vector (L2 norm)
     double norm() const {
          double res = 0.0;
          for (int i=0; i<size; i++)  {
               res += data[i]*data[i];
          }
          return sqrt(res);
     }

     //! Calculate squared length of vector (L2 norm squared)
     double norm_sq() const {
          double res = 0.0;
          for (int i=0; i<size; i++)  {
               res += data[i]*data[i];
          }
          return res;
     }
 

     //! Normalize vector
     Vector_nD normalize() const {
          double norm = this->norm(); 
          if (norm != 0) {
               return *this/norm;
          } else {
               return *this;
          }
     }

     //! Return internal data array
     double *get_array() {
          return this->data;
     }

     //! Return Vector_nD as std::vector
     std::vector<double> get_vector() {
          return (std::vector<double>(this->data, this->data+this->size));
     }

     //! Create vector initialized with range of numbers
     friend Vector_nD range(double first_element, double last_element, double step_size);

     //! entrywise multiplication
     friend Vector_nD entrywise_multiplication(const Vector_nD &v1, const Vector_nD &v2);

     //! Return maximum element in vector
     friend double max(const Vector_nD &v);

     //! Return minimum element in vector
     friend double min(const Vector_nD &v);

     //! Return sum of elements
     friend double sum(const Vector_nD &v);

     //! Return average value of elements
     friend double mean(const Vector_nD &v);
     
     //! Return log Vector_nD of input Vector_nD
     friend Vector_nD log(const Vector_nD &v);

     //! Return new vector containing cumulative sum of the elements so far
     friend Vector_nD cum_sum(const Vector_nD &v);

     //! Save vector to file
     friend inline void save_vector(std::string filename, const Vector_nD &A, bool append);

     //! Apply specified function to all elements in vector
     friend Vector_nD fmap(double(*f)(double), const Vector_nD &v);
     
     //! Overload output operator
     friend std::ostream& operator<<(std::ostream& o, const Vector_nD &v) {
          o << "(";
          for (int i=0; i<v.size; i++) {
               if (i>0) {
                    o << ", ";
               }
               o << v.data[i];
          }
          o << ")";
          return o;
     }
};


//! Create vector initialized with range of numbers
inline Vector_nD range(double first_element, double last_element, double step_size) {
     int size = (int)std::fabs(ceil((last_element-first_element)/step_size));
     Vector_nD res(size);
     int direction;
     if (first_element < last_element) {
	  direction = 1;
     } else {
	  direction = -1;
     }
     for (int i=0; i<size; i++) {
	  res.data[i] = first_element + step_size * i * direction;
     }
     return res;
}
     
inline Vector_nD entrywise_multiplication(const Vector_nD &v1, const Vector_nD &v2) {
     Vector_nD result(v1);
     for(int i=0;i<v1.size;i++) {
          result[i]=v1[i]*v2[i];
     }
     return result;
}

//! Return maximum element in vector
inline double max(const Vector_nD &v) {
     double max_value = v.data[0];
     for (int i=1; i<v.size; i++) {
	  if (v.data[i] > max_value) {
	       max_value = v.data[i];
	  }
     }
     return max_value;
}

//! Return minimum element in vector
inline double min(const Vector_nD &v) {
     double min_value = v.data[0];
     for (int i=1; i<v.size; i++) {
	  if (v.data[i] < min_value) {
	       min_value = v.data[i];
	  }
     }
     return min_value;
}
     
//! Return sum of elements
inline double sum(const Vector_nD &v) {
     double sum = 0.0;
     for (int i=0; i<v.size; i++) {
	  sum += v.data[i];
     }
     return sum;
}

//! Return log vector_nd of input vector_nd
inline Vector_nD log(const Vector_nD &v) {
     Vector_nD res(v);
     for (int i=0; i<v.size; i++) {
          res.data[i] = std::log(v.data[i]);
     }
     return res;
}

//! Return average value of elements
inline double mean(const Vector_nD &v) {
     double sum = 0.0;
     for (int i=0; i<v.size; i++) {
	  sum += v.data[i];
     }
     return sum/v.size;
}

//! Saves Vector_nD (A) to file (filename) 
//!
//! \param filename output filename of matrix
//! \param A Vector_nD;
//! \param append Whether to append to file instead of overwriting it.
inline void save_vector(std::string filename, const Vector_nD &A, bool append=true) {
     int j=0;
     std::ofstream file_stream;
     if(append) {
          file_stream.open(filename.c_str(),std::ios::app);
     } else {
          file_stream.open(filename.c_str());
     }
     assert(file_stream);

     for(;j<A.size;j++) {
          file_stream << A.data[j]  << " ";
     }
     file_stream << std::endl;
     
     file_stream.close();
}


//! Return new vector containing cumulative sum of the elements so far
inline Vector_nD cum_sum(const Vector_nD &v) {
     Vector_nD res(v.size);
     double sum = 0.0;
     for (int i=0; i<v.size; i++) {
	  sum += v.data[i];
	  res.data[i] = sum;
     }
     return res;
}


//! Apply specified function to all elements in vector
inline Vector_nD fmap(double(*f)(double), const Vector_nD &v) {
     Vector_nD res(v.size);
     for (int i=0; i<v.size; i++) {
	  res.data[i] = (*f)(v.data[i]);
     }
     return res;
}
     
}

#endif
