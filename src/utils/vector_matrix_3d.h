// vectorMatrix_3D.h --- 3D vector and matrix classes
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


#ifndef VECTORMATRIX_3D_H
#define VECTORMATRIX_3D_H

#define SGN(x) (((x)<0) ? -1 : 1)

#include "utils.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace phaistos {

class Matrix_3D;
class Vector_3D;

// Forward declaration
Matrix_3D& svd_rotation_matrix(Matrix_3D& Sigma, Matrix_3D& Gamma, Vector_3D& singular_values);          

//! 3D Vector class
class Vector_3D {
private:

     //! Internal data
     double data[3];

public:
     //! Whether vector has been initialized
     bool initialized;

     //! Constructor (default)
     Vector_3D() {
          initialized=false;
     }

     //! Constructor - initialize based on x,y,z components.
     Vector_3D(const double x, const double y, const double z) {
          data[0] = x;
          data[1] = y;
          data[2] = z;
          initialized=true;
     }

     //! Constructor - initialize based on array of components
     Vector_3D(const double array[3]) {
          data[0] = array[0];
          data[1] = array[1];
          data[2] = array[2];
          initialized=true;          
     }

     //! Constructor - initialize based on vector of components
     Vector_3D(const std::vector<double> &v) {
          data[0] = v[0];
          data[1] = v[1];
          data[2] = v[2];
          initialized=true;          
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     Vector_3D(const Vector_3D &other) {
          initialized = other.initialized;
          data[0] = other.data[0];
          data[1] = other.data[1];
          data[2] = other.data[2];
     }

     //! Initialize all elements to zero
     void init_to_zero() {
          data[0] = 0.0;
          data[1] = 0.0;
          data[2] = 0.0;
          initialized=true;
     }     

     //! Return internal data array
     double *get_array() {
          return this->data;
     }
     
     //! Get size (always 3)
     int size() {
          return 3;
     }
     
     //! Overload [] indexing operator (const)
     double operator[](const int index) const {
          //assert(index<3);
          return data[index];
     }

     //! Overload indexing operator (non-const)
     double& operator[](const int index) {
          //assert(index<3);
          return data[index];
     }

     //! Overload assignment operator
     const Vector_3D& operator=(const Vector_3D &v) {
          data[0] = v.data[0];
          data[1] = v.data[1];
          data[2] = v.data[2];
          this->initialized = v.initialized;
          return *this;
     }

     //! Overload assignment operator (array assignment)
     const Vector_3D& operator=(const double array[3]) {
          data[0] = array[0];
          data[1] = array[1];
          data[2] = array[2];
          this->initialized = true;
          return *this;
     }

     //! Overload += operator
     const Vector_3D& operator+=(const Vector_3D& v2) {
          data[0] += v2.data[0];
          data[1] += v2.data[1];
          data[2] += v2.data[2];
          return *this;
     }

     //! Overload += operator (with scalar value)
     const Vector_3D& operator+=(const double value) {
          data[0] += value;
          data[1] += value;
          data[2] += value;
          return *this;
     }

     //! Overload -= operator
     const Vector_3D& operator-=(const Vector_3D& v2) {
          data[0] -= v2.data[0];
          data[1] -= v2.data[1];
          data[2] -= v2.data[2];
          return *this;
     }

     //! Overload -= operator (with scalar value)
     const Vector_3D& operator-=(const double value) {
          data[0] -= value;
          data[1] -= value;
          data[2] -= value;
          return *this;
     }

     //! Overload *= operator (with scalar value)
     const Vector_3D& operator*=(const double value) {
          data[0] *= value;
          data[1] *= value;
          data[2] *= value;
          return *this;
     }

     //! Overload /= operator (with scalar value
     const Vector_3D& operator/=(const double value) {
          data[0] /= value;
          data[1] /= value;
          data[2] /= value;
          return *this;
     }

     //! Overload == operator
     bool operator==(const Vector_3D& v2) {
          return (
               (data[0] == v2.data[0]) &&
               (data[1] == v2.data[1]) &&
               (data[2] == v2.data[2])
            );
     }

     //! Overload != operator
     bool operator!=(const Vector_3D& v2) {
          return !( *this == v2 );
     }

     //! Check if vector is a null-vector
     bool is_null() const {
          return (data[0]==0 &&
                  data[1]==0 &&
                  data[2]==0);
     }

     //! Calculate L2 norm (length of vector)
     double norm() const {
          return sqrt(data[0]*data[0] +
                      data[1]*data[1] +
                      data[2]*data[2]);
     }

     //! Calculate squared L2 norm (squared length of vector)
     double norm_squared() const {
          return (data[0]*data[0] +
                  data[1]*data[1] +
                  data[2]*data[2]);
     }

     //! Normalize vector
     Vector_3D normalize() const {
          double norm = this->norm(); 
          if (norm != 0) {
               return *this/norm;
          } else {
               return *this;
          }
     }

     //! Overload output operator
     friend std::ostream& operator<<(std::ostream& Out, const Vector_3D &v) {
          Out << std::setiosflags( std::ios::fixed );
          Out << std::setprecision(4);
          if (v.initialized)
               Out << "[" << v.data[0] << ", " << v.data[1]<< ", " << v.data[2] << "]";
          else
               Out << "[" << UNINITIALIZED << ", " << UNINITIALIZED << ", " << UNINITIALIZED << "]";
          return Out;
     }

     //! Overload + operator (Vector_3D + Vector_3D)
     friend Vector_3D operator+(const Vector_3D &v1, const Vector_3D& v2);

     //! Overload + operator (Vector_3D + scalar)
     friend Vector_3D operator+(const Vector_3D &v1, const double value);

     //! Overload - operator (Vector_3D - Vector_3D)
     friend Vector_3D operator-(const Vector_3D &v1, const Vector_3D& v2);

     //! Overload - operator (Vector_3D - scalar)     
     friend Vector_3D operator-(const Vector_3D &v1, const double value);

     //! Overload - operator (negate)          
     friend Vector_3D operator-(const Vector_3D &v1);

     //! Overload * operator (dot product)
     friend double operator*(const Vector_3D &v1, const Vector_3D& v2);

     //! Overload * operator (Vector_3D * scalar)     
     friend Vector_3D operator*(const Vector_3D &v1, const double value);

     //! Overload * operator (scalar * Vector_3D)     
     friend Vector_3D operator*(const double value, const Vector_3D &v1);

     //! Overload / operator (Vector_3D / scalar)
     friend Vector_3D operator/(const Vector_3D &v1, const double value);

     //! Overload ^ operator (Outer product)
     friend Matrix_3D operator^(const Vector_3D &v1,const Vector_3D &v2);

     //! Overload % operator (Cross product)
     friend Vector_3D operator%(const Vector_3D &v1,const Vector_3D &v2);

     //! Cross product
     friend Vector_3D cross_product(const Vector_3D& v1, const Vector_3D& v2);

     //! Absolute value of elements
     friend Vector_3D fabs(const Vector_3D& v);

     //! Sum of elements
     friend double sum(const Vector_3D& v);
     
     //! Calculate angle between two normalized vectors
     friend double vector_angle_norm(const Vector_3D & v1,const Vector_3D &v2);

     //! Calculate angle between two vectors
     friend double vector_angle(const Vector_3D& v1, const Vector_3D &v2);

     //! Calculate angle defined by 3 position vectors
     friend double calc_angle(const Vector_3D &v1, 
                              const Vector_3D &v2, 
                              const Vector_3D &v3);

     //! Calculate dihedral angle from two binormals (must be normalized)
     friend double dihedral_from_binormal(const Vector_3D &b1,
                                          const Vector_3D &b2,
                                          const Vector_3D &connection);

     //! Calculate dihedral angle from 3 vectors
     friend double calc_dihedral(const Vector_3D &v12,const Vector_3D &v23,const Vector_3D &v34);

     //! Calculate dihedral from 4 position vectors
     friend double calc_dihedral(const Vector_3D &v1, const Vector_3D &v2, 
                                 const Vector_3D &v3, const Vector_3D &v4);

     //! Convert spherical coordinate vector(r,theta,phi) to cartesian coordinate set. 
     //! Angular range: 0<theta<2pi, 0<phi<pi
     friend Vector_3D spherical_to_cartesian(Vector_3D& v);

     //! Convert cartesian coordinate set to spherical coordinate set vector(r,theta,phi)
     friend Vector_3D cartesian_to_spherical(Vector_3D& v);

     //! Return a vector containing zeros
     friend Vector_3D null_vector();

     //! Calculate position D based on three existing positions A, B, C
     //! and the angles ACD and DCB
     //! direction==1 gives solution in direction of normal CAxCB.
     //! direction==-1 gives solution in opposite direction.
     friend Vector_3D calc_pos_from_3_pos_and_2_angles(Vector_3D &A,
                                                       Vector_3D &C,
                                                       Vector_3D &B,
                                                       double angle1,
                                                       double angle2,
                                                       double length_CD,
                                                       bool direction=1);

     //! Calculate position D based on three existing positions A, B, C
     //! and the angles ACD and DCB - - geometric solution - slightly slower
     friend Vector_3D calc_pos_from_3_pos_and_2_angles_alt(Vector_3D &B,
                                                           Vector_3D &A,
                                                           Vector_3D &C,
                                                           double alpha,
                                                           double beta,
                                                           double l,
                                                           bool choice=0);

     //! Local Iterator typedef
     typedef std::iterator<std::random_access_iterator_tag, Vector_3D > Iterator;
};


//! Overload + operator (Vector_3D + Vector_3D)
inline Vector_3D operator+(const Vector_3D &v1, const Vector_3D& v2) {
     return Vector_3D(v1.data[0] + v2.data[0],
                      v1.data[1] + v2.data[1],
                      v1.data[2] + v2.data[2]);
}
     
//! Overload + operator (Vector_3D + scalar)
inline Vector_3D operator+(const Vector_3D &v1, const double value) {
     return Vector_3D(v1.data[0] + value,
                      v1.data[1] + value,
                      v1.data[2] + value);
}

//! Overload - operator (Vector_3D - Vector_3D)
inline Vector_3D operator-(const Vector_3D &v1, const Vector_3D& v2) {
     return Vector_3D(v1.data[0] - v2.data[0],
                      v1.data[1] - v2.data[1],
                      v1.data[2] - v2.data[2]);
}

//! Overload - operator (Vector_3D - scalar)     
inline Vector_3D operator-(const Vector_3D &v1, const double value) {
     return Vector_3D(v1.data[0] - value,
                      v1.data[1] - value,
                      v1.data[2] - value);
}

//! Overload - operator (negate)          
inline Vector_3D operator-(const Vector_3D &v1) {
     return Vector_3D(-v1.data[0],
                      -v1.data[1],
                      -v1.data[2]);
}

//! Overload * operator (dot product)
inline double operator*(const Vector_3D &v1, const Vector_3D& v2) {
     return v1.data[0] * v2.data[0] + v1.data[1] * v2.data[1] + v1.data[2] * v2.data[2];
}

//! Overload * operator (Vector_3D * scalar)     
inline Vector_3D operator*(const Vector_3D &v1, const double value) {
     return Vector_3D(v1.data[0] * value,
                      v1.data[1] * value,
                      v1.data[2] * value);
}

//! Overload * operator (scalar * Vector_3D)     
inline Vector_3D operator*(const double value, const Vector_3D &v1) {
     return Vector_3D(v1.data[0] * value,
                      v1.data[1] * value,
                      v1.data[2] * value);
}

//! Overload / operator (Vector_3D / scalar)
inline Vector_3D operator/(const Vector_3D &v1, const double value) {
     return Vector_3D(v1.data[0] / value,
                      v1.data[1] / value,
                      v1.data[2] / value);
}

//! Cross product
inline Vector_3D cross_product(const Vector_3D& v1, const Vector_3D& v2) {
     return Vector_3D((v1.data[1]*v2.data[2])-(v1.data[2]*v2.data[1]),
                      (v1.data[2]*v2.data[0])-(v1.data[0]*v2.data[2]),
                      (v1.data[0]*v2.data[1])-(v1.data[1]*v2.data[0]));
}

//! Overload % operator (Cross product)
inline Vector_3D operator%(const Vector_3D &v1,const Vector_3D &v2) {
     return cross_product(v1,v2);
}


//! Absolute value of elements
inline Vector_3D fabs(const Vector_3D &v) {
     return Vector_3D(std::fabs(v.data[0]),
                      std::fabs(v.data[1]),
                      std::fabs(v.data[2]));
}

//! Sum of elements
inline double sum(const Vector_3D &v) {
     return v.data[0] + v.data[1] + v.data[2];
}

//! Calculate angle between two normalized vectors
inline double vector_angle_norm(const Vector_3D & v1,const Vector_3D &v2) {
     double product = v1* v2;
     if (product>1.0) 
          product=1.0;
     else if (product<-1.0) 
          product=-1.0;
     return acos(product);
}

//! Calculate angle between two vectors
inline double vector_angle(const Vector_3D& v1, const Vector_3D &v2) {
     if (v1.is_null() || v2.is_null()) return 0;
     return vector_angle_norm(v1.normalize(),v2.normalize());
}


//! Calculate angle defined by 3 position vectors
inline double calc_angle(const Vector_3D &v1, const Vector_3D &v2, 
                        const Vector_3D &v3) {
     return vector_angle(v1 - v2, v3 - v2);
}
     

//! Calculate dihedral angle from two binormals (must be normalized)
inline double dihedral_from_binormal(const Vector_3D &b1,
                                     const Vector_3D &b2, 
                                     const Vector_3D &connection) {
     double angle = vector_angle_norm(b1, b2);
     // Determine sign of angle
     Vector_3D orientation=b1%b2;
     if (orientation*connection<0) angle=-angle;
     return angle;
}


//! Calculate dihedral angle from 3 vectors
inline double calc_dihedral(const Vector_3D &v12,const Vector_3D &v23,const Vector_3D &v34) {
     Vector_3D cross1 = cross_product(v12, v23); //binormal at v2
     Vector_3D cross2 = cross_product(v23, v34); //binormal at v3
     double angle=vector_angle(cross1,cross2);
     if (cross1*v34<0) angle=-angle;
     return angle;      
}

//! Calculate dihedral from 4 position vectors
inline double calc_dihedral(const Vector_3D &v1, const Vector_3D &v2, 
                           const Vector_3D &v3, const Vector_3D &v4) {
     return calc_dihedral(v2-v1,v3-v2,v4-v3); 
}
     
//! Convert spherical coordinate vector(r,theta,phi) to cartesian coordinate set. 
//! Angular range: 0<theta<2pi, 0<phi<pi
inline Vector_3D spherical_to_cartesian(Vector_3D& v) {
     double r = v[0];
     double theta = v[1];
     double phi = v[2];
     double x = r * cos(theta) * sin(phi);
     double y = r * sin(theta) * sin(phi);
     double z = r * cos(phi);
     return Vector_3D(x,y,z);
}

//! Convert cartesian coordinate set to spherical coordinate set vector(r,theta,phi)
inline Vector_3D cartesian_to_spherical(Vector_3D& v) {
     double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
     double theta = atan2(v[1], v[0]);
     double phi = acos(v[2]/r);
     return Vector_3D(r,theta,phi);
}

//! Return a vector containing zeros
inline Vector_3D null_vector() {
     return Vector_3D(0,0,0);
}

//! Calculate position D based on three existing positions A, B, C
//! and the angles ACD and DCB
//! direction==1 gives solution in direction of normal CAxCB.
//! direction==-1 gives solution in opposite direction.
inline Vector_3D calc_pos_from_3_pos_and_2_angles(Vector_3D &A,
                                                  Vector_3D &C,
                                                  Vector_3D &B,
                                                  double angle1,
                                                  double angle2,
                                                  double length_CD,
                                                  bool direction) {

     Vector_3D CA = A-C;
     Vector_3D CB = B-C;
     double CAlength = CA.norm();
     double CBlength = CB.norm();

     // Save results for speed
     double cosAngle1 = cos(angle1);
     double cosAngle2 = cos(angle2);
     double factor = (CA[0]*CB[1] - CA[1]*CB[0]);

     // Coefficients
     double a = ((CA[2]*CB[1] - CA[1]*CB[2]) * (CA[2]*CB[1] - CA[1]*CB[2]) +
                 (CA[0]*CB[2] - CA[2]*CB[0]) * (CA[0]*CB[2] - CA[2]*CB[0]) +
                 factor*factor);

     double b = 2*length_CD* ((CA[1]*CA[2]*CB[1] - CA[1]*CA[1]*CB[2] -
                              CA[0]*CA[0]*CB[2] + CA[0]*CA[2]*CB[0])*CBlength*cosAngle2 -
                             (CA[2]*CB[1]*CB[1] - CA[1]*CB[1]*CB[2] -
                              CA[0]*CB[0]*CB[2] + CA[2]*CB[0]*CB[0])*CAlength*cosAngle1);

     double c = length_CD*length_CD*((CA[1]*CBlength*cosAngle2 - CB[1]*CAlength*cosAngle1)*
                                   (CA[1]*CBlength*cosAngle2 - CB[1]*CAlength*cosAngle1)+
                                   (CA[0]*CBlength*cosAngle2 - CB[0]*CAlength*cosAngle1)*
                                   (CA[0]*CBlength*cosAngle2 - CB[0]*CAlength*cosAngle1) -
                                   factor*factor);
     // Discriminant
     double d = b*b-4*a*c;

     Vector_3D res;
     if (d>1e-11) {     // This method seems to be numerically unstable for small d
          double CDX1 = -((CA[2]*CB[1] - CA[1]*CB[2])/factor);
          double CDX2 = length_CD*((CA[1]*CBlength*cosAngle2 - CB[1]*CAlength*cosAngle1)/
                                  factor);
          double CDY1 = -((CA[0]*CB[2] - CA[2]*CB[0])/factor);
          double CDY2 = length_CD*((CA[0]*CBlength*cosAngle2 - CB[0]*CAlength*cosAngle1)/
                                  factor);

          double sqrtD = sqrt(d);
          double q = -0.5*(b+SGN(b)*sqrtD);
          res[2] = c/q;
          res[0] = CDX1*res[2] - CDX2;
          res[1] = CDY1*res[2] + CDY2;

          // Use normal vector to determine which solution to pick
          Vector_3D normal = (A-C)%(B-C);
          if ((res*normal)*direction < 0) {
               res[2] = q/a;
               res[0] = CDX1*res[2] - CDX2;
               res[1] = CDY1*res[2] + CDY2;
          }

          // Translate CD vector to D vector
          res = res+C;
               
     } else {   // Try alternative method
          #if DEBUGLEVEL > 1
          std::cerr << "Warning: calcPosFrom3PosAnd2Angles - Falling back to geometric method (discriminant=" << d << ")\n";
          #endif
          res = calc_pos_from_3_pos_and_2_angles_alt(A,
                                                     C,
                                                     B,
                                                     angle1,
                                                     angle2,
                                                     length_CD);
     }
          
     return res;
}

//! Calculate position D based on three existing positions A, B, C
//! and the angles ACD and DCB - - geometric solution - slightly slower
inline Vector_3D calc_pos_from_3_pos_and_2_angles_alt(Vector_3D &B,
                                                      Vector_3D &A,
                                                      Vector_3D &C,
                                                      double alpha,
                                                      double beta,
                                                      double l,
                                                      bool choice) {

     double gamma = calc_angle(B, A, C);
     double Dx = l*cos(alpha);
     double Dy = l*(cos(beta)-cos(alpha)*cos(gamma))/sin(gamma);
     double Dz = sqrt(l*l-Dx*Dx-Dy*Dy);

     if (std::isnan(Dz))
          return Vector_3D();
     
     double Dab   = Dx-cos(gamma)/sin(gamma)*Dy;
     double Dac   = Dy/sin(gamma);
     double Dabac = Dz;

     Vector_3D AB    = (B-A).normalize();
     Vector_3D AC    = (C-A).normalize();
     Vector_3D ABxAC = (AB%AC).normalize();

     Vector_3D AD;
     if (choice == 0) {
          AD = AB*Dab + AC*Dac + ABxAC*Dabac;
     } else {
          AD = AB*Dab + AC*Dac - ABxAC*Dabac;
     }

     return A+AD;     
}



//! Matrix class
class Matrix_3D {
private:

     //! Row vectors     
     Vector_3D data[3];

public:

     //! Constructor (default)
     Matrix_3D() {}

     //! Constructor - Initialize based on three vectors (optionally inserted as columns)
     Matrix_3D(const Vector_3D& v1, const Vector_3D& v2, const Vector_3D& v3, const bool as_columns=false) {
          if (!as_columns) {
               data[0] = v1;
               data[1] = v2;
               data[2] = v3;
          } else {
               data[0] = Vector_3D(v1[0], v2[0], v3[0]);
               data[1] = Vector_3D(v1[1], v2[1], v3[1]);
               data[2] = Vector_3D(v1[2], v2[2], v3[2]);
          }
     }
     
     //! Constructor - Initialize based on 3x3 array
     Matrix_3D(const double matrix_array[][3]) {
          data[0] = Vector_3D(matrix_array[0]);
          data[1] = Vector_3D(matrix_array[1]);
          data[2] = Vector_3D(matrix_array[2]);       
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     Matrix_3D(const Matrix_3D &other) {
          data[0] = Vector_3D(other.data[0]);
          data[1] = Vector_3D(other.data[1]);
          data[2] = Vector_3D(other.data[2]);
     }

     //! Initialize to zero
     void init_to_zero() {
          data[0].init_to_zero();
          data[1].init_to_zero();
          data[2].init_to_zero();
     }
          
     //! Overload [] indexing operator (const)
     const Vector_3D operator[](const int index) const {
          assert(index<3);
          return data[index];
     }

     //! Overload [] indexing operator (non const)
     Vector_3D& operator[](const int index) {
          assert(index<3);
          return data[index];
     }

     //! Overload () indexing operator using two indices (const)
     double operator()(const int index1, const int index2) const {
          // assert(index1<3);
          return data[index1][index2];
     }

     //! Overload () indexing operator using two indices (non const)
     double& operator()(const int index1, const int index2) {
          // assert(index1<3);
          return data[index1][index2];
     }

     //! Overload + operator (Matrix_3D + Matrix_3D)
     Matrix_3D operator+(const Matrix_3D& m2) const{
          return Matrix_3D(Vector_3D(data[0]+m2.row_vector(0)),Vector_3D(data[1]+m2.row_vector(1)),Vector_3D(data[2]+m2.row_vector(2)));
          
     }
   
     //! Overload += operator
     const Matrix_3D& operator+=(const Matrix_3D& m2) {
          data[0] += m2.data[0];
          data[1] += m2.data[1];
          data[2] += m2.data[2];
          return *this;
     }
     
     //! Overload - operator (Matrix_3D - Matrix_3D)
     Matrix_3D operator-(const Matrix_3D& m2) const{
          return Matrix_3D(Vector_3D(data[0]-m2.row_vector(0)),
                           Vector_3D(data[1]-m2.row_vector(1)),
                           Vector_3D(data[2]-m2.row_vector(2)));
          
     }
   
     //! Overload -= operator
     const Matrix_3D& operator-=(const Matrix_3D& m2) {
          data[0] -= m2.data[0];
          data[1] -= m2.data[1];
          data[2] -= m2.data[2];
          return *this;
     }
     
     //! Overload * operator (Matrix_3D * scalar value)
     Matrix_3D operator*(const double val) const {
          return Matrix_3D(Vector_3D(val*data[0]), Vector_3D(val*data[1]), Vector_3D(val*data[2]));
     }   
     
     //! Overload *= operator (Matrix_3D * scalar value)
     Matrix_3D operator*=(const double val) {
          data[0] *= val;
          data[1] *= val;
          data[2] *= val;
          return *this;
     }   
     
     //! Overload * operator (matrix multiplication)
     Matrix_3D operator*(const Matrix_3D& m2) const {
          return Matrix_3D(Vector_3D(data[0]*m2.col_vector(0),
                                     data[0]*m2.col_vector(1),
                                     data[0]*m2.col_vector(2)),
                           Vector_3D(data[1]*m2.col_vector(0),
                                     data[1]*m2.col_vector(1),
                                     data[1]*m2.col_vector(2)),
                           Vector_3D(data[2]*m2.col_vector(0),
                                     data[2]*m2.col_vector(1),
                                     data[2]*m2.col_vector(2)));
     }
   
     

     //! Overload * operator (with scalar value)
     Vector_3D operator*(const Vector_3D& v) const {
          return Vector_3D(data[0]*v,
                           data[1]*v,
                           data[2]*v);
     }
     
     
     //! Return copy of row as vector
     //!
     //! \param index Row index
     Vector_3D row_vector(const int index) const {
          return Vector_3D(data[index][0], data[index][1], data[index][2]);
     }

     //! Return copy of column as vector
     //!
     //! \param index Row index
     Vector_3D col_vector(const int index) const {
          return Vector_3D(data[0][index], data[1][index], data[2][index]);
     }

     //! Calculate determinant
     double determinant() const {
          return 
               data[0][0]*(data[1][1]*data[2][2] - data[2][1]*data[1][2]) -
               data[0][1]*(data[1][0]*data[2][2] - data[2][0]*data[1][2]) +
               data[0][2]*(data[1][0]*data[2][1] - data[2][0]*data[1][1]) ;
     }
     
     //! Fill 3x3 array with values from matrix
     void as_array(double matrix[][3]) {
          for (int i=0; i<3; i++) {
               for (int j=0; j<3; j++) {
                    matrix[i][j] = data[i][j];
               }
          }
     }

     //! Overload output operator
     friend std::ostream& operator<<(std::ostream& o, const Matrix_3D &m) {
          o << "[" << m[0] << ",\n " << m[1]<< ",\n " << m[2] << "]\n";
          return o;
     }

     //! Return transpose matrix
     friend Matrix_3D transpose(const Matrix_3D &m);

     //! Return identity matrix
     friend Matrix_3D identity_matrix();

     //! Return matrix containing zeros
     friend Matrix_3D null_matrix();

     //! Construct orthonormal base vectors and insert them as rows in matrix
     friend Matrix_3D orthonormal_basis_matrix(const Vector_3D& v1, const Vector_3D& v2);

     //! Return rotation matrix for rotation around x-axis
     friend Matrix_3D rotation_matrix_x(const double angle);

     //! Return rotation matrix for rotation around y-axis
     friend Matrix_3D rotation_matrix_y(const double angle);

     //! Return rotation matrix for rotation around z-axis
     friend Matrix_3D rotation_matrix_z(const double angle);
     
     //! Return rotation matrix for rotation around unit vector
     //! NOTE: r must be normalized!     
     friend Matrix_3D rotation_matrix_from_dihedral(const Vector_3D &r,double dihedral);

     //! Absolute value of all entries
     friend Matrix_3D fabs(const Matrix_3D& v);

     //! Sum of all entries
     friend double sum(const Matrix_3D& v);
     
     
     //! Matrix multiple a 3xN matrix with a Nx3 matrix (both represented as arrays of 3D-vectors)
     template <typename ITERATORTYPE>     
     friend Matrix_3D vector_list_matrix_multiplication(ITERATORTYPE vector_iterator1_start, ITERATORTYPE vector_iterator1_end,
                                                        Vector_3D &translation1,
                                                        ITERATORTYPE vector_iterator2_start, ITERATORTYPE vector_iterator2_end,
                                                        Vector_3D &translation2);

     //! Calculate superimposition rotation matrix for two vectors (iterators) of Vector_3D of equal length
     template <typename ITERATORTYPE>     
     friend void calc_super_imposition_rotation_matrix(ITERATORTYPE vector_iterator1_start, ITERATORTYPE vector_iterator1_end,
                                                       ITERATORTYPE vector_iterator2_start, ITERATORTYPE vector_iterator2_end,
                                                       Vector_3D& cm1, Vector_3D& cm2, Matrix_3D& rotation_matrix, Vector_3D& singular_values);
};


//! Return transposed matrix
inline Matrix_3D transpose(const Matrix_3D &m) {
     return Matrix_3D(Vector_3D(m.data[0][0], m.data[1][0], m.data[2][0]),
                      Vector_3D(m.data[0][1], m.data[1][1], m.data[2][1]),
                      Vector_3D(m.data[0][2], m.data[1][2], m.data[2][2]));        
}

//! Return identity matrix
inline Matrix_3D identity_matrix() {
     return Matrix_3D(Vector_3D(1,0,0),
                      Vector_3D(0,1,0),
                      Vector_3D(0,0,1));
}

//! Return matrix containing zeros
inline Matrix_3D null_matrix() {
     return Matrix_3D(Vector_3D(0,0,0),
                      Vector_3D(0,0,0),
                      Vector_3D(0,0,0));
}

//! Construct orthonormal base vectors and insert them as rows in matrix
inline Matrix_3D orthonormal_basis_matrix(const Vector_3D& v1, const Vector_3D& v2) {
     Vector_3D X = v1.normalize();
     Vector_3D Z = cross_product(v2,X).normalize();
     Vector_3D Y = cross_product(X,Z);

     // Vector_3D X = v1.normalize();
     // Vector_3D Y = v2.normalize();
     // Y = (Y-(X*(X*Y))).normalize();
     
     // Vector_3D Z = cross_product(X,Y);

     bool as_columns = true;
     return Matrix_3D(X,Y,Z,as_columns);
} 

//! Return rotation matrix for rotation around x-axis
inline Matrix_3D rotation_matrix_x(const double angle) {
     return Matrix_3D(Vector_3D(1,          0,          0),
                      Vector_3D(0, cos(angle), -sin(angle)),
                      Vector_3D(0,sin(angle), cos(angle)));          
}

//! Return rotation matrix for rotation around y-axis
inline Matrix_3D rotation_matrix_y(const double angle) {
     return Matrix_3D(Vector_3D( cos(angle), 0, sin(angle)),
                      Vector_3D(          0, 1,           0),
                      Vector_3D(-sin(angle), 0, cos(angle)));
}

//! Return rotation matrix for rotation around z-axis
inline Matrix_3D rotation_matrix_z(const double angle) {
     return Matrix_3D(Vector_3D( cos(angle),-sin(angle), 0),
                      Vector_3D( sin(angle), cos(angle), 0),
                      Vector_3D(          0,          0, 1));
}

//! Absolute value of all entries
inline Matrix_3D fabs(const Matrix_3D &m) {
     return Matrix_3D(fabs(m[0]),
                      fabs(m[1]),
                      fabs(m[2]));
}

//! Sum of all entries
inline double sum(const Matrix_3D &m) {
     return sum(m[0]) + sum(m[1]) + sum(m[2]);
}

//! Return rotation matrix for rotation around unit vector
//! NOTE: r must be normalized!     
inline Matrix_3D rotation_matrix_from_dihedral(const Vector_3D &r,double dihedral) {
     double cosa=cos(dihedral),sina=sin(dihedral);
     double dcosa=1.0-cosa;
     double sinx=sina*r[0],siny=sina*r[1],sinz=sina*r[2];
     return Matrix_3D(Vector_3D(r[0]*r[0]*dcosa+cosa,r[0]*r[1]*dcosa-sinz,r[0]*r[2]*dcosa+siny),
                      Vector_3D(r[1]*r[0]*dcosa+sinz,r[1]*r[1]*dcosa+cosa,r[1]*r[2]*dcosa-sinx),
                      Vector_3D(r[2]*r[0]*dcosa-siny,r[2]*r[1]*dcosa+sinx,r[2]*r[2]*dcosa+cosa));

}


//! Matrix multiple a 3xN matrix with a Nx3 matrix (both represented as arrays (iterators) of 3D-vectors).
//! Uses new style iterators (containing knowledge of their own endpoint)
template <typename ITERATORTYPE>     
inline Matrix_3D vector_list_matrix_multiplication(ITERATORTYPE it1,
                                                   Vector_3D &translation1,
                                                   ITERATORTYPE it2,
                                                   Vector_3D &translation2) {

     Matrix_3D res = Matrix_3D(Vector_3D(0,0,0),
                               Vector_3D(0,0,0),
                               Vector_3D(0,0,0));

     for (int k=0; k<3; k++) {
          for (int j=0; j<3; j++) {
               ITERATORTYPE it1_tmp(it1);
               ITERATORTYPE it2_tmp(it2);
               for (; !it1_tmp.end() && !it2_tmp.end(); ++it1_tmp, ++it2_tmp) {
                    res[j][k] += ((*it1_tmp)[j] - translation1[j]) * ((*it2_tmp)[k] - translation2[k]);
               }
          }
     }
     return res;
}


//! Matrix multiple a 3xN matrix with a Nx3 matrix (both represented as arrays (iterators) of 3D-vectors)
template <typename ITERATORTYPE>     
inline Matrix_3D vector_list_matrix_multiplication(ITERATORTYPE vector_iterator1_start, ITERATORTYPE vector_iterator1_end,
                                                   Vector_3D &translation1,
                                                   ITERATORTYPE vector_iterator2_start, ITERATORTYPE vector_iterator2_end,
                                                   Vector_3D &translation2) {

     Matrix_3D res = Matrix_3D(Vector_3D(0,0,0),
                               Vector_3D(0,0,0),
                               Vector_3D(0,0,0));

     for (int k=0; k<3; k++) {
          for (int j=0; j<3; j++) {
               for (ITERATORTYPE it1 = vector_iterator1_start, it2 = vector_iterator2_start;
                    it1 != vector_iterator1_end && it2 != vector_iterator2_end;
                    ++it1, ++it2) {
                    res[j][k] += ((*it1)[j] - translation1[j]) * ((*it2)[k] - translation2[k]);
               }
          }
     }
     return res;
}


//! Calculate center of mass within range
template <typename ITERATOR_TYPE>     
Vector_3D center_of_mass(ITERATOR_TYPE it) {
     Vector_3D res = Vector_3D(0,0,0);
     int counter = 0;
     for (; !it.end(); ++it) {
          res += *it;
          counter++;
     }
     res /= counter;
     return res;
}

//! Calculate center of mass within range
template <typename ITERATOR_TYPE>     
Vector_3D center_of_mass(ITERATOR_TYPE begin, ITERATOR_TYPE end) {
     Vector_3D res = Vector_3D(0,0,0);
     int counter = 0;
     for (ITERATOR_TYPE it=begin; it != end; ++it) {
          res += *it;
          counter++;
     }
     res /= counter;
     return res;
}


//! Calculate coveriance matrix over an iterator of positions
template <typename ITERATOR_TYPE>     
Matrix_3D covariance_matrix(ITERATOR_TYPE it) {

     Matrix_3D cov;
     cov.init_to_zero();

     Vector_3D cm = center_of_mass(it);

     for (; !it.end(); ++it) {
          Vector_3D vec = (*it-cm);
          cov += vec^vec;
     }
     return cov;
}


//! Calculate coveriance matrix over an iterator of positions
template <typename ITERATOR_TYPE>     
Matrix_3D covariance_matrix(ITERATOR_TYPE begin, ITERATOR_TYPE end) {

     Matrix_3D cov;
     cov.init_to_zero();

     Vector_3D cm = center_of_mass(begin, end);

     for (ITERATOR_TYPE it=begin; it != end; ++it) {
          Vector_3D vec = (*it-cm);
          cov += vec^vec;
     }
     return cov;
}


//! Calculate superimposition rotation matrix for two vectors (iterators) of Vector_3D of equal length
template <typename ITERATOR_TYPE>     
inline void calc_superimpose_rotation_matrix(ITERATOR_TYPE it1, ITERATOR_TYPE it2,
                                             Vector_3D& cm1, Vector_3D& cm2, Matrix_3D& rotation_matrix, Vector_3D& singular_values) {

     // Calculate center of mass
     cm1 = center_of_mass(it1);
     cm2 = center_of_mass(it2);

     // Calculate correlation matrix
     Matrix_3D cm_matrix = vector_list_matrix_multiplication(it2, cm2,
                                                            it1, cm1);

     // SVD
     svd_rotation_matrix(cm_matrix, rotation_matrix, singular_values);
}     


//! Calculate superimposition rotation matrix for two vectors (iterators) of Vector_3D of equal length
template <typename ITERATOR_TYPE>     
inline void calc_superimpose_rotation_matrix(ITERATOR_TYPE vector_iterator1_start, ITERATOR_TYPE vector_iterator1_end,
                                             ITERATOR_TYPE vector_iterator2_start, ITERATOR_TYPE vector_iterator2_end,
                                             Vector_3D& cm1, Vector_3D& cm2, Matrix_3D& rotation_matrix, Vector_3D& singular_values) {

     // Calculate center of mass
     cm1 = center_of_mass(vector_iterator1_start, vector_iterator1_end);
     cm2 = center_of_mass(vector_iterator2_start, vector_iterator2_end);

     // Calculate correlation matrix
     Matrix_3D cm_matrix = vector_list_matrix_multiplication(vector_iterator2_start, vector_iterator2_end, cm2,
                                                            vector_iterator1_start, vector_iterator1_end, cm1);

     // SVD
     svd_rotation_matrix(cm_matrix, rotation_matrix, singular_values);
}     

//! Overload ^ operator (Outer product)
inline Matrix_3D operator^(const Vector_3D &v1,const Vector_3D &v2){
     Matrix_3D M;
     for (int i = 0 ; i<3; i++){
          for (int j = 0; j < 3; j++){
               M[i][j] = v1[i] * v2[j];
          }
     }
     return M;
}

}

#endif
