// math.h --- numerical constants
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

#ifndef MATH_H
#define MATH_H

namespace phaistos {

//! Numerical constants
template <typename REAL_TYPE>
class Math {
public:

     //! Cutoff when values are considered zero
     static const REAL_TYPE zero_tolerance;

     //! Square root of 2
     static const REAL_TYPE sqrt_2;

     //! Inverse square root of 2
     static const REAL_TYPE inv_sqrt_2;

     //! Inverse square root of 2*PI
     static const REAL_TYPE sqrt_2_pi;

     //! Square
     static REAL_TYPE sqr(REAL_TYPE value) {
          return value*value;
     }
};

}

#endif
