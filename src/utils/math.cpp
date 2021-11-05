// math.cpp --- numerical constants
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

#include <math.h>
#include "math.h"

namespace phaistos {

//@{ 
//! Float implementation
template<> const float Math<float>::zero_tolerance = 1e-06f;
template<> const float Math<float>::sqrt_2 = (float)(sqrt(2.0)); 
template<> const float Math<float>::sqrt_2_pi = (float)sqrt(2.0*M_PI);
template<> const float Math<float>::inv_sqrt_2 = 1.0f/Math<float>::sqrt_2; 
//@}

//@{
//! Double implementation
template<> const double Math<double>::zero_tolerance = 1e-08;
template<> const double Math<double>::sqrt_2 = sqrt(2.0); 
template<> const double Math<double>::sqrt_2_pi = (double)sqrt(2.0*M_PI);
template<> const double Math<double>::inv_sqrt_2 = 1.0/Math<double>::sqrt_2;
//@}

}
