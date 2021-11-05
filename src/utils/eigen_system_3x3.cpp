// eigen_system_3x3.cpp --- Efficient calculation of eigen vectors using analytical solution
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

#include "eigen_system_3x3.h"

namespace phaistos {

const double EigenSystem3x3::inv_3 = 1.0/3.0;
const double EigenSystem3x3::sqrt_3 = sqrt(3.0);

}
