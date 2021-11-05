// profasi_energy_term_common.h --- Common code for PROFASI energy terms
// Copyright (C) 2010 Pengfei Tian, Wouter Boomsma
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

// Implementation of the PROFASI force field
// An effective all-atom potential for proteins, Irback, Mitternacht, Mohanty, PMC Biophysics 2009

#ifndef PROFASI_ENERGY_TERM_COMMON_H
#define PROFASI_ENERGY_TERM_COMMON_H

namespace phaistos {

//! Base class for Profasi energy terms
class EnergyTermCommonProfasi {
protected:
     // In "An effective all-atom potential for proteins, Irback,
     // Mitternacht, Mohanty, PMC Biophysics 2009" a simulation with
     // the profasi force field produces a heat capacity maximum at
     // RT=0.4722 eu (profasi units). Setting this temperature equal
     // to the melting temperature (315K), we have:
     // 1eu = 315/0.4722 * 1.9858/1000 kCal/mol
     const static double profasi_energy_in_kcal_per_mol;
};

//! conversion factor between eu and kcal/mol
const double EnergyTermCommonProfasi::profasi_energy_in_kcal_per_mol = (315.0/0.4722)*1.9858/1000.0;

}

#endif
