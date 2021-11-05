// xtc_chain.h --- Input/Output chains in XTC trajectory format
// Copyright (C) 2008-2012 Wouter Boomsma
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

#ifndef XTC_CHAIN_H
#define XTC_CHAIN_H

extern "C" {
#include "xdrfile.h"
#include "xdrfile_xtc.h"
}

#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "protein/chain.h"

namespace phaistos {

//! Write chain in XTC format.
//! \param chain Molecule chain
//! \param xdr_file XDR File object 
//! \param step Iteration index
//! \param time time into simulation
//! \param prec precision
template <typename CHAIN_TYPE>
void xtc_write_chain(const CHAIN_TYPE &chain, XDRFILE *xdr_file, int step, float time, float prec=1000) {

     std::vector<Vector_3D> coordinate_vector;
     for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(chain); !it.end(); ++it) {
          coordinate_vector.push_back(it->position);
     }

     unsigned int n_atoms = coordinate_vector.size();
     rvec *coordinates = new rvec[n_atoms];

     for (unsigned int i=0; i<n_atoms; ++i) {
          for(unsigned j=0; j<DIM; j++) {
               // convert to nanometer
               coordinates[i][j] = 0.1*coordinate_vector[i][j];
          }               
     }

     matrix box = {{0,0,0},
                   {0,0,0},
                   {0,0,0}};

     write_xtc(xdr_file, 
               n_atoms,
               step,
               time,
               box,
               coordinates,
               prec);

     delete[] coordinates;
}

//! Transfer coordinates to chain
//! \param chain Molecule chain
//! \param atom coordinates
template <typename CHAIN_TYPE>
void xtc_coordinates_to_chain(CHAIN_TYPE &chain, rvec *coordinates) {
     int atom_index = 0;
     for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(chain); !it.end(); ++it,++atom_index) {
          it->position = Vector_3D(coordinates[atom_index][0],
                                   coordinates[atom_index][1],
                                   coordinates[atom_index][2])*10; // convert to Aangstroem
     }
     chain.update_angles(); 
}

//! Read chain in XTC format.
//! \param chain Molecule chain
//! \param xdr_file XDR File object 
//! \param step Iteration index
//! \param time time into simulation
//! \param prec precision
//! \param n_atoms pointer to variable containing the number of atoms (to avoid recalculation with repeated calls to this function)
//! \param coordinates Optional atom coordinates
//! \param dry_run Do not do assignment to chain
template <typename CHAIN_TYPE>
bool xtc_read_chain(CHAIN_TYPE &chain, XDRFILE *xdr_file, int &step, float &time, float &prec, int *n_atoms=NULL, boost::optional<rvec **> coordinates=boost::optional<rvec**>(), bool dry_run=false) {

     int n_atoms_internal = 0;
     if (n_atoms==NULL) {
          n_atoms = &n_atoms_internal;
     }

     if ((*n_atoms) == 0) {
          for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(chain); !it.end(); ++it) {
               (*n_atoms)++;
          }
     }


     bool coordinates_owner = false;
     if (!coordinates) {
          coordinates_owner = true;
     }
     if (!coordinates || **coordinates == NULL) {
          **coordinates = new rvec[*n_atoms];
     }

     matrix box;
     int success = !read_xtc(xdr_file, *n_atoms, &step, &time, box, **coordinates, &prec);

     if (success && !dry_run) {
          xtc_coordinates_to_chain(chain, **coordinates);
     }
     
     if (coordinates_owner)
          delete[] (**coordinates);

     return success;
}


//! Check whether trajectory is consistent with chain
//! \param chain Molecule chain
//! \param xdr_filename Filename of XDR trajectory file  
template <typename CHAIN_TYPE>
void xtc_chain_check_consistency(CHAIN_TYPE &chain, std::string &xdr_filename) {
     int n_atoms_chain = 0;
     for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(chain); !it.end(); ++it) {
          n_atoms_chain++;
     }

     int n_atoms_xdr_file;
     if ((read_xtc_natoms((char *)xdr_filename.c_str(), &n_atoms_xdr_file)) != exdrOK) {
          std::cerr << "Error in determining number of atoms in XTC file\n";
          return;
     }
     
     if (n_atoms_chain != n_atoms_xdr_file) {
          std::cerr << "Warning: Number of atoms in chain and trajectory do not match\n";
     }
}

}

#endif
