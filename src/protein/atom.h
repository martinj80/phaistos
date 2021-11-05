// atom.h --- Atom class
// Copyright (C) 2006-2011 Wouter Boomsma
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


#ifndef ATOM_H
#define ATOM_H

#include "utils/vector_matrix_3d.h"
#include "definitions.h"

namespace phaistos {

class Residue;

//! Atom class.

//! Atoms have a double representation. Each atom is parameterized
//! using an angle, a dihedral and a bondlength (internal coordinates),
//! but also contains a position. These representations are both maintained, 
//! and can be checked to be consistent using check_consistency(). 
class Atom {
public:

     //! Type of atom: N,CA,C,...
     definitions::AtomEnum atom_type;

     //! Residue in which atom is located
     Residue *residue;

     //! Index in enclosing resiue
     int index;

     //! Additional identifier used by some energy functions
     int biotype;               
     
     //! Bond angle for atom i defined by atoms: i-1, i, i+1.
     double *angle;		

     //! dihedral for atom i defined by atoms: i-2, i-1, i, i+1.
     double *dihedral;  

     //! length of bond to previous atom in chain        
     double bond_length;

     //! 3D coordinate
     Vector_3D position;

     //! Mass of atom
     double mass;

     //! Whether the atom belongs to the residue sidechain
     bool is_sidechain_atom;

     //! indicates that the atom belongs to the backbone
     bool is_backbone_atom; 	

     //! Whether angle and dihedral reflect the position of the atom itself. 
     //! Normally, the atom angle at position i is defined by the
     //! position of atoms i-1, i and i+1, and the dihedral at position
     //! i is defined by positions i-2, i-1, i and i+1. In some cases
     //! (e.g. in sidechains) it is more convenient to let the atom
     //! keep track of the angle and dihedral necessary to position
     //! itself, i.e., the angle of atom i is defined by position i-2,
     //! i-1, i, and the dihedral by i-3, i-2, i-1, i. In this latter
     //! case, self_maintained_positioning=true     
     bool self_maintained_positioning;                                      

     //! Whether the atom is responsible for cleaning up dihedral pointer
     bool owner_of_dihedral;

     //! Whether the atom is responsible for cleaning up angle pointer
     bool owner_of_angle;
     
     //! Neighbour cache of first 3 neighbours in various iteration modes
     Atom *neighbour_atom[7][definitions::ITERATE_ENUM_SIZE];

     //! Vector of all covalent neighbours (atom_type, residue-offset)
     std::vector<std::pair<definitions::AtomEnum, int> > covalent_neighbours;
     

     //! Constructor - incomplete initialization - angles and position must be initialized later
     //!
     //! \param atom_type Type of atom (N,CA,C,...).
     //! \param residue Pointer to enclosing residue.
     //! \param index index in enclosing residue.
     Atom(definitions::AtomEnum atom_type, Residue *residue, int index);

     //! Constructor - incomplete initialization - angles must be initalized later
     //!
     //! \param atom_type Type of atom (N,CA,C,...).
     //! \param residue Pointer to enclosing residue.
     //! \param index index in enclosing residue.
     //! \param position 3D-coordinate
     Atom(definitions::AtomEnum atom_type, Residue *residue, int index, Vector_3D position);

     //! Constructor - incomplete initialization - position must be initalized later
     //!
     //! \param atom_type Type of atom (N,CA,C,...).
     //! \param residue Pointer to enclosing residue.
     //! \param index index in enclosing residue.
     //! \param angle bondangle.
     //! \param dihedral dihedral angle.
     Atom(definitions::AtomEnum atom_type, Residue *residue, int index, double angle, double dihedral);

     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     Atom(Atom &other);

     //! Copy constructor - with different enclosing residue.
     //!
     //! \param other Source object from which copy is made.
     //! \param residue Pointer to enclosing residue.
     Atom(Atom &other, Residue *residue);

     //! Destructor
     ~Atom();
     
     //! Initializer
     //!
     //! \param atom_type Type of atom (N,CA,C,...).
     //! \param residue Pointer to enclosing residue.
     //! \param index Index in enclosing residue.
     void init(definitions::AtomEnum atom_type, Residue *residue, int index);

     //! Initialize bond lengths
     void init_bond_length();

     //! Initialize bond angles and dihedrals
     void init_angles();

     //! Determine the atoms defining current angle
     //!
     //! \param atom1 Destination for atom1
     //! \param atom2 Destination for atom2
     //! \param atom3 Destination for atom3
     void get_angle_atoms(Atom **atom1, Atom **atom2, Atom **atom3);

     //! Determine the atoms defining current dihedral
     //!
     //! \param atom1 Destination for atom1
     //! \param atom2 Destination for atom2
     //! \param atom3 Destination for atom3
     //! \param atom4 Destination for atom4
     void get_dihedral_atoms(Atom **atom1, Atom **atom2, Atom **atom3, Atom **atom4);
     
     //! Constants necessary for determining neighbour
     //!
     //! \param offset Distance along chain to neighbour (nth neighbour)
     //! \param iteration_type Type of iteration
     //! \return Pointer to neighbouring atom
     Atom *get_neighbour_constants(int &offset, definitions::IterateEnum iteration_type);

     //! Get neighouring atom (internal method)
     //!
     //! \param offset Distance along chain to neighbour (nth neighbour)
     //! \param iteration_type Type of iteration
     //! \return Pointer to neighbouring atom
     Atom *get_neighbour_internal(int offset, definitions::IterateEnum iteration_type=definitions::ALL);

     //! Get neighouring atom (use cache if available)
     //!
     //! \param offset Distance along chain to neighbour (nth neighbour)
     //! \param iteration_type Type of iteration
     //! \return Pointer to neighbouring atom
     Atom *get_neighbour(int offset, definitions::IterateEnum iteration_type=definitions::ALL);

     //! Faster version of get_neighbour. In this version, the offset is passed
     //! along as a template parameter, which makes it possible to 
     //! specialize for offset=-3,-2,-1,1,2,3
     //!
     //! \tparam OFFSET Distance along chain to neighbour (nth neighbour)
     //! \param iteration_type Type of iteration
     //! \param direction Iteration direction along the chain
     //! \return Pointer to neighbouring atom
     template <int OFFSET>
     Atom *get_neighbour(definitions::IterateEnum iteration_type=definitions::ALL, int direction=1);

     //! Specialized version of get_neighbour for offset=-3,-2,-1,1,2,3.
     //!
     //! \tparam OFFSET Distance along chain to neighbour (nth neighbour)
     //! \param iteration_type Type of iteration
     //! \param direction Iteration direction along the chain
     //! \return Pointer to neighbouring atom
     template <int OFFSET>
     Atom *get_neighbour_cached(definitions::IterateEnum iteration_type=definitions::ALL, int direction=1);

     //! Clear the cached information on neighbouring atoms.
     void clear_neighbour_cache();

     //! Update bondangle value based on current positions.
     //! NOTE: Should only be called when the position of left and right atom neighbours are initialized
     //!
     //! \param fix_end_points If true, the angle value at the non-well-defined 
     //!        endpoint positions are set to zero.
     void update_angle(bool fix_end_points=true);

     //! Update dihedral value based on current positions.
     //! NOTE: Should only be called when the position of 2xleft and 1xright atom neighbours are initialized
     //!
     //! \param fix_end_points If true, the angle value at the non-well-defined 
     //!        endpoint positions are set to zero.
     void update_dihedral(bool fix_end_points=true);

     //! Update position value based on current angle values
     //!
     //! \param direction Iteration direction along the chain
     void update_position(int direction);

     //! Get bondangle value 
     //!
     //! \return bondangle value
     double &get_angle();

     //! Set bondangle value
     //!
     //! \param angle bondangle value
     void set_angle(double angle);

     //! Get dihedral value 
     //!
     //! \return dihedral value
     double &get_dihedral();

     //! Set dihedral value 
     //!
     //! \param dihedral dihedral value     
     void set_dihedral(double dihedral);

     //! Get the atom's mass
     //!
     //! \return mass value
     double get_mass();

     //! Apply translation operation to position
     //!
     //! \param translation translation vector
     void translate(Vector_3D translation);

     //! Apply rotation operation to position
     //!
     //! \param rotation matrix
     void rotate(Matrix_3D rotation);

     //! Check consistency between the two representations of the atom. 
     //! The chain is overspecified since it keeps track of both dihedral angles and positions.
     //! This method checks that the two representations are consistent.
     void check_consistency();

     //! Return distance along chain from atom1 to atom2 
     //!
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \return Distance along chain
     friend int chain_distance(Atom *atom1, Atom *atom2);
     
     //! Output in PDB format
     //!
     //! \param counter_offset Offset to add to atom index
     //! \param chain_number The chain id
     //! \param b_factor_string Optional string to output as b factor
     //! \return Output in pdb format
     std::string output_as_pdb(int counter_offset, const int chain_number=1, const std::string &b_factor_string="");

     //! Overload << operator for atom reference
     friend std::ostream &operator<<(std::ostream &o, Atom &a);

     //! Overload << operator for atom pointer
     friend std::ostream &operator<<(std::ostream &o, Atom *a);

};

}

#include "residue.h"

namespace phaistos {

//! Return distance along chain from atom1 to atom2 
//!
//! \param atom1 First atom
//! \param atom2 Second atom
//! \return Distance along chain
template <typename CHAIN_TYPE>
inline int chain_distance(Atom *atom1, Atom *atom2) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     Residue *res1 = atom1->residue;
     Residue *res2 = atom2->residue;
     
     // Determine whether atoms are in same residue
     if (res1 == res2) {
          return res1->chain_distance<CHAIN_TYPE>(atom1, atom2);
     } else if (res1->index < res2->index) {
          int d1 = res1->chain_distance<CHAIN_TYPE>(atom1, (*res1)[C]);
          int d2 = res2->chain_distance<CHAIN_TYPE>(atom2, (*res2)[N]);
          // The number of backbone atoms per residue is assumed to be constant
          // return (res2->index - res1->index - 1)*res1->backboneAtoms+1+d1+d2;
          return (res2->index - res1->index - 1)*(res1->iteration_range_indices[BACKBONE][Residue::END]+1)+1+d1+d2;
     } else {
          int d1 = res2->chain_distance<CHAIN_TYPE>(atom2, (*res2)[C]);
          int d2 = res1->chain_distance<CHAIN_TYPE>(atom1, (*res1)[N]);
          // The number of backbone atoms per residue is assumed to be constant
          // return (res1->index - res2->index - 1)*res1->backboneAtoms+1+d1+d2;
          return (res1->index - res2->index - 1)*(res1->iteration_range_indices[BACKBONE][Residue::END]+1)+1+d1+d2;
     }
}

//! Faster version of get_neighbour. In this version, the offset is passed
//! along as a template parameter, which makes it possible to 
//! specialize for offset=-3,-2,-1,1,2,3
//!
//! \tparam OFFSET Distance along chain to neighbour (nth neighbour)
//! \param iteration_type Type of iteration
//! \param direction Iteration direction along the chain
//! \return Pointer to neighbouring atom
template <int OFFSET>
inline Atom *Atom::get_neighbour_cached(definitions::IterateEnum iteration_type, int direction) {

     int offset = OFFSET*direction;
     if (neighbour_atom[offset+3][iteration_type] == NULL ) {
	  neighbour_atom[offset+3][iteration_type] = get_neighbour_internal(offset, iteration_type);
     }

     // Return cached value
     return neighbour_atom[OFFSET+3][iteration_type];
}

//! Faster version of get_neighbour. In this version, the offset is passed
//! along as a template parameter, which makes it possible to 
//! specialize for offset=-3,-2,-1,1,2,3
//!
//! \tparam OFFSET Distance along chain to neighbour (nth neighbour)
//! \param iteration_type Type of iteration
//! \param direction Iteration direction along the chain
//! \return Pointer to neighbouring atom
template <int OFFSET>
inline Atom *Atom::get_neighbour(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_internal(OFFSET*direction, iteration_type);
}

//! get_neighbour specialization for offset=+1
template<>
inline Atom *Atom::get_neighbour<+1>(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_cached<+1>(iteration_type, direction);
}

//! get_neighbour specialization for offset=-1
template<>
inline Atom *Atom::get_neighbour<-1>(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_cached<-1>(iteration_type, direction);
}

//! get_neighbour specialization for offset=+2
template<>
inline Atom *Atom::get_neighbour<+2>(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_cached<+2>(iteration_type, direction);
}

//! get_neighbour specialization for offset=-2
template<>
inline Atom *Atom::get_neighbour<-2>(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_cached<-2>(iteration_type, direction);
}

//! get_neighbour specialization for offset=+3
template<>
inline Atom *Atom::get_neighbour<+3>(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_cached<+3>(iteration_type, direction);
}

//! get_neighbour specialization for offset=-3
template<>
inline Atom *Atom::get_neighbour<-3>(definitions::IterateEnum iteration_type, int direction) {
     return get_neighbour_cached<-3>(iteration_type, direction);
}

}

#endif
