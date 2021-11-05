// dof_iterator.h --- Iterators over degrees-of-freedom
// Copyright (C) 2006-2010 Wouter Boomsma
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

#ifndef DOF_ITERATOR_H
#define DOF_ITERATOR_H

#include "iterator_base.h"

namespace phaistos {

//! Iterator over degrees-of-freedom
template <typename CHAIN_TYPE>
class DofIterator: public IteratorBase<double, DofIterator<CHAIN_TYPE> > {
};

//! Iterator over degrees-of-freedom.
//! Specialization for ChainFB (currently the only available implementation)
template <>
class DofIterator<ChainFB>: public IteratorBase<double, DofIterator<ChainFB> > {

protected:
     //! Which of the potential degrees of freedom are free - internally in chain
     bool dof[3][2];

     //! Which of the potential degrees of freedom are free - at N terminal
     bool dof_n_terminus[3][2];

     //! Which of the potential degrees of freedom are free - at C terminal
     bool dof_c_terminus[3][2];

     //! Current atom pointer     
     Atom *atom_current;

     //! End of iteration pointer
     Atom *atom_end;

     //! Backup: Current atom pointer     
     Atom *atom_current_backup;

     //! Current dof index (index into dof vector)
     int dof_index;

     //! End of iteration dof index (index into dof vector)
     int dof_index_end;
     
     //! Whether to include sidechain dofs
     bool include_sidechain;

     //! Sidechain degrees of freedom
     std::vector<std::pair<Atom *,definitions::AngleEnum> > *chi_atoms;

     //! Move in forward direction
     //! continues until free dof is found
     void increment() {

          // Test if we are currently iterating over a sidechain
          if (chi_atoms) {
               dof_index++;
               // Test whether we are outside chi angle bounds
               if ((dof_index-2) >= (int)chi_atoms->size()) {
                    chi_atoms = NULL;
                    dof_index = definitions::DIHEDRAL;
                    atom_current = atom_current_backup;

               // if not, simply update current atom pointer
               } else {
                    atom_current = (*this->chi_atoms)[dof_index-2].first;
               }

          // If we are currently at a dihedral, the next dof is an angle 
          } else if (dof_index == definitions::DIHEDRAL) {
               dof_index = definitions::ANGLE;

          // If we are currently at an ANGLE, move to the next atom 
          } else {

               // Find next atom along backbone
               Atom *next_atom = atom_current->get_neighbour<+1>(definitions::BACKBONE);

               // Test whether we reached the end of a residue
               // if so, check whether to start side chain iteration 
               if (include_sidechain && (!next_atom || (next_atom->residue != atom_current->residue)) && 
                   (atom_current->residue->chi_atoms.size() > 0)) {

                    // Retrieve sidechain chi angles for current residue
                    dof_index = definitions::CHI1;
                    this->chi_atoms = &atom_current->residue->chi_atoms;
                    atom_current = (*this->chi_atoms)[dof_index-2].first;
                    atom_current_backup = next_atom;

               // if no sidechain iteration, check if next_atom exists
               } else if (next_atom) {
                    atom_current = next_atom;
                    dof_index = definitions::DIHEDRAL;
               } else {
                    atom_current = NULL;
                    dof_index = definitions::DIHEDRAL;
                    return;
               }
          }

          // If current degree of freedom is not free, move on
          if (atom_current && !is_free()) {
               increment();
          }
     }

     //! Move in reverse direction
     //! continues until free dof is found
     void decrement() {

          // Test if we are currently iterating over a sidechain
          if (chi_atoms) {

               dof_index--;

               // Test whether we are outside chi angle bounds
               if (dof_index < 2) {
                    chi_atoms = NULL;
                    dof_index = definitions::ANGLE;
                    atom_current = atom_current_backup;

               // if not, simply update current atom pointer
               } else {
                    atom_current = (*this->chi_atoms)[dof_index-2].first;
               }

          // If we are currently at an angle, the next dof is a dihedral 
          } else if (dof_index == definitions::ANGLE) {
               dof_index = definitions::DIHEDRAL;

          // If we are currently at a dihedral, move to the next atom 
          // in reverse direction 
          } else {

               // Get next atom in reverse direction along backbone
               Atom *next_atom = atom_current->get_neighbour<-1>(definitions::BACKBONE);

               // Check whether next atom exists
               if (next_atom) {

                    // Test whether we reached the beginning of a residue
                    // if so, check whether to start side chain iteration 
                    if (include_sidechain &&
                        (next_atom->residue != atom_current->residue) &&
                        (next_atom->residue->chi_atoms.size() > 0)) {
                    
                         atom_current = next_atom;
                         this->chi_atoms = &atom_current->residue->chi_atoms;
                         dof_index = 2 + this->chi_atoms->size()-1;
                         atom_current = (*this->chi_atoms)[dof_index-2].first;
                         atom_current_backup = next_atom;

                    // Otherwise, just update atom_current 
                    } else {
                         atom_current = next_atom;
                         dof_index = definitions::ANGLE;
                    }

               } else {
                    atom_current = NULL;
                    dof_index = definitions::ANGLE;
                    return;
               }	  
          }

          // If current degree of freedom is not free, move on
          if (atom_current && !is_free()) {
               decrement();
          }          
     }

     //! Check whether current dof is free
     //!
     //! \return True if potential degree of freedom is free
     // Note: this function has external requirement not available at
     // this point in the code. It is included at the end of the file
     bool is_free();


public:

     //! Number of times this iterator has been incremented
     int offset;  

     //! Specific angles types
     enum AngleSelectionEnum{NONE=0,
                             N_ANGLE=1,  N_DIHEDRAL=2,
                             CA_ANGLE=4, CA_DIHEDRAL=8,
                             C_ANGLE=16, C_DIHEDRAL=32,
                             NTERM_N_ANGLE=64,   NTERM_N_DIHEDRAL=128,
                             NTERM_CA_ANGLE=256, NTERM_CA_DIHEDRAL=512,
                             NTERM_C_ANGLE=1024, NTERM_C_DIHEDRAL=2048,
                             CTERM_N_ANGLE=4096,   CTERM_N_DIHEDRAL=8192,
                             CTERM_CA_ANGLE=16384, CTERM_CA_DIHEDRAL=32768,
                             CTERM_C_ANGLE=65536, CTERM_C_DIHEDRAL=131072,
                             CHI_ANGLES=262144,
                             STANDARD_DOFS=(N_ANGLE +
                                            CA_ANGLE + CA_DIHEDRAL +
                                            C_ANGLE + C_DIHEDRAL +
                                            NTERM_CA_ANGLE +
                                            NTERM_C_ANGLE + NTERM_C_DIHEDRAL +
                                            CTERM_N_ANGLE +
                                            CTERM_CA_ANGLE + CTERM_CA_DIHEDRAL),
                             BONDANGLE_DOFS=(N_ANGLE + CA_ANGLE + C_ANGLE +
                                             NTERM_CA_ANGLE + NTERM_C_ANGLE +
                                             CTERM_N_ANGLE + CTERM_CA_ANGLE),
                             DIHEDRAL_DOFS=(CA_DIHEDRAL + C_DIHEDRAL +
                                            NTERM_C_DIHEDRAL +
                                            CTERM_CA_DIHEDRAL)};
     
     //! Overload + for AngleSelectionEnum
     friend AngleSelectionEnum operator+(AngleSelectionEnum v1, AngleSelectionEnum v2) {
          return AngleSelectionEnum((int)v1 | (int)v2);
     }

     //! Overload += for AngleSelectionEnum
     friend AngleSelectionEnum &operator+=(AngleSelectionEnum &v1, AngleSelectionEnum v2) {
          v1 = AngleSelectionEnum((int)v1 | (int)v2);
          return v1;
     }

     //! Overload - for AngleSelectionEnum
     friend AngleSelectionEnum operator-(AngleSelectionEnum v1, AngleSelectionEnum v2) {
          return AngleSelectionEnum((int)v1 ^ (int)v2);
     }     

     //! Overload -= for AngleSelectionEnum
     friend AngleSelectionEnum operator-=(AngleSelectionEnum &v1, AngleSelectionEnum v2) {
          v1 = AngleSelectionEnum((int)v1 ^ (int)v2);
          return v1;
     }     


     //! Initializer
     //!
     //! \param active_dofs Specifies with of the potential degrees of freedom are active
     void init(AngleSelectionEnum active_dofs) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          chi_atoms = NULL;

          // Keeps track of how many incremental steps have been done
          offset = 0;

          // Set degrees of freedom for internal residues
          this->dof[N][DIHEDRAL]  = active_dofs & N_DIHEDRAL;
          this->dof[N][ANGLE]     = active_dofs & N_ANGLE;
          this->dof[CA][DIHEDRAL] = active_dofs & CA_DIHEDRAL;
          this->dof[CA][ANGLE]    = active_dofs & CA_ANGLE;
          this->dof[C][DIHEDRAL]  = active_dofs & C_DIHEDRAL;
          this->dof[C][ANGLE]     = active_dofs & C_ANGLE;   

          // Set degrees of freedom for N-terminus
          this->dof_n_terminus[N][DIHEDRAL]  = active_dofs & NTERM_N_DIHEDRAL;
          this->dof_n_terminus[N][ANGLE]     = active_dofs & NTERM_N_ANGLE;
          this->dof_n_terminus[CA][DIHEDRAL] = active_dofs & NTERM_CA_DIHEDRAL;
          this->dof_n_terminus[CA][ANGLE]    = active_dofs & NTERM_CA_ANGLE;
          this->dof_n_terminus[C][DIHEDRAL]  = active_dofs & NTERM_C_DIHEDRAL;
          this->dof_n_terminus[C][ANGLE]     = active_dofs & NTERM_C_ANGLE;   
 
          // Set degrees of freedom for C-terminus
          this->dof_c_terminus[N][DIHEDRAL]  = active_dofs & CTERM_N_DIHEDRAL;
          this->dof_c_terminus[N][ANGLE]     = active_dofs & CTERM_N_ANGLE;
          this->dof_c_terminus[CA][DIHEDRAL] = active_dofs & CTERM_CA_DIHEDRAL;
          this->dof_c_terminus[CA][ANGLE]    = active_dofs & CTERM_CA_ANGLE;
          this->dof_c_terminus[C][DIHEDRAL]  = active_dofs & CTERM_C_DIHEDRAL;
          this->dof_c_terminus[C][ANGLE]     = active_dofs & CTERM_C_ANGLE;   
    
          this->include_sidechain = active_dofs & CHI_ANGLES;
     }


     //! Constructors - using atom
     //!
     //! \param atom_start Atom at which to start     
     //! \param dof_index_begin Degree of freedom index with which to begin
     //! \param active_dofs Selection of active degrees-of-freedom
     DofIterator(Atom *atom_start, int dof_index_begin = definitions::DIHEDRAL,
                 AngleSelectionEnum active_dofs=STANDARD_DOFS)
          : atom_current(atom_start),
            atom_end(NULL),
            dof_index(dof_index_begin),
            dof_index_end(-1) {

          // Set which variables are degrees of freedom
          init(active_dofs);

          // Start in sidechain if dof_index is larger than ANGLE
          if (include_sidechain && dof_index>definitions::ANGLE)
               this->chi_atoms = &atom_current->residue->chi_atoms;

          // Move forward until first degree of freedom
          if (atom_current && !is_free())
               increment();
     }
      
     //! Constructor - from chain
     //!
     //! \param chain Molecule chain
     // Note: this function has external requirement not available at
     // this point in the code. It is included at the end of the file
     DofIterator(const ChainFB &chain);

     //! Constructors - from atoms and dof indices
     //! \param atom_start Atom at which to start     
     //! \param dof_index_begin Degree of freedom index with which to begin
     //! \param atom_end Atom at which to end
     //! \param dof_index_end Degree of freedom index with which to end
     //! \param active_dofs Selection of active degrees-of-freedom
     DofIterator(Atom *atom_start, int dof_index_begin,
                 Atom *atom_end, int dof_index_end=definitions::DIHEDRAL,
                 AngleSelectionEnum active_dofs=STANDARD_DOFS)
          : atom_current(atom_start),
            atom_end(atom_end),
            dof_index(dof_index_begin),
            dof_index_end(dof_index_end) {

          // Set which variables are degrees of freedom
          init(active_dofs);
          
          // Start in sidechain if dof_index is larger than ANGLE
          if (include_sidechain && dof_index>definitions::ANGLE)
               this->chi_atoms = &atom_current->residue->chi_atoms;
          
          // Move forward until first degree of freedom
          if (atom_current && !is_free())
               increment();
     }

     //! Constructor - from residue
     //!
     //! \param r Residue pointer
     //! \param dof_index Degree of freedom index with which to begin
     //! \param active_dofs Selection of active degrees-of-freedom
     DofIterator(Residue *r, int dof_index = definitions::DIHEDRAL, 
                 AngleSelectionEnum active_dofs=STANDARD_DOFS) {

          // Start at first atom in residue
          if (r)
               atom_current = r->atoms[0];
          else
               atom_current = NULL;

          // Set which variables are degrees of freedom
          init(active_dofs);

          // Starting dof type
          this->dof_index = dof_index;

          // Start in sidechain if dof_index is larger than ANGLE
          if (include_sidechain && dof_index>definitions::ANGLE)
               this->chi_atoms = &atom_current->residue->chi_atoms;
     
          // Move forward until first degree of freedom
          if (atom_current && !is_free())
               increment();
     }
    
     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DofIterator& operator=(const DofIterator& other) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          atom_current = other.atom_current;
          dof_index = other.dof_index;
          offset = other.offset;
          include_sidechain = other.include_sidechain;

          this->dof[N][DIHEDRAL] = other.dof[N][DIHEDRAL];
          this->dof[N][ANGLE] = other.dof[N][ANGLE];
          this->dof[CA][DIHEDRAL] = other.dof[CA][DIHEDRAL];
          this->dof[CA][ANGLE] = other.dof[CA][ANGLE];
          this->dof[C][DIHEDRAL] = other.dof[C][DIHEDRAL];
          this->dof[C][ANGLE] = other.dof[C][ANGLE];
          this->dof_n_terminus[N][DIHEDRAL] = other.dof_n_terminus[N][DIHEDRAL];
          this->dof_n_terminus[N][ANGLE] = other.dof_n_terminus[N][ANGLE];
          this->dof_n_terminus[CA][DIHEDRAL] = other.dof_n_terminus[CA][DIHEDRAL];
          this->dof_n_terminus[CA][ANGLE] = other.dof_n_terminus[CA][ANGLE];
          this->dof_n_terminus[C][DIHEDRAL] = other.dof_n_terminus[C][DIHEDRAL];
          this->dof_n_terminus[C][ANGLE] = other.dof_n_terminus[C][ANGLE];
          this->dof_c_terminus[N][DIHEDRAL] = other.dof_c_terminus[N][DIHEDRAL];
          this->dof_c_terminus[N][ANGLE] = other.dof_c_terminus[N][ANGLE];
          this->dof_c_terminus[CA][DIHEDRAL] = other.dof_c_terminus[CA][DIHEDRAL];
          this->dof_c_terminus[CA][ANGLE] = other.dof_c_terminus[CA][ANGLE];
          this->dof_c_terminus[C][DIHEDRAL] = other.dof_c_terminus[C][DIHEDRAL];
          this->dof_c_terminus[C][ANGLE] = other.dof_c_terminus[C][ANGLE];

          return(*this);
     }

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(const DofIterator& other) const {
          return((atom_current == other.atom_current) && (dof_index == other.dof_index));
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     bool operator>(const DofIterator& other) const {
          if (!atom_current) {
               if (!other.atom_current) {
                    return false;
               } else {
                    return true;
               }
          } else {
               if (!other.atom_current) {
                    return false;
               } else if (atom_current->residue->index != other.atom_current->residue->index) {
                    return (atom_current->residue->index > other.atom_current->residue->index);
               } else {
                    if (atom_current->index != other.atom_current->index)
                         return (atom_current->index > other.atom_current->index);
                    else
                         return (dof_index > other.dof_index);
               }
          }
     }

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     bool operator<(const DofIterator& other) const {
          if (!atom_current) {
               if (!other.atom_current) {
                    return false;
               } else {
                    return false;
               }
          } else {
               if (!other.atom_current) {
                    return true;
               } else if (atom_current->residue->index != other.atom_current->residue->index) {
                    return (atom_current->residue->index < other.atom_current->residue->index);
               } else {
                    if (atom_current->index != other.atom_current->index)
                         return (atom_current->index < other.atom_current->index);
                    else
                         return (dof_index < other.dof_index);
               }
          }
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     DofIterator &operator++() {
          offset++;
          increment();
          return (*this);
     }

     //! Increment with value
     //!
     //! \param v Value to increment with
     //! \return Current iterator (this)
     DofIterator &operator+=(const int v) {
          for (int i=0; i<v; i++) {
               ++(*this);
          }
          return (*this);
     }

     //! Increment with value - return new iterator
     //!
     //! \param v Value to increment with
     //! \return New iterator
     DofIterator operator+(const int v) const {
          DofIterator dof_iterator(*this);
          dof_iterator += v;
          return dof_iterator;
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     DofIterator &operator--() {
          offset--;
          decrement();
          return (*this);
     }

     //! Increment with value
     //!
     //! \param v Value to increment with
     //! \return Current iterator (this)
     DofIterator &operator-=(const int v) {
          for (int i=0; i<v; i++) {
               --(*this);
          }
          return (*this);
     }

     //! Decrement with value - return new iterator
     //!
     //! \param v Value to increment with
     //! \return New iterator
     DofIterator operator-(const int v) const {
          DofIterator dof_iterator(*this);
          dof_iterator -= v;
          return dof_iterator;
     }


     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     // Implemented in different specializations above
     double &operator*() const {
          if (dof_index == definitions::ANGLE)
               return atom_current->get_angle();
          else {
               return atom_current->get_dihedral();
          }
     }

     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     bool end() const {
          return (atom_current==NULL || (atom_current == atom_end && dof_index == dof_index_end));
     }     


     //! Return type of degree of freedom (angle or dihedral)
     //!
     //! \return Type of current degree of freedom (ANGLE|DIHEDRAL|CHI1|CHI2|CHI3|CHI4)
     definitions::AngleEnum get_dof_type() {
          return definitions::AngleEnum(dof_index);
     }
     
     //! Return current atom
     //!
     //! \return Current atom pointer
     Atom *get_atom() const {
          return atom_current;
     }

     //! Return current residue
     //!
     //! \return Current residue pointer
     Residue *get_residue() const {
          if (atom_current)
               return atom_current->residue;
          else
               return NULL;
     }

     //! Return current residue
     //!
     //! \return index of current residue
     int get_residue_index() const {
          if (atom_current)
               return atom_current->residue->index;
          else
               return -1;
     }
     
     //! Return position vector corresponding to current degree of freedom
     //!
     //! \return 3D-coordinate
     Vector_3D get_position_vector() {
          return atom_current->position;
     }

     //! Return rotation vector corresponding to the current degree of freedom
     //!
     //! \return rotation axis
     Vector_3D get_rotation_vector() {
          Vector_3D ret;
          Atom *atom_prev = atom_current->get_neighbour(-1, definitions::BACKBONE);
          Atom *atom_next = atom_current->get_neighbour(+1, definitions::BACKBONE);
          if (dof_index == definitions::DIHEDRAL) {
               ret = atom_current->position - atom_prev->position;
          } else {
               Vector_3D v1 = atom_current->position - atom_prev->position;
               Vector_3D v2 = atom_next->position - atom_current->position;
               ret = v2%v1;
          }

          return ret.normalize();
     }

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, DofIterator &it) {
          const char *angle_types[6] = {"DIHEDRAL", "ANGLE", "CHI1", "CHI2", "CHI3", "CHI4"};
          o << "Type: " << angle_types[it.dof_index] << "\t    Value: " << *it << "\t    Atom: " << it.get_atom();
          return o;
     }
};

class Dof {
public:
     //! Atom pointer
     Atom *atom;

     //! Type of Degree of Freedom
     definitions::AngleEnum dof_type;

     //! Current value
     double *value;

     //! Constructor
     Dof(Atom *atom, definitions::AngleEnum dof_type, double &value)
          : atom(atom), dof_type(dof_type), value(&value) {}
};

}

// The code below requires complete types of chains

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"

namespace phaistos {

//! Constructor - from chain
//! \param chain Molecule chain
DofIterator<ChainFB>::DofIterator(const ChainFB &chain) : 
     atom_current(chain[0].get_first_atom(definitions::ALL)),
     atom_end(NULL),
     dof_index(definitions::DIHEDRAL),
     dof_index_end(definitions::DIHEDRAL) {
     
     // Set which variables are degrees of freedom
     init(STANDARD_DOFS);
     
     // Move forward until first degree of freedom
     if (atom_current && !is_free())
          increment();          
}



//! DofIterator: Check whether current potential degree-of-freedom is free
//!
//! \return True if potential degree of freedom is free.
inline bool DofIterator<ChainFB>::is_free() {
     if (dof_index<2) {
          if (atom_current->residue->index==0) {
               return (dof_n_terminus[atom_current->atom_type][dof_index]);
          } else if (atom_current->residue->index==((int)(((ResidueFB *)(atom_current->residue))->chain->size())-1)) {
               return dof_c_terminus[atom_current->atom_type][dof_index];
          } else {
               return dof[atom_current->atom_type][dof_index];
          }
     } else {
          return true;
     }
}

}

#endif
