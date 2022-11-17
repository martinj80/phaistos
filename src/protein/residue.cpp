// residue.cpp --- residue class
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

#include "residue.h"
#include "chain.h"
#include "chain_ca.h"
#include "chain_fb.h"
#include "protein/protein_ideal_values.h"

namespace phaistos {

// Import protein definitions (such as residue names)
using namespace definitions;

////////////////////////////
//// Residue base class ////
////////////////////////////

// Constructor
Residue::Residue(ResidueEnum residue_type, int index, int index_res_seq) {
     this->init();
     this->residue_type = residue_type;
     this->index = index;

     this->index_res_seq = index_res_seq;
     if (!is_initialized(index_res_seq))
          this->index_res_seq = index;          

     for (int i=0; i<ATOM_ENUM_SIZE; i++) {
          this->atom_index[i] = -1;
     }

     for (int i=0; i<ITERATE_ENUM_SIZE; ++i) {
          this->iteration_range_indices[i][0] = 0;
          this->iteration_range_indices[i][1] = 0;
     }
}

// Copy constructor
Residue::Residue(Residue &other, int index, int index_res_seq) {

     this->init();
     this->residue_type = other.residue_type;
     this->index = index;

     this->index_res_seq = index_res_seq;
     if (!is_initialized(index_res_seq))
          this->index_res_seq = index;          

     this->sidechain_status = other.sidechain_status;

     this->terminal_status = other.terminal_status;

     for (int i=0; i<ATOM_ENUM_SIZE; i++) {
          this->atom_index[i] = other.atom_index[i];
     }

     for (int i=0; i<ITERATE_ENUM_SIZE; ++i) {
          this->iteration_range_indices[i][0] = other.iteration_range_indices[i][0];
          this->iteration_range_indices[i][1] = other.iteration_range_indices[i][1];
     }

     for (Iterator it = other.atoms.begin(); it != other.atoms.end(); ++it) {
          atoms.push_back(new Atom(**it, this));
     }

}

// Destructor
Residue::~Residue() {
     for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          delete *it;
     }
     for (Residue::Iterator it = atoms_inactive.begin(); it != atoms_inactive.end(); ++it) {
          delete *it;
     }
}

// Return size
int Residue::size() const {
     return atoms.size();
}

// Overload indexing operator (using AtomEnum)
Atom *Residue::operator[](const AtomEnum atom_type) {
     return atoms.at(atom_index[atom_type]);
}

// Overload indexing operator (using index)
Atom *Residue::operator[](const unsigned int index) {
     return atoms.at(index);
}

// Determines whether atom is present
bool Residue::has_atom(const AtomEnum atom_type) {
     return (atom_index[atom_type]>=0);
}

// true if the sidechain has been initialized properly
bool Residue::has_sidechain() {
     return this->sidechain_status;
}

void Residue::init() {
     // this->backboneAtoms = 0;
     this->sidechain_status = false;
     this->atoms_length=0;
}

// define a compare function for sort_atoms
inline bool sort_atoms_cmp (const Atom *a1, const Atom *a2) {
     if (!a1)
          return false;
     if (!a2)
          return true;
     return ( (*a1).atom_type < (*a2).atom_type);
}

// Make sure that atom names adhere to the PDB standard
void Residue::standardize_atoms() {

     // Import phaistos::make_vector namespace
     using namespace vector_utils;

     if (terminal_status == NTERM && has_atom(H)) {
          // std::cout << "Renaming H in residue " << index << "(" << residue_type << ")\n";
          rename_atom(H, make_vector(H1,H2,H3),H1);
     }

     switch (residue_type) {
     case GLY:
          // Translate HA to HA2 or HA3
          if (has_atom(HA)) {
               rename_atom(HA, make_vector(HA2, HA3), HA2);
          }
          if (has_atom(HA1)) {
               rename_atom(HA1, make_vector(HA2, HA3), HA2);
          }
          break;
     case GLU:
     case MET:
     case GLN:
          // Translate HB1 to HB2 or HB3
          if (has_atom(HB1)) {
               rename_atom(HB1, make_vector(HB2, HB3), HB2);
          }
          // Translate HG1 to HG2 or HG3
          if (has_atom(HG1)) {
               rename_atom(HG1, make_vector(HG2, HG3), HG2);
          }
          break;
     case CYS:
     case LEU:
     case PHE:
     case ASN:
     case HIS:
     case SER:
     case ASP:
     case TYR:
     case TRP:

          // Translate HB1 to HB2 or HB3
          if (has_atom(HB1)) {
               rename_atom(HB1, make_vector(HB2, HB3), HB2);
          }
          // HD22 is placed by means of HD21 - so if only
          // one is present, rename it
          if (has_atom(HD22) && !has_atom(HD21)) {
               rename_atom(HD22, make_vector(HD21));
          }

          // Remove HE2 hydrogens (later replaced with HD1)
          if (residue_type==HIS && has_atom(HE2)) {
               remove_atom(HE2, true);
          }

          break;
     case ILE:
          if (has_atom(HG11)) {
               rename_atom(HG11, make_vector(HG12, HG13), HG12);
          }
          if (has_atom(HD1)) {
               rename_atom(HD1, make_vector(HD11), HD11);
          }
          if (has_atom(HD2)) {
               rename_atom(HD2, make_vector(HD12), HD12);
          }
          if (has_atom(HD3)) {
               rename_atom(HD3, make_vector(HD13), HD13);
          }
          break;
     case PRO:
     case ARG:
          // Translate HB1 to HB2 or HB3
          if (has_atom(HB1)) {
               rename_atom(HB1, make_vector(HB2, HB3), HB2);
          }
          // Translate HG1 to HG2 or HG3
          if (has_atom(HG1)) {
               rename_atom(HG1, make_vector(HG2, HG3), HG2);
          }
          // Translate HD1 to HD2 or HD3
          if (has_atom(HD1)) {
               rename_atom(HD1, make_vector(HD2, HD3), HD2);
          }
          break;
     case LYS:
          // Translate HB1 to HB2 or HB3
          if (has_atom(HB1)) {
               rename_atom(HB1, make_vector(HB2, HB3), HB2);
          }
          // Translate HG1 to HG2 or HG3
          if (has_atom(HG1)) {
               rename_atom(HG1, make_vector(HG2, HG3), HG2);
          }
          // Translate HD1 to HD2 or HD3
          if (has_atom(HD1)) {
               rename_atom(HD1, make_vector(HD2, HD3), HD2);
          }
          // Translate HE1 to HE2 or HE3
          if (has_atom(HE1)) {
               rename_atom(HE1, make_vector(HE2, HE3), HE2);
          }
          break;          
          // Translate HB1 to HB2 or HB3
          if (has_atom(HB1)) {
               rename_atom(HB1, make_vector(HB2, HB3), HB2);
          }
          break;

     //Modified by MJ: Possibly add more error checking for SEP, TPO, PTR
     case SEP:
         break;
     case TPO:
         break;
     case PTR:
         break;
     default:
          break;
     }
}


// Rename an atom to any available name in possible_replacements
void Residue::rename_atom(definitions::AtomEnum atom_type, 
                          std::vector<definitions::AtomEnum> possible_replacements, 
                          definitions::AtomEnum preferred_type) {

     int index = atom_index[atom_type];
     Atom *atom = (*this)[atom_type];

     // Remove atom index
     atom_index[atom_type] = -1;
     atom->atom_type = XX;

     // First try preferred type
     if (preferred_type != XX && !has_atom(preferred_type)) {

          atom->atom_type = preferred_type;
          atom_index[preferred_type] = index;

          // Initialize bond lengths and angles
          atom->init_bond_length();
          atom->init_angles();

     } else {
     
          // Change name to any available name in possibleReplacements
          for (unsigned int i=0; i<possible_replacements.size(); i++) {
               if (!has_atom(possible_replacements[i])) {

                    // If preferredType is set, move atom at preferredIndex to
                    // available AtomEnum, and atom_type to the preferredType.
                    // This is to allow translations of the type:
                    // HB1, HB2 -> HB2, HB3 instead of HB3, HB2
                    if (preferred_type != XX) {
                         int index2 = atom_index[preferred_type];
                         Atom *atom2 = atoms[index2];

                         // First move atom with preferred atomtype to new atom type
                         atom2->atom_type = possible_replacements[i];
                         atom_index[possible_replacements[i]] = index2;
                         atom2->init_bond_length();
                         atom2->init_angles();
                              
                         // Then translate atom_type to preferredType
                         atom->atom_type = preferred_type;
                         atom_index[preferred_type] = index;
                         atom->init_bond_length();
                         atom->init_angles();

                    } else {
                    
                         atom->atom_type = possible_replacements[i];
                         atom_index[possible_replacements[i]] = index;

                         // Initialize bond lengths and angles
                         atom->init_bond_length();
                         atom->init_angles();
                    }
                    break;
               }
          }
     }

     // If no possible replacement was found, delete the atom
     if (atom->atom_type == XX) {
          assert(false);
          delete atoms[index];
          atoms[index] = NULL;
     } else {
          // Otherwise, reinitialize atom
          atom->init(atom->atom_type, this, index);
          atom->init_bond_length();
     }

}


// Sort the atoms vector according to the standard pdb order and
// remove NULL entries
void Residue::sort_atoms() {

     // Make sure all atoms names are correct according to PDB standard
     standardize_atoms();
     
     // sort the atoms vector 
     std::sort(this->atoms.begin(), this->atoms.end(), sort_atoms_cmp);

     // Update the atom_index vector which keeps track of the positions
     // of the atoms in the atoms vector. Also set the internal index
     // of each atom, and the pointers to where the backbone atoms end
     // and where the sidechain atoms begin
     for (int i=0; i<ITERATE_ENUM_SIZE; ++i) {
          this->iteration_range_indices[i][0] = -1;
          this->iteration_range_indices[i][1] = -1;
     }

     // this->sidechainAtomsBegin = -1;
     // this->sidechainAtomsEnd = 0;
     for ( unsigned int i=0;i<this->atoms.size(); i++) {

          // The deleted atoms are placed last by the sorting method
          if (atoms[i] == NULL) {
               atoms.resize(i);
               break;
          }
          
          // update index
          this->atom_index[(*(this->atoms.at(i))).atom_type] = i;

          // update internal atom index
          (this->atoms.at(i))->index = i;

          // set the pointers to the first and last sidechain atom
          // Trick: casting -1 to an unsigned int turns it a large number
          // therefore - we only have to make a single check
          if ( this->atoms.at(i)->is_sidechain_atom && (i < (unsigned int)this->iteration_range_indices[SC_ONLY][BEGIN])) {
               // sidechainAtomsBegin = i;
               this->iteration_range_indices[SC_ONLY][BEGIN] = i;
          }
          // if ( !this->atoms.at(i)->is_sidechain_atom && i > (unsigned int)sidechainAtomsBegin) {
          // std::cout << "\t" << i << " " << this->atoms.at(i) << " " << this->atoms.at(i)->is_sidechain_atom << " " << ((int)i > this->iteration_range_indices[SC_ONLY][END]) << "\n";
          if ( this->atoms.at(i)->is_sidechain_atom && (int)i > this->iteration_range_indices[SC_ONLY][END]) {
               // sidechainAtomsEnd = i;
               this->iteration_range_indices[SC_ONLY][END] = i;
          }

          // set the pointers to the first and last backbone atom
          if ( this->atoms.at(i)->is_backbone_atom && (i < (unsigned int)iteration_range_indices[BACKBONE][BEGIN])) {
               this->iteration_range_indices[BACKBONE][BEGIN] = i;
          }
          if ( this->atoms.at(i)->is_backbone_atom && (int)i > (int)iteration_range_indices[BACKBONE][END]) {
               this->iteration_range_indices[BACKBONE][END] = i;
          }          
     }
     if (this->iteration_range_indices[SC_ONLY][END] == -1) {
          // sidechainAtomsEnd = this->size()-1;
          this->iteration_range_indices[SC_ONLY][END] = this->size()-1;
     }
     if (this->iteration_range_indices[SC_ONLY][BEGIN] == -1) {
          // sidechainAtomsBegin = this->size()-1;
          this->iteration_range_indices[SC_ONLY][BEGIN] = this->size()-1;
     }

     this->iteration_range_indices[ALL][BEGIN] = 0;
     this->iteration_range_indices[ALL][END] = this->size()-1;
     this->iteration_range_indices[CA_ONLY][BEGIN] = this->atom_index[CA];
     this->iteration_range_indices[CA_ONLY][END] = this->atom_index[CA];

     // Since the order of atoms has been modified, clear all neighbour information
     clear_neighbour_cache();
     // std::cout << index << ": " << this->iteration_range_indices[SC_ONLY][BEGIN] << " " << this->iteration_range_indices[SC_ONLY][END] << "\n";
}

// Sort the atoms_inactive vector according to the standard pdb order
// and remove NULL entries
void Residue::sort_atoms_inactive() {
     // Sort the atoms_inactive vector
     std::sort(this->atoms_inactive.begin(), this->atoms_inactive.end(), sort_atoms_cmp);
     for ( unsigned int i=0;i<this->atoms_inactive.size(); i++) {
          // The deleted atoms are placed last by the sorting method
          if (atoms_inactive[i] == NULL) {
               atoms_inactive.resize(i);
               break;
          }
          
          // update index
          this->atom_index[this->atoms_inactive.at(i)->atom_type] = (i+2)*-1;
          // update internal atom index
          (this->atoms_inactive.at(i))->index = (i+2)*-1;          
     }            
}

// clear all the atoms cached neighbour-information 
void Residue::clear_neighbour_cache() {
     // std::cerr << "\n\nCLEARING NEIGHBOUR CACHE\n\n";
     for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->clear_neighbour_cache();
     }
     // Clear chain distance vector
     atom_chain_distances.clear();
}

// Update angles based on surrounding positions
void Residue::update_angles(bool fix_end_points) {
     for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->update_angle(fix_end_points);
          (*it)->update_dihedral(fix_end_points);
     }
}

// Initialize bond lenghts 
void Residue::init_bond_lengths() {
     for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->init_bond_length();
     }
}

// Initialize bond angle contants 
void Residue::init_angles() {
     for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->init_angles();
     }
}

// Return vector of atom positions
std::vector<Vector_3D> Residue::get_positions() {
     std::vector<Vector_3D> positions;
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          positions.push_back((*it)->position);
     }
     return positions;
}

// Apply translation to all atoms
void Residue::translate(Vector_3D translation) {
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->translate(translation);
     }
}

// Apply rotation to all atoms
void Residue::rotate(Matrix_3D rotation) {
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->rotate(rotation);
     }
}

// Remove atom with atom_type from the chain
void Residue::remove_atom(AtomEnum atom_type, bool execute_sort_atoms) {
     if (this->has_atom(atom_type)) {
          int index = atom_index[atom_type];
          atom_index[atom_type] = -1;
          delete atoms[index];
          atoms[index] = NULL;

          // Sort the atoms vector
          // Sorting must always be done after a deletion/insertion of atoms. However, if several atoms
          // are removed, it is only necessary to sort once. It can therefore be convenient to turn the
          // sort_atoms option off, and do it manually afterwards.
          if (execute_sort_atoms) {
               this->sort_atoms();
          }       
     }
}

// Moves an atom from the atoms_inactive vector to the atoms vector
bool Residue::activate_atom(AtomEnum atom_type, bool execute_sort_atoms) {

     int index = atom_index[atom_type];

     // Check whether atom_type is in an inactive state
     if (index >= 0) {
          // this atom is already active
          return true;
     }
     if (index == -1) {
          // this Atom has never been set.
          return false;
     }

     // Inactive atoms are represented by negative indices starting at -2
     index = (index+2)*-1;

     // Remove atom from inactive atoms vector
     Atom *atom = atoms_inactive[index];
     atoms_inactive[index] = NULL;

     // Add atom to atoms vector
     atoms.push_back(atom);

     // Sort the atoms vector
     // Sorting must always be done after a deletion/insertion of atoms. However, if several atoms
     // are removed, it is only necessary to sort once. It can therefore be convenient to turn the
     // sort_atoms option off, and do it manually afterwards.
     if (execute_sort_atoms) {
          this->sort_atoms();
     }
     return true;
}


// Moves an atom from the atoms vector to the atoms_inactive vector
void Residue::deactivate_atom(AtomEnum atom_type, bool execute_sort_atoms) {
     if (this->has_atom(atom_type)) {

          int index = atom_index[atom_type];
          Atom *atom = atoms[index];

          // Remove atom from atom vector
          atoms[index] = NULL;

          // Add atom to inactive atom vector
          atoms_inactive.push_back(atom);

          // Set atom_index. Inactive atoms are represented by negative indices starting at -2
          index = atoms_inactive.size()-1;
          atom->index = -1*index-2;
          atom_index[atom_type] = atom->index;

          // Sort the atoms vector
          // Sorting must always be done after a deletion/insertion of atoms. However, if several atoms
          // are removed, it is only necessary to sort once. It can therefore be convenient to turn the
          // sort_atoms option off, and do it manually afterwards.
          if (execute_sort_atoms) {
               this->sort_atoms();
          }
     }
}


// Removes a selection of atoms from the residue
void Residue::remove_atoms(AtomSelectionEnum atom_selection, bool execute_sort_atoms) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Remove backbone hydrogens
     if (atom_selection & BACKBONE_H_ATOMS) {
          remove_atom(H,   false);
          remove_atom(H1, false);
          remove_atom(H2, false);
          remove_atom(H3, false);
     }

     // Remove non-backbone hydrogens
     if (atom_selection & NON_BACKBONE_H_ATOMS) {
          for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
               if ((*it) && ((*it)->atom_type>=HA && (*it)->atom_type<=HZ3)) {
                    remove_atom((*it)->atom_type, false);
               }
          }       
     }

     // Remove O atoms
     if (atom_selection & BACKBONE_O_ATOMS) {
          remove_atom(O, false);
          remove_atom(OXT, false);
     }

     // Remove Sidechain atoms
     if (atom_selection & SIDECHAIN_ATOMS) {
          for (Residue::ReverseIterator it = atoms.rbegin(); it != atoms.rend(); ++it) {
               // CB is not disabled together with the other sidechain atoms
               if ((*it) && (*it)->is_sidechain_atom && (*it)->atom_type!=CB) {
                    remove_atom((*it)->atom_type, false);
               }
          }
          sidechain_status = false;
     }

     // Remove CB atom
     if (atom_selection & CB_ATOMS) {
          remove_atom(CB, execute_sort_atoms);
     }
     
     // Sort atom vector
     if (execute_sort_atoms)
          this->sort_atoms();
     
}

// Activates a selection of atoms
void Residue::activate_atoms(AtomSelectionEnum atom_selection, bool execute_sort_atoms) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Activate O atoms
     if (atom_selection & BACKBONE_O_ATOMS) {
          activate_atom(O, false);
          activate_atom(OXT, false);
     }     

     // Activate CB atom
     if (atom_selection & CB_ATOMS) {
          activate_atom(CB, false);
     }

     // Activate backbone hydrogens
     if (atom_selection & BACKBONE_H_ATOMS) {
          activate_atom(H,   false);
          activate_atom(H1, false);
          activate_atom(H2, false);
          activate_atom(H3, false);
     }

     // Activate Sidechain atoms
     if (atom_selection & SIDECHAIN_ATOMS) {
          for (Residue::ReverseIterator it = atoms_inactive.rbegin(); it != atoms_inactive.rend(); ++it) {
               if ((*it) && (*it)->is_sidechain_atom) {
                    activate_atom((*it)->atom_type, false);
                    sidechain_status = true;
               }
          }
     }     
          
     // Activate non-backbone hydrogens
     if (atom_selection & NON_BACKBONE_H_ATOMS) {
          for (Residue::Iterator it = atoms_inactive.begin(); it != atoms_inactive.end(); ++it) {
               if ((*it) && ((*it)->atom_type>=HA && (*it)->atom_type<=HZ3)) {
                    activate_atom((*it)->atom_type, false);
               }
          }       
     }

     if (execute_sort_atoms) {
          this->sort_atoms();
          this->sort_atoms_inactive();
     }
}


// Deactivates a selection of atoms
void Residue::deactivate_atoms(AtomSelectionEnum atom_selection, bool execute_sort_atoms) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Deactivate backbone hydrogens
     if (atom_selection & BACKBONE_H_ATOMS) {
          deactivate_atom(H,   false);
          deactivate_atom(H1, false);
          deactivate_atom(H2, false);
          deactivate_atom(H3, false);
     }

     // Deactivate non-backbone hydrogens
     if (atom_selection & NON_BACKBONE_H_ATOMS) {
          for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
               if ((*it) && ((*it)->atom_type>=HA && (*it)->atom_type<=HZ3)) {
                    deactivate_atom((*it)->atom_type, false);
               }
          }       
     }

     // Deactivate Sidechain atoms
     if (atom_selection & SIDECHAIN_ATOMS) {
          for (Residue::ReverseIterator it = atoms.rbegin(); it != atoms.rend(); ++it) {
               // CB is not disabled together with the other sidechain atoms
               if ((*it) && (*it)->is_sidechain_atom && (*it)->atom_type!=CB) {
                    deactivate_atom((*it)->atom_type, false);
               }
          }
          sidechain_status = false;
     }     
          
     // Deactivate O atoms
     if (atom_selection & BACKBONE_O_ATOMS) {
          deactivate_atom(O, false);
          deactivate_atom(OXT, false);
     }     

     // Deactivate CB atom
     if (atom_selection & CB_ATOMS) {
          deactivate_atom(CB, execute_sort_atoms);
     }


     // Sort atom vector
     if (execute_sort_atoms) {
          this->sort_atoms();
     }
}

// Add a Cbeta atom
void Residue::add_CB(AtomSelectionEnum atom_selection, bool execute_sort_atoms) {
     if (atom_selection & CB_ATOMS) {
          if (residue_type!=GLY && !has_atom(CB)) {
               atoms.push_back(new Atom(CB, this, atoms.size()));
          }
          init_bond_lengths();
          init_angles();

          if (execute_sort_atoms)
               this->sort_atoms();         
     }
}

// Add full-atom sidechain to residue
void Residue::add_sidechain(AtomSelectionEnum atom_selection, const std::vector<double> &chi_angles, bool execute_sort_atoms) {

     // Add CB if not already present
     if (residue_type!=GLY && !has_atom(CB)) {
          add_CB(atom_selection, false);
     }

     // Add sidechain atoms
     if (atom_selection & SIDECHAIN_ATOMS) {          
          double chi[4];
          for (unsigned int i=0; i<4; i++) {
               if (i<chi_angles.size())
                    chi[i] = chi_angles[i];
               else
                    chi[i] = UNINITIALIZED;
          }
     
          switch (residue_type) {
          case ALA:
               assert(chi_angles.size()==0);
               break;
          case CYS:
               if (!has_atom(SG)) {
                    atoms.push_back(new Atom(SG, this, atoms.size(), UNINITIALIZED, chi[0]));
               }
               break;
          case ASP:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(OD1))
                    atoms.push_back(new Atom(OD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(OD2))
                    atoms.push_back(new Atom(OD2, this, atoms.size()));
               break;
          case GLU:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD))
                    atoms.push_back(new Atom(CD, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(OE1))
                    atoms.push_back(new Atom(OE1, this, atoms.size(), UNINITIALIZED, chi[2]));
               if (!has_atom(OE2))
                    atoms.push_back(new Atom(OE2, this, atoms.size()));
               break;
          case PHE:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD1))
                    atoms.push_back(new Atom(CD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CD2))
                    atoms.push_back(new Atom(CD2, this, atoms.size()));
               if (!has_atom(CE1))
                    atoms.push_back(new Atom(CE1, this, atoms.size()));
               if (!has_atom(CE2))
                    atoms.push_back(new Atom(CE2, this, atoms.size()));
               if (!has_atom(CZ))
                    atoms.push_back(new Atom(CZ, this, atoms.size()));
               break;
          case GLY:
               break;
          case HIS:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(ND1))
                    atoms.push_back(new Atom(ND1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CD2))
                    atoms.push_back(new Atom(CD2, this, atoms.size()));
               if (!has_atom(CE1))
                    atoms.push_back(new Atom(CE1, this, atoms.size()));
               if (!has_atom(NE2))
                    atoms.push_back(new Atom(NE2, this, atoms.size()));
               break;
          case ILE:
               if (!has_atom(CG1))
                    atoms.push_back(new Atom(CG1, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CG2))
                    atoms.push_back(new Atom(CG2, this, atoms.size()));
               if (!has_atom(CD1))
                    atoms.push_back(new Atom(CD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               break;
          case LYS:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD))
                    atoms.push_back(new Atom(CD, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CE))
                    atoms.push_back(new Atom(CE, this, atoms.size(), UNINITIALIZED, chi[2]));
               if (!has_atom(NZ))
                    atoms.push_back(new Atom(NZ, this, atoms.size(), UNINITIALIZED, chi[3]));
               break;
          case LEU:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD1))
                    atoms.push_back(new Atom(CD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CD2))
                    atoms.push_back(new Atom(CD2, this, atoms.size()));
               break;
          case MET:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(SD))
                    atoms.push_back(new Atom(SD, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CE))
                    atoms.push_back(new Atom(CE, this, atoms.size(), UNINITIALIZED, chi[2]));
               break;
          case ASN:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(OD1))
                    atoms.push_back(new Atom(OD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(ND2))
                    atoms.push_back(new Atom(ND2, this, atoms.size()));
               break;
          case PRO:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD))
                    atoms.push_back(new Atom(CD, this, atoms.size(), UNINITIALIZED, chi[1]));
               break;
          case GLN:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD))
                    atoms.push_back(new Atom(CD, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(OE1))
                    atoms.push_back(new Atom(OE1, this, atoms.size(), UNINITIALIZED, chi[2]));
               if (!has_atom(NE2))
                    atoms.push_back(new Atom(NE2, this, atoms.size()));
               break;
          case ARG:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD))
                    atoms.push_back(new Atom(CD, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(NE))
                    atoms.push_back(new Atom(NE, this, atoms.size(), UNINITIALIZED, chi[2]));
               if (!has_atom(CZ))
                    atoms.push_back(new Atom(CZ, this, atoms.size(), UNINITIALIZED, chi[3]));
               if (!has_atom(NH1))
                    atoms.push_back(new Atom(NH1, this, atoms.size()));
               if (!has_atom(NH2))
                    atoms.push_back(new Atom(NH2, this, atoms.size()));
               break;
          case SER:
               if (!has_atom(OG))
                    atoms.push_back(new Atom(OG, this, atoms.size(), UNINITIALIZED, chi[0]));     
               break;
          case THR:
               if (!has_atom(OG1))
                    atoms.push_back(new Atom(OG1, this, atoms.size(), UNINITIALIZED, chi[0]));    
               if (!has_atom(CG2))
                    atoms.push_back(new Atom(CG2, this, atoms.size()));   
               break;
          case VAL:
               if (!has_atom(CG1))
                    atoms.push_back(new Atom(CG1, this, atoms.size(), UNINITIALIZED, chi[0]));    
               if (!has_atom(CG2))
                    atoms.push_back(new Atom(CG2, this, atoms.size()));   
               break;
          case TRP:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD1))
                    atoms.push_back(new Atom(CD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CD2))
                    atoms.push_back(new Atom(CD2, this, atoms.size()));
               if (!has_atom(NE1))
                    atoms.push_back(new Atom(NE1, this, atoms.size()));
               if (!has_atom(CE2))
                    atoms.push_back(new Atom(CE2, this, atoms.size()));
               if (!has_atom(CE3))
                    atoms.push_back(new Atom(CE3, this, atoms.size()));
               if (!has_atom(CZ2))
                    atoms.push_back(new Atom(CZ2, this, atoms.size()));
               if (!has_atom(CZ3))
                    atoms.push_back(new Atom(CZ3, this, atoms.size()));
               if (!has_atom(CH2))
                    atoms.push_back(new Atom(CH2, this, atoms.size()));
               break;
          case TYR:
               if (!has_atom(CG))
                    atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
               if (!has_atom(CD1))
                    atoms.push_back(new Atom(CD1, this, atoms.size(), UNINITIALIZED, chi[1]));
               if (!has_atom(CD2))
                    atoms.push_back(new Atom(CD2, this, atoms.size()));
               if (!has_atom(CE1))
                    atoms.push_back(new Atom(CE1, this, atoms.size()));
               if (!has_atom(CE2))
                    atoms.push_back(new Atom(CE2, this, atoms.size()));
               if (!has_atom(CZ))
                    atoms.push_back(new Atom(CZ, this, atoms.size()));
               if (!has_atom(OH))
                    atoms.push_back(new Atom(OH, this, atoms.size()));
               break;

          //Modified by MJ: Maybe also add P, O1P-O3P for SEP, TPO, PTR?
          case SEP:
              if (!has_atom(OG))
                  atoms.push_back(new Atom(OG, this, atoms.size(), UNINITIALIZED, chi[0]));
              break;

          case TPO:
              if (!has_atom(OG1))
                  atoms.push_back(new Atom(OG1, this, atoms.size(), UNINITIALIZED, chi[0]));
              if (!has_atom(CG2))
                  atoms.push_back(new Atom(CG2, this, atoms.size()));
              break;

          case PTR:

              if (!has_atom(CG))
                  atoms.push_back(new Atom(CG, this, atoms.size(), UNINITIALIZED, chi[0]));
              if (!has_atom(CD1))
                  atoms.push_back(new Atom(CD1, this, atoms.size(), UNINITIALIZED, chi[1]));
              if (!has_atom(CD2))
                  atoms.push_back(new Atom(CD2, this, atoms.size()));
              if (!has_atom(CE1))
                  atoms.push_back(new Atom(CE1, this, atoms.size()));
              if (!has_atom(CE2))
                  atoms.push_back(new Atom(CE2, this, atoms.size()));
              if (!has_atom(CZ))
                  atoms.push_back(new Atom(CZ, this, atoms.size()));
              if (!has_atom(OH))
                  atoms.push_back(new Atom(OH, this, atoms.size()));
              break;
          default:
               break;
          }
          init_bond_lengths();
          init_angles();
          sidechain_status = true;

     }

     if (execute_sort_atoms)
          this->sort_atoms();

}


// Add pseudo-sidechain atom to the residue
void Residue::add_pseudo_sidechain(AtomSelectionEnum atom_selection, const std::vector<double> &bond_length_and_angles, bool execute_sort_atoms) {
     // make sure not to add double atoms
     if (has_atom(PS))
          return;

     double bond_length = UNINITIALIZED;
     double angle = UNINITIALIZED;
     double dihedral = UNINITIALIZED;
     if (bond_length_and_angles.size() == 3) {
          bond_length = bond_length_and_angles[0];
          angle = bond_length_and_angles[1];
          dihedral = bond_length_and_angles[2];
     }

     // Create PS (pseudo-sidechain) atom
     Atom *atom_PS = new Atom(PS, this, atoms.size(), angle, dihedral);

     atom_PS->init_bond_length();
     if (is_initialized(bond_length))
          atom_PS->bond_length = bond_length;

     atoms.push_back(atom_PS);
     atom_PS->init_angles();

     // Sort atom vector
     if (execute_sort_atoms)
          this->sort_atoms();     
}

// calculates the position of the side chain center of mass
Vector_3D Residue::calc_sidechain_center_of_mass() {
     double x=0., y=0., z=0.;
     double ws=0;
     //
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          if ((*it)->is_sidechain_atom && (*it)->atom_type != CB && !is_atom_XH((*it)->atom_type) ) {
               double w=(*it)->get_mass();
               ws+=w;
               x+= w*(*it)->position[0];
               y+= w*(*it)->position[1];
               z+= w*(*it)->position[2];
          }
     }
     if (ws == 0)
          return Vector_3D(0,0,0);
     Vector_3D cm ((x/ws),(y/ws),(z/ws));
     return cm;
}

// calculates the position of the backbone center of mass
Vector_3D Residue::calc_backbone_center_of_mass(){
     double x=0., y=0., z=0.;
     double ws=0;
     //
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          if (!((*it)->is_sidechain_atom) && !is_atom_XH((*it)->atom_type) ) {
               double w=(*it)->get_mass();
               ws+=w;
               x+= w*(*it)->position[0];
               y+= w*(*it)->position[1];
               z+= w*(*it)->position[2];
          }
     }
     if (ws == 0)
          return Vector_3D(0,0,0);
     Vector_3D cm ((x/ws),(y/ws),(z/ws));
     return cm;
}


// calculates the position of the residue center of mass
Vector_3D Residue::calc_center_of_mass(){
     double x=0., y=0., z=0.;
     double ws=0;
     //
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          if (!is_atom_XH((*it)->atom_type)) {
               double w=(*it)->get_mass();
               ws+=w;
               x+= w*(*it)->position[0];
               y+= w*(*it)->position[1];
               z+= w*(*it)->position[2];
          }
     }
     if (ws == 0)
          return Vector_3D(0,0,0);
     Vector_3D cm ((x/ws),(y/ws),(z/ws));
     return cm;
}

// switching sidechain representation
void Residue::toggle_sidechain() {
     // lets try to find out what we need to do
     // do we have a sidechain here?
     if (this->has_atom(PS)) {
          // we have a Pseudo sidechain
          // .. lets go full atom
          this->toggle_sidechain_to_full_atom();

     } else if (this->has_sidechain()) {
          // there is a full atom sidechain
          // lets switch to PS
          this->toggle_sidechain_to_PS();

     } else {
          // per default .. switch to full atom
          this->toggle_sidechain_to_full_atom();
     }
}

// switching to pseudo sidechains
void Residue::toggle_sidechain_to_PS() {
     //
     // calculate the center of mass
     Vector_3D cm = this->calc_sidechain_center_of_mass();
     this->update_angles();
     //
     // deactivate the FullAtom Sidechain
     this->deactivate_atoms(SIDECHAIN_ATOMS, true);
     //
     // try to activate the PS atom .. otherwise create one.
     if (!this->activate_atom(PS, true)) {
          // lets create the atom first
          this->add_pseudo_sidechain(PSEUDO_SIDECHAIN_ATOMS);
          // lets update the position here already once,
          // so that everything get initialized properly
          this->atoms[this->atom_index[PS]]->update_position(0);
     }
     this->sort_atoms();
     //
     // now set the new position
     if (this->residue_type != ALA && this->residue_type != GLY) {
          // set it to its new position
          this->atoms[this->atom_index[PS]]->position = cm;
     }
     //
     // now update all the atoms since we removed
     // some atoms
     this->update_angles();
     // .. and set the proper bondlength
     if (this->residue_type != ALA && this->residue_type != GLY)
          this->atoms[this->atom_index[PS]]->init_bond_length();
     //
     // now update the positions again
     this->update_positions_non_backbone();
}

//
// switching to pseudo sidechains
void Residue::toggle_sidechain_to_full_atom() {

     this->activate_atoms(SIDECHAIN_ATOMS, true);
     // just to be sure all the atoms are in place
     this->add_sidechain(SIDECHAIN_ATOMS);
     //
     // then deactivate the pseudo sidechain as well.
     this->deactivate_atom(PS, true);
     //
     // and update the positions properly
     this->update_positions_backbone(1);
     this->update_positions_non_backbone();
}


// Return pointer to the first atom, given an iteration mode
Atom *Residue::get_first_atom(IterateEnum iteration_mode) {
     return atoms[iteration_range_indices[iteration_mode][BEGIN]];
}

// Return pointer to the last atom, given an iteration mode
Atom *Residue::get_last_atom(IterateEnum iteration_mode) {
     return atoms[iteration_range_indices[iteration_mode][END]];
}


// The chain is overspecified since it keeps track of both dihedral angles and positions.
// This method checks that the two representations are consistent
// The check is only done if (debuglevel > 0)
void Residue::check_consistency() {

#if DEBUGLEVEL > 0
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          (*it)->check_consistency();
     }
#endif
}

// Output to stream
void Residue::output(std::ostream &o) {
     o << "Residue[" << index << "]: " << residue_type << "\n";
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          o << "\t" << *it << "\n";
     }
}

// Output in PDB format
std::string Residue::output_as_pdb(int counterOffset, const int chain_number, std::string *b_factor_string) {

     std::map<AtomEnum,std::string> b_factor_map;
     if ((*b_factor_string) != "") {

          int comma_pos=0;
          do {
               comma_pos = b_factor_string->find_first_of(",");
               std::string substring = b_factor_string->substr(0, comma_pos);

               boost::trim(substring);

               if (substring.size() == 0)
                    break;

               int colon_pos = substring.find_first_of(":");
               std::string tag = substring.substr(0, colon_pos);
               std::string value = substring.substr(colon_pos+1, std::string::npos);

               std::string residue_string = tag;
               std::string atom_string = "*";
               size_t bracket_open_pos = tag.find_first_of("[");
               if (bracket_open_pos != std::string::npos) {
                    size_t bracket_close_pos = tag.find_first_of("]");
                    residue_string = tag.substr(0,bracket_open_pos);
                    atom_string = tag.substr(bracket_open_pos+1, (bracket_close_pos-bracket_open_pos));
               }

               int residue_index = boost::lexical_cast<int>(residue_string);
               if (residue_index == this->index)
                    *b_factor_string = b_factor_string->substr(comma_pos+1, std::string::npos);
               else
                    break;

               if (atom_string == "*")
                    b_factor_map.insert(std::make_pair(XX, value));
               else
                    b_factor_map.insert(std::make_pair(string_to_atom(atom_string), value));
          } while (comma_pos != (int)(std::string::npos));          
     }

     std::string output;
     // this->sort_atoms();
     for (Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          std::string b_factor_value = "";
          if (b_factor_string) {
               b_factor_value = map_lookup(b_factor_map, (*it)->atom_type, std::string(""));
               if (b_factor_value == "") {
                    b_factor_value = map_lookup(b_factor_map, XX, std::string(""));
               }
          }
          output += (*it)->output_as_pdb(counterOffset, chain_number, b_factor_value);
     }

     return output;
}

// Output
std::ostream & operator<<(std::ostream &o, Residue &r) {
     r.output(o);
     return o;
}



// Update backbone-atom positions
// The update procedure is slit in two (first backbone - then the rest), because some
// of the non backbone atoms are dependent on atoms outside the boundary of the current
// residue (for instance the H and C atoms)
void Residue::update_positions_backbone(int direction) {
     if (direction >=0) {
          Residue::Iterator begin = atoms.begin();
          // Residue::Iterator end = begin+this->backboneAtoms;
          Residue::Iterator end = begin+iteration_range_indices[BACKBONE][Residue::END]+1;
          for (Residue::Iterator it = begin; it != end; ++it) {
               (*it)->update_position(direction);
          }
     } else {
          // Residue::ReverseIterator begin = atoms.rbegin()+(atoms.size()-backboneAtoms);
          Residue::ReverseIterator begin = atoms.rbegin()+(atoms.size()-(iteration_range_indices[BACKBONE][Residue::END]+1));
          Residue::ReverseIterator end = atoms.rend();
          for (Residue::ReverseIterator it = begin; it != end; ++it) {
               (*it)->update_position(direction);
          }
     }
}


// Update non-backbone-atom positions
// The update procedure is split in two (first backbone - then the rest), because some
// of the non backbone atoms are dependent on atoms outside the boundary of the current
// residue (for instance the H and C atoms)
void Residue::update_positions_non_backbone() {
     // Residue::Iterator begin = atoms.begin()+backboneAtoms;
     Residue::Iterator begin = atoms.begin()+iteration_range_indices[BACKBONE][Residue::END]+1;
     Residue::Iterator end = atoms.end();
     for (Residue::Iterator it = begin; it != end; ++it) {
          (*it)->update_position(+1);
     }
}


// Set ideal bond lengths and bond angles, and update positions
void Residue::idealize(AtomSelectionEnum atom_selection) {

     // Get current degrees of freedom
     std::vector<double> dofs_backbone = get_angles();
     std::vector<double> dofs_sidechain = get_sidechain_dof_values(SIDECHAIN_ATOMS);
     std::vector<double> dofs_pseudosidechain = get_sidechain_dof_values(PSEUDO_SIDECHAIN_ATOMS);

     // Test which atoms should be idealized
     for (Residue::Iterator it = atoms.begin(); it != atoms.end(); ++it) {
          if (((atom_selection == ALL_ATOMS)) ||
              ((atom_selection & BACKBONE_ATOMS) && (*it)->is_backbone_atom) ||
              ((atom_selection & SIDECHAIN_ATOMS) && (*it)->is_sidechain_atom) ||
              ((atom_selection & PSEUDO_SIDECHAIN_ATOMS) && (*it)->atom_type==PS) ||
              ((atom_selection & CB_ATOMS) && (*it)->atom_type==CB) ||
              ((atom_selection & BACKBONE_O_ATOMS) && ((*it)->atom_type==O || (*it)->atom_type==OXT)) ||
              ((atom_selection & BACKBONE_H_ATOMS) && ((*it)->atom_type==H || (*it)->atom_type==H1 || (*it)->atom_type==H2 || (*it)->atom_type==H3)) ||
              ((atom_selection & NON_BACKBONE_H_ATOMS) && ((*it)->atom_type>=HA && (*it)->atom_type<=HZ3))) {

               (*it)->position.initialized = false;
               (*it)->set_angle(UNINITIALIZED);
               (*it)->set_dihedral(UNINITIALIZED);
          }
     }

     // Set degrees of freedom again
     set_angles(dofs_backbone);
     set_sidechain_dof_values(dofs_sidechain, SIDECHAIN_ATOMS);
     set_sidechain_dof_values(dofs_pseudosidechain, PSEUDO_SIDECHAIN_ATOMS);

     // Set ideal bond lengths
     init_bond_lengths();

     // Initialize angles
     init_angles();

}



// Set degree-of-freedom values of side chain.
// two modes: SIDECHAIN_ATOMS -> dofValueVector contains chi angles
//            PSEUDO_SIDECHAIN_ATOMS -> dofValueVector contains (bondlength, angle, dihedral)
void Residue::set_sidechain_dof_values(std::vector<double> &dof_value_vector,
                                    AtomSelectionEnum mode) {

     if (mode == PSEUDO_SIDECHAIN_ATOMS || !has_sidechain()) {
          if (this->has_atom(PS)) {

               // Do not update PS position for GLY and ALA
               if (residue_type==GLY || residue_type==ALA)
                    return;
               
               assert(dof_value_vector.size() == 3);
               double bond_length = dof_value_vector[0];
               double angle = dof_value_vector[1];
               double dihedral = dof_value_vector[2];
               
               Atom *atom_PS = (*this)[PS];
               atom_PS->bond_length = bond_length;
               atom_PS->set_angle(angle);
               atom_PS->set_dihedral(dihedral);
          }
     } else {
          assert(chi_atoms.size() >= dof_value_vector.size());
          
          // Set new chi values
          for (unsigned int i=0; i<dof_value_vector.size(); i++) {
               chi_atoms[i].first->get_dihedral() = dof_value_vector[i];
          }
     }
}



// Get degree-of-freedom values of side chain.
// two modes: SIDECHAIN_ATOMS -> dofValueVector contains chi angles
//            PSEUDO_SIDECHAIN_ATOMS -> dofValueVector contains (bondlength, angle, dihedral)
std::vector<double> Residue::get_sidechain_dof_values(AtomSelectionEnum mode) {

     std::vector<double> dof_value_vector;

     if (mode == PSEUDO_SIDECHAIN_ATOMS || (mode==ALL_ATOMS && !has_sidechain())) {

          if (this->has_atom(PS)) {
               dof_value_vector.resize(3);
               Atom *atom_PS = (*this)[PS];
               dof_value_vector[0] = atom_PS->bond_length;
               dof_value_vector[1] = atom_PS->get_angle();
               dof_value_vector[2] = atom_PS->get_dihedral();
          }       
     } else {
     
          for (unsigned int i=0; i<chi_atoms.size(); i++) {
               dof_value_vector.push_back(chi_atoms[i].first->get_dihedral());
          }
     }
          
     return dof_value_vector;
}

// Set values of bondangles in residue which can be modified freely 
// without modifying the backbone
void Residue::set_minor_dof_values(const std::vector<double> &dof_value_vector) {
     for (unsigned int i=0; i<minor_dof_atoms.size(); i++) {
          AngleEnum angle_type = minor_dof_atoms[i].second;
          if (angle_type == ANGLE)
               minor_dof_atoms[i].first->get_angle() = dof_value_vector[i];
          else
               minor_dof_atoms[i].first->get_dihedral() = dof_value_vector[i];
     }
}

// Get values of bondangles in residue which can be modified freely 
// without modifying the backbone
std::vector<double> Residue::get_minor_dof_values() {
     std::vector<double> dof_value_vector;
     for (unsigned int i=0; i<minor_dof_atoms.size(); ++i) {
          AngleEnum angle_type = minor_dof_atoms[i].second;
          if (angle_type == ANGLE)          
               dof_value_vector.push_back(minor_dof_atoms[i].first->get_angle());
          else
               dof_value_vector.push_back(minor_dof_atoms[i].first->get_dihedral());
     }
     return dof_value_vector;
}

// A residue keeps track of pointers to the atoms corresponding to
// its sidechain chi angles for easy reference. This method sets up such a pointer
void Residue::set_chi_atom(Atom *atom, int index) {

     if (index < 0) {
          chi_atoms.push_back(std::make_pair(atom,DIHEDRAL));
     } else {
          if ((int)chi_atoms.size() <= index)
               chi_atoms.resize(index+1);
          chi_atoms[index] = std::make_pair(atom,DIHEDRAL);
     }
}

// A residue keeps track of pointers to the atoms corresponding to
// the atoms which have bond angle dofs that can be modified without
// effecting the backbone
void Residue::set_minor_dof_atom(Atom *atom, AngleEnum atom_type, int index) {
     if (index < 0) {
          minor_dof_atoms.push_back(std::make_pair(atom,atom_type));
     } else {
          if ((int)minor_dof_atoms.size() <= index)
               minor_dof_atoms.resize(index+1);
          minor_dof_atoms[index] = std::make_pair(atom,atom_type);
     }
}



/////////////////////////////////////
//// ResidueCA  - Calpha backbone////
/////////////////////////////////////

// Constructor - incomplete initialization
ResidueCA::ResidueCA(ResidueEnum residue_type,
                     Chain<ResidueCA> *chain, 
                     int index,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = (ChainCA*)chain;
}

// Constructor - incomplete initialization - angles must be initalized later
ResidueCA::ResidueCA(ResidueEnum residue_type, ChainCA *chain, int index,
                     Vector_3D position_CA,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {
     this->chain = chain;
     init();
     this->atoms.push_back(new Atom(AtomEnum(CA), this, 0, position_CA));
     init_bond_lengths();
}

// Constructor - incomplete initialization - angles must be initalized later
ResidueCA::ResidueCA(ResidueEnum residue_type, ChainCA *chain, int index,
                     std::vector<std::pair<definitions::AtomEnum, Vector_3D> > &positions,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = chain;
     init();

     for (unsigned int i=0; i<positions.size(); i++) {
          this->atoms.push_back(new Atom(positions[i].first, this, i, positions[i].second));
     }
     init_bond_lengths();
     this->sort_atoms();     
}

// Constructor - incomplete initialization - atom positions must be initalized later
ResidueCA::ResidueCA(ResidueEnum residue_type,
                     ChainCA *chain, 
                     int index, 
                     double theta, double tau,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = chain;
     init();
     this->atoms.push_back(new Atom(AtomEnum(CA), this, 0, theta, tau));
     init_bond_lengths();
}

// Copy constructor
ResidueCA::ResidueCA(ResidueCA &r, ChainCA *chain, 
                     int index,
                     int index_res_seq)
     : Residue(r, index, index_res_seq) {

     this->chain = chain;
     // init();
}

// Initializer
void ResidueCA::init() {
     // backboneAtoms = 1;
     // init_bond_lengths();

     terminal_status = INTERNAL;
     if (index==0)
          terminal_status = NTERM;
     else if (index==(chain->size()-1))
          terminal_status = CTERM;

     this->sort_atoms();

     atoms_length = size();

}

// Destructor
ResidueCA::~ResidueCA(){}


// Add a selection of atoms to the residue
void ResidueCA::add_atoms(AtomSelectionEnum atom_selection, const std::vector<double> &sidechain_dof_values, bool execute_sort_atoms) {

     // Only add atoms to non-terminal residues (positioning not well-defined at terminals)
     if (terminal_status == INTERNAL) {
          add_CB(atom_selection, false);

          if (atom_selection & PSEUDO_SIDECHAIN_ATOMS) {
               add_pseudo_sidechain(atom_selection, sidechain_dof_values, false);
          }
     
          if (execute_sort_atoms)
               this->sort_atoms();
     }
}


// Set new theta value
void ResidueCA::set_theta(double theta) {
     (*this)[CA]->set_angle(theta);
}

// Return theta value
double ResidueCA::get_theta() {
     return (*this)[CA]->get_angle();
}

// Set new tau value
void ResidueCA::set_tau(double tau) {
     (*this)[CA]->set_dihedral(tau);
}

// Return tau value
double ResidueCA::get_tau() {
     return (*this)[CA]->get_dihedral();
}

// Return vector of theta, tau angles
// The includeBondAngles option has no effect here
std::vector<double> ResidueCA::get_angles(bool include_omega, bool include_bond_angles) {
     std::vector<double> angles;
     angles.push_back(get_theta());
     angles.push_back(get_tau());
     return angles;
}

// Set new theta, tau pair
void ResidueCA::set_angles(std::vector<double> angles) {
     set_theta(angles[0]);
     set_tau(angles[1]);
}


// Get neighbouring residue
ResidueCA *ResidueCA::get_neighbour(int offset) {
     int new_index = index + offset;

     if (new_index >= 0 && new_index < chain->size()) {
          return &(*chain)[new_index];
     } else {
          return NULL;
     }
}


// Output to stream
void ResidueCA::output(std::ostream &o) {
     Residue::output(o);
     o << "\t(Theta, Tau) = (" << std::setw(7) << get_theta() << ", " << std::setw(7) << get_tau() << ")";
}

// Output
std::ostream & operator<<(std::ostream &o, ResidueCA &r) {
     r.output(o);
     return o;
}




///////////////////////////////////////////
//// ResidueFB - Full backbone residue ////
///////////////////////////////////////////

// Constructor - incomplete initialization - angles must be initalized later
ResidueFB::ResidueFB(ResidueEnum residue_type,
                     Chain<ResidueFB> *chain, 
                     int index,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = (ChainFB*)chain;
     init();
}

// Constructor - incomplete initialization - angles must be initalized later
ResidueFB::ResidueFB(ResidueEnum residue_type,
                     ChainFB *chain, 
                     int index, 
                     Vector_3D position_N,
                     Vector_3D position_CA, 
                     Vector_3D position_C,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = chain;
     this->atoms.push_back(new Atom(AtomEnum(N), this, 0, position_N));
     this->atoms.push_back(new Atom(AtomEnum(CA), this, 1, position_CA));
     this->atoms.push_back(new Atom(AtomEnum(C), this, 2, position_C));
     init_bond_lengths();
}

// Constructor - incomplete initialization - angles must be initalized later
ResidueFB::ResidueFB(ResidueEnum residue_type,
                     ChainFB *chain, 
                     int index, 
                     std::vector<std::pair<definitions::AtomEnum, Vector_3D> > &positions,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = chain;
     init();

     // Move all hydrogens to the end of the array - the hydrogens assume that the backbone atoms are in place
     std::vector<std::pair<AtomEnum, Vector_3D> > hydrogen_positions;
     for (unsigned int i=0; i<positions.size();) {
          if (is_atom_XH(positions[i].first)) {
               hydrogen_positions.push_back(positions[i]);
               positions.erase(positions.begin()+i);
          } else {
               ++i;
          }
     }     

     // Add hydrogens to the end of original vector
     positions.insert(positions.end(), hydrogen_positions.begin(), hydrogen_positions.end() );

     for (unsigned int i=0; i<positions.size(); i++) {
          this->atoms.push_back(new Atom(positions[i].first, this, i, positions[i].second));
     }
     this->sort_atoms();
     init_bond_lengths();
     
}

// Constructor - incomplete initialization - atom positions must be initalized later
ResidueFB::ResidueFB(ResidueEnum residue_type,
                     ChainFB *chain, 
                     int index, 
                     double phi, double psi, double omega,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     this->chain = chain;
     init();

     double angle_N = bond_angle_constants(C, N, CA, residue_type, NULL, terminal_status);
     double angle_CA = bond_angle_constants(N, CA, C, residue_type, NULL, terminal_status);
     double angle_C = bond_angle_constants(CA, C, N, residue_type, NULL, terminal_status);

     this->atoms.push_back(new Atom(AtomEnum(N), this, 0, angle_N, omega));
     this->atoms.push_back(new Atom(AtomEnum(CA), this, 1, angle_CA, phi));
     this->atoms.push_back(new Atom(AtomEnum(C), this, 2, angle_C, psi));
     init_bond_lengths();
}

// Constructor - incomplete initialization - atom positions must be initalized later
ResidueFB::ResidueFB(ResidueEnum residue_type,
                     ChainFB *chain, 
                     int index, 
                     const std::vector<double> &angles,
                     int index_res_seq)
     : Residue(residue_type, index, index_res_seq) {

     assert(angles.size() > 1);

     this->chain = chain;
     init();
     
     double phi = angles[0];
     double psi = angles[1];

     double omega = M_PI;

     double angle_N = bond_angle_constants(C, N, CA, residue_type, NULL, terminal_status);
     double angle_CA = bond_angle_constants(N, CA, C, residue_type, NULL, terminal_status);
     double angle_C = bond_angle_constants(CA, C, N, residue_type, NULL, terminal_status);

     // Check if there are additional angles in the vector
     if (angles.size() > 2)
          omega = angles[2];
     else if (angles.size() > 3)
          angle_N = angles[3];
     else if (angles.size() > 4)
          angle_CA = angles[4];
     else if (angles.size() > 5)
          angle_C = angles[5];

     this->atoms.push_back(new Atom(AtomEnum(N), this, 0, angle_N, omega));
     this->atoms.push_back(new Atom(AtomEnum(CA), this, 1, angle_CA, phi));
     this->atoms.push_back(new Atom(AtomEnum(C), this, 2, angle_C, psi));
     init_bond_lengths();
}

// Destructor
ResidueFB::~ResidueFB() {}

// Copy constructor
ResidueFB::ResidueFB(ResidueFB &r, 
                     ChainFB *chain, 
                     int index, 
                     int index_res_seq)
     : Residue(r, index, index_res_seq) {

     this->chain = chain;
     // init();
}

// Initializer
void ResidueFB::init() {
     terminal_status = INTERNAL;
     if (index==0)
          terminal_status = NTERM;
     else if (index==(chain->size()-1))
          terminal_status = CTERM;

     this->sort_atoms();

     atoms_length = size();
}

// Add backbone oxygen atom(s) to residue
void ResidueFB::add_oxygens(AtomSelectionEnum atom_selection, bool execute_sort_atoms) {

     if (atom_selection & BACKBONE_O_ATOMS) {

          if (terminal_status == CTERM) {
               if (!has_atom(O)) {
                    Atom *atomO = new Atom(O, this, atoms.size());
                    atoms.push_back(atomO);
               }
               if (!has_atom(OXT)) {
                    Atom *atomOXT = new Atom(OXT, this, atoms.size());
                    atoms.push_back(atomOXT);
               }
          } else {
               if (!has_atom(O)) {
                    Atom *atomO = new Atom(O, this, atoms.size());
                    atoms.push_back(atomO);
               }
          }

          init_bond_lengths();
          init_angles();

          if (execute_sort_atoms)
               this->sort_atoms();         
     }
}


// Add hydrogen atom(s) to residue
// atom_selection determines whether to add backbone H-atoms, sidechain H-atoms, or both
void ResidueFB::add_hydrogens(AtomSelectionEnum atom_selection, bool execute_sort_atoms) {

     // Backbone H-atoms
     if (atom_selection & BACKBONE_H_ATOMS) {

          if (terminal_status == NTERM) {

               if (!has_atom(H1)) {
                    Atom *atom1H = new Atom(H1, this, atoms.size());
                    atoms.push_back(atom1H);
               }

               if (!has_atom(H2)) {
                    Atom *atom2H = new Atom(H2, this, atoms.size());
                    atoms.push_back(atom2H);
               }

               if (residue_type!=PRO) {
                    if (!has_atom(H3)) {
                         Atom *atom3H = new Atom(H3, this, atoms.size());
                         atoms.push_back(atom3H);
                    }
               }
          } else {
               if (residue_type!=PRO) {   
                    if (!has_atom(H)) {
                         Atom *atomH = new Atom(H, this, atoms.size());
                         atoms.push_back(atomH);
                    }
               }
          }
     }

     // Sidechain H-atoms
     if (atom_selection & NON_BACKBONE_H_ATOMS) {

          if (!has_atom(HA) && residue_type!=GLY)
               atoms.push_back(new Atom(HA, this, atoms.size()));

          if (has_sidechain()) {          
               // Residue specific side chain hydrogens
               switch (residue_type) {
               case ALA:
                    if (!has_atom(HB1))
                         atoms.push_back(new Atom(HB1, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));
                    break;
               case CYS:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));
                    if (!has_atom(HG))
                         atoms.push_back(new Atom(HG, this, atoms.size()));
                    break;
               case ASP:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    break;
               case GLU:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG2))
                         atoms.push_back(new Atom(HG2, this, atoms.size()));
                    if (!has_atom(HG3))
                         atoms.push_back(new Atom(HG3, this, atoms.size()));           
                    break;
               case PHE:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HD1))
                         atoms.push_back(new Atom(HD1, this, atoms.size()));
                    if (!has_atom(HD2))
                         atoms.push_back(new Atom(HD2, this, atoms.size()));           
                    if (!has_atom(HE1))
                         atoms.push_back(new Atom(HE1, this, atoms.size()));
                    if (!has_atom(HE2))
                         atoms.push_back(new Atom(HE2, this, atoms.size()));           
                    if (!has_atom(HZ))
                         atoms.push_back(new Atom(HZ, this, atoms.size()));
                    break;
               case GLY:
                    if (!has_atom(HA2))
                         atoms.push_back(new Atom(HA2, this, atoms.size()));
                    if (!has_atom(HA3))
                         atoms.push_back(new Atom(HA3, this, atoms.size()));
                    break;
               case HIS:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));
                    if (!has_atom(HD2))
                         atoms.push_back(new Atom(HD2, this, atoms.size()));           
                    if (!has_atom(HE1))
                         atoms.push_back(new Atom(HE1, this, atoms.size()));
                    // At pH 7, HIS has either an HE2 or an HD1 
                    if (!has_atom(HE2) && !has_atom(HD1)) {
                         // atoms.push_back(new Atom(HE2, this, atoms.size()));
                         atoms.push_back(new Atom(HD1, this, atoms.size()));
                    }
                    break;
               case ILE:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB))
                         atoms.push_back(new Atom(HB, this, atoms.size()));
                    if (!has_atom(HG12))
                         atoms.push_back(new Atom(HG12, this, atoms.size()));
                    if (!has_atom(HG13))
                         atoms.push_back(new Atom(HG13, this, atoms.size()));
                    if (!has_atom(HG21))
                         atoms.push_back(new Atom(HG21, this, atoms.size()));
                    if (!has_atom(HG22))
                         atoms.push_back(new Atom(HG22, this, atoms.size()));
                    if (!has_atom(HG23))
                         atoms.push_back(new Atom(HG23, this, atoms.size()));
                    if (!has_atom(HD11))
                         atoms.push_back(new Atom(HD11, this, atoms.size()));
                    if (!has_atom(HD12))
                         atoms.push_back(new Atom(HD12, this, atoms.size()));
                    if (!has_atom(HD13))
                         atoms.push_back(new Atom(HD13, this, atoms.size()));
                    break;
               case LYS:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG2))
                         atoms.push_back(new Atom(HG2, this, atoms.size()));
                    if (!has_atom(HG3))
                         atoms.push_back(new Atom(HG3, this, atoms.size()));           
                    if (!has_atom(HD2))
                         atoms.push_back(new Atom(HD2, this, atoms.size()));
                    if (!has_atom(HD3))
                         atoms.push_back(new Atom(HD3, this, atoms.size()));           
                    if (!has_atom(HE2))
                         atoms.push_back(new Atom(HE2, this, atoms.size()));
                    if (!has_atom(HE3))
                         atoms.push_back(new Atom(HE3, this, atoms.size()));           
                    if (!has_atom(HZ1))
                         atoms.push_back(new Atom(HZ1, this, atoms.size()));           
                    if (!has_atom(HZ2))
                         atoms.push_back(new Atom(HZ2, this, atoms.size()));
                    if (!has_atom(HZ3))
                         atoms.push_back(new Atom(HZ3, this, atoms.size()));           
                    break;
               case LEU:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG))
                         atoms.push_back(new Atom(HG, this, atoms.size()));
                    if (!has_atom(HD11))
                         atoms.push_back(new Atom(HD11, this, atoms.size()));
                    if (!has_atom(HD12))
                         atoms.push_back(new Atom(HD12, this, atoms.size()));
                    if (!has_atom(HD13))
                         atoms.push_back(new Atom(HD13, this, atoms.size()));
                    if (!has_atom(HD21))
                         atoms.push_back(new Atom(HD21, this, atoms.size()));
                    if (!has_atom(HD22))
                         atoms.push_back(new Atom(HD22, this, atoms.size()));
                    if (!has_atom(HD23))
                         atoms.push_back(new Atom(HD23, this, atoms.size()));
                    break;
               case MET:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG2))
                         atoms.push_back(new Atom(HG2, this, atoms.size()));
                    if (!has_atom(HG3))
                         atoms.push_back(new Atom(HG3, this, atoms.size()));           
                    if (!has_atom(HE1))
                         atoms.push_back(new Atom(HE1, this, atoms.size()));
                    if (!has_atom(HE2))
                         atoms.push_back(new Atom(HE2, this, atoms.size()));
                    if (!has_atom(HE3))
                         atoms.push_back(new Atom(HE3, this, atoms.size()));           
                    break;
               case ASN:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HD21))
                         atoms.push_back(new Atom(HD21, this, atoms.size()));
                    if (!has_atom(HD22))
                         atoms.push_back(new Atom(HD22, this, atoms.size()));
                    break;
               case PRO:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG2))
                         atoms.push_back(new Atom(HG2, this, atoms.size()));
                    if (!has_atom(HG3))
                         atoms.push_back(new Atom(HG3, this, atoms.size()));           
                    if (!has_atom(HD2))
                         atoms.push_back(new Atom(HD2, this, atoms.size()));
                    if (!has_atom(HD3))
                         atoms.push_back(new Atom(HD3, this, atoms.size()));           
                    break;
               case GLN:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG2))
                         atoms.push_back(new Atom(HG2, this, atoms.size()));
                    if (!has_atom(HG3))
                         atoms.push_back(new Atom(HG3, this, atoms.size()));           
                    if (!has_atom(HE21))
                         atoms.push_back(new Atom(HE21, this, atoms.size()));          
                    if (!has_atom(HE22))
                         atoms.push_back(new Atom(HE22, this, atoms.size()));          
                    break;
               case ARG:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG2))
                         atoms.push_back(new Atom(HG2, this, atoms.size()));
                    if (!has_atom(HG3))
                         atoms.push_back(new Atom(HG3, this, atoms.size()));           
                    if (!has_atom(HD2))
                         atoms.push_back(new Atom(HD2, this, atoms.size()));
                    if (!has_atom(HD3))
                         atoms.push_back(new Atom(HD3, this, atoms.size()));           
                    if (!has_atom(HE))
                         atoms.push_back(new Atom(HE, this, atoms.size()));            
                    if (!has_atom(HH11))
                         atoms.push_back(new Atom(HH11, this, atoms.size()));          
                    if (!has_atom(HH12))
                         atoms.push_back(new Atom(HH12, this, atoms.size()));          
                    if (!has_atom(HH21))
                         atoms.push_back(new Atom(HH21, this, atoms.size()));          
                    if (!has_atom(HH22))
                         atoms.push_back(new Atom(HH22, this, atoms.size()));          
                    break;
               case SER:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HG))
                         atoms.push_back(new Atom(HG, this, atoms.size()));            
                    break;
               case THR:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB))
                         atoms.push_back(new Atom(HB, this, atoms.size()));
                    if (!has_atom(HG1))
                         atoms.push_back(new Atom(HG1, this, atoms.size()));
                    if (!has_atom(HG21))
                         atoms.push_back(new Atom(HG21, this, atoms.size()));
                    if (!has_atom(HG22))
                         atoms.push_back(new Atom(HG22, this, atoms.size()));
                    if (!has_atom(HG23))
                         atoms.push_back(new Atom(HG23, this, atoms.size()));
                    break;
               case VAL:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB))
                         atoms.push_back(new Atom(HB, this, atoms.size()));
                    if (!has_atom(HG11))
                         atoms.push_back(new Atom(HG11, this, atoms.size()));
                    if (!has_atom(HG12))
                         atoms.push_back(new Atom(HG12, this, atoms.size()));
                    if (!has_atom(HG13))
                         atoms.push_back(new Atom(HG13, this, atoms.size()));
                    if (!has_atom(HG21))
                         atoms.push_back(new Atom(HG21, this, atoms.size()));
                    if (!has_atom(HG22))
                         atoms.push_back(new Atom(HG22, this, atoms.size()));
                    if (!has_atom(HG23))
                         atoms.push_back(new Atom(HG23, this, atoms.size()));
                    break;
               case TRP:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HD1))
                         atoms.push_back(new Atom(HD1, this, atoms.size()));
                    if (!has_atom(HE1))
                         atoms.push_back(new Atom(HE1, this, atoms.size()));
                    if (!has_atom(HE3))
                         atoms.push_back(new Atom(HE3, this, atoms.size()));
                    if (!has_atom(HZ2))
                         atoms.push_back(new Atom(HZ2, this, atoms.size()));
                    if (!has_atom(HZ3))
                         atoms.push_back(new Atom(HZ3, this, atoms.size()));           
                    if (!has_atom(HH2))
                         atoms.push_back(new Atom(HH2, this, atoms.size()));           
                    break;
               case TYR:
                    if (!has_atom(HA))
                         atoms.push_back(new Atom(HA, this, atoms.size()));
                    if (!has_atom(HB2))
                         atoms.push_back(new Atom(HB2, this, atoms.size()));
                    if (!has_atom(HB3))
                         atoms.push_back(new Atom(HB3, this, atoms.size()));           
                    if (!has_atom(HD1))
                         atoms.push_back(new Atom(HD1, this, atoms.size()));
                    if (!has_atom(HD2))
                         atoms.push_back(new Atom(HD2, this, atoms.size()));
                    if (!has_atom(HE1))
                         atoms.push_back(new Atom(HE1, this, atoms.size()));
                    if (!has_atom(HE2))
                         atoms.push_back(new Atom(HE2, this, atoms.size()));
                    if (!has_atom(HH))
                         atoms.push_back(new Atom(HH, this, atoms.size()));            
                    break;

               //Modified by MJ> Maybe add more error checking
               case SEP:
                   if (!has_atom(HA))
                       atoms.push_back(new Atom(HA, this, atoms.size()));
                   if (!has_atom(HB2))
                       atoms.push_back(new Atom(HB2, this, atoms.size()));
                   if (!has_atom(HB3))
                       atoms.push_back(new Atom(HB3, this, atoms.size()));
                   break;
               case TPO:
                   if (!has_atom(HA))
                       atoms.push_back(new Atom(HA, this, atoms.size()));
                   if (!has_atom(HB))
                       atoms.push_back(new Atom(HB, this, atoms.size()));
                   if (!has_atom(HG21))
                       atoms.push_back(new Atom(HG21, this, atoms.size()));
                   if (!has_atom(HG22))
                       atoms.push_back(new Atom(HG22, this, atoms.size()));
                   if (!has_atom(HG23))
                       atoms.push_back(new Atom(HG23, this, atoms.size()));
                   break;
               case PTR:
                   if (!has_atom(HA))
                       atoms.push_back(new Atom(HA, this, atoms.size()));
                   if (!has_atom(HB2))
                       atoms.push_back(new Atom(HB2, this, atoms.size()));
                   if (!has_atom(HB3))
                       atoms.push_back(new Atom(HB3, this, atoms.size()));
                   if (!has_atom(HD1))
                       atoms.push_back(new Atom(HD1, this, atoms.size()));
                   if (!has_atom(HD2))
                       atoms.push_back(new Atom(HD2, this, atoms.size()));
                   if (!has_atom(HE1))
                       atoms.push_back(new Atom(HE1, this, atoms.size()));
                   if (!has_atom(HE2))
                       atoms.push_back(new Atom(HE2, this, atoms.size()));
                   break;
               default:
                    break;
               }
          }
     }
     
     init_bond_lengths();
     init_angles();

     if (execute_sort_atoms)
          this->sort_atoms();      
}

// Add a selection of atoms to the residue
void ResidueFB::add_atoms(AtomSelectionEnum atom_selection, const std::vector<double> &dofValues, bool execute_sort_atoms) {
     add_oxygens(atom_selection, false);
     add_CB(atom_selection, false);
     
     if (atom_selection & SIDECHAIN_ATOMS) {
          add_sidechain(atom_selection, dofValues, false);
     } else if (atom_selection & PSEUDO_SIDECHAIN_ATOMS) {
          add_pseudo_sidechain(atom_selection, dofValues, false);
     }
     
     add_hydrogens(atom_selection, false);

     if (execute_sort_atoms)
          this->sort_atoms();
}



// Set new phi dihedral
void ResidueFB::set_phi(double phi) {
     (*this)[CA]->set_dihedral(phi);
}

// Return phi dihedral
double &ResidueFB::get_phi() {
     return (*this)[CA]->get_dihedral();
}

// Set new psi dihedral
void ResidueFB::set_psi(double psi) {
     (*this)[C]->set_dihedral(psi);
}

// Return psi dihedral
double &ResidueFB::get_psi() {
     return (*this)[C]->get_dihedral();
}

// Set new omega dihedral
void ResidueFB::set_omega(double omega) {
     (*this)[N]->set_dihedral(omega);
}

// Return omega dihedral
double &ResidueFB::get_omega() {
     return (*this)[N]->get_dihedral();
}


// Set new angles in residue
void ResidueFB::set_angles(std::vector<double> angles) {
     set_phi(angles[0]);
     set_psi(angles[1]);

     if (angles.size() > 2)
          set_omega(angles[2]);
     else
          return;

     if (angles.size() > 3)
          // C-N-CA angle
          (*this)[N]->set_angle(angles[3]);
     else
          return;

     if (angles.size() > 4)
          // N-CA-C angle
          (*this)[CA]->set_angle(angles[4]);
     else
          return;

     if (angles.size() > 5)
          // CA-C-N angle
          (*this)[C]->set_angle(angles[5]);
     else
          return;
}

// Return current psi psi dihedrals as a vector (plus optional bond angles)
std::vector<double> ResidueFB::get_angles(bool include_omega, bool include_bond_angles) {
     std::vector<double> angles;
     angles.push_back(get_phi());
     angles.push_back(get_psi());
     if(include_omega)
             angles.push_back(get_omega());

     if (include_bond_angles) {
          angles.push_back((*this)[N]->get_angle());
          angles.push_back((*this)[CA]->get_angle());
          angles.push_back((*this)[C]->get_angle());
     }
     return angles;
}

// Get neighbouring residue
ResidueFB *ResidueFB::get_neighbour(int offset) {
     int new_index = index + offset;

     if (new_index >= 0 && new_index < chain->size())
          return &(*chain)[new_index];
     else {
          return NULL;
     }
}

// Returns residue name - potentially dependent on residue 
// protonation state
std::string ResidueFB::get_protonation_dependent_name() {

     if (residue_type == HIS) {
          if (has_atom(HD1) && has_atom(HE2))
               return "HISH";
          else if (has_atom(HD1))
               return "HISD";
          else if (has_atom(HE2))
               return "HISE";
          else
               return "HIS";
     }
     else 
          return boost::lexical_cast<std::string>(residue_type);
}

// Output to stream
void ResidueFB::output(std::ostream &o) {
     Residue::output(o);
     o << "\t(Phi, Psi) = (" << std::setw(7) << get_phi() << ", " << std::setw(7) << get_psi() << ")\n";
     if (has_sidechain())
         o << "\t(chi) = " ;
     else
         o << "\t(d,Omega,Phi) = " ;
     o << this->get_sidechain_dof_values() << "";
}

// Output
std::ostream & operator<<(std::ostream &o, ResidueFB &r) {
     r.output(o);
     return o;
}

}
