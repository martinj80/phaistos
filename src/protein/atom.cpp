// atom.cpp --- Atom class
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


#include "atom.h"
#include "chain.h"
#include "residue.h"
#include "protein/protein_ideal_values.h"
#include "definitions.h"

namespace phaistos {

// Import protein definitions (such as residue names)
using namespace definitions;

// Initializer
void Atom::init(AtomEnum atom_type, Residue *residue, int index) {
     this->atom_type = atom_type;
     this->residue = residue;
     this->index = index;
     this->residue->atom_index[atom_type] = index;

     this->biotype = -1;
     this->bond_length = UNINITIALIZED;
     this->mass = UNINITIALIZED;
     
     // Create dihedral value pointer
     if (!this->dihedral)
          this->dihedral = new double(UNINITIALIZED);
     
     // Create angle value pointer
     if (!this->angle)
          this->angle = new double(UNINITIALIZED);

     // Usually, angle and dihedral pointers are owned by the atom,
     // which is therefore responsible for cleaning it up. In the N
     // terminus, the H-atom dihedral points to the already existing
     // phi-dihedral, which should therefore not be deleted on the
     // deletion of the H-atom. Likewise for the O-atom in the C-terminus
     this->owner_of_dihedral = true;
     this->owner_of_angle = true;
     
     // If this is the N-terminus, let phi determine position of H
     if (atom_type == H1 || atom_type == H) {
          if (residue->terminal_status == NTERM) {
               if (dihedral && dihedral != &(((ResidueFB*)residue)->get_phi()))
                    delete dihedral;
               dihedral = &(((ResidueFB*)residue)->get_phi());
               this->owner_of_dihedral = false;
          }
     }

     // If this is the C-terminus, let last psi determine position of O
     if (atom_type==O) {
          if (residue->terminal_status == CTERM) {
               if (dihedral && dihedral != &(((ResidueFB*)residue)->get_psi()))
                    delete dihedral;
               dihedral = &(((ResidueFB*)residue)->get_psi());
               this->owner_of_dihedral = false;
          }
     }

     // Determine whether atom is a sidechain or backbone atom
     is_backbone_atom = is_backbone_atom_type(atom_type);
     is_sidechain_atom = is_sidechain_atom_type(atom_type);
     
     // Non-backbone atom keeps track of neighbouring angles themselves
     // see header file for details
     self_maintained_positioning = false;
     if (!is_backbone_atom) {
          self_maintained_positioning = true;
     }

     // Update #backboneAtoms counter in residue
     if (is_backbone_atom) {
          // if ((index+1) > residue->backboneAtoms) {
          //      residue->backboneAtoms = index+1;
          // }
          if (index > residue->iteration_range_indices[BACKBONE][Residue::END]) {
               residue->iteration_range_indices[BACKBONE][Residue::END] = index;
          }
     }
     
     // Cache for quick navigation of immediate neighbours
     for (int i=0; i<definitions::ITERATE_ENUM_SIZE; i++) {
          neighbour_atom[0][i] = NULL;
          neighbour_atom[1][i] = NULL;
          neighbour_atom[2][i] = NULL;
          neighbour_atom[3][i] = NULL;
          neighbour_atom[4][i] = NULL;
          neighbour_atom[5][i] = NULL;
          neighbour_atom[6][i] = NULL;
     }


     covalent_neighbours.clear();
     switch (atom_type) {
     case N:
          if (residue->terminal_status != NTERM)
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(C, -1));
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
          if (residue->terminal_status == NTERM) {
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(H1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(H2, 0));
               if (residue->residue_type!=PRO)
                    covalent_neighbours.push_back(std::pair<AtomEnum, int>(H3, 0));
               else
                    covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
          } else {
               if (residue->residue_type==PRO)
                    covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
          }
          // The H is always added, even in NTERM residues - currently, in Phaistos
          // the beginning of a chain is marked NTERM even if there are no H1 and H2
          // atoms present (this might change in the future)
          if (residue->residue_type!=PRO)
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(H,  0));
          this->mass = definitions::atom_n_weight;
          break;
     case CA:
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(N,  0));
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(C,  0));
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(HA, 0));
          if (residue->residue_type==GLY) {
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HA2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HA3, 0));
          }
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(PS,  0));
          this->mass = definitions::atom_c_weight;
          break;
     case C:
          if (residue->terminal_status != CTERM)
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(N,  +1));
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(O,  0));
          if (residue->terminal_status == CTERM)
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OXT,  0));
          this->mass = definitions::atom_c_weight;
          break;
     case PS:
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
          this->mass = 1;
          break;
     case O:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(C,  0));
          this->mass = definitions::atom_o_weight;
          break;
     case OXT:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(C,  0));
          this->mass = definitions::atom_o_weight;
          break;
     case H:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(N,  0));
          this->mass = definitions::atom_h_weight;
          break;
     case H1:
          residue->set_minor_dof_atom(this, ANGLE);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(N,  0));
          this->mass = definitions::atom_h_weight;
          break;
     case H2:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(N,  0));
          this->mass = definitions::atom_h_weight;
          break;
     case H3:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(N,  0));
          this->mass = definitions::atom_h_weight;
          break;
     case HA:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
          this->mass = definitions::atom_h_weight;
          break;
     case HA2:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
          this->mass = definitions::atom_h_weight;
          break;
     case HA3:
          residue->set_minor_dof_atom(this, ANGLE);
          residue->set_minor_dof_atom(this, DIHEDRAL);
          covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
          this->mass = definitions::atom_h_weight;
          break;
     default:
          break;
     }
     
     // Residue specific initialization
     switch (residue->residue_type) {
     case ALA:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB1:
               // Chi 1 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          default:
               break;
          }
          break;
     case CYS:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(SG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case SG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG, 0));
               this->mass = definitions::atom_s_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG:
               // Chi 2 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(SG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          default:
               break;
          }
          break;
     case ASP:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1 
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case OD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD1, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case OD2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD1:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OD1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               // Chi 3 (Hydrogen degree-of-freedom) - assumes that HD1 is not also present
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OD2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case GLU:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OE1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case OE1:
               // Chi 3
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case OE2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HE1:
               // Chi 4 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OE1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE2:
               // Chi 4 (Hydrogen degree-of-freedom) - assumes that HE1 is not also present
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OE2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case PHE:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD1, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE1:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE2:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CZ:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HZ, 0));
               this->mass = definitions::atom_c_weight;
               break;          
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HZ:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case GLY:
          break;
     case HIS:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(ND1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case ND1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD1, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case CD2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));         
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE1:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(ND1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE2,  0));        
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case NE2:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1,  0));        
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2,  0));
               this->mass = definitions::atom_n_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(ND1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case ILE:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG1:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG12, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG13,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG21, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG22, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG23, 0));
               this->mass = definitions::atom_c_weight;
               break;          
          case CD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1,   0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD11, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD12, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD13, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG12:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG13:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG21:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG23:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD11:
               // Chi 4 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD12:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD13:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case LYS:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE:
               // Chi 3
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NZ,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case NZ:
               // Chi 4
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HZ1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HZ2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HZ3, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HE3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HZ1:
               // Chi 5 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 4);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NZ,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HZ2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NZ,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HZ3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NZ,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case LEU:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD11, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD12, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD13, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD21, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD22, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD23, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD11:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD12:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD13:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD21:
               // Chi 4 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD23:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case MET:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(SD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case SD:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               this->mass = definitions::atom_s_weight;
               break;
          case CE:
               // Chi 3
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(SD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1, 0));                 
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2, 0));                 
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE1:
               // Chi 4 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HE2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HE3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case ASN:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OD1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(ND2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case OD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case ND2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD21, 0));        
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD22, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD21:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(ND2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(ND2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case PRO:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(N, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case GLN:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OE1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE2,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case OE1:
               // Chi 3
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case NE2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE21, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE22, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE21:
               // Chi 4 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HE22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case ARG:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD3, 0));                 
               this->mass = definitions::atom_c_weight;
               break;
          case NE:
               // Chi 3
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case CZ:
               // Chi 4
               residue->set_chi_atom(this, 3);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE, 0));          
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NH1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NH2,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case NH1:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH11, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH12, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case NH2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH21, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH22, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HD3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HH11:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NH1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HH12:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NH1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HH21:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NH2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HH22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NH2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case SER:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case OG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG,  0));
               this->mass = definitions::atom_o_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG:
               // Chi 2 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;

     //added by MJ
     case SEP:
         switch (atom_type) {
         case CB:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
             this->mass = definitions::atom_c_weight;
             break;
         case OG:
             // Chi 1
             residue->set_chi_atom(this, 0);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;
         
         case P:
             // Chi 2
             residue->set_chi_atom(this, 1);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(O1P, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(O2P, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(O3P, 0));
             this->mass = definitions::atom_p_weight;
             break;

         case O1P:
             //// Chi 3
             residue->set_chi_atom(this, 2);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);
             /*residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);*/

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;

         case O2P:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;

         case O3P:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;

         case HB2:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
             this->mass = definitions::atom_h_weight;
             break;

         case HB3:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
             this->mass = definitions::atom_h_weight;
             break;

         default:
             break;

         }
         break;


     case THR:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case OG1:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG1,  0));
               this->mass = definitions::atom_o_weight;
               break;
          case CG2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG21, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG22, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG23, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB:
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               this->mass = definitions::atom_h_weight;
               break;          
          case HG1:
               // Chi 2 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG21:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG23:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;

     case TPO:
         switch (atom_type) {
         case CB:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG1, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB, 0));
             this->mass = definitions::atom_c_weight;
             break;
         case OG1:
             // Chi 1
             residue->set_chi_atom(this, 0);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG1, 0));
             this->mass = definitions::atom_o_weight;
             break;
         case CG2:
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG21, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG22, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG23, 0));
             this->mass = definitions::atom_c_weight;
             break;
         case HB:
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             this->mass = definitions::atom_h_weight;
             break;

         case P:
             // Chi 2
             residue->set_chi_atom(this, 1);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(OG1, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(O1P, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(O2P, 0));
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(O3P, 0));
             this->mass = definitions::atom_p_weight;
             break;

         case O1P:
             //// Chi 4
             residue->set_chi_atom(this, 3);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);
             /*residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);*/

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;

         case O2P:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;

         case O3P:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
             this->mass = definitions::atom_o_weight;
             break;

         case HG21:
             // Chi 3 (Hydrogen degree-of-freedom)
             residue->set_chi_atom(this, 2);
             residue->sidechain_status = true;
             residue->set_minor_dof_atom(this, ANGLE);

             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
             this->mass = definitions::atom_h_weight;
             break;
         case HG22:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
             this->mass = definitions::atom_h_weight;
             break;
         case HG23:
             residue->set_minor_dof_atom(this, ANGLE);
             residue->set_minor_dof_atom(this, DIHEDRAL);
             covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
             this->mass = definitions::atom_h_weight;
             break;
         default:
             break;
         }
         break;

     case VAL:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB,  0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG1:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG11, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG12, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG13, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG21, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG22, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HG23, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG11:
               // Chi 2 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG12:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG13:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HG21:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG22:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HG23:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case TRP:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD1, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case NE1:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1, 0));
               this->mass = definitions::atom_n_weight;
               break;
          case CE2:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE3:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ3, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CZ2:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CH2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HZ2, 0));         
               this->mass = definitions::atom_c_weight;
               break;
          case CZ3:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE3, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CH2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HZ3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CH2:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ3, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(NE1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE3,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HZ2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HZ3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ3,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HH2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CH2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
     case TYR:
          switch (atom_type) {
          case CB:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CG:
               // Chi 1
               residue->set_chi_atom(this, 0);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD1:
               // Chi 2
               residue->set_chi_atom(this, 1);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD1, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CD2:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE1:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CE2:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case CZ:
               residue->sidechain_status=true;

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OH, 0));                                          
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1,  0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
               this->mass = definitions::atom_c_weight;
               break;
          case OH:
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH, 0));
               this->mass = definitions::atom_o_weight;
               break;
          case HB2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;
          case HB3:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HD2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE1:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HE2:
               residue->set_minor_dof_atom(this, ANGLE);
               residue->set_minor_dof_atom(this, DIHEDRAL);
               covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          case HH:
               // Chi 3 (Hydrogen degree-of-freedom)
               residue->set_chi_atom(this, 2);
               residue->sidechain_status=true;
               residue->set_minor_dof_atom(this, ANGLE);

               covalent_neighbours.push_back(std::pair<AtomEnum, int>(OH,  0));
               this->mass = definitions::atom_h_weight;
               break;          
          default:
               break;
          }
          break;
          //Added by MJ:
     case PTR:
          switch (atom_type) {
              case CB:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CA, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB2, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HB3, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case CG:
                  // Chi 1
                  residue->set_chi_atom(this, 0);
                  residue->sidechain_status = true;
                  residue->set_minor_dof_atom(this, ANGLE);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case CD1:
                  // Chi 2
                  residue->set_chi_atom(this, 1);
                  residue->sidechain_status = true;
                  residue->set_minor_dof_atom(this, ANGLE);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD1, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case CD2:
                  residue->sidechain_status = true;
                  residue->set_minor_dof_atom(this, DIHEDRAL);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CG, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HD2, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case CE1:
                  residue->sidechain_status = true;

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE1, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case CE2:
                  residue->sidechain_status = true;

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HE2, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case CZ:
                  residue->sidechain_status = true;

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(OH, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
                  this->mass = definitions::atom_c_weight;
                  break;
              case OH:
                  residue->sidechain_status = true;
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CZ, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(HH, 0));
                  this->mass = definitions::atom_o_weight;
                  break;
              case HB2:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
                  this->mass = definitions::atom_h_weight;
                  break;
              case HB3:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CB, 0));
                  this->mass = definitions::atom_h_weight;
                  break;
              case HD1:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD1, 0));
                  this->mass = definitions::atom_h_weight;
                  break;
              case HD2:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CD2, 0));
                  this->mass = definitions::atom_h_weight;
                  break;
              case HE1:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE1, 0));
                  this->mass = definitions::atom_h_weight;
                  break;
              case HE2:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(CE2, 0));
                  this->mass = definitions::atom_h_weight;
                  break;

              case P:
                  // Chi 3
                  residue->set_chi_atom(this, 2);
                  residue->sidechain_status = true;
                  residue->set_minor_dof_atom(this, ANGLE);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(OH, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(O1P, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(O2P, 0));
                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(O3P, 0));
                  this->mass = definitions::atom_p_weight;
                  break;

              case O1P:
                  // Chi 4
                  residue->set_chi_atom(this, 3);
                  residue->sidechain_status = true;
                  residue->set_minor_dof_atom(this, ANGLE);
                  /*residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);*/

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
                  this->mass = definitions::atom_o_weight;
                  break;

              case O2P:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
                  this->mass = definitions::atom_o_weight;
                  break;

              case O3P:
                  residue->set_minor_dof_atom(this, ANGLE);
                  residue->set_minor_dof_atom(this, DIHEDRAL);

                  covalent_neighbours.push_back(std::pair<AtomEnum, int>(P, 0));
                  this->mass = definitions::atom_o_weight;
                  break;

              default:
                  break;
              }
              break;

     default:
          break;
     }
}

// returns the mass of the atom
double Atom::get_mass() {
     return this->mass;
}

// Determine the 3 atoms defining current angle
void Atom::get_angle_atoms(Atom **atom1, Atom **atom2, Atom **atom3) {
     if (this->self_maintained_positioning) {
          *atom1 = get_neighbour(-2, POSITIONING);
          *atom2 = get_neighbour(-1, POSITIONING);
          *atom3 = this;
     } else {
          *atom1 = get_neighbour(-1, POSITIONING);
          *atom2 = this;
          *atom3 = get_neighbour(+1, POSITIONING);
     }
}

// Determine the 4 atoms defining current dihedral
void Atom::get_dihedral_atoms(Atom **atom1, Atom **atom2,
                              Atom **atom3, Atom **atom4) {
     if (this->self_maintained_positioning) {
          *atom1 = get_neighbour(-3, POSITIONING);
          *atom2 = get_neighbour(-2, POSITIONING);
          *atom3 = get_neighbour(-1, POSITIONING);
          *atom4 = this;
     } else {
          *atom1 = get_neighbour(-2, POSITIONING);
          *atom2 = get_neighbour(-1, POSITIONING);
          *atom3 = this;
          *atom4 = get_neighbour(+1, POSITIONING);
     }
}

// Initialize bond lengths either from neighbouring positions (if
// initialized), or from idealized values
void Atom::init_bond_length() {
     Atom *prev = get_neighbour(-1, BACKBONE);
     if (prev) {
          // Use positions to determine bond length if positions are available
          if (prev->position.initialized && this->position.initialized) {
               this->bond_length = (this->position - prev->position).norm();
          // otherwise, use default bond lengths
          } else {
               this->bond_length = bond_length_constants(prev->atom_type, this->atom_type, prev->residue->residue_type);
          }
     }
}


// Set angle and dihedral values that have not been initialized
void Atom::init_angles() {

     // Check if angle has been initialized
     if (!is_initialized(get_angle())) {
          
          Atom *atom1, *atom2, *atom3;
          get_angle_atoms(&atom1, &atom2, &atom3);

          if (atom1 && atom2 && atom3) {
               this->set_angle(bond_angle_constants(atom1->atom_type,
                                                    atom2->atom_type,
                                                    atom3->atom_type,
                                                    residue->residue_type,
                                                    NULL,
                                                    residue->terminal_status));
          }
          
          // if (this->self_maintained_positioning) {
          //      Atom *prevPrev = get_neighbour(-2, POSITIONING);
          //      Atom *prev     = get_neighbour(-1, POSITIONING);
          //      Atom *current  = this;
          //      if (prevPrev && prev)
          //           this->set_angle(bond_angle_constants(prevPrev->atom_type,
          //                                          prev->atom_type,
          //                                          current->atom_type,
          //                                          residue->residue_type,
          //                                          NULL,
          //                                          residue->terminal_status));
          // } else {
          //      Atom *prev    = get_neighbour(-1, POSITIONING);
          //      Atom *current = this;
          //      Atom *next    = get_neighbour(+1, POSITIONING);
          //      if (prev && next)
          //           this->set_angle(bond_angle_constants(prev->atom_type,
          //                                          current->atom_type,
          //                                          next->atom_type,
          //                                          residue->residue_type,
          //                                          NULL,
          //                                          residue->terminal_status));
          // }
     }

     // Check if dihedral has been initialized
     if (!is_initialized(get_dihedral())) {

          Atom *atom1, *atom2, *atom3, *atom4;
          get_dihedral_atoms(&atom1, &atom2, &atom3, &atom4);

          if (atom1 && atom2 && atom3 && atom4) {
               this->set_dihedral(dihedral_constants(atom1->atom_type,
                                                   atom2->atom_type,
                                                   atom3->atom_type,
                                                   atom4->atom_type,
                                                   residue->residue_type));
          }
          
          // // Check if atom keeps track of prev dihedral itself
          // if (this->self_maintained_positioning) {
          //      Atom *prevPrevPrev = get_neighbour(-3, POSITIONING);
          //      Atom *prevPrev     = get_neighbour(-2, POSITIONING);
          //      Atom *prev         = get_neighbour(-1, POSITIONING);
          //      Atom *current      = this;
          //      if (prevPrevPrev && prevPrev && prev)
          //           this->set_dihedral(dihedral_constants(prevPrevPrev->atom_type,
          //                                            prevPrev->atom_type,
          //                                            prev->atom_type,
          //                                            current->atom_type,
          //                                            residue->residue_type));

          // } else {
          //      Atom *prevPrev = get_neighbour(-2, POSITIONING);
          //      Atom *prev     = get_neighbour(-1, POSITIONING);
          //      Atom *current  = this;
          //      Atom *next     = get_neighbour(+1, POSITIONING);
          //      if (prevPrev && prev && next)
          //           this->set_dihedral(dihedral_constants(prevPrev->atom_type,
          //                                            prev->atom_type,
          //                                            current->atom_type,
          //                                            next->atom_type,
          //                                            residue->residue_type));
          // }
     }
}


// Constructor - incomplete initialization
Atom::Atom(AtomEnum atomType, Residue *residue, int index):angle(NULL), dihedral(NULL) {
     init(atomType, residue, index);
}


// Constructor - incomplete initialization - angles must be initalized later
Atom::Atom(AtomEnum atomType, Residue *residue, int index, Vector_3D position)
     : angle(NULL), dihedral(NULL), position(position) {
     init(atomType, residue, index);
}

// Constructor - incomplete initialization - position must be initalized later
Atom::Atom(AtomEnum atomType, Residue *residue, int index,
           double angle, double dihedral):angle(NULL), dihedral(NULL) {
     init(atomType, residue, index);

     set_angle(angle);
     set_dihedral(dihedral);
}

// Copy constructor
Atom::Atom(Atom &other):angle(NULL), dihedral(NULL) {

     init(other.atom_type, other.residue, other.index);

     this->bond_length = other.bond_length;
     this->position = other.position;
     set_angle(other.get_angle());
     set_dihedral(other.get_dihedral());

     this->biotype = other.biotype;
     
     this->is_sidechain_atom = other.is_sidechain_atom;
     this->is_backbone_atom = other.is_backbone_atom;
}

// Copy constructor
Atom::Atom(Atom &other, Residue *residue):angle(NULL), dihedral(NULL) {

     init(other.atom_type, residue, other.index);

     this->bond_length = other.bond_length;
     this->position = other.position;
     set_angle(other.get_angle());
     set_dihedral(other.get_dihedral());

     this->biotype = other.biotype;
     
     this->is_sidechain_atom = other.is_sidechain_atom;
     this->is_backbone_atom = other.is_backbone_atom;
}

// Destructor
Atom::~Atom() {
     if (this->angle && this->owner_of_angle)
          delete this->angle;

     if (this->dihedral && this->owner_of_dihedral)
          delete this->dihedral;
}

// atom-neighbour constants. 
// This method is slow, but the values will be cached
// by the get_neighbour method
Atom *Atom::get_neighbour_constants(int &offset, definitions::IterateEnum iteration_type) {

     if (offset==0) {
          return this;
     }

     if (iteration_type == CA_ONLY) {
          // Ensure that CA iteration starts at a CA atom
          if (atom_type != CA) {
               if (offset>0 && (atom_type==N || atom_type==H)) {
                    offset -= 1;
               } else if ((offset<0) && (!(atom_type==N || atom_type==H))) {
                    offset += 1;
               }
               return (*residue)[CA];
          }

     } else if (iteration_type == SC_ONLY) {

          // Ensure that SC_ONLY iteration starts at a sidechain atom
          if (!this->is_sidechain_atom) {
               if (offset>0) {
                    // For atoms before sidechain, go to CB - this counts as a step
                    // if (atom_type==N || atom_type==CA || atom_type==H || atom_type==H1 || atom_type== H2 || atom_type== H3) {
                    if (atom_type==N || atom_type==CA) {
                         offset -= 1;
                         return residue->atoms[residue->iteration_range_indices[SC_ONLY][Residue::BEGIN]];
                    // Backbone hydrogens should simply be skipped 
                    } else if (atom_type==H || atom_type==H1 || atom_type== H2 || atom_type== H3) {
                         int nextOffset = +1;
                         return get_neighbour_constants(nextOffset, ALL);                         
                    // For atoms after sidechain, go to last sidechain atom, the first step
                    // made hereafter will take you to the sidechain of the next residue
                    } else {
                         return residue->atoms[residue->iteration_range_indices[SC_ONLY][Residue::END]-1];
                    }
               } else {
                    // For atoms before the sidechain, go to CB, the first step
                    // made hereafter will take you to the sidechain of the previous residue
                    if (atom_type==N || atom_type==CA || atom_type==H || atom_type==H1 || atom_type== H2 || atom_type== H3) {
                         return residue->atoms[residue->iteration_range_indices[SC_ONLY][Residue::BEGIN]];
                    // For backbone atoms after sidechain, go to the last sidechain atom
                    // this counts as a step
                    } else if (atom_type==C || atom_type==O || atom_type==OXT) {
                         offset += 1;
                         return residue->atoms[residue->iteration_range_indices[SC_ONLY][Residue::END]-1];
                    // For sidechain hydrogens, go to nearest sidechain atom, this counts as a step 
                    } else if (atom_type>=HA && atom_type<=HZ3) {
                         int prev_offset = -1;
                         return get_neighbour_constants(prev_offset, BACKBONE);                    
                    }
               }
          }       
     }


     // The code below is only relevant for POSITIONING and BACKBONE iterators
     if (iteration_type != POSITIONING && iteration_type != BACKBONE) {
          return NULL;
     }
     
     // Residue independent constants
     switch (atom_type) {
     case PS:           // PS
     case CB:           // CB
          if (iteration_type==POSITIONING) {
               switch (offset) {
               case -1: offset=+1; return (*residue)[CA]; 
               case -2: offset=-1; return (*residue)[CA]; 
               case -3: offset=0; return (*residue)[CA]; 
               default: offset=0; return NULL; 
               }
          } else {
               if (offset<0) offset+=1;
               return (*residue)[CA];
               break;
          }
          break;
     case O: {          // O
          
          Residue *next_res = residue->get_neighbour(+1);

          if (iteration_type==POSITIONING) {
               if (next_res) {
                    switch (offset) {
                    case -1: offset=+1; return (*residue)[C]; 
                    case -2: offset=-1; return (*residue)[C]; 
                    case -3: offset=0; return (*residue)[C]; 
                    default: offset=0; return NULL; 
                    }
               } else { // C-terminal
                    switch (offset) {
                    case -1: offset=0; return (*residue)[C]; 
                    case -2: offset=-1; return (*residue)[C]; 
                    case -3: offset=-2; return (*residue)[C]; 
                    default: offset=0; return NULL; 
                    }
               }
          } else {
               if (offset < 0) offset += 1;
               return (*residue)[C];
          }
          break;
     }
     case OXT:                 // C-terminal
          if (iteration_type==POSITIONING) {
               switch (offset) {
               case -1: offset=0; return (*residue)[C]; 
               case -2: offset=-1; return (*residue)[C]; 
               case -3: offset=0; return (*residue)[O]; 
               default: offset=0; return NULL;
               }
          } else {
               if (offset < 0) offset += 1;
               return (*residue)[C];
          }
          break;
     case H:            // H
     case H1: {         // H - N-terminal
          Residue *prevRes = residue->get_neighbour(-1);
          if (iteration_type==POSITIONING) {
               if (prevRes) {
                    switch (offset) {
                    case -1: offset=+1; return (*residue)[N]; 
                    case -2: offset=-1; return (*residue)[N]; 
                    case -3: offset=0; return (*residue)[N];
                    default: offset=0; return NULL;
                    }
               } else {         // N-terminal
                    switch (offset) {
                    case -1: offset=0; return (*residue)[N]; 
                    case -2: offset=+1; return (*residue)[N]; 
                    case -3: offset=+2; return (*residue)[N];
                    default: offset=0; return NULL;
                    }
               }
          } else {
               if (offset < 0) offset += 1;
               return (*residue)[N];
          }       
          break;
     }
     case H2:
     case H3:
          if (iteration_type==POSITIONING) {
               switch (offset) {
               case -1: offset=0; return (*residue)[N]; 
               case -2: offset=0; return (*residue)[CA]; 
               case -3: offset=0; return (*residue)[H1];
               default: offset=0; return NULL;
               }
          } else {
               if (offset < 0) offset += 1;
               return (*residue)[N];
          }
          break;
     case HA:           // HA
     case HA2:
     case HA3:
          if (iteration_type==POSITIONING) {
               switch (offset) {
               case -1: offset=+1; return (*residue)[CA]; 
               case -2: offset=-1; return (*residue)[CA]; 
               case -3: offset=0; return (*residue)[CA];
               default: offset=0; return NULL;
               }
          } else {
               if (offset < 0) offset += 1;
               return (*residue)[CA];
          }
          break;
     default: break;
     }


     // Residue dependent constants
     switch (residue->residue_type) {
     case ALA:
          switch (atom_type) {
          case HB1:
               if (iteration_type==POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA];
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type==POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[HB1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }
               }
               break;
          default: return NULL;
          }
          break;
     case CYS:
          switch (atom_type) {
          case SG:
               if (iteration_type==POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[SG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB];
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[SG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[SG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case ASP:
          switch (atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case OD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA];
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case OD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[OD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OD1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OD1];
                    case -2: offset=0; return (*residue)[CG];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OD2];
                    case -2: offset=0; return (*residue)[CG];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          default: return NULL; 
          }
          break;
     case GLU:
          switch (atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case OE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case OE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[OE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG2:
          case HG3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OE1]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OE1];
                    case -2: offset=0; return (*residue)[CD];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OE2]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OE2];
                    case -2: offset=0; return (*residue)[CD];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          default: return NULL; 
          }
          break;
     case PHE:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    // case -3: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CZ:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CG];               
                    case -4: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1];
                    case -2: offset=0; return (*residue)[CD1];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CE2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2];
                    case -2: offset=0; return (*residue)[CD2];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HZ:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2]; 
                    case -2: offset=0; return (*residue)[CE1]; 
                    case -3: offset=0; return (*residue)[CZ]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ];
                    case -2: offset=0; return (*residue)[CE1];
                    case -3: offset=0; return (*residue)[CD1];
                    case -4: offset=0; return (*residue)[CG];
                    case -5: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case HIS:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case ND1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[ND1]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case NE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[ND1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[ND1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2]; 
                    case -2: offset=0; return (*residue)[ND1]; 
                    case -3: offset=0; return (*residue)[CE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1];
                    case -2: offset=0; return (*residue)[ND1];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[NE2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2];
                    case -2: offset=0; return (*residue)[CD2];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case ILE:
          switch(atom_type) {
          case CG1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CG2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CG1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }                               
               }
               break;
          case CD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG12:
          case HG13:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG22:
          case HG23:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[HG21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HD11:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG1]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG1];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD12:
          case HD13:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG1]; 
                    case -3: offset=0; return (*residue)[HD11]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG1];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case LYS:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case NZ:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE];
                    case -2: offset=0; return (*residue)[CD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG2:
          case HG3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD2:
          case HD3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE2:
          case HE3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NZ]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CE]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE];
                    case -2: offset=0; return (*residue)[CD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HZ1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NZ]; 
                    case -2: offset=0; return (*residue)[CE]; 
                    case -3: offset=0; return (*residue)[CD]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NZ];
                    case -2: offset=0; return (*residue)[CE];
                    case -3: offset=0; return (*residue)[CD];
                    case -4: offset=0; return (*residue)[CG];
                    case -5: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HZ2:
          case HZ3:            
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NZ]; 
                    case -2: offset=0; return (*residue)[CE]; 
                    case -3: offset=0; return (*residue)[HZ1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NZ];
                    case -2: offset=0; return (*residue)[CE];
                    case -3: offset=0; return (*residue)[CD];
                    case -4: offset=0; return (*residue)[CG];
                    case -5: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }
               }
               break;
          default: return NULL; 
          }
          break;
     case LEU:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD11:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD12:
          case HD13:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[HD11]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD22:
          case HD23:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[HD21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case MET:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case SD:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[SD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[SD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG2:
          case HG3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[SD]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE]; 
                    case -2: offset=0; return (*residue)[SD]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE];
                    case -2: offset=0; return (*residue)[SD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE2:
          case HE3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE]; 
                    case -2: offset=0; return (*residue)[SD]; 
                    case -3: offset=0; return (*residue)[HE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE];
                    case -2: offset=0; return (*residue)[SD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case ASN:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case OD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case ND2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[OD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD22:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[HD21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[ND2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case PRO:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case CD:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG2:
          case HG3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD2:
          case HD3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=-1; return (*residue)[CA]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case GLN:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case OE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case NE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[OE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG2:
          case HG3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2];
                    case -2: offset=0; return (*residue)[CD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HE22:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[HE21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE2];
                    case -2: offset=0; return (*residue)[CD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;          
          default: return NULL; 
          }
          break;
     case ARG:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case NE:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case CZ:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[CG];               
                    case -4: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }
               }                    
               break;
          case NH1:
          case NH2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[NE]; 
                    case -3: offset=0; return (*residue)[CD]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[NE]; 
                    case -3: offset=0; return (*residue)[CD]; 
                    case -4: offset=0; return (*residue)[CG];               
                    case -5: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HG2:
          case HG3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD2:
          case HD3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CD]; 
                    case -3: offset=0; return (*residue)[NE]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE];
                    case -2: offset=0; return (*residue)[CD];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HH11:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH1]; 
                    case -2: offset=0; return (*residue)[CZ]; 
                    case -3: offset=0; return (*residue)[NE]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH1];
                    case -2: offset=0; return (*residue)[CZ];
                    case -3: offset=0; return (*residue)[NE];
                    case -4: offset=0; return (*residue)[CD];
                    case -5: offset=0; return (*residue)[CG];
                    case -6: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=7;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HH12:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH1]; 
                    case -2: offset=0; return (*residue)[CZ]; 
                    case -3: offset=0; return (*residue)[HH11]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH1];
                    case -2: offset=0; return (*residue)[CZ];
                    case -3: offset=0; return (*residue)[NE];
                    case -4: offset=0; return (*residue)[CD];
                    case -5: offset=0; return (*residue)[CG];
                    case -6: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=7;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HH21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH2]; 
                    case -2: offset=0; return (*residue)[CZ]; 
                    case -3: offset=0; return (*residue)[NE]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH2];
                    case -2: offset=0; return (*residue)[CZ];
                    case -3: offset=0; return (*residue)[NE];
                    case -4: offset=0; return (*residue)[CD];
                    case -5: offset=0; return (*residue)[CG];
                    case -6: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=7;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HH22:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH2]; 
                    case -2: offset=0; return (*residue)[CZ]; 
                    case -3: offset=0; return (*residue)[HH21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NH2];
                    case -2: offset=0; return (*residue)[CZ];
                    case -3: offset=0; return (*residue)[NE];
                    case -4: offset=0; return (*residue)[CD];
                    case -5: offset=0; return (*residue)[CG];
                    case -6: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=7;
                         return (*residue)[CA];
                    }               
               }
               break;          
          default: return NULL; 
          }
          break;
     case SER:
          switch(atom_type) {
          case OG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          default: return NULL; 
          }
          break;               
     case THR:
          switch (atom_type) {
          case OG1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CG2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[OG1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OG1]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OG1]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OG1];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG22:
          case HG23:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[HG21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;
     case VAL:
          switch(atom_type) {
          case CG1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CG2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CG1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HB:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HG11:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HG12:
          case HG13:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[HG11]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG1];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HG21:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HG22:
          case HG23:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[HG21]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG2];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;          
          default: return NULL; 
          }
          break;
     case TRP:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case CD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }                    
               break;
          case NE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CE2]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CZ2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CE3]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CG];               
                    case -4: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }
               }                                    
               break;
          case CZ3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE3]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CE2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE3]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CG];               
                    case -4: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }
               }                    
               break;
          case CH2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ3]; 
                    case -2: offset=0; return (*residue)[CE3]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ3]; 
                    case -2: offset=0; return (*residue)[CE3]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    case -4: offset=0; return (*residue)[CG];               
                    case -5: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[NE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[NE1];
                    case -2: offset=0; return (*residue)[CD1];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HZ2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CH2]; 
                    case -2: offset=0; return (*residue)[CE2]; 
                    case -3: offset=0; return (*residue)[CZ2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ2];
                    case -2: offset=0; return (*residue)[CE2];
                    case -3: offset=0; return (*residue)[CD2];
                    case -4: offset=0; return (*residue)[CG];
                    case -5: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }
               }
               break;          
          case HH2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ3]; 
                    case -2: offset=0; return (*residue)[CZ2]; 
                    case -3: offset=0; return (*residue)[CH2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CH2];
                    case -2: offset=0; return (*residue)[CZ2];
                    case -3: offset=0; return (*residue)[CE2];
                    case -4: offset=0; return (*residue)[CD2];
                    case -5: offset=0; return (*residue)[CG];
                    case -6: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=7;
                         return (*residue)[CA];
                    }               
               }
               break;          
          case HE3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ3]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CE3]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE3];
                    case -2: offset=0; return (*residue)[CD2];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;                  
          case HZ3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CH2]; 
                    case -2: offset=0; return (*residue)[CE3]; 
                    case -3: offset=0; return (*residue)[CZ3]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ3];
                    case -2: offset=0; return (*residue)[CE3];
                    case -3: offset=0; return (*residue)[CD2];
                    case -4: offset=0; return (*residue)[CG];
                    case -5: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }
               }
               break;                  
          default: return NULL; 
          }
          break;
     case TYR:
          switch(atom_type) {
          case CG:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=-1; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CB]; 
                    case -3: offset=0; return (*residue)[CA]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG];
                    case -2: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=3;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case CZ:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CG]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CG];               
                    case -4: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }
               }
               break;
          case OH:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CE1]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CE1]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    case -4: offset=0; return (*residue)[CG];               
                    case -5: offset=0; return (*residue)[CB];               
                    default:
                         if (offset < 0) offset+=6;
                         return (*residue)[CA];
                    }
               }
               break;
          case HB2:
          case HB3:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CG]; 
                    case -2: offset=0; return (*residue)[CA]; 
                    case -3: offset=0; return (*residue)[CB]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=2;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD1];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HD2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2]; 
                    case -2: offset=0; return (*residue)[CG]; 
                    case -3: offset=0; return (*residue)[CD2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CD2];
                    case -2: offset=0; return (*residue)[CG];
                    case -3: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=4;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE1:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CD1]; 
                    case -3: offset=0; return (*residue)[CE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE1];
                    case -2: offset=0; return (*residue)[CD1];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HE2:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CZ]; 
                    case -2: offset=0; return (*residue)[CD2]; 
                    case -3: offset=0; return (*residue)[CE2]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[CE2];
                    case -2: offset=0; return (*residue)[CD2];
                    case -3: offset=0; return (*residue)[CG];
                    case -4: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=5;
                         return (*residue)[CA];
                    }               
               }
               break;
          case HH:
               if (iteration_type == POSITIONING) {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OH]; 
                    case -2: offset=0; return (*residue)[CZ]; 
                    case -3: offset=0; return (*residue)[CE1]; 
                    default: offset=0; return NULL;
                    }
               } else {
                    switch (offset) {
                    case -1: offset=0; return (*residue)[OH];
                    case -2: offset=0; return (*residue)[CZ];
                    case -3: offset=0; return (*residue)[CE1];
                    case -4: offset=0; return (*residue)[CD1];
                    case -5: offset=0; return (*residue)[CG];
                    case -6: offset=0; return (*residue)[CB];
                    default:
                         if (offset < 0) offset+=7;
                         return (*residue)[CA];
                    }               
               }
               break;
          default: return NULL; 
          }
          break;

     //Added by MJ: case PTR, TPO, SEP
    case PTR:
        switch (atom_type) {
        case CG:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = -1; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case CD1:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG];
                case -2: offset = 0; return (*residue)[CB];
                case -3: offset = 0; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG];
                case -2: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 3;
                    return (*residue)[CA];
                }
            }
            break;
        case CD2:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG];
                case -2: offset = 0; return (*residue)[CD1];
                case -3: offset = 0; return (*residue)[CB];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG];
                case -2: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 3;
                    return (*residue)[CA];
                }
            }
            break;
        case CE1:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CD1];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CD2];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CD1];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 4;
                    return (*residue)[CA];
                }
            }
            break;
        case CE2:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CD2];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CD1];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CD2];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 4;
                    return (*residue)[CA];
                }
            }
            break;
        case CZ:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CE1];
                case -2: offset = 0; return (*residue)[CD1];
                case -3: offset = 0; return (*residue)[CG];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CE1];
                case -2: offset = 0; return (*residue)[CD1];
                case -3: offset = 0; return (*residue)[CG];
                case -4: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 5;
                    return (*residue)[CA];
                }
            }
            break;
        case OH:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CZ];
                case -2: offset = 0; return (*residue)[CE1];
                case -3: offset = 0; return (*residue)[CD1];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CZ];
                case -2: offset = 0; return (*residue)[CE1];
                case -3: offset = 0; return (*residue)[CD1];
                case -4: offset = 0; return (*residue)[CG];
                case -5: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 6;
                    return (*residue)[CA];
                }
            }
            break;
        case HB2:
        case HB3:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = 0; return (*residue)[CB];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case HD1:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CE1];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CD1];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CD1];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 4;
                    return (*residue)[CA];
                }
            }
            break;
        case HD2:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CE2];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CD2];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CD2];
                case -2: offset = 0; return (*residue)[CG];
                case -3: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 4;
                    return (*residue)[CA];
                }
            }
            break;
        case HE1:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CZ];
                case -2: offset = 0; return (*residue)[CD1];
                case -3: offset = 0; return (*residue)[CE1];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CE1];
                case -2: offset = 0; return (*residue)[CD1];
                case -3: offset = 0; return (*residue)[CG];
                case -4: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 5;
                    return (*residue)[CA];
                }
            }
            break;
        case HE2:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CZ];
                case -2: offset = 0; return (*residue)[CD2];
                case -3: offset = 0; return (*residue)[CE2];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CE2];
                case -2: offset = 0; return (*residue)[CD2];
                case -3: offset = 0; return (*residue)[CG];
                case -4: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 5;
                    return (*residue)[CA];
                }
            }
            break;
        case HH:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OH];
                case -2: offset = 0; return (*residue)[CZ];
                case -3: offset = 0; return (*residue)[CE1];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OH];
                case -2: offset = 0; return (*residue)[CZ];
                case -3: offset = 0; return (*residue)[CE1];
                case -4: offset = 0; return (*residue)[CD1];
                case -5: offset = 0; return (*residue)[CG];
                case -6: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 7;
                    return (*residue)[CA];
                }
            }
            break;
        default: return NULL;
        }
        break;
        //PTR END
    case TPO:
        switch (atom_type) {
        case OG1:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = -1; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case CG2:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = 0; return (*residue)[OG1];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case HB:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OG1];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = 0; return (*residue)[CB];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case HG1:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OG1];
                case -2: offset = 0; return (*residue)[CB];
                case -3: offset = 0; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OG1];
                case -2: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 3;
                    return (*residue)[CA];
                }
            }
            break;
        case HG21:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG2];
                case -2: offset = 0; return (*residue)[CB];
                case -3: offset = 0; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG2];
                case -2: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 3;
                    return (*residue)[CA];
                }
            }
            break;
        case HG22:
        case HG23:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG2];
                case -2: offset = 0; return (*residue)[CB];
                case -3: offset = 0; return (*residue)[HG21];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CG2];
                case -2: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 3;
                    return (*residue)[CA];
                }
            }
            break;
        default: return NULL;
        }
        break;
        //TPO END
    case SEP:
        switch (atom_type) {
        case OG:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = -1; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case HB2:
        case HB3:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OG];
                case -2: offset = 0; return (*residue)[CA];
                case -3: offset = 0; return (*residue)[CB];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 2;
                    return (*residue)[CA];
                }
            }
            break;
        case HG:
            if (iteration_type == POSITIONING) {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OG];
                case -2: offset = 0; return (*residue)[CB];
                case -3: offset = 0; return (*residue)[CA];
                default: offset = 0; return NULL;
                }
            }
            else {
                switch (offset) {
                case -1: offset = 0; return (*residue)[OG];
                case -2: offset = 0; return (*residue)[CB];
                default:
                    if (offset < 0) offset += 3;
                    return (*residue)[CA];
                }
            }
            break;
        default: return NULL;
        }
        break;
        //SEP END

     default: break;
     }
     
     return NULL;
}

// Return neighbour at given offset
// +-3 neigbours are cached when this method is called by get_neighbour()
Atom *Atom::get_neighbour_internal(int offset, definitions::IterateEnum iteration_type) {

     int actual_index = index;
     int actual_offset = offset;

     if (iteration_type==POSITIONING && (this->is_backbone_atom))
          iteration_type = BACKBONE;
     
     if (iteration_type == ALL) {
          // no need to adjust anything - we iterate over all atoms

     } else {
          // Check for neighbour constants 
          Atom *neighbour_value = get_neighbour_constants(actual_offset, iteration_type);

          if (actual_offset==0) {
               return neighbour_value;
          } else {
               if (neighbour_value)
                    actual_index = residue->atom_index[neighbour_value->atom_type];
          }                 
     }    

     int new_index = actual_index + actual_offset;

     Residue *prev_residue = NULL;
     Residue *current_residue = residue;

     int begin = 0;
     int end = current_residue->size();

     if (iteration_type == BACKBONE || iteration_type == POSITIONING) {
          // end = currentResidue->backboneAtoms;
          end = current_residue->iteration_range_indices[BACKBONE][Residue::END]+1;
     } else if (iteration_type == CA_ONLY) {
          begin = residue->atom_index[CA];
          end = begin+1;
     } else if (iteration_type == SC_ONLY) {
          begin = current_residue->iteration_range_indices[SC_ONLY][Residue::BEGIN];

          // this stunt is necessary to be able to skip over the
          // the gly residues a bit easier
          if (current_residue->residue_type != GLY) {
               end = current_residue->iteration_range_indices[SC_ONLY][Residue::END]+1 ;
          } else {
               end = current_residue->iteration_range_indices[SC_ONLY][Residue::END];
          }
     }

     // Backwards movement through chain
     while (current_residue != NULL && (new_index < begin)) {
          prev_residue = current_residue;
          current_residue = current_residue->get_neighbour(-1);
          
          if (current_residue != NULL) {

               if (iteration_type == SC_ONLY) {
                    begin = current_residue->iteration_range_indices[SC_ONLY][Residue::BEGIN];
               }               
               
               if (iteration_type==ALL) {
                    end = current_residue->size();
               } else if (iteration_type == BACKBONE || iteration_type == POSITIONING) {
                    // end = currentResidue->backboneAtoms;
                    end = current_residue->iteration_range_indices[BACKBONE][Residue::END]+1;
               } else if (iteration_type == SC_ONLY) {
                    if (current_residue->residue_type != GLY) {
                         end = current_residue->iteration_range_indices[SC_ONLY][Residue::END]+1 ;
                    } else {
                         end = current_residue->iteration_range_indices[SC_ONLY][Residue::END];
                    }
               }
               new_index += end-begin;
          }
     }

     // Forwards movement through chain
     while (current_residue != NULL && (new_index >= end)) {
          prev_residue = current_residue;
          current_residue = current_residue->get_neighbour(+1);

          if (current_residue != NULL) {
               
               if (iteration_type == SC_ONLY) {
                    begin = current_residue->iteration_range_indices[SC_ONLY][Residue::BEGIN];
               }
               
               new_index -= end-begin;

               if (iteration_type==ALL) {
                   end = current_residue->size();
               } else if (iteration_type == BACKBONE || iteration_type == POSITIONING) {
                    // end = currentResidue->backboneAtoms;
                    end = current_residue->iteration_range_indices[BACKBONE][Residue::END]+1;
               } else if (iteration_type == SC_ONLY) {
                    if (current_residue->residue_type != GLY) {
                         end = current_residue->iteration_range_indices[SC_ONLY][Residue::END]+1;
                    } else {
                         end = current_residue->iteration_range_indices[SC_ONLY][Residue::END];
                    }
               }
          }
     }

     if (new_index >= begin && new_index < end) {
          return (current_residue->atoms[new_index]);
     } else {
          return NULL;
     }

}

// Wrapper for get_neighbour_internal
// Checks whether neighbour is available in cache
Atom *Atom::get_neighbour(int offset, definitions::IterateEnum iteration_type) {

     // Detect whether number of atoms in chain has changed
     // if so, clear the cache
     if (residue->atoms_length != residue->size()) {
          residue->clear_neighbour_cache();
          residue->atoms_length = residue->size();
     }

     // Special case if offset value is not cachable
     if (offset < -3 || offset > 3) {
          return get_neighbour_internal(offset, iteration_type);
     }

     offset+=3;
     
     // If no cached neighbour value is available, call
     // get_neighbour_internal to fill the cache
     if (neighbour_atom[offset][iteration_type] == NULL ) {
          neighbour_atom[offset][iteration_type] = get_neighbour_internal(offset-3, iteration_type);
     }

     // Return cached value
     return neighbour_atom[offset][iteration_type];
}


// Delete all entries in the cache
void Atom::clear_neighbour_cache() {
     for (int i=0; i<definitions::ITERATE_ENUM_SIZE; i++) {
         neighbour_atom[0][i] = NULL;
         neighbour_atom[1][i] = NULL;
         neighbour_atom[2][i] = NULL;
         neighbour_atom[3][i] = NULL;
         neighbour_atom[4][i] = NULL;
         neighbour_atom[5][i] = NULL;
         neighbour_atom[6][i] = NULL;
    }
}


// Set angles in current atom based on surrounding atoms
// NOTE: Should only be called when the position of left and right atom neighbours are initialized
void Atom::update_angle(bool fix_end_points) {

     Atom *atom1, *atom2, *atom3;
     get_angle_atoms(&atom1, &atom2, &atom3);
     
     // Check if this atom keeps track of the previous angle itself
     if (this->self_maintained_positioning) {

          Vector_3D pos_atom1 = atom1->position; // PrevPrev
          Vector_3D pos_atom2 = atom2->position; // Prev
          Vector_3D pos_atom3 = atom3->position; // Current

          // Test if POSITIONING neighbouring atom is different
          // from the real neighbour
          // if so, correct the corresponding position vector
          Atom *prev_real = get_neighbour(-1, BACKBONE);
          if (prev_real != atom2) {
               pos_atom3 += (pos_atom2 - prev_real->position);
          }

          // Calculate angle from neighbouring positions
          this->set_angle(calc_angle(pos_atom1,
                                    pos_atom2,
                                    pos_atom3));


     // Otherwise, use default placing method     
     } else {

          if (atom1 == NULL || atom3 == NULL) {
               // If fix_end_points options is set, the angle values
               // at the non-well-defined positions are set to zero
               if (fix_end_points) {
                    this->set_angle(0.0);
               }
          } else {

               // Calculate angle from neighbouring positions
               this->set_angle(calc_angle(atom1->position,
                                          atom2->position,
                                          atom3->position));
          }
     }
}

// Set dihedral based on positions of neighbouring atoms
// NOTE: Should only be called when neighbouring position are initialized
void Atom::update_dihedral(bool fix_end_points) {

     Atom *atom1, *atom2, *atom3, *atom4;
     get_dihedral_atoms(&atom1, &atom2, &atom3, &atom4);
     
     // Check if this atom keeps track of the previous dihedral itself
     // if (this->prevAtomDihedral) {
     if (this->self_maintained_positioning) {
     
          Vector_3D pos_atom1 = atom1->position; // prevPrevPrev
          Vector_3D pos_atom2 = atom2->position; // prevPrev
          Vector_3D pos_atom3 = atom3->position; // prev
          Vector_3D pos_atom4 = atom4->position; // current

          // Test if POSITIONING neighbouring atom is different
          // from the real neighbour
          // if so, correct the corresponding position vector
          Atom *prev_real = get_neighbour(-1, BACKBONE);
          if (prev_real != atom3) {
               pos_atom4 += (pos_atom3 - prev_real->position);
          }

          // Calculate dihedral from neighbouring positions
          this->set_dihedral(calc_dihedral(pos_atom1,
                                           pos_atom2,
                                           pos_atom3,
                                           pos_atom4));
          
     // Otherwise, use default placing method 
     } else {

         if (atom1 == NULL || atom2 == NULL || atom4 == NULL) {
               // If fixEndPoints options is set, the angle values
               // at the non-well-defined positions are set to zero
              if (fix_end_points)
                   this->set_dihedral(M_PI);
         } else {
              // Calculate dihedral from neighbouring positions
              this->set_dihedral(calc_dihedral(atom1->position,
                                               atom2->position,
                                               atom3->position,
                                               atom4->position));
         }
     }
}

// Set position based on angles of neighbouring atoms
void Atom::update_position(int direction) {
     
     int dir = direction>=0 ? 1 : -1;

     Atom *prev_prev_prev = NULL;
     Atom *prev_prev = NULL;
     Atom *prev = NULL;

     if (this->self_maintained_positioning) {
          prev_prev_prev = get_neighbour(-3, POSITIONING);
          prev_prev = get_neighbour(-2 , POSITIONING);
          prev = get_neighbour(-1 , POSITIONING);
     } else {
          prev_prev_prev = get_neighbour(-3 * dir, POSITIONING);
          prev_prev = get_neighbour(-2 * dir, POSITIONING);
          prev = get_neighbour(-1 * dir, POSITIONING);
     }

     Atom *current = this;

     if (prev==NULL) {            // First atom in chain
          if (!position.initialized)
               current->position = Vector_3D(0,0,0);
     } else if (prev_prev == NULL) {
          // If the position has previously been initialized. Use the
          // old direction to place the atom again - making sure that
          // the bond length is idealized
          //
          double bondLength = direction>0 ? this->bond_length : prev->bond_length;
          if (position.initialized) {
               current->position = prev->position + ((current->position - prev->position).normalize()* bondLength);
          } else {
               current->position = Vector_3D(bondLength, 0, 0);
          }

     } else {
          
          double angle;
          double dihedral;

          // Normally, the atom angle at position i is defined by the
          // position of atoms i-1, i and i+1, and the dihedral at position
          // i is defined by positions i-2, i-1, i and i+1. In some cases
          // (e.g. in sidechains) it is more convenient to let the atom
          // keep track of the angle and dihedral necessary to position
          // itself, i.e., the angle of atom i is defined by position i-2,
          // i-1, i, and the dihedral by i-3, i-2, i-1, i. In this latter
          // case, self_maintained_positioning=true         
          if (self_maintained_positioning) {
               angle = this->get_angle();
          } else {
               angle = prev->get_angle();
          }

          if (self_maintained_positioning) {
               dihedral = this->get_dihedral();
          } else {
               dihedral = dir>0 ? *prev->dihedral : *prev_prev->dihedral;
          }
          
          double bond_length = direction>0 ? this->bond_length : prev->bond_length;
          
          // Vectors used as basis
          Vector_3D &v1 = prev->position;
          Vector_3D &v2 = prev_prev->position;
          Vector_3D v3;
          if (prev_prev_prev == NULL) {
               // If the position has previously been initialized. Use
               // plane defined by the old position and the two first
               // positions to place the atom.
               if (position.initialized) {
                    v3 = v2 + (current->position-v2);
                    dihedral = 0;
               } else {
                    v3 = (v2 - v1);
                    v3 = Vector_3D(v3[1], -v3[0], 0)+v1;
               }
          } else {
               v3 = prev_prev_prev->position;
          }

          current->position = prev_prev->position;

          Vector_3D D(bond_length * cos(M_PI-angle),
                      bond_length * cos(M_PI-dihedral) * sin(M_PI-angle),
                      bond_length * sin(M_PI-dihedral) * sin(M_PI-angle));

          Vector_3D bc = (v1-v2).normalize();
          Vector_3D n = cross_product((v1-v3), bc).normalize();
          Vector_3D nbc = cross_product(bc,n);
          Matrix_3D basis_change(bc, nbc, n, true);
          D = basis_change*D + v1;

          // Test if POSITIONING neighbouring atom is different from the real neighbour
          // if so, correct the corresponding position vector
          Atom *prev_real = get_neighbour(-1 * dir, BACKBONE);
          if (prev_real != prev) {
               // std::cout << "Adjusting position for origo offset\n";
               D -= (prev->position - prev_real->position);
          }
          current->position = D;
     }
}

// Return angle
double &Atom::get_angle() {
     return *angle;
}

// Set angle
// Also sets savedAngle if angle was not yet initialized
void Atom::set_angle(double angle) {
     *this->angle = angle;
}

// Return dihedral
double &Atom::get_dihedral() {
     return *dihedral;
}

// Set dihedral
// Also sets savedDihedral if dihedral was not yet initialized
void Atom::set_dihedral(double dihedral) {
     *this->dihedral = dihedral;
}

// // return position
// Vector_3D &Atom::getPosition() {
//      return position;
// }

// // Set position
// // Also sets savedPosition if position was not yet initialized
// void Atom::setPosition(Vector_3D position) {
//      this->position = position;
// }


// Translate position with a translation vector
void Atom::translate(Vector_3D translation) {
     position += translation;
}

// Rotate position by a rotation matrix
void Atom::rotate(Matrix_3D rotation) {
     position = rotation * position;
}

// Check atom consistency
// Only included in compilation when DEBUGLEVEL > 0
void Atom::check_consistency() {

#if DEBUGLEVEL > 0
     
     Atom *current_atom = this;
     Atom *prev_atom = get_neighbour(-1, POSITIONING);
     Atom *prev_prev_atom = get_neighbour(-2, POSITIONING);
     Atom *prev_prev_prev_atom = get_neighbour(-3, POSITIONING);
     Atom *next_atom = get_neighbour(+1, POSITIONING);

     Atom *prev_atom_real = get_neighbour(-1, BACKBONE);
     
     if (prev_atom) {
          Vector_3D &current_pos = current_atom->position;
          Vector_3D &prev_pos_real = prev_atom_real->position;

          double bond_length = (current_pos - prev_pos_real).norm();
          check_equality_assert(this->bond_length, bond_length);
          
          // Check if this atom keeps track of the previous angle itself
          if (this->self_maintained_positioning) {

               if (prev_atom && prev_prev_atom) {
                    Vector_3D pos_prev_prev = prev_prev_atom->position;
                    Vector_3D pos_prev      = prev_atom->position;
                    Vector_3D pos_current   = current_atom->position;

                    // Test if POSITIONING neighbouring atom is different from the real neighbour
                    // if so, correct the corresponding position vector
                    if (prev_atom_real != prev_atom) {
                         pos_current += (pos_prev - prev_atom_real->position);
                    }
                    check_equality_angle_assert(get_angle(), calc_angle(pos_prev_prev, pos_prev, pos_current));
               }

               if (prev_atom && prev_prev_atom && prev_prev_prev_atom) {
                    Vector_3D pos_prev_prev_prev = prev_prev_prev_atom->position;
                    Vector_3D pos_prev_prev      = prev_prev_atom->position;
                    Vector_3D pos_prev          = prev_atom->position;
                    Vector_3D pos_current       = current_atom->position;

                    // Test if POSITIONING neighbouring atom is different from the real neighbour
                    // if so, correct the corresponding position vector
                    if (prev_atom_real != prev_atom) {
                         pos_current += (pos_prev - prev_atom_real->position);
                    }
                    
                    check_equality_angle_assert(get_dihedral(), calc_dihedral(pos_prev_prev_prev, pos_prev_prev, pos_prev, pos_current));
               }
               
          // Otherwise, use default placing method        
          } else {

               if (prev_atom && next_atom) {

                    Vector_3D pos_prev     = prev_atom->position;
                    Vector_3D pos_current  = current_atom->position;
                    Vector_3D pos_next     = next_atom->position;

                    // checkEquality_angle(get_angle(), calc_angle(posPrev, posCurrent, posNext));
                    check_equality_angle_assert(get_angle(), calc_angle(pos_prev, pos_current, pos_next));
               }

               if (prev_atom && prev_prev_atom && next_atom) {
                    Vector_3D pos_prev_prev = prev_prev_atom->position;
                    Vector_3D pos_prev      = prev_atom->position;
                    Vector_3D pos_current   = current_atom->position;
                    Vector_3D pos_next      = next_atom->position;
                    
                    check_equality_angle_assert(get_dihedral(), calc_dihedral(pos_prev_prev, pos_prev, pos_current, pos_next));
               }
               
          }
     }
#endif
}

// Output in PDB format
std::string Atom::output_as_pdb(int counter_offset, const int chain_number, const std::string &b_factor_string) {
     char str[500];
     std::string atom_label = atom_name[atom_type];
     if (atom_label.size() < 4)
          atom_label = " " + atom_label;

     double b_factor = 0.0;
     if (b_factor_string != "") {
          b_factor = boost::lexical_cast<double>(b_factor_string);
     }

     sprintf(str,"ATOM   %4d %-4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00%6.2f\n",
             counter_offset + index + 1, atom_label.c_str(),
             residue_name[residue->residue_type], static_cast<char>(static_cast<int>('A')+chain_number-1),
             residue->index_res_seq, position[0], position[1], position[2], b_factor);
     
     std::string output = std::string(str);
     
     return output;
}

// Output for atom reference
std::ostream & operator<<(std::ostream &o, Atom &a) {
     o << &a;
     return o;
}

// Output for atom pointer
std::ostream & operator<<(std::ostream &o, Atom *a) {
     if (a != NULL)
          o << "Atom[" << a->index << "][" << a->atom_type << "]: (" << a->position << ") in Residue " << a->residue->residue_type << "(" << a->residue->index << ")";
     else
          o << "Atom[NULL]";
     return o;
}

}

