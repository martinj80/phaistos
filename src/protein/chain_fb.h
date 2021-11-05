// chain_fb.h --- Full-atom backbone chain (inherits from Chain)
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

#include "chain.h"

#ifndef CHAIN_FB_H
#define CHAIN_FB_H

namespace phaistos {

//! Full atom chain class
class ChainFB: public Chain<ResidueFB> {
private:

     //! Initializer based on ProteinData object
     //!
     //! \param data ProteinData parser object
     //! \param atom_selection Optional specification of additional atoms to be added
     //! \param chain_index If data contains multiple chains, specify chain_index to select which one to use.
     void init_from_protein_data(ProteinData &data,
                                 definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS +
                                                                                definitions::BACKBONE_H_ATOMS +
                                                                                definitions::CB_ATOMS),
                                 int chain_index=-1) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          if (data.n_polypeptides() > 1) {
               // If chain index was specified, check if it is within range - otherwise issue warning
               if (!((chain_index >= 0) && (chain_index < data.n_polypeptides()))) {
                    fprintf(stderr, "Chain contains multiple chains or chain breaks. Only using first chain segment\n");
                    chain_index = 0;
               }
          } else {
               chain_index = 0;
          }

          bool toggle_PS = false;

          if ((atom_selection & PSEUDO_SIDECHAIN_ATOMS) && !(atom_selection & SIDECHAIN_ATOMS)) {
               // in order to calculate the ps atom position from the
               // full atom positions, we need to read those in first.
               toggle_PS = true;
               atom_selection = atom_selection+SIDECHAIN_ATOMS;
          }

          // Set size of residues array
          residues.resize(data[chain_index].size());
          
          for (unsigned int i=0; i<data[chain_index].size(); i++) {

               if (data[chain_index][i].size() >=3) {
                    ResidueFB *residue;

                    std::vector<std::pair<AtomEnum, Vector_3D> > positions;
                    std::vector<std::string> alt_locs;

                    // Keep track of duplicates
                    std::string active_alt_loc = "";
                    std::vector<int> added_atoms(ATOM_ENUM_SIZE, -1);

                    std::string aa_type_str = *data[chain_index][i][0]->residue;
                    ResidueEnum aa_type = ResidueEnum(str_to_aa(aa_type_str));

                    for (unsigned int k=0; k<data[chain_index][i].size(); k++) {
                         AtomEnum atom_type = string_to_atom( (*data[chain_index][i][k]->atom), aa_type_str );

                         if (atom_type == XX)
                              continue;

                         if (atom_type == O && !(atom_selection & BACKBONE_O_ATOMS))
                              continue;

                         if ((atom_type == H || atom_type==H1 || atom_type==H2 || atom_type==H3) &&
                             !(atom_selection & BACKBONE_H_ATOMS))
                              continue;

                         if (atom_type == CB && !(atom_selection & CB_ATOMS)) 
                              continue;

                         if (is_sidechain_atom_type(atom_type) && (atom_type != CB) && !(atom_selection & SIDECHAIN_ATOMS) && !(atom_selection & PSEUDO_SIDECHAIN_ATOMS && atom_type == PS) )  {
                              continue;
                         }

                         // if (atom_type>H3 && atom_type<=HZ3 && !(atom_selection & NON_BACKBONE_H_ATOMS))
                         if (atom_type>HA3 && atom_type<=HZ3 && !(atom_selection & NON_BACKBONE_H_ATOMS))
                              continue;
                         
                         alt_locs.push_back(*data[chain_index][i][k]->altloc);
                         positions.push_back(std::pair<AtomEnum, Vector_3D>(atom_type,
                                                                            Vector_3D(data[chain_index][i][k]->coordinates)));

                         // Check for duplicates
                         if (added_atoms[atom_type] == -1) {
                              added_atoms[atom_type] = positions.size()-1;
                         } else {
                              active_alt_loc = alt_locs[added_atoms[atom_type]];
                              std::cerr << "WARNING: Residue " << i << " has multiple occurences of same atomtype - using altloc: " << active_alt_loc << "\n";
                              continue;
                         }
                         
                    }


                    // Remove multiple occurances of same atom_type
                    std::vector<std::pair<AtomEnum, Vector_3D> > positions_unique;
                    for (unsigned int k=0; k<positions.size(); k++) {
                    
                         // If there is an active altloc, only accept entries with
                         // empty or matchin alt_locs
                         if (active_alt_loc!="" && alt_locs[k]!="_" && alt_locs[k] != active_alt_loc) {
                              continue;
                         }
                         positions_unique.push_back(positions[k]);
                    }

                    residue = new ResidueFB(aa_type, this, i, positions_unique, *data[chain_index][i][0]->resseq);

                    residues[i] = residue;
               } else {
                    if (i==0) {
                         fprintf(stderr, "Incomplete residue at position 0. Quitting.\n");
                         exit(1);
                    }
                    fprintf(stderr, "Error in ChainFB: Encountered incomplete residue at position %d. Using residues 0-%d\n", i, i-1);
                    break;
               }
          }

          if (toggle_PS) {
               // and remove all the sidechains atoms again
               // and switch to the PS representation
               atom_selection = atom_selection-SIDECHAIN_ATOMS;
               this->toggle_sidechains_to_PS() ;
          }

          update_angles();
     }


public:

     //! Constructor
     ChainFB() : Chain<ResidueFB>() {}

     //! Destructor
     ~ChainFB() {}

     //! Constructor - Initialize chain based on list of positions
     //!
     //! \param positions Vector of 3D-coordinates
     //! \param sequence Optional amino acid vector
     //! \param res_seq vector of residue indices (RESSEQ entry in PDB)
     ChainFB(const std::vector<Vector_3D> positions,
             const std::vector<int> sequence=std::vector<int>(),
             const std::vector<int> res_seq=std::vector<int>()) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          int chain_length = positions.size()/3;

          // Set size of residues array
          residues.resize(chain_length);
          
          for (int i=0; i<chain_length; i++) {
               ResidueFB *residue;
               Vector_3D position_N  = positions[i*3 + 0];
               Vector_3D position_CA = positions[i*3 + 1];
               Vector_3D position_C  = positions[i*3 + 2];

               ResidueEnum aa;
               if (sequence.size() != 0)
                    aa = ResidueEnum(sequence[i]);
               else
                    aa = ResidueEnum(ALA);

               int res_seq_index;
               if (res_seq.size() != 0)
                    res_seq_index = res_seq[i];
               else
                    res_seq_index = i;

               residue = new ResidueFB(aa, this, i, position_N, position_CA, position_C, res_seq_index);

               residues[i] = residue;
          }

          // Update angles based on positions
          update_angles();

          // Add non-backbone atoms
          add_atoms(BACKBONE_O_ATOMS + BACKBONE_H_ATOMS + CB_ATOMS);
          update_positions();
     }

     //! Constructor - Initialize chain based on list of angles.
     //!
     //! \param angles Vector of angles
     //! \param sequence Optional Amino acid vector
     //! \param res_seq vector of residue indices (RESSEQ entry in PDB)
     //! \param atom_selection Optional specification of additional atoms to be added
     //! \param sidechain_dof_values Optionally set values for sidechain degrees of freedom
     ChainFB(const std::vector<std::vector<double> > &angles,
             std::vector<int> sequence=std::vector<int>(),
             const std::vector<int> res_seq=std::vector<int>(),
             definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS +
                                                            definitions::BACKBONE_H_ATOMS +
                                                            definitions::CB_ATOMS),
             const std::vector<std::vector<double> > &sidechain_dof_values=std::vector<std::vector<double> >()) {
          
          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Set size of residues array
          residues.resize(angles.size());
          
          for (unsigned int i=0; i<angles.size(); i++) {
               ResidueFB *residue;

               ResidueEnum aa;
               if (sequence.size() != 0)
                    aa = ResidueEnum(sequence[i]);
               else
                    aa = ResidueEnum(ALA);

               int res_seq_index;
               if (res_seq.size() != 0)
                    res_seq_index = res_seq[i];
               else
                    res_seq_index = i;

               residue = new ResidueFB(aa, this, i, angles[i], res_seq_index);
               residues[i] = residue;
          }

          add_atoms(atom_selection, sidechain_dof_values);
          
          update_positions();
     }



     //! Constructor - initialize based on PDB file
     //!
     //! \param pdb_filename Name of PDB file
     //! \param atom_selection Optional specification of additional atoms to be added
     ChainFB(std::string pdb_filename,
             definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS + 
                                                            definitions::BACKBONE_H_ATOMS +
                                                            definitions::CB_ATOMS)) {
          
          ProteinData protein_data = read_pdb_input((char *)pdb_filename.data());
          init_from_protein_data(protein_data, atom_selection);
          this->name = pdb_filename;
     }

     //! Constructor - initialize based on ProteinData object
     //!
     //! \param data ProteinData parser object
     //! \param atom_selection Optional specification of additional atoms to be added
     //! \param chain_index If data contains multiple chains, specify chain_index to select which one to use.
     ChainFB(ProteinData &data,
             definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS + 
                                                            definitions::BACKBONE_H_ATOMS +
                                                            definitions::CB_ATOMS),
             int chain_index=-1) {
          init_from_protein_data(data, atom_selection, chain_index);
     }

     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     ChainFB(const ChainFB &other)
          : Chain<ResidueFB>() {

          // Set size of residues array
          residues.resize(other.residues.size());
          
          for (int i=0; i<other.size(); i++) {
               ResidueFB *residue = new ResidueFB(*other.residues[i], this, i, other.residues[i]->index_res_seq);
               residues[i] = residue;
          }
     }


     //! Copy constructor - sub-chain
     //!
     //! \param other Source object from which copy is made.
     //! \param start_index Where in the sequence to start
     //! \param end_index Where in the sequence to end
     ChainFB(const ChainFB &other, int start_index, int end_index)
          : Chain<ResidueFB>() {

          if (start_index < 0)
               start_index=0;
          if (end_index < 0)
               end_index = other.size();

          // Set size of residues array
          residues.resize(end_index-start_index);
          
          for (int i=start_index; i<end_index; i++){
               ResidueFB *residue = new ResidueFB(*other.residues[i], this, i-start_index);
               residues[i-start_index] = residue;
          }
     }

     //! Clone - virtual function that acts as copy constructor
     ChainFB *clone() {
          return new ChainFB(*this);
     }
};

}

#endif
