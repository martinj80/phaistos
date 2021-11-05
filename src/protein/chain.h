// chain.h --- protein chain base class
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


#include <boost/algorithm/string/replace.hpp>

#include "residue.h"
#include "atom.h"
#include "protein/pdb_input.h"
#include "utils/utils.h"
#include "iterators/atom_iterator.h"
#include "iterators/residue_iterator.h"

#ifndef CHAIN_H
#define CHAIN_H

namespace phaistos {

// Forward declarations
namespace chaintree {
class RSS;
template<typename CHAIN_TYPE, typename BV_TYPE> class ChainTree;
}
template<typename RESIDUE_TYPE> class ResidueIterator;

//! Chain class.
template<typename RESIDUE_TYPE>
class Chain {
protected:

     //! Vector of residue objects
     std::vector<RESIDUE_TYPE *> residues;

     //! Chain name
     std::string name;

     //! Residue iterator
     typedef typename std::vector<RESIDUE_TYPE *>::iterator Iterator;

     //! Residue iterator - reverse
     typedef typename std::vector<RESIDUE_TYPE *>::reverse_iterator ReverseIterator;

     //! Return iteration pointing to beginning of residue vector.
     //!
     //! \return begin iterator
     Iterator begin() {
          return residues.begin();
     }

     //! Return iteration pointing to end of residue vector.
     //!
     //! \return  end iterator
     Iterator end() {
          return residues.end();
     }

     //! Return iteration pointing to end of residue vector (for reverse iteration).
     //!
     //! \return begin iterator (reverse: points to end of vector)
     ReverseIterator rbegin() {
          return residues.rbegin();
     }

     //! Return iteration pointing to beginning of residue vector (for reverse iteration).
     //!
     //! \return end iterator (reverse: points to beginning of vector)
     ReverseIterator rend() {
          return residues.rend();
     }

public:

     //! Allow direct access to ResidueIterator
     friend class ::phaistos::ResidueIterator<typename RESIDUE_TYPE::ChainType>;

     //! Time stamp
     long int time_stamp;

     //! Local version of Residue type.
     //! This is to avoid passing along a RESIDUE_TYPE parameter in situations where 
     //! CHAINTYPE is available.
     typedef RESIDUE_TYPE Residue;

     //! Local version of ChainTree type.
     typedef chaintree::ChainTree<Chain<RESIDUE_TYPE> , chaintree::RSS> ChainTree;

     //! Chain tree pointer (activated by get_chain_tree method)
     ChainTree *chain_tree;

     //! Constructor
     Chain() {
          chain_tree = NULL;
          time_stamp = 0;
          name = "";
     }

     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     // Note: this function has external requirement not available at
     // this point in the code. It is included at the end of the file
     Chain(const Chain &other);

     //! Destructor
     virtual ~Chain() {
          for (Iterator it = begin(); it != end(); ++it) {
               delete (*it);
          }
          if (chain_tree)
               delete chain_tree;
     }

     //! Clone - virtual function that acts as copy constructor
     virtual Chain<RESIDUE_TYPE> *clone() = 0;

     //! Return size (number of residues)
     //!
     //! \return number of residues
     int size() const {
          return residues.size();
     }

     //! Get amino acid sequence - either as
     //! vector of int or ResidueEnums
     //!
     //! \tparam RESIDUEENUM_OR_INT Whether to return int or ResidueEnum values
     //! \return amino acid sequence
     template<typename RESIDUEENUM_OR_INT>
     std::vector<RESIDUEENUM_OR_INT> get_sequence() {
          std::vector<RESIDUEENUM_OR_INT> sequence;
          for (Iterator it = begin(); it != end(); ++it) {
               sequence.push_back((*it)->residue_type);
          }
          return sequence;
     }

     //! Get current angles
     //!
     //! \param include_omega Whether omega angles should be included in the output vector
     //! \param include_bond_angles Whether bond angles should be included in the output vector.
     //! \return vector of angles.
     std::vector<std::vector<double> > get_angles(bool include_omega = false,
                                                  bool include_bond_angles = false) {
          std::vector<std::vector<double> > angles;
          for (Iterator it = begin(); it != end(); ++it) {
               angles.push_back((*it)->get_angles(include_omega, include_bond_angles));
          }
          return angles;
     }

     //! Set angles
     //!
     //! \param angles Vector of angles for each residue in range
     //! \param start_index Where in the sequence to start
     //! \param end_index Where in the sequence to end
     void set_angles(std::vector<std::vector<double> > angles, int start_index = 0, int end_index = -1) {
          if (end_index == -1)
               end_index = size() + 1;
          std::vector<std::vector<double> >::iterator it_angles = angles.begin() + start_index;
          for (Iterator it_res = get_iterator_at_index(start_index);
               it_res != get_iterator_at_index(end_index); (++it_res, ++it_angles)) {
               (*it_res)->set_angles(*it_angles);
          }
     }

     //! Set the chain name
     //!
     //! \param name New name
     void set_name(std::string name) {
          this->name = name;
     }

     //! Get the chain name
     //!
     //! \return Chain name
     std::string get_name() const {
          return this->name;
     }

     //! Update angle values based on current positions.
     //!
     //! \param fix_end_points If true, the angle value at the non-well-defined 
     //!        endpoint positions are set to zero.
     void update_angles(bool fix_end_points = true) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->update_angles(fix_end_points);
          }
     }

     //! Get side chain degree-of-freedom values.
     //! 
     //! \param mode SIDECHAIN_ATOMS -> dofValueVector contains chi angles
     //!              PSEUDO_SIDECHAIN_ATOMS -> dofValueVector contains (bondlength, angle, dihedral)
     //! \param start_index Where in the sequence to start
     //! \param end_index Where in the sequence to end
     //! \return vector of degree-of-freedom values for each residue in range
     std::vector<std::vector<double> > get_sidechain_dof_values(definitions::AtomSelectionEnum mode = definitions::SIDECHAIN_ATOMS, 
                                                                int start_index = 0, int end_index = -1) {
          if (end_index == -1)
               end_index = size() + 1;

          std::vector<std::vector<double> > dof_value_vector;

          for (Iterator it = get_iterator_at_index(start_index); it != get_iterator_at_index(end_index); ++it) {
               dof_value_vector.push_back((*it)->get_sidechain_dof_values(mode));
          }
          return dof_value_vector;
     }

     //! Set side chain degree-of-freedom values
     //! 
     //! \param dof_value_vector Vector of input values
     //! \param mode SIDECHAIN_ATOMS -> dofValueVector contains chi angles
     //!             PSEUDO_SIDECHAIN_ATOMS -> dofValueVector contains (bondlength, angle, dihedral)
     //! \param start_index Where in the sequence to start
     //! \param end_index Where in the sequence to end
     void set_sidechain_dof_values(std::vector<std::vector<double> > &dof_value_vector,
                                   definitions::AtomSelectionEnum mode = definitions::SIDECHAIN_ATOMS,
                                   int start_index = 0,
                                   int end_index = -1) {
          if (end_index == -1)
               end_index = size() + 1;

          int i = 0;
          for (Iterator it = get_iterator_at_index(start_index); it != get_iterator_at_index(end_index); ++it) {
               (*it)->set_sidechain_dof_values(dof_value_vector[i], mode);
               i++;
          }
     }

     //! Find the update direction that involvest fewest residues
     //!
     //! \param modified_angles_start_index First residue that has been modified
     //! \param modified_angles_end_index Last residue that has been modified
     //! \param fixed_direction Overrides default behaviour by directly setting a direction.
     //! \return direction
     int find_shortest_direction(int modified_angles_start_index,
                                 int modified_angles_end_index,
                                 int fixed_direction = 0) {
          int direction = 1;
          int begin_offset = modified_angles_start_index;
          int end_offset = size() - 1 - modified_angles_end_index;

          // Set direction
          if (fixed_direction != 0)
               direction = fixed_direction;
          else if (begin_offset < end_offset)
               direction = -1;

          return direction;
     }

     //! Return vector of position values
     //!
     //! \tparam ITERATION_TYPE How to iterate over the chain when collection positions.
     //! \return vector of 3D-coordinates
     template<definitions::IterateEnum ITERATION_TYPE>
     std::vector<Vector_3D> get_positions() {
          std::vector<Vector_3D> positions;
          for (AtomIterator<typename RESIDUE_TYPE::ChainType, ITERATION_TYPE, Vector_3D> it(*(typename RESIDUE_TYPE::ChainType*) this);
               !it.end(); ++it) {
               positions.push_back(*it);
          }
          return positions;
     }

     //! Set position values
     //! 
     //! \param positions Vector of 3D-coordinates
     void set_positions(std::vector<Vector_3D> positions) {
          std::vector<Vector_3D>::iterator it_pos = positions.begin();
          for (Iterator it_res = begin(); it_res != end() && it_pos != positions.end(); ++it_res) {
               (*it_res)->set_positions(positions, it_pos); // it_pos is incremented by residue
          }
     }

     //! Update position values based on internal degrees of freedom - range specified
     //!
     //! \param modified_angles_start_index First residue that has been modified
     //! \param modified_angles_end_index Last residue that has been modified
     //! \param fixed_direction Overrides automatic detection of direction in which to update.
     //! \return Direction
     int update_positions(int modified_angles_start_index, 
                          int modified_angles_end_index, 
                          int fixed_direction = 0) {

          int begin_offset = modified_angles_start_index;
          int end_offset = size() - 1 - modified_angles_end_index;
          int direction = find_shortest_direction(modified_angles_start_index, modified_angles_end_index, fixed_direction);

          // The update procedure is slit in two (first backbone - then the rest), because some
          // of the non backbone atoms are dependent on atoms outside the boundary of the current
          // residue (for instance the H and C atoms)
          if (direction > 0) {
               for (Iterator it = begin() + begin_offset; it != end(); ++it) {
                    (*it)->update_positions_backbone(direction);
               }
               for (Iterator it = begin() + begin_offset; it != end(); ++it) {
                    (*it)->update_positions_non_backbone();
               }
          } else {
               for (ReverseIterator it = rbegin() + end_offset; it != rend(); ++it) {
                    (*it)->update_positions_backbone(direction);
               }
               for (ReverseIterator it = rbegin() + end_offset; it != rend(); ++it) {
                    (*it)->update_positions_non_backbone();
               }
          }
          return direction;
     }

     //! Update position values based on internal degrees of freedom
     void update_positions() {
          update_positions(0, size() - 1);
     }

     //! Update position values based on internal degrees of freedom - backbone only
     //!
     //! \param modified_angles_start_index First residue that has been modified
     //! \param modified_angles_end_index Last residue that has been modified
     //! \param fixed_direction Overrides automatic detection of direction in which to update.
     //! \return Direction
     int update_positions_backbone(int modified_angles_start_index, int modified_angles_end_index, int fixed_direction = 0) {

          int begin_offset = modified_angles_start_index;
          int end_offset = size() - 1 - modified_angles_end_index;
          int direction = find_shortest_direction(modified_angles_start_index, modified_angles_end_index, fixed_direction);

          // The update procedure is slit in two (first backbone - then the rest), because some
          // of the non backbone atoms are dependent on atoms outside the boundary of the current
          // residue (for instance the H and C atoms)
          if (direction > 0) {
               for (Iterator it = begin() + begin_offset; it != end(); ++it) {
                    (*it)->update_positions_backbone(direction);
               }
          } else {
               for (ReverseIterator it = rbegin() + end_offset; it != rend(); ++it) {
                    (*it)->update_positions_backbone(direction);
               }
          }
          return direction;
     }

     //! Updates positions from start to end (not included), and only these. Direction is always +1
     //! and the call will invariably introduce chain breakage which must later be closed.
     //!
     //! \param start_index Where in the sequence to start
     //! \param end_index Where in the sequence to end
     void update_positions_backbone_segment(int start_index, int end_index) {
          int direction = 1;
          Iterator begin = get_iterator_at_index(start_index);
          Iterator end = get_iterator_at_index(end_index);
          for (Iterator it = begin; it != end; ++it) {
               (*it)->update_positions_backbone(direction);
          }
     }

     //! Add atoms to chain
     //!
     //! \param atom_selection Which atoms to add
     //! \param sidechain_dof_values Optionally set values for sidechain degrees of freedom
     void add_atoms(definitions::AtomSelectionEnum atom_selection = (definitions::BACKBONE_O_ATOMS + 
                                                                     definitions::BACKBONE_H_ATOMS +
                                                                     definitions::CB_ATOMS),
                    const std::vector<std::vector<double> > &sidechain_dof_values = std::vector<std::vector<double> >()) {

          for (Iterator it = begin(); it != end(); ++it) {
               if ((*it)->index < (int) sidechain_dof_values.size()) {
                    (*it)->add_atoms(atom_selection, sidechain_dof_values[(*it)->index]);
               } else {
                    (*it)->add_atoms(atom_selection);
               }
          }
          update_positions();
     }

     //! Remove atoms from chain
     //!
     //! \param atom_selection Which atoms to remove
     void remove_atoms(definitions::AtomSelectionEnum atom_selection = (definitions::BACKBONE_O_ATOMS + 
                                                                        definitions::BACKBONE_H_ATOMS +
                                                                        definitions::CB_ATOMS +
                                                                        definitions::SIDECHAIN_ATOMS +
                                                                        definitions::NON_BACKBONE_H_ATOMS)) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->remove_atoms(atom_selection);
          }
     }

     //! Activate atoms that were previously deactivated
     //!
     //! \param atom_selection Which atoms to activate
     void activate_atoms(definitions::AtomSelectionEnum atom_selection = (definitions::BACKBONE_O_ATOMS + 
                                                                          definitions::BACKBONE_H_ATOMS +
                                                                          definitions::CB_ATOMS)) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->activate_atoms(atom_selection);
          }
     }

     //! Deactivate atoms in chain
     //!
     //! \param atom_selection Which atoms to deactivate
     void deactivate_atoms(definitions::AtomSelectionEnum atom_selection = (definitions::BACKBONE_O_ATOMS + 
                                                                            definitions::BACKBONE_H_ATOMS +
                                                                            definitions::CB_ATOMS +
                                                                            definitions::SIDECHAIN_ATOMS +
                                                                            definitions::NON_BACKBONE_H_ATOMS)) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->deactivate_atoms(atom_selection);
          }
     }

     //! Idealize chain - set ideal bond lengths and bond angles, and update positions
     //!
     //! \param atom_selection Which atoms to idealize     
     void idealize(definitions::AtomSelectionEnum atom_selection = definitions::ALL_ATOMS) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->idealize(atom_selection);
          }
          update_positions();
     }

     //! Apply translation to all residues
     //!
     //! \param translation translation vector
     void translate(Vector_3D translation) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->translate(translation);
          }
     }

     //! Apply rotation to all residues
     //!
     //! \param rotation rotation vector
     void rotate(Matrix_3D rotation) {
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->rotate(rotation);
          }
     }

     //! Overload indexing operator
     //!
     //! \param index Residue index
     //! \return Residue reference
     RESIDUE_TYPE &operator[](const int index) const {
          return (*residues.at(index));
     }

     //! Overload () indexing operator
     //! \param index Residue index
     //! \return Residue reference
     RESIDUE_TYPE &operator()(const int index) const {
          return (*residues.at(index));
     }

     //! Overload () operator to directly access atom
     //! \param index Residue index
     //! \param atom_type Type of atom
     //! \return Atom reference
     Atom *operator()(const int index, const definitions::AtomEnum atom_type) const {
          if (index < 0 || index >= (int) size()) {
               return NULL;
          } else {
               RESIDUE_TYPE *res = residues.at(index);
               if (res->has_atom(atom_type)) {
                    return (*residues.at(index))[atom_type];
               } else {
                    return NULL;
               }
          }
     }

     //! Return chaintree for the current chain.
     //! If there is no chaintree associated with the chain, create it
     //!
     //! \return Pointer to chaintree object
     // Note: this function has external requirement not available at
     // this point in the code. It is included at the end of the file
     ChainTree *get_chain_tree();

     //! Remove chaintree from chain.
     // Note: this function has external requirement not available at
     // this point in the code. It is included at the end of the file
     void remove_chain_tree();

     //! Toggle between sidechain states. 
     //! Switches between the different sidechain representations
     //! trying to automatically determine the right move and
     //! leave the chain in a consistent state.
     void toggle_sidechains() {
          for (Iterator it = this->begin(); it != this->end(); ++it) {
               (*it)->toggle_sidechain();
          }
     }

     //! Switches from a full atom sidechain representation to
     //! pseudo sidechains. 
     void toggle_sidechains_to_PS() {
          for (Iterator it = this->begin(); it != this->end(); ++it) {
               (*it)->toggle_sidechain_to_PS();
          }
          // make sure everything is in good order
          this->update_positions();
          this->update_angles();
          this->check_consistency();
     }

     //! Switches from a pseudo-sidechain representation to
     //! a full atom sidechain representation. 
     void toggle_sidechains_to_full_atom() {
          for (Iterator it = this->begin(); it != this->end(); ++it) {
               (*it)->toggle_sidechain_to_full_atom();
          }
          // make sure everything is in good order
          this->update_positions();
          this->update_angles();
          this->check_consistency();
     }


     //! Superimpose this chain upon other chain
     //!
     //! \param other Target chain on which to superimpose onto
     //! \param begin_offset Offset from beginning of chain
     //! \param end_offset Offset from end of chain
     void superimpose_onto_chain(Chain &other, int begin_offset = 0, int end_offset = 0) {

          Matrix_3D rotation_matrix;
          Vector_3D singular_values;
          Vector_3D cm1, cm2;
          strong_assert(this->size() == other.size());
          if (begin_offset == 0 && end_offset == 0) {
               calc_superimpose_rotation_matrix(AtomIterator<Chain<RESIDUE_TYPE>, definitions::CA_ONLY, Vector_3D>(*this),
                                                AtomIterator<Chain<RESIDUE_TYPE>, definitions::CA_ONLY, Vector_3D>(other),
                                                cm1, cm2, rotation_matrix, singular_values);
          } else {
               calc_superimpose_rotation_matrix(AtomIterator<Chain<RESIDUE_TYPE>, definitions::CA_ONLY, Vector_3D>(*this, begin_offset,
                                                                                                                   (this->size() - end_offset)),
                                                AtomIterator<Chain<RESIDUE_TYPE> , definitions::CA_ONLY,Vector_3D>(other, begin_offset,
                                                                                                                   (other.size() - end_offset)),
                                                cm1, cm2, rotation_matrix, singular_values);
          }

          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->translate(-cm1);
               (*it)->rotate(rotation_matrix);
               (*it)->translate(cm2);
          }
     }

     //! Output in PDB format
     //! 
     //! \param superimpose_chain Optionally specify another chain to superimpose onto
     //! \param begin_offset Offset from beginning of chain
     //! \param end_offset Offset from end of chain
     //! \param model_index Model number (for pdb format) 
     //! \param header Optional remark added to header
     //! \param b_factor_string Optional vector string to output as b factors
     //! \return pdb output string
     std::string output_as_pdb(Chain *superimpose_chain = NULL,
                               int begin_offset = 0, int end_offset = 0,
                               const int model_index = 1, 
                               std::string header = "",
                               std::string b_factor_string = "") {
          std::string output = "";
          int chain_number = 1;
          char date[10];
          time_t tt = time(NULL);
          char *transformed_time;
          transformed_time = asctime(localtime(&tt));
          sprintf(date, "%.2s-%.3s-%.2s", &(transformed_time[8]), &(transformed_time[4]), &(transformed_time[22]));
          if (model_index == 1) {
               output += std::string("HEADER    Sampled Structure                          ") + std::string(date) + std::string("         XXXX\n");
          }
          if (header != "") {
               output += std::string("REMARK 0\n");
               output += std::string("REMARK 0 Created by Phaistos\n");
               // output += std::string("REMARK 0 ") + replace(strip(header), "\n", "\nREMARK 0 ");
               //output += std::string("REMARK 0 ") + boost::replace_all_copy(boost::trim_copy(header), "\n", "\nREMARK 0 ");
               output += std::string("REMARK 0 \n");
               std::vector<std::string> header_lines;
               boost::split(header_lines, header, boost::is_any_of("\n"));
               for (unsigned int j = 0; j < header_lines.size(); j++) {
                   if (header_lines[j].find("SHEET") == 0 or header_lines[j].find("HELIX") == 0) {
                       output += header_lines[j] + std::string("\n");
                   } else {
                       output += std::string("REMARK 0 ") + header_lines[j] + std::string("\n");
                   }
               }
               output += std::string("REMARK 0\n");
          }
          output += "MODEL        " + stringify(model_index) + "\n";

          Chain<RESIDUE_TYPE> *tmp_chain = this->clone();
          // Move chain to center of mass if no superimpose chain is specified
          if (superimpose_chain == NULL) {
               // Vector_3D CM = center_of_mass(atomBeginCA(), atomEndCA());
               Vector_3D cm = center_of_mass(AtomIterator<Chain<RESIDUE_TYPE>, definitions::CA_ONLY, Vector_3D> (*this));
               tmp_chain->translate(-cm);
          } else {
               if (superimpose_chain != this)
                    tmp_chain->superimpose_onto_chain(*superimpose_chain, begin_offset, end_offset);
          }

          int counter = 0;
          for (Iterator it = tmp_chain->begin(); it != tmp_chain->end(); ++it) {

               output += (*it)->output_as_pdb(counter, chain_number, &b_factor_string);
               counter += (*it)->size();
          }
          delete tmp_chain;
          return output;
     }

     //! Output to PDB file
     //! 
     //! \param filename Output filename
     //! \param superimpose_chain Optionally specify another chain to superimpose onto
     //! \param begin_offset Offset from beginning of chain
     //! \param end_offset Offset from end of chain
     //! \param model_index Chain id (for pdb format) 
     //! \param header Optional remark added to header
     //! \param b_factor_string Optional vector string to output as b factors
     //! \return pdb output string
     void output_as_pdb_file(std::string filename,
                             Chain *superimpose_chain = NULL,
                             int begin_offset = 0, int end_offset = 0, const int model_index = 1,
                             std::string header = "",
                             std::string b_factor_string = "") {
          std::ofstream pdb_file;
          pdb_file.open(filename.c_str());
          pdb_file << output_as_pdb(superimpose_chain, begin_offset, end_offset, model_index, header, b_factor_string);
          pdb_file.close();
     }

     //! Output to PDB file
     //! 
     //! \param filename Output filename
     //! \param superimpose_chain Optionally specify another chain to superimpose onto
     //! \param header Optional remark added to header
     //! \param b_factor_string Optional vector string to output as b factors
     //! \return pdb output string
     void output_as_pdb_file(std::string filename,
                             Chain *superimpose_chain, std::string header,
                             std::string b_factor_string = "") {
          std::ofstream pdb_file;
          pdb_file.open(filename.c_str());
          pdb_file << output_as_pdb(superimpose_chain, 0, 0, 1, header, b_factor_string);
          pdb_file.close();
     }


     //! Check chain consistency (if debuglevel>0)
     void check_consistency() {

#if DEBUGLEVEL > 0
          for (Iterator it = begin(); it != end(); ++it) {
               (*it)->check_consistency();
          }
          if (chain_tree) {
               chain_tree->check_consistency();
          }
#endif
     }

     //! Get iterator for given position
     //!
     //! \param index Residue index
     //! \return Residue iterator
     Iterator get_iterator_at_index(int index) {
          if (index >= (int) residues.size()) {
               return residues.end();
          } else if (index < 0) {
               return residues.begin();
          } else {
               return residues.begin() + index;
          }
     }

     //! Output to stream
     //!
     //! \param o Output stream
     void output(std::ostream &o) {
          for (Iterator it = begin(); it != end(); ++it) {
               o << **it << "\n";
          }
     }

     //! Overload << operator
     friend std::ostream & operator<<(std::ostream &o, Chain<RESIDUE_TYPE> &c) {
          c.output(o);
          return o;
     }

};

}

// Member functions using external objects that could not
// be included at the top of the file (since Chain was not
// defined there)

#include "chaintree.h"


namespace phaistos {

// Copy constructor
template<typename RESIDUE_TYPE>
Chain<RESIDUE_TYPE>::Chain(const Chain<RESIDUE_TYPE> &other)
     : name(other.name) {
     if (other.chain_tree) {
          this->chain_tree = new typename Chain<RESIDUE_TYPE>::ChainTree(this);
     } else {
          this->chain_tree = NULL;
     }
     time_stamp = other.time_stamp;
}

// Get chaintree pointer
// if there is no chaintree associated with the chain, create it
template<typename RESIDUE_TYPE>
typename Chain<RESIDUE_TYPE>::ChainTree *Chain<RESIDUE_TYPE>::get_chain_tree() {
     if (!chain_tree) {
          chain_tree = new typename Chain<RESIDUE_TYPE>::ChainTree(this);
     }
     return chain_tree;
}

// Remove chaintree
template<typename RESIDUE_TYPE>
void Chain<RESIDUE_TYPE>::remove_chain_tree() {
     delete chain_tree;
     chain_tree = NULL;
}

}

#endif
