// clash_simple.h --- Fast clash detection using chaintree
// Copyright (C) 2008-2011 Wouter Boomsma
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


#ifndef CLASH_FAST_H
#define CLASH_FAST_H

#include "../energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace phaistos {

//! Clash detection energy term. 
//! Uses the chaintree structure to quickly detect clashing atoms.
template <typename CHAIN_TYPE>
class TermClashFast: public EnergyTermCommon<TermClashFast<CHAIN_TYPE>, CHAIN_TYPE> {
protected:

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Chaintree iterator settings - described in detail in chaintree.h
     typename chaintree::PairIterator<CHAIN_TYPE,NodeType,NodeType>::Settings it_settings;

private:
     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermClashFast<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;     

     // Matrix of cutoff distances
     std::vector<std::vector<double> > max_interaction_distance_matrix;

     // Lookup table for indices into max_interaction_distance_matrix
     std::vector<std::vector<int> > distance_matrix_index_lookup;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! Specifies whether energy should consider only pairs modified by the last
          //! move. If a clashfree structure is maintained at all times, this is sufficient
          //! to detect all clashes
          bool only_modified_pairs;

          //! Specifies whether energy should work in clash/non-clash mode and return
          //! infinity/0 (true) or count all the clashes and return the number of clashes (false)
          bool boolean_mode;

          //! Minimum distance along chain (measured in number of
          //! residues) before pair is taken into account by the energy
          int minimum_residue_distance;

          //! Clash distance constant - hydrogens
          double clash_distance_H;

          //! Clash distance constant - Nitrogen-Oxygen
          double clash_distance_NO;

          //! Clash distance constant - Pseudo-atoms
          double clash_distance_PS;

          //! Clash distance constant - sulfur
          double clash_distance_SG;

          //! Clash distance constant - Any other pair
          double clash_distance_any_pair;

          //! Constructor. Defines default values for settings object.
          Settings(bool only_modified_pairs=true,
                   bool boolean_mode=true,
                   int minimum_residue_distance=2,
                   double clash_distance_H=1.5,
                   double clash_distance_NO=2.3,
                   double clash_distance_any_pair=2.3,
                   double clash_distance_PS=1.9,
                   double clash_distance_SG=1.8)
               : only_modified_pairs(only_modified_pairs),
                 boolean_mode(boolean_mode),
                 minimum_residue_distance(minimum_residue_distance),
                 clash_distance_H(clash_distance_H),
                 clash_distance_NO(clash_distance_NO),
                 clash_distance_PS(clash_distance_PS),
                 clash_distance_SG(clash_distance_SG),
                 clash_distance_any_pair(clash_distance_any_pair) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "only-modified-pairs:" << settings.only_modified_pairs << "\n";
               o << "boolean-mode:" << settings.boolean_mode << "\n";
               o << "minimum-residue-distance:" << settings.minimum_residue_distance << "\n";
               o << "clash-distance-h:" << settings.clash_distance_H << "\n";
               o << "clash-distance-no:" << settings.clash_distance_NO << "\n";
               o << "clash-distance-any-pair:" << settings.clash_distance_any_pair << "\n";
               o << "clash-distance-ps:" << settings.clash_distance_PS << "\n";
               o << "clash-distance-sg:" << settings.clash_distance_SG << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }                    
     } settings;  //!< Local settings object 


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermClashFast(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "clash-fast", settings, random_number_engine),
            settings(settings) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          enum MatrixIndex {INDEX_N, INDEX_O, INDEX_XH, INDEX_PS, INDEX_SG, INDEX_XX, MATRIX_INDEX_SIZE};

          std::vector<AtomEnum> selected_atoms(MATRIX_INDEX_SIZE);
          selected_atoms[INDEX_N] = N;
          selected_atoms[INDEX_O] = O;
          selected_atoms[INDEX_XH] = XH;
          selected_atoms[INDEX_PS] = PS;
          selected_atoms[INDEX_SG] = SG;
          selected_atoms[INDEX_XX] = XX;

          for (int i=0; i<chain->size(); ++i) {
               distance_matrix_index_lookup.push_back(std::vector<int>((*chain)[i].size(), -1));
               for (int j=0; j<(*chain)[i].size(); ++j) {
                    Atom *atom = (*chain)[i][j];

                    for (unsigned int k=0; k<selected_atoms.size() && distance_matrix_index_lookup[i][j]==-1; ++k) {
                         AtomEnum atom_type = selected_atoms[k];
                         if (atom->atom_type == atom_type) {
                              distance_matrix_index_lookup[i][j] = k;
                         } else if (is_atom_wildcard(atom_type)) {
                              if (atom_type == XX) {
                                   distance_matrix_index_lookup[i][j] = k;                                        
                              } else {
                                   for (int l=0; l<atom_type_wildcards_size; l++) {
                                        if (atom->atom_type == atom_type_wildcards[atom_type-XS][l]) {
                                             distance_matrix_index_lookup[i][j] = k;                                        
                                        }
                                   }
                              }
                         } 
                    }
               }
          }

          max_interaction_distance_matrix = std::vector<std::vector<double> >(selected_atoms.size(),
                                                                              std::vector<double>(selected_atoms.size(),
                                                                                                  settings.clash_distance_any_pair));
          max_interaction_distance_matrix[INDEX_N][INDEX_O] = settings.clash_distance_NO;
          max_interaction_distance_matrix[INDEX_O][INDEX_N] = settings.clash_distance_NO;
          
          max_interaction_distance_matrix[INDEX_SG][INDEX_SG] = settings.clash_distance_SG;
          
          // PS settings
          for (unsigned int j=0; j<max_interaction_distance_matrix[3].size(); j++) {
               max_interaction_distance_matrix[INDEX_PS][j] = settings.clash_distance_PS;
               max_interaction_distance_matrix[j][INDEX_PS] = settings.clash_distance_PS;
          }     
          // Hydrogen settings
          for (unsigned int j=0; j<max_interaction_distance_matrix[2].size(); j++) {
               max_interaction_distance_matrix[INDEX_XH][j] = settings.clash_distance_H;
               max_interaction_distance_matrix[j][INDEX_XH] = settings.clash_distance_H;
          }

          double distance_cutoff = 0.0;
          for (unsigned int i=0; i<max_interaction_distance_matrix.size(); i++) {
               for (unsigned int j=0; j<max_interaction_distance_matrix[i].size(); j++) {
                    if (max_interaction_distance_matrix[i][j] > distance_cutoff) {
                         distance_cutoff = max_interaction_distance_matrix[i][j];
                    }
               }
          }

          // Square all entries in max_interaction_distance_matrix
          for (unsigned int i=0; i<max_interaction_distance_matrix.size(); i++) {
               for (unsigned int j=0; j<max_interaction_distance_matrix[i].size(); j++) {
                    max_interaction_distance_matrix[i][j] *= max_interaction_distance_matrix[i][j];
               }
          }          

          it_settings = typename chaintree::PairIterator<CHAIN_TYPE,NodeType,NodeType>::Settings(distance_cutoff,
                                                                                                 settings.only_modified_pairs,
                                                                                                 settings.minimum_residue_distance);
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermClashFast(const TermClashFast &other, 
                   RandomNumberEngine *random_number_engine,
                   int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            it_settings(other.it_settings),
            max_interaction_distance_matrix(other.max_interaction_distance_matrix),
            distance_matrix_index_lookup(other.distance_matrix_index_lookup),
            settings(other.settings) {
     }     


     //! Check whether pair is within clashing distance
     //! \param atom1 First atom
     //! \param atom2 Second atom
     bool detect_clash(const Atom *atom1, const Atom *atom2) {
          double distance_squared = (atom1->position - 
                                     atom2->position).norm_squared();
          int index1 = distance_matrix_index_lookup[atom1->residue->index][atom1->index];
          int index2 = distance_matrix_index_lookup[atom2->residue->index][atom2->index];
          double distance_cutoff_squared = max_interaction_distance_matrix[index1][index2];

          return (distance_squared < distance_cutoff_squared);
     }

     //! Check whether pair is within clashing distance
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param distance_squared_pointer Allows the user to provide a variable in which to write the calculated distance
     bool detect_clash(const Atom *atom1, const Atom *atom2, double *distance_squared_pointer) {
          double distance_squared = (atom1->position - 
                                     atom2->position).norm_squared();
          *distance_squared_pointer = distance_squared;
          int index1 = distance_matrix_index_lookup[atom1->residue->index][atom1->index];
          int index2 = distance_matrix_index_lookup[atom2->residue->index][atom2->index];
          double distance_cutoff_squared = max_interaction_distance_matrix[index1][index2];

          return (distance_squared < distance_cutoff_squared);
     }

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     inline double evaluate(MoveInfo *move_info=NULL) {

          if (settings.boolean_mode) {

               // Iterator over node pairs
               for (typename chaintree::PairIterator<CHAIN_TYPE, NodeType, NodeType> it(*this->chain, it_settings); !it.end(); ++it) {

                    // Check if nodes are identical
                    if (it->first == it->second) {

                         // Iterate over all atoms in first node
                         for (unsigned int i=0; i<it->first->size(); ++i) {
                              Atom *atom1 = (*it->first)[i];
                              int index1 = distance_matrix_index_lookup[atom1->residue->index][atom1->index];

                              const std::vector<double> &distance_matrix_row = max_interaction_distance_matrix[index1];

                              // Iterate over all atoms in second node
                              for (unsigned int j=i+1; j<it->second->size(); ++j) {
                                   Atom *atom2 = (*it->second)[j];
                                   if (chain_distance<CHAIN_TYPE>(atom1, atom2)>1) { //ADDED BY SOLSSON JUL 28 2011
                                        int index2 = distance_matrix_index_lookup[atom2->residue->index][atom2->index];
                                        double distance_cutoff_squared = distance_matrix_row[index2];
                                        double distance_squared = (atom1->position - 
                                                                   atom2->position).norm_squared();
                                        // if (distance_squared < distance_cutoff*distance_cutoff) {
                                        if (distance_squared < distance_cutoff_squared) {
                                        // if (detect_clash(atom1, atom2)) {
                                             return  std::numeric_limits<double>::infinity();
                                        }
                                   }
                              }
                         }
                    } else {
                         // Iterate over all atoms in first node
                         for (unsigned int i=0; i<it->first->size(); ++i) {
                              Atom *atom1 = (*it->first)[i];
                              int index1 = distance_matrix_index_lookup[atom1->residue->index][atom1->index];

                              const std::vector<double> &distance_matrix_row = max_interaction_distance_matrix[index1];                              

                              // Iterate over all atoms in second node
                              for (unsigned int j=0; j<it->second->size(); ++j) {
                                   Atom *atom2 = (*it->second)[j];
                                   int index2 = distance_matrix_index_lookup[atom2->residue->index][atom2->index];
                                   double distance_cutoff_squared = distance_matrix_row[index2];
                                   double distance_squared = (atom1->position - 
                                                              atom2->position).norm_squared();
                                   if (distance_squared < distance_cutoff_squared) {
                                   // if (detect_clash(atom1, atom2)) {
                                        return  std::numeric_limits<double>::infinity();
                                   }
                              }
                         }
                    }
               }               
               return 0.0;

          } else {

               int clash_counter = 0;
               // Iterator over node pairs
               for (typename chaintree::PairIterator<CHAIN_TYPE, NodeType, NodeType> it(*this->chain, it_settings); !it.end(); ++it) {
                    
                    // Check if nodes are identical
                    if (it->first == it->second) {

                         // Iterate over all atoms in first node
                         for (unsigned int i=0; i<it->first->size(); ++i) {

                              // Iterate over all atoms in second node
                              for (unsigned int j=i+1; j<it->second->size(); ++j) {
                                   Atom *atom1 = (*it->first)[i];
                                   Atom *atom2 = (*it->second)[j];
                                        if (detect_clash(atom1, atom2)) {
                                             ++clash_counter;
                                        }
                                   
                              }
                         }
                    } else {
                         // Iterate over all atoms in first node
                         for (unsigned int i=0; i<it->first->size(); ++i) {

                              // Iterate over all atoms in second node
                              for (unsigned int j=0; j<it->second->size(); ++j) {
                                   Atom *atom1 = (*it->first)[i];
                                   Atom *atom2 = (*it->second)[j];
                                   if (detect_clash(atom1, atom2)) {
                                        ++clash_counter;
                                   }
                              }
                         }
                    }
               }               

               return (double)clash_counter;
          }
     }
};

}

#endif
