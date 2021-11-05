// move_fixed_structure.h --- Initializing move - changes chain to fixed reference structure
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


#ifndef MOVE_FIXED_STRUCTURE_H
#define MOVE_FIXED_STRUCTURE_H

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"

namespace phaistos {

//! Initializing move: changes chain to fixed reference structure
template <typename CHAIN_TYPE>
class MoveFixedStructure: public MoveCommon<MoveFixedStructure<CHAIN_TYPE>, CHAIN_TYPE> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveFixedStructure<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;     

     //! Reference chain to which current chain is modified at each call of apply()
     CHAIN_TYPE reference_chain;

     //! update direction
     int direction;
     
public:

     //! Use base class Settings class
     typedef typename Move<CHAIN_TYPE>::Settings Settings;

     //! Local settings object
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param reference_chain Reference chain to which current chain is modified at each call of apply()
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param settings Local Settings object
     MoveFixedStructure(CHAIN_TYPE *chain, 
                        CHAIN_TYPE &reference_chain,
                        RandomNumberEngine *random_number_engine = &random_global,
                        const Settings &settings = Settings(-1,-1)): 
          MoveCommon(chain, "fixed-structure", 
                     settings, 
                     random_number_engine),
          reference_chain(reference_chain),
          settings(settings) {}

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveFixedStructure(const MoveFixedStructure<CHAIN_TYPE> &other):
          MoveCommon(other),
          reference_chain(other.reference_chain),
          direction(other.direction),
          settings(other.settings) {}

     
     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveFixedStructure(const MoveFixedStructure<CHAIN_TYPE> &other, 
                        RandomNumberEngine *random_number_engine, 
                        int thread_index,
                        CHAIN_TYPE *chain):
          MoveCommon(other,random_number_engine, chain),
          reference_chain(other.reference_chain),
          direction(other.direction),
          settings(other.settings) {}


     //! Apply move
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether move was succesful
     bool apply(int start_index=-1, int end_index=-1) {

          // Call base class apply method
          Move<CHAIN_TYPE>::apply(start_index, end_index);

	  if (this->chain_backup) {
	       delete this->chain_backup;
	  }

	  // Set start and end index if necessary
	  if(end_index < 0 || end_index > this->chain->size())
	       end_index = this->chain->size();
	  
	  if(start_index < 0 || start_index > this->chain->size())
	       start_index = 0;	  
	  
	  // Find direction of position update
	  this->direction = this->chain->find_shortest_direction(start_index, end_index-1);

	  // Set range of modified indices
	  if (direction > 0) {
               this->move_info = (new MoveInfo())->add_info(std::make_pair(start_index, this->chain->size()),
                                                            std::make_pair(start_index, this->chain->size()));
	  } else {
               this->move_info = (new MoveInfo())->add_info(std::make_pair(0, end_index),
                                                            std::make_pair(0, end_index));
	  }	  
	  
	  // Make a copy of the modified portion of the chain. This is
	  // used if the move is rejected
	  this->chain_backup = new CHAIN_TYPE(*this->chain,
                                              this->move_info->modified_positions[0].first,
					      this->move_info->modified_positions[0].second);


          for(ResidueIterator<CHAIN_TYPE> res_it(*(this->chain)); !res_it.end(); ++res_it) {
               typename CHAIN_TYPE::Residue *res = &*res_it;
               int res_index = res->index;
               Residue *res_reference = &reference_chain[res_index];
               if (res->residue_type != res_reference->residue_type) {
                    std::cerr << "Error - move_fixed_structure: reference chain is not compatible with chain. Residue mismatch at index " << res_index << ".\n";
                    assert(false);
               }
               for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(*res_it); !it.end(); ++it) {
                    if (!res_reference->has_atom(it->atom_type)) {
                         std::cerr << "Error - move_fixed_structure: reference chain is not compatible with chain. Atom " << it->atom_type << " missing in residue " << res_index << "in reference chain.\n";
                         assert(false);
                    }
                    Atom *atom_reference = (*res_reference)[it->atom_type];
                    it->set_angle(atom_reference->get_angle());
                    it->set_dihedral(atom_reference->get_dihedral());
                    it->position = atom_reference->position;               
               }
          }


          // Initialize bond lengths
	  // for (deprecated::AtomIteratorAll it1=chain_begin; it1!=chain_end; ++it1) {
	  for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(*this->chain, start_index, end_index); !it.end(); ++it) {
               it->init_bond_length();
          }          
	  
	  return true;
     }

     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {
	  return 0.0;
     }
     
     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();

	  delete this->chain_backup;
	  this->chain_backup = NULL;
     }

     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();

	  assert(this->chain_backup != NULL);

          AtomIterator<CHAIN_TYPE,
                       definitions::ALL> it_current(*this->chain, 
                                                    this->move_info->modified_positions[0].first,
                                                    this->move_info->modified_positions[0].second);
          AtomIterator<CHAIN_TYPE,
                       definitions::ALL> it_backup(*this->chain_backup);

          for (; !it_current.end() && !it_backup.end(); ++it_current,++it_backup) {
               it_current->set_angle(it_backup->get_angle());
               it_current->set_dihedral(it_backup->get_dihedral());
               it_current->position = it_backup->position;
          }


	  delete this->chain_backup;
	  this->chain_backup = NULL;
	  
     }
};

}

#endif
