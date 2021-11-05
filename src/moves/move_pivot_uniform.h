// move_pivot_uniform.h --- Move: Simple pivot move sampling changes to phi,psi angles from a uniform distribution
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


#ifndef MOVE_PIVOT_UNIFORM_H
#define MOVE_PIVOT_UNIFORM_H

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"
#include "move.h"

namespace phaistos {

//! Simple pivot move sampling changes to phi,psi angles from a uniform distribution
template <typename CHAIN_TYPE>
class MovePivotUniform: public MoveCommon<MovePivotUniform<CHAIN_TYPE>,
                                          CHAIN_TYPE> {

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MovePivotUniform<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;

public:

     //! Local Settings class.
     const class Settings: public Move<ChainFB>::Settings {
     public:

          //! Whether only to resample a single dof in each iteration (selected randomly)
          bool single_dof_only;

          //! Whether to skip prolines phi angles (modifiying the proline phi angle introduces an improper torsion change)
          bool skip_proline_phi;          

          //! Maximum change in angle value
          double max_delta;
        
          //! Constructor - set default move_length to 1
          Settings(bool single_dof_only=true,
                   bool skip_proline_phi=true,
                   double max_delta=M_PI)
               : Move<ChainFB>::Settings(1,1),
                 single_dof_only(single_dof_only),
                 skip_proline_phi(skip_proline_phi),
                 max_delta(max_delta) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               o << "single-dof-only:" << settings.single_dof_only << "\n";               
               o << "skip-proline-phi:" << settings.skip_proline_phi << "\n";               
               o << "max-delta:" << settings.max_delta << "\n";               
               return o;
          }                    
     } settings;    //!< Local settings object
     

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MovePivotUniform(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                      RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "pivot-uniform", settings, random_number_engine), 
            settings(settings) {
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     MovePivotUniform(const MovePivotUniform &other)
          : MoveCommon(other),
            settings(other.settings) {
     }
     
     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MovePivotUniform(const MovePivotUniform &other, 
                      RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            settings(other.settings) {
     }

     //! Apply move
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether move was succesful
     bool apply(int start_index=-1, int end_index=-1) {
          
          // Call base class apply method
          Move<CHAIN_TYPE>::apply(start_index, end_index);

	  if (this->chain_backup) {
	       delete this->chain_backup;
	       this->chain_backup=NULL;
	  }

	  // If start and end index are not set, choose random values
	  // within the length range specified by the constructor
	  if (start_index == end_index) {
               this->random_move_range(&start_index, &end_index);
	  }

	  // Remember which angles are modified
          std::pair<int,int> modified_angles = std::make_pair(start_index, end_index);

	  // Find direction of position update
	  // update_positions takes range of modified angles - end_index
	  // itself was not modified
	  int direction = this->chain->find_shortest_direction(start_index, end_index-1);

	  // Set range of modified indices
          std::pair<int,int> modified_positions;
	  if (direction > 0) {
               modified_positions = std::make_pair(start_index, this->chain->size());
	  } else {
               modified_positions = std::make_pair(0, end_index);
	  }
          this->move_info = (new MoveInfo(definitions::NON_LOCAL))->add_info(modified_angles,
                                                                             modified_positions);

          // Make a copy of the modified portion of the chain. This is
          // used if the move is rejected
          this->chain_backup = new CHAIN_TYPE(*this->chain,
                                              modified_positions.first,
                                              modified_positions.second);
          
          for (int i=modified_angles.first; i<modified_angles.second; ++i) {

               typename CHAIN_TYPE::Residue &res = (*this->chain)[i];

               std::vector<double> dihedrals = (*this->chain)[i].get_angles();

               if (settings.single_dof_only) {
                    int index = 0;
                    if (((res.residue_type == definitions::PRO) && settings.skip_proline_phi) ||
                        (this->random_generator_uniform_01() < 0.5)) {
                         index = 1;
                    }
                    dihedrals[index] += this->random_generator_uniform_01()*2*settings.max_delta - settings.max_delta;
                    dihedrals[index] = fmod(dihedrals[index]+3*M_PI,(2*M_PI))-M_PI;
               } else {
                    if (!((res.residue_type == definitions::PRO) && settings.skip_proline_phi)) {
                         dihedrals[0] += this->random_generator_uniform_01()*2*settings.max_delta - settings.max_delta;
                         dihedrals[0] = fmod(dihedrals[0]+3*M_PI,(2*M_PI))-M_PI;
                    }
                    dihedrals[1] += this->random_generator_uniform_01()*2*settings.max_delta - settings.max_delta;
                    dihedrals[1] = fmod(dihedrals[1]+3*M_PI,(2*M_PI))-M_PI;
               }
               (*this->chain)[i].set_angles(dihedrals);
          }

          this->chain->update_positions(start_index, end_index-1, direction);

	  return this->move_info->success;
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
          
          AtomIterator<CHAIN_TYPE,definitions::ALL> it_current(*this->chain, 
                                                               this->move_info->modified_positions[0].first,
                                                               this->move_info->modified_positions[0].second);
          AtomIterator<CHAIN_TYPE,definitions::ALL> it_backup(*this->chain_backup);

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
