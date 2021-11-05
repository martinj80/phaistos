// move_dbn.h --- Move: Resample angles using DBN
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


#ifndef MOVE_DBN_H
#define MOVE_DBN_H

#include "move.h"

#include "backbone_dbn.h"
#include "protein/chain_fb.h"
#include "protein/chain_ca.h"

namespace phaistos {

//! Local enum defining different types of sampling
enum MoveBackboneDBNResampleMode{RESAMPLE_ALL=0, RESAMPLE_HIDDEN_ONLY, RESAMPLE_ANGLES_ONLY, RESAMPLE_MODE_SIZE};

//! Names for each of the enum values
static const std::string MoveDbnResampleModeNames[] = {"resample-all", "resample-hidden-only", "resample-angles-only"};

//! Input MoveBackboneDBNResampleMode from stream
inline std::istream &operator>>(std::istream &input, MoveBackboneDBNResampleMode &rm) {
     std::string raw_string;
     input >> raw_string;

     for (unsigned int i=0; i<RESAMPLE_MODE_SIZE; ++i) {
          if (raw_string == MoveDbnResampleModeNames[i]) {
               rm = MoveBackboneDBNResampleMode(i);
          }
     }
     return input;
}

//! Output MoveBackboneDBNResampleMode to stream
inline std::ostream &operator<<(std::ostream &o, const MoveBackboneDBNResampleMode &rm) {
     o << MoveDbnResampleModeNames[static_cast<unsigned int>(rm)];
     return o;
}


//! Update backbone using dihedral angles sampled from a dynamic Bayesian network.
template <typename CHAIN_TYPE, typename DBN_TYPE>
class MoveBackboneDBN: public MoveCommon<MoveBackboneDBN<CHAIN_TYPE,DBN_TYPE>,
                                                         CHAIN_TYPE> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveBackboneDBN<CHAIN_TYPE, DBN_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Pointer to dynamic Bayesian Network model
     DBN_TYPE *dbn;

     //! update direction
     int direction;

public:

     //! Local Settings class.
     const class Settings: public Move<CHAIN_TYPE>::Settings {
     public:

          //! Specifies what is resampled when calling apply (RESAMPLE_ALL, RESAMPLE_HIDDEN_ONLY, RESAMPLE_ANGLES_ONLY)
          MoveBackboneDBNResampleMode resample_mode;

          //! If false, the DBN bias is divided out
          bool implicit_energy;

          //! Size of window used (on each side) when bringing the dbn back to consistency 
          //! (negative number => resample entire sequence)
          int dbn_consistency_window_size;

          //! Size of window used (on each side) when calculating the dbn bias.
          //! The bias for the move is P(X)/P(X'), where X is the old angle sequence and X' is the
          //! new. If only the angles in the interval [i, j] are changed, the bias can be approximated
          //! by P(X[i-w..j+w])/P(X'[i-w..j+w]), where w in the window size.
          //! A good value for the window size is >7 and a negative window size, means that the full bias
          //! should be calculate.
          int dbn_bias_window_size;

          //! Constructor
          Settings(MoveBackboneDBNResampleMode resample_mode=RESAMPLE_ALL, 
                   bool implicit_energy=true,
                   int dbn_consistency_window_size=10,
                   int dbn_bias_window_size=10)
               : Move<CHAIN_TYPE>::Settings(1,1),
                 resample_mode(resample_mode),
                 implicit_energy(implicit_energy),
                 dbn_consistency_window_size(dbn_consistency_window_size),
                 dbn_bias_window_size(dbn_bias_window_size) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "resample-mode:" << settings.resample_mode << "\n";
               o << "implicit-energy:" << settings.implicit_energy << "\n";
               o << "dbn-consistency-window-size:" << settings.dbn_consistency_window_size << "\n";
               o << "dbn-bias-window-size:" << settings.dbn_bias_window_size << "\n";
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                              
     } settings;    //!< Local settings object

     //! Constructor
     //! \param chain Molecule chain
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveBackboneDBN(CHAIN_TYPE *chain, DBN_TYPE *dbn,
              const Settings &settings=Settings(),
              RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "dbn", settings, random_number_engine),
            dbn(dbn),
            settings(settings) {

          // When using a DBN move, the angle node is flagged as fixed
          // and is then unfixed every time a move is made. This ensures
          // that when the dbn is resynchronized after a change has been 
          // made to the angles outside the DBN, the hidden nodes will be 
          // resampled based on the new angle values
          dbn->template set_emission_state<typename DBN_TYPE::ALL_ANGLE_NODES>(true);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveBackboneDBN(const MoveBackboneDBN<CHAIN_TYPE,DBN_TYPE> &other)
          : MoveCommon(other),
            dbn(other.dbn),
            settings(other.settings) {
     }
     
     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveBackboneDBN(const MoveBackboneDBN<CHAIN_TYPE,DBN_TYPE> &other, 
              RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            dbn(&other.dbn->get_copy(thread_index)),
            settings(other.settings) {

          // When using a DBN move, the angle node is flagged as fixed
          // and is then unfixed every time a move is made. This ensures
          // that when the dbn is resynchronized after a change has been 
          // made to the angles outside the DBN, the hidden nodes will be 
          // resampled based on the new angle values
          dbn->template set_emission_state<typename DBN_TYPE::ALL_ANGLE_NODES>(true);
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
	  this->direction = this->chain->find_shortest_direction(start_index, end_index-1);

	  // Set range of modified indices
          std::pair<int,int> modified_positions;
	  if (direction > 0) {
               modified_positions = std::make_pair(start_index, this->chain->size());
	  } else {
               modified_positions = std::make_pair(0, end_index);
	  }
          this->move_info = (new MoveInfoDbn<DBN_TYPE>(dbn, definitions::NON_LOCAL))->
               add_info(modified_angles,
                        modified_positions);

          // Make sure that DBN is consistent
          dbn->enforce_consistency(settings.dbn_consistency_window_size);

	  if(settings.resample_mode == RESAMPLE_HIDDEN_ONLY) {

	       // If only hidden states are updated, we don't need to backup the chain
	       this->chain_backup = NULL;

	       // Resample
	       dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->sample(start_index, end_index);

	  } else {

	       // Make a copy of the modified portion of the chain. This is
	       // used if the move is rejected
	       this->chain_backup = new CHAIN_TYPE(*this->chain,
                                                   this->move_info->modified_positions[0].first,
                                                   this->move_info->modified_positions[0].second);

	       // Save fixed/unfixed state of angleNode
	       // bool angle_emission_state = dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->fixed;
               
               // Set emission status of angle node to false (so that angles can be modified)
               bool value = false, set_uninitialized = false;
               dbn->template set_emission_state<typename DBN_TYPE::ALL_ANGLE_NODES>(value, set_uninitialized);

	       // Resample
	       if (settings.resample_mode == RESAMPLE_ANGLES_ONLY) {
		    dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->sample(start_index, end_index);
		    // dbn->template sample<typename DBN_TYPE::ALL_ANGLE_NODES>(start_index, end_index);
	       } else {
		    dbn->sample(start_index, end_index);
	       }

               // Reset angles emission state
               value = true;
               dbn->template set_emission_state<typename DBN_TYPE::ALL_ANGLE_NODES>(value, set_uninitialized);

	       // Update positions in chain
	       std::vector<std::vector<double> > angles = 
                    dbn->template get_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES>(start_index,
                                                                                          end_index);

	       for (int i=start_index; i<end_index; i++) {
		    (*this->chain)[i].set_angles(angles[i-start_index]);
	       }

	       this->chain->update_positions_backbone(start_index, end_index-1, direction);

	       finalize(this->chain);
	  }

	  return this->move_info->success;
     }


     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {
	  double bias = 0.0;
	  if (settings.implicit_energy)
	       bias = 0.0;
	  else {

	       if (settings.resample_mode == RESAMPLE_HIDDEN_ONLY) {

		    // Calculate likelihood before move
		    bool use_backup = true;
		    double old_LL = dbn->get_log_likelihood_conditional(this->move_info->modified_angles[0].first,
                                                                        this->move_info->modified_angles[0].second, use_backup);
		    double new_LL = dbn->get_log_likelihood_conditional(this->move_info->modified_angles[0].first,
                                                                        this->move_info->modified_angles[0].second);

		    bias = -(new_LL - old_LL);

	       } else {

		    // Calculate likelihood before the move and after the move
		    bool use_backup = true;
		    double old_LL, new_LL;

		    if (settings.resample_mode == RESAMPLE_ANGLES_ONLY) {
		         old_LL = dbn->get_log_likelihood_conditional(this->move_info->modified_angles[0].first, 
                                                                      this->move_info->modified_angles[0].second, use_backup);
	                 new_LL = dbn->get_log_likelihood_conditional(this->move_info->modified_angles[0].first, 
                                                                      this->move_info->modified_angles[0].second);
	            } else {
	                 int start, end;

	                 if (settings.dbn_bias_window_size >= 0) {
	                      start = std::max(0,                   this->move_info->modified_angles[0].first  - settings.dbn_bias_window_size);
	                      end   = std::min(this->chain->size(), this->move_info->modified_angles[0].second + settings.dbn_bias_window_size);
	                 } else {
	                      start = -1;
	                      end   = -1;
	                 }

	                 old_LL = dbn->get_log_likelihood(start, end, use_backup);
	                 new_LL = dbn->get_log_likelihood(start, end);
	            }
                    //std::cout << new_LL << " "  <<  old_LL << " " << -(new_LL - old_LL) << "\n";
                    //std::cout << "MOVE: " << new_LL << "\n" << *dbn << "\n";

		    bias = -(new_LL - old_LL);
	       }
	  }
	  return bias;
     }

     //! Finalize move - calculates new positions for CB and O atoms
     void finalize(CHAIN_TYPE *chain) {
	  for (int i=this->move_info->modified_positions[0].first; i<this->move_info->modified_positions[0].second; i++) {
	       (*chain)[i].update_positions_non_backbone();
	  }
     }

     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();
          
	  delete this->chain_backup;
	  this->chain_backup = NULL;
	  
	  this->dbn->accept(this->move_info->modified_angles_start,
                            this->move_info->modified_angles_end);
     }


     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();
          
	  // We don't have to reconstruct the old chain if we only resampled hidden values
	  if (settings.resample_mode != RESAMPLE_HIDDEN_ONLY) {
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
               
	  }
	  delete this->chain_backup;
	  this->chain_backup = NULL;
	  this->dbn->reject(this->move_info->modified_angles_start,
                            this->move_info->modified_angles_end);
     }

     //! Specifies how move reacts to the knowledge that another moves has been executed
     //! In the case of a move_DBN, any modification to backbone angles must
     //! be registered with the DBN, so it remains consistent
     void notify(MoveInfo *move_info) {
          if (!move_info || move_info->modified_angles.empty())
               return;
          if (move_info->move_type != definitions::SIDECHAIN) {
               MoveInfoDbn<DBN_TYPE> *move_info_dbn = dynamic_cast<MoveInfoDbn<DBN_TYPE> *>(move_info);
               if (!move_info_dbn) {

                    std::vector<std::vector<double> > angles;
                    bool include_omega = true;
                    for (int i=move_info->modified_angles_start; i<move_info->modified_angles_end; ++i) {
                         angles.push_back((*this->chain)[i].get_angles(include_omega));
                    }

                    dbn->template set_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES>(angles, move_info->modified_angles_start);

                    this->dbn->register_inconsistency(move_info->modified_angles_start,
                                                      move_info->modified_angles_end);
               }
          }
     }
};

}

#endif
