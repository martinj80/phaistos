// term_backbone_dbn.h --- BackboneDBN model energy
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


#ifndef TERM_BACKBONE_DBN_H
#define TERM_BACKBONE_DBN_H

#include "../energy/energy_term.h"

namespace phaistos {

//! BackboneDBN energy term. Local energy term of backbone dihedral angles.
template<typename CHAIN_TYPE, typename DBN_TYPE>
class TermBackboneDBN : public EnergyTermCommon<TermBackboneDBN<CHAIN_TYPE,DBN_TYPE>, CHAIN_TYPE> {
private:

     typedef phaistos::EnergyTermCommon<TermBackboneDBN<CHAIN_TYPE,DBN_TYPE>,CHAIN_TYPE> EnergyTermCommon;

     //! BackboneDbn model object
     DBN_TYPE *dbn;

     double energy_value;
     double energy_value_backup;
     
     //! Likelihood of old state minus new state - used for calculating move bias
     double ll_move_bias;

     std::vector<std::vector<double> >angle_cache;
     int start_index;
     int end_index;

public:

     //! Local settings class.     
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:
          
          //! Whether to always force update of all angles in the DBN from the chain
          bool always_full_update;

          //! Whether to update the angles in the DBN from the chain when necessary
          bool enable_dbn_update;

          //! Size of window used (on each side) when calculating the dbn energy.
          //! A good value for the window size is >7 and a negative window size, means that the full bias
          //! will be calculate.
          int window_size;
          
          //! Whether to divide out move bias corresponding to the energy
          bool eliminate_move_bias;

          //! Constructor
          Settings(bool always_full_update=false, 
                   bool enable_dbn_update=true,
                   int window_size=-1,
                   bool eliminate_move_bias = false) 
               : always_full_update(always_full_update), 
                 enable_dbn_update(enable_dbn_update),
                 window_size(window_size),
                 eliminate_move_bias(eliminate_move_bias) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "always-full-update:" << settings.always_full_update << "\n";
               o << "enable-dbn-update:" << settings.enable_dbn_update << "\n";
               o << "window-size:" << settings.window_size << "\n";
               o << "eliminate-move-bias:" << settings.eliminate_move_bias << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings> (settings);
               return o;
          }
     } settings;  //!< Local settings object 

     //! Constructor.
     //! \param chain Molecule chain
     //! \param dbn Dynamic Bayesian Network model
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermBackboneDBN(CHAIN_TYPE *chain, DBN_TYPE *dbn, 
                  const Settings &settings=Settings(),
                  RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "backbone-dbn", settings, random_number_engine), 
            dbn(new DBN_TYPE(*dbn)), 
            energy_value_backup(UNINITIALIZED),
            ll_move_bias(0.0),
            start_index(0), end_index(chain->size()),
            settings(settings) {}


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain     
     TermBackboneDBN(const TermBackboneDBN &other, 
                  RandomNumberEngine *random_number_engine,
                  int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain), 
            // dbn(&other.dbn->get_copy(thread_index)), 
            dbn(new DBN_TYPE(*other.dbn)), 
            energy_value_backup(UNINITIALIZED),
            ll_move_bias(0.0),
            start_index(0), end_index(chain->size()),
            settings(other.settings) {
     }

     //! Destructor
     ~TermBackboneDBN() {
          delete dbn;
     }

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {

          energy_value = 0.0;

          // Use cached version if last move was a sidechain move
          if (!settings.always_full_update &&
              is_initialized(energy_value_backup) &&
              move_info &&
              move_info->move_type == definitions::SIDECHAIN) {

               energy_value = energy_value_backup;

          // Otherwise calculate energy
          } else {
                        
               // Save fixed/unfixed state of angleNode
               bool angle_emission_state = dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->fixed;

               // Default range is the entire chain
               this->start_index = 0;
               this->end_index = this->chain->size();

               // Check whether to update dbn from chain
               if (settings.enable_dbn_update) {

                    // If available, use information about the range of the move
                    if (!settings.always_full_update && move_info && is_initialized(move_info->modified_angles_start)) {
                         start_index = move_info->modified_angles_start;
                         end_index = move_info->modified_angles_end;
                    }

                    // Backup angle values
                    this->angle_cache = this->dbn->template get_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES> (start_index, 
                                                                                                                     end_index);
                    // Insert angle values from chain into dbn
                    std::vector<std::vector<double> > angles_tmp;
                    bool include_omega = true;
                    for (int i = start_index; i < end_index; i++) {
                         std::vector<double> angles = ((*this->chain)[i].get_angles(include_omega));
                         angles_tmp.push_back(angles);
                    }
                    bool fix_emission=false;
                    bool set_observed=false;
                    this->dbn->template set_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES>(angles_tmp, 
                                                                                                start_index,
                                                                                                fix_emission,
                                                                                                set_observed);
               }

               // Set emission status of angle node to true
               bool value = true, set_uninitialized = false;
               dbn->template set_emission_state<typename DBN_TYPE::ALL_ANGLE_NODES>(value, set_uninitialized);

               // Calculate likelihood
               double ll=0.0;
               if (!move_info || !is_initialized(energy_value_backup) || settings.window_size < 0) {
                    ll = dbn->get_log_likelihood();
               } else {
                    int start_index = std::max(0,                   move_info->modified_angles_start - settings.window_size);
                    int end_index   = std::min(this->chain->size(), move_info->modified_angles_end + settings.window_size);
                         
                    bool use_backup = true;
                    double old_ll = dbn->get_log_likelihood(start_index, end_index, use_backup);
                    double new_ll = dbn->get_log_likelihood(start_index, end_index);
                    ll = -energy_value_backup - old_ll + new_ll;
               }

               // Revert to old emission state of angleNode
               value = angle_emission_state; set_uninitialized = false;
               dbn->template set_emission_state<typename DBN_TYPE::ALL_ANGLE_NODES>(value, set_uninitialized);

               // Return the negative log likelihood
               energy_value = -ll;
          }
          if (settings.eliminate_move_bias && is_initialized(energy_value_backup)) {
              this->ll_move_bias = -energy_value_backup+energy_value;
          }
          return energy_value;
     }
     
     //! Get the log bias of the move
     //!
     //! Only used if eliminate-move-bias is set to true. Also, the bias is only
     //! calculated if the move corresponds to the energy.
     //!
     double get_log_bias(MoveInfo *move_info = NULL) {
          MoveInfoDbn<DBN_TYPE> * move_info_dbn = dynamic_cast<MoveInfoDbn<DBN_TYPE> *>(move_info);
          //Is eliminate-move-bias set to true and does the move correspond to the energy 
          if (settings.eliminate_move_bias && move_info_dbn != NULL) {
               return this->ll_move_bias;
          }
          return .0;
     }

     void reject () {
          if (!angle_cache.empty()) {
               bool fix_emission=false;
               bool set_observed=false;
               dbn->template set_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES>(angle_cache, 
                                                                                     start_index,
                                                                                     fix_emission,
                                                                                     set_observed);
               angle_cache.clear();
          }
     }

     void accept () {
          energy_value_backup = energy_value;
          angle_cache.clear();

          dbn->accept(start_index,
                      end_index);
     }

};

}

#endif
