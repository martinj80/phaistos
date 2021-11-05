// move_sidechain_rotamer.h --- Move: Resample chi-angles using rotamer library
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


#ifndef MOVE_SIDECHAIN_ROTAMER_H
#define MOVE_SIDECHAIN_ROTAMER_H

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"
#include "move.h"

namespace phaistos {

//! Resample chi-angles using rotamer library
template <typename CHAIN_TYPE, typename ROTAMER_LIB_TYPE>
class MoveSidechainRotamer: public MoveCommon<MoveSidechainRotamer<CHAIN_TYPE,
                                                                   ROTAMER_LIB_TYPE>,
                                      CHAIN_TYPE> {
private:
     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveSidechainRotamer<CHAIN_TYPE,
                                                         ROTAMER_LIB_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Rotamer library object
     ROTAMER_LIB_TYPE rotamer_library;

     //! Last accepted chi angles (for a range of residues)
     std::vector<std::vector<double> > chi_angle_backup;

     //! Current rotamer index (for a range of residues)
     std::vector<int> rotamer_index;

     //! Last rotamer index used (for a range of residues)
     std::vector<int> rotamer_index_prev;

public:
     
     //! Local Settings class.
     const class Settings: public Move<ChainFB>::Settings {
     public:
          
          //! Specifies whether bias of move should be included
          //! when get_log_bias is called. If set to true, the move will act
          //! as an implicit energy term in the simulation.
          bool implicit_energy;
          
          //! Makes it possible to flatten/sharpen the rotamer distributions
          //! by scaling sigma
          double sigma_scale_factor;

          //! Frequency of moves where rotamer state is resampled
          double rotamer_state_resample_frequency;

          //! Whether hydrogen chi angles should be resampled uniformly
          bool sample_hydrogen_chis;

          //! Whether to skip prolines (prolines introduce a change in bond length which must be taken into account by the forcefield)
          bool skip_proline;

          //! Constructor
          Settings(bool implicit_energy=true,
                   double sigma_scale_factor=1.0,
                   double rotamer_state_resample_frequency=1.0,
                   bool sample_hydrogen_chis=true,
                   bool skip_proline=true)
               : Move<ChainFB>::Settings(1,1),
                 implicit_energy(implicit_energy),
                 sigma_scale_factor(sigma_scale_factor),
                 rotamer_state_resample_frequency(rotamer_state_resample_frequency),
                 sample_hydrogen_chis(sample_hydrogen_chis),
                 skip_proline(skip_proline) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "implicit-energy:" << settings.implicit_energy << "\n";
               o << "sigma-scale-factor:" << settings.sigma_scale_factor << "\n";
               o << "rotamer-state-resample-frequency:" << settings.rotamer_state_resample_frequency << "\n";
               o << "sample-hydrogen-chis:" << settings.sample_hydrogen_chis << "\n";
               o << "skip-proline:" << settings.skip_proline << "\n";
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                    
     } settings;    //!< Local settings object
     

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveSidechainRotamer(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                   RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "sc-rotamer", settings, random_number_engine), 
            rotamer_library(settings.sigma_scale_factor, random_number_engine),
            rotamer_index(chain->size(),-1), rotamer_index_prev(chain->size(),-1),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveSidechainRotamer(const MoveSidechainRotamer &other)
          : MoveCommon(other),
            rotamer_library(other.rotamer_library),
            chi_angle_backup(other.chi_angle_backup),
            rotamer_index(other.rotamer_index), 
            rotamer_index_prev(other.rotamer_index_prev),
            settings(other.settings) {

     }


     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveSidechainRotamer(const MoveSidechainRotamer &other,
                   RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            rotamer_library(other.rotamer_library, random_number_engine),
            chi_angle_backup(other.chi_angle_backup),
            rotamer_index(other.rotamer_index), 
            rotamer_index_prev(other.rotamer_index_prev),
            settings(other.settings) {
     }


     //! Apply move
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether move was succesful
     bool apply(int start_index=-1, int end_index=-1) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Call base class apply method
          Move<CHAIN_TYPE>::apply(start_index, end_index);
                    
	  // If start and end index are not set, choose random values
	  // within the length range specified by the constructor
	  if (start_index == end_index) {
               this->random_move_range(&start_index, &end_index);
          }

          // Check whether there are any dofs in region
          bool dofs_found = false;
          for (int i=start_index; i<end_index; i++) {
               typename CHAIN_TYPE::Residue &res = (*this->chain)[i];
               if ((res.residue_type != GLY) &&
                   !(settings.skip_proline && res.residue_type == PRO)) {
                    dofs_found = true;
                    break;
               }
          }

	  // Remember which angles are modified
          this->move_info = (new MoveInfo(SIDECHAIN))->add_info(std::make_pair(start_index, end_index),
                                                                std::make_pair(start_index, end_index));

          // If there are no degrees of freedom here, mark the move as failed
          if (!dofs_found) {
               return (this->move_info->success = false);
          }
          
          // Clear backup
          chi_angle_backup.clear();

          for (int i=start_index; i<end_index; i++) {

               // Get residue
               typename CHAIN_TYPE::Residue *res = &(*this->chain)[i];

               // Save old values
               chi_angle_backup.push_back(res->get_sidechain_dof_values(SIDECHAIN_ATOMS));

               // Save old rotamer index
               rotamer_index_prev[i] = rotamer_index[i];

               std::vector<double> chi_angles;
               if ((rotamer_index[i] == -1) || 
                   (this->random_generator_uniform_01() < settings.rotamer_state_resample_frequency)) {
                    // Sample chi angles from rotamer library
                    chi_angles = rotamer_library.sample(res->residue_type, 
                                                      &rotamer_index[i]);
               } else {
                    // Sample chi angles from rotamer library
                    chi_angles = rotamer_library.sample(res->residue_type, 
                                                      rotamer_index[i]);
               }

               // Optionally resample hydrogen chi angles uniformly
               if (settings.sample_hydrogen_chis) {
                    while(chi_angles.size() < res->chi_atoms.size()) {
                         chi_angles.push_back(2*M_PI*this->random_generator_uniform_01() - M_PI);
                    }
               }

               // Set values in residue
               res->set_sidechain_dof_values(chi_angles, SIDECHAIN_ATOMS);

               // Update positions
               (*this->chain)[i].update_positions_non_backbone();
          }

          return this->move_info->success;
     }

     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {
	  if (settings.implicit_energy) {
               return 0.0;
          } else {
               double old_LL = 0.0;
               double new_LL = 0.0;
               for (int i=this->move_info->modified_angles_start; i<this->move_info->modified_angles_end; i++) {

                    // Get residue
                    typename CHAIN_TYPE::Residue *res = &(*this->chain)[i];
                    
                    old_LL += rotamer_library.get_log_likelihood(res->residue_type,
                                                                chi_angle_backup[i-this->move_info->modified_angles_start]);
                    new_LL += rotamer_library.get_log_likelihood(res->residue_type,
                                                                res->get_sidechain_dof_values(definitions::SIDECHAIN_ATOMS));
               }
               return -(new_LL - old_LL);
          }
     }

     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();
          
          // Clear saved chi angles
          chi_angle_backup.clear();
     }


     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();

          // If modified_angles is empty, do nothing
          if (!this->move_info->success) {
               return;
          }
                    
          // Restore saved chi angles
          for (int i=this->move_info->modified_angles_start; i<this->move_info->modified_angles_end; i++) {

               // Revert to old rotamer state
               rotamer_index[i] = rotamer_index_prev[i];
               
               // Get residue
               typename CHAIN_TYPE::Residue *res = &(*this->chain)[i];

               // Set values in residue
               res->set_sidechain_dof_values(chi_angle_backup[i-this->move_info->modified_angles_start],
                                          definitions::SIDECHAIN_ATOMS);

               // Update positions
               res->update_positions_non_backbone();
          }
     }
     
};
     
}

#endif
