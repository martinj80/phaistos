// move_sidechain_uniform.h --- Move: Simple sidechain moves sampling chi angles from a uniform distribution
// Copyright (C) 2006-2011 Wouter Boomsma
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


#ifndef MOVE_SIDECHAIN_UNIFORM_H
#define MOVE_SIDECHAIN_UNIFORM_H

namespace phaistos {

//! Simple pivot move sampling phi,psi angles from a uniform distribution
template <typename CHAIN_TYPE>
class MoveSidechainUniform: public MoveCommon<MoveSidechainUniform<CHAIN_TYPE>,
                                              CHAIN_TYPE> {

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveSidechainUniform<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Last accepted chi angles (for a range of residues)
     std::vector<std::vector<double> > chi_angle_backup;

public:

     //! Local Settings class.
     const class Settings: public Move<ChainFB>::Settings {
     public:

          //! Whether only to resample a single dof in each iteration (selected randomly)
          bool single_dof_only;

          //! Whether to skip prolines (prolines introduce a change in bond length which must be taken into account by the forcefield)
          bool skip_proline;
        
          //! Constructor - set default move_length to 1
          Settings(bool single_dof_only=true,
                   bool skip_proline=true)
               : Move<ChainFB>::Settings(1,1),
                 single_dof_only(single_dof_only),
                 skip_proline(skip_proline) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               o << "single-dof-only:" << settings.single_dof_only << "\n";               
               o << "skip-proline:" << settings.skip_proline << "\n";               
               return o;
          }                    
     } settings;    //!< Local settings object

     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveSidechainUniform(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                   RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "sc-uniform", settings, random_number_engine), 
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveSidechainUniform(const MoveSidechainUniform &other)
          : MoveCommon(other),
            chi_angle_backup(other.chi_angle_backup),
            settings(other.settings) {
     }

     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveSidechainUniform(const MoveSidechainUniform &other,
                   RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            chi_angle_backup(other.chi_angle_backup),
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
               std::vector<double> chi_angles = res->get_sidechain_dof_values(SIDECHAIN_ATOMS);
               chi_angle_backup.push_back(chi_angles);

               if (chi_angles.size() == 0)
                    continue;

               if (settings.single_dof_only) {

                    // Sample index in range
                    boost::variate_generator<RandomNumberEngine&,
                         boost::uniform_int<> >
                         random_generator_uniform_int(*this->random_number_engine,
                                                      boost::uniform_int<>(0,
                                                                           chi_angles.size()-1));
                    int random_index = random_generator_uniform_int();
                    chi_angles[random_index] = 2*M_PI*this->random_generator_uniform_01() - M_PI;
                    
               } else {
                    for (unsigned int j=0; j<chi_angles.size(); ++j) {
                         chi_angles[j] = 2*M_PI*this->random_generator_uniform_01() - M_PI;
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
          return 0.0;
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
