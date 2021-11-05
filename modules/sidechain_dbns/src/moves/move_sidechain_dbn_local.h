// move_sidechain_dbn_local.h --- Move: Resample sidechain angles using DBN - local version
// Copyright (C) 2008-2009 Tim Harder, Wouter Boomsma
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


#ifndef MOVE_SIDECHAIN_DBN_LOCAL_H
#define MOVE_SIDECHAIN_DBN_LOCAL_H

#include "protein/chain_fb.h"
#include "moves/move.h"
#include "mocapy.h"
#include "models/basilisk/basilisk_dbn.h"
#include "models/compas/compas_dbn.h"
#include <typeinfo>

namespace phaistos {

template<typename DBN_TYPE>

//! Sidechain DBN Move - local
//!
//! This move uses the side chain dbn model to propose small
//! changes to the side chain conformation. This is accomplished
//! by sampling the hidden node states from the current
//! angles and subsequently sampling a set of new angles from the
//! this fixed hidden node sequence. Most samples will remain in
//! the same rotameric states as the current conformation.
class MoveSidechainDBNLocal: public MoveCommon<MoveSidechainDBNLocal<DBN_TYPE> , ChainFB> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveSidechainDBNLocal<DBN_TYPE> , ChainFB> MoveCommon;

     //! Dynamic Bayesian Network object (can be both a Compas or a Basilisk DBN)
     DBN_TYPE *dbn;

     //! Backup of the current angles to allow a rollback     
     std::vector<std::vector<double> > chi_backup;
     
     //! Likelihood of old state - before making move
     double pre_ll;

     //! Reference
     ChainFB *reference;

public:

     //! Local Settings class.
     const class Settings: public Move<ChainFB>::Settings {

          //! Default model path for basilisk
          static std::string set_model_filename_default(BasiliskDBN *dbn=NULL) {return "basilisk.dbn";}

          //! Default model path for compas
          static std::string set_model_filename_default(CompasDBN *dbn=NULL) {return "compas.dbn";}

     public:
          
          //! Path to model file
          std::string model_filename;

          //! Directory containing Mocapy models          
          std::string mocapy_dbn_dir;

          //! The default settings for this move is to run in biased
          //! mode, where the move contributes with an implicit energy.
          //! Turn on this settings to divide out the bias an regain an unbiased simulation.
          bool implicit_energy;

          //! Do not use backbone information          
          bool ignore_bb;
          
          //! Whether the angles are extracted from a fixed reference chain, rather than the current chain
          bool fixed_reference_chain;


          //! Constructor
          Settings(std::string model_filename=set_model_filename_default((DBN_TYPE *)NULL),
                   std::string mocapy_dbn_dir="../data/mocapy-dbns",
                   bool implicit_energy = true, 
                   bool ignore_bb = false, 
                   bool fixed_reference_chain=false)
               : Move<ChainFB>::Settings(1,1), // modify only a single sidechain at a time
                 model_filename(model_filename),
                 mocapy_dbn_dir(mocapy_dbn_dir),
                 implicit_energy(implicit_energy), 
                 ignore_bb(ignore_bb), 
                 fixed_reference_chain(fixed_reference_chain) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "model-filename:" << settings.model_filename << "\n";
               o << "mocapy-dbn-dir:" << settings.mocapy_dbn_dir << "\n";
               o << "implicit-energy:" << settings.implicit_energy << "\n";
               o << "ignore-bb:" << settings.ignore_bb << "\n";
               o << "fixed-reference-chain:" << settings.fixed_reference_chain << "\n";
               o << static_cast<const typename Move<ChainFB>::Settings &> (settings);
               return o;
          }
     } settings;    //!< Local settings object


     //! Constructor.
     //! \param chain Molecule chain
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveSidechainDBNLocal(ChainFB *chain, DBN_TYPE *dbn, const Settings &settings = Settings(), 
                    RandomNumberEngine *random_number_engine = &random_global) 
          : MoveCommon(chain, "sc-"+dbn->name+"-local", settings, random_number_engine),
            dbn(dbn),
            pre_ll(-1.0),
            reference(NULL),
            settings(settings) {

          if (settings.fixed_reference_chain) {
               this->reference = new ChainFB(*(this->chain));
          } else {
               this->reference = this->chain;
          }
     }

     //! Copy constructor
     //!
     //! \param other instance of a MoveScDbn to create a copy from
     MoveSidechainDBNLocal(const MoveSidechainDBNLocal<DBN_TYPE> &other) :
          MoveCommon(other), dbn(other.dbn), pre_ll(other.pre_ll), 
          reference(other.reference), 
          settings(other.settings) {
     }

     //! Copy constructor
     //!
     //! Copy constructor to be used in a multithread environment,
     //! allows setting the correct random number engine as well
     //! as choosing the correct instance of the moveSc itself
     //!
     //! \param other instance of a MoveScDbn to create a copy from
     //! \param random_number_engine Random number enginge collection, allowing thread safe runs
     //! \param thread_index
     //! \param chain Molecule chain
     MoveSidechainDBNLocal(const MoveSidechainDBNLocal<DBN_TYPE> &other, 
                    RandomNumberEngine *random_number_engine, 
                    int thread_index, ChainFB *chain) :
          MoveCommon(other, random_number_engine, chain), 
          dbn(&other.dbn->get_copy(thread_index)), 
          pre_ll(other.pre_ll), 
          reference(other.reference), 
          settings(other.settings) {
     }

     //! Apply move
     //!
     //! This method does the actual modification of the conformation.
     //! Depending on the setting, a backbone dependent or independent
     //! sample will be drawn from the according model after fixing the
     //! the hidden node sequence to the one sampled conditioned on the
     //! the current sidechain angles.
     //!
     //! \param start_index residue index of the first modified residue
     //! \param end_index residue index of the last residue to be modified
     bool apply(int start_index = -1, int end_index = -1) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Call base class apply method
          Move<ChainFB>::apply(start_index, end_index);

          // If start and end index are not set, choose random values
          // within the length range specified by the constructor
          if (start_index == end_index) {
               this->random_move_range(&start_index, &end_index);
          }

          // Remember which residues we modified (for moveInfo object)
          this->move_info = (new MoveInfo(SIDECHAIN))->add_info(std::make_pair(start_index, end_index),
                                                                std::make_pair(start_index, end_index));

          // get all the chi angles we need to get back
          // to the old state afterwards ..
          this->chi_backup = this->chain->get_sidechain_dof_values(this->dbn->get_atom_selection(), start_index, end_index);

          // if the bias shall be correct, we need
          // to calculate the previous probablity
          if (!settings.implicit_energy) {
               this->pre_ll = 0;
               for (int k = start_index; k < end_index; k++) {
                    ResidueFB *t = &(*this->chain)[k];
                    if (settings.ignore_bb) {
                         this->pre_ll += this->dbn->calc_ll((*t).residue_type, 
                                                            UNINITIALIZED, UNINITIALIZED, 
                                                            (*t).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    } else {
                         this->pre_ll += this->dbn->calc_ll((*t).residue_type, 
                                                            (*t).get_phi(), (*t).get_psi(), 
                                                            (*t).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    }
               }
          }

          std::vector<std::vector<double> > new_angles;
          std::vector<double> chis_old;
          for (int i = start_index; i < end_index; i++) {
               ResidueFB *t = &(*this->reference)[i];
               ResidueEnum res = (*t).residue_type;
               chis_old = (*t).get_sidechain_dof_values(this->dbn->get_atom_selection());
               if (settings.ignore_bb) {
                    new_angles.push_back(this->dbn->get_local_angle_sequence_for_residue(res, chis_old));
               } else {
                    new_angles.push_back(this->dbn->get_local_angle_sequence_for_residue(res, chis_old, t->get_phi(), t->get_psi()));
               }
          }

          // put them in the chain
          if (new_angles.size())
               this->chain->set_sidechain_dof_values(new_angles, this->dbn->get_atom_selection(), start_index, end_index);

          // Update positions
          finalize();

          // check that the move produced consisted
          // atom positions
          this->move_info->success = check_proline_rings();

          return this->move_info->success;
     }

     //! finalizing the move
     //!
     //! After the move, which modified the angles of the sidechain,
     //! it is necessary to update the structure to reflect the modifications
     //! in the conformation as well.
     void finalize() {
          for (int i = this->move_info->modified_angles[0].first; i < this->move_info->modified_angles[0].second; i++) {
               (*(this->chain))[i].update_positions_non_backbone();
          }
     }

     //! Accept last move
     //!
     //! Invalidate the backup cache.
     void accept() {

          // Call base class accept method
          Move<ChainFB>::accept();

          // free the backup again
          this->chi_backup.clear();
     }

     //! Reject last move
     //!
     //! Rollback to the previous conformation, since the move was
     //! rejected
     void reject() {

          // Call base class accept method
          Move<ChainFB>::reject();

          if (this->chi_backup.size())
               this->chain->set_sidechain_dof_values(this->chi_backup, this->dbn->get_atom_selection(),
                                                     this->move_info->modified_angles_start,
                                                     this->move_info->modified_angles_end);
          // Update positions
          finalize();

          // free the backup again
          this->chi_backup.clear();
     }


     //! Get the log bias of the move
     //!
     //! This method returns the exact bias introduced by the proposed move.
     //! In most use case it is actually a desired behaviour to implicitly
     //! including the move bias into the sampling to implicitly also include
     //! the probability distribution encoded in the models into the sampling.
     //! In some cases however, this is not wanted or the probablility should
     //! be included explicitly. In such case it is necessary to explicitly
     //! remove the bias.
     double get_log_bias() {

          // If we are using the dbn as an implicit energy, return 0.0
          if (settings.implicit_energy) {
               return 0.0;
          // otherwise, divide out the proposal
          } else {
               double bias = 0.;

               std::pair<int, int> tmp = this->move_info->modified_angles.back();

               for (int k = tmp.first; k < tmp.second; k++) {
                    ResidueFB *t = &(*this->chain)[k];
                    if (settings.ignore_bb) {
                         bias += this->dbn->calc_ll((*t).residue_type, -
                                                    UNINITIALIZED, UNINITIALIZED, 
                                                    (*t).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    } else {
                         bias += this->dbn->calc_ll((*t).residue_type, 
                                                    (*t).get_phi(), (*t).get_psi(), 
                                                    (*t).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    }
               }
               return (this->pre_ll - bias);
          }
     }

     //! Check the resulting conformation for proline ring breaks
     //!
     //! In order to ensure a proper ring closure for newly sampled
     //! proline residues, we perform a check after proposing new
     //! conformations. Here the ring closure is check and the move
     //! rejected if an invalid structure was proposed.
     bool check_proline_rings() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          bool success = true;

          for (int i = this->move_info->modified_angles[0].first; i < this->move_info->modified_angles[0].second; i++) {
               ResidueEnum res = (*(this->chain))[i].residue_type;
               if (res == PRO) {
                    Vector_3D cd = (*(this->chain))[i][CD]->position;
                    Vector_3D n = (*(this->chain))[i][N]->position;
                    double dist = (cd - n).norm();
                    if (dist < 1.45 || dist > 1.49) {
                         success = false;
                    }
               }
          }
          return success;
     }

};

}

#endif
