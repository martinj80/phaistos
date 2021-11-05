// move_sidechain_dbn_multi.h --- Move: Resample sidechain angles using DBN - multi version
// Copyright (C) 2008-2009 Tim Harder
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


#ifndef MOVE_SIDECHAIN_DBN_MULTI_H
#define MOVE_SIDECHAIN_DBN_MULTI_H

#include "protein/chain_fb.h"
#include "moves/move.h"
#include "mocapy.h"
#include "models/basilisk/basilisk_dbn.h"
#include "models/compas/compas_dbn.h"
#include <typeinfo>

namespace phaistos {

//! Sidechain DBN Move - multi
//!
//! This move allows to modify multiple, disjunct side chains in a
//! a single move. This can be particularly interesting since it
//! effectively allows side chains to swap places in a densely packed
//! protein core.
template<typename DBN_TYPE>
class MoveSidechainDBNMulti: public MoveCommon<MoveSidechainDBNMulti<DBN_TYPE>, ChainFB> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveSidechainDBNMulti<DBN_TYPE> , ChainFB> MoveCommon;

     //! Dynamic Bayesian Network object (can be both a Compas or a Basilisk DBN)     
     DBN_TYPE *dbn;

     //! Vector of start indices
     std::vector<int> start_residues;

     //! Vector of end indices
     std::vector<int> end_residues;

     //! Backup of the current angles to allow a rollback
     std::vector<std::vector<std::vector<double> > > chi_backups;

     //! Likelihood of old state - before making move     
     double pre_ll;

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

          //! Number of moves
          int move_count;

          //! Constructor
          Settings(std::string model_filename=set_model_filename_default((DBN_TYPE *)NULL),
                   std::string mocapy_dbn_dir="../data/mocapy-dbns",
                   bool implicit_energy = true, 
                   bool ignore_bb = false, 
                   int move_count = 3)
               : Move<ChainFB>::Settings(1,1), // modify only a single sidechain at a time
                 model_filename(model_filename),
                 mocapy_dbn_dir(mocapy_dbn_dir),
                 implicit_energy(implicit_energy), 
                 ignore_bb(ignore_bb), 
                 move_count(move_count) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "model-filename:" << settings.model_filename << "\n";
               o << "mocapy-dbn-dir:" << settings.mocapy_dbn_dir << "\n";
               o << "implicit-energy:" << settings.implicit_energy << "\n";
               o << "ignore-bb:" << settings.ignore_bb << "\n";
               o << "move-count:" << settings.move_count << "\n";
               o << static_cast<const typename Move<ChainFB>::Settings &> (settings);
               return o;
          }
     } settings;    //!< Local settings object


     //! Constructor.
     //! \param chain Molecule chain
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveSidechainDBNMulti(ChainFB *chain, DBN_TYPE *dbn, 
                    const Settings &settings = Settings(), 
                    RandomNumberEngine *random_number_engine = &random_global) 
          : MoveCommon(chain, "sc-"+dbn->name+"-multi", settings, random_number_engine), 
            dbn(dbn),
            pre_ll(-1.0),
            settings(settings) {
     }

     //! Copy constructor
     //!
     //! \param other instance of a MoveScDbn to create a copy from
     MoveSidechainDBNMulti(const MoveSidechainDBNMulti<DBN_TYPE> &other) :
          MoveCommon(other), dbn(other.dbn), pre_ll(other.pre_ll), settings(other.settings) {
     }

     //! Copy constructor
     //!
     //! Copy constructor to be used in a multithread environment,
     //! allows setting the correct random number engine as well
     //! as choosing the correct instance of the moveSc itself
     //!
     //! \param other instance of a MoveScDbn to create a copy from
     //! \param random_number_engine Random number enginge collection, allowing thread safe runs
     //! \param thread_index Which thread the copy will run in
     //! \param chain pointer to the current chain object
     MoveSidechainDBNMulti(const MoveSidechainDBNMulti<DBN_TYPE> &other, 
                    RandomNumberEngine *random_number_engine, 
                    int thread_index, ChainFB *chain) 
          : MoveCommon(other, random_number_engine, chain),
            dbn(&other.dbn->get_copy(thread_index)), 
            pre_ll(other.pre_ll), 
            settings(other.settings) {
     }

     //! Get/Set the range
     //!
     //! this Method allows to get a number of random ranges
     //! across the protein chain. It also sets the values
     //! as the start/end residues for the next move .. allowing
     //! to get the residues before applying the move, in order
     //! to get some local energies or similar.
     //!
     //! \param start_index vector of residue indeces of the first modified residue
     //! \param end_index vector of residue indeces of the last residue to be modified
     void get_set_range(std::vector<int> &start_index, std::vector<int> &end_index) {
          int length = this->random_move_length();
          if (length > this->chain->size())
               length = this->chain->size();

          int start = 0;
          int end = 0;

          if (this->start_residues.size() > 0 || this->end_residues.size() > 0) {
               std::cerr << "Residues vectors not empty .. aborting\n";
               return;
          }

          bool doublette = false;
          for (int i = 0; i < settings.move_count; i++) {
               this->random_move_range(length, 0, this->chain->size(), &start, &end);
               doublette = false;

               for (int j = 0; j < (int) start_index.size(); j++) {
                    if (start_index[j] == start) {
                         // if the position is already in the
                         // vector, redo that step.
                         i--;
                         doublette = true;
                         break;
                    }
               }

               if (doublette)
                    continue;

               start_index.push_back(start);
               end_index.push_back(end);
               // saving those for later
               this->start_residues.push_back(start);
               this->end_residues.push_back(end);
          }
     }

     //! Apply move
     //!
     //! This method does the actual modification of the conformation.
     //! Depending on the setting, a backbone dependent or independent
     //! sample will be drawn from the according model.
     //!
     //! \param start_index residue index of the first modified residue
     //! \param end_index residue index of the last residue to be modified
     bool apply(int start_index = -1, int end_index = -1) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Call base class apply method
          Move<ChainFB>::apply(start_index, end_index);

          this->pre_ll = 0.;

          // we need to register the modified residues
          // for the chaintree to update properly
          // NOTE: THIS SHOULD BE CHANGED - IT IS NOW POSSIBLE
          // TO REGISTER A VECTOR OF (START,END) PAIRS
          this->move_info = (new MoveInfo(SIDECHAIN))->add_info(std::make_pair(0, this->chain->size()),
                                                                std::make_pair(0, this->chain->size()));

          // if we haven't selected any specific move range yet
          // lets do so now
          if ((int) this->start_residues.size() < settings.move_count || (int) this->end_residues.size() < settings.move_count) {
               this->start_residues.clear();
               this->end_residues.clear();
               std::vector<int> tmp_start_index;
               std::vector<int> tmp_end_index;
               this->get_set_range(tmp_start_index, tmp_end_index);
          }
          this->pre_ll = 0;

          for (int i = 0; i < settings.move_count; i++) {

               start_index = this->start_residues[i];
               end_index = this->end_residues[i];

               // get all the chi angles we need to get back
               // to the old state afterwards ..
               this->chi_backups.push_back(this->chain->get_sidechain_dof_values(this->dbn->get_atom_selection(), 
                                                                                 start_index, end_index));

               // if the bias shall be correct, we need
               // to calculate the previous probablity
               if (!settings.implicit_energy) {
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
               //
               for (int i = start_index; i < end_index; i++) {
                    ResidueFB *t = &(*this->chain)[i];
                    ResidueEnum res = (*t).residue_type;
                    if (settings.ignore_bb) {
                         new_angles.push_back(this->dbn->get_angle_sequence_for_residue_with_bb(res));
                    } else {
                         new_angles.push_back(this->dbn->get_angle_sequence_for_residue_with_bb(res, t->get_phi(), t->get_psi()));
                    }
               }

               if (new_angles.size()) {
                    this->chain->set_sidechain_dof_values(new_angles, this->dbn->get_atom_selection(), start_index, end_index);
               }
          }

          finalize();

          // check that the move produced consistent
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
          for (int j = 0; j < settings.move_count; j++) {
               for (int i = this->start_residues[j]; i < this->end_residues[j]; i++) {
                    (*(this->chain))[i].update_positions_non_backbone();
               }
          }
     }

     // int getStartResidue(int movecount) {
     //      return this->start_residues[movecount];
     // }
     // int getEndResidue(int movecount) {
     //      return this->end_residues[movecount];
     // }
     // std::vector<int> getStartResidues() {
     //      return this->start_residues;
     // }
     // std::vector<int> getEndResidues() {
     //      return this->end_residues;
     // }

     //! Accept last move
     //!
     //! Invalidate the backup cache.
     void accept() {
          //
          // Call base class accept method
          Move<ChainFB>::accept();
          // free the backup again
          this->chi_backups.clear();
          this->start_residues.clear();
          this->end_residues.clear();
     }

     //! Reject last move
     //!
     //! Rollback to the previous conformation, since the move was
     //! rejected
     void reject() {
          // Call base class accept method
          Move<ChainFB>::reject();
          // we have to do that backwards, because otherwise
          // it may introduce errors if a residue is present
          // more than once
          if (not this->chi_backups.size())
               return;

          for (int j = settings.move_count - 1; j >= 0; j--) {
               this->chain->set_sidechain_dof_values(this->chi_backups[j], this->dbn->get_atom_selection(), 
                                                     this->start_residues[j], this->end_residues[j]);
          }
          finalize();

          // free the backup again
          this->chi_backups.clear();
          this->start_residues.clear();
          this->end_residues.clear();
     }

     //! get the log bias of the move
     //!
     //! This method returns the exact bias introduced by the proposed move.
     //! In most use case it is actually a desired behaviour to implicitly
     //! including the move bias into the sampling to implicitly also include
     //! the probability distribution encoded in the models into the sampling.
     //! In some cases however, this is not wanted or the probablility should
     //! be included explicitly  .. in such case it is necessary to explicitly
     //! remove the bias.
     double get_log_bias() {
          // If we are using the dbn as an implicit energy, return 0.0
          if (settings.implicit_energy) {
               return 0.0;
               // otherwise, divide out the proposal
          } else {
               double ll = 0;
               for (int i = 0; i < settings.move_count; i++) {
                    int start_index = this->start_residues[i];
                    int end_index = this->end_residues[i];
                    for (int k = start_index; k < end_index; k++) {
                         ResidueFB *t = &(*this->chain)[k];
                         if (settings.ignore_bb) {
                              ll += this->dbn->calc_ll((*t).residue_type, 
                                                       UNINITIALIZED, UNINITIALIZED, 
                                                       (*t).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                         } else {
                              ll += this->dbn->calc_ll((*t).residue_type, 
                                                       (*t).get_phi(), (*t).get_psi(), 
                                                       (*t).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                         }
                    }
               }
               return (this->pre_ll - ll);
          }
     }

     //! Check the resulting conformation for proline ring breaks
     //!
     //! In order to ensure a proper ring closure for newly sampled
     //! proline residues, we perform a check after proposing new
     //! conformations. Here the ring closure is checked and the move
     //! rejected if an invalid structure was proposed.
     bool check_proline_rings() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          bool success = true;

          for (int j = 0; j < settings.move_count; j++) {
               for (int i = this->start_residues[j]; i < this->end_residues[j]; i++) {
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
          }
          return success;
     }

};

}

#endif
