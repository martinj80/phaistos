// sidechain_dbn.h  --- explicit energy term for the sidechain dbn models
// Copyright (C) 2010 Tim Harder
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


#ifndef SIDECHAIN_CACHED_H
#define SIDECHAIN_CACHED_H

#include "energy/energy_term.h"
#include "models/basilisk/basilisk_dbn.h"
#include "models/compas/compas_dbn.h"
#include "moves/move_sidechain_dbn.h"

namespace phaistos {

//! Basilisk/Compas energy term
//!
//! Calculating the correct energy of all the sidechain rotamers
//! in the chain.
template<typename CHAIN_TYPE, typename DBN_TYPE>
class TermSidechainDBN: public EnergyTermCommon<TermSidechainDBN<CHAIN_TYPE, DBN_TYPE> , CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermSidechainDBN<CHAIN_TYPE, DBN_TYPE> , CHAIN_TYPE> EnergyTermCommon;

     //! Bayesian network model
     DBN_TYPE *dbn;

     //! Cache values
     std::vector<double> cache;

     //! Previous cache values (for rejections)
     std::vector<double> roll_back_cache;

     //! Total likelihood value
     double ll_cache;

     //! Total likelihood value (for rejections)
     double ll_old;
     
     //! Likelihood of old state minus new state - used for calculating move bias
     double ll_move_bias;

     //! Defines name for Basilisk model
     static std::string get_name(BasiliskDBN *dbn) {
          return "basilisk";
     }

     //! Defines name for Compas model
     static std::string get_name(CompasDBN *dbn) {
          return "compas";
     }
     
     //! Defines type for Basilisk model
     static int get_type(BasiliskDBN *dbn) {
          return definitions::BASILISK;
     }

     //! Defines type for Compas model
     static int get_type(CompasDBN *dbn) {
          return definitions::COMPAS;
     }

public:
     //! Local settings class.
     //! The settings object allows to  generate a general settings
     //! object that can be used to set all parameters correctly.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
          //@{
          //! Internal functions selecting the default model filename
          static std::string set_model_filename_default() {return set_model_filename_default_inner((DBN_TYPE *)NULL);}
          static std::string set_model_filename_default_inner(BasiliskDBN *dbn=NULL) {return "basilisk.dbn";}
          static std::string set_model_filename_default_inner(CompasDBN *dbn=NULL) {return "compas.dbn";}          
          //@}
     public:

          //! Path to file containing model parameters
          std::string model_filename;

          //! Directory containing mocapy files
          std::string mocapy_dbn_dir;

          //! Whether to use the cache
          bool use_cache;

          //! Whether to ignore the protein backbone
          bool ignore_bb;
          
          //! Whether to divide out move bias corresponding to the energy
          bool eliminate_move_bias;
          
          //! Constructor
          Settings(std::string model_filename=set_model_filename_default(),
                   std::string mocapy_dbn_dir="../data/mocapy-dbns",
                   bool use_cache = true, 
                   bool ignore_bb = false,
                   bool eliminate_move_bias = false) 
               : model_filename(model_filename),
                 mocapy_dbn_dir(mocapy_dbn_dir),
                 use_cache(use_cache), 
                 ignore_bb(ignore_bb),
                 eliminate_move_bias(eliminate_move_bias) {
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "model-filename:" << settings.model_filename << "\n";
               o << "mocapy-dbn-dir:" << settings.mocapy_dbn_dir << "\n";
               o << "use-cache:" << settings.use_cache << "\n";
               o << "ignore-bb:" << settings.ignore_bb << "\n";
               o << "eliminate-move-bias:" << settings.eliminate_move_bias << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings> (settings);
               return o;
          }
     } settings;

     //! Constructor
     //!
     //! \param chain Chain pointer
     //! \param dbn Sidechain DBN pointer
     //! \param settings a SidechainsEnergy settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermSidechainDBN(CHAIN_TYPE *chain, DBN_TYPE *dbn, 
                         const Settings &settings = Settings(),
                         RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, get_name(dbn), settings, random_number_engine), 
            settings(settings) {

          this->dbn = dbn;
          this->ll_old = 0;
          this->ll_move_bias = 0;

          this->cache.resize(chain->size(), 0.);
          this->roll_back_cache.resize(chain->size(), 0.);
     }

     //! Copy constructor
     //!
     //! Create a copy from an existing object
     //!
     //! \param other Sidechain energy term instance
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index id of the current thread (required to select the correct random number engine)
     //! \param chain chain pointer
     TermSidechainDBN(const TermSidechainDBN &other, 
                         RandomNumberEngine *random_number_engine,
                         int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain), 
            dbn(&other.dbn->get_copy(thread_index)), 
            ll_old(other.ll_old), 
            settings(other.settings) {

          this->ll_old = 0;
          // empty and resize our vectors
          this->cache.clear();
          this->roll_back_cache.clear();
          this->cache.resize(chain->size(), 0.);
          this->roll_back_cache.resize(chain->size(), 0.);
     }

     //! Evaluate energy term.
     //!
     //! This method is going to be called for every evaluation of the energy
     //! term. This specific instance is cached for performance optimization.
     //! The caching can be disabled via the settings object.
     //!
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info = NULL) {
          double ll = 0.;

          if (settings.use_cache && settings.ignore_bb && move_info != NULL && 
              move_info->move_type != definitions::SIDECHAIN && this->ll_old != 0) {
               // if possible return a cached value for backbone move
               // (works only in backbone independent mode!!)
               return this->ll_old;
          }

          this->ll_cache = this->ll_old;

          if (not settings.use_cache || move_info == NULL || this->ll_old == 0) {
               // run a full update over all residues
               // take that chance to rebuild the cache
               this->cache.clear();

               // run over all residues of the chain and sum the individual
               // contributions .. including a check whether to use the
               // backbone dependent or independent variant.
               for (ResidueIterator<CHAIN_TYPE> res1(*(this->chain)); !res1.end(); ++res1) {
                    double tmp = 0;
                    if (settings.ignore_bb) {
                         tmp = this->dbn->calc_ll((*res1).residue_type, UNINITIALIZED, UNINITIALIZED, 
                                                  (*res1).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    } else {
                         tmp = this->dbn->calc_ll((*res1).residue_type, (*res1).get_phi(), (*res1).get_psi(), 
                                                  (*res1).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    }
                    this->roll_back_cache[res1->index] = this->cache[res1->index];
                    this->cache[res1->index] = tmp;
                    ll += tmp;
               }

          } else {
               //  run a local update only (requires detailed move information)
               int start = move_info->modified_angles_start;
               int end = move_info->modified_angles_end;
               ll = -this->ll_old;

               // only iterate over the subset of residues that actually were
               // modified in the preceding move.
               for (int i = start; i < end; i++) {
                    // ridiculous shortcut.
                    ResidueFB *res = &(*this->chain)[i];

                    double tmp = 0;
                    if (settings.ignore_bb) {
                         tmp = this->dbn->calc_ll((*res).residue_type, UNINITIALIZED, UNINITIALIZED, 
                                                  (*res).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    } else {
                         tmp = this->dbn->calc_ll((*res).residue_type, (*res).get_phi(), (*res).get_psi(), 
                                                  (*res).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    }
                    // deduct the old contribution
                    ll -= this->cache[i];
                    // add the new contrib.
                    ll += tmp;
                    // and remember what we did
                    this->roll_back_cache[i] = this->cache[i];
                    this->cache[i] = tmp;

               }

          }
          
          this->ll_move_bias = this->ll_old+ll;
          this->ll_old = -ll;
          return -ll;
     }
     
     //! Get the log bias of the move
     //!
     //! Only used if eliminate-move-bias is set to true. Also, the bias is only
     //! calculated if the move corresponds to the energy.
     //!
     //! This is useful when using backbone dependent sidechain-moves, 
     //! since the energy must be added explicitly to preserve detailed balence.
     //! get_log_bias calculates the bias introduced by the proposed sidechain move
     //! more speed efficiantly, than using the implicit-energy setting.
     //!
     double get_log_bias(MoveInfo *move_info = NULL) {
          MoveInfoSidechainDbn<DBN_TYPE> * move_info_sidechain_dbn = dynamic_cast<MoveInfoSidechainDbn<DBN_TYPE> *>(move_info);
          //Is eliminate-move-bias set to true and does the move correspond to the energy 
          if (settings.eliminate_move_bias && move_info_sidechain_dbn != NULL && this->get_type(this->dbn) == move_info_sidechain_dbn->get_specific_move_type(move_info_sidechain_dbn->dbn)) {
               return -this->ll_move_bias;
          }
          return .0;
     }

     //! reject method
     //!
     //! This method is going to be called every time a move is rejected
     //! by the MCMC method. In this case it is necessary to reset the
     //! the cache properly.
     void reject() {
          for (int i = 0; i < (int) this->roll_back_cache.size(); i++) {
               if (this->roll_back_cache[i] != 0) {
                    this->cache[i] = this->roll_back_cache[i];
                    this->roll_back_cache[i] = 0;
               }
          }
          if (this->ll_cache != 0)
               this->ll_old = this->ll_cache;
          this->ll_cache =0;
     }

     //! accept method
     //!
     //! This method is being called every time a move is accepted by the
     //! MCMC method. Here we need to invalidate the cache.
     void accept() {
          for (int i = 0; i < (int) this->roll_back_cache.size(); i++) {
               this->roll_back_cache[i] = 0;
          }
          this->ll_cache =0;
     }
};

}

#endif
