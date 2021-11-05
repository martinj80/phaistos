// move_sidechain_dbn.h --- Move: Resample sidechain angles using DBN
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


#ifndef MOVE_SIDECHAIN_DBN_H
#define MOVE_SIDECHAIN_DBN_H

#include "protein/chain_fb.h"
#include "moves/move.h"
#include "mocapy.h"
#include "models/basilisk/basilisk_dbn.h"
#include "models/compas/compas_dbn.h"
#include <typeinfo>

namespace phaistos {

namespace definitions {
//! Specific Sidechain Move types
enum SpecificSidechainMoveTypeEnum {BASILISK=11223, COMPAS=2432};
}

template<typename DBN_TYPE>
class MoveInfoSidechainDbn: public MoveInfo {
public:
     DBN_TYPE *dbn;

     static int get_specific_move_type(BasiliskDBN *dbn) {
          return definitions::BASILISK;
     }

     static int get_specific_move_type(CompasDBN *dbn) {
          return definitions::COMPAS;
     }

     MoveInfoSidechainDbn(DBN_TYPE * dbn, definitions::MoveTypeEnum move_type=definitions::NON_LOCAL)
               : MoveInfo(move_type), dbn(dbn) {}
};


//! Sidechain DBN Move
//!
//! This class allows to use our probabilistic models of the
//! amino acid side chain geometry. Using this move, new
//! set of angles will be sampled directly from those models
//! (either for a full atom conformation or one using pseudo
//! side chain beads).
template<typename DBN_TYPE>
class MoveSidechainDBN: public MoveCommon<MoveSidechainDBN<DBN_TYPE>, ChainFB> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveSidechainDBN<DBN_TYPE> , ChainFB> MoveCommon;

     //! Dynamic Bayesian Network object (can be both a Compas or a Basilisk DBN)
     DBN_TYPE *dbn;

     //! Backup of the current angles to allow a rollback
     std::vector<std::vector<double> > chi_backup;

     //! Likelihood of old state - before making move
     double pre_ll;
     
     //! Gaussian random number generator
     boost::variate_generator<RandomNumberEngine&,
                              boost::normal_distribution<> > random_generator_gaussian;

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

          //! Whether to reject moves that produce prolines with broken rings
          bool reject_broken_prolines;

          //! Whether to include hydrogen chi values
          bool sample_hydrogen_chis;

          //! Whether to sample hydrogen chis from a normal distribution (default=uniform)
          bool sample_hydrogen_chis_normal;

          //! The standard deviation of used when sample_hydrogen_chis_normal=true
          double sample_hydrogen_chis_sigma;

          //! Constructor
          Settings(std::string model_filename=set_model_filename_default((DBN_TYPE *)NULL),
                   std::string mocapy_dbn_dir="../data/mocapy-dbns",
                   bool implicit_energy = true, 
                   bool ignore_bb = false, 
                   bool reject_broken_prolines = true, 
                   bool sample_hydrogen_chis = false, 
                   bool sample_hydrogen_chis_normal = true, 
                   double sample_hydrogen_chis_sigma=0.2/180.0*M_PI) :
               Move<ChainFB>::Settings(1,1), // modify only a single sidechain at a time
               model_filename(model_filename),
               mocapy_dbn_dir(mocapy_dbn_dir),
               implicit_energy(implicit_energy), 
               ignore_bb(ignore_bb),
               reject_broken_prolines(reject_broken_prolines), sample_hydrogen_chis(sample_hydrogen_chis), 
               sample_hydrogen_chis_normal(sample_hydrogen_chis_normal), 
               sample_hydrogen_chis_sigma(sample_hydrogen_chis_sigma) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "model-filename:" << settings.model_filename << "\n";
               o << "mocapy-dbn-dir:" << settings.mocapy_dbn_dir << "\n";
               o << "implicit-energy:" << settings.implicit_energy << "\n";
               o << "ignore-bb:" << settings.ignore_bb << "\n";
               o << "reject-broken-prolines:" << settings.reject_broken_prolines << "\n";
               o << "sample-hydrogen-chis:" << settings.sample_hydrogen_chis << "\n";
               o << "sample-hydrogen-chis-normal:" << settings.sample_hydrogen_chis_normal << "\n";
               o << "sample-hydrogen-chis-sigma:" << settings.sample_hydrogen_chis_sigma << "\n";
               o << static_cast<const typename Move<ChainFB>::Settings &> (settings);
               return o;
          }
     } settings;    //!< Local settings object


     //! Constructor.
     //! \param chain Molecule chain
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveSidechainDBN(ChainFB *chain, DBN_TYPE *dbn, const Settings &settings = Settings(), 
               RandomNumberEngine *random_number_engine = &random_global) 
          : MoveCommon(chain, "sc-"+dbn->name, settings, random_number_engine), 
            dbn(dbn),
            pre_ll(-1.0),
            random_generator_gaussian(*random_number_engine,
                                      boost::normal_distribution<>()), 
            settings(settings) {
     }

     //! Copy constructor
     //!
     //! \param other instance of a MoveSidechainDBN to create a copy from
     MoveSidechainDBN(const MoveSidechainDBN<DBN_TYPE> &other) 
          : MoveCommon(other), dbn(other.dbn), pre_ll(other.pre_ll), 
            random_generator_gaussian(*other.random_generator_gaussian,
                                      boost::normal_distribution<>()), settings(other.settings) {
     }

     //! Copy constructor
     //!
     //! Copy constructor to be used in a multithread environment,
     //! allows setting the correct random number engine as well
     //! as choosing the correct instance of the moveSc itself
     //!
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveSidechainDBN(const MoveSidechainDBN<DBN_TYPE> &other, 
               RandomNumberEngine *random_number_engine, 
               int thread_index, 
               ChainFB *chain) 
          : MoveCommon(other, random_number_engine, chain), 
            dbn(&other.dbn->get_copy(thread_index)), 
            pre_ll(other.pre_ll), 
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()), settings(other.settings) {
     }

     //! Apply move
     //!
     //! This method does the actual modification of the conformation.
     //! Depending on the setting, a backbone dependent or independent
     //! sample will be drawn from the according model.
     //!
     //! \param start_index residue index of the first modified residue
     //! \param end_index residue index of the last residue to be modified
     bool apply(int start_index = 0, int end_index = 0) {

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
          this->move_info = (new MoveInfoSidechainDbn<DBN_TYPE>(this->dbn, SIDECHAIN))->add_info(std::make_pair(start_index, end_index),
                                                                                      std::make_pair(start_index, end_index));

          // get all the chi angles we need to get back
          // to the old state afterwards ..
          this->chi_backup = this->chain->get_sidechain_dof_values(this->dbn->get_atom_selection(), 
                                                                   start_index, end_index);

          // if the bias shall be correct, we need
          // to calculate the previous probability
          if (!settings.implicit_energy) {
               this->pre_ll = 0;
               for (int k = start_index; k < end_index; k++) {
                    ResidueFB *res = &(*this->chain)[k];
                    if (settings.ignore_bb) {
                         this->pre_ll += this->dbn->calc_ll((*res).residue_type, 
                                                            UNINITIALIZED, UNINITIALIZED, 
                                                            (*res).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    } else {
                         this->pre_ll += this->dbn->calc_ll((*res).residue_type, 
                                                            (*res).get_phi(), (*res).get_psi(), 
                                                            (*res).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                    }
               }
          }

          std::vector<std::vector<double> > new_angles;
          double random_number=0.0, delta=0.0;
          int k;

          for (int i = start_index; i < end_index; i++) {
               ResidueFB *t = &(*this->chain)[i];
               ResidueEnum res = (*t).residue_type;
               typename ChainFB::Residue *resp = &(*this->chain)[i];
               if (settings.ignore_bb) {
                    new_angles.push_back(this->dbn->get_angle_sequence_for_residue_with_bb(res));
               } else {
                    new_angles.push_back(this->dbn->get_angle_sequence_for_residue_with_bb(res, t->get_phi(), t->get_psi()));
               }
               k=new_angles[i-start_index].size();
               delta=0.0;
               random_number=0.0;
               if (settings.sample_hydrogen_chis) {
                    while(new_angles[i-start_index].size() < resp->chi_atoms.size()) {
                         if(settings.sample_hydrogen_chis_normal) {
                              delta = random_generator_gaussian()*settings.sample_hydrogen_chis_sigma;
                              random_number = fmod(resp->chi_atoms[k].first->get_dihedral() + 3*M_PI + delta, 2*M_PI) - M_PI;
                              k++;
                         } else {
                              random_number=2*M_PI*this->random_generator_uniform_01() - M_PI;
                         }
                         new_angles[i-start_index].push_back(random_number);
                    }
               }

          }

          // put them in the chain
          if (new_angles.size())
               this->chain->set_sidechain_dof_values(new_angles, this->dbn->get_atom_selection(), start_index, end_index);

          // Update positions
          finalize();

          // Check that proline rings were not broken
          if (settings.reject_broken_prolines) {
               this->move_info->success = check_proline_rings();
          } else {
               this->move_info->success = true;
          }

          return this->move_info->success;
     }

     //! Finalize the move: update positions
     //!
     //! After the move, which modified the angles of the sidechain,
     //! it is necessary to update the structure to reflect the modifications
     //! in the conformation as well.
     void finalize() {
          for (int i = this->move_info->modified_angles[0].first; i < this->move_info->modified_angles[0].second; i++) {
               (*this->chain)[i].update_positions_non_backbone();
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

          if (this->chi_backup.size()) {
               this->chain->set_sidechain_dof_values(this->chi_backup, this->dbn->get_atom_selection(),
                                                     this->move_info->modified_angles_start,
                                                     this->move_info->modified_angles_end);

               // Update positions
               finalize();

               // free the backup again
               this->chi_backup.clear();
          }
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
               if (settings.move_length_max == 1) {
                    // if we only have a single rotamer changed
                    // then we can actually ask the dbn directly
                    double bias = this->dbn->calc_ll();
                    return (this->pre_ll - bias);
               } else {
                    // otherwise we need to do it a little more complicated
                    double bias = 0.;

                    std::pair<int,int> tmp = this->move_info->modified_angles.back();

                    for (int k = tmp.first; k < tmp.second; k++) {
                         ResidueFB *res = &(*this->chain)[k];
                         if (settings.ignore_bb) {
                              bias += this->dbn->calc_ll((*res).residue_type, 
                                                         UNINITIALIZED, UNINITIALIZED, 
                                                         (*res).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                         } else {
                              bias += this->dbn->calc_ll((*res).residue_type, 
                                                         (*res).get_phi(), (*res).get_psi(), 
                                                         (*res).get_sidechain_dof_values(this->dbn->get_atom_selection()));
                         }
                    }
                    return (this->pre_ll - bias);
               }
          }
          return 0;
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

          for (int i = this->move_info->modified_angles[0].first; i < this->move_info->modified_angles[0].second; i++) {
               ResidueEnum res = (*this->chain)[i].residue_type;
               if (res == PRO && (*this->chain)[i].has_atom(CD)) {
                    Vector_3D cd = (*this->chain)[i][CD]->position;
                    Vector_3D n = (*this->chain)[i][N]->position;
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
