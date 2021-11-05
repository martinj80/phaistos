// move_pivot_local.h --- Move: Pivot move using Gaussian distributions around current state
// Copyright (C) 2006-2010 Jes Frellsen, Wouter Boomsma
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


#ifndef MOVE_PIVOT_LOCAL_H
#define MOVE_PIVOT_LOCAL_H

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"
#include "move.h"

namespace phaistos {

//! Pivot move using Gaussian distributions around current state
template <typename CHAIN_TYPE, typename PRIOR_TYPE=MovePrior<BondAnglePriorUninformative, 
                                                             DihedralPriorUninformative> >
class MovePivotLocal: public MoveCommon<MovePivotLocal<CHAIN_TYPE,PRIOR_TYPE>,
                                        CHAIN_TYPE> {

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MovePivotLocal<CHAIN_TYPE,PRIOR_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Factory for construction of prior objects
     PRIOR_TYPE prior_factory;

     //! Current prior value
     PRIOR_TYPE *prior;
     
     //! Random number generator - Gaussian
     boost::variate_generator<RandomNumberEngine&, 
                              boost::normal_distribution<> > random_generator_gaussian;

     //@{
     //! Degree-of-freedom iterators
     DofIterator<ChainFB> *begin;
     DofIterator<ChainFB> *backup_begin;
     DofIterator<ChainFB> *end;
     DofIterator<ChainFB> *backup_end;
     //@}

     //@{
     //! Selection of which angles to include in degree-of-freedom iterators
     DofIterator<ChainFB>::AngleSelectionEnum begin_dofs; 
     DofIterator<ChainFB>::AngleSelectionEnum end_dofs; 
     DofIterator<ChainFB>::AngleSelectionEnum backup_begin_dofs; 
     DofIterator<ChainFB>::AngleSelectionEnum backup_end_dofs; 
     definitions::AngleEnum end_dof_type;
     definitions::AtomEnum end_atom;
     //@}

     //@{
     //! Matrices (saved here for get_log_bias)
     Matrix Lt;
     Matrix mu_tilde;
     Matrix W_tilde;
     Matrix dphi;
     //@}

     //! Construct the necessary matrices
     void construct_all(int dof_count, Matrix &Lt, Matrix &W_tilde, Matrix &mu_tilde) {

          // Assumption: a prior that does not induce a bias does not shift the mean
          // If it does not shift the mean, the calculation is simpler
          if (!prior->induces_bias) {

               W_tilde = Matrix(dof_count, dof_count, 0.0);
               mu_tilde = Matrix(dof_count, 1);

               // Construct prior covariance matrix
               prior->construct_multivariate_gaussian(W_tilde, mu_tilde, *begin, *end);
               prior->add_constraint(W_tilde, mu_tilde, *begin, *end, 
                                    settings.std_dev_bond_angle,
                                    settings.std_dev_phi_psi,
                                    settings.std_dev_omega);

          } else  {

               // Construct covariance matrix and mean
               Matrix W(dof_count, dof_count, 0.0);
               Matrix mu_0(dof_count, 1);
          
               // Construct prior covariance matrix
               prior->construct_multivariate_gaussian(W, mu_0, *begin, *end);
               prior->add_constraint(W, mu_0, *begin, *end, 
                                    settings.std_dev_bond_angle,
                                    settings.std_dev_phi_psi,
                                    settings.std_dev_omega);

               W_tilde = W;
               mu_tilde = Matrix(dof_count, 1);

               // Inverse of modified coveriance matrix
               Matrix W_tilde_inverse = W_tilde.inverse();
               
               // Update mu_tilde
               for (int i=0; i<dof_count; i++) {
                    mu_tilde(i,0) = 0.0;
                    for (int j=0; j<dof_count; j++) {
                         for (int k=0; k<dof_count; k++) {
                              mu_tilde(i,0) += W_tilde_inverse(i,j)*W(j,k)*mu_0(k,0);
                         }
                    }
               }
          }

          // Calculate "square root" of matrix W_tilde
          Lt = W_tilde.cholesky('U');
     }

     //! Calculate bias for prerotation move (a->b)
     //! \return log-probability
     double calculate_bias_forward() {
          int size = mu_tilde.n_row;

          double log_det_L = 0.0;
          double sum = 0.0;
          for (int i=0; i<size; i++) {
               log_det_L += std::log(Lt(i,i));
               sum += 0.5*(dphi(i,0) - mu_tilde(i,0)) * W_tilde(i,i) * (dphi(i,0) - mu_tilde(i,0));
               for (int j=i+1; j<size; j++) {
                    sum += (dphi(i,0) - mu_tilde(i,0)) * W_tilde(i,j) * (dphi(j,0) - mu_tilde(j,0));
               }
          }
          // Calculate biasing probability a->b
          return log_det_L - sum;
     }

     //! Calculate bias for reverse move (b->a)
     //! NOTE: Assumes that apply has been called
     //! \return log-probability
     double calculate_bias_reverse() {

          int size = mu_tilde.n_row;

          // Calculate Lt'
          Matrix Lt_prime;
          Matrix W_tilde_prime;
          Matrix mu_tilde_prime;
          construct_all(size, Lt_prime, W_tilde_prime, mu_tilde_prime);
               
          double log_det_L_prime = 0.0;
          double sum = 0.0;
          for (int i=0; i<size; i++) {
               log_det_L_prime += std::log(Lt_prime(i,i));
               sum += 0.5*((-dphi(i,0)) - mu_tilde_prime(i,0)) *
                    W_tilde_prime(i,i) * ((-dphi(i,0)) - mu_tilde_prime(i,0));
               for (int j=i+1; j<size; j++) {
                    sum += ((-dphi(i,0)) - mu_tilde_prime(i,0)) *
                         W_tilde_prime(i,j) * ((-dphi(j,0)) - mu_tilde_prime(j,0));
               }
          }

          // Calculate biasing probability b->a
          return log_det_L_prime - sum;
     }
     


public:

     //! Local Settings class.
     const class Settings: public Move<ChainFB>::Settings, 
                           public PRIOR_TYPE::Settings {
     public:
        
          //! Whether to sample bond_angle degrees of freedom
          bool sample_bond_angle_dofs;

          //! Whether to sample phi,psi degrees of freedom
          bool sample_phi_psi_dofs;

          //! Whether to sample omega degrees of freedom          
          bool sample_omega_dofs;

          //! Standard deviation of bond_angle degrees of freedom (degr.)
          double std_dev_bond_angle;

          //! Standard deviation of phi,psi degrees of freedom (degr.)
          double std_dev_phi_psi;

          //! Standard deviation of omega degrees of freedom (degr.)
          double std_dev_omega;
          
          //! Constructor
          Settings(bool sample_bond_angle_dofs = true,
                   bool sample_phi_psi_dofs = true,
                   bool sample_omega_dofs = false,
                   double std_dev_bond_angle=0.8,
                   double std_dev_phi_psi=4.0,
                   double std_dev_omega=0.8)
               : Move<ChainFB>::Settings(1,1),
                 sample_bond_angle_dofs(sample_bond_angle_dofs),
                 sample_phi_psi_dofs(sample_phi_psi_dofs),
                 sample_omega_dofs(sample_omega_dofs),
                 std_dev_bond_angle(std_dev_bond_angle),
                 std_dev_phi_psi(std_dev_phi_psi),
                 std_dev_omega(std_dev_omega) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "sample-bond-angle_dofs:" << settings.sample_bond_angle_dofs << "\n";
               o << "sample-phi-psi-dofs:" << settings.sample_phi_psi_dofs << "\n";
               o << "sample-omega-dofs:" << settings.sample_omega_dofs << "\n";
               o << "std-dev-bond-angle:" << settings.std_dev_bond_angle << "\n";
               o << "std-dev-phi-psi:" << settings.std_dev_phi_psi << "\n";
               o << "std-dev-omega:" << settings.std_dev_omega << "\n";
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               o << static_cast<const typename PRIOR_TYPE::Settings &>(settings);
               return o;
          }                    
     } settings;    //!< Local settings object
     

     //! Initialization
     void init() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          begin = NULL;     
          backup_begin = NULL;
          end = NULL;
          backup_end = NULL;

          begin_dofs        = DofIterator<ChainFB>::NONE;
          end_dofs          = DofIterator<ChainFB>::NONE;
          backup_begin_dofs = DofIterator<ChainFB>::NONE;
          backup_end_dofs   = DofIterator<ChainFB>::NONE;
          end_atom = CA;
          
          if (settings.sample_bond_angle_dofs) {
               begin_dofs        += DofIterator<ChainFB>::BONDANGLE_DOFS;
               end_dofs          += DofIterator<ChainFB>::BONDANGLE_DOFS;
               backup_begin_dofs += DofIterator<ChainFB>::BONDANGLE_DOFS;
               backup_end_dofs   += DofIterator<ChainFB>::BONDANGLE_DOFS;
               end_dof_type = definitions::ANGLE;
          }
          if (settings.sample_phi_psi_dofs) {
               begin_dofs        += DofIterator<ChainFB>::DIHEDRAL_DOFS;
               end_dofs          += DofIterator<ChainFB>::DIHEDRAL_DOFS;
               backup_begin_dofs += DofIterator<ChainFB>::DIHEDRAL_DOFS;
               backup_end_dofs   += DofIterator<ChainFB>::DIHEDRAL_DOFS;
               end_dof_type = definitions::DIHEDRAL;
          }
          if (settings.sample_omega_dofs) {

               // Test if omega dof is supported by prior
               if (prior_factory.dof_supported(DofIterator<ChainFB>::N_DIHEDRAL) &&
                   prior_factory.dof_supported(DofIterator<ChainFB>::CTERM_N_DIHEDRAL)) {

                    begin_dofs        += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
                    end_dofs          += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
                    backup_begin_dofs += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
                    backup_end_dofs   += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
               }
          }
     }
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MovePivotLocal(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "pivot-local", settings, random_number_engine), 
            prior_factory(settings),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings) {
          
          init();
     }

     //! Constructor - specific for DBN priors
     //! \param chain Molecule chain
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     template <typename DBN_TYPE>
     MovePivotLocal(CHAIN_TYPE *chain, DBN_TYPE *dbn, const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "pivot-local", settings, random_number_engine), 
            prior_factory(dbn, settings),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings) {
          
          init();
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     MovePivotLocal(const MovePivotLocal &other)
          : MoveCommon(other),
            prior_factory(other.settings),
            random_generator_gaussian(*other.random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(other.settings) {
          
          init();
     }
     
     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MovePivotLocal(const MovePivotLocal &other, 
                    RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            prior_factory(other.prior_factory, thread_index),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(other.settings) {
          
          init();
          
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
               const std::pair<int,int> &region = this->random_move_region();
               int length = this->random_move_length(region);
               this->random_move_range(region, length,
                                       -(length-2), this->chain->size()+length-2,
                                       &start_index, &end_index);
          }
          
          // Determine move type
          definitions::TerminalEnum move_region = definitions::INTERNAL;
          if (start_index<0) {
               move_region = definitions::NTERM;
               start_index = 0;
          } else if (end_index > this->chain->size()-1) {
               move_region = definitions::CTERM;
               end_index = this->chain->size();               
          }

	  // Find direction of position update
	  // update_positions takes range of modified angles - end_index
	  // itself was not modified
          int direction;
          if (move_region == definitions::INTERNAL) {
               direction = this->chain->find_shortest_direction(start_index, end_index-1);
          } else {
               direction = ((move_region==definitions::NTERM)?-1:1);
          } 

	  // Set range of modified indices
          std::pair<int,int> modified_positions;
	  if (direction > 0) {
               modified_positions = std::make_pair(start_index, this->chain->size());
	  } else {
               modified_positions = std::make_pair(0, end_index);
	  }

          // Make a copy of the prior - prior initialization
          // might be done in the copy constructor - this ensures
          // that it is executed
          prior = new PRIOR_TYPE(prior_factory, start_index, end_index);


          // std::cout << start_index << " " << end_index << " " << modified_positions << "\n";


          // Set which positions and angles are updated
          // The postrotation determines whether the move is local or not
	  // Remember which angles are modified
          this->move_info = prior->create_move_info(std::make_pair(start_index, end_index), modified_positions, NON_LOCAL);
          //this->move_info = (new MoveInfo(NON_LOCAL))->add_info(std::make_pair(start_index, end_index),
          //                                                      modified_positions);
     
          // Create temporary chain (backup in case of reject)
          this->chain_backup = new CHAIN_TYPE(*this->chain,
                                              this->move_info->modified_positions[0].first,
                                              this->move_info->modified_positions[0].second);

          // Construct Degree-of-freedom iterators
          begin = new DofIterator<ChainFB>((*this->chain)(start_index, N),
                                           definitions::DIHEDRAL,
                                           begin_dofs);
          backup_begin = new DofIterator<ChainFB>((*this->chain_backup)(0, N),
                                                  definitions::DIHEDRAL,
                                                  backup_begin_dofs);
          end = new DofIterator<ChainFB>((*this->chain)(end_index, N),
                                         end_dof_type,
                                         end_dofs);
          backup_end = new DofIterator<ChainFB>((Atom*)NULL,
                                                end_dof_type,
                                                backup_end_dofs);

          // Count degrees of freedom
          int dof_count = 0;
          for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {
               dof_count++;               
          }

          // Special case: only diagonal elements and no induces bias
          if (!prior->off_diagonal_elements && !prior->induces_bias) {

               DiagonalMatrix sigma_inv(dof_count, 0.0);
               Matrix mu_0(dof_count, 1);

               // Construct prior covariance matrix
               prior->construct_multivariate_gaussian(sigma_inv, mu_0, *begin, *end);
               prior->add_constraint(sigma_inv, mu_0, *begin, *end, 
                                    settings.std_dev_bond_angle,
                                    settings.std_dev_phi_psi,
                                    settings.std_dev_omega);

               dphi = Matrix(dof_count, 1);
               for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {
                    dphi(it.offset,0) = (1.0/sqrt(sigma_inv(it.offset,it.offset)))*random_generator_gaussian();
               }

          } else {

               // Construct Lt matrix
               construct_all(dof_count, this->Lt, this->W_tilde, this->mu_tilde);
          
               // Create 0,1 gaussian distributed values
               Matrix dchi(dof_count, 1);
               for (int i=0; i<dof_count; i++) {
                    dchi(i,0) = random_generator_gaussian();
               }

               // Determine changes to angular degrees of freedom
               dphi = solve_ax_b_triangular(Lt, dchi, 'U');

               // Assumption: a prior that does not induce a bias does not shift the mean
               if (prior->induces_bias) {
                    // Shift dphi based on mu_tilde
                    for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {
                         dphi(it.offset,0) += mu_tilde(it.offset,0);               
                    }
               }
          }

          // Construct new conformation
          for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {

               *it += dphi(it.offset,0);
               *it = fmod(*it+3*M_PI,(2*M_PI))-M_PI;

               // Take absolute value if current dof is an angle
               if (it.get_dof_type() == definitions::ANGLE) {
                    *it = std::fabs(*it);
               }

          }

          this->chain->update_positions(start_index,
                                        end_index-1,
                                        direction);               
          
          if (this->move_info->success)
               prior->update_state(this->chain);

          return this->move_info->success;
     }


     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();
          
          prior->accept();

          // Cleanup backup chain
          delete this->chain_backup;
          this->chain_backup = NULL;

          // Delete prior
          delete prior;

          // Cleanup dof-iterators
          delete begin;
          delete end;
          delete backup_begin;
          delete backup_end;
     }


     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();
          
          prior->reject();

          AtomIterator<ChainFB,definitions::ALL> it_current(*this->chain, 
                                                            this->move_info->modified_positions[0].first,
                                                            this->move_info->modified_positions[0].second);
          AtomIterator<ChainFB,definitions::ALL> it_backup(*this->chain_backup);
          for (; !it_current.end() && !it_backup.end(); ++it_current,++it_backup) {
               assert(it_current->atom_type == it_backup->atom_type);
	       it_current->set_angle(it_backup->get_angle());
	       it_current->set_dihedral(it_backup->get_dihedral());
	       it_current->position = it_backup->position;
	  }

          // Cleanup backup chain
          delete this->chain_backup;
          this->chain_backup = NULL;

          // Delete prior
          delete prior;

          // Cleanup dof-iterators
          delete begin;
          delete end;
          delete backup_begin;
          delete backup_end;
     }

     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {

          double bias = 0.0;

          if (!prior->induces_bias) {
               return bias;
          } else {

               double log_prerotation_bias_forward = calculate_bias_forward();
               double log_prerotation_bias_reverse = calculate_bias_reverse();
               bias += log_prerotation_bias_reverse - log_prerotation_bias_forward;
               
               // Bias from prior distribution
               bias += prior->get_log_bias(*begin, *end, *backup_begin, *backup_end);
          }
          return bias;
     }

};

}

#endif
