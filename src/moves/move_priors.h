// move_priors.h --- Priors used by various moves
// Copyright (C) 2006-2010 Wouter Boomsma, Sandro Bottaro, Jesper Ferkinghoff-Borg
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



#ifndef MOVE_PRIORS
#define MOVE_PRIORS

#include "backbone_dbn.h"
#include "protein/protein_ideal_values.h"

namespace phaistos {
     
//! Engh-Huber bond angle move prior
class BondAnglePriorEnghHuber {

     //! Calculate a single entry of the bias (implicit energy) introduced by the prior - internal method (see get_log_bias())
     //! \param atom Atom pointer
     //! \param dof_type Type of degree of freedom
     //! \return log-probability
     inline double calculate_log_bias_entry(Atom *atom,
                                            const definitions::AngleEnum &dof_type,
                                            const double &value) const {

          if (dof_type == definitions::ANGLE) {
               Atom *next_atom = atom->get_neighbour(+1, definitions::BACKBONE);
               Atom *prev_atom = atom->get_neighbour(-1, definitions::BACKBONE);
               std::pair<double,double> param = get_parameters(prev_atom, atom, next_atom);
               double mean = param.first;
               double std_dev = param.second;
               double angle = value;

               // Gaussian density
               return -0.5*(Math<double>::sqr(angle-mean)/Math<double>::sqr(std_dev)) - std::log(std_dev*sqrt(2.0*M_PI));
          }
          return 0;
     }

     //! Calculate the bias (implicit energy) introduced by the prior - internal method (see get_log_bias())
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \return log-probability
     double get_log_bias_internal(const DofIterator<ChainFB> &begin, 
                                  const DofIterator<ChainFB> &end) const {

          double bias = 0.0;

          if (settings.implicit_bond_angle_energy) {
                    
               for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
                    bias += calculate_log_bias_entry(it.get_atom(), it.get_dof_type(), *it);
               }
          }
          return bias;
     }

     //! Calculate the bias (implicit energy) introduced by the prior - internal method (see get_log_bias())
     //! \param dofs Vector of degrees of freedom
     //! \return log-probability
     double get_log_bias_internal(const std::vector<Dof> &dofs) const {

          double bias = 0.0;

          if (settings.implicit_bond_angle_energy) {
                    
               for (unsigned int i=0; i<dofs.size(); ++i) {
                    bias += calculate_log_bias_entry(dofs[i].atom, dofs[i].dof_type, *dofs[i].value);
               }
          }
          return bias;
     }


public:

     //! Whether the prior induces a bias
     const static bool induces_bias = true;

     //! Whether the prior has off-diagonal elements
     const static bool off_diagonal_elements = false;

     //! Local Settings class
     const class Settings {
     public:
          //! Whether to use this prior as an implicit energy. If false
          //! the bias introduced by the prior will be "divided out" in
          //! the acceptance ratio
          bool implicit_bond_angle_energy;

          //! Constructor
          Settings(bool implicit_bond_angle_energy=true):
               implicit_bond_angle_energy(implicit_bond_angle_energy) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "implicit-bond-angle-energy:" << settings.implicit_bond_angle_energy << "\n";
               return o;
          }                    
     } settings;    //!< Local settings object

          
     //! String identifier
     std::string id;

          
     //! Constructor
     //! \param settings Local settings object
     BondAnglePriorEnghHuber(const Settings &settings):
          settings(settings), id("EnghHuber") {}

     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     BondAnglePriorEnghHuber(const BondAnglePriorEnghHuber &other, 
                             int start_index, int end_index):
          settings(other.settings), id(other.id){}


     //! Construct a multivariate gaussian prior (covariance matrix and mean)
     //! \param sigma Target covariance matrix
     //! \param mu Target mean vector (1 column matrix)
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     template <typename MATRIX_TYPE>
     void construct_multivariate_gaussian(MATRIX_TYPE &sigma, Matrix &mu,
                                          const DofIterator<ChainFB> &start_dof, 
                                          const DofIterator<ChainFB> &end_dof) const {
                                                  
          // Iterate over all bond-angle degrees-of-freedom
          for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
               
               if (it.get_dof_type() == definitions::ANGLE) {

                    // Get mean and stddev
                    Atom *atom = it.get_atom();
                    Atom *next_atom = atom->get_neighbour(+1, definitions::BACKBONE);
                    Atom *prev_atom = atom->get_neighbour(-1, definitions::BACKBONE);
                    std::pair<double,double> param = get_parameters(prev_atom, atom, next_atom);
                    double mean = param.first;
                    double std_dev = param.second;

                    sigma(it.offset, it.offset) = 1.0/(Math<double>::sqr(std_dev));
                    mu(it.offset,0) = mean - *it;

               }
          }
     }

     //! Constrain the covariance matrix
     //! NOTE: this method should be called on an unmodified covariance matrix
     //!       produced by construct_multivariate_gaussian()
     //! \param sigma covariance matrix
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     //! \param std_dev_bond_angle Defines size of contraint. If infinite, this function is a no-op.
     //! \return Boolean specifying whether constraint was added
     template <typename MATRIX_TYPE>
     bool add_constraint(MATRIX_TYPE &sigma, 
                         const DofIterator<ChainFB> &start_dof, 
                         const DofIterator<ChainFB> &end_dof, 
                         double std_dev_bond_angle) const {

          bool constraint_added = false;

          // Only constraint if std_dev_bond_angle is finite
          if (std::isfinite(std_dev_bond_angle)) {
                    
               double total_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(std_dev_bond_angle));

               constraint_added = true;

               for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
                    if (it.get_dof_type() == definitions::ANGLE) {

                         // Back-calculate std_dev from covariance matrix entry
                         // (NOTE: this assumes that add constraint is called on
                         // an unmodified covariance matrix)                              
                         double std_dev = sqrt(1.0/sigma(it.offset, it.offset));
                              
                         if (std_dev_bond_angle < std_dev*(180/M_PI)) {
                              sigma(it.offset, it.offset) = total_constraint;
                         } 
                    }
               }
          }

          return constraint_added;
     }


     //! Update internal state based on newly modified chain state:
     //! Insert new angles into dbn
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void update_state(CHAIN_TYPE *chain) {}


          
     //! Return angle parameters for atom triplet
     //! \return (mean, stddev) pair
     std::pair<double,double> get_parameters(Atom *atom1, Atom *atom2, Atom *atom3) const {
          double std_dev;
          double mean = bond_angle_constants(atom1->atom_type, atom2->atom_type,
                                           atom3->atom_type, atom2->residue->residue_type,
                                           &std_dev);
          return std::make_pair(mean, std_dev);
     }

          
     //! Calculate the bias (implicit energy) introduced by the prior.
     //! If this prior is to be used as an implicit energy
     //! we include the likelihood of the angles as a bias term
     //! note: the we are not dividing out a proposal term here -
     //! this is already done by the prerotation itself. We are simply
     //! putting in the energy corresponding to the prior in the proposal, 
     //! to ensure that even when simulating without energies, the 
     //! bond_angle is kept within reasonable bounds.
     //!
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \param backup_begin Degree-of-freedom iterator start - backup chain
     //! \param backup_end Degree-of-freedom iterator end - backup chain
     //! \return log-probability
     double get_log_bias(const DofIterator<ChainFB> &begin, 
                         const DofIterator<ChainFB> &end,
                         const DofIterator<ChainFB> &backup_begin, 
                         const DofIterator<ChainFB> &backup_end) const {

          double bond_angle_LL_after = get_log_bias_internal(begin, end);
          double bond_angle_LL_before = get_log_bias_internal(backup_begin, backup_end);
               
          return bond_angle_LL_after - bond_angle_LL_before;
     }

     //! Calculate the bias (implicit energy) introduced by the prior.
     //! If this prior is to be used as an implicit energy
     //! we include the likelihood of the angles as a bias term
     //! note: the we are not dividing out a proposal term here -
     //! this is already done by the prerotation itself. We are simply
     //! putting in the energy corresponding to the prior in the proposal, 
     //! to ensure that even when simulating without energies, the 
     //! bond_angle is kept within reasonable bounds.
     //!
     //! \param dofs Vector of degrees of freedom
     //! \param dofs_backup Vector of degrees of freedom - backup chain
     //! \return log-probability
     double get_log_bias(const std::vector<Dof> &dofs, 
                         const std::vector<Dof> &dofs_backup) const {

          double bond_angle_LL_after = get_log_bias_internal(dofs);
          double bond_angle_LL_before = get_log_bias_internal(dofs_backup);
               
          return bond_angle_LL_after - bond_angle_LL_before;
     }

     //! Accept last move
     void accept() {}

     //! Reject last move
     void reject() {}          

     // Whether a particular dof is supported by this prior
     bool dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
          return true;
     }
};


//! Bond angle move prior - uninformative
class BondAnglePriorUninformative {
public:

     //! Whether this prior induces a bias
     const static bool induces_bias = false;

     //! Whether this prior has off-diagonal elements
     const static bool off_diagonal_elements = false;

     //! Local Settings class.
     const class Settings {
     public:

          // Constructor
          Settings() {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               return o;
          }                    

     } settings;    //!< Local settings object


     //! String identifier
     std::string id;
	  
     //! Constructor
     //! \param settings Local settings object
     BondAnglePriorUninformative(const Settings &settings):
          settings(settings), id("bond-angle-uninformative") {}


     //! Copy Constructor - thread specified (multi-threading)
     //! \param other Source object from which copy is made
     //! \param thread_index Which thread the copy will run in
     BondAnglePriorUninformative(const BondAnglePriorUninformative &other, int thread_index):
          settings(other.settings), id(other.id) {}
	  

     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     BondAnglePriorUninformative(const BondAnglePriorUninformative &other,
                                  int start_index, int end_index):
          settings(other.settings), id(other.id) {}


     //! Construct a multivariate gaussian prior (covariance matrix and mean)
     //! \param sigma Target covariance matrix
     //! \param mu Target mean vector (1 column matrix)
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     template <typename MATRIX_TYPE>
     void construct_multivariate_gaussian(MATRIX_TYPE &sigma, Matrix &mu,
                                          const DofIterator<ChainFB> &start_dof,
                                          const DofIterator<ChainFB> &end_dof) const {
	    
          for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {

               if (it.get_dof_type() == definitions::ANGLE) {

                    assert(it.offset < sigma.n_row);
                    sigma(it.offset, it.offset) = 1.0;
                    mu(it.offset,0) = 0.0;
		   
               }
          }
	    
     }

          
     //! Update internal state based on newly modified chain state: do nothing
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void update_state(CHAIN_TYPE *chain) {}


     //! Constrain the covariance matrix
     //! NOTE: this method should be called on an unmodified covariance matrix
     //!       produced by construct_multivariate_gaussian()
     //! \param sigma covariance matrix
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     //! \param std_dev_bond_angle Defines size of contraint. If infinite, this function is a no-op.
     //! \return Boolean specifying whether constraint was added
     template <typename MATRIX_TYPE>
     bool add_constraint(MATRIX_TYPE &sigma, 
                         const DofIterator<ChainFB> &start_dof, const DofIterator<ChainFB> &end_dof,
                         double std_dev_bond_angle) const {
          bool constraint_added = false;
  
          // Only constraint if std_dev_bond_angle is finite
          if (std::isfinite(std_dev_bond_angle)) {
               constraint_added = true;
               double total_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(std_dev_bond_angle));
                    
               for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {

                    if (it.get_dof_type() == definitions::ANGLE) {
                         sigma(it.offset, it.offset) += total_constraint;
                    }
               }
          } 
          return constraint_added;

     }

          
     //! Calculate the bias (implicit energy) introduced by the prior.
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \param backup_begin Degree-of-freedom iterator start - backup chain
     //! \param backup_end Degree-of-freedom iterator end - backup chain
     //! \return log-probability
     double get_log_bias(const DofIterator<ChainFB> &begin, const DofIterator<ChainFB> &end,
                         const DofIterator<ChainFB> &backup_begin, const DofIterator<ChainFB> &backup_end) const {
          return 0.0;
     }

     //! Calculate the bias (implicit energy) introduced by the prior.
     //! \param dofs Vector of degrees of freedom
     //! \param dofs_backup Vector of degrees of freedom - backup chain
     //! \return log-probability
     double get_log_bias(const std::vector<Dof> &dofs,
                         const std::vector<Dof> &dofs_backup) const {
          return 0.0;
     }


     //! Accept last move
     void accept() {}

     //! Reject last move
     void reject() {}

     // Whether a particular dof is supported by this prior
     bool dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
          return true;
     }
};


               
//! Dihedral angle move prior - uninformative
class DihedralPriorUninformative {
public:

     //! Whether this prior induces a bias
     const static bool induces_bias = false;

     //! Whether this prior has off-diagonal elements
     const static bool off_diagonal_elements = false;

     //! Local Settings class.
     const class Settings {
     public:

          //! Constructor
          Settings() {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               return o;
          }                    
     } settings;    //!< Local settings object

     //! String identifier
     std::string id;
          
     //! Constructor
     //! \param settings Local settings object
     DihedralPriorUninformative(const Settings &settings):
          settings(settings), id("dihedral_uninformative") {}


     //! Copy Constructor - thread specified (multi-threading)
     //! \param other Source object from which copy is made
     //! \param thread_index Which thread the copy will run in
     DihedralPriorUninformative(const DihedralPriorUninformative &other, int thread_index):
          settings(other.settings), id(other.id) {}


     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     DihedralPriorUninformative(const DihedralPriorUninformative &other,
                                 int start_index, int end_index):
          settings(other.settings), id(other.id) {}


     //! Construct a multivariate gaussian prior (covariance matrix and mean)
     //! \param sigma Target covariance matrix
     //! \param mu Target mean vector (1 column matrix)
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     template <typename MATRIX_TYPE>
     void construct_multivariate_gaussian(MATRIX_TYPE &sigma, Matrix &mu,
                                          const DofIterator<ChainFB> &start_dof,
                                          const DofIterator<ChainFB> &end_dof) const {

          for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
               
               if (it.get_dof_type() == definitions::DIHEDRAL) {
               
                    sigma(it.offset, it.offset) = 1.0;
                    mu(it.offset,0) = 0.0;

               }
          }
     }

          
     //! Update internal state based on newly modified chain state: do nothing
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void update_state(CHAIN_TYPE *chain) {}


     //! Constrain the covariance matrix
     //! NOTE: this method should be called on an unmodified covariance matrix
     //!       produced by construct_multivariate_gaussian()
     //! \param sigma covariance matrix
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     //! \param std_dev_phi_psi Defines size of phi,psi dihedral contraint. If infinite, no constraint is added.
     //! \param std_dev_omega Defines size of omega dihedral contraint. If infinite, no constraint is added
     //! \return Boolean specifying whether constraint was added
     template <typename MATRIX_TYPE>
     bool add_constraint(MATRIX_TYPE &sigma, 
                         const DofIterator<ChainFB> &start_dof, const DofIterator<ChainFB> &end_dof,
                         double std_dev_phi_psi, double std_dev_omega) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          bool constraint_added = false;

          // Only constraint if std_dev_phi_psi is finite
          if (std::isfinite(std_dev_phi_psi)) {

               constraint_added = true;
               double phi_psi_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(std_dev_phi_psi));
               double omega_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(std_dev_omega));
                    
               for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
                    if (it.get_dof_type() == definitions::DIHEDRAL) {
                         if(it.get_atom()->atom_type == N) {
                              sigma(it.offset, it.offset) += omega_constraint;
                         } else { 
                              sigma(it.offset, it.offset) += phi_psi_constraint;
                         }
                    }
               }
          }
          return constraint_added;
     }

          
     //! Calculate the bias (implicit energy) introduced by the prior.
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \param backup_begin Degree-of-freedom iterator start - backup chain
     //! \param backup_end Degree-of-freedom iterator end - backup chain
     //! \return log-probability
     double get_log_bias(const DofIterator<ChainFB> &begin, const DofIterator<ChainFB> &end,
                         const DofIterator<ChainFB> &backup_begin, const DofIterator<ChainFB> &backup_end) const {
          return 0.0;
     }

     //! Calculate the bias (implicit energy) introduced by the prior.
     //! \param dofs Vector of degrees of freedom
     //! \param dofs_backup Vector of degrees of freedom - backup chain
     //! \return log-probability
     double get_log_bias(const std::vector<Dof> &dofs,
                         const std::vector<Dof> &dofs_backup) const {
          return 0.0;
     }

     //! Accept last move
     void accept() {}

     //! Reject last move
     void reject() {}

     // Whether a particular dof is supported by this prior
     bool dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
          return true;
     }
};

//! Dihedral angle prior - Using Dynamic Bayesian Network
template <typename DBN_TYPE=TorusDBN>
class DihedralPriorDbn {
public:

     //! Whether this prior induces a bias
     const static bool induces_bias = true;

     //! Whether this prior has off-diagonal elements
     const static bool off_diagonal_elements = true;

     //! Local Settings class.
     const class Settings {
     public:

          //! Whether to use this prior as an implicit energy. If false
          //! the bias introduced by the prior will be "divided out" in
          //! the acceptance ratio. If true, the bias will also be divided out,
          //! and replaced by the correct DBN contribution (since the proposal
          //! is only an approximation to the real distribution)
          bool implicit_dihedral_energy;

          //! Whether to resample the hidden node sequence of the DBN
          //! whenever a prior of this type is constructed
          bool resample_hidden_nodes;

          //! Size of window used (to each side) when bringing the dbn back to consistency 
          //! (negative number => resample entire sequence)
          int dbn_consistency_window_size;

          //! Size of window used when calculating bias. Approximates 
          //! the move bias as P(X[i-w,j+w])/P(X'[i-w,j+w]), where w is 
          //! the window size and [i,j] is the interval where angles have been
          //! changed. A good value for the window size is >7, and
          //! a negative window size means that the full bias is used.
          int dbn_bias_window_size;

          //! Constructor
          Settings(bool implicit_dihedral_energy=true,
                   bool resample_hidden_nodes=true,
                   int dbn_consistency_window_size=10,
                   int dbn_bias_window_size=10):
               implicit_dihedral_energy(implicit_dihedral_energy),
               resample_hidden_nodes(resample_hidden_nodes),
               dbn_consistency_window_size(dbn_consistency_window_size),
               dbn_bias_window_size(dbn_bias_window_size) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "implicit-dihedral-energy:" << settings.implicit_dihedral_energy << "\n";
               o << "resample-hidden-nodes:" << settings.resample_hidden_nodes << "\n";
               o << "dbn-consistency-window-size:" << settings.dbn_consistency_window_size << "\n";
               o << "dbn-bias-window-size:" << settings.dbn_bias_window_size << "\n";
               return o;
          }                    

     } settings;    //!< Local settings object

     //! String identifier
     std::string id;
          
     //! Pointer to the Dynamic Bayesian Network model used
     DBN_TYPE *dbn;

     //! Range of prior - start index
     int start_index;

     //! Range of prior - end index
     int end_index;
          
     //! P(h|x) probability (saved as attribute because it is used for bias calculation)
     double prob_h_x;

     //! DBN parameters - phi,psi kappa values;
     std::vector<std::vector<double> > kappas;

     //! DBN parameters - phi,psi mean values;
     std::vector<std::vector<double> > means;

     //! DBN parameters - omega kappa values;
     std::vector<double> omega_kappas;

     //! DBN parameters - omega mean values;
     std::vector<double> omega_means;


     //! Set omega parameters - only relevant for certain models
     //! \param dbn Pointer to dynamic Bayesian network model.
     void initialize_omega_shape_parameters(DBN_TYPE *dbn);
          
     //! Constructor
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local settings object
     DihedralPriorDbn(DBN_TYPE *dbn, const Settings &settings):
          settings(settings), id("DBN"), dbn(dbn) {

          kappas = dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->get_shape_parameters();
          means = dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->get_means();

          initialize_omega_shape_parameters(dbn);
     }


     //! Copy Constructor - thread specified (multi-threading)
     //! \param other Source object from which copy is made
     //! \param thread_index Which thread the copy will run in
     DihedralPriorDbn(const DihedralPriorDbn &other, int thread_index):
          settings(other.settings), id(other.id), dbn(&other.dbn->get_copy(thread_index)),
          start_index(start_index), end_index(end_index),
          kappas(other.kappas), means(other.means),
          omega_kappas(other.omega_kappas), omega_means(other.omega_means) {
     }


     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     DihedralPriorDbn(const DihedralPriorDbn &other, int start_index, int end_index)
          : settings(other.settings), id(other.id), dbn(other.dbn),
            start_index(start_index), end_index(end_index),
            kappas(other.kappas), means(other.means),
            omega_kappas(other.omega_kappas), omega_means(other.omega_means) {

          // Make sure that DBN is consistent
          dbn->enforce_consistency(settings.dbn_consistency_window_size);

          // Resample hidden nodes based on current angles
          // prob_h_x = dbn->hiddenNode->resample(start_index, end_index);
          dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->sample(start_index, end_index);

          prob_h_x = 0.0;
          if (!settings.implicit_dihedral_energy) {
               bool sampling = false;
               if (settings.dbn_bias_window_size < 0) {
                    prob_h_x = dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()
                         ->sample(-1, -1, sampling);
               } else {
                    prob_h_x = dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()
                         ->sample(start_index-settings.dbn_bias_window_size, 
                                  end_index  +settings.dbn_bias_window_size, sampling);
               }
          }

          if (settings.resample_hidden_nodes) {
               // Always accept
               dbn->accept();
          } else {
               // Reject the new hidden node sequence (we just need the hDistribution_forward)
               dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->reject();
          }
     }

          
     //! Update internal state based on newly modified chain state: 
     //! Insert new angles into dbn
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void update_state(CHAIN_TYPE *chain) {

          std::vector<std::vector<double> > angles;
          bool include_omega = true;
          for (int i=start_index; i<end_index; i++) {
               angles.push_back((*chain)[i].get_angles(include_omega));
          }
          dbn->template set_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES>(angles, start_index);
     }

          
     //! Construct a multivariate gaussian prior (covariance matrix and mean)
     //! \param sigma Target covariance matrix
     //! \param mu Target mean vector (1 column matrix)
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     template <typename MATRIX_TYPE>
     void construct_multivariate_gaussian(MATRIX_TYPE &sigma, Matrix &mu,
                                          const DofIterator<ChainFB> &start_dof, const DofIterator<ChainFB> &end_dof) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          std::vector<double> means = std::vector<double>(sigma.n_col, 0.0);

          // Keep track of last phi, psi index to match up phi, psi pairs
          int last_phi_dof_index = -1;
          int last_phi_res_index = -1;
          for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
               
               if (it.get_dof_type() == definitions::DIHEDRAL) {
               
                    AtomEnum atom_type = it.get_atom()->atom_type;
                    int residue_index = it.get_residue()->index;

                    int h = dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->sequence[residue_index];

                    double k1 = kappas[h][0];
                    double k2 = kappas[h][1];
                    double k3 = kappas[h][2];

                    if (atom_type == CA) { // phi
                         last_phi_res_index = residue_index;
                         last_phi_dof_index = it.offset;
                         sigma(it.offset, it.offset) = k1-k3;
                         means[it.offset] = this->means[h][0];

                    } else if (atom_type==C) {	           // psi
                         sigma(it.offset, it.offset) = k2-k3;
                         if (last_phi_res_index == residue_index) {
                              sigma(last_phi_dof_index, it.offset) = k3;
                              sigma(it.offset, last_phi_dof_index) = k3;
                         }
                         last_phi_res_index = -1;
                         last_phi_dof_index = -1;
                         means[it.offset] = this->means[h][1];

                    } else if (atom_type==N) { // omega

                         int cis = dbn->template get_node<typename DBN_TYPE::CIS_NODE>()->sequence[residue_index];
                         double kappa = omega_kappas[cis];
                         sigma(it.offset, it.offset) = kappa;
                         means[it.offset] = omega_means[cis];
                    }
                    mu(it.offset,0) = (means[it.offset] - *it);
               }
          }
     }

          
     //! Constrain the covariance matrix
     //! NOTE: this method should be called on an unmodified covariance matrix
     //!       produced by construct_multivariate_gaussian()
     //! \param sigma covariance matrix
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     //! \param std_dev_phi_psi Defines size of phi,psi dihedral contraint. If infinite, no constraint is added.
     //! \param std_dev_omega Defines size of omega dihedral contraint. If infinite, no constraint is added
     //! \return Boolean specifying whether constraint was added
     template <typename MATRIX_TYPE>
     bool add_constraint(MATRIX_TYPE &sigma,
                         const DofIterator<ChainFB> &start_dof, const DofIterator<ChainFB> &end_dof, 
                         double std_dev_phi_psi, double std_dev_omega) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          bool constraint_added = false;

          // Only constrain if std_dev_phi_psi is finite
          if (std::isfinite(std_dev_phi_psi)) {

               constraint_added = true;

               double phi_psi_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(std_dev_phi_psi));
               double omega_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(std_dev_omega));

               for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
                    if (it.get_dof_type() == definitions::DIHEDRAL) {
                         if(it.get_atom()->atom_type == N) {
                              sigma(it.offset, it.offset) += omega_constraint;
                         } else { 
                              sigma(it.offset, it.offset) += phi_psi_constraint;
                         }
                    }
               }
          }
               
          return constraint_added;
     }
          

     //! Calculate the bias (implicit energy) introduced by the prior.
     //! If this prior is to be used as an implicit energy
     //! we include the likelihood of the angles as a bias term
     //! note: the we are not dividing out a proposal term here -
     //! this is already done by the prerotation itself. We are simply
     //! putting in the energy corresponding to the prior in the proposal, 
     //! to ensure that even when simulating without energies, the 
     //! bond_angle is kept within reasonable bounds.
     //! \return log-probability
     double get_log_bias() const {
          double bias = 0.0;

          // If the DBN is used as an implicit energy
          // we include the usual LL as a bias term
          if (settings.implicit_dihedral_energy) {

               // Calculate likelihood before move
               bool use_backup = true;
               double old_LL = dbn->get_log_likelihood_conditional(start_index, end_index,
                                                                   use_backup);
               // Calculate likelihood after move
               double new_LL = dbn->get_log_likelihood_conditional(start_index, end_index);

               bias += (new_LL - old_LL);


               // If the DBN is not used an an implicit energy
               // the P(h|x) bias introduced by the prior should be divided out
          } else {
                    
               // Get P(h|x) using a forward-backtrack pass without sampling
               double prob_h_xprime = 0.0;
               bool sampling = false;
               if (settings.dbn_bias_window_size < 0) {
                    prob_h_xprime = dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->
                         sample(-1, -1, sampling);
               } else {
                    prob_h_xprime = dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->
                         sample(start_index-settings.dbn_bias_window_size, 
                                end_index  +settings.dbn_bias_window_size, sampling);
               }
                    
               // The hidden node state produced here is not used
               dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->reject();

               bias += prob_h_xprime - prob_h_x;

          }
          return bias;
     }
          
     //! Calculate the bias (implicit energy) introduced by the prior.
     //! If this prior is to be used as an implicit energy
     //! we include the likelihood of the angles as a bias term
     //! note: the we are not dividing out a proposal term here -
     //! this is already done by the prerotation itself. We are simply
     //! putting in the energy corresponding to the prior in the proposal, 
     //! to ensure that even when simulating without energies, the 
     //! bond_angle is kept within reasonable bounds.
     //!
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \param backup_begin Degree-of-freedom iterator start - backup chain
     //! \param backup_end Degree-of-freedom iterator end - backup chain
     //! \return log-probability
     double get_log_bias(const DofIterator<ChainFB> &begin, const DofIterator<ChainFB> &end,
                         const DofIterator<ChainFB> &backup_begin, const DofIterator<ChainFB> &backup_end) const {
          return get_log_bias();
     }

     //! Calculate the bias (implicit energy) introduced by the prior.
     //! If this prior is to be used as an implicit energy
     //! we include the likelihood of the angles as a bias term
     //! note: the we are not dividing out a proposal term here -
     //! this is already done by the prerotation itself. We are simply
     //! putting in the energy corresponding to the prior in the proposal, 
     //! to ensure that even when simulating without energies, the 
     //! bond_angle is kept within reasonable bounds.
     //!
     //! \param dofs Vector of degrees of freedom
     //! \param dofs_backup Vector of degrees of freedom - backup chain
     //! \return log-probability
     double get_log_bias(const std::vector<Dof> &dofs,
                         const std::vector<Dof> &dofs_backup) const {
          return get_log_bias();
     }

     //! Accept last move
     void accept() {
          // Call accept in DBN
          dbn->accept();
     }

     //! Reject last move
     void reject() {
          // Call reject in DBN
          dbn->reject();
     }          

     // Whether a particular dof is supported by this prior
     bool dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const;
};

//! Set omega parameters - only relevant for certain models
template <typename DBN_TYPE>
void DihedralPriorDbn<DBN_TYPE>::initialize_omega_shape_parameters(DBN_TYPE *dbn) {
}

//! Set omega parameters - only relevant for certain models
template <>
void DihedralPriorDbn<TorusOmegaDBN>::initialize_omega_shape_parameters(TorusOmegaDBN *dbn) {
     omega_kappas = dbn->get_node<TorusOmegaDBN::OMEGA_NODE>()->get_shape_parameters();
     omega_means = dbn->get_node<TorusOmegaDBN::OMEGA_NODE>()->get_means();          
}

template <>
void DihedralPriorDbn<TorusCsOmegaDBN>::initialize_omega_shape_parameters(TorusCsOmegaDBN *dbn) {
     omega_kappas = dbn->get_node<TorusOmegaDBN::OMEGA_NODE>()->get_shape_parameters();
     omega_means = dbn->get_node<TorusOmegaDBN::OMEGA_NODE>()->get_means();          
}

//! Whether a particular dof is supported by this prior - general version
template <typename DBN_TYPE>
bool DihedralPriorDbn<DBN_TYPE>::dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
     return true;
}

//! Whether a particular dof is supported by this prior - general version
template <>
bool DihedralPriorDbn<TorusDBN>::dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
     if (dof_type & (DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL))
          return false;
     else
          return true;
}

template <>
bool DihedralPriorDbn<TorusCsDBN>::dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
     if (dof_type & (DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL))
          return false;
     else
          return true;
}

//! Move prior. Consists of a bond angle and a dihedral prior
template <typename BONDANGLE_PRIOR_TYPE=BondAnglePriorEnghHuber,
          typename DIHEDRAL_PRIOR_TYPE=DihedralPriorUninformative>
class MovePrior {
public:

     //! Whether this prior induces a bias
     const static bool induces_bias = BONDANGLE_PRIOR_TYPE::induces_bias || DIHEDRAL_PRIOR_TYPE::induces_bias;

     //! Whether this prior has off-diagonal elements
     const static bool off_diagonal_elements = BONDANGLE_PRIOR_TYPE::off_diagonal_elements || DIHEDRAL_PRIOR_TYPE::off_diagonal_elements;

     //! Local Settings class - inherits settings from dihedralprior and bond-angle prior
     const class Settings: public DIHEDRAL_PRIOR_TYPE::Settings,
                           public BONDANGLE_PRIOR_TYPE::Settings {
     public:
          //! Constructor
          Settings()
               : DIHEDRAL_PRIOR_TYPE::Settings(),
                 BONDANGLE_PRIOR_TYPE::Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename DIHEDRAL_PRIOR_TYPE::Settings &>(settings);
               o << static_cast<const typename BONDANGLE_PRIOR_TYPE::Settings &>(settings);
               return o;
          }                    
     } settings;    //!< Local settings object
     
     //! Dihedral prior
     DIHEDRAL_PRIOR_TYPE dihedral_prior;

     //! Bond-angle prior
     BONDANGLE_PRIOR_TYPE bond_angle_prior;

     //! String identifier
     std::string id;

     //! Constructor
     //! \param settings Local settings object
     MovePrior(const Settings &settings=Settings()):
          settings(settings),
          dihedral_prior(settings), bond_angle_prior(settings),
          id(dihedral_prior.id+","+bond_angle_prior.id) {}

          
     //! Constructor (used for construction with DBN - see MovePriorDbn class)
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local settings object
     template <typename DBN_TYPE>          
     MovePrior(DBN_TYPE *dbn, const Settings &settings):
          settings(settings),
          dihedral_prior(dbn,settings), bond_angle_prior(settings),
          id(dihedral_prior.id+","+bond_angle_prior.id) {}

          
     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     MovePrior(const MovePrior &other, int start_index, int end_index):
          settings(other.settings),
          dihedral_prior(other.dihedral_prior, start_index, end_index),
          bond_angle_prior(other.bond_angle_prior, start_index, end_index),
          id(other.id) {}

     //! Copy Constructor - thread specified (multi-threading)
     //! \param other Source object from which copy is made
     //! \param thread_index Which thread the copy will run in
     MovePrior(const MovePrior &other, int thread_index):
          settings(other.settings),
          dihedral_prior(other.dihedral_prior, thread_index),
          bond_angle_prior(other.bond_angle_prior),
          id(other.id) {}

     //! Construct a multivariate gaussian prior (covariance matrix and mean)
     //! \param sigma Target covariance matrix
     //! \param mu Target mean vector (1 column matrix)
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     template <typename MATRIX_TYPE>
     void construct_multivariate_gaussian(MATRIX_TYPE &sigma, Matrix &mu,
                                          const DofIterator<ChainFB> &start_dof, const DofIterator<ChainFB> &end_dof) const {
          bond_angle_prior.construct_multivariate_gaussian(sigma, mu, start_dof, end_dof);
          dihedral_prior.construct_multivariate_gaussian(sigma, mu, start_dof, end_dof);
     }


     //! Constrain the covariance matrix          
     //! Calls the add_constraint of the bond_angle and dihedral priors, and
     //! subsequently recalculates the means (mu).
     //! \param sigma covariance matrix
     //! \param mu mean vector (1 column matrix)
     //! \param start_dof Degree-of-freedom iterator start
     //! \param end_dof Degree-of-freedom iterator end
     //! \param std_dev_bond_angle Defines size of contraint. If infinite, this function is a no-op.
     //! \param std_dev_phi_psi Defines size of phi,psi dihedral contraint. If infinite, no constraint is added.
     //! \param std_dev_omega Defines size of omega dihedral contraint. If infinite, no constraint is added
     template <typename MATRIX_TYPE>
     void add_constraint(MATRIX_TYPE &sigma, Matrix &mu,
                         const DofIterator<ChainFB> &start_dof, 
                         const DofIterator<ChainFB> &end_dof, 
                         double std_dev_bond_angle, 
                         double std_dev_phi_psi,
                         double std_dev_omega) const {

          // Assumption: a prior that does not induce a bias does not shift the mean
          if (!this->induces_bias) {

               bond_angle_prior.add_constraint(sigma, 
                                              start_dof, end_dof, 
                                              std_dev_bond_angle);
               dihedral_prior.add_constraint(sigma,
                                             start_dof, end_dof, 
                                             std_dev_phi_psi,
                                             std_dev_omega);

          } else {
               
               MATRIX_TYPE sigma_constrained(sigma);
               bool constraint_added_bond_angle = 
                    bond_angle_prior.add_constraint(sigma_constrained, 
                                                   start_dof, end_dof, 
                                                   std_dev_bond_angle);
               bool constraint_added_dihedral = 
                    dihedral_prior.add_constraint(sigma_constrained,
                                                  start_dof, end_dof, 
                                                  std_dev_phi_psi,
                                                  std_dev_omega);

               // If constraints were added, mu must be updated
               if (constraint_added_bond_angle || constraint_added_dihedral) {

                    int size = sigma.n_row;
                    Matrix mu_constrained = Matrix(size, 1);
               
                    // Inverse of modified coveriance matrix
                    MATRIX_TYPE sigma_constrained_inverse = sigma_constrained.inverse();
               
                    // Update mu_tilde
                    for (int i=0; i<size; i++) {
                         mu_constrained(i,0) = 0.0;
                         for (int j=0; j<size; j++) {
                              for (int k=0; k<size; k++) {
                                   mu_constrained(i,0) += (sigma_constrained_inverse(i,j)*
                                                           sigma(j,k)*mu(k,0));
                              }
                         }
                    }
                    mu = mu_constrained;
                    sigma = sigma_constrained;
               }
          }
     }


     //! Update of internal state based on newly modified chain state
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void update_state(CHAIN_TYPE *chain) {
          dihedral_prior.update_state(chain);
          bond_angle_prior.update_state(chain);
     }


     //! Calculate the bias (implicit energy) introduced by the prior.
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \param backup_begin Degree-of-freedom iterator start - backup chain
     //! \param backup_end Degree-of-freedom iterator end - backup chain
     //! \return log-probability
     double get_log_bias(const DofIterator<ChainFB> &begin, const DofIterator<ChainFB> &end,
                         const DofIterator<ChainFB> &backup_begin, const DofIterator<ChainFB> &backup_end) const {

          double bias = 0.0;

          // Bias due to dihedral prior
          bias += dihedral_prior.get_log_bias(begin, end, backup_begin, backup_end);

          // Bias due to bond angle prior
          bias += bond_angle_prior.get_log_bias(begin, end, backup_begin, backup_end);

          return bias;
     }

     //! Calculate the bias (implicit energy) introduced by the prior.
     //! \param dofs Vector of degrees of freedom
     //! \param dofs_backup Vector of degrees of freedom - backup chain
     //! \return log-probability
     double get_log_bias(const std::vector<Dof> &dofs,
                         const std::vector<Dof> &dofs_backup) const {

          double bias = 0.0;

          // Bias due to dihedral prior
          bias += dihedral_prior.get_log_bias(dofs, dofs_backup);

          // Bias due to bond angle prior
          bias += bond_angle_prior.get_log_bias(dofs, dofs_backup);

          return bias;
     }
          
     //! Accept last move
     void accept() {
          // Call accept in dihedral prior
          dihedral_prior.accept();

          // Call accept in bond angle prior
          bond_angle_prior.accept();
     }

          
     //! Reject last move
     void reject() {
          // Call reject in dihedral prior
          dihedral_prior.reject();

          // Call reject in bond angle prior
          bond_angle_prior.reject();
     }                    

     //! Defines the type of move info for a move using this prior
     //! \param modified_angle_range [start_index,end_index] pair specifying for which residues backbone angles have changed
     //! \param modified_position_range [start_index,end_index] pair specifying for which residues positions have changed
     //! \param move_type type of move
     MoveInfo *create_move_info(const std::pair<int,int> &modified_angle_range, const std::pair<int,int> &modified_position_range, definitions::MoveTypeEnum move_type) {
          return (new MoveInfo(move_type))->add_info(modified_angle_range,
                                                     modified_position_range);
     }

     //! Specifies how a move using this prior reacts to the knowledge that another moves has been executed
     //! Default case: do nothing
     //! \param move_info Object with information about last
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void notify(MoveInfo *move_info, const CHAIN_TYPE &chain) {
     }

     // Whether a particular dof is supported by this prior
     bool dof_supported(const DofIterator<ChainFB>::AngleSelectionEnum &dof_type) const {
          return dihedral_prior.dof_supported(dof_type) && bond_angle_prior.dof_supported(dof_type);
     }

};


     
//! Move prior (DBN). Consists of a bond angle and a dihedral prior
template <typename BONDANGLE_PRIOR_TYPE=BondAnglePriorEnghHuber, typename DBN_TYPE=TorusDBN>
class MovePriorDbn: public MovePrior<BONDANGLE_PRIOR_TYPE, DihedralPriorDbn<DBN_TYPE> > {
public:

     //! Whether this prior induces a bias
     const static bool induces_bias = true;
          
     //! Whether this prior has off-diagonal elements
     const static bool off_diagonal_elements = true;

     //! Local Settings settings
     const class Settings: public MovePrior<BONDANGLE_PRIOR_TYPE,
                                            DihedralPriorDbn<DBN_TYPE> >::Settings {} settings; //!< Local settings object

     //! Constructor
     //! \param dbn Pointer to dynamic Bayesian network model.
     //! \param settings Local settings object
     MovePriorDbn(DBN_TYPE *dbn, const Settings &settings):
          MovePrior<BONDANGLE_PRIOR_TYPE, DihedralPriorDbn<DBN_TYPE> >(dbn, settings),
          settings(settings) {}

     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     MovePriorDbn(const MovePriorDbn &other, int start_index, int end_index):
          MovePrior<BONDANGLE_PRIOR_TYPE, DihedralPriorDbn<DBN_TYPE> >(other, start_index, end_index),
          settings(other.settings) {}

     //! Copy Constructor - thread specified (multi-threading)
     //! \param other Source object from which copy is made
     //! \param thread_index Which thread the copy will run in
     MovePriorDbn(const MovePriorDbn &other, int thread_index):
          MovePrior<BONDANGLE_PRIOR_TYPE, DihedralPriorDbn<DBN_TYPE> >(other, thread_index),
          settings(other.settings) {}

     //! Defines the type of move info for a move using this prior
     //! \param modified_angle_range [start_index,end_index] pair specifying for which residues backbone angles have changed
     //! \param modified_position_range [start_index,end_index] pair specifying for which residues positions have changed
     //! \param move_type type of move
     MoveInfo *create_move_info(const std::pair<int,int> &modified_angle_range, const std::pair<int,int> &modified_position_range, definitions::MoveTypeEnum move_type) {
          return (new MoveInfoDbn<DBN_TYPE>(this->dihedral_prior.dbn, move_type))->add_info(modified_angle_range,
                                                                                            modified_position_range);
     }

     //! Specifies how a move using this prior reacts to the knowledge that another moves has been executed
     //! DBN case: notify DBN
     //! \param move_info Object with information about last
     //! \param chain Molecule chain object
     template <typename CHAIN_TYPE>
     void notify(MoveInfo *move_info, const CHAIN_TYPE &chain) {
          if (!move_info || move_info->modified_angles.empty())
               return;
          if (move_info->move_type != definitions::SIDECHAIN) {
               MoveInfoDbn<DBN_TYPE> *move_info_dbn = dynamic_cast<MoveInfoDbn<DBN_TYPE> *>(move_info);
               if (!move_info_dbn) {

                    std::vector<std::vector<double> > angles;
                    bool include_omega = true;
                    for (int i=move_info->modified_angles_start; i<move_info->modified_angles_end; ++i) {
                         angles.push_back((*chain)[i].get_angles(include_omega));
                    }

                    this->dihedral_prior.dbn->template set_sequence_vector<typename DBN_TYPE::ALL_ANGLE_NODES>(angles, move_info->modified_angles_start);

                    this->dihedral_prior.dbn->register_inconsistency(move_info->modified_angles_start,
                                                                     move_info->modified_angles_end);
               }
          }               
     }
};

}

#endif
