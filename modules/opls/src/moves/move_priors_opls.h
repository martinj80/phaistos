// move_priors_opls.h --- Priors used by various moves - OPLS version
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


#ifndef MOVE_PRIORS_OPLS
#define MOVE_PRIORS_OPLS

#include "energy/OPLS/OPLS_anglebend.h"

//! OPLS bond angle prior
class BondAnglePriorOpls {

     //! Calculate a single entry of the bias (implicit energy) introduced by the prior - internal method (see get_log_bias())
     //! \param atom Atom pointer
     //! \param dof_type Type of degree of freedom
     //! \return log-probability
     inline double calculate_log_bias_entry(Atom *atom,
                                            const definitions::AngleEnum &dof_type,
                                            const double &value) const {

          if (dof_type() == DofIterator<ChainFB>::ANGLE) {
               Atom *next_atom = atom->get_neighbour(+1, BACKBONE);
               Atom *prev_atom = atom->get_neighbour(-1, BACKBONE);

               std::pair<double,double> param = getParameters(prev_atom, atom, next_atom);

               double mean = param.first;
               double std_dev = param.second;
                              
               double angle = *it;

               // Gaussian density
               return -0.5*(sqr(angle-mean)/sqr(std_dev)) - log(std_dev*sqrt(2.0*M_PI));
          }
          return 0;
     }

     //! Calculate the bias (implicit energy) that was introduced by the prior - internal method (see get_log_bias())
     //! \param begin Degree-of-freedom iterator start
     //! \param end Degree-of-freedom iterator end
     //! \return log-probability     
     double get_log_bias_internal(const DofIterator<ChainFB> &begin, 
                                  const DofIterator<ChainFB> &end) const {

          double bias = 0.0;

          // If this prior is to be used as an implicit energy
          // we include the likelihood of the angles as a bias term
          // NOTE: THE WE ARE NOT DIVIDING OUT A PROPOSAL TERM HERE -
          // THIS IS ALREADY DONE BY THE PREROTATION ITSELF. WE ARE SIMPLY
          // PUTTING IN THE ENERGY CORRESPONDING TO THE PROPOSAL, TO ENSURE
          // THAT EVEN WHEN SIMULATING WITHOUT ENERGIES, THE BONDANGLE IS
          // KEPT WITHIN REASONABLE BOUNDS
          // THIS IS TO MAKE VALID THE ASSUMPTION THAT PER DEFAULT IN PHAISTOS
          // WHEN SAMPLING WITHOUT ENERGIES, DIHEDRALS ARE DISTRIBUTED ACCORDING
          // TO BACKBONEDBN AND BONDANGLES ACCORDING TO ENGH-HUBER
          if (settings.implicit_bondangle_energy) {
                    
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

     //! Whether this prior induces a bias
     const static bool induces_bias = true;

     //! Whether this prior has off-diagonal elements
     const static bool off_diagonal_elements = false;

     //! Local Settings class
     const class Settings {
     public:

          //! Whether to use this prior as an implicit energy. If false
          //! the bias introduced by the prior will be "divided out" in
          //! the acceptance ratio
          bool implicit_bondangle_energy;

          //! Constructor
          Settings(bool implicit_bondangle_energy=true,
                   double max_bondangle=0.8, 
                   double postrotational_scale_bond=1.0):
               implicit_bondangle_energy(implicit_bondangle_energy) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               std::cout << "implicit_bondangle_energy: " << settings.implicit_bondangle_energy << "\n";
               return o;
          }               
     } settings;    //!< Local settings object

          
     //! String identifier
     std::string id;

          
     //! Tinker Parameter Lookup
     AnglebendParameters *oplsAngleParameters;

     //! Whether this object owns the Tinker parameters
     bool oplsAngleParametersOwner;

     //! Constructor
     //! \param settings Local settings object
     BondAnglePriorOpls(const Settings &settings):
          settings(settings), id("OPLS") {
          oplsAngleParameters = new AnglebendParameters();
          oplsAngleParametersOwner = true;
     }
          
     //! Copy constructor - factory method
     //! \param other Source object from which copy is made
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     BondAnglePriorOpls(const BondAnglePriorOpls &other, int start_index, int end_index):
          settings(other.settings), id(other.id),
          oplsAngleParameters(other.oplsAngleParameters),
          oplsAngleParametersOwner(false) {}

     //! Destructor
     ~BondAnglePriorOpls() {
          if (oplsAngleParametersOwner)
               delete oplsAngleParameters;
     }
          
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
               if (it.get_dof_type() == DofIterator<ChainFB>::ANGLE) {

                    // Get mean and stddev
                    Atom *atom = it.get_atom();
                    Atom *next_atom = atom->get_neighbour(+1, BACKBONE);
                    Atom *prev_atom = atom->get_neighbour(-1, BACKBONE);
                    std::pair<double,double> param = get_parameters(prev_atom, atom, next_atom);
                    double mean = param.first;
                    double std_dev = param.second;

                    // Set covariance matrix entry
                    sigma(it.offset, it.offset) = 1.0/(sqr(std_dev));

                    // Set mean entry
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
                    
               double total_constraint = sqr(180.0/M_PI)*(1.0/sqr(std_dev_bond_angle));	      
                    
               constraint_added = true;
               for (DofIterator<ChainFB> it=start_dof; it!=end_dof; ++it) {
                    if (it.get_dof_type() == DofIterator<ChainFB>::ANGLE) {
		  
                         // Back-calculate std_dev from covariance matrix entry
                         // (NOTE: this assumes that add constraint is called on
                         // an unmodified covariance matrix)
                         double std_dev = sqrt(1.0/sigma(it.offset, it.offset));
		  
                         if (std_dev_bond_angle < std_dev*(180/M_PI)) {
                              sigma(it.offset, it.offset) = total_constraint;
                              // sqr(180.0/M_PI)*(1.0/sqr(settings.std_dev_bond_angle));
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
	      
          AnglebendParameters::Parameter param =
               oplsAngleParameters->get(atom1, atom2, atom3);
          return std::make_pair(param.eqAng, sqrt(0.6/(2.0*param.k)));
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

          return (get_log_bias_internal(begin, end) -
                  get_log_bias_internal(backup_begin, backup_end));
     }

          
     //! Accept last move
     void accept() {}

     //! Reject last move
     void reject() {}          
          
};

#endif
