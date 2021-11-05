// move_crankshaft.h --- Crankshaft move: Rotate atoms around an axis defined by two Calpha atoms
//                       References: [1] Betancourt (2005) Efficient monte carlo trial moves for 
//                                       polypeptide simulations, J Chem Phys, 123, 174905
//                                   [2] Smith, Kortemme (2008) Backrub-like backbone simulation
//                                       recapitulates natural protein conformational variability and 
//                                       improves mutatant side-chain prediction, JMB, 380, 742-756
// Copyright (C) 2006-2010 Wouter Boomsma & Sandro Bottaro
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


#ifndef MOVE_CRANKSHAFT_H
#define MOVE_CRANKSHAFT_H

#include "move_priors.h"

namespace phaistos {

//! Crankshaft (backrub) move
template <typename CHAIN_TYPE>
class MoveCrankshaft: public MoveCommon<MoveCrankshaft<CHAIN_TYPE>, CHAIN_TYPE> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveCrankshaft<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Make moveStatistics a friend                                                                                                                                                                             
     friend class MoveStatistics_localmove<CHAIN_TYPE, MoveCrankshaft>;

     //! Local enum: move regions
     definitions::TerminalEnum move_region;   // INTERNAL, CTERM, NTERM

     //! Rotation angle around axis
     double current_rotation_angle;

     //! Size of allowed rotation interval
     double length_rotation_interval;

     //! Functionality associated with an end move
     class EndMove {
     public:

          //@{
          //! Degree-of-freedom iterators
          DofIterator<ChainFB> *begin;
          DofIterator<ChainFB> *end;
          DofIterator<ChainFB> *backup_begin;
          DofIterator<ChainFB> *backup_end;     
          //@}


          //! Number of degrees-of-freedom in this endmove
          int dof_count;

          //! Length of move (number of residues)
          int move_length;

          //! Use uninformative move prior
          MovePrior<BondAnglePriorUninformative, 
                    DihedralPriorUninformative> prior;

          //! Random number engine from which random number generators can be created.
          RandomNumberEngine *random_number_engine;

          //! Random number generator - uniform
          boost::variate_generator<RandomNumberEngine&, 
                                   boost::uniform_real<> > random_generator_uniform_01;

          //! Random number generator - Gaussian
          boost::variate_generator<RandomNumberEngine&, 
                                   boost::normal_distribution<> > random_generator_gaussian;

          //! Use settings object from outer class
          const typename MoveCrankshaft<CHAIN_TYPE>::Settings &settings;

          //! Constructor
          //! \param settings Settings object
          //! \param chain Molecule chain
          //! \param move_region whether we are at the N or the C terminal
          //! \param start_index Start index
          //! \param start_index End index
          //! \param constrain_bond_angle If true, only dihedrals are sampled
          //! \param random_number_engine Object from which random number generators can be created.
          EndMove(const typename MoveCrankshaft<CHAIN_TYPE>::Settings &settings,
                  CHAIN_TYPE *chain, CHAIN_TYPE *chain_backup,
                  definitions::TerminalEnum move_region, int start_index, int end_index,
                  bool constrain_bond_angle,
                  RandomNumberEngine *random_number_engine=&random_global)
               : move_length(end_index-start_index),
                 random_number_engine(random_number_engine),
                 random_generator_uniform_01(*random_number_engine, 
                                             boost::uniform_real<>(0,1)),
                 random_generator_gaussian(*random_number_engine, 
                                           boost::normal_distribution<>()),
                 settings(settings) {


               // Import protein definitions (such as residue names)
               using namespace definitions;

               // If the the angle bond variation is constrained the endmoves 
               // resample only dihedral angles, to avoid the bond angles to assume
               // values outside the allowed range
               DofIterator<ChainFB>::AngleSelectionEnum begin_dofs        = DofIterator<ChainFB>::DIHEDRAL_DOFS;
               DofIterator<ChainFB>::AngleSelectionEnum end_dofs          = DofIterator<ChainFB>::DIHEDRAL_DOFS;
               DofIterator<ChainFB>::AngleSelectionEnum backup_begin_dofs = DofIterator<ChainFB>::DIHEDRAL_DOFS;
               DofIterator<ChainFB>::AngleSelectionEnum backup_end_dofs   = DofIterator<ChainFB>::DIHEDRAL_DOFS;
               AtomEnum end_atom = CA;
               AngleEnum end_dof_type = definitions::DIHEDRAL;

               if(!constrain_bond_angle){
                    begin_dofs        = DofIterator<ChainFB>::STANDARD_DOFS;
                    end_dofs          = DofIterator<ChainFB>::STANDARD_DOFS;
                    backup_begin_dofs = DofIterator<ChainFB>::STANDARD_DOFS;
                    backup_end_dofs   = DofIterator<ChainFB>::STANDARD_DOFS;
                    end_dof_type = definitions::ANGLE;
               }


               if (move_region == NTERM) {
                    begin_dofs        += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
                    backup_begin_dofs += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
                    
               } else if (move_region == CTERM) {
                    begin_dofs        += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;
                    backup_begin_dofs +=(DofIterator<ChainFB>::CTERM_C_DIHEDRAL +
                                         DofIterator<ChainFB>::NTERM_CA_DIHEDRAL);
                    end_dofs          += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;
                    backup_end_dofs   += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;
                    end_atom = C;
                    if(!constrain_bond_angle){
                         begin_dofs += DofIterator<ChainFB>::CTERM_C_ANGLE;
                         backup_begin_dofs += DofIterator<ChainFB>::CTERM_C_ANGLE;
                         end_dofs += DofIterator<ChainFB>::CTERM_C_ANGLE;
                         backup_end_dofs += DofIterator<ChainFB>::CTERM_C_ANGLE;

                    }
               }

               // Construct Degree-of-freedom iterators
               begin        = new DofIterator<ChainFB>((*chain)(start_index, CA),
                                                       definitions::DIHEDRAL,
                                                       begin_dofs);
               backup_begin = new DofIterator<ChainFB>((*chain_backup)(0, CA),
                                                       definitions::DIHEDRAL,
                                                       backup_begin_dofs);
               end          = new DofIterator<ChainFB>((*chain)(end_index-1, end_atom),
                                                       end_dof_type,
                                                       end_dofs);
               backup_end   = new DofIterator<ChainFB>((*chain_backup)(end_index-1-start_index,
                                                                       end_atom),
                                                       end_dof_type,
                                                       backup_end_dofs);
                              
               // Calculate number of degrees of freedom in prerotation
               dof_count = 0;
               for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {
                    dof_count++;               
               }

               Matrix Lt;
               Matrix W_tilde;
               Matrix mu_tilde;
               construct_move(Lt, W_tilde, mu_tilde);
               
               // Generate Gaussian random values
               Matrix dchi;
               generate_dchi(dchi);

               // Determine changes to angular degrees of freedom
               Matrix dphi;
               dphi = solve_ax_b_triangular(Lt, dchi, 'U');
        
               // Shift dphi based on mu_tilde
               for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {
                    dphi(it.offset,0) += mu_tilde(it.offset,0);
               
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

               // Update chain positions
               int forced_direction = (move_region==NTERM)?-1:1;
               chain->update_positions_backbone(begin->get_residue()->index,
                                                end->get_residue()->index+1,
                                                forced_direction);               

               // Calculate new positions non backbone atoms
               for (int i=start_index; i<end_index; i++) {
                    (*chain)[i].update_positions_non_backbone();
               }
          }

          //! Destructor
          ~EndMove() {
               delete begin;
               delete backup_begin;
               delete end;
               delete backup_end;
          }

          //! Call the relevant methods for construction of an end move
          //! Matrices Lt, W_tilde and mu_tilde are calculated by this method.
          void construct_move(Matrix &Lt, Matrix &W_tilde, Matrix &mu_tilde) const {

               W_tilde = Matrix(dof_count, dof_count, 0.0);
               mu_tilde = Matrix(dof_count, 1);

               // Construct prior covariance matrix
               prior.construct_multivariate_gaussian(W_tilde, mu_tilde, *begin, *end);
               prior.add_constraint(W_tilde, mu_tilde, *begin, *end, 
                                    settings.endmove_std_dev_bond_angle,
                                    settings.endmove_std_dev_phi_psi,
                                    settings.endmove_std_dev_omega);

               // Calculate "square root" of matrix W_tilde
               Lt = W_tilde.cholesky('U');
          }

          //! Generate a vector of random normally distributed samples
          //! \param dchi Matrix in which to store samples 
          void generate_dchi(Matrix &dchi) {

               int size = dof_count;

               dchi = Matrix(size, 1);
               for (int i=0; i<size; i++) {
                    dchi(i,0) = random_generator_gaussian();
               }
          }

          //! Calculate the bias (implicit energy ratio) that was introduced by the prerotation          
          //! \return log probability
          double get_log_bias() const {

               double bias = 0.0;
               return bias;
          }

     };


     //! Pointer to end_move object
     EndMove *end_move;
     
public:

     //! Local settings class.
     const class Settings: public Move<CHAIN_TYPE>::Settings {
     public:

          //! Whether only to do internal moves (no moves at the N and C terminus)
          bool only_internal_moves;

          //! Whether generate rotations such that the new values for bond angles belong to the interval ideal_bond_angle +/- bond_angle_tolerance
          bool constrain_bond_angle;

          //! Tolerance of bond angles - internal moves
          double bond_angle_tolerance;

          //! Ideal value of bond angles - internal moves
          double ideal_bond_angle;

          //! Allowed variations of bond angles - end moves
          double endmove_std_dev_bond_angle;

          //! Allowed variations of phi/psi angles - end moves
          double endmove_std_dev_phi_psi;

          //! Allowed variations of omega angles - end moves
          double endmove_std_dev_omega;          

          //! Constructor
          Settings(bool only_internal_moves=false, 
                   bool constrain_bond_angle=false, 
                   double bond_angle_tolerance =15.0,
                   double ideal_bond_angle = 111.1,
                   double endmove_std_dev_bond_angle=0.8,
                   double endmove_std_dev_phi_psi=4.0,
                   double endmove_std_dev_omega=0.8)
               : Move<CHAIN_TYPE>::Settings(2,12),
                 only_internal_moves(only_internal_moves),
                 constrain_bond_angle(constrain_bond_angle),
                 bond_angle_tolerance(bond_angle_tolerance),
                 ideal_bond_angle(ideal_bond_angle),
                 endmove_std_dev_bond_angle(endmove_std_dev_bond_angle),
                 endmove_std_dev_phi_psi(endmove_std_dev_phi_psi),
                 endmove_std_dev_omega(endmove_std_dev_omega) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "only-internal-moves:" << settings.only_internal_moves << "\n";
               o << "constrain-bond-angle:" << settings.constrain_bond_angle << "\n";
               o << "bond-angle-tolerance:" << settings.bond_angle_tolerance << "\n";
               o << "ideal-bond-angle: " << settings.ideal_bond_angle << "\n";
               o << "endmove-std-dev-bond-angle:" << settings.endmove_std_dev_bond_angle << "\n";
               o << "endmove-std-dev-phi_psi:" << settings.endmove_std_dev_phi_psi << "\n";
               o << "endmove-std-dev-omega:" << settings.endmove_std_dev_omega << "\n";
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                    
     } settings;    //!< Local settings object

     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveCrankshaft(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "crankshaft", settings, random_number_engine),
            end_move(NULL),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveCrankshaft(const MoveCrankshaft &other)
          : MoveCommon(other),
            end_move(NULL),
            settings(other.settings) {
     }


     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveCrankshaft(const MoveCrankshaft &other, 
                    RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            end_move(NULL),
            settings(other.settings) {
     }

     //! Apply move
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether move was succesful
     bool apply(int start_index=-1, int end_index=-1) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          bool success = true;
          // Call base class apply method
          Move<CHAIN_TYPE>::apply(start_index, end_index);

          // Move range from -(length-2) to size+length-2
          if (start_index == end_index) {
               const std::pair<int,int> &region = this->random_move_region();
               int length = this->random_move_length(region);
               if (settings.only_internal_moves) {
                    this->random_move_range(region, length, 1, this->chain->size()-1, &start_index, &end_index);
               } else {
                    this->random_move_range(region, length,
                                            -(length-2), this->chain->size()+length-2,
                                            &start_index, &end_index);
               }
          }

          if (settings.debug) {
               std::cout << start_index << " " << end_index << "\n";
          }

	  // Determine whether to do loop closure or not
	  move_region = definitions::INTERNAL;
	  if (start_index<1) {
               move_region = definitions::NTERM;
	       start_index = 0;
	  }
	  if (end_index > this->chain->size()-1) {
               move_region = definitions::CTERM;
	       end_index = this->chain->size();
	  }

	  // Remember which angles are modified
          this->move_info = (new MoveInfo(LOCAL))->add_info(std::make_pair(start_index, end_index),
                                                            std::make_pair(start_index, end_index));

	  // Create temporary chain (backup in case of reject)
	  if (this->chain_backup) 
	       delete this->chain_backup;
	  this->chain_backup = new CHAIN_TYPE(*this->chain, 
                                              start_index, end_index);

          if (move_region == INTERNAL) {

               // Reject the move if the starting or ending residue is proline
               typename CHAIN_TYPE::Residue &res_start = (*this->chain)[start_index];
               typename CHAIN_TYPE::Residue &res_end = (*this->chain)[end_index];
               if ((res_start.residue_type == PRO) || (res_end.residue_type == PRO))
                    return (this->move_info->success = false);

               // Define relevant variables
 
               Atom *begin_atom = (*this->chain)(start_index,CA);
               Atom *begin_atom_neighbour = begin_atom->get_neighbour(+1, BACKBONE);
               Atom *begin_atom_previous = begin_atom->get_neighbour(-1, BACKBONE);
               Atom *end_atom = (*this->chain)(end_index-1,CA);
               Atom *end_atom_neighbour = end_atom->get_neighbour(+1, BACKBONE);
               Atom *end_atom_previous = end_atom->get_neighbour(-1, BACKBONE);
               
               int move_length = end_index - start_index;

               Vector_3D origin = begin_atom->position;
               Vector_3D v_1_begin = (begin_atom->position - begin_atom_previous->position).normalize();
               Vector_3D v_2_begin = (begin_atom_neighbour->position - begin_atom->position).normalize();
               Vector_3D v_1_end = (end_atom->position - end_atom_previous->position).normalize();
               Vector_3D v_2_end = (end_atom_neighbour->position - end_atom->position).normalize();
               

               // Define rotation axis
               Vector_3D axis = (end_atom->position - begin_atom->position).normalize();
               
               // vector containining minimum and maximum value for rotation
               std::vector<double> range;

               // return the maximum allowed rotation
               success = get_rotation_size(&range, axis, v_1_begin, v_2_begin, v_1_end, v_2_end, move_length);
               if (!success) {
                    return (this->move_info->success = false);
               }
               
               // Generate random rotation in range
               this->length_rotation_interval = std::fabs(range[0] - range[1]);
               double random_number = this->random_generator_uniform_01();

               current_rotation_angle = (random_number)*(std::fabs(range[0] - range[1])) - std::fabs(range[0]);


               // Perform the rotation
               Matrix_3D rotation_matrix = rotation_matrix_from_dihedral(axis,current_rotation_angle);

               // Create iterator
               AtomIterator<CHAIN_TYPE,ALL> it(*this->chain, 
                                               start_index, CA,
                                               end_index-1, CA);

               // Move all atoms in range
               for (++it; !it.end(); ++it) {
                    Vector_3D new_position = it->position;
                    new_position -= origin;
                    new_position = rotation_matrix*new_position;
                    new_position += origin;
                    it->position = new_position;
               }

               // Update degrees of freedom
               begin_atom->update_angle();
               begin_atom->update_dihedral();
               begin_atom_neighbour->update_dihedral();
               end_atom->update_angle();
               end_atom->update_dihedral();
               if (end_atom_neighbour) {
                    end_atom_neighbour->update_dihedral();
               }

               // The sidechain at the endpoints are effected
               if ((*this->chain)(start_index).has_sidechain() || (*this->chain)(start_index).has_atom(PS)) {
                    (*this->chain)(start_index).update_positions_non_backbone();
               }
               if ((*this->chain)(end_index-1).has_sidechain() || (*this->chain)(start_index).has_atom(PS)) {
                    (*this->chain)(end_index-1).update_positions_non_backbone();
               }


               // calculate bond angle after rotation (for testing purposes)
               if (settings.constrain_bond_angle){
                    Vector_3D v_begin_new = (begin_atom->position - begin_atom_previous->position).normalize();
                    Vector_3D u_begin_new = (begin_atom_neighbour->position - begin_atom->position).normalize();
                    Vector_3D v_end_new = (end_atom->position - end_atom_previous->position).normalize();
                    Vector_3D u_end_new = (end_atom_neighbour->position - end_atom->position).normalize();
                    double alpha_new_begin = (180.0/M_PI)*acos(-u_begin_new*v_begin_new);
                    double alpha_new_end = (180.0/M_PI)*acos(-u_end_new*v_end_new);
//                    double alpha_begin = (180.0/M_PI)*acos(-v_1_begin*v_2_begin);
//                    double alpha_end = (180.0/M_PI)*acos(-v_1_end*v_2_end);

                    if ((alpha_new_begin > settings.ideal_bond_angle+ settings.bond_angle_tolerance) || (alpha_new_begin < settings.ideal_bond_angle - settings.bond_angle_tolerance) ||
                        (alpha_new_end > settings.ideal_bond_angle + settings.bond_angle_tolerance) || (alpha_new_end < settings.ideal_bond_angle - settings.bond_angle_tolerance)){
//                         std::cerr << " WARNING - I am proposing a crankshaft move that violates the angle constraint - the move will be rejected " <<"\n";
//                         std::cerr << " boundaries : "<<  settings.ideal_bond_angle-  settings.bond_angle_tolerance << " " <<  settings.ideal_bond_angle+ settings.bond_angle_tolerance << "\n";
//                    std::cerr << " Proposed angle (start): " << alpha_new_begin << " old angle: " << alpha_begin<<" Proposed angle (end): "<< alpha_new_end << " old angle: " << alpha_end<<" " << current_rotation_angle <<"\n";                
                         return (this->move_info->success = false);
                    }
               }
               
          } else {
               end_move = new EndMove(settings, 
                                      this->chain, this->chain_backup,
                                      move_region, start_index, end_index,
                                      settings.constrain_bond_angle,
                                      this->random_number_engine);
          }

          return success;
     }


     //! Determine the interval from which the rotation angle tau is drawn 
     bool get_rotation_size(std::vector<double> *range, Vector_3D axis, 
                            Vector_3D v_1_begin, Vector_3D v_2_begin,
                            Vector_3D v_1_end,Vector_3D v_2_end, int move_size) {
          

          // Set tolerance on bond angle variation 
          double tolerance_bond_angle = settings.bond_angle_tolerance/180.0*M_PI;
                  
          //calculate maximum allowed rotation (Backrub-style)
          double backrub_constraint = calculate_tau_global(move_size);
          double minimum = -1.0*backrub_constraint;
          double maximum = backrub_constraint;

          if(settings.constrain_bond_angle){
               
               // Define vectors containing intervals
               std::vector<double> tau_1;
               std::vector<double> tau_2;

               // The fixed vector is passed to the function before the rotating vector
               calculate_tau_bond_angle(axis, v_1_begin, v_2_begin, tolerance_bond_angle,&tau_1);
               calculate_tau_bond_angle(axis, v_2_end, v_1_end, tolerance_bond_angle,&tau_2);

               double lower_bound = -M_PI/2.0;
               double upper_bound = M_PI/2.0;

               // find the intersection of the bond angle intervals
                if((tau_1[1] < tau_2[0]) || tau_2[1] < tau_1[0]){
                    return false;

               } else {
                    if (tau_2[0]> tau_1[0])
                         lower_bound = tau_2[0];
                    else 
                         lower_bound = tau_1[0];

                    if (tau_2[1] < tau_1[1])
                         upper_bound = tau_2[1];
                    else 
                         upper_bound = tau_1[1];
               }

               // intersection with the backrub interval
               if((upper_bound < minimum) || lower_bound > maximum){
                    return false;

               } else {
                    if(upper_bound < maximum)
                         maximum = upper_bound;
                    if(lower_bound > minimum)
                         minimum = lower_bound;
               }
               if((maximum-minimum) < 1.0e-15){
                    return false;}
          }

      
          range->push_back(minimum);
          range->push_back(maximum);
          return true;
          
     }

     //! Calculate the maximum rotation angle tau based on [2] (see header)
     double calculate_tau_global(int move_size){

          // Backrub formula
          double tau = 23.0 - (double)move_size;
          if(move_size == 2){
               tau = 40.0;
          }
          tau = (tau/180.0)*M_PI;
          return tau;
     }


     //! Calculate the rotation interval for tau that satisfies the restriction on bond angles 
     bool calculate_tau_bond_angle(Vector_3D axis, Vector_3D fixed_vector, 
                                  Vector_3D rotating_vector, double tolerance, 
                                  std::vector<double> *vec){

          // define relevant quantities (Betancourt convention)
          Vector_3D v = fixed_vector;
          Vector_3D u = rotating_vector;
          
          double cos_sigma_u = u*axis;
          double cos_sigma_v = v*axis;
          double sin_sigma_v_sq = 1.0 - (cos_sigma_v*cos_sigma_v);
          double sin_sigma_u_sq = 1.0 - (cos_sigma_u*cos_sigma_u);          

          // Avoid numerical issues with sqrt
          if (std::fabs(sin_sigma_v_sq < 1.0e-04) || std::fabs(sin_sigma_v_sq < 1.0e-04) ){
               vec->push_back(-M_PI/2.0);
               vec->push_back(M_PI/2.0);
               return true;
          }
          double sin_sigma_u = sqrt(sin_sigma_u_sq);
          double sin_sigma_v = sqrt(sin_sigma_v_sq);

          long double cos_alpha  = -1.0*u*v;
          
          // set range of allowed bond angles (exclude boundaries)
          double alpha_ideal = settings.ideal_bond_angle*(M_PI/180.0);
          long double max_value = alpha_ideal + tolerance;
          long double min_value = alpha_ideal - tolerance;
          long double cos_alpha_new_plus =  cos(max_value);
          long double cos_alpha_new_minus =  cos(min_value);
          
          if(acos(cos_alpha) > max_value  || acos(cos_alpha) < min_value ){
               std::cerr << "WARNING: angle bend value out of range - current value: " << acos(cos_alpha)*(180.0/M_PI) << "\n";
               std::cerr << "Should be between: " << (180.0/M_PI)* max_value << " and: " <<  min_value*(180.0/M_PI) << "\n";
               std::cerr << "If you are using only crankshaft moves this angle will never be resampled." << "\n";

               return false;
          }
          
          double sign = (v*(axis%u))/std::fabs(v*(axis%u));
          long double arg = (cos_alpha +(cos_sigma_u*cos_sigma_v))/(sin_sigma_u*sin_sigma_v);               
          long double tau_u = sign*acos(arg);

          if ((u*v - (cos_sigma_u*cos_sigma_v - (sin_sigma_u*sin_sigma_v)*cos(tau_u))) > 0.01)          {
               std::cerr << " ERROR: Sanity check in Cranshaft move failed " << "\n";
          }
          
          long double arg_new_plus = (cos_alpha_new_plus +(cos_sigma_u*cos_sigma_v))/(sin_sigma_u *sin_sigma_v);
          long double arg_new_minus = (cos_alpha_new_minus +(cos_sigma_u*cos_sigma_v))/(sin_sigma_u *sin_sigma_v);
          double tau_lower = -M_PI/2.0;
          double tau_upper = M_PI/2.0;

          if((std::fabs(arg_new_plus) < 1.0) && (std::fabs(arg_new_plus) > 0.0) && (std::fabs(arg_new_minus) < 1.0) && (std::fabs(arg_new_minus) > 0.0)){
               long double tau_1 = sign*acos(arg_new_plus) - tau_u;
               long double tau_2 = sign*acos(arg_new_minus) - tau_u;
               if(tau_1 < tau_2){
                    tau_upper = tau_2;
                    tau_lower = tau_1;
               } else {
                    tau_upper = tau_1;
                    tau_lower = tau_2;
               }
          } else {
               if((std::fabs(arg_new_plus) < 1.0) && (std::fabs(arg_new_plus) > 0.0)){
                    long double tau_1 = sign*acos(arg_new_plus) - tau_u;
                    if(tau_1 < 0.0)
                         tau_lower = tau_1;
                    else 
                         tau_upper = tau_1;
               }
               if((std::fabs(arg_new_minus) < 1.0) && (std::fabs(arg_new_minus) > 0.0)){
                    long double tau_1 = sign*acos(arg_new_minus) - tau_u;
                    if(tau_1 < 0.0)
                         tau_lower = tau_1;
                    else 
                         tau_upper = tau_1;
               }

          }
          if ((std::fabs(tau_lower) < 1.0e-03))
               tau_lower = 0.0;
          if ((std::fabs(tau_upper) < 1.0e-03))
               tau_upper = 0.0;

          vec->push_back(tau_lower);
          vec->push_back(tau_upper);
          
          return true;
     }
     
     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          if (end_move) {
               return end_move->get_log_bias();
          }else {
               if(settings.constrain_bond_angle){
                    int start_index = this->move_info->modified_positions[0].first;
                    int end_index = this->move_info->modified_positions[0].second;
                    Atom *begin_atom = (*this->chain)(start_index,CA);
                    Atom *begin_atom_neighbour = begin_atom->get_neighbour(+1, ALL);
                    Atom *begin_atom_previous = begin_atom->get_neighbour(-1, ALL);
                    Atom *end_atom = (*this->chain)(end_index-1,CA);
                    Atom *end_atom_neighbour = end_atom->get_neighbour(+1, ALL);
                    Atom *end_atom_previous = end_atom->get_neighbour(-1, ALL);
                    int move_length = end_index- start_index;
                    Vector_3D v_begin = (begin_atom->position - begin_atom_previous->position).normalize();
                    Vector_3D u_begin = (begin_atom_neighbour->position - begin_atom->position).normalize();
                    Vector_3D v_end = (end_atom->position - end_atom_previous->position).normalize();
                    Vector_3D u_end = (end_atom_neighbour->position - end_atom->position).normalize();
                    Vector_3D axis = (begin_atom->position - end_atom->position).normalize();
                    std::vector<double> range;
                    if (!(get_rotation_size(&range, axis, u_begin, v_begin, u_end, v_end, move_length))){
                         this->move_info->success = false;
                         return 1.0e20;
                    }
                    double length_rotation_interval_prime = std::fabs(range[0] - range[1]);
                    return std::log(this->length_rotation_interval) - std::log(length_rotation_interval_prime);
               }
               else
                    return 0.0;
          }
     }

     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();
          
	  // std::cout << "accept\n";
	  delete this->chain_backup;
	  this->chain_backup = NULL;

          if (end_move) {
               delete end_move;
               end_move = NULL;
          }
     }

     //! Reject last move
     void reject() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();
                              
	  // Update positions and angles in the current chain based on the test chain
	  assert(this->chain_backup != NULL);

          AtomIterator<ChainFB,ALL> it_current(*this->chain, 
                                               this->move_info->modified_positions[0].first, N,
                                               this->move_info->modified_positions[0].second, CA);
          AtomIterator<ChainFB,ALL> it_backup(*this->chain_backup);
          for (; !it_current.end() && !it_backup.end(); ++it_current,++it_backup) {
	       it_current->set_angle(it_backup->get_angle());
	       it_current->set_dihedral(it_backup->get_dihedral());
	       it_current->position = it_backup->position;
	  }

          if (end_move) {
               delete end_move;
               end_move = NULL;
          }
     }     
};

}

#endif
