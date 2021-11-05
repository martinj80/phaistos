// move_crisp.h --- Concerted Rotation Integrating Structural Priors
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


#include <utility>
#include "utils/matrix_vec.h"
#include "moves/move.h"
#include "moves/move_priors.h"
#include <typeinfo>

#ifndef MOVE_CRISP
#define MOVE_CRISP


namespace phaistos {

// Name space for CRISP internals
namespace crisp {

//////////////////////////////////////
//        PREROTATION OBJECTS       //
//////////////////////////////////////

//! Prerotation class
template <typename PRIOR_TYPE=MovePrior<> >
class PreRotation {

     //! Random number engine from which random number generators can be created.
     RandomNumberEngine *random_number_engine;

     //! Random number generator  - uniform
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > random_generator_uniform_01;
     //! Random number generator  - Gaussian
     boost::variate_generator<RandomNumberEngine&, 
                              boost::normal_distribution<> > random_generator_gaussian;

public:

     //! Local settings object
     const class Settings: public PRIOR_TYPE::Settings {
     public:

          //! Number of potential dofs modified by prerotation: -1=ALL
          unsigned int prerotation_active_dofs;     
          
          //! Parameter controlling the amplitude of bond angles variation
          double std_dev_bond_angle;                 

          //! Parameter controlling the amplitude of phi and psi  angles variation
          double std_dev_phi_psi;                   

          //! Parameter controlling the amplitude of omega angles variation
          double std_dev_omega;                     

          //! Constructor
          Settings(unsigned int prerotation_active_dofs=(unsigned int)-1,
                   double std_dev_bond_angle=0.8,
                   double std_dev_phi_psi=4.0,
                   double std_dev_omega=0.8):
               prerotation_active_dofs(prerotation_active_dofs),
               std_dev_bond_angle(std_dev_bond_angle),
               std_dev_phi_psi(std_dev_phi_psi),
               std_dev_omega(std_dev_omega){
          }
  
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "prerotation-active-dofs:" << settings.prerotation_active_dofs << "\n";
               o << "std-dev-bond-angle:" << settings.std_dev_bond_angle << "\n";
               o << "std-dev-phi_psi:" << settings.std_dev_phi_psi << "\n";
               o << "std-dev-omega:" << settings.std_dev_omega << "\n";
               o << static_cast<const typename PRIOR_TYPE::Settings &>(settings);
               return o;
          }                    

     } settings;    //!< Local settings object
     
     //! Prior Type
     PRIOR_TYPE prior;

     //! String identification of Prior
     std::string id;
     
     //! Move region (INTERNAL,CTERM,NTERM)
     definitions::TerminalEnum move_region;
     
     //! Number of degrees of freedom effected by prerotation
     unsigned int prerotation_dof_count;
     
     //! Number of degrees of freedom effected by postrotation (=6)
     const static int postrotation_dof_count = 6;

     //! Chain object
     ChainFB *chain;

     //! length of move (number of residues)
     int move_length;

     //@{
     //! Degree-of-freedom iterators          
     const DofIterator<ChainFB> &begin;
     const DofIterator<ChainFB> &breakpoint;
     const DofIterator<ChainFB> &end;
     const DofIterator<ChainFB> &backup_begin;
     const DofIterator<ChainFB> &backup_breakpoint;
     const DofIterator<ChainFB> &backup_end;

     DofIterator<ChainFB> begin_inactive;
     DofIterator<ChainFB> backup_begin_inactive;
     //@}

     //! W tilde matrix (saved here for get_log_bias)
     Matrix W_tilde;

     //! Choleski decomposition of W tilde
     Matrix Lt;

     //! Collection of random numbers
     Matrix dphi;
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param id String identification of Prior
     PreRotation(ChainFB *chain,
                 const Settings &settings,
                 RandomNumberEngine *random_number_engine=&random_global,
                 std::string id="PreRotation")
          : random_number_engine(random_number_engine),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings),
            prior(settings),
            id(id+"<"+prior.id+">"),
            chain(chain),
            begin((Atom *)NULL),
            breakpoint((Atom *)NULL),
            end((Atom *)NULL),
            backup_begin((Atom *)NULL),
            backup_breakpoint((Atom *)NULL),
            backup_end((Atom *)NULL),
            begin_inactive((Atom *)NULL),
            backup_begin_inactive((Atom *)NULL) {
     }

     
     //! Constructor for DBN-based preRotation
     //! \param chain Molecule chain
     //! \param dbn DBN object
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param id String identification of Prior
     template <typename DBN_TYPE>
     PreRotation(ChainFB *chain, DBN_TYPE *dbn, const Settings &settings,
                 RandomNumberEngine *random_number_engine=&random_global,
                 std::string id="PreRotation")
          : random_number_engine(random_number_engine),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings),
            prior(dbn, settings),
            id(id+"<"+prior.id+">"),
            chain(chain),
            begin((Atom *)NULL),
            breakpoint((Atom *)NULL),
            end((Atom *)NULL),
            backup_begin((Atom *)NULL),
            backup_breakpoint((Atom *)NULL),
            backup_end((Atom *)NULL),
            begin_inactive((Atom *)NULL),
            backup_begin_inactive((Atom *)NULL) {
     }

     
     //! Copy constructor, with specification of iterators
     //! This allows prerotation objects to be created from one initial
     //! prerotation object (the prerotationFactory)
     //! \param other Prerotation object from which the copy is made
     //! \param move_region Whether we are at the N or the C terminal
     //! \param move_length Length of move (number of residues)
     //! \param begin Begin dof iterator
     //! \param breakpoint Breakpoint dof iterator
     //! \param end End dof iterator
     //! \param backup_begin Begin dof iterator in backup chain 
     //! \param backup_breakpoint Breakpoint dof iterator in backup chain 
     //! \param backup_end End dof iterator in backup chain 
     PreRotation(const PreRotation &other, definitions::TerminalEnum move_region,
                 int move_length,
                 const DofIterator<ChainFB> &begin,
                 const DofIterator<ChainFB> &breakpoint,
                 const DofIterator<ChainFB> &end,
                 const DofIterator<ChainFB> &backup_begin,
                 const DofIterator<ChainFB> &backup_breakpoint,
                 const DofIterator<ChainFB> &backup_end)
          : random_number_engine(other.random_number_engine),
            random_generator_uniform_01(other.random_generator_uniform_01),
            random_generator_gaussian(other.random_generator_gaussian),
            settings(other.settings), 
            prior(other.prior, begin.get_residue()->index, end.get_residue()->index+1),
            id(other.id), move_region(move_region), chain(other.chain), move_length(move_length),
            begin(begin), breakpoint(breakpoint), end(end),
            backup_begin(backup_begin), backup_breakpoint(backup_breakpoint), backup_end(backup_end),
            begin_inactive(begin),
            backup_begin_inactive(backup_begin) {

          // If the prerotation_active_dofs parameter is set, 
          // the prerotation doesn't iterate all the way to breakpoint
          prerotation_dof_count = 0;

          if (move_region == definitions::INTERNAL) {
               for (; begin_inactive!=breakpoint && prerotation_dof_count < settings.prerotation_active_dofs; 
                    ++begin_inactive,++backup_begin_inactive,++prerotation_dof_count) {
                    //std::cout << begin_inactive << "\n";
               }
          } else {
               for (; begin_inactive!=breakpoint; ++begin_inactive,prerotation_dof_count++) {}               
          }
     }


     //! Copy constructor - thread specified (multi-threading)
     //! \param other Prerotation object from which the copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index of the thread
     //! \param chain Molecule chain
     PreRotation(const PreRotation &other, 
                 RandomNumberEngine *random_number_engine,
                 int thread_index,
                 ChainFB *chain)
          : random_number_engine(random_number_engine),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(other.settings), 
            prior(other.prior, thread_index),
            id(other.id), 
            chain(chain),
            move_length(other.move_length),
            begin((Atom *)NULL),
            breakpoint((Atom *)NULL),
            end((Atom *)NULL),
            backup_begin((Atom *)NULL),
            backup_breakpoint((Atom *)NULL),
            backup_end((Atom *)NULL),
            begin_inactive(begin),
            backup_begin_inactive(backup_begin) {
     }


     //! Destructor
     virtual ~PreRotation() {}

     
     //! Execute prerotation move
     //! \return Boolean indicating whether move was successful
     bool apply() {
          if (move_region == definitions::INTERNAL) {
               if(!apply_internalmove())
                    return false;
          } else {
               apply_endmove();
          }
          
          return true;
     }

     
     //! Execute internal move
     //! \return Boolean indicating whether it was possible to perform the move
     bool apply_internalmove() {

          // Construct Lt matrix
          if(!construct_all(this->Lt, this->W_tilde))
               return false;

          // Generate Gaussian random values
          Matrix dchi;
          generate_dchi(dchi);

          // Determine changes to angular degrees of freedom
          dphi = solve_ax_b_triangular(Lt, dchi, 'U');
          
          // Construct new conformation
          for (DofIterator<ChainFB> it=begin; it!=begin_inactive; ++it) {
               *it += dphi(it.offset,0);
               *it = fmod(*it+3*M_PI,(2*M_PI))-M_PI;

               // Take absolute value if current dof is an angle
               if (it.get_dof_type() == definitions::ANGLE) {
                    *it = std::fabs(*it);
               }
          }

          // Update chain positions
          chain->update_positions_backbone_segment(begin.get_residue()->index,
                                                   breakpoint.get_residue()->index+1);
          return true;
     }

     
     //! Execute constrained end-move 
     void apply_endmove() {

          // Construct Lt matrix
          construct_endmove(this->Lt, this->W_tilde);

          // Generate Gaussian random values
          Matrix dchi;
          generate_dchi(dchi);

          // Determine changes to angular degrees of freedom
          dphi = solve_ax_b_triangular(Lt, dchi, 'U');
   
          // Construct new conformation
          for (DofIterator<ChainFB> it=begin; it!=breakpoint; ++it) {
               *it += dphi(it.offset,0);
               *it = fmod(*it+3*M_PI,(2*M_PI))-M_PI;

               // Take absolute value if current dof is an angle
               if (it.get_dof_type() == definitions::ANGLE) {
                    *it = std::fabs(*it);
               }
          }

          // Update chain positions
          int forced_direction = (move_region==definitions::NTERM)?-1:1;
          chain->update_positions_backbone(begin.get_residue()->index,
                                         breakpoint.get_residue()->index+1,
                                         forced_direction);

     }

     
     //! Create and return move info object 
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \param move_type type of move
     MoveInfo *create_move_info(int start_index, int end_index, 
                                definitions::MoveTypeEnum move_type) {
          return prior.create_move_info(std::make_pair(start_index, end_index),std::make_pair(start_index, end_index), move_type);
     }

     //! Specifies how move reacts to the knowledge that another moves has been executed
     //! \param move_info Move information
     void notify(MoveInfo *move_info) {
          prior.notify(move_info, chain);
     }
     
     //! Update of internal state based on newly modified chain state
     //! \param chain Molecule chain
     template <typename CHAIN_TYPE>
     void update_state(CHAIN_TYPE *chain) {
          prior.update_state(chain);
     }

     //! Generate a vector of random normally distributed samples
     //! \param dchi List of random numbers
     void generate_dchi(Matrix &dchi) {

          int size = prerotation_dof_count;

          dchi = Matrix(size, 1);
          for (int i=0; i<size; i++) {
               dchi(i,0) = random_generator_gaussian();
          }
     }

     
     //! Calculate the bias (implicit energy ratio) that was introduced by the prerotation 
     //! \return Log probability
     double get_log_bias() const {

          // Basic move Bias
          double bias = 0.0;

          if (move_region == definitions::INTERNAL) {
               double log_prerotation_bias_forward = calculate_bias_forward();
               double log_prerotation_bias_reverse = calculate_bias_reverse();
               bias += log_prerotation_bias_reverse - log_prerotation_bias_forward;
          }
                                   
          // Bias from prior distribution
          bias += prior.get_log_bias(begin, end, backup_begin, backup_end);

          return bias;
     }

     
     //! Calculate bias for prerotation move (a->b)
     //! \return Log bias introduced by the proposed move
     double calculate_bias_forward() const {
          int size = prerotation_dof_count;

          double logdetL = 0.0;
          double sum = 0.0;
          for (int i=0; i<size; i++) {
               logdetL += std::log(Lt(i,i));
               sum += 0.5*(dphi(i,0)) * W_tilde(i,i) * (dphi(i,0));
               for (int j=i+1; j<size; j++) {
                    sum += (dphi(i,0)) * W_tilde(i,j) * (dphi(j,0));
               }
          }
          // Calculate biasing probability a->b
          return logdetL - sum;
     }

     //! Calculate bias for reverse move (b->a)
     //! \return Log bias introduced by the reverse move
     // NOTE: Assumes that apply has been called
     double calculate_bias_reverse() const {

          int size = prerotation_dof_count;

          // Calculate Lt'
          Matrix Lt_prime;
          Matrix W_tilde_prime;
          if (move_region == definitions::INTERNAL) {
               if(!construct_all(Lt_prime, W_tilde_prime))
                    return std::numeric_limits<double>::infinity();
          } else {
               construct_endmove(Lt_prime, W_tilde_prime);
          }
          
          double logdetL_prime = 0.0;
          double sum = 0.0;
          for (int i=0; i<size; i++) {
               logdetL_prime += std::log(Lt_prime(i,i));
               sum += 0.5*((-dphi(i,0))) *
                    W_tilde_prime(i,i) * ((-dphi(i,0)));
               for (int j=i+1; j<size; j++) {
                    sum += (-dphi(i,0)) *
                         W_tilde_prime(i,j) * (-dphi(j,0));
               }
          }

          // Calculate biasing probability b->a
          return logdetL_prime - sum;
     }

     //! Accept last move
     void accept() {
          // Call accept in prior
          prior.accept();
     }

     //! Reject last move
     void reject() {
          // Call reject in prior
          prior.reject();
     }
     

     //! Construct all the matrices for internal move
     //! \param Lt Choleski decomposition of W tilde
     //! \param W_tilde Matrix containing structural priors
     //! \return Boolean indicating whether the construction of the matrices and the decompoistion were successful 
     bool construct_all(Matrix &Lt, Matrix &W_tilde) const {
          
          int size = this->prerotation_dof_count;
          
          // Construct diagonal covariance matrix
          Matrix W(size, size, 0.0);
          
          double phi_psi_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_phi_psi));
          double omega_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_omega));
          double bond_angle_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_bond_angle));
                                                         
          for (DofIterator<ChainFB> it=begin; it!=begin_inactive; ++it) {
               if (it.get_dof_type() == definitions::DIHEDRAL) {
                    if(it.get_atom()->atom_type == definitions::N) {
                         W(it.offset, it.offset) += omega_constraint;
                    } else {
                         W(it.offset, it.offset) += phi_psi_constraint;
                    }
               } else {
                    W(it.offset, it.offset) += bond_angle_constraint;
               }
          }
          
          
          if(!construct_M_tilde(W))
               return false;

          W_tilde = W;         
                       
          // Calculate "square root" of matrix W_tilde
          try {
               Lt = W_tilde.cholesky('U');
          } catch (Matrix::MatrixError &e) {
               return false;
          }

          return true;
     }
     
     //! Construct M matrix, which combines the chain closure 
     //! information with the constraint on bond and omega angles
     //! \param W_tilde Matrix containing structural priors
     //! \return Boolean indicating whether he construction of the matrices was successful
     bool construct_M_tilde(Matrix &W_tilde) const {

          //construct matrices
          int size_pre = this->prerotation_dof_count;
          int size_post = this->postrotation_dof_count;

          //Construct S matrix
          Matrix S(size_post, size_pre, 0.0);
          if(!construct_S(S))
               return false;

          //declare, initialize and construct W_post for postrotational DOFs

          Matrix W_post(size_post, size_post, 0.0);

          double phi_psi_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_phi_psi));
          double omega_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_omega));
          double bond_angle_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_bond_angle));
          
          for (DofIterator<ChainFB> it=breakpoint; it!=end; ++it) {
               if (it.get_dof_type() == definitions::DIHEDRAL) {
                    if(it.get_atom()->atom_type == definitions::N) {
                         W_post(it.offset, it.offset) += omega_constraint;
                    } else {
                         W_post(it.offset, it.offset) += phi_psi_constraint;
                    }
               } else {
                    W_post(it.offset, it.offset) += bond_angle_constraint;
               }
          }
          
          //declaration and initialization of M_tilde
          Matrix M_tilde(size_pre, size_pre, 0.0);

          // Construct S(t)W_postS and S(t)W_post, we will need them to build M_tilde
          // and to find the new means mu_tilde respectively
          Matrix SW_post = S.transpose()*W_post;
          Matrix SW_postS(size_pre, size_pre);
          SW_postS = SW_post*S;
          
          Matrix SW(size_pre, size_pre+size_post, 0.0);
          M_tilde = W_tilde+SW_postS;
          for (int i = 0; i<size_pre; i++){
               for (int j = 0; j<size_pre; j++) {
                    SW(i,j) = W_tilde(i,j);
               }
          }
          
          for (int i = 0; i<size_pre; i++) {
               for (int j = size_pre; j<size_pre+size_post; j++){
                    SW(i,j) = SW_post(i,j-size_pre);
               }
          }
          
          W_tilde = M_tilde;

          return true;
     }


     //! construct S matrix, as the product of the two special matrices Z and Lambda.
     //! \param S Matrix containing the first order approximation to chain closure 
     //! \return Boolean indicating whether the construction of S was successful
     bool construct_S(Matrix &S) const {

          //construct Z
          MatrixVec Z(this->postrotation_dof_count, this->postrotation_dof_count-2);
          if(!construct_Z(Z))
               return false;

          // Construct Lambda
          MatrixVec Lambda(this->postrotation_dof_count-2, this->prerotation_dof_count);
          construct_lambda(Lambda);

          S = Z*Lambda;
          return true;
     }


     //! Construct Z matrix
     //! \param Z Matrix containing the derivatives of the  postrotated DOFs chi with respect to the bond vectors p
     //! \return Boolean indicating whether the construction of Z was successful
     bool construct_Z(MatrixVec &Z) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          int dchi_size = this->postrotation_dof_count-2;                    

          // vector of indices for where in the Z matrix
          // to place the dchi value
          // int indices[this->postrotation_dof_count] = {1,1,2,2,3,4};
          int indices[6] = {1,1,2,2,3,4};

          // Iterate over all post-rotational degrees of freedom
          for (DofIterator<ChainFB> it=this->breakpoint; (it!=this->end &&
                                                 it.offset < this->postrotation_dof_count); ++it) {

               //std::cout << "postrot: " << it << "\n";
               // Initialize with null vectors
               std::vector<Vector_3D> dchi(dchi_size);
               for (int i=0; i<dchi_size; i++)
                    dchi[i] = null_vector();
          
               if (it.get_dof_type() == definitions::DIHEDRAL) {

                    // Get relevant atoms
                    Atom *prevPrev = it.get_atom()->get_neighbour(-2, BACKBONE);
                    Atom *prev = it.get_atom()->get_neighbour(-1, BACKBONE);
                    Atom *current = it.get_atom();
                    Atom *next = it.get_atom()->get_neighbour(+1, BACKBONE);

                    // Get relevant bond vectors
                    Vector_3D pa = prev->position - prevPrev->position;
                    Vector_3D pb = current->position - prev->position;
                    Vector_3D pc = next->position - current->position;

                    // Calculate dchi values
                    get_dchi_dp_dihedral(pa, pb, pc, dchi, *it, indices[it.offset]);

               } else {

                    // Get relevant atoms
                    Atom *prev = it.get_atom()->get_neighbour(-1, BACKBONE);
                    Atom *current = it.get_atom();
                    Atom *next = it.get_atom()->get_neighbour(+1, BACKBONE);

                    // Get relevant bond vectors
                    Vector_3D pa = current->position - prev->position;
                    Vector_3D pb = next->position - current->position;

                    // // Check wether pa and pb are parallel
                    // double epsilon = 1.0e-03;
                    // double current_angle = *(current->angle);
                    // if((fabs(current_angle-M_PI) < epsilon) || (fabs(current_angle) < epsilon)){
                       
                    //      return false;
                    // }
                    
                    // Calculate dchi values
                    get_dchi_dp_bond_angle(pa, pb, dchi, *it, indices[it.offset]);
               }

               // Place dchi values in Z matrix
               for (int i=0; i<dchi_size; i++) {
                    Z(it.offset, i) = dchi[i];
               }
          }
          return true;
     }
     

     //! Call the relevant methods for construction of an end move
     //! \param Lt Choleski decomposition of W tilde
     //! \param W_tilde Matrix containing structural priors
     void construct_endmove(Matrix &Lt, Matrix &W_tilde) const {

          int size = prerotation_dof_count;
          W_tilde = Matrix(size, size, 0.0);

          double phi_psi_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_phi_psi));
          double omega_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_omega));
          double bond_angle_constraint = Math<double>::sqr(180.0/M_PI)*(1.0/Math<double>::sqr(settings.std_dev_bond_angle));
                                                         
          for (DofIterator<ChainFB> it=begin; it!=begin_inactive; ++it) {
               if (it.get_dof_type() == definitions::DIHEDRAL) {
                    if(it.get_atom()->atom_type == definitions::N) {
                         W_tilde(it.offset, it.offset) += omega_constraint;
                    } else {
                         W_tilde(it.offset, it.offset) += phi_psi_constraint;
                    }
               } else {
                    W_tilde(it.offset, it.offset) += bond_angle_constraint;
               }
          }

          // Calculate "square root" of matrix W_tilde
          Lt = W_tilde.cholesky('U');
     }

     //! Calculate the derivatives of the postrotated DOFs chi with
     //! respect to the bond vectors p for the bond angles. 
     //! \param p_i Bond vector at position i
     //! \param p_i_plus Bond vector at position i+1
     //! \param dchi Dchi/dp vector
     //! \param theta Bond angle value
     //! \param index index in Z matrix
     void get_dchi_dp_bond_angle(Vector_3D &p_i,Vector_3D &p_i_plus,
                                std::vector<Vector_3D> &dchi,
                                double theta, int index) const {

          int dchi_size = this->postrotation_dof_count-2;

          // Where to place dchi values in the dchi array
          int index_chi = index;
          int index_chi_plus = index+1;
          
          if (index_chi < dchi_size) {
               double prefactor = (+1.0)/(sin(theta)*p_i.norm()*p_i_plus.norm());
               
               dchi[index_chi] = prefactor*p_i_plus;

               if (index_chi_plus < dchi_size) {
                    dchi[index_chi_plus] = prefactor*p_i;
               }
          }
     }
     
     //! Calculate the derivatives of the postrotated DOFs dchi with
     //! respect to the bond vectors p for the dihedral angles. 
     //! \param p_i_minus Bond vector at position i -1 
     //! \param p_i Bond vector at position i
     //! \param p_i_plus Bond vector at position i+1
     //! \param dchi Dchi/dp vector
     //! \param chi dihedral angle value
     //! \param index index in Z matrix
     void get_dchi_dp_dihedral(Vector_3D &p_i_minus,Vector_3D &p_i, Vector_3D &p_i_plus,
                               std::vector<Vector_3D> &dchi,
                               double chi, int index) const {

          int dchi_size = this->postrotation_dof_count-2;
          
          //define binormal vectors
          Vector_3D b_minus = p_i_minus%p_i;
          double b_minus_length = b_minus.norm();
          b_minus = b_minus.normalize();
          Vector_3D b_plus = p_i%p_i_plus;
          double b_plus_length = b_plus.norm();
          b_plus = b_plus.normalize();

          // Where to place dchi values in the dchi array
          int index_chi_minus = index-1;
          int index_chi = index;
          int index_chi_plus = index+1;

          // Two solutions: sine and cosine variants
          if (std::fabs(sin(chi)) > ((1.0)/(sqrt(2.0)))){

               // Test if index_chi_minus is within bounds of dchi array
               if (index_chi_minus < dchi_size) {
               
                    double prefactor_1 = (-1.0)/(sin(chi)*b_minus_length);

                    dchi[index_chi_minus] = (prefactor_1*((p_i%b_plus)-
                                                          ((cos(chi))*p_i%b_minus)));

                    // Test if index_chi is within bounds of dchi array
                    if (index_chi < dchi_size) {
                         double prefactor_2 = (-1.0)/(sin(chi)*b_plus_length);
                         dchi[index_chi] = (((-1.0*prefactor_1)*((p_i_minus%b_plus)-
                                                                 (cos(chi)*(p_i_minus%b_minus))))+
                                            (prefactor_2*((p_i_plus%b_minus)-(cos(chi))*
                                                          (p_i_plus%b_plus))));

                         
                         // Test if index_chi_plus is within bounds of dchi array
                         if (index_chi_plus < dchi_size) {
                              dchi[index_chi_plus] = (((-1.0)*prefactor_2)*
                                                      ((p_i%b_minus)-(cos(chi))*(p_i%b_plus)));
                         }
                    }
               }
               
          } else {

               // Test if index_chi_minus is within bounds of dchi array
               if (index_chi_minus < dchi_size) {
                    
                    dchi[index_chi_minus] = ((1.0)/(cos(chi)*b_minus_length)*
                                             ((p_i.norm()*b_plus)+(sin(chi)*(b_minus%p_i))));
                    
                    // Test if index_chi is within bounds of dchi array
                    if (index_chi < dchi_size) {

                         Vector_3D a = b_minus%b_plus;
                         Vector_3D b = (b_plus%p_i)%p_i_minus + p_i.norm()*sin(chi)*(p_i_minus%b_minus);
                         Vector_3D c = (b_minus%p_i)%p_i_plus - p_i.norm()*sin(chi)*(p_i_plus%b_plus);
                         dchi[index_chi] = ((1.0/(p_i.norm()*cos(chi)))*
                                            (a + b/b_minus_length + c/b_plus_length));

                         // Test if index_chi_plus is within bounds of dchi array
                         if (index_chi_plus < dchi_size) {

                              dchi[index_chi_plus] = ((1.0/(cos(chi)*b_plus_length))*
                                                      ((p_i.norm()*b_minus)-(sin(chi)*(b_plus%p_i))));
                         }
                    }
               }
          }
     }

     //! Construct Lambda matrix
     //! \param lambda Matrix containing dp/dphi, i.e. the derivatives of the bond vectors p with respect to the prerotated DOFs phi
     void construct_lambda(MatrixVec & lambda) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Define relevant vectors and matrices
          Residue &residue_breakpoint = *this->breakpoint.get_residue();
          Residue &prev_residue_breakpoint = (*this->breakpoint.get_residue()->get_neighbour(-1));

          Vector_3D p_0 = (residue_breakpoint[N]->position -
                           prev_residue_breakpoint[C]->position);
          Vector_3D p_1 = (residue_breakpoint[CA]->position -
                           residue_breakpoint[N]->position);
          Vector_3D C_alpha_pos = residue_breakpoint[CA]->position;

          // Construct PI_2 and PI_3
          Matrix_3D PI_2 = null_matrix();
          Matrix_3D PI_3 = null_matrix();
          construct_PI(PI_2,PI_3);

          // Fill in lambda matrix
          for (DofIterator<ChainFB> it=this->begin; it!=this->begin_inactive; ++it) {

               Vector_3D rot_vector = it.get_rotation_vector();
               lambda(0,it.offset) = rot_vector%p_0;
               if ((unsigned int)it.offset == this->prerotation_dof_count-1) {
                    lambda(0,it.offset) = null_vector();
               }
               lambda(1,it.offset) = rot_vector%p_1;
               lambda(2,it.offset) = PI_2*(rot_vector%(C_alpha_pos - it.get_position_vector()));
               lambda(3,it.offset) = PI_3*(rot_vector%(C_alpha_pos - it.get_position_vector()));
          }
          
     }

     //! Construct PI matrices
     //! \param PI_2 Transformation matrix dp_2 = PI_2*dr_1
     //! \param PI_3 Transformation matrix dp_3 = PI_3*dr_1
     virtual void construct_PI(Matrix_3D &PI_2, Matrix_3D &PI_3) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          Residue &residue_breakpoint = *this->breakpoint.get_residue();
          Residue &residue_end = *this->end.get_residue();

          // define the relevant vectors
          Vector_3D p2 = (residue_breakpoint[C]->position -
                          residue_breakpoint[CA]->position);
          Vector_3D p3 = (residue_end[N]->position -
                          residue_breakpoint[C]->position);
          double len_p3_sq = p3.norm()*p3.norm();
          Vector_3D p4 = (residue_end[CA]->position -
                          residue_end[N]->position);
          double len_p4_sq = p4.norm()*p4.norm();
          Vector_3D delta = (residue_end[N]->position -
                             residue_breakpoint[CA]->position);
          double len_delta = delta.norm();
     
          Vector_3D e1 = delta/len_delta;
          Vector_3D e2 = (p4-(p4*e1)*e1).normalize();
          Vector_3D e3 = e1%e2;

          double p_21 = p2*e1;
          double p_31 = p3*e1;
          double p_32 = p3*e2;
          double p_33 = p3*e3;
          double p_41 = p4*e1;
          double p_42 = p4*e2;
          double p_3_p_4 = p3*p4;
     
          double cos_omega = ((p2%p3)*(p3%p4))/((p2%p3).norm() * (p3%p4).norm());
          double cos_omega_sq = cos_omega*cos_omega;
          double sin_omega_sq = 1.0 - cos_omega_sq;
          double A = (p_3_p_4*(p_31*p_3_p_4*sin_omega_sq - len_p3_sq*p_41) +
                      len_p3_sq*len_p4_sq*p_31*cos_omega_sq);
          double B = len_p3_sq*(p_41*len_p3_sq - p_3_p_4*p_31);
          double s = (p_31*p_42-p_32*p_41)/p_42;
     
          double denom = (len_p3_sq*(p_3_p_4*cos_omega_sq - p_31*p_41) +
                          p_3_p_4*p_31*p_31*sin_omega_sq);
          A = -A/denom;
          double A_tilde = (p_21/p_42)*(A-p_41);
          B = -B/denom;
          double B_tilde = B - s;
          double F = 0.0;
          double G = 0.0;
          // Check p_33 to avoid numerical issues around p_33 = 0
          if(std::fabs(p_33) > 1e-05){
               F = -((p_32*A_tilde + p_21*p_31)/p_33);
               G = -((p_32*B_tilde)/p_33);
          }
          PI_3 = (e1^e1)*p_21 - (e1^e2)*p_32 - (e1^e3)*p_33
               + (e2^e1)*A_tilde  + (e2^e2)*(B+(p_32*p_41)/p_42) + (e2^e3)*((p_33*p_41)/p_42)
               + (e3^e1)*F + (e3^e2)*G + (e3^e3)*s;
          
          PI_3 = PI_3*(-1.0/len_delta);
          Matrix_3D Id = identity_matrix();
          PI_2 = (Id + PI_3)*(-1.0);
     }
     
};


//////////////////////////////////////
//        POSTROTATION OBJECTS      //
//////////////////////////////////////

//! Prerotation class
//! Contains functionality and state of last postrotation
class PostRotation {
protected:
   
     //! Chain object
     ChainFB *chain;

     //@{
     //! Degree-of-freedom iterators          
     const DofIterator<ChainFB> &breakpoint;
     const DofIterator<ChainFB> &end;
     const DofIterator<ChainFB> &backup_breakpoint;
     const DofIterator<ChainFB> &backup_end;
     //@}

public:
     //! Local settings class
     const class Settings {
          //! Make output operator a friend
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {return o;}
     } settings;    //!< Local settings object

     //! \param id String identification of Postrotation
     std::string id;

     //! Moves using this postrotation are set as LOCAL
     static const definitions::MoveTypeEnum moveType = definitions::LOCAL;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Settings object
     //! \param id String identification of postrotation
     PostRotation(ChainFB *chain, const Settings &settings=Settings(),
                  std::string id="PostRotation")
          : chain(chain),
            breakpoint((Atom *)NULL), end((Atom *)NULL),
            backup_breakpoint((Atom *)NULL), backup_end((Atom *)NULL), 
            settings(settings), id(id) {}
          
     //! Copy constructor, with specification of iterators
     //! \param chain Molecule chain
     //! \param breakpoint Breakpoint dof iterator
     //! \param end End dof iterator
     //! \param backup_breakpoint Breakpoint dof iterator in backup chain 
     //! \param backup_end End dof iterator in backup chain 
     //! \param settings Settings object
     //! \param id String identification of postrotation
     PostRotation(ChainFB *chain,
                  DofIterator<ChainFB> &breakpoint, DofIterator<ChainFB> &end,
                  DofIterator<ChainFB> &backup_breakpoint, DofIterator<ChainFB> &backup_end,
                  const Settings &settings=Settings(),
                  std::string id="PostRotation")
          : chain(chain), 
            breakpoint(breakpoint), end(end), 
            backup_breakpoint(backup_breakpoint), backup_end(backup_end),
            settings(settings), id(id) {
     }

     //! Copy constructor
     //! \param other Postrotation object from which the copy is made
     //! \param breakpoint Breakpoint dof iterator
     //! \param end End dof iterator
     //! \param backup_breakpoint Breakpoint dof iterator in backup chain 
     //! \param backup_end End dof iterator in backup chain 
     PostRotation(PostRotation &other,
                  DofIterator<ChainFB> &breakpoint, DofIterator<ChainFB> &end,
                  DofIterator<ChainFB> &backup_breakpoint, DofIterator<ChainFB> &backup_end)
          : chain(other.chain), 
            breakpoint(breakpoint), end(end), 
            backup_breakpoint(backup_breakpoint), backup_end(backup_end),
            settings(other.settings),
            id(other.id){}

     // Copy Constructor - thread specified (multi-threading)
     //! \param other Postrotation object from which the copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index of the thread
     //! \param chain Molecule chain
     PostRotation(const PostRotation &other,                                      
                  RandomNumberEngine *random_number_engine,
                  int thread_index,
                  ChainFB *chain)
          : chain(chain), 
            breakpoint(other.breakpoint), end(other.end), 
            backup_breakpoint(other.backup_breakpoint), backup_end(other.backup_end),
            settings(other.settings),
            id(other.id){}

     //! Destructor
     virtual ~PostRotation(){}
     
     //! Execute postrotation 
     //! \return Boolean indicating whether chain closure was succesful
     bool apply() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          Residue &residue_breakpoint = *breakpoint.get_residue();
          Residue &residue_end = *end.get_residue();

          Residue &backup_residue_breakpoint = *backup_breakpoint.get_residue();
          Residue &backup_residue_end = *backup_end.get_residue();
          DofIterator<ChainFB> local_begin = breakpoint;
          DofIterator<ChainFB> local_end = end;
          Vector_3D p2 = backup_residue_breakpoint[C]->position - backup_residue_breakpoint[CA]->position;
          double len_p2 = p2.norm();
          
          Vector_3D p3 = backup_residue_end[N]->position - backup_residue_breakpoint[C]->position;
          double len_p3 = p3.norm();
          Vector_3D delta = (residue_end[N]->position - residue_breakpoint[CA]->position);
          double len_delta = delta.norm();

          // Check whether solution exists
          if (len_delta > len_p2 + len_p3|| (len_delta < std::fabs(len_p2-len_p3))){
               //std::cout << "distance check fails" << "\n";
               return false;
          }
          // Define vectors and scalar products needed for finding the solution
          Vector_3D p4 = backup_residue_end[CA]->position - backup_residue_end[N]->position;
          double len_p4 = p4.norm();

          Vector_3D e1 = delta/len_delta;
          Vector_3D e2 = (p4-((p4*e1)*e1)).normalize();
          Vector_3D e3 = e1%e2;
          Vector_3D delta_old = (backup_residue_end[N]->position - backup_residue_breakpoint[CA]->position);
          double len_delta_old = delta_old.norm();
          Vector_3D e1_old = delta_old/len_delta_old;
          double p_3_1 = (len_delta*len_delta  + len_p3*len_p3 - len_p2*len_p2)/(2.0*len_delta);
         
          double cos_omega = ((p2%p3)*(p3%p4))/((p2%p3).norm() * (p3%p4).norm());
          double sin_sq_omega = 1.0 - (cos_omega*cos_omega);
          double p_3_perp_squared = len_p3*len_p3 - p_3_1*p_3_1;
          double p_4_1 = p4*e1;
          double p_4_2 = p4*e2;
          double determinant = len_p3*len_p3*p_4_2*p_4_2 - len_p4*len_p4*p_3_perp_squared*sin_sq_omega;

          // Check the sign of the Delta
          if(determinant < 0.0)          return false; 
          Vector_3D e2_old = (p4-((p4*e1_old)*e1_old)).normalize();
          Vector_3D e3_old = e1_old%e2_old;

          double a = p_3_1*p_4_1*len_p3*len_p3;
          double b = len_p3*std::fabs(cos_omega)*sqrt(p_3_perp_squared*determinant);
          double d = len_p3*len_p3 - p_3_perp_squared*sin_sq_omega;
          double p_3_4_plus =  (a+b)/d;
          double p_3_4_minus = (a-b)/d;
          double p_3_2_minus = (p_3_4_minus - p_3_1*p_4_1)/p_4_2;
          double p_3_2 = (p_3_4_plus - p_3_1*p_4_1)/p_4_2;
          if(std::fabs(p_3_2_minus-p3*e2_old) < std::fabs(p_3_2-p3*e2_old))
               p_3_2 = p_3_2_minus;
          double p_3_3_sq = len_p3*len_p3 - p_3_1*p_3_1 - p_3_2*p_3_2;
          if(p_3_3_sq < 1e-10){
               p_3_3_sq = 0.0;
          }

          double p_3_3 =  sqrt(p_3_3_sq);
          if(p3*e3_old < 0)
               p_3_3 *= -1.0;

          Vector_3D p_3_new = p_3_1*e1 + p_3_2*e2 + p_3_3*e3;
                  
          // Reposition the C atom
          Vector_3D Cpos = residue_end[N]->position - p_3_new;
          for (int i=0; i<3; ++i) {
               assert(std::isfinite(Cpos[i]));
          }
          residue_breakpoint[C]->position = Cpos;

          // Check the invariance of the bond length and of the omega angle
          Vector_3D p_3_current = residue_end[N]->position - residue_breakpoint[C]->position;
          Vector_3D p_2_current = residue_breakpoint[C]->position - residue_breakpoint[CA]->position;
          double cos_omega_new = ((p_2_current%p_3_current)*(p_3_current%p4))/((p_2_current%p_3_current).norm() * (p_3_current%p4).norm());
         
          // If the sign is  wrong, reject the solution
          if(std::fabs(cos_omega-cos_omega_new) > std::fabs(cos_omega)/1000.0){
               return false;
          }

          assert(std::fabs(len_p3 - (residue_end[N]->position - residue_breakpoint[C]->position).norm()) < len_p3/100.0);
          assert(std::fabs(len_p2 -(residue_breakpoint[C]->position - residue_breakpoint[CA]->position).norm()) < len_p2/100.0);
          
          // Update angles
          residue_breakpoint[CA]->update_dihedral();
          residue_breakpoint[CA]->update_angle();
          residue_breakpoint[C]->update_dihedral();
          residue_breakpoint[C]->update_angle();
          residue_end[N]->update_angle();
          residue_end[CA]->update_dihedral();

          // Enforce Omega assumption
          // (if it wasn't exactly planar before, it is now)
          residue_end[N]->update_dihedral();
          
          return true;
     }

     //! Calculate the log of the Jacobian introduced by postrotation
     //! \return Log jacobian 
     double get_log_bias() {
          double log_Jacobian_before = calculate_log_Jacobian_before();
          double log_Jacobian_after = calculate_log_Jacobian_after();

          return log_Jacobian_after - log_Jacobian_before;
     }
     
     //! Calculate log-Jacobian on current structure
     //! \return Log Jacobian for proposed move
     double calculate_log_Jacobian_before() {
          return calculate_log_Jacobian(backup_breakpoint, backup_end);
     }

     //! Calculate log-Jacobian after move
     //! \return Log-Jacobian for reverse move
     double calculate_log_Jacobian_after() {
          return calculate_log_Jacobian(breakpoint, end);
     }
     
     //! Calculate logarithm of Jacobian determinant
     //! \param breakpoint Breakpoint dof iterator
     //! \param end End dof iterator
     //! \return Log-Jacobian 
     double calculate_log_Jacobian(const DofIterator<ChainFB> &breakpoint, const DofIterator<ChainFB> &end) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          Matrix A(5,5);

          Vector_3D s = (*end.get_atom()->residue)[N]->position;
          Vector_3D u = (*end.get_atom()->residue)[CA]->position - s;   //u.normalize() not needed.

          // The Omega_4 dihedral and Alpha_4 angle are excluded from end
          DofIterator<ChainFB> local_begin = breakpoint;
          DofIterator<ChainFB> local_end = end-1;
          for (DofIterator<ChainFB> it=local_begin; it!=local_end; ++it) {
               Vector_3D rotation = it.get_rotation_vector();
               Vector_3D ds_dphi_i = rotation % (s - it.get_position_vector());
               Vector_3D du_dphi_i = rotation % u;

               int i = it.offset;
               for (int j=0; j<3; j++) {
                    A(j,i) = ds_dphi_i[j];
               }
               for (int j=3; j<5; j++) {
                    A(j,i) = du_dphi_i[j-3];
               }	       
          }

          // Check for possible linearity between rows in A
          // If true is returned the row number "null_vector_index" needs  to be replaced 
          int null_vector_index = -1;
          bool linearity = linearity_check(A.get_row(3), A.get_row(4), &null_vector_index);
	  
          //correcting the A matrix if the linearity check returns true
          if (linearity) {
               for (DofIterator<ChainFB> it=local_begin; it!=local_end; ++it) {
                    Vector_3D rotation = it.get_rotation_vector();
                    Vector_3D du_dphi_i = rotation % u;
                    int i = it.offset;
                    A(null_vector_index+3,i) = du_dphi_i[2];
               }
          }
			  
          return -std::log(std::fabs(A.determinant())); 
     }

     //! Check whether du_dphi vectors are independent or not.
     //! \param du_dphi_e1 Vector du_dphi dot e1 in Jacobian matrix
     //! \param du_dphi_e2 Vector du_dphi dot e2 in Jacobian matrix
     //! \param vector_index vector index 
     //! \return Boolean indicating whether du_dphi_e1 and du_dphi_e2 are linearly dependent
     bool linearity_check(const Vector_nD &du_dphi_e1, const Vector_nD &du_dphi_e2, int* vector_index) {
          double epsilon_length = 1E-6;
          double epsilon_product = 1E-3;

          if (du_dphi_e1*du_dphi_e1 < epsilon_length){
               *vector_index = 0;
               return true;
          } 
          if (du_dphi_e2*du_dphi_e2 < epsilon_length){
               *vector_index = 1;
               return true;
          }
          double scalar_product = (du_dphi_e1*du_dphi_e2)/((du_dphi_e1.norm())*(du_dphi_e2.norm()));
          if (std::fabs(scalar_product-1.0)<epsilon_product) {
               *vector_index = 1;
               return true;
          }
          return false;
     }

     //! Accept last move
     void accept() {
     }

     //! Reject last move
     void reject() {
     }     
};
    
}


//! CRISP move class
template <typename CHAIN_TYPE,
          typename PREROTATION_TYPE=crisp::PreRotation<>,
          typename POSTROTATION_TYPE=crisp::PostRotation>
class MoveCRISP: public MoveCommon<MoveCRISP<CHAIN_TYPE,
                                               PREROTATION_TYPE,
                                               POSTROTATION_TYPE>,
                                    CHAIN_TYPE> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveCRISP<CHAIN_TYPE,
                                               PREROTATION_TYPE,
                                               POSTROTATION_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Pre-rotation template
     PREROTATION_TYPE preRotationFactory;

     //! Post-rotation template
     POSTROTATION_TYPE postRotationFactory;

     //! Current pre rotation object
     PREROTATION_TYPE *preRotation;

     //! Current pre rotation object
     POSTROTATION_TYPE *postRotation;

     //@{
     //! Degree-of-freedom iterators          
     DofIterator<ChainFB> *begin;
     DofIterator<ChainFB> *breakpoint;
     DofIterator<ChainFB> *end;
     DofIterator<ChainFB> *backup_begin;
     DofIterator<ChainFB> *backup_breakpoint;
     DofIterator<ChainFB> *backup_end;     
     //@}

     //! Local enum: move regions
     definitions::TerminalEnum move_region;

     //! Make moveStatistics a friend
     friend class MoveStatisticsLocalMove<CHAIN_TYPE, MoveCRISP>;
     
public:

     //! Move settings
     const class Settings: public Move<CHAIN_TYPE>::Settings, 
                           public PREROTATION_TYPE::Settings,
                           public POSTROTATION_TYPE::Settings {
     public:

          //! Whether to perform only internal moves
          bool only_internal_moves;
          
          //! Specifies the range of bond_angles that are excluded from the phase-space
          //! in order to avoid gimbal lock around 0 and pi
          //! The range is [cutoff, pi-cutoff]
          double bond_angle_cutoff;

          //! Whether omega is resampled during prerotation
          bool sample_omega;

          //! Whether bond angles are resampled during prerotation
          bool sample_bond_angle;
          
          //! Constructor
          Settings(bool only_internal_moves=false,
                   double bond_angle_cutoff = 0.1,
                   double sample_omega = false,
                   double sample_bond_angle = true)
               : Move<CHAIN_TYPE>::Settings(5,5),
                 only_internal_moves(only_internal_moves),
                 bond_angle_cutoff(bond_angle_cutoff),
                 sample_omega(sample_omega),
                 sample_bond_angle(sample_bond_angle){
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "only-internal-moves:" << settings.only_internal_moves << "\n";
               o << "bond-angle-cutoff:" << settings.bond_angle_cutoff << "\n";
               o << "sample-omega:" << settings.sample_omega << "\n";
               o << "sample-bond-angle:" << settings.sample_bond_angle << "\n";
               o << static_cast<const typename PREROTATION_TYPE::Settings &>(settings);
               o << static_cast<const typename POSTROTATION_TYPE::Settings &>(settings);
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                    
     } settings;


     //! Initializer
     void init() {

          // Set name
          this->id = "crisp<"+preRotationFactory.id+","+postRotationFactory.id+">";

          // // The postrotation determines whether the move is local or not
          // this->moveType = postRotationFactory.moveType;

          // Statistics object
          delete this->statistics;  // Remove existing statistics object
          this->statistics = new MoveStatisticsLocalMove<CHAIN_TYPE, MoveCRISP>(this);          
     }

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveCRISP(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "crisp", settings, random_number_engine),
            preRotationFactory(chain, settings, random_number_engine), 
            postRotationFactory(chain, settings),
            preRotation(NULL), postRotation(NULL),
            settings(settings) {

          init();
     }

     //! Constructor for DBN-based preRotation
     //! \param chain Molecule chain
     //! \param dbn DBN object
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     template <typename DBN_TYPE>
     MoveCRISP(CHAIN_TYPE *chain, DBN_TYPE *dbn, const Settings &settings = Settings(),
               RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "crisp", settings, random_number_engine),
            preRotationFactory(chain, dbn, settings, random_number_engine), 
            postRotationFactory(chain, settings),
            preRotation(NULL), postRotation(NULL),
            settings(settings) {

          // When using a CRISP move with DBN prior, the angle node is
          // flagged as fixed and is then unfixed every time a move is
          // made. This ensures that when the dbn is resynchronized
          // after a change has been made to the angles outside the
          // DBN, the hidden nodes will be resampled based on the new
          // angle values
          dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->fixed = true;

          init();
     }

     //! Copy constructor
     //! \param other MoveCRISP object from which the copy is made
     MoveCRISP(const MoveCRISP<CHAIN_TYPE,PREROTATION_TYPE,POSTROTATION_TYPE> &other)
          : MoveCommon(other),
            preRotationFactory(other.preRotationFactory),
            postRotationFactory(other.postRotationFactory),
            preRotation(NULL), postRotation(NULL),
            settings(other.settings)  {
          init();
     }

     
     //! Copy constructor - thread specified (multi-threading)
     //! \param other MoveCRISP object from which the copy is made     
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index of the thread
     //! \param chain Molecule chain
     MoveCRISP(const MoveCRISP<CHAIN_TYPE,PREROTATION_TYPE,POSTROTATION_TYPE> &other, 
               RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            preRotationFactory(other.preRotationFactory, random_number_engine, thread_index, chain),
            postRotationFactory(other.postRotationFactory, random_number_engine, thread_index, chain),
            preRotation(NULL), postRotation(NULL),
            settings(other.settings) {
          init();
     }


     //! Apply local move
     //! \param end_index End index in sequence
     //! \param start_index Start index in sequence
     //! \return Boolean indicating whether the move is successful
     // start_index<0: Resample beginning of chain without loop closure
     // end_index>chain.size: Resample end of chain without loop closure
     // start_index=end_index: Calculate random range
     bool apply(int start_index=0, int end_index=0) {

          // Import crisp namespace
          using namespace crisp;
          
          // Import protein definitions (such as residue names)
          using namespace definitions;

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

          // std::cout << start_index << " " << end_index << "\n";

          // Determine move type
          move_region = INTERNAL;
          if (start_index<0) {
               move_region = NTERM;
               start_index = 0;
          } else if (end_index > this->chain->size()-1) {
               move_region = CTERM;
               end_index = this->chain->size();               
          }

          // Set which positions and angles are updated
          // The postrotation determines whether the move is local or not
          this->move_info = preRotationFactory.create_move_info(start_index, end_index,
                                                                postRotationFactory.moveType);

          // Create temporary chain (backup in case of reject)
          if (this->chain_backup)
               delete this->chain_backup;
          this->chain_backup = new CHAIN_TYPE(*this->chain,
                                              start_index, end_index);

          // Setup the DofIterators
          int breakpoint_index = end_index-2;
          DofIterator<ChainFB>::AngleSelectionEnum begin_dofs             = DofIterator<ChainFB>::STANDARD_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum breakpoint_dofs        = DofIterator<ChainFB>::STANDARD_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum end_dofs               = DofIterator<ChainFB>::STANDARD_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum backup_begin_dofs      = DofIterator<ChainFB>::STANDARD_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum backup_breakpoint_dofs = DofIterator<ChainFB>::STANDARD_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum backup_end_dofs        = DofIterator<ChainFB>::STANDARD_DOFS;
          

          AtomEnum breakpoint_atom = CA;
          AtomEnum end_atom = CA;
          AngleEnum breakpoint_dof_type = DIHEDRAL;
          AngleEnum end_dof_type = ANGLE;

          if(settings.sample_omega) {
               begin_dofs  += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
               backup_begin_dofs  += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
          }

          if(!settings.sample_bond_angle) {
               begin_dofs -= DofIterator<ChainFB>::BONDANGLE_DOFS;
               backup_begin_dofs -= DofIterator<ChainFB>::BONDANGLE_DOFS;
               end_dof_type = DIHEDRAL;
          }


          switch(move_region) {
          case INTERNAL:
               if (start_index != 0)
                    backup_begin_dofs += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
               break;
          case NTERM:
               begin_dofs        += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
               backup_begin_dofs += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
               breakpoint_index = end_index-1;
               breakpoint_dof_type = ANGLE;
               
               if(!settings.sample_bond_angle) {
                    breakpoint_atom = C;
                    breakpoint_dof_type = DIHEDRAL;
               }
               
               break;
          case CTERM:
               begin_dofs        += DofIterator<ChainFB>::CTERM_C_DIHEDRAL + DofIterator<ChainFB>::CTERM_C_ANGLE;
               backup_begin_dofs +=(DofIterator<ChainFB>::CTERM_C_DIHEDRAL + DofIterator<ChainFB>::CTERM_C_ANGLE +
                                    DofIterator<ChainFB>::NTERM_CA_DIHEDRAL);
               breakpoint_dofs        += DofIterator<ChainFB>::CTERM_C_DIHEDRAL + DofIterator<ChainFB>::CTERM_C_ANGLE;
               backup_breakpoint_dofs += DofIterator<ChainFB>::CTERM_C_DIHEDRAL + DofIterator<ChainFB>::CTERM_C_ANGLE;
               end_dofs               += DofIterator<ChainFB>::CTERM_C_DIHEDRAL + DofIterator<ChainFB>::CTERM_C_ANGLE;
               backup_end_dofs        += DofIterator<ChainFB>::CTERM_C_DIHEDRAL + DofIterator<ChainFB>::CTERM_C_ANGLE;
               breakpoint_atom = C;
               end_atom = C;
               breakpoint_dof_type = ANGLE;
               if(!settings.sample_bond_angle) {
                    breakpoint_dof_type = DIHEDRAL;
                    begin_dofs        -= DofIterator<ChainFB>::CTERM_C_ANGLE;
                    backup_begin_dofs -= DofIterator<ChainFB>::CTERM_C_ANGLE;
                                        
                    breakpoint_dofs        -=  DofIterator<ChainFB>::CTERM_C_ANGLE;
                    backup_breakpoint_dofs -=  DofIterator<ChainFB>::CTERM_C_ANGLE;
                    end_dofs               -=  DofIterator<ChainFB>::CTERM_C_ANGLE;
                    backup_end_dofs        -=  DofIterator<ChainFB>::CTERM_C_ANGLE;
               }

               breakpoint_index = end_index-1;
               break;
          }
          
          

          // Construct Degree-of-freedom iterators
          begin        = new DofIterator<ChainFB>((*this->chain)(start_index, CA),
                                         DIHEDRAL,
                                         begin_dofs);
          backup_begin = new DofIterator<ChainFB>((*this->chain_backup)(0, CA),
                                         DIHEDRAL,
                                         backup_begin_dofs);

          breakpoint        = new DofIterator<ChainFB>((*this->chain)(breakpoint_index, breakpoint_atom),
                                              breakpoint_dof_type,
                                              breakpoint_dofs);
          backup_breakpoint = new DofIterator<ChainFB>((*this->chain_backup)(breakpoint_index-start_index,
                                                                    breakpoint_atom),
                                              breakpoint_dof_type,
                                              backup_breakpoint_dofs);
               
          end        = new DofIterator<ChainFB>((*this->chain)(end_index-1, end_atom),
                                       end_dof_type,
                                       end_dofs);
          backup_end = new DofIterator<ChainFB>((*this->chain_backup)(end_index-1-start_index, end_atom),
                                       end_dof_type,
                                       backup_end_dofs);
          

          // Initialize prerotation
          preRotation = new PREROTATION_TYPE(preRotationFactory, move_region,
                                             end_index-start_index,
                                             *begin, *breakpoint, *end,
                                             *backup_begin, *backup_breakpoint, *backup_end);

          // Initialize postrotation
          if(move_region==INTERNAL) {
               postRotation = new POSTROTATION_TYPE(postRotationFactory, 
                                                    *breakpoint, *end, *backup_breakpoint, *backup_end);
          }

          // Apply prerotation
          if(!preRotation->apply()){
               return (this->move_info->success = false);
          }
          // Optionally apply postrotation
          if (move_region == INTERNAL)
               this->move_info->success = postRotation->apply();



          if (this->move_info->success){

               // Ensure that no bond angles are in illegal range
               for (DofIterator<ChainFB> it=*begin; it!=*end; ++it) {
                    
                    if (it.get_dof_type() == ANGLE) {
                         // Check whether angle is close to 0 or pi
                         if((std::fabs(*it-M_PI) < settings.bond_angle_cutoff) || (std::fabs(*it) < settings.bond_angle_cutoff)){
                              return (this->move_info->success = false);
                         }
                    }
               }

               // Calculate new positions non backbone atoms
               for (int i=start_index; i<end_index; i++) {
                    (*this->chain)[i].update_positions_non_backbone();
               }

               // Update state of preRotation object (e.g. DBN object)
               preRotation->update_state(this->chain);
          } 

          return this->move_info->success;
     }
     
     
     //! Calculate the log-bias introduced by the move
     //! \return Log bias
     double get_log_bias() {

          // Import crisp namespace
          using namespace crisp;
          
          double bias = 0.0;

          // Always evaluate bias of preRotation
          bias += preRotation->get_log_bias();

          // Optionally evaluate bias of postrotation
          if (move_region == definitions::INTERNAL) {
               bias += postRotation->get_log_bias();
          }

          return bias;
     }

     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();

          // Call accept in preRotation (update state in prior)
          preRotation->accept();
          
          // Call accept in postRotation
          if (move_region == definitions::INTERNAL) {
               postRotation->accept();
          }
          
          // Cleanup backup chain
          delete this->chain_backup;
          this->chain_backup = NULL;

          // Cleanup pre and post rotation objects
          delete preRotation;
          preRotation = NULL;
          if (postRotation) {
               delete postRotation;
               postRotation = NULL;
          }

          // Cleanup dof-iterators
          delete begin;
          delete breakpoint;
          delete end;
          delete backup_begin;
          delete backup_breakpoint;
          delete backup_end;
     }     

     //! Reject last move
     void reject() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Call base class reject method
          Move<CHAIN_TYPE>::reject();

          // Update positions and angles in the current chain based on the backup chain
          assert(this->chain_backup != NULL);

          AtomIterator<ChainFB,ALL> it_current(*this->chain, 
                                               this->move_info->modified_positions[0].first,
                                               this->move_info->modified_positions[0].second);
          AtomIterator<ChainFB,ALL> it_backup(*this->chain_backup);
          for (; !it_current.end() && !it_backup.end(); ++it_current,++it_backup) {
	       it_current->set_angle(it_backup->get_angle());
	       it_current->set_dihedral(it_backup->get_dihedral());
	       it_current->position = it_backup->position;
	  }

          // Cleanup backup chain
          delete this->chain_backup;
          this->chain_backup = NULL;

          // Call reject in preRotation (update state in prior)
          preRotation->reject();
                    
          // Call reject in postRotation
          if (move_region == INTERNAL) {
               postRotation->reject();
          }

          // Cleanup pre and post rotation objects
          delete preRotation;
          preRotation = NULL;
          if (postRotation) {
               delete postRotation;
               postRotation = NULL;
          }

          // Cleanup dof-iterators
          delete begin;
          delete breakpoint;
          delete end;
          delete backup_begin;
          delete backup_breakpoint;
          delete backup_end;

     }          


     //! Specifies how move reacts to the knowledge that another moves has been executed
     void notify(MoveInfo *move_info) {
          preRotationFactory.notify(move_info);
     }
     
};

}

#endif
