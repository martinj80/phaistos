// move_bgs.h --- Biased Gaussian Step move
//                Reference: Monte Carlo update for chain molecules: biased Gaussian steps in torsional space
//                           G. Favrin, A. Irbäck and F. Sjunnesson, Journal of Chemical Physics, 2001
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

#ifndef MOVE_BGS
#define MOVE_BGS


namespace phaistos {

// Name space for BGS internals
namespace bgs {

//! Gaussian Step Class
template <typename PRIOR_TYPE=MovePrior<> >
class GaussianStep{

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

          //! Value of parameter a (global constraint)
          double constraint_a; 

          //! Value of parameter b (locality constraint)
          double constraint_b; 

          //! Scaling factor for omega angles
          double omega_scaling;       

          //! Scaling factor for bond angles
          double bond_angle_scaling;   

          //! Constructor
          Settings(double constraint_a=300,
                   double constraint_b=10,
                   double omega_scaling = 8,
                   double bond_angle_scaling = 8):
               constraint_a(constraint_a),
               constraint_b(constraint_b),
               omega_scaling(omega_scaling),
               bond_angle_scaling(bond_angle_scaling){ 
          }
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "constraint-a:" << settings.constraint_a << "\n";
               o << "constraint-b:" << settings.constraint_b << "\n";
               o << "omega-scaling:" << settings.omega_scaling << "\n";
               o << "bond-angle-scaling:" << settings.bond_angle_scaling << "\n";
               o << static_cast<const typename PRIOR_TYPE::Settings &>(settings);
               return o;
          }                    

     } settings;    //!< Local settings object

     //! Prior type
     PRIOR_TYPE prior;

     //! String identification of Prior
     std::string id;
     
     //! Move region (INTERNAL,CTERM,NTERM)
     definitions::TerminalEnum move_region;
     
     //! Number of degrees of freedom effected by the bgs move
     unsigned int dof_count;
     
     //! Chain object
     ChainFB *chain;

     //! length of move (number of residues)
     int move_length;

     //! Vector of the degrees of freedom to be modified
     const std::vector<Dof> &dofs;

     //! Vector of the degrees of freedom in the backup chain
     const std::vector<Dof> &dofs_backup;

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
     GaussianStep(ChainFB *chain,
                  const Settings &settings,
                  RandomNumberEngine *random_number_engine=&random_global,
                  std::string id="GaussianStep")
          : random_number_engine(random_number_engine),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings),
            prior(settings),
            id(id+"<"+prior.id+">"),
            chain(chain),
            dofs(std::vector<Dof>()),
            dofs_backup(std::vector<Dof>()) {
     }
     
     
     //! Constructor for DBN-based Gaussian Step
     //! \param chain Molecule chain
     //! \param dbn DBN object
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param id String identification of Prior
     template <typename DBN_TYPE>
     GaussianStep(ChainFB *chain, DBN_TYPE *dbn, const Settings &settings,
                  RandomNumberEngine *random_number_engine=&random_global,
                  std::string id="bgs_step")
          : random_number_engine(random_number_engine),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings),
            prior(dbn, settings),
            id(id+"<"+prior.id+">"),
            chain(chain),
            dofs(std::vector<Dof>()),
            dofs_backup(std::vector<Dof>()) {
     }

     
     //! Copy constructor, with specification of iterators
     //! \param other Gaussian Step object from which the copy is made
     //! \param move_region Whether we are at the N or the C terminal
     //! \param move_length Length of move (number of residues)
     //! \param dofs Vector of degree-of-freedom values
     //! \param dofs_backup Vector of degree-of-freedom backup values
     GaussianStep(const GaussianStep &other, definitions::TerminalEnum move_region,
                  int move_length,
                  const std::vector<Dof> &dofs,
                  const std::vector<Dof> &dofs_backup)
          : random_number_engine(other.random_number_engine),
            random_generator_uniform_01(other.random_generator_uniform_01),
            random_generator_gaussian(other.random_generator_gaussian),
            settings(other.settings), 
            prior(other.prior, dofs.front().atom->residue->index, dofs.back().atom->residue->index+1),
            id(other.id), move_region(move_region), chain(other.chain), move_length(move_length),
            dofs(dofs),
            dofs_backup(dofs_backup) {
          
          dof_count = dofs.size();
     }

     //! Copy constructor - thread specified (multi-threading)
     //! \param other Gaussian Step object from which the copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index of the thread
     //! \param chain Molecule chain
     GaussianStep(const  GaussianStep &other, 
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
            dofs(std::vector<Dof>()),
            dofs_backup(std::vector<Dof>()) {
     }


     //! Destructor
     virtual ~GaussianStep() {}

     //! Execute  GaussianStep move
     //! \return Boolean indicating whether move was successful
     bool apply() {
          if(!apply_internalmove())
               return false;

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
          for (unsigned int i=0; i<dofs.size(); ++i) {
               assert(dofs[i].atom != NULL);
               double &value = *dofs[i].value;
               value += dphi(i,0);
               value = fmod(value+3*M_PI,(2*M_PI))-M_PI;
               
               if (dofs[i].dof_type == definitions::ANGLE) {
                    value = std::fabs(value);
               }
          }

          // Update chain positions
          chain->update_positions_backbone_segment(dofs[0].atom->residue->index,
                                                   this->chain->size());
          return true;
     }

     //! Create and return move info object 
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \param move_type type of move
     MoveInfo *create_move_info(int start_index, int end_index, definitions::MoveTypeEnum move_type) {
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

          int size = dof_count;

          dchi = Matrix(size, 1);
          for (int i=0; i<size; i++) {
               dchi(i,0) = random_generator_gaussian();
          }
     }

     
     //! Calculate the bias (implicit energy ratio) that was introduced by the  GaussianStep          
     //! \return Log probability
     double get_log_bias() const {

          // Basic move Bias
          double bias = 0.0;

          double log_gaussianstep_bias_forward = calculate_bias_forward();
          double log_gaussianstep_bias_reverse = calculate_bias_reverse();
          bias += log_gaussianstep_bias_reverse - log_gaussianstep_bias_forward;

          // Bias from prior distribution
          bias += prior.get_log_bias(dofs, dofs_backup);

          return bias;
     }

     
     //! Calculate bias for gaussianstep move (a->b)
     //! \return Log bias introduced by the proposed move
     double calculate_bias_forward() const {
          int size = dof_count;

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

          int size = dof_count;

          // Calculate Lt'
          Matrix Lt_prime;
          Matrix W_tilde_prime;
          if(!construct_all(Lt_prime, W_tilde_prime))
               return std::numeric_limits<double>::infinity();
                        
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
          Residue &residue_end = *dofs.back().atom->residue;

          // Always update positions in the forward direction (due to assymetry in CRA moves)
          int fixed_direction = 1;
          chain->update_positions(residue_end.index,
                                  chain->size(),
                                  fixed_direction);

          // Call reject in prior
          prior.reject();
     }
     
     
     //! Construct all the matrices for internal move
     //! \param Lt Choleski decomposition of W tilde
     //! \param W_tilde Matrix containing structural priors
     //! \return Boolean indicating whether the construction of the matrices and the decompoistion were successful 
     bool construct_all(Matrix &Lt, Matrix &W_tilde) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          unsigned int size = dof_count;

          // Construct covariance matrix (identity)
          Matrix W(size, size, 0.0);
          for (unsigned int i=0; i<size; ++i) {
               W(i,i) = 1.0;
          }

          // Construct I matrix
          std::vector<Matrix> I;
          std::vector<std::vector<Vector_3D> > da_dphi;
          construct_I(I, da_dphi);

          // Add d-constraint
          for (unsigned int i=0; i<I.size(); i++) {
               W += I[i]*settings.constraint_b;
          }

          W = W*(settings.constraint_a);

          W_tilde = W;
          // Calculate "square root" of matrix W_tilde
          try {
               Lt = W.cholesky('U');
          } catch (Matrix::MatrixError &e) {
               std::cout<< "cholesky decomposition failed"<< "\n";
               return false;
          }

          // Constrain bond and omega angles
          for (unsigned int i=0; i<size; ++i) {
               for (unsigned int j=i; j<size; ++j) {
                    if (dofs[j].dof_type == definitions::DIHEDRAL) {
                         if (dofs[j].atom->atom_type == N) {
                              Lt(i,j) *= settings.omega_scaling;
                         }
                    } else {
                         Lt(i,j) *= settings.bond_angle_scaling;
                    }
               }
          }

          W_tilde = (Lt.transpose())*Lt;

          return true;
     }

     //! Construct Constraint matrices
     //! \param I Constraint matrix
     //! \param da_dphi Vector of derivatives d_endpoint/d_angles
     void construct_I(std::vector<Matrix> &I, std::vector<std::vector<Vector_3D> > &da_dphi) const {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Constrain locality for CA, C and O atoms
          std::vector<Vector_3D> pos;
          pos.push_back((*chain)[dofs.back().atom->residue->index][CA]->position);
          pos.push_back((*chain)[dofs.back().atom->residue->index][C]->position);
          pos.push_back((*chain)[dofs.back().atom->residue->index][O]->position);

          int size =  dof_count;

          for (unsigned int index=0; index<pos.size(); index++) {

               I.push_back(Matrix(size, size));
               da_dphi.push_back(std::vector<Vector_3D>(size));

               // Calculate da_dphi
               for (unsigned int i=0; i<dofs.size(); ++i) {
                    Vector_3D rotation_vector;
                    definitions::AngleEnum dof_type = dofs[i].dof_type;
                    Atom *atom_current = dofs[i].atom;
                    Atom *atom_prev = atom_current->get_neighbour(-1, definitions::BACKBONE);
                    Atom *atom_next = atom_current->get_neighbour(+1, definitions::BACKBONE);
                    if (dof_type == definitions::DIHEDRAL) {
                         rotation_vector = (atom_current->position - atom_prev->position).normalize();
                    } else {
                         Vector_3D v1 = atom_current->position - atom_prev->position;
                         Vector_3D v2 = atom_next->position - atom_current->position;
                         rotation_vector = (v2%v1).normalize();
                    }

                    // The differential of a position vector around a
                    // rotation vector can be written as a cross product
                    // between the two vectors
                    da_dphi[index][i] = rotation_vector % (pos[index]-atom_current->position);
               }

               // I = da_dphi^2
               for (int i=0; i<size; i++) {
                    for (int j=0; j<size; j++) {
                         I[index](i,j)=da_dphi[index][i]*da_dphi[index][j];
                    }
               }
          }
     }
     
};

}


//! BGS move class
template <typename CHAIN_TYPE,
          typename BGS_TYPE=bgs::GaussianStep<> >
class MoveBGS: public MoveCommon<MoveBGS<CHAIN_TYPE,BGS_TYPE>,
                                  CHAIN_TYPE> {
private:
     
     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveBGS<CHAIN_TYPE,
                                    BGS_TYPE >,CHAIN_TYPE> MoveCommon;
     //! BGS template
     BGS_TYPE BGSFactory;

     //! Current gaussian steps object
     BGS_TYPE *BGSRotation;

     //! Degree-of-freedom vector
     std::vector<Dof> dofs;

     //! Degree-of-freedom vector - backup
     std::vector<Dof> dofs_backup;

     //! Local enum: move regions
     definitions::TerminalEnum move_region;

     //! Make moveStatistics a friend
     friend class MoveStatisticsLocalMove<CHAIN_TYPE, MoveBGS>;
     
public:

     //! Move settings
     const class Settings: public Move<CHAIN_TYPE>::Settings, 
                           public BGS_TYPE::Settings {
     public:

          //! Whether to perform only internal moves
          bool only_internal_moves;
       
          //! Specifies the range of bondangles that are excluded from the phase-space
          //! in order to avoid gimbal lock around 0 and pi
          //! The range is [cutoff, pi-cutoff]
          double bond_angle_cutoff;

          //! Whether omega is resampled during prerotation
          bool sample_omega;

          //! Whether bond angles are resampled during prerotation
          bool sample_bond_angle;

          //! Whether to skip prolines phi angles (modifiying the proline phi angle introduces an improper torsion change)          
          bool skip_proline_phi;

          //! Constructor
          Settings(bool only_internal_moves=false,
                   double bond_angle_cutoff = 0.1,
                   double sample_omega = false,
                   double sample_bond_angle = false,
                   bool skip_proline_phi = false)
               : Move<CHAIN_TYPE>::Settings(4,4),
                 only_internal_moves(only_internal_moves),
                 bond_angle_cutoff(bond_angle_cutoff),
                 sample_omega(sample_omega),
                 sample_bond_angle(sample_bond_angle),
                 skip_proline_phi(skip_proline_phi) {
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "only-internal-moves:" << settings.only_internal_moves << "\n";
               o << "bond-angle-cutoff:" << settings.bond_angle_cutoff << "\n";
               o << "sample-omega:" << settings.sample_omega << "\n";
               o << "sample-bond-angle:" << settings.sample_bond_angle << "\n";
               o << "skip-proline-phi:" << settings.skip_proline_phi << "\n";
               o << static_cast<const typename BGS_TYPE::Settings &>(settings);
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                    
     } settings;


     //! Initializer
     void init() {

          // Set name
          this->id = "semilocal<"+BGSFactory.id+">";
     }

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveBGS(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "semilocal_bgs", settings, random_number_engine),
            BGSFactory(chain, settings, random_number_engine), 
            BGSRotation(NULL),
            settings(settings) {

          init();
     }

     //! Constructor for DBN-based preRotation
     //! \param chain Molecule chain
     //! \param dbn DBN object
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     template <typename DBN_TYPE>
     MoveBGS(CHAIN_TYPE *chain, DBN_TYPE *dbn, const Settings &settings = Settings(),
                RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "semilocal_bgs", settings, random_number_engine),
            BGSFactory(chain, dbn, settings, random_number_engine), 
            BGSRotation(NULL),
            settings(settings) {

          // When using a BGS move with DBN prior, the angle node is
          // flagged as fixed and is then unfixed every time a move is
          // made. This ensures that when the dbn is resynchronized
          // after a change has been made to the angles outside the
          // DBN, the hidden nodes will be resampled based on the new
          // angle values
          dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->fixed = true;

          init();
     }

     //! Copy constructor
     //! \param other Move_CRISP object from which the copy is made
     MoveBGS(const MoveBGS<CHAIN_TYPE,BGS_TYPE> &other)
          : MoveCommon(other),
            BGSFactory(other.BGSFactory),
            BGSRotation(NULL), 
            settings(other.settings)  {
          init();
     }

     //! Copy constructor - thread specified (multi-threading)
     //! \param other Move_CRISP object from which the copy is made     
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index of the thread
     //! \param chain Molecule chain
     MoveBGS(const MoveBGS<CHAIN_TYPE,BGS_TYPE> &other, 
                RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            BGSFactory(other.BGSFactory, random_number_engine, thread_index, chain),
            BGSRotation(NULL), 
            settings(other.settings) {
          init();
     }

     //! Apply local move
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether the move is successful
     bool apply(int start_index=0, int end_index=0) {

          // Import bgs namespace
          using namespace bgs;
          
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
          this->move_info = BGSFactory.create_move_info(start_index, end_index,NON_LOCAL);

          // Create temporary chain (backup in case of reject)
          if (this->chain_backup)
               delete this->chain_backup;
          this->chain_backup = new CHAIN_TYPE(*this->chain,
                                              start_index, end_index);

          // Setup the DofIterators

          DofIterator<ChainFB>::AngleSelectionEnum begin_dofs             = DofIterator<ChainFB>::DIHEDRAL_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum end_dofs               = DofIterator<ChainFB>::DIHEDRAL_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum backup_begin_dofs      = DofIterator<ChainFB>::DIHEDRAL_DOFS;
          DofIterator<ChainFB>::AngleSelectionEnum backup_end_dofs        = DofIterator<ChainFB>::DIHEDRAL_DOFS;
          

          if(settings.sample_omega) {
               begin_dofs  += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
               backup_begin_dofs  += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
               end_dofs  += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
               backup_end_dofs  += DofIterator<ChainFB>::N_DIHEDRAL + DofIterator<ChainFB>::CTERM_N_DIHEDRAL;
          }

          if(settings.sample_bond_angle) {
               begin_dofs += DofIterator<ChainFB>::BONDANGLE_DOFS;
               backup_begin_dofs += DofIterator<ChainFB>::BONDANGLE_DOFS;
               end_dofs += DofIterator<ChainFB>::BONDANGLE_DOFS;
               backup_end_dofs += DofIterator<ChainFB>::BONDANGLE_DOFS;
          }


          if(move_region==CTERM){
               begin_dofs        += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;
               backup_begin_dofs += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;
               end_dofs               += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;
               backup_end_dofs        += DofIterator<ChainFB>::CTERM_C_DIHEDRAL;

          } else {
               if (move_region==NTERM){
                    begin_dofs        += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
                    backup_begin_dofs += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
                    end_dofs               += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
                    backup_end_dofs        += DofIterator<ChainFB>::NTERM_CA_DIHEDRAL;
               }
          }

          // Construct Degree-of-freedom iterators
          DofIterator<ChainFB> begin((*this->chain)(start_index, CA),
                                     definitions::DIHEDRAL,
                                     begin_dofs);
          DofIterator<ChainFB> backup_begin((*this->chain_backup)(0, CA),
                                            definitions::DIHEDRAL,
                                            backup_begin_dofs);

               
          DofIterator<ChainFB> end((*this->chain)(end_index-1,C),
                                   definitions::DIHEDRAL,
                                   end_dofs);
          DofIterator<ChainFB> backup_end((*this->chain_backup)(end_index-1-start_index, C),
                                          definitions::DIHEDRAL,
                                          backup_end_dofs);
          
          dofs.clear();
          for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
               // Optionally skip proline phi angles
               if (settings.skip_proline_phi &&
                   (it.get_residue()->residue_type == definitions::PRO && 
                    it.get_atom()->atom_type == definitions::CA)) {
                    continue;
               }
               dofs.push_back(Dof(it.get_atom(),
                                  it.get_dof_type(),
                                  *it));
          }          

          dofs_backup.clear();
          for (DofIterator<ChainFB> it=backup_begin; it!=backup_end; ++it) {
               // Optionally skip proline phi angles
               if (settings.skip_proline_phi &&
                   (it.get_residue()->residue_type == definitions::PRO && 
                    it.get_atom()->atom_type == definitions::CA)) {
                    continue;
               }

               dofs_backup.push_back(Dof(it.get_atom(),
                                         it.get_dof_type(),
                                         *it));
          }          

          // Initialize move 
          BGSRotation = new BGS_TYPE(BGSFactory, move_region,
                                     end_index-start_index,
                                     dofs, dofs_backup);

          // Apply prerotation
          if(!BGSRotation->apply()){
               return (this->move_info->success = false);
          }
          
          if (this->move_info->success){
               // std::cout << "move successful \n";
               // Ensure that no bond angles are in illegal range
               for (unsigned int i=0; i<dofs.size(); ++i) {
                    //std::cout << it << "\n";
                    if (dofs[i].dof_type == definitions::ANGLE) {
                         // Check whether angle is close to 0 or pi
                         if((std::fabs(*dofs[i].value-M_PI) < settings.bond_angle_cutoff) || (std::fabs(*dofs[i].value) < settings.bond_angle_cutoff)){
                              return (this->move_info->success = false);
                         }
                    }
               }

               // Rigidly update the rest of the chain

               int fixed_direction = 1;
               this->chain->update_positions(dofs.back().atom->residue->index,
                                             this->chain->size(),
                                             fixed_direction);          
              
               // Calculate new positions non backbone atoms
               for (int i=start_index; i<this->chain->size(); i++) {
                    (*this->chain)[i].update_positions_non_backbone();
               }

               // Update state of preRotation object (e.g. DBN object)
               BGSRotation->update_state(this->chain);
          } 

          return this->move_info->success;
     }
     
     
     //! Calculate the log-bias introduced by the move
     //! \return Log bias
     double get_log_bias() {

          // Import crisp namespace
          using namespace bgs;
          
          double bias = 0.0;

          // Always evaluate bias of BGSRotation
          bias += BGSRotation->get_log_bias();

          return bias;
     }

     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();

          // Call accept in BGSRotation (update state in prior)
          BGSRotation->accept();
        
          
          // Cleanup backup chain
          delete this->chain_backup;
          this->chain_backup = NULL;

          // Cleanup pre and post rotation objects
          delete BGSRotation;
          BGSRotation = NULL;
     }     

     //! Reject last move
     void reject() {

          // Call base class reject method
          Move<CHAIN_TYPE>::reject();

          // Update positions and angles in the current chain based on the backup chain
          assert(this->chain_backup != NULL);

          AtomIterator<ChainFB,definitions::ALL> it_current(*this->chain, 
                                                            this->move_info->modified_positions[0].first,
                                                            this->move_info->modified_positions[0].second);
          AtomIterator<ChainFB,definitions::ALL> it_backup(*this->chain_backup);
          for (; !it_current.end() && !it_backup.end(); ++it_current,++it_backup) {
	       it_current->set_angle(it_backup->get_angle());
	       it_current->set_dihedral(it_backup->get_dihedral());
	       it_current->position = it_backup->position;
	  }

          // Cleanup backup chain
          delete this->chain_backup;
          this->chain_backup = NULL;

          // Call reject in BGSRotation (update state in prior)
          BGSRotation->reject();
                    
          // Cleanup pre and post rotation objects
          delete BGSRotation;
          BGSRotation = NULL;

     }          


     //! Specifies how move reacts to the knowledge that another moves has been executed
     void notify(MoveInfo *move_info) {
          BGSFactory.notify(move_info);
     }
     
};

}

#endif
