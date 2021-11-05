// move_cra.h --- Concerted rotation (with bond-angle modifications) move.
//                See: UlmSchneider, Jorgensen, J. Chem. Phys, 2003
//                Includes new strategy for post-ration solution selection
// Copyright (C) 2006-2008 Wouter Boomsma, Sandro Bottaro, Jesper Ferkinghoff-Borg
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


#ifndef MOVE_CRA
#define MOVE_CRA

#include "utils/random.h"
#include "../src/utils/matrix.h"
#include "moves/move_crisp.h"

namespace phaistos {

// Name space for CRISP internals
namespace cra {

//! Scale and tolerance parameters for the CRA algorithm
class Parameters {
public:
     
     //! specific scale factors for individual degrees of freedom.
     double scale[3][2];

     //! the tolerance of change for the postrotational part
     double tolerance[3][2];    

     //! parameter c1:overall scale -large value favours small changes on all  degrees of freedom in the prerotation
     double scale_common;  
    
     //! parameter c2:large value favours small changes of end-point
     double scale_local;       

     //! Constructor
     Parameters() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Indices
          int ANGLE = definitions::ANGLE;
          int DIHEDRAL = definitions::DIHEDRAL;

	  this->scale_common = 300.0;
	  this->scale_local = 8.0;

	  // Angle scale factors
	  scale[N][ANGLE] = 20.0;
	  scale[CA][ANGLE] = 20.0;
	  scale[C][ANGLE] = 20.0;

	  // Dihedral scale factors
	  scale[N][DIHEDRAL] = 40.0; // Omega
	  scale[CA][DIHEDRAL] = 1.0; // Phi
	  scale[C][DIHEDRAL] = 1.0;  // Psi

	  double c = M_PI/180.0;
	  // Angle tolerance
	  tolerance[N][ANGLE] = c*20.0;
	  tolerance[CA][ANGLE] = c*20.0;
	  tolerance[C][ANGLE] = c*20.0;

	  // Dihedral tolerance
	  tolerance[N][DIHEDRAL] = c*50.0;  // Omega
	  tolerance[CA][DIHEDRAL] = c*50.0; // Phi
	  tolerance[C][DIHEDRAL] = c*50.0;	 // Psi
     }
};


//! Prerotation class
// Contains functionality and state of last prerotation
class PreRotation {

     
     //! Random number generator  - uniform
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > random_generator_uniform_01;

     //! Random number generator  - Gaussian
     boost::variate_generator<RandomNumberEngine&, 
                              boost::normal_distribution<> > random_generator_gaussian;

public:

     //! CRA parameters
     Parameters &parameters;
    
     //! Chain object
     ChainFB *chain;

     //@{
     //! Degree-of-freedom iterators          
     DofIterator<ChainFB> &begin;
     DofIterator<ChainFB> &end;
     DofIterator<ChainFB> &breakpoint;
     //@}

     //! Number of degrees of freedom effected by prerotation
     int prerotation_dof_count;

     //! Choleski decomposition of the inverce covariance matrix (saved here for get_log_bias)
     Matrix Lt;

     //! Collection of random numbers
     Matrix dchi;

     //! Collection of angular variations 
     Matrix dphi;
     
     //! Constructor
     //! \param parameters CRA parameters 
     //! \param chain Molecule chain
     //! \param begin Begin dof iterator
     //! \param breakpoint Breakpoint dof iterator
     //! \param end End dof iterator
     //! \param random_number_engine Object from which random number generators can be created.
     PreRotation(Parameters &parameters, ChainFB *chain,
                 DofIterator<ChainFB> &begin, DofIterator<ChainFB> &breakpoint, DofIterator<ChainFB> &end,
                 RandomNumberEngine *random_number_engine=&random_global)
          : random_generator_uniform_01(*random_number_engine, boost::uniform_real<>(0,1)),
            random_generator_gaussian(*random_number_engine, boost::normal_distribution<>(0,(1.0/sqrt(2.0)))),
            parameters(parameters), begin(begin), end(end), breakpoint(breakpoint){
	  
	  this->chain = chain;

	  // Calculate number of degrees of freedom in prerotation
	  prerotation_dof_count = 0;
	  // for (deprecated::DofIterator it=begin; it!=breakpoint; ++it) {
	  for (DofIterator<ChainFB> it=begin; it!=breakpoint; ++it) {
	       prerotation_dof_count++;
	  }
     }

     //! Execute prerotation
     void apply() {

	  // Construct Lt matrix
	  construct_Lt(this->Lt);

	  // Generate Gaussian random values
	  generate_dchi(this->Lt.n_col);
	  
	  // Determine changes to angular degrees of freedom
	  dphi = solve_ax_b_triangular(Lt, dchi, 'U');

	  // Construct new conformation
	  // for (deprecated::DofIterator it=begin; it!=breakpoint; ++it) {
	  for (DofIterator<ChainFB> it=begin; it!=breakpoint; ++it) {
	       *it += dphi(it.offset, 0);
	       *it = fmod(*it+3*M_PI,(2*M_PI))-M_PI;
	       if (it.get_dof_type() == definitions::ANGLE) {
		    *it = std::fabs(*it);
	       }
	  }

	  // Update chain
	  chain->update_positions_backbone_segment(begin.get_residue()->index,
	        				   breakpoint.get_residue()->index+1);
     }


     //! Execute prerotation for end of chain
     void apply_endmove() {

	  // Calculate number of degrees of freedom
	  int size = 0;
	  // for (deprecated::DofIterator it=begin; it!=end; ++it)
	  for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
	       size++;
          }

	  // Only fill out diagonal in Lt
	  this->Lt = Matrix(size, size);
	  this->Lt.init_to_zero();
	  // for (deprecated::DofIterator it=begin; it!=end; ++it) {
	  for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
	       Lt(it.offset, it.offset) = sqrt((parameters.scale_common*
						parameters.scale[it.get_atom()->atom_type][it.get_dof_type()]));
	  }

	  // Generate Gaussian random values
	  generate_dchi(size);

	  // Determine changes to angular degrees of freedom
	  dphi = solve_ax_b_triangular(Lt, dchi, 'U');

	  // Construct new conformation
	  // for (deprecated::DofIterator it=begin; it!=end; ++it) {
	  for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
	       *it += dphi(it.offset,0);
	       *it = fmod(*it+3*M_PI,(2*M_PI))-M_PI;
	       if (it.get_dof_type() == definitions::ANGLE) {
		    *it = std::fabs(*it);
	       }
	  }

	  // Update chain positions
	  chain->update_positions_backbone(begin.get_residue()->index,
                                         end.get_residue()->index);

     }


     //! Construct the matrix Lt
     //! \param Lt Matrix Lt
     //! \param use_I1 Constrain CA position at breakpoint
     //! \param use_I2 Constrain N position at breakpoint
     //! \param use_I3 Constrain C position at breakpoint
     void construct_Lt(Matrix &Lt, bool use_I1=true, bool use_I2=false, bool use_I3=false) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

	  // Number of degrees of freedom in prerotation
	  int size = prerotation_dof_count;

	  std::vector<Vector_3D> pos;
	  if (use_I1)
	       pos.push_back((*chain)[breakpoint.get_residue()->index][CA]->position);
	  if (use_I2)
	       pos.push_back((*chain)[breakpoint.get_residue()->index][N]->position);
	  if (use_I3)
	       pos.push_back((*chain)[breakpoint.get_residue()->index][C]->position);
	  
	  std::vector<Matrix> I;
	  for (unsigned int index=0; index<pos.size(); index++) {

	       I.push_back(Matrix(size, size));
	       
	       Vector_3D *da_dphi = new Vector_3D[size];
	       // for (deprecated::DofIterator it=begin; it!=breakpoint; ++it) {
	       for (DofIterator<ChainFB> it=begin; it!=breakpoint; ++it) {
	  	    // The differential of a position vector around a
	  	    // rotation vector can be written as a cross product
	  	    // between the two vectors
	  	    da_dphi[it.offset] = it.get_rotation_vector() % (pos[index]-it.get_position_vector());
	       }

	       // I = da_dphi^2
	       for (int i=0; i<size; i++) {
	  	    for (int j=i; j<size; j++) {
	  		 I[index](i,j)=da_dphi[i]*da_dphi[j];
	  	    }
	       }
	       delete[] da_dphi;
	  }

	  // J: Linear transformation between normally distributed values
	  Matrix J(size,size);
	  for (int i=0; i<size; i++) {
	       double I_val = 0.0;
	       for (unsigned int index=0; index<I.size(); index++) {
		    I_val += I[index](i,i);
	       }
	       J(i,i)=parameters.scale_common*(1+parameters.scale_local*(I_val));
	       // J(i,i)=parameters.scale_common*(1+parameters.scale_local*I(i,i));
	       for (int j=i+1; j<size; j++) {
		    double I_val = 0.0;
		    for (unsigned int index=0; index<I.size(); index++) {
			 I_val += I[index](i,j);
		    }
 		    J(i,j) = parameters.scale_common*parameters.scale_local*I_val;
	       }
	  }

	  // Calculate "square root" of matrix J
	  Lt = J.cholesky('U');

	  // Scale to enforce different behaviour between bond angles and dihedrals
	  // for (deprecated::DofIterator it1=begin; it1!=breakpoint; ++it1) {
	  for (DofIterator<ChainFB> it1=begin; it1!=breakpoint; ++it1) {
	       // for (deprecated::DofIterator it2=begin; it2!=breakpoint; ++it2) {
	       for (DofIterator<ChainFB> it2=begin; it2!=breakpoint; ++it2) {
	  	    Lt(it1.offset, it2.offset) = (Lt(it1.offset, it2.offset) *
	  					  parameters.scale[it2.get_atom()->atom_type][it2.get_dof_type()]);
	       }
	  }
     }

     //! Generate a vector of random normally distributed samples
     //! \param size size of the vector
     void generate_dchi(int size) {

	  this->dchi = Matrix(size, 1);
	  for (int i=0; i<size; i++) {
	       this->dchi(i,0) = random_generator_gaussian();
	  }	  
     }

     //! Calculate bias for prerotation move (a->b)
     //! \return Log bias introduced by the proposed move
     double calculate_bias() {
	  double logdetL = 0.0;
	  double d_squared = 0.0;
	  for (int i=0;i<this->Lt.n_col;i++) {
	       logdetL += std::log(Lt(i,i));
	       d_squared += Math<double>::sqr(dchi(i,0));
	  }

	  // Calculate biasing probability a->b
	  return logdetL - d_squared;
     }     

     //! Calculate bias for reverse move (b->a)
     //! NOTE: Assumes that apply has been called
     //! \return Log bias introduced by the reverse move
     double calculate_bias_reverse() {

	  // Calculate Lt'
	  Matrix Lt_prime;
	  construct_Lt(Lt_prime);

	  // Calculate Chi'
	  Vector_nD dchi_prime = Lt_prime * dphi.get_column(0);

	  double logdetL_prime = 0.0;
	  double d_squared_prime = 0.0;
	  for (int i=0; i<prerotation_dof_count; i++) {
	       logdetL_prime += std::log(Lt_prime(i,i));	       
	       d_squared_prime += Math<double>::sqr(dchi_prime[i]);
	  }

	  // Calculate biasing probability a->b
	  return logdetL_prime - d_squared_prime;
     }	  

     //! Accept last move
     void accept() {
     }

     //! Reject last move
     void reject() {
     }     
};

}



//! Move CRA class
template <typename CHAIN_TYPE>
class MoveCRA: public MoveCommon<MoveCRA<CHAIN_TYPE>, CHAIN_TYPE> {
private:

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveCRA<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Pre-rotation object
     cra::PreRotation *preRotation;
     
     //! Post-rotation object (from crisp namespace)
     crisp::PostRotation *postRotation;

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
     friend class MoveStatisticsLocalMove<CHAIN_TYPE, MoveCRA>;
     
public:

     //! CRA parameters
     cra::Parameters parameters;

     //! Move Settings
     const class Settings: public Move<CHAIN_TYPE>::Settings {
     public:

          //! Whether bond angle bias should be divided out
          bool implicit_energy;

          //! Constructor
          Settings(bool implicit_energy=false)
               : Move<CHAIN_TYPE>::Settings(5,5),
                 implicit_energy(implicit_energy) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "implicit-energy:" << settings.implicit_energy << "\n";
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                    
     } settings;
     

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings CRA Settings
     //! \param random_number_engine Object from which random number generators can be created.
     // Note: If resample_hidden_nodes is set to true, the relevant
     // part of the hidden node sequence will be resampled prior to
     // loop closure. This move is ALWAYS accepted and will
     // not be rejected even if reject is called subsequently.
     MoveCRA(ChainFB *chain, const Settings &settings = Settings(),
             RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "cra", settings, random_number_engine),
            settings(settings) {

          init();
     }

     //! Initializer
     void init() {

          // Statistics object
          delete this->statistics;  // Remove existing statistics object
          this->statistics = new MoveStatisticsLocalMove<CHAIN_TYPE, MoveCRA>(this);                    
     }

     //! Copy constructor
     //! \param other Move CRA from which the copy is made
     MoveCRA(const MoveCRA<CHAIN_TYPE> &other)
          : MoveCommon(other),
            settings(other.settings) {
          init();
     }

     
     //! Copy constructor
     //! \param other Move CRA from which the copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain
     MoveCRA(const MoveCRA<CHAIN_TYPE> &other, 
             RandomNumberEngine *random_number_engine, int thread_index,
              CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
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

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Call base class apply method
          Move<CHAIN_TYPE>::apply(start_index, end_index);          
          
	  // Move range from -(length-1) to size+length-1
	  if (start_index == end_index) {
               const std::pair<int,int> &region = this->random_move_region();
               int length = this->random_move_length(region);
	       this->random_move_range(region, length, -(length-2), this->chain->size()+length-2, &start_index, &end_index);

	       // // only local moves (no endmoves)
               //this->random_move_range(length, 1, chain->size()-1, &start_index, &end_index);
	  }

	  

	  // Determine whether to do loop closure or not
	  move_region = INTERNAL;
	  if (start_index<1) {
               move_region = NTERM;
	       start_index = 0;
	  }
	  if (end_index > this->chain->size()-1) {
               move_region = CTERM;
	       end_index = this->chain->size();
	  }

          // Set which positions are updated
          this->move_info = (new MoveInfo(LOCAL))->add_info(std::make_pair(start_index, end_index),
                                                            std::make_pair(start_index, end_index));

	  // Create temporary chain (backup in case of reject)
	  if (this->chain_backup) 
	       delete this->chain_backup;
	  this->chain_backup = new CHAIN_TYPE(*this->chain, 
                                              start_index, end_index);

	  int breakpoint_index = end_index-2;
	  // If no loop closure is done, the breakpoint iterator will be equal to the end iterator
	  if (move_region != INTERNAL) {
	       breakpoint_index = end_index-1;
	  }


          // Indices
          int ANGLE = definitions::ANGLE;
          int DIHEDRAL = definitions::DIHEDRAL;

          if(move_region == CTERM){

               begin = new DofIterator<ChainFB>((*this->chain)(start_index, CA), DIHEDRAL,
                                       (DofIterator<ChainFB>::STANDARD_DOFS +
                                        DofIterator<ChainFB>::CTERM_C_DIHEDRAL + 
                                        DofIterator<ChainFB>::CTERM_C_ANGLE));
               backup_begin = new DofIterator<ChainFB>((*this->chain_backup)(0, CA), DIHEDRAL,
                                              (DofIterator<ChainFB>::STANDARD_DOFS + 
                                               DofIterator<ChainFB>::NTERM_CA_DIHEDRAL +
                                               DofIterator<ChainFB>::CTERM_C_DIHEDRAL + 
                                               DofIterator<ChainFB>::CTERM_C_ANGLE));

               breakpoint = new DofIterator<ChainFB>((*this->chain)(breakpoint_index, C),ANGLE,
                                            (DofIterator<ChainFB>::STANDARD_DOFS +
                                             DofIterator<ChainFB>::CTERM_C_DIHEDRAL + 
                                             DofIterator<ChainFB>::CTERM_C_ANGLE));
               backup_breakpoint = new DofIterator<ChainFB>((*this->chain_backup)(breakpoint_index-start_index, C),
                                                   ANGLE,
                                                   (DofIterator<ChainFB>::STANDARD_DOFS +
                                                    DofIterator<ChainFB>::CTERM_C_DIHEDRAL + 
                                                    DofIterator<ChainFB>::CTERM_C_ANGLE));

               end = new DofIterator<ChainFB>((*this->chain)(end_index-1, C), ANGLE,
                                     (DofIterator<ChainFB>::STANDARD_DOFS +
                                      DofIterator<ChainFB>::CTERM_C_DIHEDRAL + 
                                      DofIterator<ChainFB>::CTERM_C_ANGLE));
               backup_end = new DofIterator<ChainFB>((*this->chain_backup)(end_index-1-start_index, C), ANGLE,
                                            (DofIterator<ChainFB>::STANDARD_DOFS +
                                             DofIterator<ChainFB>::CTERM_C_DIHEDRAL + 
                                             DofIterator<ChainFB>::CTERM_C_ANGLE));

          } else if (move_region == NTERM){

               begin = new DofIterator<ChainFB>((*this->chain)(start_index, CA), DIHEDRAL,
                                                (DofIterator<ChainFB>::STANDARD_DOFS + 
                                                 DofIterator<ChainFB>::NTERM_CA_DIHEDRAL));
               backup_begin = new DofIterator<ChainFB>((*this->chain_backup)(0, CA), DIHEDRAL,
                                                       (DofIterator<ChainFB>::STANDARD_DOFS + 
                                                        DofIterator<ChainFB>::NTERM_CA_DIHEDRAL));

               breakpoint = new DofIterator<ChainFB>((*this->chain)(breakpoint_index, CA), ANGLE);
               backup_breakpoint = new DofIterator<ChainFB>((*this->chain_backup)(breakpoint_index-start_index, CA), ANGLE);

               end = new DofIterator<ChainFB>((*this->chain)(end_index-1, CA), ANGLE);
               backup_end = new DofIterator<ChainFB>((*this->chain_backup)(end_index-1-start_index, CA), ANGLE);

          } else {
               begin = new DofIterator<ChainFB>((*this->chain)(start_index, CA));
               if (start_index == 0) {
                    backup_begin = new DofIterator<ChainFB>((*this->chain_backup)(0, CA), DIHEDRAL,
                                                   DofIterator<ChainFB>::STANDARD_DOFS);
               } else {
                    backup_begin = new DofIterator<ChainFB>((*this->chain_backup)(0, CA), DIHEDRAL,
                                                            (DofIterator<ChainFB>::STANDARD_DOFS + 
                                                             DofIterator<ChainFB>::NTERM_CA_DIHEDRAL));
               }

               breakpoint = new DofIterator<ChainFB>((*this->chain)(breakpoint_index, CA));
               backup_breakpoint = new DofIterator<ChainFB>((*this->chain_backup)(breakpoint_index-start_index, CA));

               end = new DofIterator<ChainFB>((*this->chain)(end_index-1, CA), ANGLE);
               backup_end = new DofIterator<ChainFB>((*this->chain_backup)(end_index-1-start_index, CA), ANGLE);
          }


	  // Initialize prerotation
	  preRotation = new cra::PreRotation(parameters, this->chain, *begin, *breakpoint, *end, this->random_number_engine);

	  // Initialize postrotation
	  postRotation = new crisp::PostRotation(this->chain, *breakpoint, *end, *backup_breakpoint, *backup_end);

	  if (move_region == INTERNAL) {
	       preRotation->apply();
	       this->move_info->success = postRotation->apply();

	  } else {
               preRotation->apply_endmove();
	  }

	  if (this->move_info->success)
	       finalize(this->chain);

	  return this->move_info->success;
     }


     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return Log bias
     double get_log_bias() {

	  double bias = 0.0;
	  
	  // Determine whether last move was a loop closure or an end move
	  if (move_region == definitions::INTERNAL) {

	       double log_Jacobian_before = postRotation->calculate_log_Jacobian_before();
	       double log_Jacobian_after = postRotation->calculate_log_Jacobian_after();

	       double log_prerotation_bias = preRotation->calculate_bias();
	       double log_prerotation_bias_reverse = preRotation->calculate_bias_reverse();

	       bias += (log_Jacobian_after + log_prerotation_bias_reverse -
			log_Jacobian_before - log_prerotation_bias);

	  // } else {

          //      // The distribution over bondangles should follow a sine distribution
          //      // We compensate for the bias here
          //      DofIterator<ChainFB> it=*begin;
          //      DofIterator<ChainFB> it_backup=*backup_begin;
          //      for (; it!=*end && it_backup!=*backup_end; ++it,++it_backup) {
                    
          //           if ((it.get_atom()->atom_type != it_backup.get_atom()->atom_type) || 
          //               (it.get_atom()->residue->residue_type != it_backup.get_atom()->residue->residue_type)) {
          //                assert(false);
          //           }

          //           if (it.get_dof_type() == DofIterator<ChainFB>::ANGLE) {
          //                bias += log(sin(*it)) - log(sin(*it_backup));
          //           }
          //      }               
          }

	  // Include a bias directing the bond angles to their native
	  // values. This will act as an implicit energy term.
	  if (settings.implicit_energy) {

	       double angle_LL_after = get_log_likelihood_bond_angles(*begin, *end);
	       double angle_LL_before = get_log_likelihood_bond_angles(*backup_begin, *backup_end);

	       bias += angle_LL_after - angle_LL_before;
	  }
	  
	  return bias;
     }
     

     //! Calculate the likelihood of the current bond angles
     //! The bond angle energies are Gaussians around their ideal Engh-Huber values
     //! \param begin Begin dof iterator
     //! \param end End dof iterator
     //! \return Log likelyhood bond angles
     double get_log_likelihood_bond_angles(DofIterator<ChainFB> &begin, DofIterator<ChainFB> &end) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

	  double log_probability = 0.0;

	  for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
	       if (it.get_dof_type() == definitions::ANGLE) {
		    Atom *atom=it.get_atom();
		    Atom *nextAtom = atom->get_neighbour(+1, BACKBONE);
		    Atom *prevAtom = atom->get_neighbour(-1, BACKBONE);
		    double stdDev = 0.0;
		    double mean=bond_angle_constants(prevAtom->atom_type,
						   atom->atom_type,
						   nextAtom->atom_type,
						   it.get_residue()->residue_type, &stdDev);

		    double angle = *it;

		    // Accept with higher tolerance - for ergodicity reasons
		    // stdDev *= 3.0;
		    
		    // Gaussian density 
		    log_probability += -0.5*(Math<double>::sqr(angle-mean)/Math<double>::sqr(stdDev)) - std::log(stdDev*sqrt(2.0*M_PI));
	       }
	  }

	  return log_probability;	  
     }
     
     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();
          
	  // std::cout << "accept\n";
	  delete this->chain_backup;
	  this->chain_backup = NULL;

          // Call accept in preRotation (update state in prior)
          preRotation->accept();
          
          // Call accept in postRotation
          if (move_region == definitions::INTERNAL) {
               postRotation->accept();
          }          

	  delete preRotation;
	  delete postRotation;
	  delete begin;
	  delete breakpoint;
	  delete end;
	  delete backup_begin;
	  delete backup_breakpoint;
	  delete backup_end;

     }


     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();
                              
	  // std::cout << "reject\n";
	  // Update positions and angles in the current chain based on the test chain
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

          // Call reject in preRotation (update state in prior)
          preRotation->reject();
                    
          // Call reject in postRotation
          if (move_region == definitions::INTERNAL) {
               postRotation->reject();
          }

	  delete this->chain_backup;
	  this->chain_backup = NULL;
	  
	  delete preRotation;
	  delete postRotation;
	  delete begin;
	  delete breakpoint;
	  delete end;
	  delete backup_begin;
	  delete backup_breakpoint;
	  delete backup_end;

     }


     //! Calculate new positions for non backbone atoms (FB chain)
     //! \param chain Molecule chain
     void finalize(ChainFB *chain) {
	  for (int i=this->move_info->modified_positions[0].first; i<this->move_info->modified_positions[0].second; i++) {
	       (*chain)[i].update_positions_non_backbone();	       
	  }
     }

     //! Calculate new positions for non backbone atoms (CA chain)
     //! \param chain Molecule chain
     void finalize(ChainCA *chain) {}
     
};

}

#endif
