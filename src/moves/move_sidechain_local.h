// move_sidechain_local.h --- Move: Local updates of sidechain chi angles using constraints on hydrogen-bond donors and acceptors
// Copyright (C) 2006-2010 Wouter Boomsma
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


#ifndef MOVE_SIDECHAIN_LOCAL_H
#define MOVE_SIDECHAIN_LOCAL_H

#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"
#include "move.h"

namespace phaistos {

//! Local enum defining possible move types
enum MoveSidechainLocalMode {SIMPLE,SCALE_SIGMA,SELECT_DOFS,CONSTRAIN_ALL_ENDPOINTS,CONSTRAIN_ONE_ENDPOINT,MODE_SIZE};

//! Names associated with local enum
static const std::string MoveScLocalModeNames[] = {"simple","scale-sigma","select-dofs","constrain-all-endpoints","constrain-one-endpoint"};

//! Input MoveSidechainLocalMode from stream
inline std::istream &operator>>(std::istream &input, MoveSidechainLocalMode &rm) {
     std::string raw_string;
     input >> raw_string;

     for (unsigned int i=0; i<MODE_SIZE; ++i) {
          if (raw_string == MoveScLocalModeNames[i]) {
               rm = MoveSidechainLocalMode(i);
          }
     }
     return input;
}

//! Output MoveSidechainLocalMode to stream
inline std::ostream &operator<<(std::ostream &o, const MoveSidechainLocalMode &rm) {
     o << MoveScLocalModeNames[static_cast<unsigned int>(rm)];
     return o;
}




//! Local updates of sidechain chi angles using constraints on hydrogen-bond donors and acceptors
template <typename CHAIN_TYPE>
class MoveSidechainLocal: public MoveCommon<MoveSidechainLocal<CHAIN_TYPE>,
                                                               CHAIN_TYPE> {

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveSidechainLocal<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;

     //! Container for backed up major dofs (mainly chi angles)
     std::vector<std::vector<double> > major_dof_value_backup;

     //! Container for backed up minor dofs (mainly bond angles) 
     std::vector<std::vector<double> > minor_dof_value_backup;

     //! Gaussian random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::normal_distribution<> > random_generator_gaussian;


     //@{
     //! Matrices used by constrained endpoint move
     std::vector<Matrix> Lt;
     std::vector<Matrix> dphi;
     std::vector<Matrix> precision_matrix;
     std::vector<int> constrained_atom_index;
     //@}


     //! Calculate the derivative of the endpoint with respect to the current dof
     //! \param res Residue pointer
     //! \param (atom-type, dof-type pair)
     //! \param constrained_atom Pointer to endpoint atom
     //! \return derivative (vector)
     Vector_3D calculate_da_dphi(typename CHAIN_TYPE::Residue *res,
                                 const std::pair<definitions::AtomEnum,definitions::AngleEnum> &dof,
                                 const Atom *constrained_atom) {
          Atom *dof_atom = (*res)[dof.first];
          Atom *dof_atom_prev = dof_atom->get_neighbour<-1>(definitions::POSITIONING);
          Atom *dof_atom_prev_prev = dof_atom->get_neighbour<-2>(definitions::POSITIONING);

          Vector_3D rotation_vector;
          if (dof.second == definitions::DIHEDRAL)
               rotation_vector = (dof_atom_prev->position - 
                                  dof_atom_prev_prev->position).normalize();
          else {
               Vector_3D v1 = dof_atom_prev->position - dof_atom_prev_prev->position;
               Atom *dof_atom_prev_real = dof_atom->get_neighbour<-1>(definitions::BACKBONE);
               Vector_3D v2 = dof_atom->position - dof_atom_prev_real->position;
               rotation_vector = (v2%v1).normalize();
          }
                    
          Vector_3D arm = (constrained_atom->position - dof_atom_prev->position);

          // The differential of a position vector around a
          // rotation vector can be written as a cross product
          // between the two vectors
          return rotation_vector % arm;
     }


     //! Computer the precision matrix (inverse covariance matrix)
     //! \param res Residue pointer
     //! \param precision_matrix Target inverse covariance matrix
     //! \param selected_constrained_atom_index Optionally specify a particular constraint atom (rather than picking one randomly)
     //! \return whether constraint was imposed (this is false if there are no constrained atoms in a residue)
     bool construct_precision_matrix(typename CHAIN_TYPE::Residue *res,
                                     Matrix &precision_matrix,
                                     int selected_constrained_atom_index=-1) {
          
          // Import protein definitions (such as residue names)
          using namespace definitions;

          std::vector<Matrix> I;
          std::vector<std::vector<Vector_3D> > da_dphi;

          // Check if there are constrained atoms for this residue type
          if (settings.constrained_atoms.count(res->residue_type)) {

               // Retrieve constrained_atom information
               const std::vector<std::pair<AtomEnum,std::pair<std::vector<
               std::pair<AtomEnum,AngleEnum> >,std::vector<std::pair<AtomEnum,AngleEnum> > > > > &data = 
                    settings.constrained_atoms.find(res->residue_type)->second;

               int size = precision_matrix.n_row;

               int start_index = 0;
               int end_index = data.size();

               // In mode=CONSTRAIN_ONE_ENDPOINT, we pick a random constraint atom
               if (settings.mode == CONSTRAIN_ONE_ENDPOINT) {
                    if (selected_constrained_atom_index == -1) {
                         boost::variate_generator<RandomNumberEngine&,
                              boost::uniform_int<> >
                              random_generator_uniform_int(*this->random_number_engine,
                                                           boost::uniform_int<>(0,
                                                                                end_index-1));
                         selected_constrained_atom_index = random_generator_uniform_int();
                         
                         // Save index - so the same index can be used in the reverse step
                         constrained_atom_index.push_back(selected_constrained_atom_index);
                    }

                    start_index = selected_constrained_atom_index;
                    end_index = start_index+1;
               } else {
                    constrained_atom_index.push_back(-1);
               }

               // Relates the index that the atom has in the residue to the index that it has 
               // in the precision matrix - both for angles and dihedrals
               std::vector<std::vector<int> > dof_indices(res->atoms.size(), std::vector<int>(2,-1));

               // Initialize precision matrix diagonal using sigmas
               int offset = 0;
               if (settings.sample_major_dofs) {
                    for (unsigned int j=0; j<res->chi_atoms.size(); ++j) {
                         precision_matrix(offset+j,offset+j) = 1.0/(Math<double>::sqr(settings.sigma_major_dofs));

                         // Register the index in precision_matrix of all dofs
                         AtomEnum atom_type = res->chi_atoms[j].first->atom_type;
                         dof_indices[res->atom_index[atom_type]][DIHEDRAL] = offset+(int)j;
                    }
                    offset += res->chi_atoms.size();
               }
               if (settings.sample_minor_dofs) {
                    for (unsigned int j=0; j<res->minor_dof_atoms.size(); ++j) {
                         // if (offset+j >= precision_matrix.n_row) {
                         //      std::cout << offset << " " << j << " " << precision_matrix.n_row << "\n";
                         //      assert(false);
                         // }
                         precision_matrix(offset+j,offset+j) = 1.0/(Math<double>::sqr(settings.sigma_minor_dofs));

                         // Register the index in precision_matrix of all dofs
                         AtomEnum atom_type = res->minor_dof_atoms[j].first->atom_type;
                         AngleEnum angle_type = res->minor_dof_atoms[j].second;
                         dof_indices[res->atom_index[atom_type]][angle_type] = offset+(int)j;
                    }
               }
               
               // Iterate over constrained atoms
               for (int index=start_index; index<end_index; ++index) {

                    // Construct matrices
                    I.push_back(Matrix(size, size, 0.0));
                    da_dphi.push_back(std::vector<Vector_3D>(size,Vector_3D(0,0,0)));

                    // Extract information about constrained atom and constrained chi angles from data structure
                    Atom *constrained_atom = (*res)[data[index].first];
                    
                    // Iterate over constrained chi angles
                    if (settings.sample_major_dofs) {
                         const std::vector<std::pair<AtomEnum,AngleEnum> > &constrained_major_dofs = data[index].second.first;
                         for (unsigned int i=0; i<constrained_major_dofs.size(); ++i) {
                              
                              // Find index in precision matrix
                              AtomEnum atom_type = constrained_major_dofs[i].first;
                              AngleEnum angle_type = constrained_major_dofs[i].second;                              
                              int precision_matrix_index = dof_indices[res->atom_index[atom_type]][angle_type];

                              // Calculate the derivative of the endpoint with respect to the current dof
                              da_dphi.back()[precision_matrix_index] = 
                                   calculate_da_dphi(res, constrained_major_dofs[i], constrained_atom);
                              
                         }
                    }

                    // Iterate over constrained chi angles
                    if (settings.sample_minor_dofs) {
                         const std::vector<std::pair<AtomEnum,AngleEnum> > &constrained_minor_dofs = data[index].second.second;
                         for (unsigned int i=0; i<constrained_minor_dofs.size(); ++i) {

                              // Find index in precision matrix
                              AtomEnum atom_type = constrained_minor_dofs[i].first;
                              AngleEnum angle_type = constrained_minor_dofs[i].second;                              
                              int precision_matrix_index = dof_indices[res->atom_index[atom_type]][angle_type];

                              // Calculate the derivative of the endpoint with respect to the current dof
                              da_dphi.back()[precision_matrix_index] = 
                                   calculate_da_dphi(res, constrained_minor_dofs[i], constrained_atom);
                         }
                    }

                    // I = da_dphi^2
                    for (int i=0; i<size; i++) {
                         for (int j=0; j<size; j++) {
                              I.back()(i,j)=da_dphi.back()[i]*da_dphi.back()[j];
                         }
                    }

                    if (settings.debug)
                         std::cout << "I: " << I << "\n";
               }

               // Add constraints
               for (unsigned int j=0; j<I.size(); ++j) {
                    precision_matrix += (I[j] * (settings.lagrange_multiplier/(double)I.size()));
               }

               return true;
          } else {
               constrained_atom_index.push_back(-1);
               return false;
          }
     }


     //! Sampling in case were no endpoint constrained is requested
     //! If mode==SCALE_SIGMA, the standard deviation is scaled so that
     //! dofs further away from CA have larger variance. If mode==SELECT_DOFS
     //! the dofs are chosen at random according to a weight given by their
     //! distance to the CA atom.
     //!
     //! \param res Residue pointer
     //! \param min_distance_to_CA The minimum distance to CA among all dofs in the current residue
     //! \param max_distance_to_CA The maximum distance to CA among all dofs in the current residue
     //! \param dof_atoms vector of [atom-pointer,dof-type] pairs
     //! \param standard deviation
     //! \param selection_cutoff_value Random sample to compare against - enforcing the weights on the different dofs
     //! \param dof_values Target vector of dof values
     void sample_unconstrained_endpoint(typename CHAIN_TYPE::Residue *res,
                                        int min_distance_to_CA, 
                                        int max_distance_to_CA,
                                        const std::vector<std::pair<Atom*,definitions::AngleEnum> > &dof_atoms, 
                                        double sigma,
                                        double selection_cutoff_value,
                                        std::vector<double> &dof_values) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          for (unsigned int j=0; j<dof_atoms.size(); ++j) {

               // Calculate distance from current chi angle atom to CA
               int distance_to_CA = chain_distance<CHAIN_TYPE>(dof_atoms[j].first, (*res)[CA]);

               // Calculate weight - used both by SCALE_SIGMA (for weighting) and SELECT_DOFS (for selection)
               double weight = (1+distance_to_CA-min_distance_to_CA)/(double)(1+max_distance_to_CA-min_distance_to_CA);

               // Sample stepsize
               double delta = random_generator_gaussian()*sigma;
               if (settings.mode == SCALE_SIGMA) {
                    // Delta is distributed with stddev=sigma - so we simply scale
                    delta *= weight;
                    dof_values[j] = fmod(dof_values[j] + 3*M_PI + delta, 2*M_PI) - M_PI;
               } else {
                    if (selection_cutoff_value < weight) {
                         dof_values[j] = fmod(dof_values[j] + 3*M_PI + delta, 2*M_PI) - M_PI;
                    } 
               }

               // Take absolute value if current dof is an angle
               if (dof_atoms[j].second == ANGLE) {
                    dof_values[j] = std::fabs(dof_values[j]);
               }
          }
     }


     //! Calculate minimum and maximum distance to CA along the chain for all dofs
     //! \param res Residue pointer
     //! \param dof_atoms vector of [atom-pointer,dof-type] pairs
     //! \param min_distance_to_CA Target variable in which to write minimum distance
     //! \param max_distance_to_CA Target variable in which to write maximum distance
     void get_min_max_chain_distance(typename CHAIN_TYPE::Residue *res,
                               std::vector<std::pair<Atom*,definitions::AngleEnum> > &dof_atoms,
                               int &min_distance_to_CA, int &max_distance_to_CA) {
          
          for (unsigned int i=0; i<dof_atoms.size(); ++i) {
               int distance = chain_distance<CHAIN_TYPE>(dof_atoms[i].first, (*res)[definitions::CA]);
               if (distance<min_distance_to_CA)
                    min_distance_to_CA = distance;
               if (distance > max_distance_to_CA)
                    max_distance_to_CA = distance;
          }
     }


     //! Calculate bias for reverse move using the constrained end mode
     //! \return log-probability
     double calculate_log_bias_forward() {

          double log_bias = 0.0;
          for (int index=this->move_info->modified_angles_start; 
               index<this->move_info->modified_angles_end; ++index) {

               int i = index-this->move_info->modified_angles_start; 

               double log_det_L = 0.0;
               double sum = 0.0;
               int size = Lt[i].n_row;
               for (int j=0; j<size; ++j) {
                    log_det_L += std::log(Lt[i](j,j));
                    sum += 0.5*dphi[i](j,0) * precision_matrix[i](j,j) * dphi[i](j,0);
                    for (int k=j+1; k<size; k++) {
                         sum += dphi[i](j,0) * precision_matrix[i](j,k) * dphi[i](k,0);
                    }
               }

               log_bias += log_det_L - sum;
          }


          return log_bias;
     }

     //! Calculate bias for reverse move using the constrained end mode
     //! \return log-probability
     double calculate_log_bias_reverse() {

          double log_bias = 0.0;
          for (int index=this->move_info->modified_angles_start; 
               index<this->move_info->modified_angles_end; ++index) {

               // Get residue
               typename CHAIN_TYPE::Residue *res = &(*this->chain)[index];

               int i = index-this->move_info->modified_angles_start; 

               int size = 0;
               if (settings.sample_major_dofs)
                    size += res->chi_atoms.size();
               if (settings.sample_minor_dofs)
                    size += res->minor_dof_atoms.size();
               Matrix precision_matrix_prime(size, size, 0.0);
               if (construct_precision_matrix(res, precision_matrix_prime, constrained_atom_index[i])) {

                    // Calculate "square root" of matrix W_tilde
                    Matrix Lt_prime(precision_matrix_prime.cholesky('U'));

                    double log_det_L = 0.0;
                    double sum = 0.0;
                    int size = Lt_prime.n_row;
                    for (int j=0; j<size; ++j) {
                         log_det_L += std::log(Lt_prime(j,j));
                         sum += 0.5*(-dphi[i](j,0)) * precision_matrix_prime(j,j) * (-dphi[i](j,0));
                         for (int k=j+1; k<size; k++) {
                              sum += (-dphi[i](j,0)) * precision_matrix_prime(j,k) * (-dphi[i](k,0));
                         }
                    }
                    log_bias += log_det_L - sum;
               }
          }
          return log_bias;
     }



public:

     //! Local Settings class.
     const class Settings: public Move<ChainFB>::Settings {
     public:

          //! MoveSidechainLocalMode: SIMPLE|SCALE_SIGMA|SELECT_DOFS|CONSTRAIN_ALL_ENDPOINTS|CONSTRAIN_ONE_ENDPOINT
          MoveSidechainLocalMode mode;

          //! Default standard deviation of proposals - major dofs
          double sigma_major_dofs;

          //! Default standard deviation of proposals - minor dofs
          double sigma_minor_dofs;

          //! Whether to resample major dofs (mainly chi angles)
          bool sample_major_dofs;

          //! Whether to resample minor dofs (mainly bondangles)
          bool sample_minor_dofs;

          //! Lagrange multiplier used in CONSTRAIN_ENDPOINT mode
          double lagrange_multiplier;

          //! Whether to skip prolines (prolines introduce a change in bond length which must be taken into account by the forcefield)
          bool skip_proline;          

          //! Which atoms should be constrained in the different sidechains
          std::map<definitions::ResidueEnum, 
                   std::vector<
                        std::pair<
                             definitions::AtomEnum,
                             std::pair<
                                  std::vector<
                                       std::pair<definitions::AtomEnum,
                                                 definitions::AngleEnum> >,
                                  std::vector<
                                       std::pair<definitions::AtomEnum,
                                                 definitions::AngleEnum> > > > > > constrained_atoms;

          //! Constructor
          Settings(MoveSidechainLocalMode mode=SIMPLE,
                   double sigma_major_dofs=2.0/180.0*M_PI,
                   double sigma_minor_dofs=0.2/180.0*M_PI,
                   bool sample_major_dofs=true,
                   bool sample_minor_dofs=false,
                   double lagrange_multiplier=200,
                   bool skip_proline=true)
               : Move<ChainFB>::Settings(1,1),
                 mode(mode),
                 sigma_major_dofs(sigma_major_dofs),
                 sigma_minor_dofs(sigma_minor_dofs),
                 sample_major_dofs(sample_major_dofs),
                 sample_minor_dofs(sample_minor_dofs),
                 lagrange_multiplier(lagrange_multiplier),
                 skip_proline(skip_proline) {

               using namespace std;
               using namespace vector_utils;
               using namespace boost::assign;

               // Import protein definitions (such as residue names)
               using namespace definitions;

               // Set which atoms are constrained for each residue type
               // format: (ResidueType, [(constrained-atom1, 
               //                          ([(major_dof_atom1, major_dof_type1),
               //                            (major_dof_atom2, major_dof_type2),
               //                            ...], 
               //                           [(minor_dof_atom1, minor_dof_type1),
               //                            (minor_dof_atom2, minor_dof_type2),
               //                            ...])),
               //                        (constrained-atom2, 
               //                          ...
               constrained_atoms.insert(
                    make_pair(
                         CYS,
                         make_vector(
                              make_pair(
                                   SG,                           // restrained endpoint atom
                                   make_pair(
                                        make_vector(make_pair(SG,DIHEDRAL)),   // effected dihedrals
                                        make_vector(make_pair(CB,ANGLE),       // effected minor dofs
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(SG,ANGLE))))))); 
               constrained_atoms.insert(
                    make_pair(
                         ASP,
                         make_vector(
                              make_pair(
                                   OD1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(OD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(OD1,ANGLE),
                                                    make_pair(OD2,DIHEDRAL)))),
                              make_pair(
                                   OD2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(OD1,DIHEDRAL)), //OD1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(OD2,ANGLE),
                                                    make_pair(OD2,DIHEDRAL))))))); 
               constrained_atoms.insert(
                    make_pair(
                         GLU,
                         make_vector(
                              make_pair(
                                   OE1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(OE1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(OE1,ANGLE),
                                                    make_pair(OE2,DIHEDRAL)))),
                              make_pair(
                                   OE2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(OE1,DIHEDRAL)), //OE1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(OE2,ANGLE),
                                                    make_pair(OE2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         PHE,
                         make_vector(
                              make_pair(
                                   CZ, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD1,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         HIS,
                         make_vector(
                              make_pair(
                                   NE2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(ND1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(ND1,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         ILE,
                         make_vector(
                              make_pair(
                                   CD1, 
                                   make_pair(
                                        make_vector(make_pair(CG1,DIHEDRAL),
                                                    make_pair(CD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG1,ANGLE),
                                                    make_pair(CG2,DIHEDRAL),
                                                    make_pair(CD1,ANGLE)))),
                              make_pair(
                                   CG2, 
                                   make_pair(
                                        make_vector(make_pair(CG1,DIHEDRAL)),//CG11 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG2,ANGLE),
                                                    make_pair(CG2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         LYS,
                         make_vector(
                              make_pair(
                                   NZ, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(CE,DIHEDRAL),
                                                    make_pair(NZ,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(CE,ANGLE),
                                                    make_pair(NZ,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         LEU,
                         make_vector(
                              make_pair(
                                   CD1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD1,ANGLE),
                                                    make_pair(CD2,DIHEDRAL)))),
                              make_pair(
                                   CD2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD1,DIHEDRAL)),//CD1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD2,ANGLE),
                                                    make_pair(CD2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         MET,
                         make_vector(
                              make_pair(
                                   SD, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(SD,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(SD,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         ASN,
                         make_vector(
                              make_pair(
                                   OD1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(OD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(OD1,ANGLE),
                                                    make_pair(ND2,DIHEDRAL)))),
                              make_pair(
                                   ND2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(OD1,DIHEDRAL)),//OD1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(ND2,ANGLE),
                                                    make_pair(ND2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         PRO,
                         make_vector(
                              make_pair(
                                   CD, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         GLN,
                         make_vector(
                              make_pair(
                                   OE1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(OE1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(OE1,ANGLE),
                                                    make_pair(NE2,DIHEDRAL)))),
                              make_pair(
                                   NE2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(OE1,DIHEDRAL)),//OE1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(NE2,ANGLE),
                                                    make_pair(NE2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         ARG,
                         make_vector(
                              make_pair(
                                   NE,  
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(NE,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(NE,ANGLE)))),
                              make_pair(
                                   NH1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(NE,DIHEDRAL),
                                                    make_pair(CZ,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(NE,ANGLE),
                                                    make_pair(CZ,ANGLE),
                                                    make_pair(NH1,ANGLE),
                                                    make_pair(NH1,DIHEDRAL)))),
                              make_pair(
                                   NH2, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD,DIHEDRAL),
                                                    make_pair(NE,DIHEDRAL),
                                                    make_pair(CZ,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD,ANGLE),
                                                    make_pair(NE,ANGLE),
                                                    make_pair(CZ,ANGLE),
                                                    make_pair(NH2,ANGLE),
                                                    make_pair(NH2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         SER,
                         make_vector(
                              make_pair(
                                   OG,
                                   make_pair(
                                        make_vector(make_pair(OG,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(OG,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         THR,
                         make_vector(
                              make_pair(
                                   OG1, 
                                   make_pair(
                                        make_vector(make_pair(OG1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(OG1,ANGLE),
                                                    make_pair(CG2,DIHEDRAL)))),
                              make_pair(
                                   CG2, 
                                   make_pair(
                                        make_vector(make_pair(OG1,DIHEDRAL)),//OG1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG2,ANGLE),
                                                    make_pair(CG2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         VAL,
                         make_vector(
                              make_pair(
                                   CG1, 
                                   make_pair(
                                        make_vector(make_pair(CG1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG1,ANGLE),
                                                    make_pair(CG2,DIHEDRAL)))),
                              make_pair(
                                   CG2, 
                                   make_pair(
                                        make_vector(make_pair(CG1,DIHEDRAL)),//CG1 is the chi angle
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG2,ANGLE),
                                                    make_pair(CG2,DIHEDRAL)))))));
               constrained_atoms.insert(
                    make_pair(
                         TRP,
                         make_vector(
                              make_pair(
                                   NE1, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD1,ANGLE)))))));
                              // make_pair(
                              //      CZ3, 
                              //      make_pair(
                              //           make_vector(make_pair(CG,DIHEDRAL),
                              //                       make_pair(CD1,DIHEDRAL)),
                              //           make_vector(make_pair(CB,ANGLE),
                              //                       make_pair(CB,DIHEDRAL),
                              //                       make_pair(CG,ANGLE)))))));
               constrained_atoms.insert(
                    make_pair(
                         TYR,
                         make_vector(
                              make_pair(
                                   OH, 
                                   make_pair(
                                        make_vector(make_pair(CG,DIHEDRAL),
                                                    make_pair(CD1,DIHEDRAL)),
                                        make_vector(make_pair(CB,ANGLE),
                                                    make_pair(CB,DIHEDRAL),
                                                    make_pair(CG,ANGLE),
                                                    make_pair(CD1,ANGLE),
                                                    make_pair(OH,ANGLE)))))));

          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "mode:" << settings.mode << "\n";
               o << "sigma-major-dofs:" << settings.sigma_major_dofs << "\n";
               o << "sigma-minor-dofs:" << settings.sigma_minor_dofs << "\n";
               o << "sample-major-dofs:" << settings.sample_major_dofs << "\n";
               o << "sample-minor-dofs:" << settings.sample_minor_dofs << "\n";
               o << "lagrange-multiplier:" << settings.lagrange_multiplier << "\n";
               o << "skip-proline:" << settings.skip_proline << "\n";
               o << static_cast<const typename Move<CHAIN_TYPE>::Settings &>(settings);
               return o;
          }                    
     } settings;    //!< Local settings object

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveSidechainLocal(CHAIN_TYPE *chain, const Settings &settings = Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "sc-local", settings, random_number_engine), 
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveSidechainLocal(const MoveSidechainLocal &other)
          : MoveCommon(other),
            random_generator_gaussian(*other.random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(other.settings) {

     }

     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveSidechainLocal(const MoveSidechainLocal &other, 
                 RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            random_generator_gaussian(*random_number_engine, 
                                      boost::normal_distribution<>()),
            settings(other.settings) {
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
               this->random_move_range(&start_index, &end_index);
          }

          // Check whether there are any dofs in region
          bool dofs_found = false;
          for (int i=start_index; i<end_index; i++) {
               typename CHAIN_TYPE::Residue &res = (*this->chain)[i];
               if (!(settings.skip_proline && res.residue_type == definitions::PRO)) {
                    dofs_found = true;
                    break;
               }
          }

	  // Remember which angles are modified
          this->move_info = (new MoveInfo(SIDECHAIN))->add_info(std::make_pair(start_index, end_index),
                                                                std::make_pair(start_index, end_index));

          // If there are no degrees of freedom here, mark the move as failed
          if (!dofs_found) {
               return (this->move_info->success = false);
          }
          
          // Clear backup
          major_dof_value_backup.clear();
          minor_dof_value_backup.clear();

          // Clear matrices
          Lt.clear();
          dphi.clear();
          precision_matrix.clear();
          constrained_atom_index.clear();

          // Iterate over residue positions
          for (int i=start_index; i<end_index; i++) {

               // Get residue
               typename CHAIN_TYPE::Residue *res = &(*this->chain)[i];

               // Save old values
               major_dof_value_backup.push_back(res->get_sidechain_dof_values(SIDECHAIN_ATOMS));
               minor_dof_value_backup.push_back(res->get_minor_dof_values());

               // Copy old values to new values
               std::vector<double> major_dof_values(major_dof_value_backup.back());
               std::vector<double> minor_dof_values(minor_dof_value_backup.back());

               if (settings.mode == SCALE_SIGMA || settings.mode == SELECT_DOFS) {

                    // Get maximum chain-distance from CA to any chi angle
                    int min_distance_to_CA = std::numeric_limits<int>::max();
                    int max_distance_to_CA = 0;

                    // Calculate minimum and maximum distance along chain for all dofs
                    if (settings.sample_major_dofs)
                         get_min_max_chain_distance(res, res->chi_atoms, min_distance_to_CA, max_distance_to_CA);
                    if (settings.sample_minor_dofs)
                         get_min_max_chain_distance(res, res->minor_dof_atoms, min_distance_to_CA, max_distance_to_CA);

                    // Draw random number once (used by SELECT_DOFS mode
                    double random_number = this->random_generator_uniform_01();

                    if (settings.sample_major_dofs)
                         sample_unconstrained_endpoint(res, min_distance_to_CA, max_distance_to_CA,
                                                       res->chi_atoms, settings.sigma_major_dofs, random_number, major_dof_values);
                    
                    if (settings.sample_minor_dofs)
                         sample_unconstrained_endpoint(res, min_distance_to_CA, max_distance_to_CA,
                                                       res->minor_dof_atoms, settings.sigma_minor_dofs, random_number, minor_dof_values);

               } else {

                    if (settings.mode == CONSTRAIN_ALL_ENDPOINTS || settings.mode == CONSTRAIN_ONE_ENDPOINT) {

                         int size = 0;
                         if (settings.sample_major_dofs)
                              size += res->chi_atoms.size();
                         if (settings.sample_minor_dofs)
                              size += res->minor_dof_atoms.size();

                         // Create precision matrix
                         precision_matrix.push_back(Matrix(size,size,0.0));
                         if (construct_precision_matrix(res, precision_matrix.back())) {

                              // Calculate "square root" of matrix W_tilde
                              Lt.push_back(precision_matrix.back().cholesky('U'));
                              
                              // Generate 0,1 distributed values
                              Matrix dchi(size, 1);
                              for (int i=0; i<size; i++) {
                                   dchi(i,0) = random_generator_gaussian();
                              }

                              // Determine changes to angular degrees of freedom
                              dphi.push_back(solve_ax_b_triangular(Lt.back(), dchi, 'U'));

                              // Update dof value vectors
                              int offset = 0;
                              if (settings.sample_major_dofs) {
                                   for (unsigned int j=0; j<res->chi_atoms.size(); ++j) {
                                        major_dof_values[j] = fmod(major_dof_values[j] + 3*M_PI + dphi.back()(j+offset,0), 2*M_PI) - M_PI;
                                   }
                                   offset += res->chi_atoms.size();
                              }
                              if (settings.sample_minor_dofs) {
                                   for (unsigned int j=0; j<res->minor_dof_atoms.size(); ++j) {
                                        minor_dof_values[j] = fmod(minor_dof_values[j] + 3*M_PI + dphi.back()(j+offset,0), 2*M_PI) - M_PI;

                                        // Take absolute value if current dof is an angle
                                        if (res->minor_dof_atoms[j].second == ANGLE) {
                                             minor_dof_values[j] = std::fabs(minor_dof_values[j]);
                                        }
                                   }
                              }


                         // If no constraints were found for this residue type, simply insert empty matrices.
                         // This is done to preserve the indexing (for when we calculate the bias)
                         } else {
                              precision_matrix.push_back(Matrix(0,0));
                              Lt.push_back(Matrix(0,0));
                              dphi.push_back(Matrix(0,0));
                         }
                    }

                    // Execute a basic move if no constraints were found - or if settings.mode is none of the above
                    if (precision_matrix.empty() || precision_matrix.back().n_row == 0) {
                         if (settings.sample_major_dofs) {
                              for (unsigned int j=0; j<major_dof_values.size(); ++j) {
                                   double delta = random_generator_gaussian() * settings.sigma_major_dofs;
                                   major_dof_values[j] = fmod(major_dof_values[j] + 3*M_PI + delta, 2*M_PI) - M_PI;
                              }
                         }
                         if (settings.sample_minor_dofs) {
                              for (unsigned int j=0; j<res->minor_dof_atoms.size(); ++j) {
                                   double delta = random_generator_gaussian() * settings.sigma_minor_dofs;
                                   minor_dof_values[j] = fmod(minor_dof_values[j] + 3*M_PI + delta, 2*M_PI) - M_PI;

                                   // Take absolute value if current dof is an angle
                                   if (res->minor_dof_atoms[j].second == ANGLE) {
                                        minor_dof_values[j] = std::fabs(minor_dof_values[j]);
                                   }
                              }
                         }
                    }
               }
               
               // Set values in residue
               if (!major_dof_values.empty())
                    res->set_sidechain_dof_values(major_dof_values, SIDECHAIN_ATOMS);
               if (!minor_dof_values.empty())
                    res->set_minor_dof_values(minor_dof_values);

               // Update positions
               (*this->chain)[i].update_positions_non_backbone();
          }

          return this->move_info->success;
     }

     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();
          
          // Clear saved angles
          major_dof_value_backup.clear();
          minor_dof_value_backup.clear();
     }


     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();

          // If modified_angles is empty, do nothing
          if (!this->move_info->success) {
               return;
          }
                    
          // Restore saved chi angles
          for (int i=this->move_info->modified_angles_start; i<this->move_info->modified_angles_end; i++) {

               // Get residue
               typename CHAIN_TYPE::Residue *res = &(*this->chain)[i];

               if (settings.sample_major_dofs) {
                    // Set values in residue
                    res->set_sidechain_dof_values(major_dof_value_backup[i-this->move_info->modified_angles_start],
                                               definitions::SIDECHAIN_ATOMS);
               }
               if (settings.sample_minor_dofs) {
                    // Set values in residue
                    res->set_minor_dof_values(minor_dof_value_backup[i-this->move_info->modified_angles_start]);
               }

               // Update positions
               res->update_positions_non_backbone();
          }
     }


     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {

          // Only the constrain end point setting introduces a constraint
          if (settings.mode == CONSTRAIN_ALL_ENDPOINTS || settings.mode == CONSTRAIN_ONE_ENDPOINT) {

               double log_bias_forward = calculate_log_bias_forward();
               double log_bias_reverse = calculate_log_bias_reverse();

               return (log_bias_reverse - log_bias_forward);
               
          } else {
               return 0.0;
          }
     }
};

}

#endif
