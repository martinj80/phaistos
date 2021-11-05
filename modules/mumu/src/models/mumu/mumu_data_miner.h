// mumu_data_miner.h --- 
// Copyright (C) 2012-2013 Kristoffer Enøe Johansson
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


#ifndef MUMU_TEST14_DATA_MINER_H
#define MUMU_TEST14_DATA_MINER_H

#include <cmath>
#include <iostream>
#include <cstdio>
#include <vector>

#include "mocapy.h"

#include "models/visible_volume/sphere.h"
#include "models/mumu/chemical_group_definitions.h"
#include "models/mumu/mumu_inf.h"

#include "protein/chain_fb.h"
#include "utils/vector_matrix_3d.h"

/*
 *  Test14 Model on raw9 data:          H 
 *                                /  /  |  \  \
 *                              AA  N  MN  PP  VV
 *
 *   N: Number of contacts in [0,23]. Counts higher than 23 is set to 23.
 *      24D cathegorial.
 *
 *  MN: Non-local structure modeled by a 13 cathegory multinomial:
 *      Chemical group contacts (min_count 2 with 106 sphere points).
 *      Neighboring residues are included and
 *      own backbone is excluded (incl terminal COO and NH3!).
 *        0    1    2   3   4      5    6    7    8    9  10   11  12
 *      NCO  NH3  COO  OH  SH  C2+C3  CH3  ARO  IMI  GUA  SC  CCN  BB
 *
 *  PP: Phi and psi. Torus
 *
 *  VV: 3D gaussian of uppointing window log volume. 
 *
 */

class MumuDataMinerTest14 {

private:

     phaistos::ChainFB *chain;
     int chain_length;
     phaistos::Vector_3D n_term_vec, c_term_vec;
     Sphere *sphere;
     std::vector<phaistos::Vector_3D> bb_nco, contact_center;
     std::vector<std::vector<std::pair<ChemicalGroupEnum, phaistos::Vector_3D> > > sc_groups;
     mocapy::Sequence data_point;

     std::vector<std::set<int> > neighborhood;
     std::vector<bool> neighborhood_init;
          
     std::vector<std::pair<ChemicalGroupEnum, phaistos::Vector_3D> > find_chemical_groups(phaistos::ResidueFB *res) {
          std::vector<std::pair<ChemicalGroupEnum, phaistos::Vector_3D> > ret(0);
          phaistos::definitions::ResidueEnum rt = res->residue_type;
          int n_groups = sc_group_number[rt];
          for (int g=0; g<n_groups; g++) {
               ChemicalGroupEnum chemical_group = sc_group[rt][g];
               phaistos::Vector_3D pos(0.0, 0.0, 0.0);
               int n_atoms = group_atom_number[rt][g];
               for (int a=0; a<n_atoms; a++) {
                    phaistos::definitions::AtomEnum atom_type = group_atom[rt][g][a];
                    if ( ! res->has_atom(atom_type) ) {
                         printf("Missing %s in %s %d",atom_name[atom_type],residue_name[rt],res->index_res_seq);
                         assert(false);
                    }
                    pos += (*res)[atom_type]->position;
               }
               pos /= 1.0*n_atoms;
               ret.push_back( std::pair<ChemicalGroupEnum, phaistos::Vector_3D>(chemical_group, pos) );
          }
          return ret;
     }
     
     //! Calculate the position of CB
     phaistos::Vector_3D position_cb(phaistos::Vector_3D c_pos, phaistos::Vector_3D n_pos, phaistos::Vector_3D ca_pos) {
          double bond_length = 1.53;           //enghHuberConstants.h L71
          double angle = 89.83*M_PI/180.0;    //enghHuberConstants.h L660
          double dihedral = 52.25*M_PI/180.0; //enghHuberConstants.h L1262
          phaistos::Vector_3D D(bond_length * cos(M_PI-angle),
                                bond_length * cos(M_PI-dihedral) * sin(M_PI-angle),
                                bond_length * sin(M_PI-dihedral) * sin(M_PI-angle));
          phaistos::Vector_3D bc = (c_pos-n_pos).normalize();
          phaistos::Vector_3D no = cross_product((c_pos-ca_pos), bc).normalize();
          phaistos::Vector_3D nbc = cross_product(bc,no);
          phaistos::Matrix_3D basis_change(bc, nbc, no, true);
          D = basis_change*D + c_pos;
          D -= (c_pos - ca_pos);
          return D;
     }
     
     phaistos::Vector_3D find_cb(phaistos::ResidueFB *res) {
          if ( res->has_atom(CB) ) {
               return (*res)[CB]->position;
          } else {
               if ( ! res->has_atom(N) ) {
                    printf("Missing %s in %s %d",atom_name[N],residue_name[res->residue_type],res->index_res_seq);          
                    assert(false);
               }
               if ( ! res->has_atom(C) ) {
                    printf("Missing %s in %s %d",atom_name[C],residue_name[res->residue_type],res->index_res_seq);          
                    assert(false);
               }
               if ( ! res->has_atom(CA) ) {
                    printf("Missing %s in %s %d",atom_name[CA],residue_name[res->residue_type],res->index_res_seq);          
                    assert(false);
               }
               phaistos::Vector_3D cb_pos = position_cb((*res)[C]->position, (*res)[N]->position, (*res)[CA]->position);
               return cb_pos;
          }
     }
     
     phaistos::Vector_3D find_bb_nco(phaistos::ResidueFB *res1, phaistos::ResidueFB *res2) {
          if ( ! res1->has_atom(C) ) {
               printf("Missing %s in %s %d",atom_name[C],residue_name[res1->residue_type],res1->index_res_seq);          
               assert(false);
          }
          if ( ! res2->has_atom(N) ) {
               printf("Missing %s in %s %d",atom_name[N],residue_name[res2->residue_type],res2->index_res_seq);          
               assert(false);
          }
          phaistos::Vector_3D pos = ((*res1)[C]->position + (*res2)[N]->position )/2.0;
          return pos;
     }
     
     phaistos::Vector_3D find_n_term(phaistos::ChainFB *chain) {
          phaistos::ResidueFB *res = &(*chain)[0];
          if ( ! res->has_atom(N) ) {
               printf("Missing %s in %s %d",atom_name[N],residue_name[res->residue_type],res->index_res_seq);          
               assert(false);
          }
          return (*res)[N]->position;
     }
     
     phaistos::Vector_3D find_c_term(phaistos::ChainFB *chain) {
          phaistos::ResidueFB *res = &(*chain)[chain->size()-1];
          if (res->has_atom(O) && res->has_atom(OXT)) {
               return ((*res)[O]->position + (*res)[OXT]->position)/2.0;
          } else {
               if ( ! res->has_atom(C) ) {
                    printf("Missing %s in %s %d",atom_name[C],residue_name[res->residue_type],res->index_res_seq);          
                    assert(false);
               }
               return (*res)[C]->position;
          }
     }
     
public:
     
     //! Default dbn file
     std::string default_dbn_filename;
     std::string default_reference_dbn_filename;

     //! Missing masks
     mocapy::MDArray<mocapy::eMISMASK> mismask_joint, mismask_margi;

     //! Index of hidden node in mismask
     int hidden_index;
     
     //! Default Constructor
     MumuDataMinerTest14() {
          sphere = NULL;
     };

     //! Copy Constructor - uesr must init with new chain
     MumuDataMinerTest14(MumuDataMinerTest14 &other) {};

     //! Dectructor
     ~MumuDataMinerTest14() {
          if (sphere) {
               delete sphere;
               sphere = NULL;
          }
     }
     
     //! Initialize
     void init(phaistos::ChainFB *chain) {
          // default dbn file
          this->default_dbn_filename = "mumu.dbn";
          this->default_reference_dbn_filename = "not-trained-yet";

          // init data point
          this->data_point.set_shape(1,MUMU_FVEC_SIZE);
          this->data_point.set(0,MUMU_FVEC_H,0);

          // init missing mask
          mismask_joint.set_shape(1, MUMU_MASK_SIZE);
          mismask_margi.set_shape(1, MUMU_MASK_SIZE);

          // index of hidden node in mismask
          this->hidden_index = MUMU_MASK_H;     
          
          // Flag mismasks
          // LL(AA|WV,PP,CONT) = LL(AA,VV,PP,CONT) - LL(VV,PP,CONT)
          std::vector<int> v = mocapy::vec(MOCAPY_ALL, MOCAPY_ALL);
          mismask_joint.set_wildcard(v, mocapy::MOCAPY_OBSERVED);
          mismask_margi.set_wildcard(v, mocapy::MOCAPY_OBSERVED);
          v = mocapy::vec(MOCAPY_ALL, MUMU_MASK_H);
          mismask_joint.set_wildcard(v, mocapy::MOCAPY_HIDDEN);
          mismask_margi.set_wildcard(v, mocapy::MOCAPY_HIDDEN);
          v = mocapy::vec(MOCAPY_ALL, MUMU_MASK_AA);
          mismask_margi.set_wildcard(v, mocapy::MOCAPY_MISSING);

          // // test without contacts node
          // v = mocapy::vec(MOCAPY_ALL, MUMU_MASK_CONTACTS);
          // mismask_joint.set_wildcard(v, mocapy::MOCAPY_MISSING);
          // mismask_margi.set_wildcard(v, mocapy::MOCAPY_MISSING);

          // save chain
          this->chain = chain;
          this->chain_length = chain->size();

          neighborhood.resize(chain_length);
          neighborhood_init.assign(chain_length, false);
          
          // construct sphere object
          double sphere_radius = 11.0;
          int min_sphere_points = 100;
          sphere = new Sphere(sphere_radius, min_sphere_points);

          // nuke vectors
          bb_nco = std::vector<phaistos::Vector_3D>(0);
          contact_center = std::vector<phaistos::Vector_3D>(0);
          sc_groups = std::vector<std::vector<std::pair<ChemicalGroupEnum, phaistos::Vector_3D> > >(0);

          // find positions and chemical type of everything
          for (int j=0; j<chain_length; j++) {
               phaistos::ResidueFB *res = &(*chain)[j];
               // find chemical group centers and type for side chain
               sc_groups.push_back( find_chemical_groups(res) );
               // use C beta as contact center
               contact_center.push_back( find_cb(res) );
               // N terminal has no peptide bond
               if (j+1 < chain_length)
                    bb_nco.push_back( find_bb_nco(res, &(*chain)[j+1]) );
          }
          
          // Get position of terminal chemical groups
          n_term_vec = find_n_term(chain);
          c_term_vec = find_c_term(chain);
     };

     // Update position of everything relevant to last move
     void update(phaistos::MoveInfo *move_info=NULL) {
          for (int j=0; j<chain_length; j++) {
               phaistos::ResidueFB *res = &(*chain)[j];
               // find chemical group centers and type for side chain
               sc_groups[j] = find_chemical_groups(res);
               // use C beta as contact center
               contact_center[j] = find_cb(res);
               // N terminal has no peptide bond
               if (j+1 < chain_length)
                    bb_nco[j] = find_bb_nco(res, &(*chain)[j+1]);
          }
          
          // Get position of terminal chemical groups
          n_term_vec = find_n_term(chain);
          c_term_vec = find_c_term(chain);
     };

     std::set<int> get_neighborhood(int resi) {
          if (neighborhood_init[resi])
               return neighborhood[resi];
          else
               assert(false); //neighborhood not found for res
          
          return std::set<int>();
     };
     
     // Make a data point and return a pointer to it
     mocapy::Sequence* make_featurevector(int index) {
          double gf = 1.2;
          phaistos::ResidueFB *res = &(*chain)[index];

          neighborhood[index].clear();
          neighborhood_init[index] = true;

          // ==== ADD NEIGHBOURS IN SPHERE MODEL ==== //
          for (int j=0; j<chain_length; j++) {
               int aa_type = (*chain)[j].residue_type;               
               if (j == 0) {
                    // shadow from N terminal amide
                    phaistos::Vector_3D pos = n_term_vec - contact_center[index];
                    sphere->shade(pos.get_array(), chemical_group_radius[NH3]*gf, j, -2, aa_type, NH3);
                    // shadow from backbone carboxamin
                    pos = bb_nco[j] - contact_center[index];
                    sphere->shade(pos.get_array(), chemical_group_radius[BB]*gf, j, -1, aa_type, BB);
               }
               else if (j == chain_length-1) {
                    // shadow from C terminal carboxyl
                    phaistos::Vector_3D pos = c_term_vec - contact_center[index];
                    sphere->shade(pos.get_array(), chemical_group_radius[COO]*gf, j, -3, aa_type, COO);
               } else {
                    // shadow from backbone carboxamin
                    phaistos::Vector_3D pos = bb_nco[j] - contact_center[index];
                    sphere->shade(pos.get_array(), chemical_group_radius[BB]*gf, j, -1, aa_type, BB);
               }
               // no shadow from own side chain
               if (index==j)
                    continue;
               // shadow from side chain chemical groups
               int sc_groups_size = sc_groups[j].size();
               for (int g=0; g<sc_groups_size; g++) {
                    ChemicalGroupEnum group_type = sc_groups[j][g].first;
                    phaistos::Vector_3D pos = sc_groups[j][g].second - contact_center[index];
                    sphere->shade(pos.get_array(), chemical_group_radius[group_type]*gf, j, g, aa_type, group_type);
               }
          }

          this->data_point.set(0, MUMU_FVEC_AA, res->residue_type);

          // ==== CONTACTS MULTINOMIAL ==== //
          std::vector<std::vector<std::vector<int> > > contacts = sphere->get_contacts();

          if (contacts.size() == 0) {
               for (int i=0; i<13; i++) 
                    this->data_point.set(0, i+MUMU_FVEC_CONTACTS, 0);
               this->data_point.set(0, MUMU_FVEC_N, 0);
               this->mismask_margi.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_MISSING);
               this->mismask_joint.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_MISSING);
          } else {
               // make multinomial
               std::vector<int> multinomial(CHEMICAL_GROUP_ENUM_SIZE,0);
               const int c_size = contacts[0].size();
               const int min_counts = 2;
               for (int c=0; c<c_size; c++) {
                    if (contacts[0][c][2] >= min_counts) {
                         // do not include own backbone
                         if (contacts[0][c][0] == index)
                              continue;
                         // convert atom index to group_type
                         ChemicalGroupEnum group_type = NONE;
                         if (contacts[0][c][1] == -1)
                              group_type = BB;
                         else if (contacts[0][c][1] == -2)
                              group_type = NH3;
                         else if (contacts[0][c][1] == -3)
                              group_type = COO;
                         else {
                              ResidueEnum res_type = (*chain)[contacts[0][c][0]].residue_type;
                              group_type = sc_group[res_type][contacts[0][c][1]];
                         }
                         // check if type is still None
                         assert(group_type >= 0);
                         // count in final multinomial
                         multinomial[group_type] += 1;

                         // registre res id in neighbor set
                         neighborhood[index].insert(contacts[0][c][0]);
                    }
               }

               // put multinomail in data point
               int offset = 0;
               int total_counts = 0;
               for (int i=0; i<13; i++) {
                    int counts = multinomial[i+offset];
                    if (i == C2) {
                         counts += multinomial[C3];
                         offset = 1;
                    }
                    total_counts += counts;
                    this->data_point.set(0, i+MUMU_FVEC_CONTACTS, counts);
               }

               if (total_counts > 23)
                    total_counts = 23;

               if (total_counts == 0) {
                    this->mismask_margi.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_MISSING);
                    this->mismask_joint.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_MISSING);
               } else {
                    this->mismask_margi.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_OBSERVED);
                    this->mismask_joint.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_OBSERVED);
               }
               
               this->data_point.set(0, MUMU_FVEC_N, total_counts);
               // this->mismask_margi.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_OBSERVED);
               // this->mismask_joint.set(0, MUMU_MASK_CONTACTS, mocapy::MOCAPY_OBSERVED);
          }

          // ==== VISIBLE WINDOW VOLUME ==== //
          // calculate C alpha - contact center vector
          phaistos::Vector_3D ca_cc_vec3d = contact_center[index] - (*res)[CA]->position;
          double *ca_cc_array = ca_cc_vec3d.get_array();
          std::vector<double> ca_cc_vec(ca_cc_array, ca_cc_array+3);

          // calculate C alpha - C vector
          phaistos::Vector_3D ca_c_vec3d = (*res)[C]->position - (*res)[CA]->position;
          double *ca_c_array = ca_c_vec3d.get_array();
          std::vector<double> ca_c_vec(ca_c_array, ca_c_array+3);

          // Orient sphere
          sphere->assign_tetrahedral_window(ca_cc_vec, ca_c_vec);
          
          // Calculate visible volume seen through each window
          std::vector<double> wv = sphere->calc_visible_volume();

          // Put in data
          this->data_point.set(0, MUMU_FVEC_VV0, log(wv[0]));
          this->data_point.set(0, MUMU_FVEC_VV1, log(wv[1]));
          this->data_point.set(0, MUMU_FVEC_VV2, log(wv[2]));

          // ==== PHI PSI ==== //
          if (index != 0 && index != chain_length-1) {
               this->mismask_margi.set(0, MUMU_MASK_PP, mocapy::MOCAPY_OBSERVED);
               this->mismask_joint.set(0, MUMU_MASK_PP, mocapy::MOCAPY_OBSERVED);
               this->data_point.set(0, MUMU_FVEC_PHI, res->get_phi());
               this->data_point.set(0, MUMU_FVEC_PSI, res->get_psi());
          } else {
               this->mismask_margi.set(0, MUMU_MASK_PP, mocapy::MOCAPY_MISSING);
               this->mismask_joint.set(0, MUMU_MASK_PP, mocapy::MOCAPY_MISSING);
               this->data_point.set(0, MUMU_FVEC_PHI, 99.0);
               this->data_point.set(0, MUMU_FVEC_PSI, 99.0);
          }
          
          // std::cout<<mismask_margi<<std::endl;
          // std::cout<<mismask_joint<<std::endl;
          // std::cout<<contacts<<std::endl;
          // std::cout<<*res<<":\n"<<data_point<<std::endl;
          // assert(false);
          
          sphere->reset();
          
          return &this->data_point;
     };
};

#endif
