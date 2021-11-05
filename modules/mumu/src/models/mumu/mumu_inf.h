// mumu_inf.h --- 
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


#ifndef MUMU_TEST14_INF_H
#define MUMU_TEST14_INF_H

#include "mocapy.h"

// MUMU featurevector indices
#define MUMU_FVEC_SIZE 21
#define MUMU_FVEC_H 0
#define MUMU_FVEC_AA 1
#define MUMU_FVEC_N 2
#define MUMU_FVEC_CONTACTS 3
#define MUMU_FVEC_PHI 16
#define MUMU_FVEC_PSI 17
#define MUMU_FVEC_VV0 18
#define MUMU_FVEC_VV1 19
#define MUMU_FVEC_VV2 20

#define MUMU_CONTACTS_SIZE 13

// MUMU mask indices
#define MUMU_MASK_SIZE 6
#define MUMU_MASK_H 0
#define MUMU_MASK_AA 1
#define MUMU_MASK_N 2
#define MUMU_MASK_CONTACTS 3
#define MUMU_MASK_PP 4
#define MUMU_MASK_VV 5

class MumuTest14Inf {

     mocapy::DiscreteDensities *aa_densities;
     mocapy::DiscreteDensities *n_densities;
     mocapy::MultinomialDensities *contacts_densities;
     mocapy::VonMises2dDensities *pp_densities;
     mocapy::GaussianDensities *vv_densities;

public:

     // hidden node parameters
     mocapy::MDArray<double> h_param;
     int h_size;

     MumuTest14Inf(mocapy::DBN *dbn) {
          
          // Get nodes from dbn. First slice should be identical to the rest so just use this.
          std::vector<mocapy::Node*> nodes = dbn->getNodes0();
          
          // Get hidden node parameters and set h_size
          this->h_param = nodes[MUMU_MASK_H]->get_parameters()[0];
          this->h_size = this->h_param.size();
          
          // Get individual nodes and density objects
          // Amino acid type
          mocapy::DiscreteNode* aa_node = dynamic_cast<mocapy::DiscreteNode*>(nodes[MUMU_MASK_AA]);
          this->aa_densities = aa_node->get_densities();

          // Total contacts
          mocapy::DiscreteNode* n_node = dynamic_cast<mocapy::DiscreteNode*>(nodes[MUMU_MASK_N]);
          this->n_densities = n_node->get_densities();

          // Contacts
          mocapy::MultinomialNode* contacts_node = dynamic_cast<mocapy::MultinomialNode*>(nodes[MUMU_MASK_CONTACTS]);
          this->contacts_densities = contacts_node->get_densities();

          // Phi psi
          mocapy::VonMises2dNode* pp_node = dynamic_cast<mocapy::VonMises2dNode*>(nodes[MUMU_MASK_PP]);
          this->pp_densities = pp_node->get_densities();

          // Visible volume
          mocapy::GaussianNode* vv_node = dynamic_cast<mocapy::GaussianNode*>(nodes[MUMU_MASK_VV]);
          this->vv_densities = vv_node->get_densities();
     }

     //! Evaluate the log likelihood of an amino acid given the environment
     void calc_log_lik(std::vector<double> &data_point,
                         double *return_ll_joint,
                         double *return_ll_env,
                         double *return_ll_aa = NULL) {

          assert(data_point.size() == MUMU_FVEC_SIZE);

          // Craft "ptv" vectors (parent,data_1,data_2,...)
          std::vector<double> aa_ptv = mocapy::vec(0.0, data_point[MUMU_FVEC_AA]);
          std::vector<double> n_ptv = mocapy::vec(0.0, data_point[MUMU_FVEC_N]);
          std::vector<double> pp_ptv = mocapy::vec(0.0,
                                                   data_point[MUMU_FVEC_PHI],
                                                   data_point[MUMU_FVEC_PSI]);
          std::vector<double> vv_ptv = mocapy::vec(0.0,
                                               data_point[MUMU_FVEC_VV0],
                                               data_point[MUMU_FVEC_VV1],
                                               data_point[MUMU_FVEC_VV2]);
          // Get HSE MN densities of data vector
          std::vector<double> contacts_ptv(MUMU_CONTACTS_SIZE+1, 0.0);
          for (int i=0; i<MUMU_CONTACTS_SIZE; i++) {
               contacts_ptv[i+1] = data_point[MUMU_FVEC_CONTACTS+i];
          }

          // Check for phi psi
          bool pp_missing = false;
          if (data_point[MUMU_FVEC_PHI] > 90.0 || data_point[MUMU_FVEC_PSI] > 90.0)
               pp_missing = true;

          // Get LL for each node and set maxLL
          std::vector<double> ll_h(h_size,0.0);
          std::vector<double> ll_aa(h_size,0.0);
          std::vector<double> ll_env(h_size,0.0);

          const double lower_ll_lim = -1000; //discrete node zero is 1e-50
          long double log_lik;
          
          for (int h=0; h<h_size; h++) {          

               // log weight of component h
               ll_h[h] = log( h_param[h] );
               
               // Amino acid log probability component of h
               aa_ptv[0] = h+0.001;
               ll_aa[h] = aa_densities->get_lik(aa_ptv, true);
               
               // Sum environment log component of h
               n_ptv[0] = h+0.001;
               ll_env[h] = n_densities->get_lik(n_ptv, true);

               if (data_point[MUMU_FVEC_N] > 0.01) {
                    contacts_ptv[0] = h+0.001;
                    ll_env[h] += contacts_densities->get_lik(contacts_ptv, true);
               }

               if (! pp_missing) {
                    pp_ptv[0] = h+0.001;
                    log_lik = pp_densities->get_lik(pp_ptv, true);
                    if (log_lik < lower_ll_lim) {
                         // std::cout<<"# PP hit lower limit: "<<log_lik<<std::endl;
                         log_lik = lower_ll_lim;
                    }
                    ll_env[h] += log_lik;
               }

               vv_ptv[0] = h+0.001;
               log_lik = vv_densities->get_lik(vv_ptv, true);
               if (log_lik < lower_ll_lim) {
                    // std::cout<<"# VV hit lower limit: "<<log_lik<<std::endl;
                    log_lik = lower_ll_lim;
               }
               ll_env[h] += log_lik;
               
          }
     
          long double LL_joint=0.0; // P(AA,N,C,VV,PP)
          long double LL_env=0.0; // P(N,C,VV,PP)
          long double LL_aa=0.0; // P(AA)

          for (int h=0; h<h_size; h++) {
               if (std::fabs(ll_aa[h] + ll_env[h] + ll_h[h]) > 3000) {
                    std::cout<<"# MUMU ERROR: log lik for h="<<h<<" is to large for sum in non-log space: "<<
                         ll_aa[h] + ll_env[h] + ll_h[h]<<std::endl;
                    // assert(false); //overload in mumu
               }

               LL_joint += expl(ll_aa[h] + ll_env[h] + ll_h[h]);
               LL_env += expl(ll_env[h] + ll_h[h]);
               LL_aa += expl(ll_aa[h] + ll_h[h]);
          }
          
          LL_joint = log(LL_joint);
          LL_env = log(LL_env);
          
          *return_ll_joint = (double) LL_joint;
          *return_ll_env = (double) LL_env;
          if (return_ll_aa) {
               LL_aa = log(LL_aa);
               *return_ll_aa = (double) LL_aa;
          }
     }

     //! Evaluate the log likelihood of a single residue given the environment
     //! log[ L(AA|N,C,VV,PP) ] = LL(AA,N,C,VV,PP)-LL(N,C,VV,PP)
     double calc_aa_ll(std::vector<double> &data_point) {
          double ll_joint, ll_env;
          calc_log_lik(data_point, &ll_joint, &ll_env);
          return ll_joint-ll_env;
     }
          
     //! Evaluate the potential of mean force of a single residue
     //! -log[ L(AA,N,C,VV,PP)/{L(AA)*L(N,C,VV,PP)} ] = LL(AA,N,C,VV,PP)-LL(N,C,VV,PP)
     double calc_pmf(std::vector<double> &data_point) {
          double ll_joint, ll_env, ll_aa;
          calc_log_lik(data_point, &ll_joint, &ll_env, &ll_aa);
          return ll_aa+ll_env-ll_joint;
     }
          
     //! Evaluate the log likelihood of a single residue
     //! LL(AA,N,C,VV,PP)
     double calc_ll(std::vector<double> &data_point) {
          double ll_joint, ll_env;
          calc_log_lik(data_point, &ll_joint, &ll_env);
          return ll_joint;
     }
          
     //! Evaluate the log likelihood of a mocapy Sequence of residues
     //! sum[ LL(AA,N,C,VV,PP) ]
     double calc_ll(mocapy::Sequence &data) {
          int dim = data.get_shape().size();
          assert(dim > 0);

          if (dim == 1) {
               double ll_joint, ll_env;
               calc_log_lik(data.get_values(), &ll_joint, &ll_env);
               return ll_joint;

          } else if (dim == 2) {
               double ll_sum = 0.0;
               double ll_joint, ll_env;
               int n = data.get_shape()[0];
               for (int i=0; i<n; i++) {
                    calc_log_lik(data.get_view(i).get_values(), &ll_joint, &ll_env);
                    ll_sum += ll_joint;
               }
               return ll_sum;
          }
          
          assert(false); // data should have dimension 1 or 2
          return 0.0;
     }
};

#endif
