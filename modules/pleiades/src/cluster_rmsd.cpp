// cluster_rmsd.cpp --- cluster class using RMSD distances
// Copyright (C) 2011  Tim Harder, Thomas Hamelryck, Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#include <stdio.h>
#include <string>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <vector>
#include <iostream>	
#include <fstream>
#include <algorithm>
#include <cstring>

#include "energy/term_rmsd.h"

#include "utils/vector_matrix_3d.h"
#include "protein/chain_ca.h"
#include "protein/chain_fb.h"
#include "protein/chain.h"

#include "cluster_rmsd.h"

namespace phaistos {

// For a given chain, extract all CA positions into a vector of 3D-Vectors.
std::vector<Vector_3D> get_ca_positions( ChainCA chain ) {
     std::vector<Vector_3D> return_me;
     AtomIterator<ChainCA,definitions::CA_ONLY,Vector_3D> it(chain);
     for(; !it.end(); ++it) {
          return_me.push_back((*it));
     }
     return return_me;
}


// Constructor
ClusterRmsd::ClusterRmsd(int id, std::vector<ChainCA> *members, std::vector<std::string> *chain_names, std::vector<double> *weights, int residue_start, int residue_end, bool debug) {
     this->n_members = 0;
     this->members = members;
     this->chain_names = chain_names;
     this->weights = weights;
     this->member_ids = new std::vector<int>;
     this->id = id;
     this->debug = debug;
     this->residue_start = residue_start;
     this->residue_end = residue_end;
     this->is_initialized = false;

     if (this->debug) {
          std::cout << "D> created clusterRMSD (id:" << this->id << ")\n";
          std::cout.flush();
     }

}

// Copy constructor
ClusterRmsd::ClusterRmsd(const ClusterRmsd &other) {

     this->member_ids = new std::vector<int>;
     for (int i=0; i<(int)other.member_ids->size(); i++) {
          this->member_ids->push_back((*(other.member_ids))[i]);
     }

     this->debug = other.debug;
     this->members = other.members;
     this->chain_names = other.chain_names;
     this->weights = other.weights;
     this->mean = other.mean;
     this->n_members = other.n_members;
     this->id = other.id;
     this->residue_start = other.residue_start;
     this->residue_end = other.residue_end;
     this->is_initialized = other.is_initialized;

     if (this->debug) {
          std::cout << "D> copied clusterRMSD (id:"<< other.id << " and is now id: " << this->id << ")\n";
          std::cout.flush();
     }
}

// Desctructor
ClusterRmsd::~ClusterRmsd() {

     if (this->debug) {
          std::cout << "D> killed clusterRMSD (id:" << this->id << ")\n";
          std::cout.flush();
     }

     //if (this->member_ids != NULL)
     //     delete this->member_ids;
     this->member_ids = NULL;
     this->members = NULL;
     this->mean = NULL;
     this->chain_names = NULL;
}


// simple getters and setters
void ClusterRmsd::set_id(int new_id) {
     this->id = new_id;
}
int ClusterRmsd::get_id() const {
     return this->id;
}
int ClusterRmsd::size() const {
     return this->n_members;
}


// release all the members
void ClusterRmsd::free_members() {
     this->member_ids->clear();
     this->n_members = 0;
}


// get all the member IDs 
std::vector<int> *ClusterRmsd::get_member_ids() const {
     return this->member_ids;
}

// returns the mean 
ChainCA *ClusterRmsd::get_mean() {
     return this->mean;
}

// Returns the structure with the lowest distance the 
// the mean structure 
ChainCA *ClusterRmsd::get_median(const std::vector<double> &energies) {
     int median = -1;

     double weight_sum = 0.0;
     for (int i = 0; i < this->n_members; i++) {
          weight_sum += get_weight(i);
     }

     // If energies are given, calculate a weighted average and
     // standard deviation
     double energies_mean = 0.0;
     double energies_std_dev = 0.0;
     if (!energies.empty()) {
          for (int i = 0; i < this->n_members; i++) {
               double energy = energies[(*this->member_ids)[i]];
               double weight = get_weight(i);

               energies_mean += weight*energy;
               energies_std_dev += weight * (energy*energy);
          }
          energies_mean /= weight_sum;
          energies_std_dev = sqrt((energies_std_dev/weight_sum) - energies_mean*energies_mean);
     } 


     double min_distance = std::numeric_limits<double>::infinity();
     for (int i = 0; i < this->n_members; i++) {
          double distance = calc_rmsd<definitions::CA_ONLY>( ((*this->members)[(*this->member_ids)[i]]) , *(this->mean),this->residue_start,this->residue_end);

          if (!energies.empty()) {
               double energy = energies[(*this->member_ids)[i]];
               // distance += std::fabs(energy-energies_mean);

               // Skip structures with energies more than 0.5 standard deviation away from mean
               if (std::fabs(energy-energies_mean) > (0.5*energies_std_dev))
                    continue;
          }

          if (distance < min_distance) {
               min_distance = distance;
               median = (*this->member_ids)[i];
          }
     }

     return &(*this->members)[median];
}

// Add a new chain to the cluster
void ClusterRmsd::add_member(int new_member_id, bool recalc_mean) {

     this->member_ids->push_back(new_member_id);
     this->n_members++;

     if (this->debug)
          std::cout << "D> Add member to cluster " << this->id << " with chains ID : " << new_member_id << "\n";

     if (!this->is_initialized && this->n_members == 1) {
          // we are the first
          this->mean = new ChainCA( ((*this->members)[new_member_id]) );
          this->is_initialized = true;
          return;
     }

     if (recalc_mean) {
          //////
          // at this point we need to calculate the
          // new mean or centeroid structure
          (*this->members)[new_member_id].superimpose_onto_chain(*(this->mean));
          std::vector<Vector_3D> pos1;
          std::vector<Vector_3D> pos2;
          pos1 = get_ca_positions(((*this->members)[new_member_id]));
          pos2 = get_ca_positions(*(this->mean));

          int i = 0;

          AtomIterator<ChainCA,definitions::CA_ONLY,Vector_3D> it(*(this->mean));
          for(; !it.end(); ++it) {
               Vector_3D res_vec;
               res_vec = ((pos1[i] - pos2[i]) / (this->n_members)) + pos2[i];
               (*it)=res_vec;
               i++;
          }
     }
}

// Recalculate the mean structure from all members
void ClusterRmsd::calc_mean() {
     if (this->n_members < 1) {
          return;
     }
     if (this->n_members == 1) {
          this->mean = new ChainCA(((*this->members)[(*this->member_ids)[0]]) );
          return;
     }

     double weight_sum = 0.0;
     for (int i = 0; i < this->n_members; i++) {
          weight_sum += get_weight(i);
     }
     for (int i = 0; i < this->n_members; i++) {
          (*this->members)[(*this->member_ids)[i]].superimpose_onto_chain(*(this->mean));
          std::vector<Vector_3D> pos1;
          std::vector<Vector_3D> pos2;
          pos1 = get_ca_positions(((*this->members)[(*this->member_ids)[i]]));
          pos2 = get_ca_positions(*this->mean);

          int j = 0;
          AtomIterator<ChainCA,definitions::CA_ONLY,Vector_3D> it(*(this->mean));
          for(; !it.end(); ++it) {
               Vector_3D res_vec;
               res_vec = (get_weight(i)*(pos1[j] - pos2[j]) / (weight_sum)) + pos2[j];
               (*it)=res_vec;
               j++;
          }
     }

}

// Calculates the standard deviation from of the
// cluster, which is the root mean squared distance
// to the cluster centeroid
double ClusterRmsd::get_std_dev() const {
     double sd = 0.;
     double avg = 0.;
     double d, w;
     double w_sum = 0.;

     if (this->n_members == 0) {
          return 0.;
     }

     for (int i = 0; i < this->n_members; i++) {
          d = calc_rmsd<definitions::CA_ONLY>( ((*this->members)[(*this->member_ids)[i]]), *(this->mean),this->residue_start,this->residue_end);
          w = get_weight(i);
          avg += w * d;
          sd += w * (d * d);
          w_sum += w;
     }
     avg = avg / w_sum;
     sd = sqrt( (sd / w_sum )-(avg*avg) );

     return sd;
}

// Returns the sum of the weights of all members
double ClusterRmsd::get_weight_sum() const {
     double weight_sum = 0.;
     for (int i = 0; i < this->n_members; i++) {
          weight_sum += get_weight(i);
     }
     return weight_sum;
}


// Calculates cluster densitory, defined as the size of the  
// cluster over the standard deviation.
double ClusterRmsd::get_density() const {
     double dens = 0.;

     double sd = this->get_std_dev();


     // the test set contained the same points
     // points multiple times, so we need to
     // exclude too low sd.
     if (this->n_members <= 1 || sd < 0.001) {
          return 0.;
     }

     dens = (1 / sd) * this->get_weight_sum();

     return dens;
}

// Output all members of the cluster to the output stream provided.
void ClusterRmsd::dump_member_list(std::ofstream &o, ChainCA *native, int print_max, bool dump_mean, bool dump_median, const std::vector<double> &energies) {

     if (size() == 0)
          return;

     char filename[200];
     //double dist;
     this->calc_mean();

     sprintf(filename, "chain_%03d_mean.pdb", this->id);
     if (dump_mean) {
        (*this->mean).output_as_pdb_file(filename, (this->mean));
     }

     ChainCA *tmpMedian = this->get_median(energies);
     sprintf(filename, "chain_%03d_median.pdb", this->id);
     if (dump_median) {
         tmpMedian->output_as_pdb_file(filename, (this->mean));
     }

     char tmp[1000];
     if (native != NULL) {
          sprintf(tmp, "# Median Structure : rmsdToMean %10.3f\trmsdToNative %10.3f\t%s", calc_rmsd<definitions::CA_ONLY>(*tmpMedian, *(this->mean),this->residue_start,this->residue_end),calc_rmsd<definitions::CA_ONLY>(*tmpMedian, *(native),this->residue_start,this->residue_end), (*tmpMedian).get_name().c_str());
     } else  {
          sprintf(tmp, "# Median Structure : rmsdToMean %10.3f\t%s", calc_rmsd<definitions::CA_ONLY>(*tmpMedian, *(this->mean),this->residue_start,this->residue_end), (*tmpMedian).get_name().c_str());
     }

     o << tmp << std::endl;

     int len = print_max;
     if (len < 0 || len > this->n_members)
          len = this->n_members;

     for (int i = 0; i < len; i++) {
          if (native != NULL) {
               sprintf( tmp, "median: %10.3f\t native: %10.3f\t %s\t%d\t%d\t%f", calc_rmsd<definitions::CA_ONLY>( (*this->members)[(*this->member_ids)[i]], *(this->mean),this->residue_start,this->residue_end), calc_rmsd<definitions::CA_ONLY>( (*this->members)[(*this->member_ids)[i]], *(native),this->residue_start,this->residue_end), ((*(this->chain_names))[(*this->member_ids)[i]]).c_str(), this->id, this->size(), get_weight(i) );
               } else  {
               sprintf( tmp, "%10.3f\t%s\t%d\t%d\t%f", calc_rmsd<definitions::CA_ONLY>( (*this->members)[(*this->member_ids)[i]], *(this->mean),this->residue_start,this->residue_end), ((*(this->chain_names))[(*this->member_ids)[i]]).c_str(), this->id, this->size(), get_weight(i) );
               }

          o << tmp << std::endl;
     }
}

// Retrieve weight
double ClusterRmsd::get_weight(int index) const {
     if (weights != NULL) {
          return (*weights)[(*this->member_ids)[index]];
     } else {
          return 1.0;
     }
}

// Comparison operator for two ClusterRmsd objects. The comparision is
// made by comparing the density of the two for both contestants: a
// cluster is assumed to be "larger" when the the density is higher.
bool operator>(const ClusterRmsd &a, const ClusterRmsd &b) {
     return (a.get_density() > b.get_density());
}

}
