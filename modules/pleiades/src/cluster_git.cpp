// cluster_git.cpp --- cluster class using git distances
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

#include "cluster_git.h"
#include <vector>

#include "git.h"
#include "git_element.h"

#include "utils/vector_matrix_3d.h"
#include "protein/chain_ca.h"
#include "protein/chain_fb.h"
#include "protein/chain.h"
#include "energy/term_rmsd.h"

#include <stdio.h>
#include <string>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

namespace phaistos {


// Constructor
Cluster::Cluster(int id) {
     this->n_members = 0;
     this->id = id;
     this->members = new std::vector<GitElement>;
     this->is_initialized = false;
     this->mean = NULL;
}


// Destructor
Cluster::~Cluster() {
     return;
     if (this->mean != NULL) {
          delete this->mean;
          this->mean = NULL;
     }

     delete this->members;
     this->members = NULL;
}

// Set a new id 
void Cluster::set_id(int new_id) {
     this->id = new_id;
}

void Cluster::set_weights_of_members(const std::vector<double> &weights) {
     for (int i = 0; i < (int)(*(this->members)).size(); i++) {
          (*(this->members))[i].set_weight(weights.at((*(this->members))[i].get_index()));
     }
}

// Return the unique cluster id
int Cluster::get_id() {
     return this->id;
}

// Returns the size of the cluster (number of members)
int Cluster::get_size() const {
     return this->n_members;
}

// Reset the cluster, empty all member arrays and reset the count.
void Cluster::free_members() {
     (*this->members).clear();
     this->n_members = 0;
}

// Get the mean vector of this cluster
GitElement Cluster::get_mean() {
     return (*this->mean);
}

// Get the Median vector of the cluster,
// which is the member that is closest to the theoretical
// mean of all members.
GitElement Cluster::get_median() {
     GitElement median;
     double min_d = 9999, d;
     if (n_members == 0) {
          return median;
     }

     for (int i = 0; i < this->n_members; i++) {
          d = get_git_distance((*this->mean).get_git(), ((*this->members)[i]).get_git());
          if (d < min_d) {
               min_d = d;
               median = (*this->members)[i];
          }
     }
     return median;
}

// Add a new member to the cluster. Per default a new mean will be calculated right away.
void Cluster::add_member(GitElement new_member, bool recalc_mean) {
     (*this->members).push_back(new_member);
     this->n_members++;

     if (!this->is_initialized && this->n_members == 1) {
          // we are the first
          this->mean = new GitElement();
          (*this->mean) = new_member;
          (*this->mean).set_name("Mean");
          this->is_initialized = true;
     }

     if (recalc_mean) {
          //std::cout << this->id << "       "<< (*this->members).size() << "\n" ;
          std::vector<double> new_mean = (*this->mean).get_git();

          for (unsigned int i = 0; i < new_mean.size(); i++) {
               new_mean[i] = new_mean[i] + ((new_member.get_git(i) - new_mean[i]) / n_members);
          }

          (*this->mean).set_git(new_mean);
     }
}

// Recalculate the mean vector from
// all the members of the cluster
bool Cluster::calc_mean() {
     if (this->n_members < 1) {
          return true;
     }
     if (this->n_members == 1) {
          std::vector<double> tmp = ((*this->members)[0]).get_git();
          (*this->mean).set_git(tmp);
          return true;
     }

     std::vector<double> old_mean = (*this->mean).get_git();
     std::vector<double> new_mean = (*this->mean).get_git();
     double tmp;
     double weight_sum;
     for (unsigned int i = 0; i < new_mean.size(); i++) {
          tmp = 0.;
          weight_sum = 0.;
          for (int j = 0; j < this->n_members; j++) {
               tmp += ((*this->members)[j]).get_git(i)*((*this->members)[j]).get_weight();
               weight_sum +=  ((*this->members)[j]).get_weight();
          }
          new_mean[i] = tmp / weight_sum;
     }
     (*this->mean).set_git(new_mean);

     if (old_mean == new_mean)
          return true;

     return false;
}

// Calculates the standard deviation from of the
// cluster, which is the root mean squared distance
// to the cluster centeroid.
double Cluster::get_std_dev() const {
     double sd = 0.;
     double avg = 0.;
     double d, w;
     double w_sum = 0.;


     if (this->n_members == 0) {
          return 0.;
     }

     for (int i = 0; i < this->n_members; i++) {
          d = get_git_distance(((*this->members)[i]).get_git(), (*this->mean).get_git());
          w = ((*this->members)[i]).get_weight();
          avg += w * d;
          sd += w * (d * d);
          w_sum += w;
     }
     avg = avg / w_sum;
     sd = sqrt( (sd / w_sum )-(avg*avg) );

     return sd;
}

// Returns the sum of the weights of all members
double Cluster::get_weight_sum() const {
     double weight_sum = 0.;
     for (int i = 0; i < this->n_members; i++) {
          weight_sum += ((*this->members)[i]).get_weight() ;
     }
     return weight_sum;
}

// Get the average distance of all members to the
// the cluster mean.
double Cluster::get_avg_dist_mean() const {
     double mean = 0.;

     if (this->n_members == 0) {
          return 0.;
     }

     for (int i = 0; i < this->n_members; i++) {
          mean += ((*this->members)[i]).get_weight() * get_git_distance(((*this->members)[i]).get_git(), (*this->mean).get_git());
     }

     mean = mean / this->get_weight_sum();
     return mean;
}

// Returns the Density of the cluster, which is the defined
// as the size of the cluster over the standard deviation
double Cluster::get_density() const {
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


//! Output all members of the cluster to the output stream provided.
void Cluster::dump_member_list(std::ofstream &o, int max_entries, ChainCA *native, std::string prefix) {

     for (int i = 0; i < this->n_members; i++) {
          ((*this->members)[i]).set_distance_to_mean((*this->mean));
     }
     std::sort((*this->members).begin(), (*this->members).end());

     int len = max_entries;
     if (max_entries < 0) {
          len = (this->n_members);
     }

     double rmsd_native_sum = 0.;
     double rmsd_cluster_sum = 0.;
     double git_cluster_sum = 0.;

     double rmsd_native_sum2 = 0.;
     double rmsd_cluster_sum2 = 0.;
     double git_cluster_sum2 = 0.;

     //
     // lets read in the centroid structure.
     char filename_centroid[1000];
     ChainCA *centroid = NULL;

     sprintf(filename_centroid, "%s/%s", prefix.c_str(), ((this->get_median()).get_name()).c_str());
     if (file_exists(filename_centroid)) {
          centroid = new ChainCA(filename_centroid);
     }
     double min_rmsd=9999999, first_rmsd=0;
     char tmp[1000];

     for (int i = 0; i < (len); i++) {

          double rmsd = 0, rmsd_centroid = 0.;
          double dist = get_git_distance(((*this->members)[i]).get_git(), (*this->mean).get_git());

          if (native != NULL || centroid != NULL) {
               char filename[200];
               sprintf(filename, "%s/%s", prefix.c_str(), (((*this->members)[i]).get_name()).c_str());

               if (file_exists(filename)) {
                    ChainCA *decoy = new ChainCA(filename);
                    if (native != NULL) {
                         rmsd = calc_rmsd<definitions::CA_ONLY> (*native, *decoy);
                    }

                    if (centroid != NULL) {
                         rmsd_centroid = calc_rmsd<definitions::CA_ONLY> (*centroid, *decoy);
                    }
                    delete decoy;
               }
          }

          if (i==0) {
               first_rmsd = rmsd;
          }
          if (rmsd<min_rmsd) {
               min_rmsd = rmsd;
          }

          rmsd_native_sum  += rmsd;
          rmsd_cluster_sum += rmsd_centroid;
          git_cluster_sum  += dist;

          rmsd_native_sum2  += (rmsd*rmsd);
          rmsd_cluster_sum2 += (rmsd_centroid*rmsd_centroid);;
          git_cluster_sum2  += (dist*dist);

          // I know this is not pretty, but way easier to format the output :-)
          sprintf(tmp, "%10.3f\t%s\t%4d\t%5d\t%.3f\t%.3f\t%e\t[", dist,
                    (((*this->members)[i]).get_name()).c_str(), this->get_id(), this->get_size(), rmsd,  rmsd_centroid , ((*this->members)[i]).get_weight());
          o << tmp ;

          // ok for the full output, lets also put the vector out
          int vecLen = ((*this->members)[i]).get_git().size();
          for (int j = 0; j < (vecLen); j++) {
               sprintf(tmp, "%.3f,", ((*this->members)[i]).get_git(j));
               o << tmp;
          }
          o << "]" << std::endl;
     }

     if (len) {
          // and lets print some statistics as well
          double avg_native   = (rmsd_native_sum  / len);
          double avg_cluster  = (rmsd_cluster_sum / len);
          double avg_git      = (git_cluster_sum  / len);

          double std_native  = sqrt((rmsd_native_sum2  / len)-(avg_native*avg_native));
          double std_cluster = sqrt((rmsd_cluster_sum2 / len)-(avg_cluster*avg_cluster));
          double std_git     = sqrt((git_cluster_sum2  / len)-(avg_git*avg_git));

          o << "#" << std::endl;
          if (rmsd_native_sum > 0) {
               sprintf(tmp, "# RMSD_native_avg_/_std_/_first_/_min_:\t%.3f\t%.3f\t%.3f\t%.3f", avg_native, std_native, first_rmsd, min_rmsd);
               o << tmp << std::endl;
          }
          if (rmsd_cluster_sum > 0) {
               sprintf(tmp, "# RMSD_cluster_avg_/_std_:\t\t%.3f\t%.3f", avg_cluster, std_cluster);
               o << tmp << std::endl;
          }
          sprintf(tmp, "# GIT_cluster_avg_/_std_:\t\t%.3f\t%.3f", avg_git, std_git) ;
          o << tmp << std::endl;
          o << "#" << std::endl;
     }
}

}
