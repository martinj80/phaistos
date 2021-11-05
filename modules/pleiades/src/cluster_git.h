// cluster_git.h --- cluster class using git distances
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

#ifndef CLUSTER_GIT_H
#define CLUSTER_GIT_H

#include <vector>
#include <fstream>
#include <iostream>

#include "git.h"
#include "git_element.h"
#include "utils/vector_matrix_3d.h"
#include "protein/chain_ca.h"
#include "protein/chain_fb.h"

namespace phaistos {

//! Represents a cluster of git vectors, providing easy access
//! to routine tasks such as finding the mean or median structure,
//! adding or removing members or printing the results to the standard out.
class Cluster {
private:

     //! Vector of member GIT elements
     std::vector<GitElement> *members;

     //! Pointer to the mean GIT element
     GitElement *mean;

     //! Number of members in the cluster
     int n_members;

     //! Unique cluster id
     int id;

     //! Whether the cluster been initialized properly
     bool is_initialized;

public:

     //! Standard constructor, requiring a unique id
     //!
     //! \param id unique cluster id
     Cluster(int id);

     //! Standard destructor
     ~Cluster();

     //! Get the mean vector of this cluster
     GitElement get_mean();

     //! Get the Median vector of the cluster,
     //! which is the member that is closest to the theoretical
     //! mean of all members.
     GitElement get_median();

     //! Returns the Density of the cluster, which is the defined
     //! as the size of the cluster over the standard deviation
     double get_density() const;

     //! Calculates the standard deviation from of the
     //! cluster, which is the root mean squared distance
     //! to the cluster centeroid.
     double get_std_dev() const;

     //! Recalculate the mean vector from
     //! all the members of the cluster
     bool calc_mean();

     //! Get the average distance of all members to the
     //! the cluster mean.
     double get_avg_dist_mean() const;

     //! Return the unique cluster id
     int get_id();

     //! Returns the size of the cluster (number of members)
     int get_size() const;

     //! Returns the sum of the weights of all members
     double get_weight_sum() const;

     //! Set a new id 
     //! \param new_id new and unique id
     void set_id(int new_id);
     
     //! Sets the member weights
     //! \param weights list of weights for all git vectors
     void set_weights_of_members(const std::vector<double> &weights);

     //! Output all members of the cluster to the output stream provided.
     //!
     //! \param o output stream to write to
     //! \param max_entries set a maximum number of members to write, negative numbers will lead to all members being written.
     //! \param native ChainCA object of the native chain, to calculate the RMSD against.
     //! \param prefix filename prefix necessary to access the decoy structure files (ie PDB files)
     void dump_member_list(std::ofstream &o, int max_entries = -1, ChainCA *native = NULL, std::string prefix = ".");

     //! Reset the cluster, empty all member arrays and reset the count.
     void free_members();

     //! Add a new member to the cluster. Per default a new mean will be calculated right away.
     //!
     //! \param new_member new GitElement to add
     //! \param recalc_mean whether to incrementally update the mean
     void add_member(GitElement new_member, bool recalc_mean = true);

     //! Comparison operator for two Cluster objects. The comparision
     //! is made by comparing the density of the two contestants: a
     //! cluster is assumed to be "larger" when the the density is
     //! higher.
     //!
     //! \param a first Cluster
     //! \param b second Cluster
     friend inline bool operator>(const Cluster &a, const Cluster &b) {
          return (a.get_density() > b.get_density());
     }
};

}
#endif

