// cluster_rmsd.h --- cluster class using RMSD distances
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

#include "protein/chain_ca.h"
#include "protein/chain.h"

#ifndef CLUSTER_RMSD_H
#define CLUSTER_RMSD_H

namespace phaistos {

//! Represents a Cluster of Protein Structures. Here the Calpha trace
//! of the protein is used for clustering.
class ClusterRmsd {
private:

     //! Vector of member ids. The stored ids are used as index
     //! for the ChainCA members vector, so that no structure
     //! needs to be duplicated.
     std::vector<int> *member_ids;

     //! Pointer to a global vector of ChainCAs. Members are
     //! accessed via their id to save some memory.
     std::vector<ChainCA> *members;

     //! Pointer to global vector of chain names. The indices should
     //! be the sames as for the members vector.
     std::vector<std::string> *chain_names;

     //! Pointer to global vector of weights. The indices should
     //! be the sames as for the members vector.
     std::vector<double> *weights;

     //! Calpha structure representing the mean.
     //! While this is the mean structure, it may not actually look
     //! like a protein or follow regular protein geometry.
     ChainCA *mean;

     //! number of members in the cluster
     int n_members;

     //! unique id of the cluster
     int id;

     //! Whether the cluster has been initialized properly
     bool is_initialized;
     
     //! Specifies the first residue used in the rmsd calculation
     int residue_start;
     
     //! Specifies the last residue used in the rmsd calculation
     int residue_end;

     //! Debug flag
     bool debug;

     //! Retrieve weight
     double get_weight(int index) const;

public:

     //! Standard constructor, requiring a unique id
     //!
     //! \param id unique cluster id
     //! \param members pointer to the member structure vector
     //! \param chain_names pointer to the chainname vector
     //! \param weights Optional vector of weights
     //! \param residue_start Optional residue range specification: start index
     //! \param residue_end Optional residue range specification: end index
     //! \param debug debug flag
     ClusterRmsd(int id, std::vector<ChainCA> *members, std::vector<std::string> *chain_names, std::vector<double> *weights=NULL, int residue_start = 0, int residue_end = -1, bool debug = false);

     //! Standard copy constructor
     //!
     //! \param other other cluster to clone
     ClusterRmsd(const ClusterRmsd &other);

     //! Standard destructor
     ~ClusterRmsd();

     //! Get the mean vector of this cluster
     ChainCA *get_mean();

     //! Get the Median vector of the cluster,
     //! which is the member that is closest to the theoretical
     //! mean of all members.
     ChainCA *get_median(const std::vector<double> &energies=std::vector<double>());

     //! Returns the sum of the weights of all Members
     double get_weight_sum() const;

     //! Returns the Density of the cluster, which is defined
     //! as the size of the cluster over the standard deviation
     double get_density() const;

     //! Calculates the standard deviation from of the
     //! cluster, which is the root mean squared distance
     //! to the cluster centeroid
     double get_std_dev() const;

     //! Recalculate the mean vector from
     //! all the members of the cluster
     void calc_mean();

     //! Returns a pointer to the member id vector
     std::vector<int> *get_member_ids() const;

     //! Returns the unique cluster id
     int get_id() const;

     //! Returns the size of the cluster (the number of members)
     int size() const;

     //! Set a new id
     //! \param new_id new and unique id
     void set_id(int new_id);

     //! Output all members of the cluster to the output stream provided.
     //!
     //! \param o output stream to write to
     //! \param native Native chain object
     //! \param print_max set a maximum number of members to write, negative numbers will lead to all members being written.
     //! \param dump_mean whether to dump calculated means
     //! \param dump_median whether to dump calculated medians
     //! \param energies Optionally specify energy observables
     void dump_member_list(std::ofstream &o, ChainCA *native, int print_max = -1, bool dump_mean = false, bool dump_median = false, const std::vector<double> &energies=std::vector<double>());

     //! Reset the cluster, empty all member arrays and reset the count.
     void free_members();

     //! Add a new member to the cluster. Per default a new mean will be calculated right away.
     //!
     //! \param new_member_id new GitElement to add
     //! \param recalc_mean whether to incrementally update the mean
     void add_member(int new_member_id, bool recalc_mean = true);

     //! Comparison operator for two ClusterRmsd objects. The comparision
     //! is made by comparing the density of the two contestants: a
     //! cluster is assumed to be "larger" when the the density is
     //! higher.
     //!
     //! \param a first Cluster
     //! \param b second Cluster
     friend bool operator>(const ClusterRmsd &a, const ClusterRmsd &b);

};

}

#endif
