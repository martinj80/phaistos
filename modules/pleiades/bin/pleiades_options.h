// pleiades_options.h --- commandline option parser for pleiades
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

#include <string>

#ifndef PLEIADES_OPTIONS_H
#define PLEIADES_OPTIONS_H

//! Simple commandline argument option parser for the Pleiades
//! GIT vector clustering.
class PleiadesOptions {
private:

     // interal use
     std::string command_line;     

public:

     //! path to the git file
     std::string git_file;

     //! verbose output switch
     bool verbose;

     //! debug output switch
     bool debug;

     //! try to guestimate a good number of k
     bool guess_k;

     //! Whether to perform k-means clustering
     bool k_means;

     //! Whether to use smart seeding (k-means++) to find initial
     //! cluster center
     bool smart_seed;

     //! Whether to perform weighted k-means clustering
     bool w_k_means;

     //! beta value used in combination with muninn
     //! weighted clustering
     double beta;

     //! Whether to write all members of a cluster into the outputfile.
     //! If FALSE, only the five members closest to the cluster
     //! centroid will be written.
     bool long_output;

     //! Defines the number of iterations to run the clustering for.
     int iterations;

     //! Specifies the number of cluster centers K
     int k;
     
     //! Evaluate clusters weighted after k-means clustering
     bool evaluate_clusters_weighted;

     //! Specifies the distance threshold use for the guestimation of K
     double guess_k_threshold;

#ifdef HAVE_MUNINNLIB

     //! The full path to the muninn log file used in the simulation.
     //! allows to reweight the different vector ie at different temperatures.
     std::string muninn_log;

#endif

     //! Specifies whether to scale the weights to the interval [0,1]
     //! Zero weights are set to the smallest none zero weight found
     //! in the process.
     bool scale_weights;

     //! File to which to write the output
     std::string out_file;

     //! Specifies the full filename of the native pdb structure.
     //! If a native structure is known for the decoys clustered,
     //! Pleiades can directly calculate the according rmsd between
     //! every decoy and the native.
     std::string rmsd_native_pdb;

     //! Specifies the path prefix to decoys clustered. The prefix
     //! plus the GIT vector name need to point to the actual decoy.
     //! If a native structure is known for the decoys clustered,
     //! Pleiades can directly calculate the according rmsd between
     //! every decoy and the native.
     std::string rmsd_decoy_prefix;

     //! Standard constructor
     //!
     //! \param argc commandline argument count
     //! \param argv commandline arguments
     PleiadesOptions(int argc, char *argv[]);

     //! Standard destructor
     ~PleiadesOptions();

     //! Specify the default values for all the
     //! arguments.
     void init();

     //! Parse the command line
     //!
     //! \param argc commandline argument count
     //! \param argv commandline arguments
     void get_command_line(int argc, char *argv[]);

     //! Write the help message to the standard output.
     //! All the class attributes will be reset to their
     //! default settings in the process.
     void print_usage();
};

#endif
