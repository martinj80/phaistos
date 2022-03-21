// pleiades_options.cpp --- commandline option parser for pleiades
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
#include "pleiades_options.h"
#include <iostream>
#include <stdio.h>
#include <cstdlib>

// Constructor
PleiadesOptions::PleiadesOptions(int argc, char *argv[]) {
     get_command_line(argc, argv);
}

// Destructor
PleiadesOptions::~PleiadesOptions() {}


// Default option values
void PleiadesOptions::init() {

     git_file = "";

     verbose = false;
     debug = false;

     guess_k = false;
     k_means = true;
     w_k_means = false;
     long_output = false;
     beta = 1.0;

     smart_seed = false;

     iterations = 100;
     k = 10;
     guess_k_threshold = 30;

#ifdef HAVE_MUNINNLIB
     muninn_log = "";
#endif

     evaluate_clusters_weighted = false;

     scale_weights = false;

     out_file = "out.cluster";

     rmsd_native_pdb = "";
     rmsd_decoy_prefix = ".";

     command_line = "";
}

// Parse the command line
void PleiadesOptions::get_command_line(int argc, char *argv[]) {
     init();

     for (int i = 0; i < argc; i++) {
          std::string arg = argv[i];
          std::string arg_name = "";

          // keyword or argument ?
          if (arg.substr(0, 2) == std::string("--")) {
               arg_name = arg.substr(2, arg.length() - 2);
          } else if (arg.substr(0, 1) == std::string("-")) {
               arg_name = arg.substr(1, arg.length() - 1);
          } else if (i != 0) {
               // we want to actually store the progname
               continue;
          }

          if (i == 0) {
               command_line = arg;
          } else if (arg_name == std::string("help") || arg_name == std::string("h")) {
               print_usage();
               exit(0);
          } else if (arg_name == std::string("git-file")) {
               i++;
               git_file = argv[i];
#ifdef HAVE_MUNINNLIB
          } else if (arg_name == std::string("muninn-log")) {
               i++;
               muninn_log = argv[i];
#endif
          } else if (arg_name == std::string("o") || arg_name == std::string("out-file")) {
               i++;
               out_file = argv[i];
          } else if (arg_name == std::string("rmsd-decoy-prefix")) {
               i++;
               rmsd_decoy_prefix = argv[i];
          } else if (arg_name == std::string("rmsd-native-pdb")) {
               i++;
               rmsd_native_pdb = argv[i];
          } else if (arg_name == std::string("verbose") || arg_name == std::string("v")) {
               verbose = true;
          } else if (arg_name == std::string("debug")) {
               debug = true;
          } else if (arg_name == std::string("guess-k")) {
               guess_k = true;
          } else if (arg_name == std::string("k-means")) {
               k_means = true;
               w_k_means = false;
          } else if (arg_name == std::string("w-k-means")) {
               w_k_means = true;
               k_means = false;
          } else if (arg_name == std::string("beta")) {
               i++;
               beta = atof(argv[i]);
          } else if (arg_name == std::string("smart-seed")) {
               smart_seed = true;
          } else if (arg_name == std::string("long-output")) {
               long_output = true;
          } else if (arg_name == std::string("scale-weights")) {
               scale_weights = true;
          } else if (arg_name == std::string("iterations")) {
               i++;
               iterations = atoi(argv[i]);
          } else if (arg_name == std::string("guess-k-threshold")) {
               i++;
               guess_k_threshold = atof(argv[i]);
          } else if (arg_name == std::string("k")) {
               i++;
               k = atoi(argv[i]);
          } else if (arg_name == std::string("evaluate-clusters-weighted")) {
               evaluate_clusters_weighted = true;
          } else {
               std::cerr << "Unknown commandline option " << arg << " will be ignored! \n";
          }
     }

}


// Write the help message to the standard output.
// All the class attributes will be reset to their
// default settings in the process.
void PleiadesOptions::print_usage() {

     // reset those values
     init();
     // generate usage informations
     std::cout << "Usage " << command_line << " [options] \n";
     std::cout << "\toptions: \n";
     std::cout << "\t--git-file\t\tgit vector file [required] \n";
     std::cout << "\t--out-file\t\tfile to store all the resulting clusters in [required] \n";
     std::cout << "\n";
     std::cout << "\t--verbose\t\tbe talky? [default=" << verbose << "] \n";
     std::cout << "\t--debug\t\t\tprint debugging information [default=" << debug << "] \n";
     std::cout << "\n";
     std::cout << "\t--guess-k\t\ttry to estimate k [default=" << guess_k << "] \n";
     std::cout << "\t--k-means\t\trun k-means clustering [default=" << k_means << "] \n";
     std::cout << "\t--w-k-means\t\trun weighted k-means clustering [default=" << w_k_means << "] \n";
     std::cout << "\t--beta\t\t\tinverse temperature for weighted clustering [default=" << beta << "] \n";
     std::cout << "\t--long-output\t\tprint full cluster lists [default=" << long_output << "] \n";
     std::cout << "\t--smart-seed\t\tsmarter initial seeding [default=" << smart_seed << "] \n";
     std::cout << "\n";
#ifdef HAVE_MUNINNLIB
     std::cout << "\t--muninn-log\t\tread in the weights from Muninn [default=" << muninn_log << "] \n";
#endif
     std::cout << "\t--scale-weights\t\tScale the weights to be in the interval [0,1[ [default=" << scale_weights << "] \n";
     std::cout << "\n";
     std::cout << "\t--rmsd-native-pdb\tNative Structure to calculate the RMSD against [default=" << rmsd_native_pdb << "] \n";
     std::cout << "\t--rmsd-decoy-prefix\tPrefix to find decoy files [default=" << rmsd_decoy_prefix << "] \n";
     std::cout << "\n";
     std::cout << "\t--iterations\t\thow many rounds of k-means? [default=" << iterations << "] \n";
     std::cout << "\t--k\t\t\thow many clusters to use in k-means [default=" << k << "] \n";
     std::cout << "\t--evaluate-clusters-weighted\tEvaluate clusters weighted after k-means clustering [default=" << evaluate_clusters_weighted << "] \n";
     std::cout << "\t--guess-k-threshold\tthreshold to use when guessing k [default=" << guess_k_threshold << "] \n";
     std::cout << "\n";
}

