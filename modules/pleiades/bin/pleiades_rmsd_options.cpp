// pleiades_rmsd_options.cpp --- commandline option parser for pleiades_rmsd
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

#include "pleiades_rmsd_options.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>

// Constructor
PleiadesRmsdOptions::PleiadesRmsdOptions(int argc, char *argv[]) {
     get_command_line(argc, argv);
}

// Destructor
PleiadesRmsdOptions::~PleiadesRmsdOptions() {} ;


// Default option values
void PleiadesRmsdOptions::init() {

     decoy_dir = "";
     pdb_list_file = "";

     native_file = "";

     verbose = false;
     debug	= false;

     guess_k 	= false;
     k_means		= true;
     w_k_means = false;
     scale_weights = false;
     beta = 1.0;

#ifdef HAVE_MUNINNLIB
     muninn_log = "";
#endif

     long_output = false;

     out_file = "out.cluster";

     iterations	= 100;
     k = 10;
     guess_k_threshold = 5;
     evaluate_clusters_weighted = false;

     residue_start = 0;
     residue_end = -1;
     
     dump_median = false;
     dump_mean = false;

     command_line = "";
     
}

// Parse the command line
void PleiadesRmsdOptions::get_command_line(int argc, char *argv[]) {
     init();

     for (int i=0; i<argc; i++) {
          std::string arg = argv[i];
          std::string arg_name = "";

          // keyword or argument ?
          if (arg.substr(0,2) == std::string("--")) {
               arg_name = arg.substr(2,arg.length()-2);
          } else if (arg.substr(0,1) == std::string("-")) {
               arg_name = arg.substr(1,arg.length()-1);
          } else if ( i!=0 ){
               // we want to actually store the progname
               continue;
          }

          if (i == 0) {
               command_line = arg;
          } else if (arg_name == std::string("help") || arg_name == std::string("h")) {
               print_usage();
               exit(0);
          } else if (arg_name == std::string("decoy-dir")) {
               i++;
               decoy_dir = argv[i];
          } else if (arg_name == std::string("pdb-list-file")) {
               i++;
               pdb_list_file = argv[i];
          } else if (arg_name == std::string("native-file")) {
               i++;
               native_file = argv[i];
          } else if (arg_name == std::string("o") || arg_name == std::string("out-file")) {
               i++;
               out_file = argv[i];
          } else if (arg_name == std::string("verbose") || arg_name == std::string("v")) {
               verbose = true;
          } else if (arg_name == std::string("debug")) {
               debug = true;
          } else if (arg_name == std::string("guess-k")) {
               guess_k = true;
          } else if (arg_name == std::string("k-means")) {
               k_means = true;
          } else if (arg_name == std::string("w-k-means")) {
               w_k_means = true;
               k_means = false;
          } else if (arg_name == std::string("scale-weights")) {
               scale_weights = true;
          } else if (arg_name == std::string("beta")) {
               i++;
               beta = atof(argv[i]);
#ifdef HAVE_MUNINNLIB
          } else if (arg_name == std::string("muninn-log")) {
               i++;
               muninn_log = argv[i];
#endif
          } else if (arg_name == std::string("long-output")) {
               long_output = true;
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
          } else if (arg_name == std::string("residue-start")) {
               i++;
               residue_start = atoi(argv[i]);
          } else if (arg_name == std::string("residue-end")) {
               i++;
               residue_end = atoi(argv[i]);
          } else if (arg_name == std::string("dump-median")) {
               dump_median = true;
          } else if (arg_name == std::string("dump-mean")) {
               dump_mean = true;
          } else {
               std::cerr << "Unknown commandline option " << arg << " will be ignored! \n";
          }
     }
}


// Write the help message to the standard output.
// All the class attributes will be reset to their
// default settings in the process.
void PleiadesRmsdOptions::print_usage() {
     // reset those values
     init();
     // generate usage informations
     std::cout << "Usage " << command_line << " [options] \n";
     std::cout << "\toptions: \n";
     std::cout << "\t--decoy-dir\t\tdecoy directory [either decoy-dir or pdb-list-file is required]\n";
     std::cout << "\t--pdb-list-file\t\tfile containing decoy filenames [either decoy-dir or pdb-list-file is required]\n";
     std::cout << "\t--out-file\t\twhere to write the output to [default="<< out_file << "] \n";
     std::cout << "\n";
     std::cout << "\t--native-file\t\tspecify the native structure [default="<< native_file << "] \n";
     std::cout << "\n";
     std::cout << "\t--verbose -v\t\tbe talky? [default=" << verbose << "] \n";
     std::cout << "\t--debug\t\t\tprint debugging information [default=" << debug << "] \n";
     std::cout << "\n";
     std::cout << "\t--guess-k\t\ttry to estimate k [default=" << guess_k << "] \n";
     std::cout << "\t--k-means\t\trun k-means clustering [default=" << k_means << "] \n";
     std::cout << "\t--w-k-means\t\trun weighted k-means clustering [default=" << w_k_means << "] \n";
     std::cout << "\t--beta\t\t\tinverse temperature for weighted clustering [default=" << beta << "] \n";
     std::cout << "\t--scale-weights\t\tScale the weights to be in the interval [0,1[ [default=" << scale_weights << "] \n";
     std::cout << "\t--long-output\t\tprint full cluster lists [default=" << long_output << "] \n";
#ifdef HAVE_MUNINNLIB
     std::cout << "\t--muninn-log\t\tread in the weights from Muninn [default=" << muninn_log << "] \n";
#endif
     std::cout << "\n";
     std::cout << "\t--iterations\t\thow many rounds of k-means? [default=" << iterations << "] \n";
     std::cout << "\t--k\t\t\thow many clusters to use in k-means [default=" << k << "] \n";
     std::cout << "\t--guess-k-threshold\tthreshold to use when guessing k [default=" << 	guess_k_threshold << "] \n";
     std::cout << "\t--evaluate-clusters-weighted\tEvaluate clusters weighted after k-means clustering [default=" << evaluate_clusters_weighted << "] \n";
     std::cout << "\t--residue-start\t\twhich start residue used in the truncated rmsd calculation [default=" << residue_start << "] \n";
     std::cout << "\t--residue-end\t\twhich end residue used in the truncated rmsd calculation [default=" << residue_end << "] \n";
     std::cout << "\t--dump-median\t\twhether to dump the median structures [default=" << dump_median << "] \n";
     std::cout << "\t--dump-mean\t\twhether to dump the mean structures [default=" << dump_mean << "] \n";
     std::cout << "\n";
}


// Print all the variables and what value is actually
// chosen. Allows to easily put print all runtime
// variable into the log file.
void PleiadesRmsdOptions::print_options() {
     // generate ouptut
     std::cout << "################################################################################################\n";
     std::cout << "# Usage " << command_line << " [options] \n";
     std::cout << "#\toptions: \n";
     std::cout << "#\tdecoy-dir\t\t=  "<< decoy_dir << " \n";
     std::cout << "#\tpdb-list-file\t\t=  "<< pdb_list_file << " \n";
     std::cout << "#\tout-file\t\t= "<< out_file << " \n";
     std::cout << "#\n";
     std::cout << "#\tnative-file\t\t= "<< native_file << " \n";
     std::cout << "#\n";
     std::cout << "#\tverbose -v\t\t= " << verbose << " \n";
     std::cout << "#\tdebug\t\t\t= " << debug << " \n";
     std::cout << "#\n";
     std::cout << "#\tguess-k\t\t= " << guess_k << " \n";
     std::cout << "#\tk-means\t\t= " << k_means << " \n";
     std::cout << "#\tw-k-means\t\t= " << w_k_means << " \n";
     std::cout << "#\tscale-weights\t\t= " << scale_weights << " \n";
     std::cout << "#\tbeta\t\t= " << beta << " \n";
#ifdef HAVE_MUNINNLIB
     std::cout << "#\tmuninn-log\t\t= " << muninn_log << " \n";
#endif
     std::cout << "#\tlong-output\t\t= " << long_output << " \n";
     std::cout << "#\n";
     std::cout << "#\titerations\t\t= " << iterations << " \n";
     std::cout << "#\tk\t\t\t= " << k << " \n";
     std::cout << "#\tguess-k-threshold\t= " << 	guess_k_threshold << " \n";
     std::cout << "#\tevaluate-clusters-weighted\t= " << evaluate_clusters_weighted << " \n";
     std::cout << "#\tresidue-start\t\t= " << residue_start << " \n";
     std::cout << "#\tresidue-end\t\t= " << residue_end << " \n";
     std::cout << "#\tdump-median\t\t= " << dump_median << " \n";
     std::cout << "#\tdump-mean\t\t= " << dump_mean << " \n";
     std::cout << "#\n";
     std::cout << "################################################################################################\n";
}
