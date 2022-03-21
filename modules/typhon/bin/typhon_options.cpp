// typhon_options.cpp --- option parser for typhon.cpp
// Copyright (C) 2011  Tim Harder, Martin Paluszewski, Thomas Hamelryck, Wouter Boomsma
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

#include <cstdlib>
#include <cstring>
#include <limits>
#include <exception>

#include "typhon_options.h"

//! Constructor
TyphonOptions::TyphonOptions(int argc, char *argv[]) {
     char * tmpRoot;
     tmpRoot = getenv("PHAISTOS_ROOT");
     if (tmpRoot != NULL && strcmp(tmpRoot, "") != 0) {
          phaistos_root = std::string(tmpRoot);
     } else {
          phaistos_root="../";
     }
     init();
     get_command_line(argc, argv);
}

//! Destructor
TyphonOptions::~TyphonOptions() {}

//! Initializer
void TyphonOptions::init() {
     //
     output_directory = "samples/";
     dump_interval = 10000;
     reinitialize_interval = 0;
     dump_git = false;
     //
     // Sequence input
     pdb_file = "";
     ss_file = "";
     restore_network = "";
     //
     // Torus parameter
     dbn_parameter_dir = std::string(phaistos_root + std::string("/data/backbone_dbn_parameters/"));
     dbn_parameter_file = "/p_TORUS_H55_BIC-2070230_FULL_SABMARK_NEW_SS.txt";
     ///
     // runtime options
     iterations = 10000000; // default to 10mio
     burnin = 1000; 
     seed = -1;
     //
     //char tmp[2000];
     //sprintf(tmp, "%s/./data/sidechain_dbns/models/basilisk/data/basilisk.dbn", phaistos_root.c_str());
     //sc_dbn_file = std::string(tmp);
     sc_dbn_file = std::string(phaistos_root + std::string("/data/mocapy_dbns/basilisk.dbn"));

     sc_move_weight = 5;
     hbond_weight = 1;
     // local sc move
     use_local_sc_move = false;
     local_sc_move_weight = 5;
     // semi local move
     use_semi_local_move = false;
     semi_local_move_weight = 1;
     // pivot move
     use_pivot_move = false;
     pivot_move_weight = 1;
     //
     // run modi
     mode_mcmc = "metropolis";
     mode_move = "crisp";
     mode_energy = "all";
     mode_analyze = false; // run analysis of the network only ..
     generate_pymol = false ;  // generate Pymol script to visualize the network

     ignore_hbonds = false;
     ignore_sc_hbonds = false;
     ignore_bbsc_hbonds = false;
     ignore_gaussian = false;
     ignore_ssbond = false;
     ignore_fixed = true;
     init_hbond_from_native = false;
     no_prune = false;
     create_phaistos_input = false;
     ignore_ss = false;

     //
     ca_cutoff = 6.;
     ca_skip = 5;
     //
     // dehydron cutoffs
     dehydron_bb_cutoff = 14;
     dehydron_bbsc_cutoff = 9;
     dehydron_sc_cutoff = 7;
     //
     // a couple std switches
     verbose = false;
     debug = false;
}

//! Parse the command line
void TyphonOptions::get_command_line(int argc, char *argv[]) {

     init();
     commandline = "";

     for (int i = 0; i < argc; i++) {
          std::string arg = argv[i];
          std::string argName = "";

          // keyword or argument ?
          if (arg.substr(0, 2) == std::string("--")) {
               argName = arg.substr(2, arg.length() - 2);
          } else if (arg.substr(0, 1) == std::string("-")) {
               argName = arg.substr(1, arg.length() - 1);
          } else if (i != 0) {
               // we want to actually store the progname
               continue;
          }

          //
          // what do we got here then ?
          if (i == 0) {
               commandline = arg;
          } else if (argName == std::string("help") || argName == std::string("h")) {
               print_usage();
               exit(0);
          } else if (argName == std::string("output-directory")) {
               i++;
               output_directory = argv[i];
          } else if (argName == std::string("dump-interval")) {
               i++;
               dump_interval = atoi(argv[i]);
          } else if (argName == std::string("dump-git")) {
               dump_git = true;
          } else if (argName == std::string("pdb-file")) {
               i++;
               pdb_file = argv[i];
          } else if (argName == std::string("restore-network-from-file") || argName == std::string("restore")) {
               i++;
               restore_network = argv[i];
          } else if (argName == std::string("ss-file")) {
               i++;
               ss_file = argv[i];
          } else if (argName == std::string("dbn-parameter-dir")) {
               i++;
               dbn_parameter_dir = argv[i];
          } else if (argName == std::string("dbn-parameter-file")) {
               i++;
               dbn_parameter_file = argv[i];
          } else if (argName == std::string("iterations")) {
               i++;
               iterations = atoi(argv[i]);
          } else if (argName == std::string("reinitialize-interval")) {
               i++;
               reinitialize_interval = atoi(argv[i]);
          } else if (argName == std::string("burnin")) {
               i++;
               burnin = atoi(argv[i]);
          } else if (argName == std::string("seed")) {
               i++;
               seed = atoi(argv[i]);
          } else if (argName == std::string("sc-dbn-file")) {
               i++;
               sc_dbn_file = argv[i];
          } else if (argName == std::string("sc-move-weight")) {
               i++;
               sc_move_weight = atoi(argv[i]);
          } else if (argName == std::string("use-local-sc-move")) {
               i++;
               use_local_sc_move = true;
          } else if (argName == std::string("local-sc-move-weight")) {
               i++;
               local_sc_move_weight = atof(argv[i]);
          } else if (argName == std::string("use-semi-local-move")) {
               i++;
               use_semi_local_move = true;
          } else if (argName == std::string("semi-local-move-weight")) {
               i++;
               semi_local_move_weight = atof(argv[i]);
          } else if (argName == std::string("use-pivot-move")) {
               i++;
               use_pivot_move = true;
          } else if (argName == std::string("pivot-move-weight")) {
               i++;
               pivot_move_weight = atof(argv[i]);
          }  else if (argName == std::string("mode-mcmc")) {
               i++;
               mode_mcmc = argv[i];
          } else if (argName == std::string("mode-move")) {
               i++;
               mode_move = argv[i];
          } else if (argName == std::string("mode-energy")) {
               i++;
               mode_energy = argv[i];
          } else if (argName == std::string("analyze") || argName == std::string("mode-analyze")) {
               mode_analyze = true;
               verbose = true;
          } else if (argName == std::string("pymol") || argName == std::string("generate-pymol")) {
               create_phaistos_input = true;
               generate_pymol = true;
          } else if (argName == std::string("hbond-weight")) {
               i++;
               hbond_weight = atoi(argv[i]);
          }  else if (argName == std::string("ignore-hbonds")) {
               ignore_hbonds = true;
          } else if (argName == std::string("ignore-sc-hbonds")) {
               ignore_sc_hbonds = true;
          } else if (argName == std::string("ignore-bbsc-hbonds")) {
               ignore_bbsc_hbonds = true;
          } else if (argName == std::string("ignore-gaussians")) {
               ignore_gaussian = true;
          } else if (argName == std::string("ignore-ssbonds")) {
               ignore_ssbond = true;
          } else if (argName == std::string("ignore-fixed")) {
               ignore_fixed = true;
          } else if (argName == std::string("enable-fixed")) {
               ignore_fixed = false;
          } else if (argName == std::string("init-hbond-from-native")) {
               init_hbond_from_native = true;
          } else if (argName == std::string("no-prune")) {
               no_prune = true;
          } else if (argName == std::string("create-phaistos-structure")) {
               create_phaistos_input = true;
          } else if (argName == std::string("ignore-ss")) {
               ignore_ss = true;
          } else if (argName == std::string("ca-cutoff")) {
               i++;
               ca_cutoff = atof(argv[i]);
          } else if (argName == std::string("ca-skip")) {
               i++;
               ca_skip = atof(argv[i]);
          } else if (argName == std::string("dehydron-bb-cutoff")) {
               i++;
               dehydron_bb_cutoff = atoi(argv[i]);
          } else if (argName == std::string("dehydron-bbsc-cutoff")) {
               i++;
               dehydron_bbsc_cutoff = atoi(argv[i]);
          } else if (argName == std::string("dehydron-sc-cutoff")) {
               i++;
               dehydron_sc_cutoff = atoi(argv[i]);
          } else if (argName == std::string("verbose")) {
               verbose = true;
          } else if (argName == std::string("debug")) {
               debug = true;
               verbose = true;
          } else {
               std::cerr
                    << "\n###############################################\nERROR: Unknown commandline option "
                    << arg << "! \n\n";
               std::cerr.flush();
               throw UnknownOption;
          }
     }

     std::cout << "\n";
}


//! Display usage information
void TyphonOptions::print_usage() {
     // reset those values .. except the debug one
     bool my_debug = debug;
     init();
     // generate usage informations
     std::cout << "Usage " << commandline << " [options] \n";
     std::cout << "\t\n";
     std::cout << "\t# common options: \n";
     std::cout << "\t--verbose\t\t\tproduce more output (recommended) [default=" << verbose << "] \n";
     std::cout << "\t--debug\t\t\t\tprint debugging information (probably not necessary) [default="	<< debug << "] \n";
     std::cout << "\t--analyze\t\t\tAnalyze and print the hbond network only, no sampling! [default=" << mode_analyze << "]  \n";
     std::cout << "\t--iterations\t\t\tHow many samples to generate in the simulation [default=" << iterations << "]  \n";
     std::cout << "\t--burnin\t\t\tHow many iterations to skip before starting to sample [default=" << burnin << "]  \n";
     std::cout << "\t--reinitialize-interval\t\tHow often to reinitialize to the native structure  [default="	<< reinitialize_interval << "]  \n";
     std::cout << "\t--seed\t\t\t\tSeed the random number generators with [default=\"seed with current time stamp aka random\"]  \n";
     std::cout << "\n";
     std::cout << "\t# input options:\n";
     std::cout << "\t--pdb-file\t\t\tWhich structure to use [required] \n";
     std::cout << "\t--ss-file\t\t\tCondition on different secondary structure [default=" << ss_file << "]  \n";
     std::cout << "\t--restore\t\t\trestore h-bond network from file [default=" << restore_network << "]  \n";
     std::cout << "\n";
     std::cout << "\t# output options:\n";
     std::cout << "\t--output-directory\t\tWhere to put all the generated structures  [default="	<< output_directory << "]  \n";
     std::cout << "\t--dump-interval\t\t\tHow often shall a structure be saved  [default="	<< dump_interval << "]  \n";
     std::cout << "\t--dump-git\t\t\tStore the GIT vectors as well? [default="	<< dump_git << "]  \n";
     std::cout << "\t--create-phaistos-structure\tdump a \"native\" structure with all atoms at the beginning of the run [default=" << create_phaistos_input << "]  \n";
     std::cout << "\t--generate-pymol\t\t\tGenerate a Python script that will visualize the network in PyMOL [default=" << generate_pymol << "]  \n";
     std::cout << "\n";
     std::cout << "\t# network options:\n";
     std::cout << "\t--no-prune\t\t\tdo not prune/optimize the hbond network  [default="	<< no_prune << "]  \n";
     std::cout << "\t--ignore-hbonds\t\t\tdisregard any hbonds found (overwrites the restore network) [default=" << ignore_hbonds << "]  \n";
     std::cout << "\t--ignore-sc-hbonds\t\tdisregard any side chain - side chain hbonds found (overwrites the restore network) [default=" << ignore_sc_hbonds << "]  \n";
     std::cout << "\t--ignore-bbsc-hbonds\t\tdisregard any backbone - side chain hbonds found (overwrites the restore network) [default=" << ignore_bbsc_hbonds << "]  \n";
     std::cout << "\t--ignore-ssbonds\t\tdisregard any ssbond found (overwrites the restore network) [default=" << ignore_ssbond << "]  \n";
     std::cout << "\t--ignore-gaussians\t\tdisregard any gaussian found (overwrites the restore network) [default=" << ignore_gaussian << "]  \n";
     std::cout << "\t--ignore-fixed\t\t\tdisregard any interactions with ligands (overwrites the restore network) [default=" << ignore_fixed << "]  \n";
     std::cout << "\t--enable-fixed\t\t\tinclude interactions with ligands (overwrites the restore network) [default=" << !ignore_fixed << "]  \n";
     std::cout << "\t--init-hbond-from-native\tper default we use parameter estimated from the top 500 set, \n\t\t\t\t\tset this flag to initialize hbond distances to the native [default=" << init_hbond_from_native << "]  \n";
     std::cout << "\t--ca-cutoff\t\t\tcutoff distance in Angstrom to be considered a Calpha contact  [default=" << ca_cutoff << "]  \n";
     std::cout << "\t--ca-skip\t\t\thow many residues to skip along the chain before considering a Calpha contact  [default=" << ca_skip << "]  \n";
     std::cout << "\n";
     std::cout << "\t# path options\n";
     std::cout << "\t--dbn-parameter-dir\t\tWhich torusDBN parameter set to use [default="	<< dbn_parameter_dir << "]  \n";
     std::cout << "\t--dbn-parameter-file\t\tWhich torusDBN parameter set to use [default=" << dbn_parameter_file << "]  \n";
     std::cout << "\t--sc-dbn-file\t\t\tWhich DBN to use for the sidechain moves [default=" << sc_dbn_file << "]  \n";
     std::cout << "\n";
     std::cout << "\t# move options\n";
     std::cout << "\t--sc-move-weight\t\tHow many sidechain moves per CRISP move [default=" << sc_move_weight << "]  \n";
     //std::cout << "\t--use-local-sc-move\t\tUse local sidechain moves, allows for more detailed dynamics [default=" << use_local_sc_move << "]  \n";
     //std::cout << "\t--local-sc-move-weight\t\tHow many sidechain moves per CRISP move [default=" << local_sc_move_weight << "]  \n";
     std::cout << "\n";
     std::cout << "\t# sampling options:\n";
     std::cout << "\t--ignore-ss\t\t\tdo not fix the secondary structure for sampling [default=" << ignore_ss << "]  \n";
     std::cout << "\n";
     std::cout << "\t# dehydron options:\n";
     std::cout << "\t--dehydron-bb-cutoff\t\tconsider bb hbonds below (d<=cutoff) as dehydrated aka weak/broken  [default=" << dehydron_bb_cutoff << "]  \n";
     std::cout << "\t--dehydron-bbsc-cutoff\t\tconsider bb-sc hbonds below (d<=cutoff) as dehydrated aka weak/broken  [default=" << dehydron_bbsc_cutoff << "]  \n";
     std::cout << "\t--dehydron-sc-cutoff\t\tconsider sc hbonds below (d<=cutoff) as dehydrated aka weak/broken  [default=" << dehydron_sc_cutoff << "]  \n";
     std::cout << "\n";
     if (my_debug) {
          std::cout << "\t# advanced options:\n";
          std::cout << "\t--mode-move\t\t\twhich moves to use [default=" << mode_move << "]  \n";
          std::cout << "\t--mode-energy\t\t\twhich energy to use [default=" << mode_energy << "]  \n";
          std::cout << "\t--mode-mcmc\t\t\twhich mcmc strategy to use [default=" << mode_mcmc << "]  \n";
          std::cout << "\t--hbond-weight\t\t\thow to weight the hbond term [default=" << hbond_weight << "]  \n";
          std::cout << "\n";
          std::cout << "\t# advanced move options:\n";
          std::cout << "\t--use-semi-local-move\t\tUse semi local backbone moves, rather radical moves [default=" << use_semi_local_move << "]  \n";
          std::cout << "\t--semi-local-move-weight\tHow many semi local moves per CRISP move [default=" << semi_local_move_weight << "]  \n";
          //std::cout << "\t--use-pivot-move\t\tUse pivot moves, very radical moves [default=" << use_pivot_move << "]  \n";
          //std::cout << "\t--pivot-move-weight\t\tHow many pivot moves per CRISP move [default=" << pivot_move_weight << "]  \n";
          std::cout << "\n";
     }
     std::cout << "\n";
     if (phaistos_root == "../") {
          std::cout << "\tWARNING: PHAISTOS_ROOT environmental variable is not set. You are fine as long you \n";
          std::cout << "\tWARNING: are running from your ./build/bin/ directory. Otherwise please set the  \n";
          std::cout << "\tWARNING: variable to match the fully qualified path to your Phaistos build directory.\n";
     }
     std::cout << "\n";
}

//! Display available options
void TyphonOptions::print_options() {
     std::cout << "##############################################################################################################\n";
     std::cout << "#\n";
     std::cout << "#\tPHAISTOS_ROOT\t\t\t = " << phaistos_root<< " \n";
     std::cout << "#\n";
     std::cout << "#\toptions used for this run: \n";
     std::cout << "#\n";
     std::cout << "#\t";
     std::cout << "#\tcommon options: \n";
     std::cout << "#\n";
     std::cout << "#\t--verbose\t\t\t = " << verbose << "\n";
     std::cout << "#\t--debug\t\t\t\t = "	<< debug << "\n";
     std::cout << "#\t--analyze\t\t\t = " << mode_analyze << " \n";
     std::cout << "#\t--iterations\t\t\t = " << iterations << " \n";
     std::cout << "#\t--burnin\t\t\t = " << burnin << " \n";
     std::cout << "#\t--reinitialize-interval\t\t\t = "	<< reinitialize_interval << " \n";
     std::cout << "#\t--seed\t\t\t\t = " << seed << " \n";
     std::cout << "#\n";
     std::cout << "#\tinput options:\n";
     std::cout << "#\t--pdb-file\t\t\t = " << pdb_file << " \n";
     std::cout << "#\t--ss-file\t\t\t = " << ss_file << " \n";
     std::cout << "#\t--restore\t\t\t = " << restore_network << " \n";
     std::cout << "#\n";
     std::cout << "#\toutput options:\n";
     std::cout << "#\t--output-directory\t\t = "	<< output_directory << " \n";
     std::cout << "#\t--dump-interval\t\t\t = "	<< dump_interval << " \n";
     std::cout << "#\t--dump-git\t\t\t = "	<< dump_git << " \n";
     std::cout << "#\t--create-phaistos-structure\t = " << create_phaistos_input << " \n";
     std::cout << "#\t--generate-pymol\t\t\t = " << generate_pymol << " \n";
     std::cout << "#\n";
     std::cout << "#\tnetwork options:\n";
     std::cout << "#\t--no-prune\t\t\t = "	<< no_prune << " \n";
     std::cout << "#\t--ignore-hbonds\t\t\t = " << ignore_hbonds << " \n";
     std::cout << "#\t--ignore-sc-hbonds\t\t = " << ignore_sc_hbonds << " \n";
     std::cout << "#\t--ignore-bbsc-hbonds\t\t = " << ignore_bbsc_hbonds << " \n";
     std::cout << "#\t--ignore-ssbonds\t\t = " << ignore_ssbond << " \n";
     std::cout << "#\t--ignore-fixed\t\t\t = " << ignore_fixed<< " \n";
     std::cout << "#\t--ignore-gaussians\t\t = " << ignore_gaussian << " \n";
     std::cout << "#\t--init-hbond-from-native\t = " << init_hbond_from_native << " \n";
     std::cout << "#\t--ca-cutoff\t\t\t = " << ca_cutoff << " \n";
     std::cout << "#\t--ca-skip\t\t\t = " << ca_skip << " \n";
     std::cout << "#\n";
     std::cout << "#\tpath options\n";
     std::cout << "#\t--dbn-parameter-dir\t\t = "	<< dbn_parameter_dir << " \n";
     std::cout << "#\t--dbn-parameter-file\t\t = " << dbn_parameter_file << " \n";
     std::cout << "#\t--sc-dbn-file\t\t\t = " << sc_dbn_file << " \n";
     std::cout << "#\n";
     std::cout << "#\tmove options\n";
     std::cout << "#\t--sc-move-weight\t\t = " << sc_move_weight << " \n";
     //std::cout << "#\t--use-local-sc-move\t\t = " << use_local_sc_move << " \n";
     std::cout << "#\t--local-move-weight\t\t = " << local_sc_move_weight << " \n";
     std::cout << "#\n";
     std::cout << "#\tSampling options:\n";
     std::cout << "#\t--ignore-ss\t\t\t = " << ignore_ss << " \n";
     std::cout << "#\n";
     std::cout << "#\tDehydron options:\n";
     std::cout << "#\t--dehydron-bb-cutoff\t\t = " << dehydron_bb_cutoff << " \n";
     std::cout << "#\t--dehydron-bbsc-cutoff\t\t = " << dehydron_bbsc_cutoff << " \n";
     std::cout << "#\t--dehydron-sc-cutoff\t\t = " << dehydron_sc_cutoff << " \n";
     std::cout << "#\n";

     if (debug) {
          std::cout << "#\n";
          std::cout << "#\tPlease only modify these options if you understand what you do. Otherwise, you\n";
          std::cout << "#\twill most probably end up with a bunch of useless files and a lot of heat.\n";
          std::cout << "#\n";
          std::cout << "#\tadvanced options:\n";
          std::cout << "#\t--mode-move\t\t\t = " << mode_move << " \n";
          std::cout << "#\t--mode-energy\t\t\t = " << mode_energy << " \n";
          std::cout << "#\t--mode-mcmc\t\t\t = " << mode_mcmc << " \n";
          std::cout << "#\t--hbond-weight\t\t\t = " << hbond_weight << " \n";
          std::cout << "#\n";
          std::cout << "#\tadvanced move options:\n";
          std::cout << "#\t--use-semi-local-move\t\t = " << use_semi_local_move << " \n";
          std::cout << "#\t--semi-local-move-weight\t  = " << semi_local_move_weight << " \n";
          //std::cout << "#\t--use-pivot-move\t\t = " << use_pivot_move << " \n";
          //std::cout << "#\t--pivot-move-weight\t\t = " << pivot_move_weight << " \n";
          std::cout << "#\n";
     }
     std::cout << "##############################################################################################################\n";

     if (phaistos_root == "../") {
          std::cout << "#\n";
          std::cout << "#\tW A R N I N G: PHAISTOS_ROOT environmental variable is not set. You are fine as long you \n";
          std::cout << "#\t               are running from your ./build/bin/ directory. Otherwise please set the  \n";
          std::cout << "#\t               variable to match the fully qualified path to your Phaistos build directory.\n";
          std::cout << "#\n";
          std::cout << "##############################################################################################################\n";
     }
}

