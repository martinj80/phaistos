// typhon_options.h --- option parser for typhon.cpp
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
#include <string>
#include <iostream>


//! Typhon command line options class
class TyphonOptions {

     //! Command line string
     std::string commandline;

public:
     
     //! Internal exception for unknown options
     class UnknownOption: public std::exception {
          virtual const char* what() const throw() {
               return "Unknown commandline option .. bailing out!\n";
          }
     } UnknownOption;

     //! Phaistos root directory
     std::string phaistos_root;

     //! Where to put all the generated structures
     std::string output_directory;

     //! How often shall a structure be saved
     int dump_interval;

     //! Store the GIT vectors as well?
     bool dump_git;

     //! Which starting structure to use
     std::string pdb_file;

     //! File used for condition on different secondary structure
     std::string ss_file;

     //! Restore h-bond network from file
     std::string restore_network;

     //! Torus DBN Parameter directory
     std::string dbn_parameter_dir;

     //! Torus DBN Parameter file
     std::string dbn_parameter_file;

     //! How many samples to generate in the simulation
     int iterations;

     //! How many iterations to skip before starting to sample
     int burnin;

     //! Seed for the random number generators
     int seed;

     //! How often to reinitialize to the native structure
     int reinitialize_interval;

     //! Which DBN to use for the sidechain moves
     std::string sc_dbn_file;

     //! Weight of sidechain move
     int sc_move_weight;

     // Whether to use local sidechain move
     bool use_local_sc_move;

     //! Weight of local sidechain move
     double local_sc_move_weight;

     // Whether to use semi local move
     bool use_semi_local_move;

     //! Weight of semi local move
     double semi_local_move_weight;

     // Whether to use pivot move
     bool use_pivot_move;

     //! Weight of pivot move
     double pivot_move_weight;

     //@{
     //! some modes .. mostly for debugging
     std::string mode_mcmc;
     std::string mode_move;
     std::string mode_energy;
     //@}

     //! Analyze and print the hbond network only, no sampling
     bool mode_analyze;

     //! Generate a Python script that will visualize the network in PyMOL
     bool generate_pymol;

     //! Weight of the the hbond term
     double hbond_weight;

     //! Disregard any hbonds found (overwrites the restore network)
     bool ignore_hbonds;

     //! Disregard any side chain - side chain hbonds found (overwrites the restore network)
     bool ignore_sc_hbonds;

     //! Disregard any backbone - side chain hbonds found (overwrites the restore network)
     bool ignore_bbsc_hbonds;

     //! Disregard any gaussian found (overwrites the restore network)
     bool ignore_gaussian;

     //! Disregard any ssbond found (overwrites the restore network)
     bool ignore_ssbond;

     //! disregard any interactions with ligands (overwrites the restore network)
     bool ignore_fixed;

     //! Per default we use parameter estimated from the top 500 set. Set this 
     //! flag to initialize hbond distances to the native structure
     bool init_hbond_from_native;

     //! Do not prune/optimize the hbond network
     bool no_prune;

     //! Dump a "native" structure with all atoms at the beginning of the run
     bool create_phaistos_input;

     //! Do not fix the secondary structure for sampling
     bool ignore_ss;

     //! Cutoff distance in Angstrom to be considered a Calpha contact
     float ca_cutoff;

     //! How many residues to skip along the chain before considering a Calpha contact
     int ca_skip;

     //@{
     //! dehyrdon cutoffs
     int dehydron_bb_cutoff;
     int dehydron_bbsc_cutoff;
     int dehydron_sc_cutoff;
     //@}

     //! Produce more output
     bool verbose;

     //! Print debugging information
     bool debug;


     //! Constructor
     TyphonOptions(int argc, char *argv[]);

     //! Destructor
     ~TyphonOptions() ;

     //! Initializer
     void init() ;

     //! Parse the command line
     void get_command_line(int argc, char *argv[]) ;

     //! Display usage information
     void print_usage() ;

     //! Display available options
     void print_options() ;

};
