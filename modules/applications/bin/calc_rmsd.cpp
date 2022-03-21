// calc_rmsd.cpp --- Calculate RMSD between two proteins
// Copyright (C) 2006-2013 Wouter Boomsma, Jan Valentin
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


#include "protein/chain_ca.h"
#include "protein/chain_fb.h"

#include "energy/term_rmsd.h"
#include "protein/iterators/atom_iterator.h"

using namespace phaistos;

//! Option class
class Options {
public:

     //! Name of first PDB file
     std::string pdb_filename1;

     //! Name of second PDB file
     std::string pdb_filename2;
     
     int residue_start;
     
     int residue_end;

     //! C-alpha only mode
     bool only_CA;

     //! Verbose mode
     bool verbose;

     //! Constructor
     Options(int argc, char *argv[]) {          
          get_command_line(argc, argv);

     }

     //! Default values
     void init() {
          pdb_filename1 = "";
          pdb_filename2 = "";
          residue_start = 0;
          residue_end = -1;
          only_CA = false;
          verbose = false;
     }

     //! Parse the command line
     void get_command_line(int argc, char *argv[]) {
          init();
          std::vector<std::string> arguments;

          for (int i=1; i<argc; i++) {
               std::string arg = argv[i];
               std::string argName = "";

               if (arg.substr(0,2) == std::string("--")) {
                    argName = arg.substr(2,arg.length()-2);
               } else if (arg.substr(0,1) == std::string("-")) {
                    argName = arg.substr(1,arg.length()-1);
               } else {
                    arguments.push_back(std::string(argv[i]));
                    continue;
               }

               if (argName == std::string("verbose") || argName == std::string("v")) {
                    this->verbose = true;
               } else if (argName == std::string("help") || argName == std::string("h")) {
                    print_usage(argv);
                    exit(0);
               } else if (argName == std::string("residue-start")) {
                    i++;
                    this->residue_start = std::atoi(argv[i]);
               } else if (argName == std::string("residue-end")) {
                    i++;
                    this->residue_end = std::atoi(argv[i]);
               } else if (argName == std::string("ca")) {
                    this->only_CA = true;
               }
          }

          if (arguments.size() != 2) {
               std::cerr << "Error: " << argv[0] << "requires two PDB file arguments.\n";
               print_usage(argv);
               exit(0);
          }
          pdb_filename1 = arguments[0];
          pdb_filename2 = arguments[1];
     }

     void print_usage(char *argv[]) {
          // generate usage informations
          std::cout << "Usage: " << argv[0] << " pdbfile1 pdbfile2 [options]\n\n";
          std::cout << "\toptions: \n";
          std::cout << "\t-h  --help          \tprint help \n";
          std::cout << "\t-v  --verbose       \tbe talky? [default; " << this->verbose << "] \n";
          std::cout << "\t--residue-start INT \t[default; INT=" << this->residue_start << "] \n";
          std::cout << "\t--residue-end INT   \t[default; INT=" << this->residue_end << "] \n";
          std::cout << "\t--ca                \tUse only CA atoms [default; " << this->only_CA << "] \n";
     }

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, Options &options) {
          o << "### Options ###\n";
          o << "PDB file 1: " << options.pdb_filename1 << "\n";
          o << "PDB file 2: " << options.pdb_filename2 << "\n";
          o << "residue-start: " << options.residue_start << "\n";
          o << "residue-end: " << options.residue_end << "\n";
          o << "ca: " << options.only_CA << "\n";
          o << "Verbose: " << options.verbose << "\n";
          return o;
     }

};


//! Main function
int main(int argc, char *argv[]) {

     Options options(argc, argv);
     
     ChainCA chain1(options.pdb_filename1);
     ChainCA chain2(options.pdb_filename2);

     if (options.verbose) {
          std::cout << options;
          printf( "# chain sizes: c1: %d, c2: %d\n", chain1.size(),chain2.size() );
     }

     // Import protein definitions (such as residue names)
     using namespace definitions;

     if (options.only_CA) {
          std::cout << calc_rmsd<definitions::CA_ONLY>(chain1, chain2, options.residue_start, options.residue_end) << "\n";
     } else {
          std::cout <<  calc_rmsd<definitions::ALL>(chain1, chain2, options.residue_start, options.residue_end) << "\n";
     }

     return 0;
}
