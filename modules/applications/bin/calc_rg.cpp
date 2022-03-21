// calcRG.cpp --- Calculate the radius of gyration for a structure
// Copyright (C) 2010 Tim Harder, Wouter Boomsma
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

#include "energy/term_rg.h"
#include "protein/chain_fb.h"
#include "protein/iterators/atom_iterator.h"
#include "protein/chain_ca.h"
#include "protein/definitions.h"
#include "utils/vector_matrix_3d.h"

using namespace phaistos;

//! Option class
class Options {
public:

     //! Name of PDB file
     std::string pdb_filename;

     //! Verbose mode     
     bool verbose;

     //! C-alpha only mode
     bool only_CA;

     //! Constructor
     Options(int argc, char *argv[]) {
          get_command_line(argc, argv);

     }

     //! Default values
     void init() {
          pdb_filename = "";
          verbose = false;
          only_CA = false;
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
               } else if (argName == std::string("ca") || argName == std::string("v")) {
                    this->only_CA = true;
               }
          }

          std::string usage = std::string("Usage ") + argv[0] + " [--ca] pdbfile1 \n\n";
          if (arguments.size() != 1) {
               std::cerr << "Error: " << argv[0] << "requires a PDB file as argument.\n" << usage;
               exit(0);
          }
          pdb_filename = arguments[0];
     }

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, Options &options) {
          o << "### Options ###\n";
          o << "PDB file: " << options.pdb_filename << "\n";
          o << "only_CA: " << options.only_CA << "\n";
          o << "verbose: " << options.verbose << "\n";
          return o;
     }

};

//! Calculate rg
template <typename CHAINTYPE, definitions::IterateEnum ITERATION_MODE>
double calc_rg_local (AtomIterator<CHAINTYPE, ITERATION_MODE, Vector_3D> it){
     double res = 0.0;
     Vector_3D CM = center_of_mass(it);
     int counter = 0;
     for (; !it.end(); ++it) {
          res += (*it-CM)*(*it-CM);
          counter++;
     }
     res=sqrt(res/(double)counter);
     return res;
}


//! Main function
int main(int argc, char *argv[]) {

     Options options(argc, argv);

     if (options.verbose) {
          std::cout << options.pdb_filename  << "\t";
     }
     if(options.only_CA){
             ChainCA chain1(options.pdb_filename);
             std::cout << calc_rg_local(AtomIterator<ChainCA, definitions::BACKBONE, Vector_3D>(chain1))  << "\n";
     } else {
             ChainFB chain1(options.pdb_filename);
             std::cout << calc_rg_local(AtomIterator<ChainFB, definitions::ALL, Vector_3D>(chain1))  << "\n";
     }
     return 0.;
}


