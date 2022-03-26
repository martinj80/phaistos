// calcQScore.cpp --- Calculate the Q Score between two proteins
//                    As described in: Assessment of CASP8 structure predictions for 
//                                     template free targets. M Ben-David, O Noivirt-Brik, 
//                                     A Paz, J Prilusky, JL Sussman, Y Levy, 2009
// Copyright (C) 2010 Tim Harder
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


#include "energy/term_rmsd.h"
#include "protein/chain_ca.h"
#include "protein/iterators/atom_iterator.h"

using namespace phaistos;

//! Option class
class Options {
public:

     //! Name of first PDB file
     std::string pdb_filename1;

     //! Name of second PDB file
     std::string pdb_filename2;

     //! Verbose mode
     bool verbose;

     //! Whether RMSD should be calculated
     bool calc_rmsd;

     //! Constructor
     Options(int argc, char *argv[]) {
          get_command_line(argc, argv);

     }

     //! Default values
     void init() {
          pdb_filename1 = "";
          pdb_filename2 = "";
          verbose = false;
          calc_rmsd = false;
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
               } else if (argName == std::string("calc_rmsd")|| argName == std::string("r")) {
                    this->calc_rmsd = true;
               }else if (argName == std::string("vr")|| argName == std::string("rv")) {
                    this->verbose = true;
                    this->calc_rmsd = true;
               }
          }

          std::string usage = std::string("Usage ") + argv[0] + " pdbfile1 pdbfile2 n\n";
          if (arguments.size() != 2) {
               std::cerr << "Error: " << argv[0] << "requires two PDB file arguments.\n" << usage;
               exit(0);
          }
          pdb_filename1 = arguments[0];
          pdb_filename2 = arguments[1];
     }

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, Options &options) {
          o << "### Options ###\n";
          o << "PDB file 1: " << options.pdb_filename1 << "\n";
          o << "PDB file 2: " << options.pdb_filename2 << "\n";
          o << "Verbose: " << options.verbose << "\n";
          o << "RMSD: " << options.calc_rmsd << "\n";
          return o;
     }

};

//! Q-score class.
//! As described in:  Assessment of CASP8 structure predictions for 
//! template free targets. M Ben-David, O Noivirt-Brik, A Paz, 
//! J Prilusky, JL Sussman, Y Levy, 2009
class QScore {

     //! Internal matrix data
     double **matrix;

     //! Length
     int len;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain object
     QScore (ChainCA &chain) {
          this->len = chain.size();

          matrix= new double *[len];
          for (int i=0; i<len; i++) {
               matrix[i] = new double[len];
          }

          for (ResidueIterator<ChainCA> res1(chain); !res1.end(); ++res1) {
               if (not (res1)->has_atom(definitions::CA))
                    continue;
               Atom *CA1 = (*res1)[definitions::CA];
               int i1=res1->index;
               // and lets go for the inner loop
               for (ResidueIterator<ChainCA> res2(chain); !res2.end(); ++res2) {
                    if (res2->index <= res1->index || not (res1)->has_atom(definitions::CA))
                         continue;
                    Atom *CA2 = (*res2)[definitions::CA];
                    int i2=res2->index;
                    matrix[i1][i2] = (CA1->position - CA2->position).norm();
               }
          }
     }

     //! Desctructor
     ~QScore () {
          for( int i = 0 ; i < this->len ; i++ )
               delete [] this->matrix[i] ;
          delete [] matrix ;
     }

     //! Return length
     int length () {
          return this->len;
     }

     //! Retrieve matrix element
     double get(int i, int j) {
          if (i>=len || j>=len)
               return 0.;
          return this->matrix[i][j];
     }

     //! Calculate score
     double get_score(QScore &target) {
          if (this->len != target.length())
               return 0;

          double qsum = 0.;
          int count =0;
          for (int i=0; i<len; i++) {
               for (int j=i+1; j<len; j++) {
                    count++;
                    qsum+=exp(-((this->matrix[i][j]-target.get(i,j))*(this->matrix[i][j]-target.get(i,j)) ));
               }
          }
          if (count!=0)
               return qsum/count;
          return 0.;
     }
};


//! Main function
int main(int argc, char *argv[]) {

     Options options(argc, argv);

     ChainCA chain1(options.pdb_filename1);
     ChainCA chain2(options.pdb_filename2);

     QScore decoy(chain1);
     QScore target(chain2);

     // Output the filenames
     if (options.verbose) {
          std::cout << options.pdb_filename1 << "\t" << options.pdb_filename2 << "\t" ;
     }
     // Output the Q score
     std::cout << decoy.get_score(target) << "\t";

     // Optionally output the RMSD
     if (options.calc_rmsd) {
          std::cout << calc_rmsd<definitions::CA_ONLY>(chain1, chain2);
     }

     std::cout << "\n";

     return 0.;
}


