// test_opls.cpp --- Test of OPLS energy class
// Copyright (C) 2008 Kristoffer Enøe Johansson, Wouter Boomsma
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

#include "protein/chain_fb.h"
#include "energy/energy.h"

#include "energy/term_opls_non_bonded.h"
#include "energy/term_gbsa.h"
#include "energy/term_opls_vdw.h"
#include "energy/term_opls_charge.h"
#include "energy/term_opls_torsion.h"
#include "energy/term_opls_angle_bend.h"
#include "energy/term_opls_bond_stretch.h"
#include "energy/term_opls_imptor.h"

#include "energy/observable.h"
#include "energy/observable_collection.h"

using namespace std;
using namespace phaistos;
using namespace phaistos::definitions;

//! Method to evaluate a PDB file using energy terms
void test_terms(ChainFB *chain, int n=1) {

     // Create Energy class 
     Energy<ChainFB> energy(chain);

     // Add terms
     // energy.add_term( new TermOplsNonBonded<ChainFB,double>(chain) );
     energy.add_term( new TermOplsNonBondedCached<double>(chain) );
     
     // energy.add_term( new TermGbsa(chain) );
     // energy.add_term( new TermGbsaCached(chain) );

     // energy.add_term( new TermOplsVdw(chain) );
     // energy.add_term( new TermOplsVdwCached(chain) );

     // energy.add_term( new TermOplsCharge(chain) );
     // energy.add_term( new TermOplsChargeCached(chain) );

     // energy.add_term( new TermOplsAngleBend(chain) );
     energy.add_term( new TermOplsAngleBendCached(chain) );

     energy.add_term( new TermOplsTorsion(chain) );

     energy.add_term( new TermOplsBondStretch(chain) );

     // energy.add_term( new TermOplsImptor(chain) );

     // Evaluate energy
     cout << "\nEvaluate terms " << n << " times\n" << endl;
     for (int i=0; i<n; i++)
          energy.evaluate();

     // Output
     cout << energy << endl;;
}



int main(int argc, char *argv[]) {

     if (argc < 2) {
          cout << "USAGE: ./test_opls <pdb-file>" <<endl;;
          exit(1);
     }

     // Create chain from PDB filename
     string pdb_filename = argv[1];
     ChainFB chain(pdb_filename, ALL_ATOMS);

     // Add atoms missing in the pdb structure
     chain.add_atoms(ALL_PHYSICAL_ATOMS);

     test_terms(&chain, 100);
}
