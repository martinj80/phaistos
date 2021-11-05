
// test_profasi.cpp --- Test of profasi energy class
// Copyright (C) 2010 Pengfei Tian
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
#include "energy/term_profasi_local_sidechain.h"
#include "energy/term_profasi_hydrophobicity.h"
#include "energy/term_profasi_hydrogen_bond.h"
#include "energy/term_profasi_sidechain_charge.h"
#include "energy/term_profasi_excluded_volume_local.h"
#include "energy/term_profasi_local.h"
#include "energy/term_profasi_excluded_volume.h"
#include "energy/term_profasi_proline_phi_torsion.h"

#include "energy/observable.h"
#include "energy/observable_collection.h"

using namespace std;
using namespace phaistos;

void test_profasi(std::string pdb_filename) {

     // Create chain from PDB filename
     ChainFB chain(pdb_filename,definitions:: ALL_ATOMS);
     // Add atoms missing in the pdb structure
     chain.add_atoms(definitions::ALL_PHYSICAL_ATOMS);

     // Create Energy class
     Energy<ChainFB> energy(&chain);

     // Add term
     //typename TermProfasiTorsion::Settings settings;
     //settings.weight = 1.0/(3.2976E-27*6.022E23 * 315);
     //energy.addTerm(new TermProfasiSidechainTorsion(&chain));
//     typename TermProfasiExcludedVolume::Settings settings_excluded_volume;
//     settings_excluded_volume.weight = 1;
//     energy.add_term(new TermProfasiExcludedVolume(&chain, settings_excluded_volume));

     TermProfasiExcludedVolumeCached::Settings settings_excluded_volume_cached;
     settings_excluded_volume_cached.weight = 1;
     energy.add_term(new TermProfasiExcludedVolumeCached(&chain, settings_excluded_volume_cached));

//     typename TermProfasiExcludedVolumeLocal::Settings settings_excluded_volume_local;
//     settings_excluded_volume_local.weight = 1;
//     energy.add_term(new TermProfasiExcludedVolumeLocal(&chain, settings_excluded_volume_local));

     TermProfasiExcludedVolumeLocalCached::Settings settings_excluded_volume_local_cached;
     settings_excluded_volume_local_cached.weight = 1;
     energy.add_term(new TermProfasiExcludedVolumeLocalCached(&chain, settings_excluded_volume_local_cached));

//     typename TermProfasiLocalSidechain::Settings settings_local_sidechain;
//     settings_local_sidechain.weight = 1;
//     energy.add_term(new TermProfasiLocalSidechain(&chain, settings_local_sidechain));
     TermProfasiLocalCached::Settings settings_local_cached;
     settings_local_cached.weight = 1;
     energy.add_term(new TermProfasiLocalCached(&chain, settings_local_cached));


     TermProfasiLocalSidechainCached::Settings settings_local_sidechain_cached;
     settings_local_sidechain_cached.weight = 1;
     energy.add_term(new TermProfasiLocalSidechainCached(&chain, settings_local_sidechain_cached));


//     typename TermProfasiHydrogenBond::Settings settings_hydrogen_bond;
//     settings_hydrogen_bond.weight = 1;
//     energy.add_term(new TermProfasiHydrogenBond(&chain, settings_hydrogen_bond));

//     typename TermProfasiHydrogenBondCached::Settings settings_hydrogen_bond_cached;
//     settings_hydrogen_bond_cached.weight = 1;
//     energy.add_term(new TermProfasiHydrogenBondCached(&chain, settings_hydrogen_bond_cached));

     TermProfasiHydrogenBondImproved::Settings settings_hydrogen_bond_improved;
     settings_hydrogen_bond_improved.weight = 1;
     energy.add_term(new TermProfasiHydrogenBondImproved(&chain, settings_hydrogen_bond_improved));

//     typename TermProfasiHydrophobicity::Settings settings_hydrophobicity;
//     settings_hydrophobicity.weight = 1;
//     energy.add_term(new TermProfasiHydrophobicity(&chain, settings_hydrophobicity));

//     typename TermProfasiHydrophobicityCached::Settings settings_hydrophobicity_cached;
//     settings_hydrophobicity_cached.weight = 1;
//     energy.add_term(new TermProfasiHydrophobicityCached(&chain, settings_hydrophobicity_cached));

     TermProfasiHydrophobicityImproved::Settings settings_hydrophobicity_improved;
     settings_hydrophobicity_improved.weight = 1;
     energy.add_term(new TermProfasiHydrophobicityImproved(&chain, settings_hydrophobicity_improved));

//     typename TermProfasiLocal::Settings settings_local;
//     settings_local.weight = 1;
//     energy.add_term(new TermProfasiLocal(&chain, settings_local));


//     typename TermProfasiSidechainCharge::Settings settings_sidechain_charge;
//     settings_sidechain_charge.weight = 1;
//     energy.add_term(new TermProfasiSidechainCharge(&chain, settings_sidechain_charge));

//     typename TermProfasiSidechainChargeCached::Settings settings_sidechain_charge_cached;
//     settings_sidechain_charge_cached.weight = 1;
//     energy.add_term(new TermProfasiSidechainChargeCached(&chain, settings_sidechain_charge_cached));


     TermProfasiSidechainChargeImproved::Settings settings_sidechain_charge_improved;
     settings_sidechain_charge_improved.weight = 1;
     energy.add_term(new TermProfasiSidechainChargeImproved(&chain, settings_sidechain_charge_improved));

     //   typename TermProfasiProlinePhiTorsion::Settings settings_proline_phi_torsion;
     // settings_proline_phi_torsion.weight = 1;
     // energy.add_term(new TermProfasiProlinePhiTorsion(&chain, settings_proline_phi_torsion));


//     typename TermProfasiHydrophobicity::Settings settings_hydrophobicity;
//     settings_hydrophobicity.weight = 1;
//     energy.addTerm(new TermProfasiHydrophobicity(&chain, settings_hydrophobicity));
//
//     typename TermProfasiSidechainCharge::Settings settings_sidechain_charge;
//     settings_sidechain_charge.weight = 1;
//     energy.addTerm(new TermProfasiSidechainCharge(&chain, settings_sidechain_charge));

     energy.evaluate();
     std::cout <<"\n\n"<< energy << "\n";
}



int main(int argc, char *argv[]) {

     if (argc < 2) {
          std::cout<<"usage: ./test_profasi <struc.pdb>"<<std::endl;
          exit(1);
     }

     std::string pdb_filename = argv[1];

     test_profasi(pdb_filename);
}
