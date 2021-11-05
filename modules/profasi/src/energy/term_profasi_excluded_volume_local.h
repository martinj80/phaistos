// profasi_excluded_volume_local.h --- PROFASI local excluded volume energy term
// Copyright (C) 2010 Pengfei Tian, Wouter Boomsma
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

// Implementation of the PROFASI force field
// An effective all-atom potential for proteins, Irback, Mitternacht, Mohanty, PMC Biophysics 2009


#ifndef PROFASI_EXCLUDED_VOLUME_LOCAL_H
#define PROFASI_EXCLUDED_VOLUME_LOCAL_H

#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"
#include "protein/iterators/atom_iterator.h"

namespace phaistos {

//! PROFASI local excluded volume energy term - base class
template<typename DERIVED_CLASS>
class TermProfasiExcludedVolumeLocalBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB>,
                                          public EnergyTermCommonProfasi {
private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

     //! Hydrogens on N-terminus
     std::vector<definitions::AtomEnum> H_NTERM;

     //! Rotate-torsion atoms (define chi angles)
     std::vector<std::vector<definitions::AtomEnum> > side_chain_rt_atoms;

     //! sigmas values for all atom types
     std::vector<double> sigmas;

public:
     
     //! Local settings class.     
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

	  //! Kappa (overall scaling)
          double k_ev;

	  //! Cutoff distance
          double r_cut;

	  //! Atom-specific distance parameters
          double sigma_S;
          double sigma_C;
          double sigma_N;
          double sigma_O;
          double sigma_H;

          //! Constructor.Defines default values for settings object
          Settings(double k_ev = 0.10,
		   double r_cut = 4.3,
		   double sigma_S = 1.77,
		   double sigma_C = 1.75,
		   double sigma_N = 1.53,
		   double sigma_O = 1.42,
		   double sigma_H = 1.00)
          :  k_ev(k_ev),
             r_cut(r_cut),
             sigma_S(sigma_S),
             sigma_C(sigma_C),
             sigma_N(sigma_N),
             sigma_O(sigma_O),
             sigma_H(sigma_H) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "k-ev:" << settings.k_ev << "\n";
               o << "r-cut:" << settings.r_cut << "\n";
               o << "sigma-S:" << settings.sigma_S << "\n";
               o << "sigma-C:" << settings.sigma_C << "\n";
               o << "sigma-N:" << settings.sigma_N << "\n";
               o << "sigma-O:" << settings.sigma_O << "\n";
               o << "sigma-H:" << settings.sigma_H << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }          

     } settings;    //!< Local settings object 


     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiExcludedVolumeLocalBase(ChainFB *chain,
                                        std::string name,
                                        const Settings &settings=Settings(),
                                        RandomNumberEngine *random_number_engine = &random_global)
	  : EnergyTermCommon(chain, name, settings, random_number_engine) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          //! Import phaistos::make_vector namespace
          using namespace vector_utils;

          side_chain_rt_atoms.resize(20);

          H_NTERM = make_vector ( H1, H2, H3, H);

          // Heavy sidechain atoms defining primary branch
          side_chain_rt_atoms[ALA] = make_vector( CA, CB);
          side_chain_rt_atoms[ARG] = make_vector( CA, CB, CG, CD, NE);
          side_chain_rt_atoms[ASN] = make_vector( CA, CB, CG);
          side_chain_rt_atoms[ASP] = make_vector( CA, CB, CG);
          side_chain_rt_atoms[CYS] = make_vector( CA, CB, SG);
          side_chain_rt_atoms[GLU] = make_vector( CA, CB, CG, CD);
          side_chain_rt_atoms[GLN] = make_vector( CA, CB, CG, CD);
          side_chain_rt_atoms[HIS] = make_vector( CA, CB, CG);
          side_chain_rt_atoms[ILE] = make_vector( CA, CB, CG1, CD1);
          side_chain_rt_atoms[LEU] = make_vector( CA, CB, CG, CD1);
          side_chain_rt_atoms[LYS] = make_vector( CA, CB, CG, CD, CE, NZ);
          side_chain_rt_atoms[MET] = make_vector( CA, CB, CG, SD, CE);
          side_chain_rt_atoms[PHE] = make_vector( CA, CB, CG);
          side_chain_rt_atoms[SER] = make_vector( CA, CB, OG);
          side_chain_rt_atoms[THR] = make_vector( CA, CB, OG1);
          side_chain_rt_atoms[TRP] = make_vector( CA, CB, CG);
          side_chain_rt_atoms[TYR] = make_vector( CA, CB, CG); 
          side_chain_rt_atoms[VAL] = make_vector( CA, CB, CG1);

	  // Assign  every atom a sigma value
	  sigmas.resize(ATOM_ENUM_SIZE);
	  for (int i=0; i<ATOM_ENUM_SIZE; ++i) {
	       if( is_atom_XC(AtomEnum(i)) ){
		    sigmas[i] = settings.sigma_C;
	       }
	       else if(is_atom_XH(AtomEnum(i)) ){
		    sigmas[i] = settings.sigma_H;
	       }
	       else if(is_atom_XN(AtomEnum(i)) ){
		    sigmas[i] = settings.sigma_N;
	       }
	       else if(is_atom_XO(AtomEnum(i)) ){
		    sigmas[i] = settings.sigma_O;
	       }
	       else if(is_atom_XS(AtomEnum(i)) ){
		    sigmas[i] = settings.sigma_S;
	       }
	  }

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiExcludedVolumeLocalBase(const TermProfasiExcludedVolumeLocalBase &other,
                                        RandomNumberEngine *random_number_engine,
                                        int thread_index, ChainFB *chain)
	  : EnergyTermCommon(other, random_number_engine, thread_index, chain),
	    H_NTERM(other.H_NTERM),
	    side_chain_rt_atoms(other.side_chain_rt_atoms),
	    sigmas(other.sigmas) {
     }

     //! Calculate the excluded volume of the given pair of atoms
     //! \param atom_i The first atom in the excluded volume interaction
     //! \param atom_j The second atom in the excluded volume interaction
     //! \return Excluded volume potential between the pair of atoms
     double localpairs_energy( Atom * atom_i, Atom * atom_j ){

          // Atom distance
          double r_ij = (atom_i->position - atom_j->position).norm();

          if (r_ij > settings.r_cut)
               return 0;

          double sigma_i = sigmas[atom_i->atom_type];
          double sigma_j = sigmas[atom_j->atom_type];
          double asa = -7*pow( (sigma_i + sigma_j) / settings.r_cut, 12);
          double bsa = 6*pow( (sigma_i + sigma_j) / settings.r_cut, 12) / settings.r_cut / settings.r_cut;

          double energy = asa + bsa*r_ij*r_ij + pow( (sigma_i + sigma_j) / r_ij, 12);
          return energy;
     }


     //! Calculate the local excluded volume of the given pair of atoms(on backbone) on the backbone
     //! \param it1 Residue iterator
     double res_localbb_energy(ResidueIterator<ChainFB> &it1) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          double E_res_locbb=0.0;

          // Define atoms which are involved the local excluded volume calculation
          Atom *atom_1_HA, *atom_1_C, *atom_1_CB, *atom_1_N, *atom_2_N, *atom_1_O ;

          atom_1_C = (*it1)[C];
          atom_1_N = (*it1)[N];
          atom_1_O = (*it1)[O];

          if(it1->residue_type == GLY){
               atom_1_HA = (*it1)[HA3];
               atom_1_CB = (*it1)[HA2];
          } else {
               atom_1_HA = (*it1)[HA];
               atom_1_CB = (*it1)[CB];
          }

          if(it1->terminal_status != CTERM){
               atom_2_N = (*(it1+1))[N];
          } else {
               atom_2_N = (*it1)[OXT];
          }

          E_res_locbb += localpairs_energy(atom_1_N, atom_1_O);
          E_res_locbb += localpairs_energy(atom_1_N, atom_2_N);

          E_res_locbb += localpairs_energy(atom_1_HA, atom_1_O);
          E_res_locbb += localpairs_energy(atom_1_HA, atom_2_N);

          E_res_locbb += localpairs_energy(atom_1_CB, atom_1_O);
          E_res_locbb += localpairs_energy(atom_1_CB, atom_2_N);

          if(it1->residue_type != PRO){

               int H_num = it1-> terminal_status == NTERM ? 3:4;

               for( int i = it1-> terminal_status == NTERM ? 0:3; i< H_num; i++){
                    Atom *atom_H = (*it1)[H_NTERM[i]];

                    E_res_locbb += localpairs_energy(atom_H, atom_1_HA);
                    E_res_locbb += localpairs_energy(atom_H, atom_1_CB);
                    E_res_locbb += localpairs_energy(atom_H, atom_1_C);

               }

               if(it1-> terminal_status != NTERM){

                    ResidueFB* it0 = it1->get_neighbour(-1);
                    Atom *atom_0_C = (*it0)[C];

                    E_res_locbb += localpairs_energy(atom_0_C, atom_1_HA);
                    E_res_locbb += localpairs_energy(atom_0_C, atom_1_CB);
                    E_res_locbb += localpairs_energy(atom_0_C, atom_1_C);
               }

          }

          return E_res_locbb;
     }

     // Calculate the local excluded volume of the given pair of atoms(in the same residue) on the backbone
     //! \param it2 Residue iterator
     double res_localrt_energy(ResidueIterator<ChainFB> &it2){

          // Import protein definitions (such as residue names)
          using namespace definitions;

          double E_res_locrt=0.0;

          int res = it2->residue_type;

          if (res==PRO || res == GLY){
               return 0;
          }
          
          // If residue is TYR, make pair of HH with CE1 and CE2
          if (res == TYR){

               Atom *atom_HH = (*it2)[HH];
               Atom *atom_CE1= (*it2)[CE1];
               Atom *atom_CE2 = (*it2)[CE2];

               E_res_locrt += localpairs_energy(atom_HH, atom_CE1);
               E_res_locrt += localpairs_energy(atom_HH, atom_CE2);
          }

          const std::vector<AtomEnum> & rt_atoms = side_chain_rt_atoms[res];

          int rt_size = rt_atoms.size();

          for(int i=0; i< (rt_size -1);i++){

               Atom *atom_1 = (*it2)[rt_atoms[i]];
               Atom *atom_2 = (*it2)[rt_atoms[i+1]];

               CovalentBondIterator<ChainFB> it_a(atom_1, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);

               for (; !it_a.end(); ++it_a) {

                    Atom *atom_pair_1 = &*it_a;

                    // Branched sidechains
                    if((res == ILE) || (res == LEU) || (res == THR) || (res == VAL) ) {///

                         if(i>0) { // If atom_1 != CA

                              Atom *atom_0 = (*it2)[rt_atoms[i-1]];

                              // Test if atom1 has a carbon neighbour that is not in the rt_atoms list
                              // if so, it defines a second branch
                              if((is_atom_XC(atom_pair_1->atom_type)) &&
                                        (atom_pair_1->atom_type != atom_0->atom_type) &&
                                        (atom_pair_1->atom_type != atom_2->atom_type )){

                                   // Iterate over the neighbours of the secondary branch carbon
                                   CovalentBondIterator<ChainFB> it_c(atom_pair_1, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
                                   for (; !it_c.end(); ++it_c){

                                        Atom *atom_pair_1_nb  = &*it_c;

                                        for (CovalentBondIterator<ChainFB> it_a(atom_1, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY); !it_a.end(); ++it_a) {
                                             Atom *atom_1_nb = &*it_a;

                                             if((atom_1_nb->atom_type != atom_pair_1->atom_type) &&
                                                       (atom_pair_1_nb -> atom_type != atom_1-> atom_type))
                                                  E_res_locrt += localpairs_energy(atom_1_nb, atom_pair_1_nb);
                                        }
                                   }
                              }
                         }
                    }

		    // Primary branch
                    for (CovalentBondIterator<ChainFB> it_b(atom_2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
                         !it_b.end(); ++it_b) {

                         Atom *atom_pair_2 = &*it_b;

                         if((atom_pair_1->atom_type != atom_2->atom_type ) &&
                                   ( atom_pair_2->atom_type != atom_1->atom_type ))
                              E_res_locrt += localpairs_energy(atom_pair_1, atom_pair_2);
                    }
               }
          }

          return E_res_locrt;
     }

};


//! PROFASI local excluded volume energy term
class TermProfasiExcludedVolumeLocal: public TermProfasiExcludedVolumeLocalBase<TermProfasiExcludedVolumeLocal> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiExcludedVolumeLocalBase<TermProfasiExcludedVolumeLocal> TermProfasiExcludedVolumeLocalBase;

public:

     //! Use same settings as base class
     typedef TermProfasiExcludedVolumeLocalBase::Settings Settings;

     //! Settings
     const Settings settings;
  
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiExcludedVolumeLocal(ChainFB *chain,
                                    const Settings &settings=Settings(),
                                    RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiExcludedVolumeLocalBase(chain, "profasi-excluded-volume-local", settings, random_number_engine),
            settings(settings){
     }

     // Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain     
     TermProfasiExcludedVolumeLocal(const TermProfasiExcludedVolumeLocal &other, 
                                    RandomNumberEngine *random_number_engine,
                                    int thread_index, ChainFB *chain)
          : TermProfasiExcludedVolumeLocalBase(other, random_number_engine, thread_index, chain),
            settings(other.settings){}

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     //! \return energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          double energy = 0.0;
          double energy_bb = 0.0;
          double energy_rt = 0.0;

          //Go through all the residues on the chain
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){

               energy_bb = this->res_localbb_energy(it);
               energy_rt = this->res_localrt_energy(it);

               energy += (energy_bb + energy_rt);
          }

          energy *= this->settings.k_ev;

          return this->profasi_energy_in_kcal_per_mol * energy;
     }

};


//! PROFASI local excluded volume energy term - cached version
class TermProfasiExcludedVolumeLocalCached: public TermProfasiExcludedVolumeLocalBase<TermProfasiExcludedVolumeLocalCached> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiExcludedVolumeLocalBase<TermProfasiExcludedVolumeLocalCached> TermProfasiExcludedVolumeLocalBase;

public:
     //! Cached iterator
     CachedIterator<ResidueIterator<ChainFB> > cached_it;

     //! Use same settings as base class
     typedef TermProfasiExcludedVolumeLocalBase::Settings Settings;

     //! Settings
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiExcludedVolumeLocalCached(ChainFB *chain,
                                          const Settings &settings=Settings(),
                                          RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiExcludedVolumeLocalBase(chain, "profasi-excluded-volume-local-cached", settings, random_number_engine),
         cached_it(*chain),
         settings(settings){}

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiExcludedVolumeLocalCached(const TermProfasiExcludedVolumeLocalCached &other, 
                                          RandomNumberEngine *random_number_engine,
                                          int thread_index, ChainFB *chain)
          : TermProfasiExcludedVolumeLocalBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings){}

     //! Evaluate chain energy
     //! \param move_info Object containing information about the last executed move
     //! \return energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {
          
          int start_index = 0;
          int end_index = this->chain->size();

          if (move_info) {
               start_index = std::max(0,move_info->modified_angles_start-1);
               end_index = std::min(this->chain->size(), move_info->modified_angles_end+1);
          }
          

          for (cached_it(*this->chain,
                         start_index,
                         end_index); !cached_it.end(); ++cached_it) {
               
               double energy_bb = this->res_localbb_energy(cached_it);
               double energy_rt = this->res_localrt_energy(cached_it);

               cached_it.register_contribution(energy_bb + energy_rt);
          }
          
          double energy = cached_it.compute_total();
          
          energy *= this->settings.k_ev;

          return this->profasi_energy_in_kcal_per_mol * energy;
     }

     //! Accept last energy evaluation
     void accept() {
          cached_it.accept();
     }

     //! Reject last energy evaluation
     void reject() {
          cached_it.reject();
     }

};

}

#endif
