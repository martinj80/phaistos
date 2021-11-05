// profasi_excluded_volume.h --- PROFASI excluded volume energy term
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


#ifndef PROFASI_EXCLUDED_VOLUME_H
#define PROFASI_EXCLUDED_VOLUME_H

#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"
#include "protein/iterators/atom_iterator.h"
#include "protein/iterators/dof_iterator.h"
#include "protein/residue.h"

namespace phaistos {

//! PROFASI excluded volume energy term - base class
template<typename DERIVED_CLASS>
class TermProfasiExcludedVolumeBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB>,
                                     public EnergyTermCommonProfasi {
private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

     //! Sigmas values for all atom types
     std::vector<double> sigmas;

public:

     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

	  //! Kappa (overall scaling)
          double k_ev;

	  //! Cutoff distance
          double r_c;

	  //! Lambda parameter
          double lambda;

          //! Cutoff distance for cached version = cutoff distance * lambda parameter
          double r_cut;

	  //! Atom-specific distance parameters
          double sigma_S;
          double sigma_C;
          double sigma_N;
          double sigma_O;
          double sigma_H;

          //! Constructor. Defines default values for settings object
          Settings( double k_ev = 0.10,
                    double r_c = 4.3,
                    double lambda = 0.75,
                    double r_cut = 3.225,
                    double sigma_S = 1.77,
                    double sigma_C = 1.75,
                    double sigma_N = 1.53,
                    double sigma_O = 1.42,
                    double sigma_H = 1.00)
                  : k_ev(k_ev),
                    r_c(r_c),
                    lambda(lambda),
                    r_cut(r_cut),
                    sigma_S(sigma_S),
                    sigma_C(sigma_C),
                    sigma_N(sigma_N),
                    sigma_O(sigma_O),
                    sigma_H(sigma_H) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "k-ev:" << settings.k_ev << "\n";
               o << "r_c:" << settings.r_c << "\n";
               o << "lambda:" << settings.lambda << "\n";
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
     TermProfasiExcludedVolumeBase(ChainFB *chain,
				   std::string name,
				   const Settings &settings=Settings(),
                                   RandomNumberEngine *random_number_engine = &random_global)
	  : EnergyTermCommon(chain, name, settings, random_number_engine) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

	  sigmas.resize(ATOM_ENUM_SIZE);

          // Assign sigma values fo evergy atom
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
     TermProfasiExcludedVolumeBase(const TermProfasiExcludedVolumeBase &other,
                                   RandomNumberEngine *random_number_engine,
			       int thread_index, ChainFB *chain)
	  : EnergyTermCommon(other, random_number_engine, thread_index, chain),
	    sigmas(other.sigmas) {}


     //! Calculate the excluded volume of the given pair of atoms
     //! \param atom_i First atom in the excluded volume interaction 
     //! \param atom_j Second atom in the excluded volume interaction
     //! \return Excluded volume potential between the pair of atoms 
     double pair_energy( Atom * atom_i, Atom * atom_j){

          // Calculate the square of distance
          double r_ij2 = (atom_i->position - atom_j->position).norm_squared();

          double r_cut2 = settings.r_c * settings.lambda * settings.r_c * settings.lambda;

          if (r_ij2 > r_cut2)
               return 0;

          double sigma_i = sigmas[atom_i->atom_type];
          double sigma_j = sigmas[atom_j->atom_type];
          double sig_ij = settings.lambda*(sigma_i + sigma_j);

          // Internal parameters
          double sr2_internal = sig_ij * sig_ij / r_cut2;
          double sr6_internal = sr2_internal * sr2_internal * sr2_internal;
          double sr12_internal = sr6_internal * sr6_internal;

          double asa = -7 * sr12_internal;
          double bsa = 6 * sr12_internal / r_cut2;

          // Internal parameters
          double r2_internal = sig_ij * sig_ij/ r_ij2;
          double r6_internal = r2_internal * r2_internal * r2_internal;

          double energy = asa + bsa*r_ij2 + r6_internal*r6_internal;

          return energy;
     }

     //! Calculate pair energy between atoms if they are on the same residue and seperated by non-constant distance
     //! \param atom_1 First atom
     //! \param atom_2 Second atom
     //! \param fixed_begin Start index of residue1's fixed region (non-hydrogen)
     //! \param fixed_end End index of residue1's fixed region (non-hydrogen)
     //! \param fixed_begin_h Start index of residue1's fixed region (hydrogen)
     //! \param fixed_end_h End index of residue1's fixed region (hydrogen)
     double pair_energy(Atom *atom_1, Atom *atom_2, int fixed_begin, int fixed_end, int fixed_begin_h, int fixed_end_h){

          int at1_in = atom_1->index;
          int at2_in = atom_2->index;

          if ((at1_in >= fixed_begin && at1_in <= fixed_end) || (at1_in >= fixed_begin_h && at1_in <= fixed_end_h)){
               if ((at2_in >= fixed_begin && at2_in <= fixed_end) || (at2_in >= fixed_begin_h && at2_in <= fixed_end_h))
                    return 0;
          }

          return pair_energy(atom_1, atom_2);
     }

     //! Calculate the excluded volume interaction of two atoms if they are separated by non-constant distance and connected by more than three convalent bonds 
     //! \param atom1 First atom in the excluded volume interaction 
     //! \param atom2 Second atom in the excluded volume interaction
     //! \return Excluded volume potential between two atoms 
     inline double calculate_contribution(Atom *atom1, Atom *atom2) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          double energy_ev = 0.0;

          Residue *res1 = atom1->residue;
          Residue *res2 = atom2->residue;

          // Residue distance on the chain
          int dif_res = res2->index - res1->index;

          int ring_fixed_range1_low;
          int ring_fixed_range1_high;
          int ring_fixed_range2_low;
          int ring_fixed_range2_high;

          if (abs(dif_res) > 1 || (chain_distance<ChainFB> (atom1, atom2) > 3)) {

               if (dif_res == 0) {

                    if (res1->residue_type == HIS) {

                         ring_fixed_range1_low = res1->atom_index[CB];
                         ring_fixed_range1_high = res1->atom_index[CE1];
                         ring_fixed_range2_low = res1->atom_index[HD1];
                         ring_fixed_range2_high = res1->atom_index[HE1];

                         if (res1->terminal_status != NTERM) {
                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         } else {
                              ring_fixed_range2_low = ring_fixed_range2_low - 1;
                              ring_fixed_range2_high = ring_fixed_range2_high - 1;

                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         }

                    } else if (res1->residue_type == PHE) {

                         ring_fixed_range1_low = res1->atom_index[CB];
                         ring_fixed_range1_high = res1->atom_index[CZ];
                         ring_fixed_range2_low = res1->atom_index[HD1];
                         ring_fixed_range2_high = res1->atom_index[HZ];

                         if (res1->terminal_status != NTERM) {
                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         } else {
                              ring_fixed_range2_low = ring_fixed_range2_low - 1;
                              ring_fixed_range2_high = ring_fixed_range2_high - 1;

                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         }

                    } else if (res1->residue_type == TRP) {

                         ring_fixed_range1_low = res1->atom_index[CB];
                         ring_fixed_range1_high = res1->atom_index[CH2];
                         ring_fixed_range2_low = res1->atom_index[HD1];
                         ring_fixed_range2_high = res1->atom_index[HZ3];

                         if (res1->terminal_status != NTERM) {
                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         } else {
                              ring_fixed_range2_low = ring_fixed_range2_low - 1;
                              ring_fixed_range2_high = ring_fixed_range2_high - 1;

                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         }

                    } else if (res1->residue_type == ARG) {

                         ring_fixed_range1_low = res1->atom_index[CD];
                         ring_fixed_range1_high = res1->atom_index[NH2];
                         ring_fixed_range2_low = res1->atom_index[HE];
                         ring_fixed_range2_high = res1->atom_index[HH22];

                         if (res1->terminal_status != NTERM) {
                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         } else {
                              ring_fixed_range2_low = ring_fixed_range2_low - 1;
                              ring_fixed_range2_high = ring_fixed_range2_high - 1;

                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         }

                    } else if (res1->residue_type == TYR) {

                         ring_fixed_range1_low = res1->atom_index[CB];
                         ring_fixed_range1_high = res1->atom_index[OH];
                         ring_fixed_range2_low = res1->atom_index[HD1];
                         ring_fixed_range2_high = res1->atom_index[HE2];

                         if (res1->terminal_status != NTERM) {
                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         } else {
                              ring_fixed_range2_low = ring_fixed_range2_low - 1;
                              ring_fixed_range2_high = ring_fixed_range2_high - 1;

                              energy_ev += pair_energy(atom1, atom2,
                                                       ring_fixed_range1_low,
                                                       ring_fixed_range1_high,
                                                       ring_fixed_range2_low,
                                                       ring_fixed_range2_high);
                         }

                    } else if (res1->residue_type == PRO) {

                         if (atom1->atom_type == O || atom2->atom_type == O) {
                              if (res1->terminal_status == INTERNAL) {
                                   energy_ev += pair_energy(atom1, atom2);
                              }
                         }
                         
                         //If PRO is C-terminal
                         if (res1-> terminal_status == CTERM) { 
                              if (atom1->atom_type == OXT || atom2->atom_type == OXT) {
                                   energy_ev += pair_energy(atom1, atom2);
                              }
                         }

                         //If PRO is N-terminal
                         if (res1->terminal_status == NTERM) {
                              if (atom1->atom_type == O) {
                                   if (atom2-> atom_type == H1
                                             || atom2-> atom_type == H2) {
                                        energy_ev += pair_energy(atom1, atom2);
                                   }
                              }

                              if (atom1->atom_type == H1 || atom1-> atom_type == H2) {
                                   if (atom2->atom_type == O) {
                                        energy_ev += pair_energy(atom1, atom2);
                                   }
                              }

                         }
                    } else {
                         energy_ev += pair_energy(atom1, atom2);
                    }

               } else if (dif_res == 1 && res2->residue_type == PRO) {

                    if (atom1->atom_type == O || atom1->atom_type == C
                              || atom1->atom_type == CA) {

                         if (atom2->atom_type == O) {
                              energy_ev += pair_energy(atom1, atom2);
                         }

                         // If residue2 is C-terminal
                         if (res2->terminal_status == CTERM) {
                              if (atom2->atom_type == OXT) {

                                   energy_ev += pair_energy(atom1, atom2);
                              }

                         }

                    }

                    if (atom1->atom_type != O && atom1->atom_type != C
                              && atom1->atom_type != CA) {
                         energy_ev += pair_energy(atom1, atom2);
                    }

               } else if (dif_res == -1 && res1->residue_type == PRO) {

                    if (atom2->atom_type == O || atom2->atom_type == C
                        || atom2->atom_type == CA) {

                         if (atom1->atom_type == O) {
                              energy_ev += pair_energy(atom1, atom2);
                         }
                         
                         // If residue1 is C-terminal
                         if (res1-> terminal_status == CTERM) {
                              if (atom1->atom_type == OXT) {
                                   energy_ev += pair_energy(atom1, atom2);
                              }

                         }

                    }

                    if (atom2->atom_type != O && atom2->atom_type != C
                        && atom2->atom_type != CA) {
                         energy_ev += pair_energy(atom1, atom2);
                    }

               } else {
                    energy_ev += pair_energy(atom1, atom2);
               }

          }
          return energy_ev;
     }
};

//! PROFASI excluded volume energy term
class TermProfasiExcludedVolume: public TermProfasiExcludedVolumeBase<TermProfasiExcludedVolume> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiExcludedVolumeBase<TermProfasiExcludedVolume> TermProfasiExcludedVolumeBase;

public:

     //! Use same settings as base class
     typedef TermProfasiExcludedVolumeBase::Settings Settings;

     //! Settings
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiExcludedVolume(ChainFB *chain,
			       const Settings &settings=Settings(),
                               RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiExcludedVolumeBase(chain, "profasi-excluded-volume", settings, random_number_engine),
            settings(settings){}


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiExcludedVolume(const TermProfasiExcludedVolume &other, 
                               RandomNumberEngine *random_number_engine,
                               int thread_index, ChainFB *chain)
          : TermProfasiExcludedVolumeBase(other, random_number_engine, thread_index, chain),
            settings(other.settings){}


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return excluded volume potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_ev = 0.0;

          // Iterate all the atom pairs on the chain
          for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {
               for(AtomIterator<ChainFB, definitions::ALL> it2(it1+1); !it2.end(); ++it2){

                    Atom *atom1 = &*it1;
                    Atom *atom2 = &*it2;
                    //Residue *res1 = atom1->residue;
                    //Residue *res2 = atom2->residue;

                    energy_ev += this->calculate_contribution(atom1, atom2);
                    //std::cout<<"atom1: "<<atom1->atom_type<<" "<<res1->residue_type<<" "<<res1->index<< " atom2: "<<atom2->atom_type<<" "<<res2->residue_type<<" "<<res2->index<<" "<< this->calculate_contribution(atom1, atom2)<<"\n";
                    //std::cout<<this->calculate_contribution(atom1, atom2)<<"\n";
               }

          }

          energy_ev *= this->settings.k_ev;
          return this->profasi_energy_in_kcal_per_mol * energy_ev;
     }


};


//! PROFASI excluded volume energy term - cached version
class TermProfasiExcludedVolumeCached: public TermProfasiExcludedVolumeBase<TermProfasiExcludedVolumeCached> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiExcludedVolumeBase<TermProfasiExcludedVolumeCached> TermProfasiExcludedVolumeBase;

public:

     //! Cached iterator (used when iterating over node pairs)
     CachedIterator<chaintree::PairIterator<ChainFB,Atom,Atom> > cached_it;

     //! Use same settings as base class
     typedef TermProfasiExcludedVolumeBase::Settings Settings;

     //! Settings
     const Settings settings;

     //! Chaintree iterator settings - described in detail in chaintree.h
     chaintree::PairIterator<ChainFB,Atom,Atom>::Settings iterator_settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiExcludedVolumeCached(ChainFB *chain,
                                     const Settings &settings=Settings(),
                                     RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiExcludedVolumeBase(chain, "profasi-excluded-volume-cached", settings, random_number_engine),
            cached_it(*chain),
            settings(settings){
          
          //! Only evaluate modified pairs (as always when caching)
          bool only_modified_pairs = true;

          iterator_settings = chaintree::PairIterator<ChainFB,Atom,Atom>::Settings(settings.r_cut,
                                                                                   only_modified_pairs);
     }


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiExcludedVolumeCached(const TermProfasiExcludedVolumeCached &other, 
                                     RandomNumberEngine *random_number_engine,
                                     int thread_index, ChainFB *chain)
          : TermProfasiExcludedVolumeBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings),
            iterator_settings(other.iterator_settings) {}


     //! Evaluate chain energy
     //! \param move_info object containing information aout last move
     //! \return excluded volume potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          for (cached_it(*this->chain, iterator_settings); !cached_it.end(); ++cached_it) {

               Atom *atom1 = cached_it->first;
               Atom *atom2 = cached_it->second;
               //Residue *res1 = atom1->residue;
               //Residue *res2 = atom2->residue;

               double contribution = this->calculate_contribution(atom1, atom2);

               cached_it.register_contribution(contribution);
               //std::cout<<"cached:: atom1: "<<atom1->atom_type<<" "<<res1->residue_type<<" "<<res1->index<< " atom2: "<<atom2->atom_type<<" "<<res2->residue_type<<" "<<res2->index<<" "<< this->calculate_contribution(atom1, atom2)<<"\n";
               //std::cout<<this->calculate_contribution(atom1, atom2)<<"\n";
          }

          double energy_ev = cached_it.compute_total();

          energy_ev *= this->settings.k_ev;

          return this->profasi_energy_in_kcal_per_mol * energy_ev;
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
