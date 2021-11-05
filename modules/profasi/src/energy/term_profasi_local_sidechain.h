// profasi_local_sidechain.h --- profasi side chain torsion angle energy term
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

// Implementation of the profasi force field
// An effective all-atom potential for proteins, Irback, Mitternacht, Mohanty, PMC Biophysics 2009


#ifndef PROFASI_LOCAL_SIDECHAIN_H
#define PROFASI_LOCAL_SIDECHAIN_H


#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"

namespace phaistos {

//! PROFASI side chain torsion energy term - base class
template <typename DERIVED_CLASS>
class TermProfasiLocalSidechainBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB>,
                                     public EnergyTermCommonProfasi {
private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

     //! (residue, chi index) -> chi class
     //! the kappa and mu vectors are indexed using the chi class
     std::vector<std::vector<int> > chi_class_map;

public:
     //! Local setting class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Kappa_local 3
          std::vector<double> kappa;
          std::vector<int> n;

          //! Constructor. Defines default values for settings object
          Settings(std::vector<double> kappa = vector_utils::make_vector(0.6, 0.3, 0.4, -0.4),
		   std::vector<int> n = vector_utils::make_vector(3, 3, 2, 2))
	       : kappa(kappa),
		 n(n) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "kappa:" << settings.kappa << "\n";
               o << "n:" << settings.n << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }          

     } settings;


     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object     
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiLocalSidechainBase(ChainFB *chain,
                                   std::string name,
                                   const Settings &settings=Settings(),
                                   RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine),
            chi_class_map(chain->size()) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          //! Import phaistos::make_vector namespace
          using namespace vector_utils;

          // Set up chi_class map
          for (unsigned int i=0; i < chi_class_map.size(); ++i) {

               ResidueEnum residue_type = (*chain)(i).residue_type;

               switch (residue_type) {
               case SER:
               case CYS:
               case THR:
               case VAL:
                    chi_class_map[i] = make_vector(1);
                    break;
               case ILE:
               case LEU:
                    chi_class_map[i] = make_vector(1,1);
                    break;
               case ASP:
               case ASN:
                    chi_class_map[i] = make_vector(1,4);
                    break;
               case HIS:
               case PHE:
               case TYR:
               case TRP:
                    chi_class_map[i] = make_vector(1,3);
                    break;
               case MET:
                    chi_class_map[i] = make_vector(1,1,2);
                    break;
               case GLU:
               case GLN:
                    chi_class_map[i] = make_vector(1,1,4);
                    break;
               case LYS:
                    chi_class_map[i] = make_vector(1,1,1,1);
                    break;
               case ARG:
                    chi_class_map[i] = make_vector(1,1,1,3);
                    break;
               default:
                    break;
               }

               // PRO not included as Profasi degree of freedom
               unsigned int chi_counter = 0;
               if (residue_type != PRO) {
                    for (unsigned int j=0; j<(*chain)(i).chi_atoms.size(); ++j) {
                         bool is_hydrogen = atom_type_XH[(*chain)(i).chi_atoms[j].first->atom_type];
                         if (!is_hydrogen) {
                              chi_counter++;
                         }
                    }
               }
               assert(chi_counter == chi_class_map[i].size());
          }
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiLocalSidechainBase(const TermProfasiLocalSidechainBase &other, 
                                   RandomNumberEngine *random_number_engine,
                                   int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            chi_class_map(other.chi_class_map) {
     }


     //! Evaluate local sidechain potential energy
     //! \param &it Current residue iterator information
     //! \return local sidechain potential energy of current residue
     double calculate_contribution(ResidueIterator<ChainFB> &it){

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          double E_sctorsion = 0.0;
          double psi = 0.0;
          double kgly = -0.15;
          double Gly_eng = 0.0;

          // Energy penalty for glycine psi values around +- 120 degree
          if (it->residue_type == GLY) {
               if (it->terminal_status == INTERNAL) {

                    psi = (*it)[C]->get_dihedral();
                    Gly_eng += kgly*(cos(psi) + 2*cos(2*psi));
                    //std::cout << (*it)[C]->get_dihedral() << "\n";
                    //std::cout << "Gly_eng: " << it->index << " = " << Gly_eng << "\n";
               }
          }

          double sum_inner = 0.0;
          
          const std::vector<int> &chi_classes = chi_class_map[it->index];

          for(unsigned int i=0; i < chi_classes.size(); i++) {

               double chi = it->chi_atoms[i].first->get_dihedral();
               double kappa = settings.kappa[chi_classes[i]-1];
               double n = settings.n[chi_classes[i]-1];

               // In order to get same result with profasi
               if(it->residue_type == ARG && i == 3){ 

                    Atom *atom_1 = (*it)[CG];
                    Atom *atom_2 = (*it)[CD];
                    Atom *atom_3 = (*it)[NE];
                    Atom *atom_4 = (*it)[HE];

                    chi = calc_dihedral(atom_1->position, atom_2->position, atom_3->position, atom_4->position );
               }

               sum_inner += kappa*cos(n*chi);
          }

          E_sctorsion = sum_inner + Gly_eng;

          return E_sctorsion;
     }

};


//! PROFASI local sidechain energy term
class TermProfasiLocalSidechain: public TermProfasiLocalSidechainBase<TermProfasiLocalSidechain> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiLocalSidechainBase<TermProfasiLocalSidechain> TermProfasiLocalSidechainBase;

public:

     //! Use same settings as base class
     typedef TermProfasiLocalSidechainBase::Settings Settings;

     //! Settings
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object     
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiLocalSidechain(ChainFB *chain,
                               const Settings &settings=Settings(),
                               RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiLocalSidechainBase(chain, "profasi-local-sidechain", settings, random_number_engine),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiLocalSidechain(const TermProfasiLocalSidechain &other, 
                               RandomNumberEngine *random_number_engine,
                               int thread_index, ChainFB *chain)
          : TermProfasiLocalSidechainBase(other, random_number_engine, thread_index, chain),
            settings(other.settings) {
     }

     //! Evaluate chain energy (normal version)
     //! \param move_info object containing information about last move
     //! \return Local sidechain potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy = 0.0;

          // Go through all the residues
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){

               energy += this->calculate_contribution(it);
          }

          return this->profasi_energy_in_kcal_per_mol * energy;
     }

};


//! PROFASI local local sidechain energy term - cached version
class TermProfasiLocalSidechainCached: public TermProfasiLocalSidechainBase<TermProfasiLocalSidechainCached> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiLocalSidechainBase<TermProfasiLocalSidechainCached> TermProfasiLocalSidechainBase;

public:

     //! Cached iterator
     CachedIterator<ResidueIterator<ChainFB> > cached_it;

     //! Use same settings as base class
     typedef TermProfasiLocalSidechainBase::Settings Settings;

     //! Local Settings object
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiLocalSidechainCached(ChainFB *chain,
                                     const Settings &settings=Settings(),
                                     RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiLocalSidechainBase(chain, "profasi-local-sidechain-cached", settings, random_number_engine),
            cached_it(*chain),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiLocalSidechainCached(const TermProfasiLocalSidechainCached &other, 
                                     RandomNumberEngine *random_number_engine,
                                     int thread_index, ChainFB *chain)
          : TermProfasiLocalSidechainBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings) {
     }

     //! Evaluate chain energy (cached version)
     //! \param move_info object containing information about last move
     //! \return Local sidechain potential energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          int start_index = 0;
          int end_index = this->chain->size();

          if (move_info) {
               start_index = std::max(0,move_info->modified_angles_start);
               end_index = std::min(this->chain->size(), move_info->modified_angles_end);
          }
          for (cached_it(*this->chain,
                              start_index,
                              end_index); !cached_it.end(); ++cached_it) {

               double contribution = this->calculate_contribution(cached_it);

               cached_it.register_contribution(contribution);
          }

          double energy = cached_it.compute_total();

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
