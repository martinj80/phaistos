// profasi_local.h --- PROFASI local energy term E_loc1: a local electrostatic
//                     interaction between adjacent peptide-bond-units(NH and CO) and 
//                     E_loc2: repulsion between  neighbouring backbone H-H and O-O.
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


#ifndef PROFASI_LOCAL_H
#define PROFASI_LOCAL_H

#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"
#include <cmath>

namespace phaistos {

//! PROFASI local energy term - base class
template<typename DERIVED_CLASS>
class TermProfasiLocalBase: public EnergyTermCommon< DERIVED_CLASS, ChainFB>,
                            public EnergyTermCommonProfasi {
private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

public:

     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Kappa_local1 (overall scaling)
          double k1_loc;

          //! Kappa_local2 (overall scaling)
          double k2_loc;

          //! Partial charge of Carbon
          double q_C;

          //! Partial charge of Oxygen
          double q_O;

          //! Partial charge of Nitrogen
          double q_N;

          //! Partial charge of Hydrogen
          double q_H;

          //! Constructor. Defines default values for settings object
          Settings(double k1_loc = 6.0,
		   double k2_loc = 1.2,
		   double q_C    = 0.42,
		   double q_O    = -0.42,
		   double q_N    = -0.2,
		   double q_H    = 0.2)
	        : k1_loc(k1_loc),
	          k2_loc(k2_loc),
                  q_C(q_C),
                  q_O(q_O),
                  q_N(q_N),
                  q_H(q_H) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "k1_loc:" << settings.k1_loc << "\n";
               o << "k2_loc:" << settings.k2_loc << "\n";
               o << "q-C:" << settings.q_C << "\n";
               o << "q-O:" << settings.q_O << "\n";
               o << "q-N:" << settings.q_N << "\n";
               o << "q-H:" << settings.q_H << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }          

     } settings;    //!< Local settings object 


     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiLocalBase(ChainFB *chain,
                          std::string name,
                          const Settings &settings=Settings(),
                          RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine) {}


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiLocalBase(const TermProfasiLocalBase &other,
                          RandomNumberEngine *random_number_engine,
                          int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {}


     //! E_loc2: OO and HH repulsion for neighboring peptide.
     //! \param it Residue iterator
     //! \return Repulsive potential energy
     double repulsion(ResidueIterator<ChainFB> &it) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          double f_vi = 0.0;
          double f_ui = 0.0;
          ResidueIterator<ChainFB> res1 = it;
          ResidueIterator<ChainFB> res2 = it + 1;
          ResidueIterator<ChainFB> res3 = it + 2;

          if (res2->residue_type != GLY) {

               Atom *atom_C1 = (*res1)[C];
               Atom *atom_C2 = (*res2)[C];
               Atom *atom_O1 = (*res1)[O];
               Atom *atom_O2 = (*res2)[O];

               double v_i = std::min(
                         (atom_O1->position - atom_C2->position).norm(),
                         (atom_C1->position - atom_O2->position).norm())
                         - (atom_O1->position - atom_O2->position).norm();
               f_vi = std::max(0.0, tanh(3* v_i ));

               if (res2->residue_type != PRO && res3->residue_type != PRO) {
                    Atom *atom_N2 = (*res2)[N];
                    Atom *atom_N3 = (*res3)[N];
                    Atom *atom_H2 = (*res2)[H];
                    Atom *atom_H3 = (*res3)[H];

                    double u_i = std::min((atom_H2->position
                              - atom_N3->position).norm(), (atom_N2->position
                              - atom_H3->position).norm()) - (atom_H2->position
                              - atom_H3->position).norm();
                    f_ui = std::max(0.0, tanh(3* u_i ));
               }
          }

          return (f_ui + f_vi);
     }

     //! E_loc1: Coulomb interactions between partial charges of neighboring peptide units
     //! \param it Residue iterator
     //! \return Coulomb potential energy
     double coulumb(ResidueIterator<ChainFB> &it){

          // Import protein definitions (such as residue names)
          using namespace definitions;

          double e = 0.0;

          ResidueIterator<ChainFB> res1 = it;
          ResidueIterator<ChainFB> res2 = it + 1;
          ResidueIterator<ChainFB> res3 = it + 2;

          Atom *atom_C1 = (*res1)[C];
          Atom *atom_C2 = (*res2)[C];
          Atom *atom_O1 = (*res1)[O];
          Atom *atom_O2 = (*res2)[O];

          double eng_C1C2 = settings.q_C * settings.q_C / (atom_C1->position - atom_C2->position).norm();
          double eng_C1O2 = settings.q_C * settings.q_O / (atom_C1->position - atom_O2->position).norm();
          double eng_O1C2 = settings.q_O * settings.q_C / (atom_O1->position - atom_C2->position).norm();
          double eng_O1O2 = settings.q_O * settings.q_O / (atom_O1->position - atom_O2->position).norm();

          e += (eng_C1C2 + eng_C1O2 + eng_O1C2 + eng_O1O2);

          if(res3->residue_type != PRO){

               Atom *atom_N3 = (*res3)[N];
               Atom *atom_H3 = (*res3)[H];

               double eng_C1N3 = settings.q_C * settings.q_N / (atom_C1->position - atom_N3->position).norm();
               double eng_C1H3 = settings.q_C * settings.q_H / (atom_C1->position - atom_H3->position).norm();
               double eng_O1N3 = settings.q_O * settings.q_N / (atom_O1->position - atom_N3->position).norm();
               double eng_O1H3 = settings.q_O * settings.q_H / (atom_O1->position - atom_H3->position).norm();

               e += (eng_C1N3 + eng_C1H3 + eng_O1N3 + eng_O1H3);

               if(res2->residue_type != PRO){

                    Atom *atom_N2 = (*res2)[N];
                    Atom *atom_H2 = (*res2)[H];

                    double eng_N2C2 = settings.q_N * settings.q_C / (atom_N2->position - atom_C2->position).norm();
                    double eng_N2O2 = settings.q_N * settings.q_O / (atom_N2->position - atom_O2->position).norm();
                    double eng_H2C2 = settings.q_H * settings.q_C / (atom_H2->position - atom_C2->position).norm();
                    double eng_H2O2 = settings.q_H * settings.q_O / (atom_H2->position - atom_O2->position).norm();
                    double eng_N2N3 = settings.q_N * settings.q_N / (atom_N2->position - atom_N3->position).norm();
                    double eng_N2H3 = settings.q_N * settings.q_H / (atom_N2->position - atom_H3->position).norm();
                    double eng_H2N3 = settings.q_H * settings.q_N / (atom_H2->position - atom_N3->position).norm();
                    double eng_H2H3 = settings.q_H * settings.q_H / (atom_H2->position - atom_H3->position).norm();

                    e += (eng_N2C2 + eng_N2O2 + eng_H2C2 + eng_H2O2 + eng_N2N3 + eng_N2H3 + eng_H2N3 + eng_H2H3);
               }
          }

          if(res3->residue_type == PRO && res2->residue_type != PRO){

               Atom *atom_N2 = (*res2)[N];
               Atom *atom_H2 = (*res2)[H];

               double eng_N2C2 = settings.q_N * settings.q_C / (atom_N2->position - atom_C2->position).norm();
               double eng_N2O2 = settings.q_N * settings.q_O / (atom_N2->position - atom_O2->position).norm();
               double eng_H2C2 = settings.q_H * settings.q_C / (atom_H2->position - atom_C2->position).norm();
               double eng_H2O2 = settings.q_H * settings.q_O / (atom_H2->position - atom_O2->position).norm();

               e += (eng_N2C2 + eng_N2O2 + eng_H2C2 + eng_H2O2);
          }

          return e;
     }

};




//! PROFASI local energy term
class TermProfasiLocal: public TermProfasiLocalBase<TermProfasiLocal> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiLocalBase<TermProfasiLocal> TermProfasiLocalBase;

public:

     // Use same settings as base class
     typedef TermProfasiLocalBase::Settings Settings;

     //! Local Settings object
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiLocal(ChainFB *chain,
                      const Settings &settings=Settings(),
                      RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiLocalBase(chain, "profasi-local", settings, random_number_engine),
            settings(settings){
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiLocal(const TermProfasiLocal &other, 
                      RandomNumberEngine *random_number_engine,
                      int thread_index, ChainFB *chain)
          : TermProfasiLocalBase(other, random_number_engine, thread_index, chain),
            settings(other.settings){
     }

     //! Evaluate chain energy
     //! \param move_info Object containing information about the last executed move
     //! \return Local energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          double energy = 0.0;

          // Go through residues
          for (ResidueIterator<ChainFB> it(*this->chain, 0, this->chain->size()-2); !it.end(); ++it) {

               double coulumb_val = this->settings.k1_loc * this->coulumb(it);
               double repulsion_val = this->settings.k2_loc * this->repulsion(it);

               energy += (coulumb_val + repulsion_val);
          }

          return this->profasi_energy_in_kcal_per_mol * energy;
     }

};


//! PROFASI local energy term - cached version
class TermProfasiLocalCached: public TermProfasiLocalBase<TermProfasiLocalCached> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiLocalBase<TermProfasiLocalCached> TermProfasiLocalBase;

public:

     //! Cached iterator
     CachedIterator<ResidueIterator<ChainFB> > cached_it;

     //! Use same settings as base class
     typedef TermProfasiLocalBase::Settings Settings;

     //! Local Settings object
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiLocalCached(ChainFB *chain,
                            const Settings &settings=Settings(),
                            RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiLocalBase(chain, "profasi-local-cached", settings, random_number_engine),
            cached_it(*chain),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiLocalCached(const TermProfasiLocalCached &other, 
                            RandomNumberEngine *random_number_engine,
                            int thread_index, ChainFB *chain)
          : TermProfasiLocalBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings) {
     }

     //! Evaluate chain energy(cached version)
     //! \param move_info Object containing information about last move
     //! \return Local energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          int start_index_full_range = 0;
          int end_index_full_range = this->chain->size()-2;

          int start_index = start_index_full_range;
          int end_index = end_index_full_range;

          if (move_info) {
               start_index = std::max(0,move_info->modified_angles_start-2);
               end_index = std::min(this->chain->size()-2, move_info->modified_angles_end);
          }

          for (cached_it(*this->chain,
                         start_index,
                         end_index,
                         start_index_full_range,
                         end_index_full_range); !cached_it.end(); ++cached_it) {

               double coulumb_val = this->settings.k1_loc * this->coulumb(cached_it);
               double repulsion_val = this->settings.k2_loc * this->repulsion(cached_it);

               cached_it.register_contribution(coulumb_val + repulsion_val);
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
