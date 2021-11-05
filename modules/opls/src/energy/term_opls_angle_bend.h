// opls_angle_bend.h --- OPLS angle-bend energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson, Wouter Boomsma
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

#ifndef OPLS_ANGLEBEND_H
#define OPLS_ANGLEBEND_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "energy/parameter_data_structures.h"
#include "energy/opls_parameters.h"

namespace phaistos {

//! Parameter container for OPLS anglebend energy term
class AnglebendParameters: public TinkerParameterBase {

public:

     //! Container class
     class Parameter {
          
     public:
          
          //! Equilibrium angle
          double eq_ang;

          //! Sprint constant
          double k;

          //! Boolean interpretation
          operator bool() {
               return (k > -0.0001);
          };

          //! Default constructor
          Parameter() : k(-1.0) {};

          // Main constructor
          Parameter(double k, double eq_ang) : eq_ang(eq_ang), k(k) {};
     };

private:

     //! Map to store parameters in
     ParameterData3D<Parameter> parameter_map;

public:
     
     //! Constructor
     AnglebendParameters() {

          // get parameters from 'angle' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls, &raw_param, "angle", 3);

          // store parameters
          for (unsigned int i = 0; i < raw_param.size(); i++) {
               std::vector<int> id = (raw_param[i].first);
               std::vector<double> param = raw_param[i].second;
               // check if parameter have already been defined
               assert( ! parameter_map(id[0], id[1], id[2]) );
               // if (parameter_map(id[0], id[1], id[2])) {
               //      std::cerr << "\nOPLS TERM ANGLEBEND - multiple definitions of parameter " << id[0] << "," << id[1] << "," << id[2] << "\n\n";
               // }
               // tinker use degrees - cute
               param[1] = deg2rad(param[1]);
               parameter_map(id[0], id[1], id[2]) = Parameter(param[0], param[1]);
          }
     }
     ;

     //! Destructor
     ~AnglebendParameters() {};

     //! Getter
     Parameter get(Atom *atom1, Atom *atom2, Atom *atom3) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          int id3 = get_param_id(atom3);
          if (id1 > id3) {
               int buf = id1;
               id1 = id3;
               id3 = buf;
          }
          Parameter param = parameter_map(id1, id2, id3);
          assert( param );
          // if (!param)
          //      std::cerr << "\nOPLS TERM ANGLEBEND - missing parameter for param id " << id1 << "," << id2 << ","
          //                << id3 << "\n\n";          
          return param;
     }
};

//! OPLS anglebend energy term - base class containing all functionality
template<typename DERIVED_CLASS>
class TermOplsAngleBendBase: public EnergyTermCommon<DERIVED_CLASS,ChainFB> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

     //! Parameters
     AnglebendParameters parameters;

     //! Number of interactions in lastevaluation
     int counter;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Whether to exclude sidechain interactions
          bool omit_sidechains;

          //! Constructor. Defines default values for settings object.
          Settings(bool omit_sidechains = false)
               : omit_sidechains(omit_sidechains) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "omit-sidechains:" << settings.omit_sidechains << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }                    
     } settings;  //!< Local settings object 

     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsAngleBendBase(ChainFB *chain, std::string name, const Settings &settings = Settings(),
                           RandomNumberEngine *random_number_engine = &random_global) :
          EnergyTermCommon(chain, name, settings, random_number_engine),
          settings(settings) {

          // Annotate the chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsAngleBendBase(const TermOplsAngleBendBase &other, 
                           RandomNumberEngine *random_number_engine,
                           int thread_index, ChainFB *chain) :
          EnergyTermCommon(other, random_number_engine, thread_index, chain), 
          parameters(other.parameters), counter(other.counter),
          settings(other.settings) {

          // Annotate the new chain with biotype
          parameters.annotate_chain(chain);
     }


     //! Evaluate anglebend energy for an atom
     //! \param atom2 Central atom
     //! \return Angle bend energy between all neighboring atoms
     inline double calc_anglebend_energy(Atom *atom2) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          double angle, energy = 0.0;
          CovalentBondIterator<ChainFB> it1(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
          for (; !it1.end(); ++it1) {
               Atom *atom1 = &*it1;
               CovalentBondIterator<ChainFB> it3(it1);
               //discard it3 = it1
               ++it3;
               for (; !it3.end(); ++it3) {
                    Atom *atom3 = &*it3;
                    if (atom1->atom_type == PS || atom2->atom_type == PS || atom3->atom_type == PS )
                         continue;
                    angle = calc_angle(atom1->position, atom2->position, atom3->position);

                    AnglebendParameters::Parameter param = parameters.get(atom1, atom2, atom3);

                    energy += calc_spring_energy(angle, param.eq_ang, param.k);
               }
          }
          return energy;
     }

     //! Evaluate the potetial energy in a spring
     //! \param x Angle
     //! \param x_eq Equilibrium angle
     //! \param k Spring constant
     //! \return Spring energy     
     inline double calc_spring_energy(double x, double x_eq, double k) {
          counter++;
          double dx = x - x_eq;
          /* if (dx > M_PI || dx < -M_PI) */
          /*      std::cout<<"dx = "<<dx<<"\n"; */
          return (k * dx * dx);
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return angle bend potential energy of the chain in the object
     double evaluate(MoveInfo *move_info = NULL) {

          double energy_sum = 0.0;
          counter = 0;

          // all angles
          int size = (this->chain)->size();
          for (int r = 0; r < size; r++) {
               Residue *res = &(*(this->chain))[r];
               int res_size = res->size();
               for (int a = 0; a < res_size; a++) {
                    Atom *atom = res->atoms[a];

                    // Atom must have more than one neighbour to contain angles
                    if (parameters.get_bond_n(atom) < 2)
                         continue;

                    if (settings.omit_sidechains && atom->is_sidechain_atom)
                         continue;

                    energy_sum += calc_anglebend_energy(atom);
               }
          }

          /* std::cout<<"Angle Bending  "<<energy_sum<<" kcal/mol  "<<counter<<" interactions\n"; */
          return energy_sum;
     }

};


//! OPLS anglebend energy term
class TermOplsAngleBend: public TermOplsAngleBendBase<TermOplsAngleBend> {

public:

     //! For convenience, define local TermOplsAngleBendBase
     typedef phaistos::TermOplsAngleBendBase<TermOplsAngleBend> TermOplsAngleBendBase;            

     //! Reuse base-class Settings object
     typedef TermOplsAngleBendBase::Settings Settings;

     //! Constructor
     TermOplsAngleBend(ChainFB *chain, const Settings &settings = Settings(),
                       RandomNumberEngine *random_number_engine = &random_global) 
          : TermOplsAngleBendBase(chain, "opls-angle-bend", settings, random_number_engine) {}

     //! Copy constructor
     TermOplsAngleBend(const TermOplsAngleBend &other, 
                       RandomNumberEngine *random_number_engine,
                       int thread_index, ChainFB *chain) 
          : TermOplsAngleBendBase(other, random_number_engine, thread_index, chain) {
     }
};


//! OPLS anglebend energy term -- Cached version
class TermOplsAngleBendCached: public TermOplsAngleBendBase<TermOplsAngleBendCached> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::TermOplsAngleBendBase<TermOplsAngleBendCached> TermOplsAngleBendBase;

     //! Cached iterator
     CachedIterator<ResidueIterator<ChainFB> > cached_it;

     //! Parameters
     AnglebendParameters parameters;
     
public:

     //! For convenience, define local standard Settings
     typedef TermOplsAngleBendBase::Settings Settings;

     //! Settings
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsAngleBendCached(ChainFB *chain, 
                             const Settings &settings = Settings(),
                             RandomNumberEngine *random_number_engine = &random_global)
          : TermOplsAngleBendBase(chain, "opls-angle-bend-cached", settings),
            cached_it(*chain),
            settings(settings) {
     }
    
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsAngleBendCached(const TermOplsAngleBendCached &other, 
                             RandomNumberEngine *random_number_engine,
                             int thread_index, ChainFB *chain)
          : TermOplsAngleBendBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings) {}

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return angle bend potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          if (move_info) {
               // Cached iterator evaluation
               for (cached_it(*this->chain, 
                              std::max(0,move_info->modified_angles_start-1), 
                              std::min(this->chain->size(), move_info->modified_angles_end+1)); !cached_it.end(); ++cached_it) {
                    
                    ChainFB::Residue &res = *cached_it;

                    int res_size = res.size();
                    double residue_sum = 0.0;
                    for (int a = 0; a < res_size; a++) {
                         Atom *atom = res.atoms[a];

                         // Atom must have more than one neighbour to contain angles
                         if (parameters.get_bond_n(atom) < 2)
                              continue;

                         if (settings.omit_sidechains && atom->is_sidechain_atom)
                              continue;

                         residue_sum += this->calc_anglebend_energy(atom);
                    }
                    cached_it.register_contribution(residue_sum);
               }
               return cached_it.compute_total();
          } else {
               return TermOplsAngleBendBase::evaluate(move_info);
          }
     }     

     //! Accept last evaluation
     void accept() {
          TermOplsAngleBendBase::accept();
          cached_it.accept();
     }

     //! Reject last evaluation
     void reject() {
          TermOplsAngleBendBase::reject();
          cached_it.reject();
     }
};

}

#endif
