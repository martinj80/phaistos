// opls_charge.h --- OPLS charge-charge interaction energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson
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

#ifndef OPLS_CHARGE_H
#define OPLS_CHARGE_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "opls_parameters.h"

namespace phaistos {

//! Parameter container for OPLS partial charge parameters
//! This class reads the OPLS parameter file, process the parameters and stores them in a
//! vector for fast lookup
class ChargeParameters: public TinkerParameterBase {

private:

     //! Vector to store parameters in
     std::vector<double> charge;

public:

     //! Constructor
     ChargeParameters() {
          // get third neighbour scale factor from parameter file header
          chg14scale = 1.0;
          std::string chg14scale_str = read_header(&parameters_opls,"chg-14-scale");
          chg14scale /= atof( chg14scale_str.c_str() );

          // get charge unit from parameter file header
          charge_unit = 332.05382; //1 e^2/Å
          std::string dielectric_str = read_header(&parameters_opls,"dielectric");
          charge_unit /= atof( dielectric_str.c_str() );

          // get parameters from 'charge' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls,&raw_param,"charge",1);

          // do not use index zero
          charge.push_back(0.0);

          // store parameters
          unsigned int id(0);
          for (unsigned int i=0; i<raw_param.size(); i++) {
               id = (raw_param[i].first)[0];
               double param = raw_param[i].second[0];
               assert(id == charge.size());
               // if (id != charge.size()) {
               //      std::cerr<<"\nOPLS TERM CHARGE - read error around charge parameter id"
               //               <<id<<"\n\n";
               // }
               charge.push_back(param);
          }
          max_charge_id = id;
     };

     //! Destructor
     ~ChargeParameters() {};

     //! Partial charge getter
     double get(Atom *atom) {
          int id = get_atom_id(atom);
          assert(id <= max_charge_id);
          return ( charge[id] );
     };

     //! Number of parameters
     int max_charge_id;

     //! Scaling factor for third neighbour given in parameter file header
     double chg14scale;

     //! Unit of the energy calculated from the dielectric constant given in parameter file header
     double charge_unit;
};


//! OPLS partial charge interaction term
class TermOplsCharge: public EnergyTermCommon<TermOplsCharge, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsCharge, ChainFB> EnergyTermCommon;     

     //! Partial charge parameters
     //! initialized on construction
     ChargeParameters parameters;

     //! Number of interactions in the last evaluation
     int counter;

     // float fast_inv_sqrt(float x);

public:

     //! Use same settings as energy base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsCharge(ChainFB *chain, 
                    const Settings &settings=Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "opls-charge", settings, random_number_engine) {

          // Annotate chain with biotype
          parameters.annotate_chain(chain);

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsCharge(const TermOplsCharge &other, 
                    RandomNumberEngine *random_number_engine,
                    int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters),
            counter(other.counter){

          // Annotate new chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Evaluate charge interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Charge of atom1
     //! \param chg2 Charge of atom2
     //! \return Charge energy for atom pair
     double calc_charge_pair_energy(Atom *atom1, Atom *atom2, const double chg1, const double chg2) {
          int d = this->parameters.get_covalent_dist<ChainFB>(atom1,atom2);
          if (d > 3) {
               return calc_coulomb_energy(atom1,atom2,chg1,chg2);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3
               return calc_coulomb_energy(atom1,atom2,chg1,chg2) * parameters.chg14scale;
          }
     }
     
     //! Evaluate a coulomb interaction between two atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Charge of atom1
     //! \param chg2 Charge of atom2
     //! \return Coulomb energy for atom pair
     double calc_coulomb_energy(Atom *atom1, Atom *atom2, const double chg1, const double chg2) {

          counter++;

          /* // saves 10s @ 500 iterations (>30% of term calc time!) */
          /* double dx = atom1->position[0] - atom2->position[0]; */
          /* double dy = atom1->position[1] - atom2->position[1]; */
          /* double dz = atom1->position[2] - atom2->position[2]; */
          /* double rSq = dx*dx + dy*dy + dz*dz; */
          /* return ( chg1*chg2*fast_inv_sqrt( (float)rSq ) ); */

          double r = (atom1->position - atom2->position).norm();
          return (chg1*chg2/r);
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
          
          double energySum=0.0;
          counter=0;

          int size = (this->chain)->size();
          for (int i=0; i<size; i++) {
               Residue *res1 = &(*(this->chain))[i];
               int res1_size = res1->size();

               // i==j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    const double chg1 = this->parameters.get(atom1);
                    if (std::fabs(chg1) < 0.001)
                         continue;
                    for (int l=k+1; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         const double chg2 = this->parameters.get(atom2);
                         energySum += calc_charge_pair_energy(atom1,atom2,chg1,chg2);
                    }
               }

               // rest
               for (int j=i+1; j<size; j++) {
                    Residue *res2 = &(*(this->chain))[j];
                    int res2_size = res2->size();

                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const double chg1 = this->parameters.get(atom1);
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const double chg2 = this->parameters.get(atom2);
                              energySum += calc_charge_pair_energy(atom1,atom2,chg1,chg2);
                         }
                    }
               }
          }
          energySum *= this->parameters.charge_unit;
          /* std::cout<<"Charge-Charge  "<<energySum<<" kcal/mol  "<<counter<<" interactions\n"; */
          return energySum;
     }

};


//! OPLS partial charge interaction term -- cached version
class TermOplsChargeCached: public EnergyTermCommon<TermOplsChargeCached, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsChargeCached, ChainFB> EnergyTermCommon;            

     //! Parameters object
     ChargeParameters parameters;

     //! Cached iterator
     CachedIterator<chaintree::PairIterator<ChainFB,Atom,Atom> > cached_it;

     //! Chaintree iterator settigns
     chaintree::PairIterator<ChainFB,Atom,Atom>::Settings iterator_settings;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {

     public:

          // Cutoff distance beyond which energies are 0
          double cutoff_distance;

          Settings(double cutoff_distance = std::numeric_limits<double>::infinity())
               : cutoff_distance(cutoff_distance) {}

          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "cutoff_distance:" << settings.cutoff_distance << "\n";
               o << static_cast<const EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }                    
     } settings;    //!< Local settings object 


     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsChargeCached(ChainFB *chain, const Settings &settings=Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "opls-charge-cached", settings, random_number_engine),
            cached_it(*chain),
            settings(settings) {

          // Only evaluate modified pairs (as always when caching)
          bool only_modified_pairs = true;

          // Annotate chain with biotype
          parameters.annotate_chain(chain);

          iterator_settings = chaintree::PairIterator<ChainFB,Atom,Atom>::Settings(settings.cutoff_distance,
                                                                                   only_modified_pairs);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsChargeCached(const TermOplsChargeCached &other, 
                          RandomNumberEngine *random_number_engine,
                          int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index,chain),
            cached_it(*chain),
            iterator_settings(other.iterator_settings) {

          // Annotate new chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Evaluate charge interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Charge of atom1
     //! \param chg2 Charge of atom2
     //! \param distance Distance between atoms
     //! \return Charge energy for atom pair
     double calc_charge_pair_energy(Atom *atom1, Atom *atom2, const double chg1, const double chg2, double distance) {

          int d = this->parameters.get_covalent_dist<ChainFB>(atom1,atom2);
          // discard atoms pairs less then 3 bonds away
          if (d > 3) {
               // for everybody else do the full calculations
               return (chg1*chg2/distance);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3
               return (chg1*chg2/distance)* this->parameters.chg14scale;
          }
     }


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return partial charge interaction energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          // Cached iterator evaluation
          for (cached_it(*this->chain, this->iterator_settings); !cached_it.end(); ++cached_it) {
               const double chg1 = this->parameters.get(cached_it->first);
               const double chg2 = this->parameters.get(cached_it->second);
               double distance = cached_it->distance;
               // double distance = (cached_it->first->position - cached_it->second->position).norm();
               double contribution = calc_charge_pair_energy(cached_it->first,cached_it->second,chg1,chg2,distance);
               cached_it.register_contribution(contribution);
          }

          double energy = cached_it.compute_total();
          energy *= this->parameters.charge_unit;

          return energy;
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



/*
// saves 10s @ 500 iterations (>30% of term calc time!)
   double dx = atom1->position[0] - atom2->position[0];
   double dy = atom1->position[1] - atom2->position[1];
   double dz = atom1->position[2] - atom2->position[2];
   double rSq = dx*dx + dy*dy + dz*dz;
   return ( chg1*chg2*fast_inv_sqrt( (float)rSq ) );

/ * // Hacky bit-manipulating code that is killed by -O2 and -O3 */
/* // Requires 32 bit IEEE 754 floats */
/* // #pragma requires gcc 4.4 or newer */
/* // from http://en.wikipedia.org/wiki/Fast_inverse_square_root */
/* #pragma optimize("",off) */
/* float TermOplsCharge::fast_inv_sqrt(float x) { */
/*           float xhalf = 0.5f*x; */
/*           int i = *(int*)&x; */
/*           i = 0x5f3759df - (i>>1); */
/*           x = *(float*)&i ; */
/*           // can be added for more accurate estimate */
/*           /\* x = x*(1.5f - xhalf*x*x); *\/ */
/*           return (x*(1.5f - xhalf*x*x)); */
/* } */
/* #pragma optimize("",on) */

}

#endif
