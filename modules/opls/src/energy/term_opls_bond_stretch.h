// opls_bond_stretch.h --- OPLS bond-stretch energy term
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

#ifndef OPLS_BONDSTRETCH_H
#define OPLS_BONDSTRETCH_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "energy/parameter_data_structures.h"
#include "energy/opls_parameters.h"

namespace phaistos {
//! Parameter container for OPLS bondstretch energy term
class BondstretchParameters: public TinkerParameterBase {

public:

     //! Container class
     class Parameter {

     public:

          //! Bond equilibrium length
          double eq_len;

          //! Bond spring constant
          double k;

          //! Boolean interpretation
          operator bool() {
               return (k > -0.0001);
          };
          
          //! Default constructor
          Parameter():k(-1.0) {};
          
          //! Main constructor
          Parameter(double k,double eq_len): eq_len(eq_len), k(k) {};
     };
     
private:
     
     //! Map to store parameters in
     ParameterData2D<Parameter> parameter_map;     
     
public:
     
     //! Constructor
     BondstretchParameters() {

          // get parameters from 'angle' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls,&raw_param,"bond",2);

          // store parameters
          for (unsigned int i=0; i<raw_param.size(); i++) {
               std::vector<int> id = (raw_param[i].first);
               std::vector<double> param = raw_param[i].second;
               // check for multiple definitions of the same parameter in the parameter file
               assert( ! parameter_map(id[0],id[1]) );
               // if ( parameter_map(id[0],id[1]) ) {
               //      std::cerr<<"\nOPLS TERM BONDSTRETCH - multiple definitions of parameter "
               //               <<id[0]<<","<<id[1]<<"\n\n";
               // }
               parameter_map(id[0],id[1]) = Parameter(param[0],param[1]);
          }
     };

     //! Destructor
     ~BondstretchParameters() {};

     //! Getter
     Parameter get(Atom *atom1, Atom *atom2) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          if (id1 > id2) {
               int buf=id1; id1=id2; id2=buf;
          }
          Parameter param = parameter_map(id1,id2);
          assert( param );
          // if (!param)
          //      std::cerr<<"\nOPLS TERM BONDSTRETCH - missing parameter for param id "
          //                    <<id1<<","<<id2<<"\n\n";
          return param;
     };
};


//! OPLS bondstretch energy term
class TermOplsBondStretch: public EnergyTermCommon<TermOplsBondStretch, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsBondStretch, ChainFB> EnergyTermCommon;

     //! Parameters
     BondstretchParameters parameters;

     //! Number of interactions calculated
     int counter;
     
public:

     // Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;     

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsBondStretch(ChainFB *chain, 
                         const Settings &settings=Settings(),
                         RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "opls-bond-stretch", settings, random_number_engine) {

          // Annotate the chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsBondStretch(const TermOplsBondStretch &other, 
                         RandomNumberEngine *random_number_engine,
                         int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters),
            counter(other.counter) {

          // Annotate the new chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Evaluate anglebend energy for a bond (atom1-atom2)
     //! \param atom2 Second atom
     //! \return Bond stretch energy to previous neighbours
     inline double calc_bondstretch_energy(Atom *atom2) {
          double energy = 0.0;
          
          CovalentBondIterator<ChainFB> it1(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
          for (; !it1.end(); ++it1) {
               Atom *atom1 = &*it1;
               
               // Only calculate the distance to covalently bonded that are before atom2 in the chain
               // (otherwise all distance would be considered twice)
               if (atom1->residue->index < atom2->residue->index ||
                   (atom1->residue->index == atom2->residue->index && atom1->index < atom2->index)) {
                    BondstretchParameters::Parameter param = parameters.get(atom1,atom2);
                    double length = (atom1->position - atom2->position).norm();
                    energy += calc_spring_energy(length, param.eq_len, param.k);
               }
          }
          
          return energy;
     }
     
     //! Evaluate the potetial energy in a spring
     //! \param x Distance
     //! \param x_eq Equilibrium distance
     //! \param k Spring constant
     //! \return Spring energy     
     inline double calc_spring_energy(double x, double x_eq, double k) {
          counter++;
          double dx = x - x_eq;
          return (k*dx*dx);
     } 
     
     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return bond stretch potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energySum=0.0;
          counter=0;
          
          //+1 to avoid getNeighbout(-1) error from N[0]
          AtomIterator<ChainFB,definitions::ALL> it(*(this->chain)); ++it;
          for (; !it.end(); ++it) {
               Atom &atom = *it;
               energySum += calc_bondstretch_energy(&atom);
          }

          /* int size = chain->size(); */
          /* Residue *res = &(*chain)[0]; */
          /* int res_size = res->size(); */
          /* for (int a=1; a<res_size; a++) { */
          /*      Atom *atom = res->atoms[a]; */
          /*      energySum += calc_torsion_energy(atom); */
          /* } */
          /* for (int r=1; r<size; r++) { */
          /*      res = &(*chain)[r]; */
          /*      res_size = res->size(); */
          /*      for (int a=0; a<res_size; a++) { */
          /*           Atom *atom = res->atoms[a]; */
          /*           energySum += calc_torsion_energy(atom); */
          /*      } */
          /* } */
          
          /* std::cout<<"Bond Stretching  "<<energySum<<" kcal/mol  "<<counter<<" interactions\n"; */
          return energySum;
     }
     
};

}

#endif
