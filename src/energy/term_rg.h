// rg.h --- Radius of gyration energy term
// Copyright (C) 2008-2011 Wouter Boomsma
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


#ifndef RG_H
#define RG_H

#include "../energy/energy_term.h"
#include "protein/iterators/atom_iterator.h"
#include "observable_base.h"

namespace phaistos {

//! Calculate radius of gyration within range specified by iterators.
//! \param begin Iterator start-point
//! \param end Iterator end-point
template <typename ITERATORTYPE>     
double calc_rg_from_iterators(ITERATORTYPE begin, ITERATORTYPE end) {
     double res = 0.0;
     Vector_3D CM = center_of_mass(begin, end);
     int counter = 0;
     for (ITERATORTYPE it1=begin; it1 != end; ++it1) {
	  res += (*it1-CM)*(*it1-CM);
	  counter++;
     }
     res=sqrt(res/(double)counter);

     return res;
}

//! Calculate radius of gyration within range specified by iterator
//! \param it Iterator (iterator contains both start and end information)
template <typename ITERATORTYPE>     
double calc_rg_from_iterators(ITERATORTYPE it) {
     double res = 0.0;
     Vector_3D CM = center_of_mass(it);
     int counter = 0;
     for (; !it.end(); ++it) {
	  res += (*it-CM)*(*it-CM);
	  counter++;
     }
     res=sqrt(res/(double)counter);

     return res;
}

//! Calculate radius of gyration
//! \param chain Molecule chain
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE>
double calc_rg(const CHAIN_TYPE &chain) {
     return calc_rg_from_iterators(AtomIterator<CHAIN_TYPE,ITERATION_MODE,Vector_3D>(chain));
}


//! Radius of gyration energy term
template <typename CHAIN_TYPE>
class TermRg: public EnergyTermCommon<TermRg<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermRg<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;               

     //! Function pointer used to select either CA or ALL iteration mode
     double (*energy_function)(const CHAIN_TYPE &chain);
     
public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! Offset. Used to compare with given rg, for instance of a native state.
          double offset_value;

          //! Whether to calculate CA-only radius of gyration
          bool ca_only;

          //! Constructor. Defines default values for settings object.
          Settings(double offset_value=0.0,
                   bool ca_only=true )
               : offset_value(offset_value),
                 ca_only(ca_only) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "offset-value:" << settings.offset_value << "\n";
               o << "ca-only:" << settings.ca_only << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }                    
     } settings;  //!< Local settings object


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermRg(CHAIN_TYPE *chain, const Settings &settings=Settings(),
            RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "rg", settings, random_number_engine),
            settings(settings) {

	  if (settings.ca_only)
	       energy_function = &calc_rg<CHAIN_TYPE,definitions::CA_ONLY>;
	  else
	       energy_function = &calc_rg<CHAIN_TYPE,definitions::ALL>;
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermRg(const TermRg &other, 
            RandomNumberEngine *random_number_engine,
            int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            energy_function(other.energy_function),
            settings(other.settings) {}

     
     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {
	  return (*energy_function)(*(this->chain))-settings.offset_value;
     }
};

}

#endif
