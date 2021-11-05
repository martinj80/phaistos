// observable.h --- Observable. Wrapper for a Phaistos Energy term.
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

#ifndef OBSERVABLE_H
#define OBSERVABLE_H

namespace phaistos {

#include "observable_base.h"

//! Observable. Wrapper for a Phaistos energy term with settings allowing
//! the user to specify how the observable should be reported.
template <typename ENERGY_TERM>
class Observable: public ENERGY_TERM, public ObservableBase {
public:

     //! Local settings class.
     const class Settings: public ENERGY_TERM::Settings, public ObservableBase::Settings {
     public:

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<typename ENERGY_TERM::Settings>(settings);
               o << static_cast<ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     

     //! Constructor.
     //! \param energy_term energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const ENERGY_TERM &energy_term, 
                const ObservableBase::Settings &settings=ObservableBase::Settings(), 
                Energy<typename ENERGY_TERM::ChainType> *reference_energy_function=NULL)
          : ENERGY_TERM(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const Observable &other, 
                RandomNumberEngine *random_number_engine,
                int thread_index, typename ENERGY_TERM::ChainType *chain, 
                Energy<typename ENERGY_TERM::ChainType> *reference_energy_function)
          : ENERGY_TERM(other, random_number_engine, thread_index, chain),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     ENERGY_TERM *clone(int thread_index=0, typename ENERGY_TERM::ChainType *chain=NULL, 
                Energy<typename ENERGY_TERM::ChainType> *reference_energy_function=NULL) {
          return new Observable<ENERGY_TERM>(*this, thread_index, chain, reference_energy_function);
     }

     //! Make observation.
     //! This default implementation simply calls evaluate_weighted, but this
     //! can be overridden by specializing observable for a specific energy term
     //! \param move_info Object containing information about the last executed move
     //! \param current_iteration Current iteration index in simulation
     //! \param register_only Whether this observation call will not produce any output (merely used to register data)
     std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {
          return boost::lexical_cast<std::string>(ENERGY_TERM::evaluate_weighted(move_info));
     }
};

}

#endif
