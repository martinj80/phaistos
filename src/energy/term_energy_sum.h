// energy_sum.h --- Wrap an entire energy into a single term
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


#ifndef ENERGY_SUM_H
#define ENERGY_SUM_H

#include "../energy/energy_term.h"
#include "observable_base.h"

namespace phaistos {

//! No-op. Functionality is in Observable<TermEnergySum>
template <typename CHAIN_TYPE>
class TermEnergySum: public EnergyTermCommon<TermEnergySum<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermEnergySum<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;               

public:
     //! Constructor.
     //! \param chain Molecule chain
     //! \param random_number_engine Object from which random number generators can be created.
     TermEnergySum(CHAIN_TYPE *chain,
                   RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "@energy-sum", typename EnergyTermCommon::Settings(), random_number_engine) {}

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermEnergySum(const TermEnergySum &other, 
                   RandomNumberEngine *random_number_engine,
                   int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {
     }     
};

//!This term wraps an entire energy function into a single term
template <typename CHAIN_TYPE>
class Observable<TermEnergySum<CHAIN_TYPE> >: public TermEnergySum<CHAIN_TYPE>, public ObservableBase {

     //! Reference to main energy function
     Energy<CHAIN_TYPE> *reference_energy_function;

public:
     
     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Whether to evaluate inner energy function. This can be set
          //! to false in cases where the energy has already been evaluated.
          bool evaluate_energy;

          //! Constructor. Defines default values for settings object.
          Settings(bool evaluate_energy=false)
               : evaluate_energy(evaluate_energy) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "evaluate-energy:" << settings.evaluate_energy << "\n";               
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }                    
     } settings;  //!< Local settings object

     //! Constructor.
     //! \param energy_term AngleHistogram energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermEnergySum<CHAIN_TYPE> &energy_term, 
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermEnergySum<CHAIN_TYPE>(energy_term),
            reference_energy_function(reference_energy_function),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const Observable &other, int thread_index, typename TermEnergySum<CHAIN_TYPE>::ChainType *chain, 
                Energy<CHAIN_TYPE> *reference_energy_function)
          : TermEnergySum<CHAIN_TYPE>(other, thread_index, chain),
            reference_energy_function(reference_energy_function),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     TermEnergySum<CHAIN_TYPE> *clone(int thread_index=0, typename TermEnergySum<CHAIN_TYPE>::ChainType *chain=NULL, 
                Energy<CHAIN_TYPE> *reference_energy_function=NULL) {
          return new Observable<TermEnergySum<CHAIN_TYPE> >(*this, thread_index, chain, reference_energy_function);
     }


     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {
          double sum = 0.0;
          if (settings.evaluate_energy)
               this->reference_energy_function->evaluate(move_info);
          sum = settings.weight*this->reference_energy_function->sum;

          return boost::lexical_cast<std::string>(sum);
     }
};

}

#endif
