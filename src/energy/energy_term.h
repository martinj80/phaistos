// energy_term.h --- Energy term base class
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


#ifndef ENERGY_TERM_H
#define ENERGY_TERM_H

#include <cstdlib>
#include "moves/move_info.h"
#include "utils/vector_matrix_3d.h"
#include "utils/settings.h"
#include "utils/random.h"
#include "observable_base.h"

#include <boost/type_traits/is_base_of.hpp>

namespace phaistos {

template <typename CHAIN_TYPE>
class Energy;

template <typename ENERGY_TERM>
class Observable;

//! Energy term. Energies in Phaistos are actually -log probabilities
//! and energy terms can therefore also represent probabilistic models.
template <typename CHAIN_TYPE>
class EnergyTerm {
protected: 

     //! Chain molecule
     CHAIN_TYPE *chain;
public:

     typedef CHAIN_TYPE ChainType;

     //! Name of term
     std::string name;

     //! Additional weight (this non-const weight can be adjusted to vary the weight during simulation)
     double extra_weight;

     //! Thread index
     int thread_index;

     //! Local settings class.
     const class Settings: public ::Settings {
     public:

          //! Weight that will be applied to a term when evaluate_weigted is called.
          //! This is how temperature is implemented.
          double weight;

          //! Level of debug information
          int debug;

          //! Constructor. Defines default values for settings object.
          Settings(double weight = 1.0,
                   int debug=0)
               : weight(weight),
                 debug(debug) {
          }

          //! Destructor
          virtual ~Settings() {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "weight:" << settings.weight << "\n";
               o << "debug:" << settings.debug << "\n";
               return o;
          }          
     } settings; //!< Local settings object 


     //! Settings used by classic energy fields.
     //! Same as normal Settings object, but with
     //! weight set to kT corresponding to 300 Kelvin
     class SettingsClassicEnergy: public Settings {
     public:

          //! 1/kT at room temperature
          static const double one_over_kT_room_temperature;

          //! Constructor
          SettingsClassicEnergy()               
               : Settings(one_over_kT_room_temperature) {}
     };


     //! Constructor.
     //! \param chain Molecule chain
     //! \param name Name of energy term
     //! \param settings Local Settings object
     EnergyTerm(CHAIN_TYPE *chain, std::string name, const Settings &settings=Settings())
          : chain(chain), name(name), extra_weight(1.0), thread_index(0), settings(settings)  {
     }


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     EnergyTerm(const EnergyTerm &other, int thread_index, CHAIN_TYPE *chain)
          : name(other.name),
            extra_weight(other.extra_weight),
            thread_index(thread_index),
            settings(other.settings) {
          if (chain==NULL) {
               this->chain = other.chain;
          } else {
               this->chain = chain;
          }
     }


     //! Destructor
     virtual ~EnergyTerm(){};


     //! Clone: Corresponds to a virtual copy constructor
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     virtual EnergyTerm *clone(RandomNumberEngine *random_number_engine=NULL, int thread_index=0, CHAIN_TYPE *chain=NULL)=0;


     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     virtual double evaluate(MoveInfo *move_info=NULL) {
          return UNINITIALIZED;
     }


     //! Evaluate energy term - weighted. Multiplies the evaluated
     //! energy with the weight associated with the term 
     //! \param move_info Object containing information about the last executed move
     virtual double evaluate_weighted(MoveInfo *move_info=NULL) {
          return get_weight() * this->evaluate(move_info);
     }

     //! Accept last energy evaluation (for caching purposes)
     virtual void accept(){}

     //! Reject last energy evaluation (for caching purposes)
     virtual void reject(){}

     //! Evaluate bias of the energy term.
     //! In cases where an energy term adds extra parameters to the
     //! state of the simulation, the resampling of these values might
     //! be associated with a bias. In such cases this method can be used
     //! to make the bias available to the caller.
     virtual double get_log_bias(MoveInfo *move_info=NULL) {return 0.0;}

     //! Display the settings object of an energy term
     virtual std::string display_settings() {
          return stringify(settings);
     }

     // //! Dump internal data. Mainly used by the iterative ratio method.
     // virtual void dump_internals(const char *filename=NULL) {}

     //! Clone current energy term to the corresponding observable
     virtual EnergyTerm *clone_to_observable(const ObservableBase::Settings &settings, 
                                             Energy<CHAIN_TYPE> *reference_energy_function)=0;

     //! Return weight (sum of constant weight from settings object and dynamic extra weight)
     virtual double get_weight() const {
          return this->settings.weight * this->extra_weight;
     }
};

// 3.2976*10^-27 kcal/K * 6.022*10^23 mol^-1 * temperature
template <typename CHAIN_TYPE>
const double EnergyTerm<CHAIN_TYPE>::SettingsClassicEnergy::one_over_kT_room_temperature = 1.0/(3.2976E-27*6.022E23 * 300);




//! Layer between user-defined energy terms and the EnergyTerm base class.
//! This class takes the derived class as a template argument, and can 
//! therefore define common methods for all energy terms that would normally
//! have to be defined explicitly in each term.
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class EnergyTermCommon: public EnergyTerm<CHAIN_TYPE> {
public:
     //! Local settings definition.
     typedef typename EnergyTerm<CHAIN_TYPE>::Settings Settings;

protected: 
     
     //! Constructor.
     //! \param chain Molecule chain
     //! \param name Name of energy term
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     EnergyTermCommon(CHAIN_TYPE *chain, std::string name, const Settings &settings,
                      RandomNumberEngine *random_number_engine)
          : EnergyTerm<CHAIN_TYPE>(chain, name, settings) {}

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     EnergyTermCommon(const EnergyTermCommon &other, 
                      RandomNumberEngine *random_number_engine,
                      int thread_index, 
                      CHAIN_TYPE *chain)
          : EnergyTerm<CHAIN_TYPE>(other, thread_index, chain){}
     

     //! Clone: Corresponds to a virtual copy constructor
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     EnergyTerm<CHAIN_TYPE> *clone(RandomNumberEngine *random_number_engine=NULL, int thread_index=0, CHAIN_TYPE *chain=NULL) {
          return new DERIVED_CLASS(dynamic_cast<DERIVED_CLASS const&>(*this), random_number_engine, thread_index, chain);
     }

     //! Display the settings object of an energy term
     virtual std::string display_settings() {
          return stringify(((DERIVED_CLASS*)(this))->settings);
     }

     //! Transform energy term into observable
     //! \param settings ObservableBase Settings object
     //! \param reference_energy_function Energy function to clone.
     virtual EnergyTerm<CHAIN_TYPE> *clone_to_observable(const ObservableBase::Settings &settings,
                                                         Energy<CHAIN_TYPE> *reference_energy_function) {
          // return new Observable<DERIVED_CLASS>(dynamic_cast<DERIVED_CLASS const&>(*this), settings, reference_energy_function);
          return clone_to_observable_inner(boost::is_base_of<ObservableBase, DERIVED_CLASS>(), settings, reference_energy_function);
     }

private:

     //! clone_to_observable functionality in case term is already an observable
     //! \param settings ObservableBase Settings object
     //! \param reference_energy_function Energy function to clone.
     EnergyTerm<CHAIN_TYPE> *clone_to_observable_inner(const boost::true_type &t, 
                                                       const ObservableBase::Settings &settings,
                                                       Energy<CHAIN_TYPE> *reference_energy_function) {
          return new DERIVED_CLASS(dynamic_cast<DERIVED_CLASS const&>(*this), settings, reference_energy_function);
     }

     //! clone_to_observable functionality in case term is not an observable
     //! \param settings ObservableBase Settings object
     //! \param reference_energy_function Energy function to clone.
     EnergyTerm<CHAIN_TYPE> *clone_to_observable_inner(const boost::false_type &t, 
                                                       const ObservableBase::Settings &settings,
                                                       Energy<CHAIN_TYPE> *reference_energy_function) {
          return new Observable<DERIVED_CLASS>(dynamic_cast<DERIVED_CLASS const&>(*this), settings, reference_energy_function);
     }
};

}
#endif
