// q_factor.h --- Measure of nativeness.
//                As described in Lindorff-Larsen, Piana, Dror, Shaw, Science, 2011
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

#ifndef Q_FACTOR_H
#define Q_FACTOR_H

#include "../energy/energy_term.h"
#include "observable.h"
#include "term_contact_map.h"


namespace phaistos {

//! Measure the nativeness of a chain
template <typename CHAIN_TYPE>
class TermQFactor: public TermContactMapBase<TermQFactor<CHAIN_TYPE>, CHAIN_TYPE> {
protected:

     //! For convenience, define local TermContactMap
     typedef phaistos::TermContactMapBase<TermQFactor<CHAIN_TYPE>,CHAIN_TYPE> TermContactMapBase;

     //! For convenience, define local Contact type
     typedef typename TermContactMapBase::Contact Contact;

public:

     //! Local settings class.
     const class Settings: public TermContactMapBase::Settings {
     public:
          
          // Whether to use contact weights to calculate a weighted version of
          // the q factor
          bool use_weights;

          //! Constructor
          Settings(bool use_weights = false)
               : use_weights(use_weights) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "use-weights: " << settings.use_weights << "\n";               
               o << (typename TermContactMapBase::Settings &)settings;
               return o;
          }

     } settings;  //!< Local settings object 

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermQFactor(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : TermContactMapBase(chain, "q-factor", settings, random_number_engine) {
     }


     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermQFactor(const TermQFactor &other, 
                 RandomNumberEngine *random_number_engine,
                 int thread_index, CHAIN_TYPE *chain)
          : TermContactMapBase(other, random_number_engine, thread_index, chain) {}


     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {

          double q = 0;

          double weight_sum = 0.0;

          for (unsigned int i=0; i<this->contact_map.size(); ++i) {
               
               Contact &contact = this->contact_map[i];

               double &reference_distance = contact.distance;

               Atom *atom1 = (*this->chain)(contact.residue_index1, contact.atom_type1);
               Atom *atom2 = (*this->chain)(contact.residue_index2, contact.atom_type2);
               
               double distance = (atom2->position - atom1->position).norm();

               double weight = 1.0;
               if (settings.use_weights)
                    weight = contact.weight;

               weight_sum += weight;

               q += weight/(1+std::exp(10*(distance - reference_distance - 1)));
          }

          q /= weight_sum;
          
          return q;
     }

};


}

#endif


