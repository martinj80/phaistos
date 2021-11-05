// observable_base.h --- Base class for observables.
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

#ifndef OBSERVABLE_BASE_H
#define OBSERVABLE_BASE_H

#include "protein/atom.h"
#include "protein/residue.h"

namespace phaistos {

// Forward declaration
class ObservableCollectionBase;

//! Base class for observables
class ObservableBase {
public: 

     class Settings {
     public:
          //! Interval (number of iterations) with which observable should be registered
          PHAISTOS_LONG_LONG register_interval;

          //! Interval (number of iterations) before observable should start registering
          PHAISTOS_LONG_LONG register_burnin;

          //! How/Where the observable should be reported:
          //!   "logfile": gather information in logfile; 
          //!   "pdb-header":output information to dumped pdb files; 
          //!   "pdb-b-factor":output information to dumped pdb files; 
          //!   "stdout/stderr/cout/cerr": Output to stdout/stderr. 
          //! Any other string is interpreted as a filename for a separate 
          //! logfile - %p is expanded to the process id, %t to the thread index
          //!           and %i is constantly replaced with the current index.
          //! It is possible to override the default style using a \#style suffix.
          //! For intance stdout#compact outputs to stdout in a compact style.
          std::string output_target;

          //! Interval (number of iterations) with which observable should be outputted
          //! This value is ignored if it is lower than register_interval
          PHAISTOS_LONG_LONG output_interval;

          //! Constructor. Defines default values for settings object.
          Settings(PHAISTOS_LONG_LONG register_interval=5000,
                   PHAISTOS_LONG_LONG register_burnin=0,
                   std::string output_target="observables_%p_%t.dat",
                   PHAISTOS_LONG_LONG output_interval=1)
               : register_interval(register_interval),
                 register_burnin(register_burnin),
                 output_target(output_target),
                 output_interval(output_interval) {}

          virtual ~Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "register-interval:" << settings.register_interval << "\n";
               o << "register-burnin:" << settings.register_burnin << "\n";
               o << "output-target:" << settings.output_target << "\n";
               o << "output-interval:" << settings.output_interval << "\n";
               return o;
          }               
     };

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false)=0;

     //! This method is called once all observables have been evaluated, allowing
     //! for communication between them
     virtual void observe_finalize(ObservableCollectionBase *observable_collection) {}

     //! Generate output annotation for residue
     static std::string vector_output_tag(const Residue &residue) {
          return boost::lexical_cast<std::string>(residue.index) + ":";
     }

     //! Generate output annotation for atom
     static std::string vector_output_tag(const Atom &atom) {
          return (boost::lexical_cast<std::string>(atom.residue->index) + "[" + 
                  boost::lexical_cast<std::string>(atom.atom_type)) + "]:";
     }
};

//! Base class for marking an energy terms as begin independent
//! of the normal output system (it controls its own IO)
class IndependentObservable {};

}

#endif
