// monte_carlo_metropolis_hastings.h --- Metropolis Hastings monte carlo simulation
// Copyright (C) 2006-2008 Wouter Boomsma
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


#ifndef MONTE_CARLO_METROPOLIS_HASTINGS_H
#define MONTE_CARLO_METROPOLIS_HASTINGS_H

#include "monte_carlo.h"

namespace phaistos {

//! Metropolis-Hastings MCMC class
template <typename CHAIN_TYPE>
class MonteCarloMetropolisHastings: public Simulation<MonteCarloMetropolisHastings<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local Simulation type
     typedef ::phaistos::Simulation<MonteCarloMetropolisHastings<CHAIN_TYPE>, CHAIN_TYPE> Simulation;     

public:

     //! Use Settings object from base class
     typedef typename Simulation::Settings Settings;

     //! Settings object
     Settings settings;

     //! Constructor
     //!
     //! \param chain Molecule chain.
     //! \param energy Energy object
     //! \param move_collection Move set
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MonteCarloMetropolisHastings(CHAIN_TYPE *chain,
                                  Energy<CHAIN_TYPE> *energy, 
                                  MoveCollection<CHAIN_TYPE> *move_collection,
                                  const Settings &settings=Settings(),
                                  RandomNumberEngine *random_number_engine= &random_global)
          : Simulation("metropolis-hastings",
                       chain, energy, move_collection, 
                       settings,
                       random_number_engine),
            settings(settings) {}

     
     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     MonteCarloMetropolisHastings(const MonteCarloMetropolisHastings &other)
          : Simulation(other),
            settings(other.settings) {}
     
     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     MonteCarloMetropolisHastings(const MonteCarloMetropolisHastings &other, 
                                  RandomNumberEngine *random_number_engine,
                                  int thread_index)
          : Simulation(other, random_number_engine, thread_index),
            settings(other.settings) {}
     

     //! Make a move using the move collection, and evaluate energy
     //! using the energy object
     //!
     //! \param steps Number of steps to make
     //! \param step_code Optional functor to evaluate in each iteration
     void move(int steps=1, 
               MonteCarloStepCode<CHAIN_TYPE> *step_code=NULL) {

          for (int i=0; i<steps; i++) {

               // This function MUST be called first in every step
               this->initialize_step();

               // Make move
               bool success = this->move_collection->apply();

               if (!success) {
                    // Reject and revert to last sample
                    this->move_collection->reject();
                    //this->energy_function->reject(); do not reject energy before evaluate
                    this->move_success=false;

               } else {

                    // Calculate energy value
                    double energy_next = this->energy_function->evaluate(this->move_collection->get_move_info());

                    double move_bias = this->move_collection->get_log_bias();
                    double energy_bias = -(energy_next - this->energy_current);

                    // In some cases, an energy term might add parameters to the state of the simulation, and resample them
                    // internally. This can give rise to a bias
                    energy_bias += this->energy_function->get_log_bias(this->move_collection->get_move_info());

                    // Acceptance criterion
                    if (this->acceptance_criterion(move_bias, 
                                                   energy_bias)) {
                    
                         this->energy_current = energy_next;

                         // Accept structure
                         this->move_collection->accept();
                         this->energy_function->accept();

                    } else {

                         // Reject and revert to last structure
                         this->move_collection->reject();
                         this->energy_function->reject();
                         this->move_success=false;
                    }

               }

               // This function MUST be called last in every step
               this->finalize_step(step_code);
          }               
     }
};

}

#endif
