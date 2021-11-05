// monte_carlo_greedy_optimization.h --- Greedy optimization
// Copyright (C) 2011 Jes Frellsen
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

#ifndef MONTE_CARLO_GREEDY_OPTIMIZATION_H
#define MONTE_CARLO_GREEDY_OPTIMIZATION_H

#include "moves/move_collection.h"
#include "energy/energy.h"
#include "monte_carlo.h"

#include <algorithm>

namespace phaistos {

//! Greedy Optimization method, where only moves that improves the current
//! energy are accepted.
template <typename CHAIN_TYPE>
class MonteCarloGreedyOptimization: public Optimization<MonteCarloGreedyOptimization<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local Optimization type
     typedef ::phaistos::Optimization<MonteCarloGreedyOptimization<CHAIN_TYPE>, CHAIN_TYPE> Optimization;

private:

public:
     //! Local settings class.
     const class Settings: public Optimization::Settings {
     public:
          //! Constructor
          Settings(int debug=0)
               : Optimization::Settings(debug) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << (typename Optimization::Settings &)settings;
               return o;
          }
     } settings;    //!< Local settings object

     //! Constructor
     //!
     //! \param chain Molecule chain.
     //! \param energy Energy object
     //! \param move_collection Move set
     //! \param iterations Number of iterations
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MonteCarloGreedyOptimization(CHAIN_TYPE *chain,
                                  Energy<CHAIN_TYPE> *energy,
                                  MoveCollection<CHAIN_TYPE> *move_collection,
                                  PHAISTOS_LONG_LONG iterations,
                                  const Settings &settings=Settings(),
                                  RandomNumberEngine *random_number_engine=&random_global)
     : Optimization("greedy-optimization", chain, energy, move_collection, settings, random_number_engine),
       settings(settings) {}

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     MonteCarloGreedyOptimization(const MonteCarloGreedyOptimization &other)
     : Optimization(other),
       settings(other.settings) {}

     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     MonteCarloGreedyOptimization(const MonteCarloGreedyOptimization &other,
                                  RandomNumberEngine *random_number_engine,
                                  int thread_index)
          : Optimization(other, random_number_engine, thread_index),
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
                         // Accept structure
                         this->energy_current = energy_next;
                         this->move_collection->accept();
                         this->energy_function->accept();
                    }
                    else {
                         // Reject and revert to last sample
                         this->move_collection->reject();
                         this->energy_function->reject();
                         this->move_success = false;
                    }
               }

               // This function MUST be called last in every step
               this->finalize_step(step_code);

               if (settings.debug > 0 && this->iteration_counter%50==0)
                    std::cout << "Energy: " << this->energy_current << std::endl;

          }
     }
};


}

#endif // MONTE_CARLO_GREEDY_OPTIMIZATION_H

