// monte_carlo_simulated_annealing.h --- Simulated annealing
// Copyright (C) 2006-2010 Mikael Borg, Wouter Boomsma, Jes Frellsen
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

#ifndef MONTE_CARLO_SIMULATED_ANNEALING_H
#define MONTE_CARLO_SIMULATED_ANNEALING_H

#include "moves/move_collection.h"
#include "energy/energy.h"

#include <algorithm>

namespace phaistos {

//! Simulated annealing MC
//! Vanilla Metropolis-Hastings with acceptance criterion 
//! \f[ \min \left[1,\exp{\left(-\Delta E/T\right)} \right] \f]
//! with (optionally) varying temperature.
//! 
//! The temperature at MC step i is
//! \f[ T_i = \mathrm{tstart}\cdot dt^i \f]
//! with
//! \f[ dt = \left( \frac{\mathrm{tend}}{\mathrm{tstart}} \right)^{\frac{1}{n}} \f]
//! where saiter is the number of MC steps for the simulated annealing cycle.
template <typename CHAIN_TYPE>
class MonteCarloSimulatedAnnealing: public Optimization<MonteCarloSimulatedAnnealing<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local Simulation type
     typedef ::phaistos::Optimization<MonteCarloSimulatedAnnealing<CHAIN_TYPE>, CHAIN_TYPE> Optimization;

private:
     //! Energy of next structure
     double energy_next;

     //! Current temperature
     double temperature;

     //! Cooling factor
     double cooling_factor;

public:
     //! Local settings class.
     const class Settings: public Optimization::Settings {
     public:

          //! Lower bound on energy          
          double energy_min;

          //! Upper bound on energy          
          double energy_max;

          //! Initial value of temperature
          double temperature_start;

          //! Final value of temperature
          double temperature_end;

          //! Constructor
          Settings(double energy_min=-std::numeric_limits<double>::infinity(),
                   double energy_max=+std::numeric_limits<double>::infinity(),
                   double temperature_start=std::numeric_limits<double>::infinity(),
                   double temperature_end=1.0,
                   int debug=2)
               : Optimization::Settings(debug),
                 energy_min(energy_min),
                 energy_max(energy_max),
                 temperature_start(temperature_start),
                 temperature_end(temperature_end) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "energy-min:" << settings.energy_min << "\n";
               o << "energy-max:" << settings.energy_max << "\n";
               o << "temperature-start:" << settings.temperature_start << "\n";
               o << "temperature-end:" << settings.temperature_end << "\n";
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
     //! \param n_threads Number of threads
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MonteCarloSimulatedAnnealing(CHAIN_TYPE *chain,
                                  Energy<CHAIN_TYPE> *energy,
                                  MoveCollection<CHAIN_TYPE> *move_collection,
                                  PHAISTOS_LONG_LONG iterations,
                                  int n_threads,
                                  const Settings &settings=Settings(),
                                  RandomNumberEngine *random_number_engine=&random_global)
     : Optimization("simulated-annealing", chain, energy, move_collection, settings, random_number_engine),
       energy_next(0.0),
       temperature(settings.temperature_start),
       settings(settings) {
          cooling_factor = pow(settings.temperature_end/settings.temperature_start, static_cast<double>(n_threads)/static_cast<double>(iterations));
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     MonteCarloSimulatedAnnealing(const MonteCarloSimulatedAnnealing &other)
     : Optimization(other),
       energy_next(other.energy_next),
       temperature(other.temperature),
       cooling_factor(other.cooling_factor),
       settings(other.settings) {}

     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     MonteCarloSimulatedAnnealing(const MonteCarloSimulatedAnnealing &other,
                                  RandomNumberEngine *random_number_engine,
                                  int thread_index)
          : Optimization(other, random_number_engine, thread_index),
            energy_next(other.energy_next),
            temperature(other.temperature),
            cooling_factor(other.cooling_factor),
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
                    energy_next = this->energy_function->evaluate(this->move_collection->get_move_info());

                    // Reject if energy falls outside legal range
                    if (energy_next >= settings.energy_max || energy_next <= settings.energy_min) {
                         // Reject and revert to last sample
                         this->move_collection->reject();
                         this->energy_function->reject();
                         this->move_success = false;
                    } else {
                         double energy_bias = (this->energy_current-energy_next)/temperature;

                         // In some cases, an energy term might add parameters to the state of the simulation, and resample them
                         // internally. This can give rise to a bias
                         energy_bias += this->energy_function->get_log_bias(this->move_collection->get_move_info());

                         // Acceptance criterion
                         if (this->acceptance_criterion(this->move_collection->get_log_bias(), energy_bias)) {
                              // Accept state
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
               }

               // This function MUST be called last in every step
               this->finalize_step(step_code);

               // Change the temperature
               temperature *= cooling_factor;

               if (this->iteration_counter%1000==0)
                    std::cout << "Temperature: " << temperature << std::endl;

          }
     }
};


}

#endif // OPTIMIZATION_SIMULATED_ANNEALING_H

