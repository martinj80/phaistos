// monte_carlo_muninn.h --- Interface to Muninn monte carlo framework
// Copyright (C) 2006-2012 Wouter Boomsma, Jes Frellsen
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


#ifndef MONTE_CARLO_MUNINN_H
#define MONTE_CARLO_MUNINN_H

#include "monte_carlo/monte_carlo.h"

#include "muninn/Factories/CGEfactory.h"
#include "muninn/UpdateSchemes/IncreaseFactorScheme.h"

#include <algorithm>
#include <iostream>

namespace phaistos {

//! MCMC class - interface to the Muninn generalized ensemble method
template <typename CHAIN_TYPE>
class MonteCarloMuninn: public Simulation<MonteCarloMuninn<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local Simulation
     typedef ::phaistos::Simulation<MonteCarloMuninn<CHAIN_TYPE>, CHAIN_TYPE> Simulation;

private:
     //! Number of updates of the weights
     PHAISTOS_LONG_LONG update_counter;

     //! Iterations of burnin
     PHAISTOS_LONG_LONG burnin_counter;

     //! Energy of next structure
     double energy_next;

     //! Secondary energy current value
     double energy_current_secondary;

     //! Secondary energy next value
     double energy_next_secondary;

     //! Update rule for the burnin counter
     void update_burnin_counter() {
          const Muninn::IncreaseFactorScheme& increase_factor_scheme = dynamic_cast<const Muninn::IncreaseFactorScheme &>(cge->get_ge().get_updatescheme());
          Muninn::Count this_max = increase_factor_scheme.get_this_max();
          burnin_counter = static_cast<PHAISTOS_LONG_LONG>(settings.burnin_fraction * static_cast<double>(this_max));
          std::cout << "Setting burnin to: " << this->burnin_counter << "\n";
     }

public:

     //! Continuous Generalized Ensemble object
     Muninn::CGE *cge;

     //! Whether the CGE is owned by this object
     bool cge_ownership;

     //! Whether weights should be recalculated (pointer shared by threads)
     bool *cge_new_weights;

     //! Secondary energy function
     Energy<CHAIN_TYPE> *energy_function_secondary;


     //! Local settings class.
     const class Settings: public Simulation::Settings, public Muninn::CGEfactory::Settings {
     public:
          //! Lower bound on energy
          double energy_min;

          //! Upper bound on energy
          double energy_max;

          //! Number of histogram updates between every reinitialization (zero or negative means no reinitialization will be done)
          //! If set to =< 0, not reinitialization will be done
          int histograms_per_reinit;

          //! The fraction of the maximum number of iterations in the current simulation round used for burn-in (in each thread)
          double burnin_fraction;
          
          //! Use the secondary energy as an additional energy in muninn - not used to estimate histograms
          bool use_energy_secondary;

          //! Constructor
          Settings(double energy_min=-std::numeric_limits<double>::infinity(),
                   double energy_max=+std::numeric_limits<double>::infinity(),
                   int histograms_per_reinit = 0,
                   double burnin_fraction = 2.0,
                   bool use_energy_secondary=false,
                   int debug=0)
               : Simulation::Settings(debug), Muninn::CGEfactory::Settings(),
                 energy_min(energy_min),
                 energy_max(energy_max),
                 histograms_per_reinit(histograms_per_reinit),
                 burnin_fraction(burnin_fraction),
                 use_energy_secondary(use_energy_secondary) {

               this->initial_beta=UNINITIALIZED;

               if (statistics_log_filename == "") {
                    pid_t pid = getpid();
                    this->statistics_log_filename = "muninn_" + boost::lexical_cast<std::string>(pid) + ".txt";
               }
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "energy-min:" << settings.energy_min << "\n";
               o << "energy-max:" << settings.energy_max << "\n";
               o << "histograms-per-reinit:" << settings.histograms_per_reinit << "\n";
               o << "burnin-fraction:" << settings.burnin_fraction << "\n";
               o << "use-energy-secondary:" << settings.use_energy_secondary << "\n";
               o << dynamic_cast<const typename Muninn::CGEfactory::Settings &>(settings);
               o << dynamic_cast<const typename Simulation::Settings &>(settings);
               return o;
          }
     } settings;    //!< Local settings object


     //! Constructor
     //!
     //! \param chain Molecule chain.
     //! \param energy Energy object
     //! \param move_collection Move set
     //! \param settings Settings object
     //! \param energy_secondary Secondary energy object - energies not used for estimating histograms
     //! \param random_number_engine Object from which random number generators can be created.
     MonteCarloMuninn(CHAIN_TYPE *chain,
                      Energy<CHAIN_TYPE> *energy,
                      MoveCollection<CHAIN_TYPE> *move_collection,
                      const Settings &settings=Settings(),
                      Energy<CHAIN_TYPE> *energy_secondary=NULL,
                      RandomNumberEngine *random_number_engine= &random_global)
          : Simulation("muninn", chain, energy, move_collection, settings, random_number_engine),
            update_counter(0),
            burnin_counter(0),
            energy_next(0.0),
            energy_current_secondary(0.0),
            energy_next_secondary(0.0),
            energy_function_secondary(energy_secondary),
            settings(settings) {

          // Make a local copy of the settings class
          Settings local_settings = settings;

          // Set init_beta to min_beta if not explicitly initialized
          if (!is_initialized(local_settings.initial_beta))
               local_settings.initial_beta = settings.min_beta;
               
          cge = Muninn::CGEfactory::new_CGE(static_cast<Muninn::CGEfactory::Settings>(local_settings));
          cge_ownership = true;
          cge_new_weights = new bool(false);

          // Set the burn-in counter
          update_burnin_counter();
     }

     
     //! Destructor
     ~MonteCarloMuninn() {
          if (cge_ownership) {
               delete cge;
               delete cge_new_weights;
          }

          if (energy_function_secondary && this->pointer_owner)
               delete energy_function_secondary;
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     MonteCarloMuninn(const MonteCarloMuninn &other)
          : Simulation(other),
            update_counter(other.update_counter),
            burnin_counter(other.burnin_counter),
            energy_next(other.energy_next),
            energy_current_secondary(other.energy_current_secondary),
            energy_next_secondary(other.energy_next_secondary),
            cge(other.cge),
            cge_ownership(false),
            cge_new_weights(other.cge_new_weights),
            energy_function_secondary(other.energy_function_secondary),
            settings(other.settings) {}

     
     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     MonteCarloMuninn(const MonteCarloMuninn &other, 
                      RandomNumberEngine *random_number_engine,
                      int thread_index)
          : Simulation(other, random_number_engine, thread_index),
            update_counter(other.update_counter),
            burnin_counter(other.burnin_counter),
            energy_next(other.energy_next),
            energy_current_secondary(other.energy_current_secondary),
            energy_next_secondary(other.energy_next_secondary),
            cge(other.cge),
            cge_ownership(false),
            cge_new_weights(other.cge_new_weights),
            settings(other.settings) {

          if (other.energy_function_secondary)
               energy_function_secondary = new Energy<CHAIN_TYPE>(*other.energy_function_secondary, this->chain, 
                                                                  random_number_engine,
                                                                  thread_index);
          else 
               energy_function_secondary = NULL;
          
     }
     
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

               if (this->burnin_counter==0) {
                    // Update the Muninn count vectors and check whether the weights should be recalculated
                    this->lock();
                    *cge_new_weights = *cge_new_weights || cge->add_observation(this->energy_current);
                    this->unlock();
               } else {
                    this->burnin_counter--;
               }

               // Make move
               bool success = this->move_collection->apply();

               if (!success) {
                    // Reject and revert to last sample
                    this->move_collection->reject();
                    //this->energy_function->reject(); do not reject energy before evaluate
                    this->move_success=false;
                    
               } else {

                    // Calculate energy value
                    energy_next = this->energy_function->evaluate(this->move_collection->get_move_info());

                    // Evaluate energy2
                    if (settings.use_energy_secondary && energy_function_secondary) {
                         energy_next_secondary = energy_function_secondary->evaluate(this->move_collection->get_move_info());
                    }

                    // Reject if energy falls outside legal range
                    if (energy_next >= settings.energy_max || energy_next <= settings.energy_min) {
                         // Reject and revert to last sample		    
                         this->move_collection->reject();
                         this->energy_function->reject();
                         this->move_success = false;

                         if (settings.use_energy_secondary && energy_function_secondary) {
                              energy_function_secondary->reject();
                         }

                         if (settings.debug >= 100)
                              std::cout << "REJECTED " << energy_next << "\n";
                         
                    } else {

                         // Get the energy_bias (this must be within a mutex, since the add_observation call might modify the CGE object)
                         this->lock();
                         double energy_bias = cge->get_lnweights(energy_next) - cge->get_lnweights(this->energy_current);
                         this->unlock();

                         // In some cases, an energy term might add parameters to the state of the simulation, and resample them
                         // internally. This can give rise to a bias
                         energy_bias += this->energy_function->get_log_bias(this->move_collection->get_move_info());

                         // If there is a secondary energyfunction, add its bias
                         if (settings.use_energy_secondary && energy_function_secondary) {
                              energy_bias -= energy_next_secondary - energy_current_secondary;

                              // In some cases, an energy term might add parameters to the state of the simulation, and resample them
                              // internally. This can give rise to a bias
                              energy_bias += energy_function_secondary->get_log_bias(this->move_collection->get_move_info());
                         }

                         // Acceptance criterion
                         if (this->acceptance_criterion(this->move_collection->get_log_bias(), energy_bias)) {

                              // Accept state
                              this->energy_current = energy_next;
                        
                              // Accept structure
                              this->move_collection->accept();
                              this->energy_function->accept();
			      
			      // Accept state and structure for the secondary energy function
                              if (settings.use_energy_secondary && energy_function_secondary) {
                                   this->energy_current_secondary = energy_next_secondary;
                                   this->energy_function_secondary->accept();
                              }

                              if (settings.debug >= 100)
                                   std::cout << "ACCEPTED\n";

                         } else {
                              
                              // Reject and revert to last structure
                              this->move_collection->reject();
                              this->energy_function->reject();
                              this->move_success=false;

			      // Reject and revert to last structure for the secondary energy function
                              if (settings.use_energy_secondary && energy_function_secondary) {
                                   this->energy_function_secondary->reject();
                              }

                              if (settings.debug >= 100)
                                   std::cout << "REJECTED\n";
                         }
                    }
               }

               // This function MUST be called last in every step
               this->finalize_step(step_code);
          }

          // Make sure that all threads are finished with all steps
          this->barrier_wait();

          // Update Muninn weights if necessary
          if (*cge_new_weights) {

               // Make sure that all threads are inside conditional before continuing
               this->barrier_wait();

               // Calculate new estimate of the free energy 
               this->lock();
               // Only one thread is allowed in here
               if (*cge_new_weights) {
                    // Find new weights
                    cge->estimate_new_weights();

                    // Set weights to false to prevent other
                    // threads from updating CGE
                    *cge_new_weights = false;
               }
               this->unlock();

               // Reinitialize chain every 'histograms_per_reinit' 4th time
               if ((settings.histograms_per_reinit > 0) && (this->update_counter % settings.histograms_per_reinit == 0)) {
                    // Note that reinitialize_structure set this->energy_current
                    std::cout<<"# DEVELOP: mcmc: move: reinitialization\n";
                    this->reinitialize_structure(settings.energy_min, settings.energy_max);
                    std::cout<<"# DEVELOP: mcmc: move: reinitialization: done\n";

                    update_burnin_counter();
               }
		    
               this->update_counter++;
          }
     }


     //! Reinitialize chain (overloads base class definition)
     //! \param energy_min Minimum energy
     //! \param energy_max Maximum energy
     void reinitialize_structure(double energy_min=-std::numeric_limits<double>::infinity(), 
                                 double energy_max=std::numeric_limits<double>::infinity()) {
       
          // Limit the arguments this->energy_min and this->energy_max
          energy_min = std::max(energy_min, settings.energy_min);
          energy_max = std::min(energy_max, settings.energy_max);

          // Call base class function
          Simulation::reinitialize_structure(energy_min, energy_max);

          // Set the secondary energy_current
          if (settings.use_energy_secondary && energy_function_secondary) {
               this->energy_current_secondary = energy_function_secondary->evaluate();
               this->energy_function_secondary->accept();
          }
     }


     //! Get energy statistics data (overloads base class definition)
     typename Energy<CHAIN_TYPE>::EnergyDataCollection get_energy_data() {
          typename Energy<CHAIN_TYPE>::EnergyDataCollection data_collection;
          data_collection.add(this->energy_function->get_data("Energy"));
          if (settings.use_energy_secondary && energy_function_secondary)
               data_collection.add(energy_function_secondary->get_data("Energy2"));
          return data_collection;
     }
};

}

#endif
