// montecarlo.h --- Base class for optimization and simulation code
// Copyright (C) 2006-2009 Wouter Boomsma
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


#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <boost/thread.hpp>

#include "utils/random.h"

#include "moves/move_collection.h"
#include "energy/energy.h"
#include "energy/term_clash_fast.h"

namespace phaistos {

// Forward declaration
template <typename CHAIN_TYPE> class MonteCarlo;

//! Functor definition for per-step functionality. 
//! Inherit from this class to define specific code
//! executed in each step of a monte carlo run
template <typename CHAIN_TYPE>
class MonteCarloStepCode {
public:
     //! Overload () operator. Execute code
     virtual void operator()(MonteCarlo<CHAIN_TYPE> *monte_carlo) {}

     //! Destructor
     virtual ~MonteCarloStepCode(){};
};

//! Statistics of Monte Carlo object
class MonteCarloStatistics {
public:
     //! Local class with summary information
     class Summary {
     private:
          //! Constant with the number of seconds/day
          const static int seconds_per_day = 60*60*24;

     public:

          //! Average time per step
          std::vector<double> time_avg;

          //! Steps done per day
          std::vector<double> steps_per_day;

          //! Constructor
          //! \param mcs MonteCarloStatistics object
          Summary(const MonteCarloStatistics &mcs) {
               time_avg = std::vector<double>(1, mcs.time_spent/((double)mcs.observations));
               steps_per_day = std::vector<double>(1, seconds_per_day/time_avg[0]);
          }

          //! Constructor - create a single summary from a vector of summaries
          //! \param summaries vector of MonteCarloStatistics objects
          Summary(const std::vector<Summary> &summaries) {
               time_avg.resize(summaries.size());
               steps_per_day.resize(summaries.size());
               for (unsigned int i=0; i<summaries.size(); i++) {
                    time_avg[i] = summaries[i].time_avg[0];
                    steps_per_day[i] = summaries[i].steps_per_day[0];
               }
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &o, const Summary &mcs) {
               if (mcs.time_avg.size()==1 )
                    o << "Average step time: " << mcs.time_avg[0] << " sec.    \tSteps/day: " << mcs.steps_per_day[0] << ".";
               else {
                    double steps_per_day_total = 0;
                    for (unsigned int i=0; i<mcs.steps_per_day.size(); i++) {
                         steps_per_day_total += mcs.steps_per_day[i];
                    }
                    o << "Average Step time: " << mcs.time_avg << " sec.    \tSteps/day: " << mcs.steps_per_day << ".    \tTotal: " << steps_per_day_total << " steps/day.";
               }
               return o;
          }          
          
     };

     //! Number of observations
     int observations;

     //! Time start mark
     double time_start;

     //! Measure of time spent
     double time_spent;

     //! Constructor
     MonteCarloStatistics()
          : observations(0),
            time_spent(0) {}

     //! Register beginning of step
     void register_step_start() {
          time_start = get_time();
     }

     //! Register end of step
     void register_step_end() {
          time_spent += get_time() - time_start;
          observations++;
     }

     //! Retrieve summary
     Summary get_summary() {
          return Summary(*this);
     }
};


//! Base class for all optimization and simulation algorithm classes
template <typename CHAIN_TYPE>
class MonteCarloBase {

public:

     //! Number of iterations done
     PHAISTOS_LONG_LONG iteration_counter;

     //! Status of last move
     bool move_success;

     //! Constructor
     MonteCarloBase()
          : iteration_counter(0) {}

     //! Destructor
     virtual ~MonteCarloBase(){};

     //! Execute move
     //! \param steps Number of steps to make
     //! \param step_code Optional functor to evaluate in each iteration
     virtual void move(int steps=1, 
                       MonteCarloStepCode<CHAIN_TYPE> *step_code=NULL)=0;

     //! Reinitialize chain
     //! \param energy_min Minimum energy
     //! \param energy_max Maximum energy     
     virtual void reinitialize_structure(double energy_min=-std::numeric_limits<double>::infinity(), 
                                         double energy_max=std::numeric_limits<double>::infinity())=0;

     //! Get energy (or energies in the case of a MonteCarloMultiThread)
     virtual std::vector<double> get_energy()=0;

     //! Get energy function (or energy functions in the case of a MonteCarloMultiThread)
     virtual std::vector<Energy<CHAIN_TYPE>*> get_energy_function()=0;

     //! Get chain (or chains in the case of a MonteCarloMultiThread)
     virtual std::vector<CHAIN_TYPE*> get_chain()=0;

     //! Get monte carlo statistics
     virtual MonteCarloStatistics::Summary get_statistics()=0;

     //! Get move statistics
     virtual typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary get_move_statistics()=0;

     //! Get energy data
     virtual typename Energy<CHAIN_TYPE>::EnergyDataCollection get_energy_data()=0;

     //! Display settings
     virtual void display_settings()=0;
};


//! Base class for single threaded monte_carlo simulations/optimizations
template <typename CHAIN_TYPE>
class MonteCarlo: public MonteCarloBase<CHAIN_TYPE> {
public:

     //! Save template argument into local type
     typedef CHAIN_TYPE Chain;

     //! Identification string
     std::string id;

     //! Molecule chain
     CHAIN_TYPE *chain;

     //! Energy function
     Energy<CHAIN_TYPE> *energy_function;

     //! Move collection
     MoveCollection<CHAIN_TYPE> *move_collection;

     //! Whether object owns chain, energy_function and move_collection pointers
     bool pointer_owner;

     //! Current energy value
     double energy_current;

     //! Random_number_engine Object from which random number generators can be created.
     RandomNumberEngine *random_number_engine;     

     //! Index of thread in which object is running
     int thread_index;

     //! Monte Carlo statistics object
     MonteCarloStatistics statistics;

     //! Local settings class.
     const class Settings: public ::Settings {
     public:

          //! Debug level
          int debug;

          //! Whether do declash chain when reinitializing
          bool declash_on_reinitialize;

          //! The number of times declashing is attempted before a complete reinitialization is done
          int maximum_declash_attempts;

          //! How often the chain is reinitialized
          int reinitialization_interval;

          //! How often is the consistency of the chain tested
          int consistency_check_interval;

          //! Constructor
          Settings(int debug=0,
                   bool declash_on_reinitialize=true, 
                   int reinitialization_interval = 0,
                   int consistency_check_interval = 10000)
               : debug(debug), 
                 declash_on_reinitialize(declash_on_reinitialize),
                 maximum_declash_attempts(20000),
                 reinitialization_interval(reinitialization_interval),
                 consistency_check_interval(consistency_check_interval) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "debug:" << settings.debug << "\n";
               o << "declash-on-reinitialize:" << settings.declash_on_reinitialize << "\n";
               o << "maximum-declash-attempts:" << settings.maximum_declash_attempts << "\n";
               o << "reinitialization-interval:" << settings.reinitialization_interval << "\n";
               o << "consistency-check-interval:" << settings.consistency_check_interval << "\n";
               return o;
          }                              
     } settings;    //!< Local settings object


     //! Constructor
     //!
     //! \param id Identification string
     //! \param chain Molecule chain.
     //! \param energy_function Energy object
     //! \param move_collection Move set
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MonteCarlo(std::string id,
                CHAIN_TYPE *chain,
                Energy<CHAIN_TYPE> *energy_function, 
                MoveCollection<CHAIN_TYPE> *move_collection,
                const Settings &settings=Settings(),
                RandomNumberEngine *random_number_engine=&random_global)
          : id(id),
            chain(chain),
            energy_function(energy_function),
            move_collection(move_collection),
            pointer_owner(false),
            energy_current(std::numeric_limits<double>::infinity()),
            random_number_engine(random_number_engine),
            thread_index(0),
            settings(settings) {}


     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     MonteCarlo(const MonteCarlo &other) {
          assert(false);
     }


     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     MonteCarlo(const MonteCarlo &other, RandomNumberEngine *random_number_engine, int thread_index)
          : MonteCarloBase<CHAIN_TYPE>(other),
            id (other.id),
            energy_current(other.energy_current), 
            random_number_engine(random_number_engine),
            thread_index(thread_index),
            settings(other.settings) {

          chain = new CHAIN_TYPE(*other.chain);                    
          energy_function = new Energy<CHAIN_TYPE>(*other.energy_function, chain, 
                                                   random_number_engine,
                                                   thread_index);
          move_collection = new MoveCollection<CHAIN_TYPE>(*other.move_collection, chain, 
                                                           random_number_engine, 
                                                           thread_index);
          pointer_owner = true;
     }

     //! Destructor
     virtual ~MonteCarlo(){
          if (pointer_owner) {
               delete chain;
               delete energy_function;
               delete move_collection;
          }
     }


     //! Initializer. Should be called as the first function in every step
     void initialize_step() {

          // Register move with statistic object
          this->statistics.register_step_start();

          //! Assume move is succesful - the derived classes can overwrite this value
          this->move_success = true;

          //! Reinitialize structure at certain intervals
          if (this->iteration_counter==0 || 
              (settings.reinitialization_interval > 0 && 
               this->iteration_counter % settings.reinitialization_interval == 0)) {
               reinitialize_structure();
          }

          //! Check consistency at certain intervals
          if (this->iteration_counter % settings.consistency_check_interval == 0) {
               chain->check_consistency();
          }
     }

     //! Should be called as the final function in every step
     virtual void finalize_step(MonteCarloStepCode<CHAIN_TYPE> *step_code=NULL) {

          // Increment iteration counters
	  this->iteration_counter++;

          // Executing user-specified code
          if (step_code)
               (*step_code)(this);

          // Register time spent
          this->statistics.register_step_end();

          if (settings.debug > 9) {
               std::cout << "Move: " << this->move_collection->get_last_move()->id << " - ";
               std::cout << (this->move_success?"accepted":"rejected") << "\n";
               std::cout << *this->move_collection->get_move_info() << "\n";
          }
     }

     //@{
     //! These functions must be present but are non-ops in the non-threaded base class
     virtual void lock() {} 
     virtual void unlock() {}
     virtual void barrier_wait() {}
     //@}


     //! Local version of TermClashFast - keeps track of a vector of interacting pairs
     //! so they can be accessed from the outside
     class TermClashFastLocal: public TermClashFast<CHAIN_TYPE> {
     public:
          
          //! Define NodeType locally for ease of reference
          typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

          //! Vector of interacting atoms
          std::vector<chaintree::PairWithDistance<Atom*,Atom*> > &interacting_atoms;

          //! Constructor
          //! \param chain Molecule chain
          //! \param interacting_atoms destination vector of interacting pairs
          //! \param only_modified_pairs Specifies whether energy should consider only pairs modified by the last
          //! move. If a clashfree structure is maintained at all times, this is sufficient
          //! to detect all clashes
          //! \param boolean_mode Specifies whether energy should work in clash/non-clash mode (true)
          //! our count all the clashes (false)
          TermClashFastLocal(CHAIN_TYPE *chain,
                             std::vector<chaintree::PairWithDistance<Atom*,Atom*> >  &interacting_atoms,
                             bool only_modified_pairs=false,
                             bool boolean_mode=false)
               : TermClashFast<CHAIN_TYPE>(chain, 
                                           typename TermClashFast<CHAIN_TYPE>::Settings(only_modified_pairs,
                                                                                        boolean_mode)),
                 interacting_atoms(interacting_atoms){}
          

          //! Evaluate energy term.
          //! \param move_info Object containing information about the last executed move
          double evaluate(MoveInfo *move_info=NULL) {
               int clash_counter = 0;
               interacting_atoms.clear();

               // Iterator over node pairs
               for (chaintree::PairIterator<CHAIN_TYPE, NodeType, NodeType> it(*this->chain, this->it_settings); !it.end(); ++it) {

                    // Check if nodes are identical
                    if (it->first == it->second) {

                         // Iterate over all atoms in first node
                         for (unsigned int i=0; i<it->first->size(); ++i) {

                              // Iterate over all atoms in second node
                              for (unsigned int j=i+1; j<it->second->size(); ++j) {
                                   Atom *atom1 = (*it->first)[i];
                                   Atom *atom2 = (*it->second)[j];
                                   double distance_squared;
                                   if (this->detect_clash(atom1, atom2, &distance_squared)) {
                                        interacting_atoms.push_back(chaintree::PairWithDistance<Atom*,Atom*>(atom1, atom2, distance_squared));
                                        ++clash_counter;
                                   }
                              }
                         }
                    } else {
                         // Iterate over all atoms in first node
                         for (unsigned int i=0; i<it->first->size(); ++i) {

                              // Iterate over all atoms in second node
                              for (unsigned int j=0; j<it->second->size(); ++j) {
                                   Atom *atom1 = (*it->first)[i];
                                   Atom *atom2 = (*it->second)[j];
                                   double distance_squared;
                                   if (this->detect_clash(atom1, atom2, &distance_squared)) {
                                        interacting_atoms.push_back(chaintree::PairWithDistance<Atom*,Atom*>(atom1, atom2, distance_squared));
                                        ++clash_counter;
                                   }
                              }
                         }
                    }
               }               

               return (double)clash_counter;
          }
     };


     //! \param energy_min Minimum energy
     //! \param energy_max Maximum energy          
     virtual void reinitialize_structure(double energy_min=-std::numeric_limits<double>::infinity(), 
                                         double energy_max=std::numeric_limits<double>::infinity()) {

          std::cout << "Monte Carlo: Reinitializing structure\n";

          // Reinitialize and intialize currentBin
          move_collection->reinitialize();

          // Evaluate normal energy
          energy_current = energy_function->evaluate();

	  
          //! vector of interacting atoms
          std::vector<chaintree::PairWithDistance<Atom*,Atom*> > clash_interacting_atoms;

          //! new proposed vector of interacting atoms
          std::vector<chaintree::PairWithDistance<Atom*,Atom*> > clash_interacting_atoms_new;

          //! Energy function
          Energy<CHAIN_TYPE> energy_function_simple(move_collection->chain);
          if (settings.declash_on_reinitialize) {
               energy_function_simple.add_term(new TermClashFastLocal(move_collection->chain, 
                                                                      clash_interacting_atoms_new));
          }
	  
          // Evaluate simple energy 
          // (move info is not sent along, since result is uncachable - the energy has just been created)
          double energy_simple = energy_function_simple.evaluate();

          // Accept energy evaluations (a reinitialization cannot be rejected)
          energy_function->accept();
          energy_function_simple.accept();
          clash_interacting_atoms = clash_interacting_atoms_new;

          // Check if energy is within legal range and is clash free according to
          // energy_simple
          if (!std::isfinite(energy_current) || 
	      ((energy_current > energy_max || energy_current < energy_min) || energy_simple>0.5)) {
	       
	       // Make sure that whenever declashing is done, we accept the first new energy
	       energy_current = std::numeric_limits<double>::infinity();

               std::cout << "# MonteCarlo: Resampling to ensure that energies are clashfree and within range\n";
	       PHAISTOS_LONG_LONG counter = 0;
               while(1) {

                    while (!move_collection->apply()) {move_collection->reject();}

                    double energy_simple_new = energy_function_simple.evaluate(move_collection->get_move_info());

                    // Reject either if the new structure has more clashes
                    // or if they both have 1 clash, and the distance of the new one
                    // is smaller than the old one (more serious clash)
                    if (energy_simple_new > energy_simple || (energy_simple > 0.5 &&
                        (energy_simple_new < 1.5 && 
                         clash_interacting_atoms_new[0].distance < clash_interacting_atoms[0].distance))) {
                         move_collection->reject();
                         energy_function_simple.reject();
                    } else {

                         // Set new energy as current energy
                         energy_simple = energy_simple_new;

                         // Set new interacting atoms as current interacting atoms
                         clash_interacting_atoms = clash_interacting_atoms_new;

                         // If we find no more clashes, evaluate normal energy
                         if (energy_simple_new < 0.5) {

			      // Evaluate full energy
                              double energy_new = energy_function->evaluate();

			      if (energy_new <= energy_current) {
				   energy_current = energy_new;

				   // This call to move_collection->accept must be placed after
				   // the call to energy_function->evaluate(), for the syncing
				   // between time stamps in chain and energy caches to work.
				   move_collection->accept();
				   energy_function_simple.accept();
				   energy_function->accept();

                                   // Break if full energy is within range
				   if ((energy_current < energy_max) && (energy_current > energy_min)) {
					break;
				   }

			      } else {
				   move_collection->reject();
				   energy_function_simple.reject();
				   energy_function->reject();
			      }

			      std::cout << "# MonteCarlo reinitialize proposal " << counter 
                                        << " (0 clashes): Evaluating full energy. Current: " 
                                        << energy_current << ". Proposed: " << energy_new 
                                        << ". Range=(" << energy_min << "," << energy_max << ")\n";

                         } else {
                              move_collection->accept();
			      energy_function_simple.accept();
                         }
                    }

		    if (counter%100==0 && energy_simple>0.5) {
		         std::cout << "# MonteCarlo reinitialize proposal "  << counter <<  ". Clashes: " <<  energy_simple << "\t    ";
                         for (unsigned int i=0; i<clash_interacting_atoms.size(); ++i) {
                              chaintree::PairWithDistance<Atom*,Atom*> &p = clash_interacting_atoms[i];
                              if (i!=0)
                                   std::cout << ", ";
                              std::cout << p.first->atom_type << "(" << p.first->residue->index << ")-" <<
                                   p.second->atom_type << "(" << p.second->residue->index << ")=" << 
                                   (p.first->position-p.second->position).norm();
                         }
                         std::cout << "\n";
		    }

		    counter++;

                    // After a certain amount of attempt, completely reinitialize to try again
                    if (counter%settings.maximum_declash_attempts==0) {
                         
                         std::cout << "# MonteCarlo: Too many attempts. Doing complete reinitialization\n";

                         // Reinitialize and intialize currentBin
                         move_collection->reinitialize();
                    
                         // Evaluate normal energy
                         energy_current = energy_function->evaluate();

                         // Evaluate simple energy
                         energy_simple = energy_function_simple.evaluate(move_collection->get_move_info());

                         // Set new interacting atoms as current interacting atoms
                         clash_interacting_atoms = clash_interacting_atoms_new;

                         // Accept energy evaluations (a reinitialization cannot be rejected)
                         energy_function->accept();
                         energy_function_simple.accept();
                    }
               }
          }
          std::cout << "# MonteCarlo: ok - energy: " << energy_current << "\n";
     }


     //! Return vector of energy values.
     //! \return The returned  vector just contains a single value - but the interface is
     //! shared with the multithreaded case, in which case a vector is returned
     std::vector<double> get_energy() {
          return std::vector<double>(1, energy_current);
     };

     //! Return vector of energy function object pointers.
     //! \return The returned  vector just contains a single value - but the interface is
     //! shared with the multithreaded case, in which case a vector is returned
     std::vector<Energy<CHAIN_TYPE>*> get_energy_function() {
          return std::vector<Energy<CHAIN_TYPE>*>(1, energy_function);
     }

     //! Return vector of chain object pointers.
     //! \return The returned  vector just contains a single value - but the interface is
     //! shared with the multithreaded case, in which case a vector is returned
     std::vector<CHAIN_TYPE*> get_chain() {
          return std::vector<CHAIN_TYPE*>(1, chain);
     }

     //! Get monte carlo statistics
     MonteCarloStatistics::Summary get_statistics() {
          return statistics.get_summary();
     }

     //! Get move statistics
     typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary get_move_statistics() {
          return move_collection->get_statistics();
     }

     //! Get energy data
     typename Energy<CHAIN_TYPE>::EnergyDataCollection get_energy_data() {
          return typename Energy<CHAIN_TYPE>::EnergyDataCollection(energy_function->get_data());
     }
};



//! Multithreaded version of MonteCarlo. Constructed using a MonteCarlo object, it
//! creates the necessary copies for all threads.
template <typename MONTE_CARLO_TYPE>
class MonteCarloMultiThread: public MonteCarloBase<typename MONTE_CARLO_TYPE::Chain> {

protected:

     //! Monte Carlo object extended with threading information
     class MonteCarloThreaded: public MONTE_CARLO_TYPE {
     public:
          //! Pointer to thread
          boost::thread *thread;

          //! Pointer to barrier          
          boost::barrier *barrier;

          //! Pointer to mutex
          boost::mutex *mutex;

          //! Constructor - constructed based on non-threaded monte carlo object
          //! \param other Source object
          //! \param random_number_engine Object from which random number generators can be created.
          //! \param thread_index Index indicating in which thread|rank the copy exists
          MonteCarloThreaded(const MONTE_CARLO_TYPE &other, 
                             RandomNumberEngine *random_number_engine,
                             int thread_index)
               : MONTE_CARLO_TYPE(other, random_number_engine, thread_index) {}


          //! Lock mutex
          void lock() { 
               mutex->lock();
          }

          //! Unlock mutex
          void unlock() {
               mutex->unlock();
          }

          //! Synchronize all threads
          void barrier_wait() {
               barrier->wait();
          }
     };

     //! Template from which the threads are created
     MONTE_CARLO_TYPE *monte_carlo_template;

     //! The n_threads MonteCarlo copies
     std::vector<MonteCarloThreaded *> threads;

     //! Number of steps executed per each call of MonteCarloMultiThread::move
     int steps_per_move;

     //! Vector of random number engines - one for each copy 
     std::vector<RandomNumberEngine *> random_number_engines;

     //! Whether object owns the random number engine
     bool random_number_engine_ownership;

public:

     //! Constructor
     //!
     //! \param monte_carlo non-threaded MonteCarlo object from which threaded version is created
     //! \param n_threads Number of threads to create
     //! \param steps_per_move Number of steps executed per each call of MonteCarloMultiThread::move
     //! \param random_number_engines vector of objects from which random number generators can be created.
     MonteCarloMultiThread(MONTE_CARLO_TYPE *monte_carlo, unsigned int n_threads, int steps_per_move,
                           const std::vector<RandomNumberEngine *> &random_number_engines = 
                           std::vector<RandomNumberEngine*>())
          : monte_carlo_template(monte_carlo), 
            steps_per_move(steps_per_move) {

          this->threads.resize(n_threads);
          this->random_number_engines.resize(n_threads);

          // If a full-sized vector was passed along, use these engines
          if (random_number_engines.size() == n_threads) {
               for (unsigned int i=0; i<n_threads; i++) {
                    this->random_number_engines[i] = random_number_engines[i];
               }
               random_number_engine_ownership = false;

          // If no random number engines were passed, generate them 
          } else {

               // Generate remaining generators from first one
               for (unsigned int i=1; i<n_threads; i++) {
                    this->random_number_engines[i] = new RandomNumberEngine();
                    this->random_number_engines[i]->seed(*random_number_engines[0]);
               }
               random_number_engine_ownership = true;
          }

          // Make n_threads number of copies of monteCarlo object
          for (unsigned int i=0; i<n_threads; i++) {
               threads[i] = new MonteCarloThreaded(*monte_carlo, random_number_engines[i], i);
          }
     }

     //! Destructor
     ~MonteCarloMultiThread() {

          // Delete template
          delete monte_carlo_template;

          // Delete all copies generated during construction
          for (unsigned int i=0; i<threads.size(); i++) {
               delete threads[i];
          }

          // Delete all random number generators generated during construction
          if (random_number_engine_ownership) {
               for (unsigned int i=1; i<threads.size(); i++) {
                    delete random_number_engines[i];
               }
          }
     }

     //! Execute move
     //! \param steps Number of steps to make
     //! \param step_code Optional functor to evaluate in each iteration
     void move(int steps=1, 
               MonteCarloStepCode<typename MONTE_CARLO_TYPE::Chain> *step_code=NULL) {

          boost::barrier barrier(threads.size());
          boost::mutex mutex;
          for (unsigned int i=0; i<threads.size(); i++) {

               // threads[i]->m_thread = new boost::thread(boost::bind(&MonteCarloThreaded::move, threads[i], _1), 
               //                                          steps*steps_per_move);
               threads[i]->thread = new boost::thread(boost::bind(&MonteCarloThreaded::move, threads[i], _1, _2), 
                                                      steps*steps_per_move, step_code);

               // The barrier and mutex objects are shared between threads
               threads[i]->barrier = &barrier;
               threads[i]->mutex = &mutex;
          }

          // Join threads
          for (unsigned int i=0; i<threads.size(); i++) {
               threads[i]->thread->join();
               delete threads[i]->thread;
          }

          // Update iteration counter
          // this->iteration_counter += steps*steps_per_move*threads.size();
          this->iteration_counter += steps*steps_per_move;
     }


     //! Reinitialize chain - for all threads
     //! \param energy_min Minimum energy
     //! \param energy_max Maximum energy     
     virtual void reinitialize_structure(double energy_min=-std::numeric_limits<double>::infinity(), 
                                         double energy_max=std::numeric_limits<double>::infinity()) {
          
          for (unsigned int i=0; i<threads.size(); i++) {
               threads[i]->reinitialize_structure();
          }          
     }

     //! Return vector of energy values.
     std::vector<double> get_energy() {
          std::vector<double> energies;
          for (unsigned int i=0; i<threads.size(); i++) {
               energies.push_back(threads[i]->get_energy()[0]);
          }
          return energies;
     };     

     //! Return vector of energy function object pointers.
     std::vector<Energy<typename MONTE_CARLO_TYPE::Chain>*> get_energy_function() {
          std::vector<Energy<typename MONTE_CARLO_TYPE::Chain>*> energy_functions;
          for (unsigned int i=0; i<threads.size(); i++) {
               energy_functions.push_back(threads[i]->get_energy_function()[0]);
          }
          return energy_functions;
     }

     //! Return vector of chain object pointers.
     std::vector<typename MONTE_CARLO_TYPE::Chain*> get_chain() {
          std::vector<typename MONTE_CARLO_TYPE::Chain*> chains;
          for (unsigned int i=0; i<threads.size(); i++) {
               chains.push_back(threads[i]->get_chain()[0]);
          }
          return chains;
     }


     //! Get monte carlo statistics
     typename MonteCarloStatistics::Summary get_statistics() {
          std::vector<MonteCarloStatistics::Summary> summaries;
          for (unsigned int i=0; i<threads.size(); i++) {
               summaries.push_back(threads[i]->get_statistics());
          }
          return typename MonteCarloStatistics::Summary(summaries);
     }

     //! Get move statistics
     typename MoveStatistics<typename MONTE_CARLO_TYPE::Chain>::MoveStatisticsSummary get_move_statistics() {
          std::vector<typename MoveStatistics<typename MONTE_CARLO_TYPE::Chain>::MoveStatisticsSummary> statistics;
          for (unsigned int i=0; i<threads.size(); i++) {
               statistics.push_back(threads[i]->move_collection->get_statistics());
          }
          return typename MoveStatistics<typename MONTE_CARLO_TYPE::Chain>::MoveStatisticsSummary(statistics);
     }

     //! Get energy data
     typename Energy<typename MONTE_CARLO_TYPE::Chain>::EnergyDataCollection get_energy_data() {
          std::vector<std::vector<typename Energy<typename MONTE_CARLO_TYPE::Chain>::EnergyData> > data;
          for (unsigned int i=0; i<threads.size(); i++) {
               typename Energy<typename MONTE_CARLO_TYPE::Chain>::EnergyDataCollection thread_data = threads[i]->get_energy_data();
               for (unsigned int j=0; j<thread_data.size(); ++j) {
                    if (j>=data.size())
                         data.resize(j+1, std::vector<typename Energy<typename MONTE_CARLO_TYPE::Chain>::EnergyData>(threads.size()));
                    data[j][i] = thread_data[j];
               }
          }
          typename Energy<typename MONTE_CARLO_TYPE::Chain>::EnergyDataCollection data_collection;
          for (unsigned int j=0; j<data.size(); ++j) {
               data_collection.add(typename Energy<typename MONTE_CARLO_TYPE::Chain>::EnergyData(data[j]));
          }
          return data_collection;
     }

     //! Display settings
     virtual void display_settings() {
          for (unsigned int i=0; i<threads.size(); i++) {
               std::cout << "THREAD: " + stringify(i) + "\n";
               threads[i]->display_settings();
          }
     }

};



//! Base class for simulation code
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class Simulation: public MonteCarlo<CHAIN_TYPE> {
public:

     //! Use Settings object from base class
     typedef typename MonteCarlo<CHAIN_TYPE>::Settings Settings;

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > random_generator_uniform_01;

     //! Constructor
     //!
     //! \param id Identification string
     //! \param chain Molecule chain.
     //! \param energy_function Energy object
     //! \param move_collection Move set
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     Simulation(std::string id,
                CHAIN_TYPE *chain,
                Energy<CHAIN_TYPE> *energy_function, 
                MoveCollection<CHAIN_TYPE> *move_collection,
                const Settings &settings=Settings(),
                RandomNumberEngine *random_number_engine= &random_global)
          : MonteCarlo<CHAIN_TYPE>(id, chain, energy_function, move_collection, 
                                   settings,
                                   random_number_engine),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)) {}

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     Simulation(const Simulation &other)
          : MonteCarlo<CHAIN_TYPE>(other),
            random_generator_uniform_01(other.random_generator_uniform_01) {}
            

     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     Simulation(const Simulation &other, 
                RandomNumberEngine *random_number_engine,
                int thread_index)
          : MonteCarlo<CHAIN_TYPE>(other, random_number_engine, thread_index),
            random_generator_uniform_01(*random_number_engine, 
                                        boost::uniform_real<>(0,1)) {}


     //! Determine whether to accept or reject based on proposal and energy biases
     //! Note, that the bias terms should be given as log-probs
     bool acceptance_criterion(double move_bias, double energy_bias) {
	  double tot_bias = move_bias + energy_bias;

	  bool accept = true;

	  if (tot_bias < 0) {
	       double d = random_generator_uniform_01();
	       if (d > std::exp(tot_bias)) {
		    accept = false;
	       }
	  }
	  return accept;
     }

     //! Display settings
     virtual void display_settings() {
          std::cout << "Simulation settings (" + this->id + "):\n";
          std::cout.flush();

          // Create output stream where everything is indented
          boost::iostreams::filtering_ostream out;
          out.push(IndentedStreamFilter(4));
          out.push(std::cout);           
          out << (((DERIVED_CLASS*)(this))->settings);
          out.flush();

          this->move_collection->display_settings();
          this->energy_function->display_settings();
     }     
};


//! Base class for optimization code
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class Optimization: public MonteCarlo<CHAIN_TYPE> {
public:
     //! Use Settings object from base class
     typedef typename MonteCarlo<CHAIN_TYPE>::Settings Settings;


     //! Constructor
     //!
     //! \param id Identification string
     //! \param chain Molecule chain.
     //! \param energy_function Energy object
     //! \param move_collection Move set
     //! \param settings Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     Optimization(std::string id,
                  CHAIN_TYPE *chain,
                  Energy<CHAIN_TYPE> *energy_function, 
                  MoveCollection<CHAIN_TYPE> *move_collection,
                  const Settings &settings=Settings(),
                  RandomNumberEngine *random_number_engine= &random_global):
          MonteCarlo<CHAIN_TYPE>(id, chain, energy_function, move_collection, 
                                 settings,
                                 random_number_engine){}

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     Optimization(const Optimization &other)
          : MonteCarlo<CHAIN_TYPE>(other) {}

     //! Copy constructor. Using different random_number_generator and with specified thread index.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index which thread the copy will run in
     Optimization(const Optimization &other,
                  RandomNumberEngine *random_number_engine,
                  int thread_index)
          : MonteCarlo<CHAIN_TYPE>(other, random_number_engine, thread_index) {}


     //! Determine whether to accept or reject based on proposal and energy biases
     //! Note, that the bias terms should be given as log-probs
     bool acceptance_criterion(double move_bias, double energy_bias) {
	  double tot_bias = move_bias + energy_bias;

	  bool accept = true;

	  if (tot_bias < 0) {
		  accept = false;
	  }
	  return accept;
     }
     
     //! Display settings
     virtual void display_settings() {
          std::cout << "Optimization settings: (" + this->id + ")\n";
          std::cout.flush();

          // Create output stream where everything is indented
          boost::iostreams::filtering_ostream out;
          out.push(IndentedStreamFilter(4));
          out.push(std::cout);           
          out << (((DERIVED_CLASS*)(this))->settings);
          out.flush();

          this->move_collection->display_settings();
          this->energy_function->display_settings();
     }     
};

}

#endif
