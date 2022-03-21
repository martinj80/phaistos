// phaistos.cpp --- Main program for the phaistos package. Runs a protein folding simulation
// Copyright (C) 2008-2010 Wouter Boomsma, Thomas Hamelryck, Mikael Borg, Jes Frellsen, Tim Harder, Kasper Stovgaard, Martin Paluszewski
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


#include "phaistos_option_definitions.h"

#include "revision.h"

//! \mainpage Phaistos
//! Phaistos is a Markov chain Monte Carlo framework for protein
//! simulations. For details, see http://www.phaistos.org

// Import Phaistos namespace
using namespace phaistos;

//! Option class. Uses boost::program_options.
//! Interface to command line options and optionally
//! a container to Settings objects
class Options: public phaistos::Options {
public:

     //! Functor containing primary option definitions - these will be parsed first.
     //! These are template-unspecific options.
     struct DefineOptionsPrimary {

          //! Functor evaluate function
          //!
          //! \param target Underlying ProgramOptionParser object to which options are added
          //! \param occurrences The number of occurences of each option (this is not used for the Primary Options)
          void operator()(ProgramOptionParser &target,
                          const ProgramOptionParser::Filter &occurrences=ProgramOptionParser::Filter(1)) const {

               using namespace boost::fusion;

               try {

                    bool hidden = true;

                    // Add supergroups
                    target.add_super_group("procedure");
                    target.add_super_group("monte carlo");
                    target.add_super_group("move");
                    target.add_super_group("energy");
                    target.add_super_group("secondary energy");
                    target.add_super_group("observable");
                    target.add_super_group("backbone-dbn");


                    // These options are not displayed in the normal option output
                    ProgramOptionParser::OptionsDescription options_hidden("Command-line-only options");
                    options_hidden.add_options()
                         ("help,h", "help message")
                         ("config-file", po::value<std::string>()->default_value(""),
                          "configuration file")
                         ;
                    target.add(options_hidden, hidden);


                    std::string data_dir = "../data";
                    // Retrieve data-dir from environment variable if available
                    char *data_dir_env = getenv("PHAISTOS_DATA_DIR");
                    char *root_env = getenv("PHAISTOS_ROOT");
                    if (data_dir_env) {
                         data_dir = std::string(data_dir_env);
                    } else if (root_env) {
                         data_dir = std::string(root_env)+"/data";
                    }

                    // General options
                    ProgramOptionParser::OptionsDescription options_general("General options");
                    options_general.add_options()
                         ("verbose", po::value<bool>()->implicit_value(1)->default_value(1),
                          "Output information during run")
                         ("debug", po::value<int>()->implicit_value(1)->default_value(0),
                          "Debug information level")
                         ("data-dir", po::value<std::string>()->default_value(data_dir),
                          "Path to Phaistos data directory")
                         ("status-interval", po::value<int>()->default_value(5),
                          "Frequency of status display")
                         ("chain-type", po::value<std::string>()->default_value("chainfb"),
#ifdef PHAISTOS_INCLUDE_CHAIN_CA
                          "Chain type (chainfb|chainca)")
#else
                          "Chain type (only chainfb is currently activated)")
#endif
                         ("seed", po::value<unsigned int>()->default_value(time(NULL)),
                          "Seed for random number generator. The default value is the current time. Remove this line from config file to use random seed.")
                         //  Timeout after M minutes
                         ("timeout-minutes", po::value<PHAISTOS_LONG_LONG>()->default_value(0),
                          "Timeout after M minutes of execution")
                         ("timeout-time", po::value<PHAISTOS_LONG_LONG>()->default_value(0),
                          "Timeout: exits at time T (seconds from Epoch)")
                         ("temperature", po::value<double>(),
                          "Temperature (used by non-probabilistic energy terms)")
                         ;

                    // Register option description types with class
                    target.add(options_general);

                    // General options
                    ProgramOptionParser::OptionsDescription options_mode("Phaistos mode (predefined combinations of options)");
                    options_mode.add_options()
                         ("mode", po::value<std::string>()->default_value(""),
                          ("Phaistos Mode (" + boost::algorithm::join(phaistos_modes, "|")+")").c_str())
                         ;
                    // Register option description types with class
                    target.add(options_mode, hidden);

                    // Input options
                    ProgramOptionParser::OptionsDescription options_input("Input options");
                    options_input.add_options()
                         ("pdb-file", po::value<std::string>()->default_value(""),
                          "Input: pdb filename")
                         ("aa-file", po::value<std::string>()->default_value(""),
                          "Input: amino acid residue sequence filename")
                         ("ss-file", po::value<std::string>()->default_value(""),
                          "Input: secondary structure sequence filename")
                         ;

                    // target.add(options_mocapy);
                    target.add(options_input);

               } catch(std::exception& e) {
                    std::cerr << "Command line parsing error: " << e.what() << "\n";
                    assert(false);
               } catch(...) {
                    exit(1);
               }
          }
     };

     //! Functor containing secondary option definitions - these will be parsed 
     //! either together with the primary or with the tertiary options.
     struct DefineOptionsSecondary {

          //! Functor evaluate function
          //!
          //! \param target Underlying ProgramOptionParser object to which options are added
          //! \param occurrences The number of occurences of each option - to allow multiply definitions of the same option
          void operator()(ProgramOptionParser &target,
                          const ProgramOptionParser::Filter &occurrences=ProgramOptionParser::Filter(1)) const {

               using namespace boost::fusion;

               try {
                    
                    bool multiple_occurrences_allowed = false;
                    bool hidden = true;
                    target.add(target.create_secondary_shorthand("Shorthand for specification of BackboneDBN model (torus|fb5|torus-cs)",
                                                                 "backbone-dbn", "backbone-dbn",
                                                                 multiple_occurrences_allowed), "shorthands", hidden);

                    std::string super_group = "backbone-dbn";
                    BackboneDBNOptions(target, occurrences, super_group);

               } catch(std::exception& e) {
                    std::cerr << "Command line parsing error: " << e.what() << "\n";
                    assert(false);
               } catch(...) {
                    exit(1);
               }
          }
     };

     //! Functor containing tertiary option definitions - these will be parsed when chain type is known
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     struct DefineOptionsTertiary {

          //! Local typedef for OptionsDescription          
          typedef ProgramOptionParser::OptionsDescription OptionsDescription;

          //! Local typedef for CommandLine
          typedef ProgramOptionParser::CommandLine CommandLine;

          //! Local typedef for OptionShorthand
          typedef ProgramOptionParser::OptionShorthand OptionShorthand;

          //! Functor evaluate function
          //!
          //! \param target Underlying ProgramOptionParser object to which options are added
          //! \param occurrences The number of occurences of each option - to allow multiply definitions of the same option
          void operator()(ProgramOptionParser &target,
                          const ProgramOptionParser::Filter &occurrences=ProgramOptionParser::Filter(1)) const {

               using namespace boost::fusion;

               target.super_group_default("monte carlo")   = "--monte-carlo metropolis-hastings";
               target.super_group_default("procedure")   = "--procedure fold";

               // Set mc-mode description based on availability of libraries
               std::vector<std::string> mc_mode_description;


               // Define shorthand options
               bool hidden = true;
               target.add(target.create_secondary_shorthand("Energy shorthand",
                                                            "energy", "energy"), "shorthands", hidden);
               target.add(target.create_secondary_shorthand("Energy shorthand (secondary energy)",
                                                            "energy2", "energy2"), "shorthands", hidden);
               target.add(target.create_secondary_shorthand("Observable",
                                                            "observable", "observable"), "shorthands", hidden);
               target.add(target.create_secondary_shorthand("Move shorthand",
                                                            "move", "move"), "shorthands", hidden);
               bool multiple_occurrences_allowed = false;
               target.add(target.create_secondary_shorthand(
                               std::string("Shorthand for specification of Phaistos procedure") +
                               std::string("(fold|compactify|comparison)"),
                               "procedure", "procedure",
                               multiple_occurrences_allowed), "shorthands", hidden);

               // Define defaults for the different modes
               ModeDefinitions mode_definitions(target, (CHAIN_TYPE*)NULL);

               // Options for protein chain
               OptionsDescription options_chain("Chain options");
               options_chain.add_options()
                    ("init-from-pdb", po::value<bool>()->implicit_value(1)->default_value(0),
                     "Whether to initialize the chain from a pdb-file")
                    ("atom-types", po::value<ProgramOptionParser::WrappedEnum<definitions::AtomSelectionEnum> >()
                     ->multitoken()->default_value(mode_definitions.atom_types_default),
                     "The types of atoms included in the chain")
                    ;
               target.add(options_chain);


               // Specifies current supergroup
               std::string super_group = "";

               try {

                    super_group = "procedure";
                    ProcedureOptions(target, occurrences, super_group, (CHAIN_TYPE*)NULL, (DBN_TYPE*)NULL);

                    super_group = "monte carlo";
                    MonteCarloOptions(target, occurrences, super_group, mc_mode_description, (CHAIN_TYPE*)NULL, (DBN_TYPE*)NULL);

//! Module general option definitions are inserted here
#include "modules/phaistos_cpp/monte_carlo_options.cpp"

//! Module general option definitions are inserted here
#include "modules/phaistos_cpp/options.cpp"


                    //////////////////
                    // MOVE OPTIONS //
                    //////////////////
                    super_group = "move";
                    MoveOptions<SettingsModifier>(target, occurrences, super_group, (CHAIN_TYPE*)NULL, (DBN_TYPE*)NULL);


                    ////////////////////
                    // ENERGY OPTIONS //
                    ////////////////////

                    static const unsigned int n_supergroups = 3;
                    std::string supergroups[n_supergroups] = {"energy", "secondary energy", "observable"};
                    std::string prefixes[n_supergroups] =    {"energy", "energy2", "observable"};

                    // Iterate over the supergroups
                    for (unsigned int i=0; i<n_supergroups; ++i) {

                         super_group = supergroups[i];
                         std::string prefix = prefixes[i];

                         EnergyOptions<SettingsModifier>(target, occurrences, super_group, prefix, 
                                                         (CHAIN_TYPE*)NULL, (DBN_TYPE*)NULL);

                    }

               } catch(std::exception& e) {
                    std::cerr << "Command line parsing error: " << e.what() << "\n";
                    assert(false);
               } catch(...) {
                    exit(1);
               }

               multiple_occurrences_allowed = false;
               target.add(target.create_secondary_shorthand("Shorthand for specification of Monte Carlo mode (" +
                                                            boost::algorithm::join(mc_mode_description, "|") + ")",
                                                            "monte-carlo", "monte-carlo",
                                                            multiple_occurrences_allowed), "shorthands", hidden);


          }

     };

     //! Main parse function
     //!
     //! \tparam CHAIN_TYPE Molecule chain type
     //! \tparam DBN_TYPE Dynamic Bayesian Network type
     //! \param include_first_primary Whether to include primary parse of first parser object
     //! \param reparsing Whether we are currently doing a re-parse.
     template <typename CHAIN_TYPE,typename DBN_TYPE>
     void parse(bool include_first_primary=false, bool reparsing=false) {
          phaistos::Options::parse<CHAIN_TYPE,DBN_TYPE,
                             DefineOptionsPrimary, 
                             DefineOptionsSecondary, 
                             DefineOptionsTertiary<CHAIN_TYPE,DBN_TYPE> >(include_first_primary, reparsing);

          if (parser_5->check_all_shorthand_uninitialized_in_super_group("procedure")) {
               std::cout << "Error: no procedure specified. \n";
               std::cout << "Please select a procedure using the --procedure option or the --mode option, \nor use --help for a list of all available options.\n\n";
               exit(1);
          }
     }


     //! Constructor
     //!
     //! \param argc Number of command line arguments
     //! \param argv Array of arguments
     Options(int argc, char* argv[])
          : phaistos::Options(argc, argv) {

          phaistos::Options::parse_initial<DefineOptionsPrimary, 
                                           DefineOptionsSecondary>();
     }
};





//! Execute a procedure=fold simulation
//!
//! \param options Command line options object
//! \param dbn Dynamic Bayesian Network object
//! \param chain Molecule chain object
//! \param monte_carlo MonteCarlo object
//! \param energy_function Energy object
//! \param energy_function_secondary Secondary energy object
//! \param observable_collection Collection of observables
//! \param move_collection Move set
//! \param random_number_engines Vector of random number engines (one for each thread)
template<typename CHAIN_TYPE, typename DBN_TYPE>
void phaistos_fold(::Options &options,
                   DBN_TYPE *dbn, CHAIN_TYPE *chain,
                   MonteCarloBase<CHAIN_TYPE> *monte_carlo,
                   Energy<CHAIN_TYPE> *energy_function,
                   Energy<CHAIN_TYPE> *energy_function_secondary,
                   ObservableCollection<CHAIN_TYPE> *observable_collection,
                   MoveCollection<CHAIN_TYPE> *move_collection,
                   std::vector<RandomNumberEngine *> *random_number_engines) {

     // Cache option values (for efficiency)
     PHAISTOS_LONG_LONG iterations = options["iterations"].as<PHAISTOS_LONG_LONG>();
     PHAISTOS_LONG_LONG timeout_minutes = options["timeout-minutes"].as<PHAISTOS_LONG_LONG>();
     PHAISTOS_LONG_LONG timeout_time = options["timeout-time"].as<PHAISTOS_LONG_LONG>();
     int status_interval = options["status-interval"].as<int>();
     int threads = options["threads"].as<int>();
     int steps_per_move = options["steps-per-move"].as<int>();
     bool verbose = options["verbose"].as<bool>();

     std::vector<Energy<CHAIN_TYPE> *> energy_function_secondary_vector(threads);
     std::vector<ObservableCollection<CHAIN_TYPE> *> observable_collection_vector(threads);
     std::vector<CHAIN_TYPE*> chains = monte_carlo->get_chain();
     if (energy_function_secondary) {
          for (int i=0; i<threads; ++i) {
               energy_function_secondary_vector[i] = new Energy<CHAIN_TYPE>(*energy_function_secondary, chains[i], (*random_number_engines)[i], i);
          }
     }
     if (observable_collection) {
          for (int i=0; i<threads; ++i) {
               observable_collection_vector[i] = new ObservableCollection<CHAIN_TYPE>(*observable_collection, 
                                                                                      chains[i], (*random_number_engines)[i], i,
                                                                                      monte_carlo->get_energy_function()[i]);
          }
     }

     //! Functor executed at the end of each step in the monte carlo run
     class StepCode: public MonteCarloStepCode<CHAIN_TYPE> {
     public:

          // std::vector<Energy<CHAIN_TYPE> *> &energy_function_secondary;
          std::vector<ObservableCollection<CHAIN_TYPE> *> &observable_collection;

          //! Constructor
          StepCode(std::vector<ObservableCollection<CHAIN_TYPE> *> &observable_collection)
               : observable_collection(observable_collection) {}

          //! Code to execute in each step
          //! Note that this code might be executed in parallel, so modifications of 
          //! shared variables must be protected by a mutex
          void operator()(MonteCarlo<CHAIN_TYPE> *monte_carlo) {

               double energy_value = monte_carlo->energy_current;


               std::string id = generate_log_filename("sample_%i_%t_%e",
                                                      energy_value,
                                                      monte_carlo->iteration_counter,
                                                      monte_carlo->thread_index);
               observable_collection[monte_carlo->thread_index]->observe(id, monte_carlo->iteration_counter);

// Module-specific code
#include "modules/phaistos_cpp/fold_step_code.cpp"


          }
     // } step_code(options, chain, energy_function_secondary_vector, observable_collection_vector);
     } step_code(observable_collection_vector);

     // Create output stream where everything is indented
     boost::iostreams::filtering_ostream out_8;
     out_8.push(IndentedStreamFilter(8));
     out_8.push(std::cout);
     boost::iostreams::filtering_ostream out_4;
     out_4.push(IndentedStreamFilter(4));
     out_4.push(std::cout);

     // Execute the first move manually to ensure that a reinitialize (and potentially
     // clash removal) is triggered before we dump the initial PDB structure.
     monte_carlo->move(1, &step_code);

     // Dump initial conformation for each thread
     if (threads > 1) {
          for (int i=0; i<threads; ++i) {
               chains[i]->output_as_pdb_file("initial_" + boost::lexical_cast<std::string>(i) + ".pdb", chains[i]);
          }
     } else {
          chains[0]->output_as_pdb_file("initial.pdb", chains[0]);
     }

     // Main loop
     time_t start_time_seconds=time(NULL);
     bool loop_finish_flag=false;
     while ( not(loop_finish_flag) && (monte_carlo->iteration_counter < iterations) ) {

          // Make MCMC move using move collection
          monte_carlo->move(1, &step_code);

          //// Add code here for non-parallel executation ////
          if (timeout_minutes > 0) {
               time_t elapsed_time_seconds=time(NULL)-start_time_seconds;
               if ( elapsed_time_seconds/60 >= static_cast<time_t>(timeout_minutes)) {
                    std::cout << "# TIMEOUT: " << elapsed_time_seconds/60 << " minutes elapsed (" << timeout_minutes << " requested). Exiting\n";
                    loop_finish_flag=true;
               }
          }
          if (timeout_time > 0) {
               time_t current_time_seconds=time(NULL);
               if ( current_time_seconds >= static_cast<time_t>(timeout_time)) {
                    printf("# TIMEOUT: time limit reached (requested stop at time %.1fs). Exiting\n",(float)current_time_seconds);
                    loop_finish_flag=true;
               }
          }
          // Output energy and move statistics
          if ( (loop_finish_flag) or (monte_carlo->iteration_counter%(status_interval*steps_per_move)==0) ) {
               if (verbose) {
                    std::string iteration_string = boost::lexical_cast<std::string>(monte_carlo->iteration_counter*threads);
                    if (threads>1)
                         iteration_string += " " + boost::lexical_cast<std::string>(std::vector<PHAISTOS_LONG_LONG>(threads, monte_carlo->iteration_counter));
                    std::cout << "\n# Iteration: " << iteration_string << "\n";
                    std::cout << "    Monte Carlo:\n";
                    std::cout.flush();
                    out_8 << monte_carlo->get_statistics() << "\n";
                    out_8.flush();
                    out_4 << monte_carlo->get_energy_data() << "\n";
                    out_4.flush();
                    if (energy_function_secondary &&
                        options["procedure-fold-energy2-evaluation-interval"].as<int>()!=0 && 
                        (monte_carlo->iteration_counter%(options["procedure-fold-energy2-evaluation-interval"].as<int>()
                                                         *options["steps-per-move"].as<int>())==0)) {
                         // Evaluate secondary energy function
                         energy_function_secondary->evaluate();
                         std::cout << "    Energy2:\n" ;
                         std::cout.flush();
                         out_8 << energy_function_secondary->get_data() << "\n";
                         out_8.flush();
                    }
                    std::cout << "    Move statistics:\n" ;
                    std::cout.flush();
                    out_8 << monte_carlo->get_move_statistics() << "\n";
                    out_8.flush();
               }
          }
     }

     for (unsigned int i=0; i<energy_function_secondary_vector.size(); ++i) {
          delete energy_function_secondary_vector[i];
     }
     for (unsigned int i=0; i<observable_collection_vector.size(); ++i) {
          delete observable_collection_vector[i];
     }
}


//! Execute a procedure=comparison simulation.
//! Compares the energies calculated by energy_function 
//! and energy_function_secondary
//!
//! \param options Command line options object
//! \param dbn Dynamic Bayesian Network object
//! \param chain Molecule chain object
//! \param monte_carlo MonteCarlo object
//! \param energy_function Energy object
//! \param energy_function_secondary Secondary energy object
//! \param move_collection Move set
template<typename CHAIN_TYPE, typename DBN_TYPE>
void phaistos_comparison(::Options &options,
                         DBN_TYPE *dbn, CHAIN_TYPE *chain,
                         MonteCarloBase<CHAIN_TYPE> *monte_carlo,
                         Energy<CHAIN_TYPE> *energy_function,
                         Energy<CHAIN_TYPE> *energy_function_secondary,
                         MoveCollection<CHAIN_TYPE> *move_collection) {

     //! Functor executed at the end of each step in the monte carlo run
     class StepCode: public MonteCarloStepCode<CHAIN_TYPE> {
     public:

          ::Options &options;
          Energy<CHAIN_TYPE> *energy_function_secondary;

          PHAISTOS_LONG_LONG n;
          double sx;
          double sy;
          double sxx;
          double syy;
          double sxy;

          PHAISTOS_LONG_LONG identical;

          //! Constructor
          StepCode(::Options &options, Energy<CHAIN_TYPE> *energy_function_secondary)
               : options(options),
                 energy_function_secondary(energy_function_secondary),
                 n(0),sx(0.0),sy(0.0),sxx(0.0),syy(0.0),sxy(0.0),identical(0) {}


          //! Code to execute in each step
          //! Note that this code might be executed in parallel, so modifications of 
          //! shared variables must be protected by a mutex
          void operator()(MonteCarlo<CHAIN_TYPE> *monte_carlo) {

               double energy_value_simulation = monte_carlo->energy_current;

               // Accept or reject based on what happened in the simulation
               if (monte_carlo->move_success) {

                    // Evaluate secondary energy function
                    double energy_value_secondary = energy_function_secondary->evaluate();
                    energy_function_secondary->accept();

                    // Update statistics
                    n++;
                    sx += energy_value_simulation;
                    sy += energy_value_secondary;
                    sxx += energy_value_simulation*energy_value_simulation;
                    syy += energy_value_secondary*energy_value_secondary;
                    sxy += energy_value_simulation*energy_value_secondary;

                    if (options["procedure-comparison-debug"].template as<int>() > 0)
                         std::cout << n << "\t" << energy_value_simulation << "\t" <<  energy_value_secondary << "\n";

                    if (fabs(energy_value_simulation-energy_value_secondary)<
                        (0.001*std::max(fabs(energy_value_simulation),fabs(energy_value_secondary)))) {
                         identical++;
                    }
               } else {
                    energy_function_secondary->reject();
               }
          }
     } step_code(options, energy_function_secondary);


     // Main loop
     while (monte_carlo->iteration_counter < options["iterations"].as<PHAISTOS_LONG_LONG>()) {

          // Make MCMC move using move collection
          monte_carlo->move(1, &step_code);


          //// Add code here for non-parallel executation ////
          double n = step_code.n;
          double sx = step_code.sx;
          double sy = step_code.sy;
          double sxx = step_code.sxx;
          double syy = step_code.syy;
          double sxy = step_code.sxy;
          double b = (n*sxy - sx*sy)/(n*sxx - sx*sx);
          double b_prime = (n*sxy - sx*sy)/(n*syy - sy*sy);
          double correlation = b*b_prime;
          double identical_percentage = 100*step_code.identical/double(n);

          //std::cout << n << " " << sx << " " <<  sy << " " << sxx << " " << sxy << "\n";

          // Output energy and move statistics
          if (monte_carlo->iteration_counter%(options["status-interval"].as<int>()*options["steps-per-move"].as<int>())==0) {
               if (options["verbose"].as<bool>()) {
                    std::cout << "\n# Iteration: " << monte_carlo->iteration_counter << ". Correlation: " << std::setprecision(10) << correlation;
                    if (correlation > 0.99) {
                         std::cout << " . Percentage of matches within 0.1% deviation: " << identical_percentage;
                    }
                    std::cout << "\n";
               }
          }

     }
}

//! Templated main function
template<typename CHAIN_TYPE, typename DBN_TYPE>
void phaistos_main(::Options &options) {

     options.parse<CHAIN_TYPE,DBN_TYPE>();

     // Print options to stdout
     if (options["verbose"].as<bool>())
          std::cout << options << "\n";

     // Initialize random number generators
     RandomNumberEngineCollection random_number_engines(options["threads"].as<int>(),
                                                        options["identical-threads"].as<bool>());


     CHAIN_TYPE *chain = NULL;
     DBN_TYPE *dbn = NULL;
     MonteCarloBase<CHAIN_TYPE> *monte_carlo = NULL;


//! Module initialization code is inserted here
#include "modules/phaistos_cpp/main_initialize.cpp"


     std::string pdb_filename = "";
     if (options.has_option("pdb-file") && options["pdb-file"].as<std::string>()!="")
          pdb_filename = options["pdb-file"].as<std::string>();

     if (pdb_filename != "" && !file_exists(pdb_filename)) {
          std::cerr << "Error: Cannot find " << pdb_filename << " PDB file. Aborting.\n";
          exit(1);
     }

     initialize_dbn(options, &dbn, pdb_filename, &(*random_number_engines));

     initialize_chain(options, &chain, dbn, pdb_filename);

     // Reseed random number generators 
     // The threads are out of sync because dbn was used in the first thread only to 
     // initialize chain. This is only important for debugging purposes with identical
     // random number generators in all seeds (identical-threads=true). 
     random_number_engines.reseed(options["identical-threads"].as<bool>());


     chain->check_consistency();


     ///////////////////
     ////// MOVES //////
     ///////////////////

     // Create move collection. Any move added to this collection
     // will be applied during simulation (using the specified weight).
     MoveCollection<CHAIN_TYPE> move_collection(chain);

     initialize_moves(options, chain, dbn, &move_collection, &(*random_number_engines));


     ////////////////////
     ////// ENERGY //////
     ////////////////////

     // Energy object. Any term added to the energy will be applied
     // during simulation. The add_term method default to a weight of
     // 1, meaning that all terms are weighted equally.
     Energy<CHAIN_TYPE> *energy = new Energy<CHAIN_TYPE>(chain);

     // Initialize energy object
     initialize_energy(options, chain, dbn, energy, &(*random_number_engines));


     // In case of mode==comparison, initialize another energy term 
     // for comparison
     Energy<CHAIN_TYPE> *energy_secondary=NULL;
     if (!options.check_super_group_uninitialized("secondary energy")) {
          energy_secondary = new Energy<CHAIN_TYPE>(chain);
          initialize_energy(options, chain, dbn, energy_secondary,&(*random_number_engines), "energy2");
     }

     // Create a collection of observables
     ObservableCollection<CHAIN_TYPE> *observable_collection = new ObservableCollection<CHAIN_TYPE>(chain);
     initialize_energy(options, chain, dbn, observable_collection, &(*random_number_engines), "observable", energy);

     // If starting in native, output Energy of native state
     if (options["init-from-pdb"].as<bool>() &&
         options["pdb-file"].as<std::string>()!="") {

          definitions::AtomSelectionEnum atom_selection = options["atom-types"].as<definitions::AtomSelectionEnum>();
          CHAIN_TYPE native(options["pdb-file"].as<std::string>(),
                            atom_selection);
          native.check_consistency();
          native.add_atoms(atom_selection);
          native.check_consistency();
          native.output_as_pdb_file("native.pdb", &native);
          Energy<CHAIN_TYPE> energy_native(*energy, &native);
          energy_native.evaluate();

          if (options["verbose"].as<bool>()) {
               std::cout << "NATIVE ENERGY\n" << energy_native << "\n";
          }
     }

     // Reseed random number generators 
     // The threads can be out of sync if constructors in moves or energy terms
     // have sampled random values. We therefore resync everything again here. 
     random_number_engines.reseed(options["identical-threads"].as<bool>());

     // Initialize Monte Carlo object
     initialize_monte_carlo(options, &monte_carlo, energy, energy_secondary, &move_collection, &*random_number_engines);

     // Display settings
     if (options["verbose"].as<bool>()) {
          monte_carlo->display_settings();
          std::cout << "DBN: \n" << *dbn << "\n";
     }
     std::cout << "\n\n";


     // Select simulation function depending on Phaistos mode
     if (options["procedure-fold"].occurrences()) {

          // Run folding procedure
          phaistos_fold(options, dbn, chain, monte_carlo, energy, energy_secondary, observable_collection, &move_collection, &(*random_number_engines));

     } else if (options["procedure-comparison"].occurrences()) {

          // Run comparison procedure
          phaistos_comparison(options, dbn, chain, monte_carlo, energy, energy_secondary, &move_collection);

     } else {

          std::cerr << "ERROR: Unknown procedure: " << options["procedure"].as<std::string>() << "\n";

     }


     // Cleanup
#include "modules/phaistos_cpp/main_finalize.cpp"

     delete chain;
     delete dbn;
     delete monte_carlo;
     delete energy;
     if (energy_secondary) {
          delete energy_secondary;
     }
     delete observable_collection;
}


//! Main function
int main(int argc, char *argv[]) {

     printf("\n");
     printf("############################################################\n");
     printf("#                 PHAISTOS v%s (rev. %s)             #\n", PHAISTOS_VERSION, SVN_REVISION);
     printf("#    A Markov Chain Monte Carlo Simulation Framework       #\n");
     printf("#                                                          #\n");
     printf("#  Please cite:                                            #\n");
     printf("#   Boomsma, Frellsen, Harder, Bottaro, Johansson, Tian    #\n");
     printf("#   Stovgaard, Andreetta, Olsson, Valentin, Antonov,       #\n");
     printf("#   Christensen, Borg, Jensen, Lindorff-Larsen,            #\n");
     printf("#   Ferkinghoff-Borg, Hamelryck, J Comput Chem. 2013       #\n");
     printf("#   doi: 10.1002/jcc.23292                                 #\n");
     printf("#                                                          #\n");
     printf("############################################################\n\n");

     // Read in command line 
     ::Options options(argc, argv);

     // Call main function either in Torus or FB5 mode
     if (options["chain-type"].as<std::string>()=="chainfb") {

          if (options["backbone-dbn-torus"].occurrences()) {
               phaistos_main<ChainFB,TorusDBN>(options);
          } else if (options["backbone-dbn-torus-cs"].occurrences()) {
               phaistos_main<ChainFB,TorusCsDBN>(options);
#ifdef PHAISTOS_INCLUDE_OMEGA_DBNS
          } else if (options["backbone-dbn-torus-omega"].occurrences()) {
               phaistos_main<ChainFB,TorusOmegaDBN>(options);
          } else if (options["backbone-dbn-torus-cs-omega"].occurrences()) {
               phaistos_main<ChainFB,TorusCsOmegaDBN>(options);
#endif
          } else {
               std::cerr << "No backbone-dbn model specified. Aborting\n";
          }

#ifdef PHAISTOS_INCLUDE_CHAIN_CA
     } else {
          if (options["backbone-dbn-fb5"].occurrences()) {
               phaistos_main<ChainCA,FB5DBN>(options);
          } else {
               std::cerr << "No backbone-dbn model specified. Aborting\n";
          }
#endif
     }

     return 0;
}
