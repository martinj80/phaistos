// evaluate_observable.cpp --- Evaluate energies/observables
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


#include "phaistos_option_definitions.h"

#include "protein/xtc_chain.h"

// Import Phaistos namespace
using namespace phaistos;

// Functor which changes the default output-target of observables to stdout
// (users will expect to see their observable be outputted to the screen)
struct SettingsModifierObservable: public SettingsModifier  {
     template <typename SETTINGS_PARENT>
     typename SETTINGS_PARENT::Settings *modify(typename SETTINGS_PARENT::Settings *settings, std::string prefix="") const {
          typename SETTINGS_PARENT::Settings *settings_new = SettingsModifier::modify<SETTINGS_PARENT>(settings, prefix);
          ObservableBase::Settings *observable_settings = dynamic_cast<ObservableBase::Settings*>(settings_new);
          if (observable_settings) {
               observable_settings->output_target = "stdout#compact";
               observable_settings->register_interval = 1;
          }
          return settings_new;
     }
};


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

                    // Add supergroups
                    target.add_super_group("energy");
                    target.add_super_group("observable");
                    target.add_super_group("backbone-dbn");


                    // These options are not displayed in the normal option output
                    ProgramOptionParser::OptionsDescription options_hidden("Command-line-only options");
                    options_hidden.add_options()
                         ("help,h", "help message")
                         ;
                    target.add(options_hidden, true);

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
                         ("verbose", po::value<bool>()->implicit_value(1)->default_value(0),
                          "Output information during run")
                         ("debug", po::value<int>()->implicit_value(1)->default_value(0),
                          "Debug information level")
                         ("data-dir", po::value<std::string>()->default_value(data_dir),
                          "Path to Phaistos data directory")
                         ("chain-type", po::value<std::string>()->default_value("chainfb"),
#ifdef PHAISTOS_INCLUDE_CHAIN_CA
                          "Chain type (chainfb|chainca)")
#else
                          "Chain type (only chainfb is currently activated)")
#endif
                         ("seed", po::value<unsigned int>()->default_value(time(NULL)),
                          "Seed for random number generator. The default value is the current time. Remove this line from config file to use random seed.")
                         ("output-directory", po::value<std::string>()->default_value("annotated_pdbs"),
                          "Output directory for sampled pdb files")
                         ("same-ensemble", po::value<bool>()->implicit_value(true)->default_value(false),
                          "Whether all pdb files are from the same ensemble. When set to true, it is assumed that all pdb files correspond to different structures for the same protein, which means that protein chain can be updated rather than reinitialized in each iteration. This mode also makes it possible to use observables to get statistics for an ensemble.")
                         ;

                    // Input options
                    ProgramOptionParser::OptionsDescription options_input("Input options");
                    options_input.add_options()
                         ("pdb-file", po::value<std::string>()->default_value(""),
                          "Input: pdb filename")
                         ("pdb-files", po::value<std::vector<std::string> >()->multitoken(), 
                          "Input: pdb filenames")
                         ("pdb-list-file", po::value<std::string>()->default_value(""), 
                          "Input: file containing list of pdb filenames")
                         ("trajectory", po::value<std::string>()->default_value(""), 
                          "Input: trajectory file")
                         ;

                    // Register option description types with class
                    target.add(options_general);
                    // target.add(options_mocapy);
                    target.add(options_input);

               } catch(std::exception& e) {
                    std::cerr << "Command line parsing error: " << e.what() << "\n";
                    exit(1);
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
                    exit(1);
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

               // Define shorthand options
               bool hidden = true;
               target.add(target.create_secondary_shorthand("Energy shorthand",
                                                            "energy", "energy"), "shorthands", hidden);
               target.add(target.create_secondary_shorthand("Observable",
                                                            "observable", "observable"), "shorthands", hidden);

               // Options for protein chain
               OptionsDescription options_chain("Chain options");
               options_chain.add_options()
                    ("atom-types", po::value<ProgramOptionParser::WrappedEnum<definitions::AtomSelectionEnum> >()
                     ->multitoken()->default_value(std::string("ALL_PHYSICAL_ATOMS")),
                     "The types of atoms included in the chain")
                    ;
               target.add(options_chain);

               // Specifies current supergroup
               std::string super_group = "";

               try {

//! Module general option definitions are inserted here
#include "modules/phaistos_cpp/options.cpp"

                    ////////////////////
                    // ENERGY OPTIONS //
                    ////////////////////

                    static const unsigned int n_supergroups = 2;
                    std::string supergroups[n_supergroups] = {"energy","observable"};
                    std::string prefixes[n_supergroups] =    {"energy","observable"};

                    // Iterate over the supergroups
                    for (unsigned int i=0; i<n_supergroups; ++i) {

                         super_group = supergroups[i];
                         std::string prefix = prefixes[i];

                         EnergyOptions<SettingsModifierObservable>(target, occurrences, super_group, 
                                                                   prefix, (CHAIN_TYPE*)NULL, (DBN_TYPE*)NULL);

                    }

                    target.super_group_default("observable")   = "--observable @energy-sum @energy-terms";

               } catch(std::exception& e) {
                    std::cerr << "Command line parsing error: " << e.what() << "\n";
                    exit(1);
               } catch(...) {
                    exit(1);
               }
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


//! Core evaluation code
template<typename CHAIN_TYPE, typename DBN_TYPE>
void evaluate(::Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
              Energy<CHAIN_TYPE> *energy, 
              ObservableCollection<CHAIN_TYPE> *observable_collection,
              std::string id,
              bool initialized,
              int iteration_number,
              bool overwrite_existing_files,
              std::string pdb_filename) {

     if (options["verbose"].as<bool>()) {
          std::cout << "DBN: \n" << *dbn << "\n";
          std::cout << "\n\n";
     }

     energy->evaluate();

     observable_collection->observe(id, iteration_number, initialized);

     if (observable_collection->has_stream(ObservableCollection<CHAIN_TYPE>::PDB_HEADER) || 
         observable_collection->has_stream(ObservableCollection<CHAIN_TYPE>::PDB_B_FACTOR)) {

          if (!initialized) {
               // Make the output directory and check if it is writable
               if (system(("mkdir -p " + options["output-directory"].as<std::string>()).c_str())) {
                    std::cerr << "ERROR: Cannot create output-directory: " + 
                         options["output-directory"].as<std::string>() << std::endl;
                    exit(1);
               }
          }

          std::string pdb_filename_short = pdb_filename;
          unsigned int slash_pos = pdb_filename.find_last_of("/");
          if (slash_pos != std::string::npos)
               pdb_filename_short = pdb_filename.substr(slash_pos+1, std::string::npos);

          std::string filename = options["output-directory"].as<std::string>() + "/" + pdb_filename_short;
          std::string observable_pdb_header = observable_collection->
               gather_string_streams(ObservableCollection<CHAIN_TYPE>::PDB_HEADER);
          std::string observable_pdb_b_factor_string = observable_collection->
               gather_string_streams(ObservableCollection<CHAIN_TYPE>::PDB_B_FACTOR);
          chain->output_as_pdb_file(filename, NULL, 
                                    observable_pdb_header,
                                    observable_pdb_b_factor_string);

                    
     }


}


//! Templated main function
template<typename CHAIN_TYPE, typename DBN_TYPE>
void evaluate_observable_main(::Options &options) {

     options.parse<CHAIN_TYPE,DBN_TYPE>();

     // Print options to stdout
     if (options["verbose"].as<bool>())
          std::cout << options << "\n";

     // Initialize random number generators
     RandomNumberEngineCollection random_number_engines;

     CHAIN_TYPE *chain = NULL;
     DBN_TYPE *dbn = NULL;

     bool same_ensemble = options["same-ensemble"].as<bool>();

     std::vector<std::string> pdb_filenames;
     if (options.has_option("pdb-files"))
          pdb_filenames = options["pdb-files"].as<std::vector<std::string> >();
     else if (options.has_option("pdb-list-file") && options["pdb-list-file"].as<std::string>()!="")
          pdb_filenames = file_to_string_vector(options["pdb-list-file"].as<std::string>());

     if (options.has_option("pdb-file") && options["pdb-file"].as<std::string>()!="")
          pdb_filenames.push_back(options["pdb-file"].as<std::string>());

     std::string trajectory_filename = "";
     if (options.has_option("trajectory") && options["trajectory"].as<std::string>()!="") {
          trajectory_filename = options["trajectory"].as<std::string>();
          same_ensemble = true;
          if (pdb_filenames.size() != 1) {
               std::cerr << "Using trajectory file: single PDB file expected (for topology)";
               exit(1);
          }
     }
         

     if (pdb_filenames.empty())
          pdb_filenames.push_back("");

     Energy<CHAIN_TYPE> *energy = NULL;
     ObservableCollection<CHAIN_TYPE> *observable_collection = NULL;

     unsigned int counter = 0;
     for (counter=0; counter<pdb_filenames.size(); ++counter) {
          
          std::string pdb_filename = pdb_filenames[counter];

          // Check if weight was specified in pdb file list
          double weight = 1.0;
          boost::trim(pdb_filename);
          std::vector<std::string> tokens;
          boost::split(tokens, pdb_filename, boost::is_any_of(" \t"), boost::token_compress_on);
          if (tokens.size() > 1) {
               pdb_filename = tokens[0];
               weight = boost::lexical_cast<double>(tokens[1].c_str());
          }

          if (pdb_filename == "")
               continue;

          // Whether to overwrite existing output files
          bool overwrite_existing_files = (counter>0);

          // If independent observables is set, allocate new objects in each iteration
          if (observable_collection==NULL || !same_ensemble) {

               if (dbn != NULL) {
                    delete dbn;
                    dbn = NULL;
               }
               initialize_dbn(options, &dbn, pdb_filename, &(*random_number_engines));

               if (chain != NULL) {
                    delete chain;
                    chain = NULL;
               }
               initialize_chain(options, &chain, dbn, pdb_filename);

               if (energy != NULL)
                    delete energy;
               // Energy object. Terms added to the energy will not be reported
               // directly - they are instead referred to using the @energy-sum
               // and @energy-terms special observable names
               energy = new Energy<CHAIN_TYPE>(chain);

               // Initialize energy object
               initialize_energy(options, chain, dbn, energy, &(*random_number_engines));

               // Create a collection of observables
               if (observable_collection)
                    delete observable_collection;
               observable_collection = new ObservableCollection<CHAIN_TYPE>(chain, overwrite_existing_files);
               initialize_energy(options, chain, dbn, observable_collection, &(*random_number_engines), "observable", energy); 

          } 
          // Otherwise, update chain with data from new pdb file
          else {
               CHAIN_TYPE *chain_pdb = NULL;
               initialize_chain(options, &chain_pdb, dbn, pdb_filename);
               MoveCollection<CHAIN_TYPE> move_collection(chain);
               MoveFixedStructure<CHAIN_TYPE> *move_fixed_structure = new MoveFixedStructure<CHAIN_TYPE>(chain, *chain_pdb);
               move_collection.set_initializer(move_fixed_structure);
               move_collection.reinitialize();
               delete chain_pdb;

          }

          for (unsigned int j=0; j<energy->terms.size(); ++j) {
               energy->terms[j]->extra_weight = weight;
          }
          for (unsigned int j=0; j<observable_collection->terms.size(); ++j) {
               observable_collection->terms[j]->extra_weight = weight;
          }

          bool initialized = (counter!=0);

          std::string id = pdb_filename;
          if (id=="")
               id = dbn->settings.initial_pdb_file;
          if (id=="")
               id = "---";
          evaluate(options, chain, dbn, energy, observable_collection, id, 
                   initialized, counter+1, overwrite_existing_files, pdb_filename);
     }

     if (trajectory_filename != "") {
          // Check consistency between chain and trajectory
          xtc_chain_check_consistency(*chain, trajectory_filename);

          int step;
          float time;
          float prec;
          XDRFILE *trajectory_file = xdrfile_open(trajectory_filename.c_str(),"r");
          while (xtc_read_chain(*chain, trajectory_file, step, time, prec)) {

               // Whether to overwrite existing output files
               bool overwrite_existing_files = (counter > 0);

               bool initialized = (counter!=0);

               std::string id = "trajectory_" + boost::lexical_cast<std::string>(step);
               evaluate(options, chain, dbn, energy, observable_collection, id, 
                        initialized, counter+1, overwrite_existing_files, trajectory_filename);
               
               ++counter;
          }
          xdrfile_close(trajectory_file);          
     }

     delete chain;
     delete dbn;
     delete energy;
     delete observable_collection;
}



//! Main function
int main(int argc, char *argv[]) {

     // Read in command line 
     ::Options options(argc, argv);

     // Call main function either in Torus or FB5 mode
     if (options["chain-type"].as<std::string>()=="chainfb") {

          if (options["backbone-dbn-torus"].occurrences()) {
               evaluate_observable_main<ChainFB,TorusDBN>(options);
          } else if (options["backbone-dbn-torus-cs"].occurrences()) {
               evaluate_observable_main<ChainFB,TorusCsDBN>(options);
#ifdef PHAISTOS_INCLUDE_OMEGA_DBNS
          } else if (options["backbone-dbn-torus-omega"].occurrences()) {
               evaluate_observable_main<ChainFB,TorusOmegaDBN>(options);
          } else if (options["backbone-dbn-torus-cs-omega"].occurrences()) {
               evaluate_observable_main<ChainFB,TorusCsOmegaDBN>(options);
#endif
          } else {
               std::cerr << "No backbone-dbn model specified. Aborting\n";
          }

#ifdef PHAISTOS_INCLUDE_CHAIN_CA
     } else {
          if (options["backbone-dbn-fb5"].occurrences()) {
               evaluate_observable_main<ChainCA,FB5DBN>(options);
          } else {
               std::cerr << "No backbone-dbn model specified. Aborting\n";
          }
#endif
     }

}
