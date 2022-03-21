// likelihood.cpp --- Calculate likelihood of structure
// Copyright (C) 2006-2010 Wouter Boomsma
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


// Disable module support in parser
#define MODULE_SUPPORT 0

#include "phaistos_option_definitions.h"

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

                    // Add supergroups
                    target.add_super_group("backbone-dbn");

                    target.super_group_default("backbone-dbn")   = "--backbone-dbn torus";

                    // These options are not displayed in the normal option output
                    ProgramOptionParser::OptionsDescription options_hidden("Command-line-only options");
                    options_hidden.add_options()
                         ("help,h", "help message")
                         ("config-file", po::value<std::string>()->default_value(""),
                          "configuration file")
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
                         ("seed", po::value<unsigned int>()->default_value(time(NULL)),
                          "Seed for random number generator. The default value is the current time. Remove this line from config file to use random seed.")
                         ("start-index", po::value<int>()->default_value(0),
                          "At which sequence position to start")
                         ("end-index", po::value<int>()->default_value(-1),
                          "At which sequence position to end")
                         ("subtract-aa-likelihood", po::value<bool>()->default_value(false),
                          "Subtract sequence log-likelihood")
                         ;

                    // Register option description types with class
                    target.add(options_general);

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

     //! Constructor
     //!
     //! \param argc Number of command line arguments
     //! \param argv Array of arguments
     Options(int argc, char* argv[])
          : phaistos::Options(argc, argv) {

          bool force_reparsing = true;
          phaistos::Options::parse_initial<DefineOptionsPrimary, 
                                           DefineOptionsSecondary>(force_reparsing);

          // Output help message if requested
          if (empty_command_line || this->has_option("help")) {
               std::cout << this->generate_help_output() << "\n";
               exit(1);
          }

     }

};


//! Templated main function.
template<typename DBN_TYPE>
void likelihood_main(::Options &options) {

     using namespace definitions;

     // Print options to stdout
     if (options["verbose"].as<bool>())
          std::cout << options << "\n";

     DBN_TYPE *dbn = NULL;
     initialize_dbn(options, &dbn);

     if (dbn->sequence_length == 0) {
          std::cerr << "Zero length DBN. Aborting\n";
          exit(1);
     }

     if (options["verbose"].as<bool>())
          std::cout << "DBN:\n" << *dbn << "\n";

     // Calculate likelihood
     double ll = dbn->get_log_likelihood(options["start-index"].as<int>(), options["end-index"].as<int>());

     // Optionally subtract sequence likelihood
     double seq_ll = 0;
     if(options["subtract-aa-likelihood"].as<bool>()) {
          // Set emission status of all nodes to false
          (*dbn).template set_emission_state<typename DBN_TYPE::ALL_CHILD_NODES>(false);
          
          // Set emission status of aa node to true
          (*dbn).template set_emission_state<typename DBN_TYPE::AA_NODE>(true);                    
          
          // Calculate likelihood
          seq_ll = dbn->get_log_likelihood(options["start-index"].as<int>(), options["end-index"].as<int>());
     }          

     std::cout << ll - seq_ll << "\n";
}


//! Main function
int main(int argc, char *argv[]) {

     // Parse command line options
     ::Options options(argc, argv);

     if (options["verbose"].as<bool>()) {
          printf("\n");
          printf("############################################################\n");
          printf("#                                                          #\n");
          printf("#          Calculate log-likelihood of input data          #\n");
          printf("#                                                          #\n");
          printf("#  This program uses TorusDBN. See                         #\n");
          printf("#     Boomsma, Mardia, Taylor, Ferkinghoff-Borg,           #\n");
          printf("#     Krogh and Hamelryck, PNAS, 2008                      #\n");
          printf("#   http://sourceforge.net/projects/phaistos               #\n");
          printf("#   http://www.phaistos.org                                #\n");
          printf("#                                                          #\n");
          printf("############################################################\n\n");
     }

     // Call main function either in Torus or FB5 mode
     if (options["backbone-dbn-fb5"].occurrences()) {
          likelihood_main<FB5DBN>(options);
     } else if (options["backbone-dbn-torus"].occurrences()) {
          likelihood_main<TorusDBN>(options);
     } else if (options["backbone-dbn-torus-omega"].occurrences()) {
          likelihood_main<TorusOmegaDBN>(options);
     } else if (options["backbone-dbn-torus-cs"].occurrences()) {
          likelihood_main<TorusCsDBN>(options);
     } else if (options["backbone-dbn-torus-cs-omega"].occurrences()) {
          likelihood_main<TorusCsOmegaDBN>(options);
     } else {
          std::cerr << "No backbone DBN model specified. Aborting\n";
          exit(1);
     }
}

