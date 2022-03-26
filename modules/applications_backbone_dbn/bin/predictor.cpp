// predictor.cpp --- Predict optimal amino acid or secondary structure labeling
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
                         ("output-mode", po::value<std::string>()->default_value("ss"),
                          "Whether to predict amino acid sequence (aa) or secondary structure sequence (ss)")
                         ("mode", po::value<std::string>()->default_value("viterbi"),
                          "Whether to use viterbi (viterbi), posterior decoding (posterior) or one-best (one-best)")
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



//! Sort vector and keep track of corresponding indices
//! \param v Input vector
//! \param v_index Output index vector
void sort_vector(std::vector<double> *v, std::vector<int> *v_index) {

     for (unsigned int i=0; i<v->size(); i++) {
	  for (unsigned int j=i; j<v->size(); j++) {
	       if ((*v)[j] > (*v)[i]) {
		    double tmp = (*v)[i];
		    (*v)[i] = (*v)[j];
		    (*v)[j] = tmp;

		    int tmp_index = (*v_index)[i];
		    (*v_index)[i] = (*v_index)[j];
		    (*v_index)[j] = tmp_index;
	       }
	  }
     }
}

//! Create output for viterbi decoding
//! \param options Program options
//! \param viterbi_path Viterbi sequence of states
//! \param viterbi_emission Distribution over output observable corresponding to viterbi sequence
//! \param labels Text labels used for output
//! \param header Header used for output
void print_viterbi_output(::Options &options, 
                          std::vector<int> viterbi_path, 
                          std::vector<std::vector<double> > viterbi_emission, 
                          const char **labels, const char *header=NULL) {

     std::string label_string = "";
     if (options["verbose"].as<bool>() && header)
	  printf("%s\n", header);
     for (unsigned int i=0; i<viterbi_path.size(); i++) {
	  
	  std::vector<int> index_vector;
	  for (unsigned int j=0; j<viterbi_emission[i].size(); j++){
	       index_vector.push_back(j);
	  }
	  
	  sort_vector(&(viterbi_emission[i]), &index_vector);

          label_string += labels[index_vector[0]];

          if (options["verbose"].as<bool>()) {
               printf("%2d\t", viterbi_path[i]+1);
               for (unsigned int j=0; j<viterbi_emission[i].size(); j++){
                    printf("%s(%.2f)   ", labels[index_vector[j]], viterbi_emission[i][j]);			 
               }
               printf("\n");
          }
     }
     std::cout << label_string << "\n";
}

//! Create output for posterior decoding
//! \param options Program options
//! \param posterior_hidden Posterior distribution over hidden nodes for each position in the sequence
//! \param posterior_emission Distribution over output variables corresponding to posterior distribution
//! \param labels Text labels used for output
//! \param header Header used for output
void print_posterior_output(::Options &options, 
                            std::vector<std::vector<double> > posterior_hidden, 
                            std::vector<std::vector<double> > posterior_emission, 
                            const char **labels, const char *header=NULL) {

     std::vector<int> max_h;
     if (options["verbose"].as<bool>()) {
	  printf("h-values (ranked)\n");
	  for (unsigned int i=0; i<posterior_hidden.size(); i++) {
	       
	       std::vector<int> index_vector;
	       for (unsigned int j=0; j<posterior_hidden[i].size(); j++){
		    index_vector.push_back(j);
	       }	  
	       sort_vector(&(posterior_hidden[i]), &index_vector);

	       for (unsigned int j=0; j<posterior_hidden[i].size(); j++){
		    printf("%2d(%.2f)   ", index_vector[j]+1, posterior_hidden[i][j]);			 
	       }
	       printf("\n");
	       
	       max_h.push_back(index_vector[0]);
	  }
	  printf("\n");
     }
     
     if (header && options["verbose"].as<bool>())
	  printf("%s\n", header);

     std::string label_string = "";
     
     for (unsigned int i=0; i<posterior_emission.size(); i++) {
	  std::vector<int> index_vector;
	  for (unsigned int j=0; j<posterior_emission[i].size(); j++){
	       index_vector.push_back(j);
	  }
	  
	  sort_vector(&(posterior_emission[i]), &index_vector);

          label_string += labels[index_vector[0]];
	       
          if (options["verbose"].as<bool>()) {
	       printf("%2d\t", max_h[i]+1);
	       for (unsigned int j=0; j<posterior_emission[i].size(); j++){
		    printf("%s(%.2f)   ", labels[index_vector[j]], posterior_emission[i][j]);			 
	       }
	       printf("\n");
	  }
     }
     std::cout << label_string << "\n";
}

//! Create output for one-best decoding
//! \param options Program options
//! \param sequence Predicted sequence
//! \param labels Text labels used for output
void print_one_best_output(::Options &options, 
                           std::vector<int> sequence, 
                           const char **labels) {

     std::string label_string = "";
     for (unsigned int i=0; i<sequence.size(); i++) {
          label_string += labels[sequence[i]];
     }
     std::cout << label_string << "\n";
}


//! Templated main function
template <typename DBN_TYPE>
void predictor_main(::Options &options) {

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

     if (options["output-mode"].as<std::string>() == "aa") {
          // Set emission status of aa node to true (fixed=false)
          dbn->template set_emission_state<typename DBN_TYPE::AA_NODE>(false);
     } else if (options["output-mode"].as<std::string>() == "ss") {
          // Set emission status of aa node to true (fixed=false)
          dbn->template set_emission_state<typename DBN_TYPE::SS_NODE>(false);
     }

     if (options["mode"].as<std::string>() == "viterbi") {
	  std::vector<int> viterbi_path = dbn->viterbi();
	  if (options["output-mode"].as<std::string>() == "ss") {
	       const char **labels = ss_name;
	       std::vector<std::vector<double> > viterbi_ss = 
                    dbn->template get_node<typename DBN_TYPE::SS_NODE>()->get_viterbi_distribution();
	       print_viterbi_output(options, viterbi_path, viterbi_ss, labels, "h-value\tSS labels (ranked)");
	  } else {
	       const char **labels = (const char **)residue_name_short;
	       std::vector<std::vector<double> > viterbi_aa = 
                    dbn->template get_node<typename DBN_TYPE::AA_NODE>()->get_viterbi_distribution();
	       print_viterbi_output(options, viterbi_path, viterbi_aa, labels, "h-value\tAmino acid labels (ranked)");
	  }
	  
     } else if (options["mode"].as<std::string>() == "posterior") {

	  std::vector<std::vector<double> > posterior_hidden = dbn->posterior();
	  if (options["output-mode"].as<std::string>() == "ss") {
	       const char **labels = ss_name;
	       std::vector<std::vector<double> > posterior_ss = 
                    dbn->template get_node<typename DBN_TYPE::SS_NODE>()->get_posterior_distribution();
	       print_posterior_output(options, posterior_hidden, posterior_ss, labels, "h-value\tSS labels (ranked)");
	  } else {
	       const char **labels = (const char **)residue_name_short;
	       std::vector<std::vector<double> > posterior_aa = 
                    dbn->template get_node<typename DBN_TYPE::AA_NODE>()->get_posterior_distribution();
	       print_posterior_output(options, posterior_hidden, posterior_aa, labels, "h-value\tAmino acid labels (ranked)");
	  }

     } else if (options["mode"].as<std::string>() == "one-best") {

          const char **labels;

	  if (options["output-mode"].as<std::string>() == "ss") {
               labels = ss_name;
               std::vector<int> ss_seq = 
                    dbn->template get_node<typename DBN_TYPE::SS_NODE>()->one_best_decoding_log();
	       print_one_best_output(options, ss_seq, labels);
          } else {
               labels = residue_name_short;
               std::vector<int> aa_seq = 
                    dbn->template get_node<typename DBN_TYPE::AA_NODE>()->one_best_decoding_log();
	       print_one_best_output(options, aa_seq, labels);
          }
     } else {
          std::cerr << "Error: Unknown mode: " << options["mode"].as<std::string>() << "\n";
     }
}


//! Main function
int main(int argc, char *argv[]) {

     // Parse command line options
     ::Options options(argc, argv);

     if (options["verbose"].as<bool>()) {
          printf("\n");
          printf("############################################################\n");
          printf("#                                                          #\n");
          printf("#   Predict amino acid or secondary structure sequence     #\n");
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
          predictor_main<FB5DBN>(options);
     } else if (options["backbone-dbn-torus"].occurrences()) {
          predictor_main<TorusDBN>(options);
     } else if (options["backbone-dbn-torus-omega"].occurrences()) {
          predictor_main<TorusOmegaDBN>(options);
     } else if (options["backbone-dbn-torus-cs"].occurrences()) {
          predictor_main<TorusCsDBN>(options);
     } else if (options["backbone-dbn-torus-cs-omega"].occurrences()) {
          predictor_main<TorusCsOmegaDBN>(options);
     } else {
          std::cerr << "No backbone DBN model specified. Aborting\n";
          exit(1);
     }

}
