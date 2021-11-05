// program_option_parser.cpp --- Utility classes for parsing command lines and config files
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


// Uncomment this line for debugging of boost::spirit parsers
// #define BOOST_SPIRIT_DEBUG

#include "program_option_parser.h"
#include <boost/regex.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "program_options_short_format_parser.h"

namespace phaistos {

// Constructor
ProgramOptionParser::ProgramOptionParser(CommandLine &command_line)
     : command_line(command_line) {

     add_super_group("");
     add_super_group("shorthands");
}


// Check whether all shorthand options of a specific OptionsDescription object are uninitialized
bool ProgramOptionParser::check_all_shorthand_uninitialized(const OptionsDescription &options) {

     const std::vector< boost::shared_ptr<boost::program_options::option_description> >& opt_vec = 
          options.options();
     bool all_uninitialized = true;
     for (unsigned int j=0; j< opt_vec.size(); j++) {
          const boost::program_options::option_description& opt = *(opt_vec[j]);
          if (option_map[opt.long_name()].value().type() == typeid(OptionShorthand) &&
              !option_map[opt.long_name()].value().empty()) {
               all_uninitialized = false;
          }
     }
     return all_uninitialized;
}


// Check whether all shorthand options in a supergroup are uninitialized
bool ProgramOptionParser::check_all_shorthand_uninitialized_in_super_group(std::string super_group_name) {

     bool all_uninitialized = true;

     int sg = options_description_super_group_map[super_group_name];
     std::vector<int> &super_group_indices = options_description_super_group_vector[sg].second;

     for (unsigned int sgi = 0; sgi < super_group_indices.size(); ++sgi) {
               
          // get index
          int i = super_group_indices[sgi];
               
          all_uninitialized = all_uninitialized && check_all_shorthand_uninitialized(*options_description_group_vector[i]);
     }
     return all_uninitialized;
}


// Uses default settings for super groups if all options in super group are uninitialized
bool ProgramOptionParser::apply_super_group_defaults() {

     bool command_line_modified = false;
     for (unsigned int sg = 0; sg < options_description_super_group_vector.size(); ++sg) {

          std::string super_group_name = options_description_super_group_vector[sg].first;

          if (check_all_shorthand_uninitialized_in_super_group(super_group_name)) {
               std::string replacement = super_group_default(super_group_name);
               if (replacement != "") {
                    command_line(command_line.raw_string + replacement + " ");
                    command_line_modified = true;
               }
          }
     }
     return command_line_modified;
}


// Main parse function
void ProgramOptionParser::parse(CommandLine &command_line, bool allow_unregistered, 
                                const std::string config_file_name) {
          
     if (allow_unregistered) {
          
          // Run basic variant of parser allowing unregistered
          // options (some of the options haven't been defined yet
          boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(command_line.argc, command_line.argv)
               .style(boost::program_options::command_line_style::default_style ^ boost::program_options::command_line_style::allow_guessing)
               .options(options_description_main).allow_unregistered().run();

          boost::program_options::store(parsed, option_map);
     } else {

          // Run basic variant of parser allowing unregistered. 
          // For any unregistered options, we check whether a prefix
          // of this option is legal. If so, we set the prefix version, 
          // and send along the rest of the option in compressed form.
          // For instance:  --energy-isd-dist-likelihood-well-alpha 1
          //            =>  --energy-isd-dist-likelihood well:[alpha:1]
          boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(command_line.argc, command_line.argv)
               .style(boost::program_options::command_line_style::default_style ^ boost::program_options::command_line_style::allow_guessing)
               .options(options_description_main).allow_unregistered().run();

          boost::property_tree::ptree ptree;

          for (unsigned int i=0; i<parsed.options.size(); ++i) {

               // Iterate over unregistered options
               if (parsed.options[i].unregistered) {

                    std::string option_name = parsed.options[i].original_tokens[0];
                    std::string remainder = "";

                    // Removing leading "--" or "-"
                    while (option_name[0] == '-')
                         option_name = option_name.substr(1);

                    // Check if a prefix of this option is registered
                    while (option_name != "" && option_name.size() > 2) {
                         size_t pos = option_name.find_last_of('-');

                         if (pos == std::string::npos) {
                              option_name = "";
                              break;
                         }

                         if (remainder.size()!=0)
                              remainder = option_name.substr(pos+1) + "." + remainder;
                         else
                              remainder = option_name.substr(pos+1);
                         option_name = option_name.substr(0,pos);

                         bool approx = false;
                         const boost::program_options::option_description *option = options_description_main.find_nothrow(option_name, approx);

                         // Exit loop if prefix option is defined 
                         if (option != NULL) {

                              // Check that the option is not a super group
                              if (options_description_super_group_map.count(option_name)==0) {
                                   break;
                              } else {
                                   option_name = "";
                              }

                         }
                    }

                    if (option_name != "") {
                         std::string option_name_original = parsed.options[i].original_tokens[0];
                         std::string option_value = parsed.options[i+1].value[0];

                         // Add to property tree
                         ptree.put(option_name + "." + remainder, option_value);

                         // Remove from command line
                         remove_option_from_command_line(option_name_original, option_value);
                    }

               }
          }

          // boost::property_tree::json_parser::write_json( std::cout, ptree.get<std::string>("well").second, false );
          boost::property_tree::ptree::const_iterator end = ptree.end();
          // Iterate over nodes at top level 
          for (boost::property_tree::ptree::const_iterator it = ptree.begin(); it != end; ++it) {
               std::stringstream sstream;
               option_parser::write_short_option_format( sstream, it->second );
               std::string option_str = "--" + it->first + " " + sstream.str();

               // Add option to command line
               command_line(command_line.raw_string + " " + option_str);

          }

          
          // Run final variant of parser 
          boost::program_options::store(boost::program_options::command_line_parser(command_line.argc, command_line.argv)
                    .style(boost::program_options::command_line_style::default_style ^ boost::program_options::command_line_style::allow_guessing)
                    .options(options_description_main).run(), 
                    option_map);
     }

     // Read config file if available
     if (config_file_name != "") {
          std::ifstream ifs(config_file_name.c_str());
          store(parse_config_file(ifs, options_description_main, allow_unregistered), option_map);
     }

     boost::program_options::notify(option_map);
}


// Parse selected super group options only
void ProgramOptionParser::parse_selected_super_groups(const CommandLine &command_line, const std::vector<std::string> &super_groups, 
                                                      const std::string config_file_name) {

     for (unsigned int s=0; s<super_groups.size(); ++s) {
          std::vector<int> &super_group_indices = 
               options_description_super_group_vector[options_description_super_group_map[super_groups[s]]].second;

          for (unsigned int sgi = 0; sgi < super_group_indices.size(); ++sgi) {
          
               // get index
               int i = super_group_indices[sgi];

               // Run basic variant of parser allowing unregistered
               // options (some of the options haven't been defined yet
               boost::program_options::store(boost::program_options::command_line_parser(command_line.argc, command_line.argv)
                         .style(boost::program_options::command_line_style::default_style ^ boost::program_options::command_line_style::allow_guessing)
                         .options(*options_description_group_vector[i]).allow_unregistered().run(), 
                         option_map);

               if (config_file_name != "") {
                    std::ifstream ifs(config_file_name.c_str());
                    store(parse_config_file(ifs, *options_description_group_vector[i], true), option_map);
               }
          }
     }
     // Register options in map
     boost::program_options::notify(option_map);          
}


// Add a super group
void ProgramOptionParser::add_super_group(std::string name) {
     options_description_super_group_vector.push_back(std::make_pair(name,std::vector<int>()));
     options_description_super_group_map[name] = options_description_super_group_vector.size()-1;
     super_group_defaults_map[name] = "";
}


// Return default setting for super group
std::string &ProgramOptionParser::super_group_default(std::string name) {
     return super_group_defaults_map[name];
}


// Add option
void ProgramOptionParser::add(boost::shared_ptr<OptionsDescription> options_description, 
                              std::string super_group, bool register_final) {

     if (!options_description)
          return;

     bool previously_defined = options_description_group_map.count(options_description->title);
     if (!previously_defined) {

          // Add to options_description vector
          options_description_group_vector.push_back(options_description);

          // Add to options_description map
          options_description_group_map.insert(make_pair(options_description->title, 
                                                         options_description_group_vector.size()-1));
          // Add group to specified supergroup
          if (options_description_super_group_map.find(super_group) != 
              options_description_super_group_map.end()) {
               options_description_super_group_vector[options_description_super_group_map[super_group]]
                    .second.push_back(options_description_group_vector.size()-1);
          }
     }

     if (register_final) {
          // Add to main options_description object
          this->options_description_main.add(*options_description);
     }
}


// Add option (using shared_ptr)
void ProgramOptionParser::add(OptionsDescription &options_description, bool hidden, std::string super_group, bool register_final) {
     boost::shared_ptr<OptionsDescription> options_description_ptr = 
          boost::shared_ptr<OptionsDescription>(new OptionsDescription(options_description, 
                                                                       hidden));
     add(options_description_ptr, super_group, register_final);
}


// Overload [] operator
const ProgramOptionParser::OptionValue ProgramOptionParser::operator[](const std::string &key) const {
     if (has_key(key))
          return ProgramOptionParser::OptionValue(key, option_map[key]);
     else 
          return ProgramOptionParser::OptionValue();
}

// Count number of elements matching key
bool ProgramOptionParser::has_key(std::string key) const {
     return option_map.count(key);
}


//! Call-back functor - called when main option is activated
struct ProgramOptionParser::SecondaryShorthandActionOnActivate {

     //! Overload () operator
     //!
     //! \param options_name Options name
     //! \param option_map Internal map of options
     //! \param options_description_main Main options_description object
     //! \param command_line Command line object
     //! \param replacement_string String to replace with
     //! \param multiple_occurrences_allowed Whether multiple occurrences of the same option are allowed
     //! \param val Parsed value
     void operator()(std::string options_name,
                     boost::program_options::variables_map &option_map,
                     boost::program_options::options_description &options_description_main,
                     CommandLine &command_line,
                     std::string replacement_string,
                     bool multiple_occurrences_allowed,
                     const OptionShorthand &val) const {

          const OptionShorthand &option = option_map[options_name].as<OptionShorthand>(); 

          // Remove relevant parts of old command line
          bool removal_succesful = false;
          unsigned int first_pos = -1;
          for (unsigned int i=0; i<option.raw_strings.size(); ++i) {
               std::string old_str = command_line.raw_string;
               std::string query = "--"+options_name+" "+option.raw_strings[i]+" ";
               std::string::size_type pos = old_str.find(query);
               if (pos<first_pos) {
                    first_pos = pos;
               }

               std::string new_str = boost::replace_all_copy(old_str, 
                                                             query,
                                                             "");
               if (old_str != new_str) {
                    command_line(new_str);
                    removal_succesful = true;
               }
          }

          // Complain if multiple values are specified where they are not allowed
          if (!multiple_occurrences_allowed && (
                   val.occurrences() > 1 ||
                   (val.occurrences() == 1 && val.occurrences(0) > 1))) {
               boost::throw_exception(ValidationError("multiple values of --" + options_name + " not allowed"));
          }

          // Parse using recursion level 0 - the individual moves will parse the rest
          unsigned int max_recursion_depth = 0;
          bool label_multiple_occurrences = false;
          bool omit_singletons = true;

          // Append new command line parts
          if (removal_succesful) {
               std::string prefix = "--";

               // If there is not replacement string, we only prefix a single dash
               // - the command_line_format function will prefix another
               if (replacement_string == "")
                    prefix = "-";
                         
               command_line(command_line.raw_string.insert(first_pos, 
                                                           option.command_line_format(prefix+replacement_string, 
                                                                                      label_multiple_occurrences,
                                                                                      omit_singletons, NULL,
                                                                                      max_recursion_depth)));
          }
     }
};


// Creates a shorthand that can be parsed to provide new shorthands
// Examples are --moves and --energies
// The optional replacement_string specifies what the original
// tag should be replaced with as prefix in the output. For instance
// "energies" will have replacement string "energy", so that
// "--energies val1 val2" becomes "--energy-val1 --energy-val2
boost::shared_ptr<ProgramOptionParser::OptionsDescription> ProgramOptionParser::create_secondary_shorthand(std::string group_name,
                                                                                                           std::string option_name,
                                                                                                           std::string replacement_string,
                                                                                                           bool multiple_occurrences_allowed) {
     // Secondary shorthands are hidden by default
     bool hidden = true;

     // Create options description object
     boost::shared_ptr<OptionsDescription> options_description(new ProgramOptionParser::OptionsDescription(group_name, hidden));

     // Create functor
     SecondaryShorthandActionOnActivate action_on_activate;

     // Add move option
     options_description->add_options()
          (option_name.c_str(), 
           boost::program_options::value<OptionShorthand>()
           ->composing()
           ->multitoken()
           ->implicit_value(OptionShorthand("1"))
           ->notifier(boost::bind<void>(action_on_activate,
                                        option_name, 
                                        boost::ref(option_map), 
                                        boost::ref(options_description_main), 
                                        boost::ref(command_line), 
                                        replacement_string,
                                        multiple_occurrences_allowed,
                                        _1)),
           ("shorthand: " + option_name).c_str());

     return options_description;
}


// Call-back functor - called when main option is activated
void ProgramOptionParser::CreateOptionsActionOnActivate::operator()(ProgramOptionParser *parser,
                                                                    std::string title,
                                                                    std::string options_name,
                                                                    boost::program_options::variables_map &option_map,
                                                                    CommandLine &command_line,
                                                                    const OptionShorthand &val) {
     const OptionShorthand &option = option_map[options_name].as<OptionShorthand>(); 

     // When activated, set option group to visible
     parser->options_description_group_vector[parser->options_description_group_map[title]]->hidden = false;

     // Register number of occurrences
     parser->option_occurrences[options_name] = val.occurrences();

     if (!(val.size()==1 && option.is_singleton(0))) {

          // Remove relevant parts of old command line
          std::vector<bool> succesful_removals(option.raw_strings.size(), false);
          unsigned int first_pos = -1;
          int occurrences = 0;
          for (unsigned int i=0; i<option.raw_strings.size(); ++i) {
               std::string old_str = command_line.raw_string;
               std::string query = "--"+options_name;

               int pos = parser->remove_option_from_command_line(options_name, option.raw_strings[i]);

               if (pos >= 0) {
                    succesful_removals[i] = true;
                    occurrences += option.is_singleton(i)?option.get_singleton(i):1;
                    if ((unsigned int)pos < first_pos) {
                         first_pos = pos;
                    }
               }
          }
               
          // Append new command line parts
          if (occurrences>0) {
               bool label_multiple_occurrences = true;
               bool omit_singletons = true;
               std::string new_command_line_str = option.command_line_format("--"+options_name, label_multiple_occurrences, omit_singletons, &succesful_removals);
               command_line(command_line.raw_string.insert(first_pos,"--"+options_name + " " + boost::lexical_cast<std::string>(occurrences) + " " + new_command_line_str));
          }
     }

}
     


//! Call-back functor - called when main option is activated
struct ProgramOptionParser::CreateOptionsReplacementActionOnActivate {

     //! Overload () operator
     //!
     //! \param parser Parser object
     //! \param title Options title
     //! \param options_name Options name
     //! \param replacement_string String to replace with
     //! \param options_description_vector Vector of shared pointers to OptionsDescription objects in group
     //! \param options_description_map (name,index) map of OptionDescription objects in group
     //! \param option_map Internal map of options
     //! \param command_line Command line object
     //! \param val Parsed value
     void operator()(ProgramOptionParser *parser,
                     std::string title,
                     std::string options_name,
                     std::string replacement_string,
                     std::vector<boost::shared_ptr<OptionsDescription> > &options_description_vector, 
                     std::map<std::string, int> &options_description_map,
                     boost::program_options::variables_map &option_map,
                     // boost::program_options::options_description &options_description_main,
                     CommandLine &command_line,
                     const OptionShorthand &val) {
                    
          if (val.size()==1 && val.is_singleton(0)) {
               int value = val.get_singleton(0);

               const OptionShorthand &option = option_map[options_name].as<OptionShorthand>(); 

               // Remove option from the command line
               int pos = parser->remove_option_from_command_line(options_name, option.raw_strings[0]);
               bool removal_succesful = false;
               if (pos > 0) {
                    removal_succesful = true;
               }

               // If removal was successful, add replacement string
               if (removal_succesful) {
                    // Append new command line parts
                    for (int i=0; i<value; ++i) {
                         OptionShorthand option_tmp(replacement_string); 
                         command_line(command_line.raw_string.insert(pos,option_tmp.command_line_format()));
                    }
               } 
          }
     }
};


// Create a simple options expansion (e.g. search & replace in command line)
boost::shared_ptr<ProgramOptionParser::OptionsDescription> ProgramOptionParser::create_option_expansion(std::string group_name,
                                                                                                        std::string options_name,
                                                                                                        std::string replacement) {


     bool previously_defined = options_description_group_map.count(group_name);

     // Create functor
     CreateOptionsReplacementActionOnActivate action_on_activate;

     boost::shared_ptr<OptionsDescription> options_description;
     if (!previously_defined) {

          bool hidden = true;
          options_description = boost::shared_ptr<OptionsDescription>(new OptionsDescription(group_name, hidden));

          // Add main option
          options_description->add_options()
               (options_name.c_str(), 
                boost::program_options::value<OptionShorthand>()
                ->composing()
                ->multitoken()
                ->implicit_value(OptionShorthand("1"))
                ->notifier(boost::bind<void>(action_on_activate,
                                             this,
                                             options_description->title,
                                             options_name,
                                             replacement,
                                             boost::ref(options_description_group_vector),
                                             boost::ref(options_description_group_map), 
                                             boost::ref(option_map), 
                                             // boost::ref(options_description_main), 
                                             boost::ref(command_line), 
                                             _1)),
                ("Activate " + options_name + " [number of occurrences]").c_str());
     } else {
          options_description = options_description_group_vector[options_description_group_map[group_name]];
     }

     return options_description;
}


// Remove substring from command line
int ProgramOptionParser::remove_option_from_command_line(std::string options_name, std::string value) {

     // Regular expression and replacement to escape query_str and value
     const boost::regex esc("([\\^\\.\\$\\|\\(\\)\\[\\]\\*\\+\\?\\/\\\\])");
     const std::string rep("\\\\\\1");

     std::string &search_string = command_line.raw_string;
     std::string replacement_string = "";
     std::string query_str = options_name;
     if (query_str.substr(0,2) != "--")
          query_str = "--"+query_str;

     query_str = regex_replace(query_str, esc, rep, 
                               boost::match_default | boost::format_sed);               

     if (value == "") {
          query_str += "\\s+(-|$)";
          replacement_string = "\\1";
     } else {

          // also escape value (since square brackets for instance occur there)
          value = regex_replace(value, esc, rep, 
                                boost::match_default | boost::format_sed);               

          query_str += "\\s+" + value + "\\s*";
     }

     boost::match_results<std::basic_string<char>::iterator> what;
     boost::regex re(query_str);
     bool success = regex_search(search_string.begin(), search_string.end(), what, re);
     int position = -1;
     if (success) {
          position = what.position();
          search_string = boost::regex_replace(search_string, re, replacement_string, boost::match_default | boost::format_first_only);
          command_line(search_string);
     }

     return position;
}


// Generate option help output
std::string ProgramOptionParser::generate_help_output() {

     std::string output = "";

     for (unsigned int sg = 0; sg < options_description_super_group_vector.size(); ++sg) {

          std::string super_group_name = options_description_super_group_vector[sg].first;
          std::vector<int> &super_group_indices = options_description_super_group_vector[sg].second;

          if (super_group_indices.size() > 0) {
               if (super_group_name == "") {
                    super_group_name = "main";
               }
               output += "\n";
               output += "#####" + std::string(super_group_name.size(), '#') + "#############\n";
               output += "#### " + boost::to_upper_copy(super_group_name)    + " OPTIONS ####\n";
               output += "#####" + std::string(super_group_name.size(), '#') + "#############\n\n";
          }

          for (unsigned int sgi = 0; sgi < super_group_indices.size(); ++sgi) {

               // get index
               int i = super_group_indices[sgi];

               output += boost::lexical_cast<std::string>(*options_description_group_vector[i]) + "\n";
          }
     }
     return output;
}

}
