// program_option_parser.h --- Utility classes for parsing command lines and config files
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

#ifndef PROGRAM_OPTION_PARSER_H
#define PROGRAM_OPTION_PARSER_H

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/bind.hpp>

#include "utils.h"
#include "settings.h"

namespace phaistos {

//! Removes quotes in string - noop for any other type.
//! Behaviour of boost::program_options regarding removal of quotes 
//! is different in different versions of boost. We therefore explicitly ensure
//! that quotes are removed.
template <typename TYPE>
inline TYPE program_option_parser_unquote_string(TYPE &value) {
     return value;
}

//! Removes quotes in string - const version
//! Behaviour of boost::program_options regarding removal of quotes 
//! is different in different versions of boost. We therefore explicitly ensure
//! that quotes are removed.
template <>
inline const std::string program_option_parser_unquote_string<const std::string>(const std::string &value) {
     return boost::trim_copy_if(value, boost::is_any_of("\"\'"));
}

//! Removes quotes in string
//! Behaviour of boost::program_options regarding removal of quotes 
//! is different in different versions of boost. We therefore explicitly ensure
//! that quotes are removed.
template <>
inline std::string program_option_parser_unquote_string<std::string>(std::string &value) {
     return boost::trim_copy_if(value, boost::is_any_of("\"\'"));
}


//! Variable which calls a function everytime it is updated
template<typename TYPE>
class CallbackVariable {
public:

     //! Default constructor
     CallbackVariable()
          : callback() {}

     //! Constructor 
     //! \param variable Variable to be wrapped
     //! \param callback Callback function
     CallbackVariable(const TYPE &variable,
                      boost::shared_ptr<SettingsCallbackBase> callback)
          : variable(variable),
            callback(callback) {}

     //! Assignment operator (with object of variable type)
     const CallbackVariable &operator=(const TYPE &other) {
          variable = other;
          if (callback)
               callback->execute();
               // ((*settings).*(callback))();
          return *this;          
     }

     //! Assignment operator (with object of same type)
     const CallbackVariable &operator=(const CallbackVariable &other) {
          variable = other.variable;
          if (this->callback == NULL) {
               // this->settings = other.settings;
               this->callback = other.callback;
          }
          if (this->callback)
               callback->execute();
               // ((*settings).*(callback))();
          return *this;          
     }

     //! String input operator
     friend std::istream &operator>>(std::istream &i, CallbackVariable &cv) {
          i >> cv.variable;
          if (cv.callback)
               cv.callback->execute();
               // (*(cv.settings).*(cv.callback))();
          return i;
     }

     //! String output operator
     friend std::ostream& operator<<(std::ostream& o, const CallbackVariable &cv) {
          o << cv.variable;
          return o;
     }

     //! Extract value
     const TYPE& value() const {
          return this->variable;
     } 
private:
     //! Variable
     TYPE variable;

     //! Callback
     boost::shared_ptr<SettingsCallbackBase> callback;
};



//! Program option parser class
class ProgramOptionParser {

public:

     //! Local Exception class
     class ValidationError: public std::logic_error {
     public:
          //! Constructor
          ValidationError(const std::string& what) 
               : logic_error(what) {}

          //! Destructor
          ~ValidationError() throw() {}
     };

     //! A simple wrapper for variable values - includes their name
     class OptionValue {

          //! Extract value as a specific type (possible involving conversion) - general case
          template <typename TYPE>
          TYPE get_value(TYPE *dummy) const {
               if (value->value().type() == typeid(ProgramOptionParser::WrappedEnum<TYPE>))
                    return (TYPE&)value->as<ProgramOptionParser::WrappedEnum<TYPE> >();
               // Behaviour of boost::program_options regarding removal of quotes 
               // is different in different versions of boost. We therefore explicitly ensure
               // that quotes are removed.
               else if (value->value().type() == typeid(std::string)) {
                    return program_option_parser_unquote_string(value->as<TYPE>());
               } else
                    return value->as<TYPE>();
          }

          //! Extract value (possible involving conversion) - string case
          std::string get_value(std::string *dummy) const {
               if (value->value().type() == typeid(std::string)) {
                    return program_option_parser_unquote_string(value->as<std::string>());
               } else {
                    return program_option_parser_unquote_string((std::string&)value->as<ProgramOptionParser::WrappedStringPointer>());                    
               }
          }

     public:
          //! option name
          std::string name;

          //! option value
          const boost::program_options::variable_value *value;

          //! Constructor (default)
          OptionValue()
               : name(""), value(NULL) {}

          //! Constructor
          //!
          //! \param name Option name
          //! \param value Option value
          OptionValue(const std::string &name, const boost::program_options::variable_value &value)
               : name(name), value(&value) {}

          //! Extract option value as specific type.
          //! Special case: if its a wrapped Enum type - cast directly to the Enum type
          //!
          //! \tparam TYPE Target type (value should be extracted as this type)
          //! \return Copy of value
          template <typename TYPE>
          TYPE as() const {
               return get_value((TYPE*)NULL);
          }

          //! Extract option value as specific type - const case
          //! Special case: if its a wrapped Enum type - cast directly to the Enum type
          //!
          //! \tparam TYPE Target type (value should be extracted as this type)
          // template <typename TYPE>
          // const TYPE &as() const {
          //      if (value->value().type() == typeid(ProgramOptionParser::WrappedEnum<TYPE>))
          //           return (TYPE&)value->as<ProgramOptionParser::WrappedEnum<TYPE> >();
          //      else
          //           return value->as<TYPE>();
          // }


          //! Return number of occurrences
          int occurrences() const {
               if (!value || value->value().empty()) {
                    return 0;
               } else if (value->value().type() == typeid(OptionShorthand)) {
                    return as<OptionShorthand>().occurrences();                    
               } else {
                    return 1;
               }
          }

          //! Return type of value
          const std::type_info& type() const {
               return value->value().type();
          }
     };


     //! Command line option type: a plain vector 
     //! This is simply a wrapper for a std::vector<TYPE>, but overrides
     //! the default parsing style
     //!
     //! \tparam TYPE Type of vector elements
     template <typename TYPE>
     class PlainVector: public std::vector<TYPE>  {
     public:

          //! Constructor (default)
          PlainVector() {}

          //! Constructor
          //!
          //! \param v Input vector
          PlainVector(const std::vector<TYPE> &v)
               : std::vector<TYPE>(v) {}

          //! Override validate function. This function is necessary
          //! for boost::program_options to use this type
          friend void validate(boost::any& v, 
                               const std::vector<std::string>& values,
                               PlainVector* target_type, int){

               try {
                    // Make sure no previous assignment was made.
                    boost::program_options::validators::check_first_occurrence(v);

               } catch(const boost::program_options::multiple_occurrences &exception) {
                    boost::throw_exception(
                         std::logic_error(
                         // boost::program_options::multiple_occurrences(
                              "Multiple occurrences of option with PlainVector type - value: " + 
                              boost::lexical_cast<std::string>(values)));                    
               }

               // Use lexical cast to transform string to value
               v = boost::any(PlainVector<TYPE>(boost::lexical_cast<std::vector<TYPE> >(values[0])));
          }          
     };


     //! Command line option type: WrappedEnum - makes input of enums possible
     //! Note: the enum itself must define an operator<< and an operator>> (in namespace std)
     //!
     //! \tparam Wrapped type
     template <typename TYPE>
     class WrappedEnum {

          //! Enum type
          TYPE val;

     public:

          //! Constructor
          //!
          //! \param val Value to be wrapped
          WrappedEnum(const TYPE &val)
               : val(val) {}

          //! Constructor - from string
          //!
          //! \param str Input string
          WrappedEnum(std::string str):
               val(static_cast<TYPE>(std::numeric_limits<int>::min())) {
               std::istringstream stream(str);
               stream >> val;

               // Check if ee.value has been changed
               if (val==std::numeric_limits<int>::min()) {
                    throw boost::program_options::invalid_option_value(str);
               }
          }

          //! Override validate function. This function is necessary
          //! for program_options to use this type
          friend void validate(boost::any& v, 
                               const std::vector<std::string>& values,
                               WrappedEnum* target_type, int){

               // Make sure no previous assignment was made.
               boost::program_options::validators::check_first_occurrence(v);

               // Call the | operator to combine elements
               for (unsigned int i=0; i<values.size(); ++i) {
                    if (v.empty())
                         v = boost::any(WrappedEnum(values[i]));
                    else {
                         v = boost::any(
                                  WrappedEnum(
                                      static_cast<TYPE>(
                                         static_cast<TYPE>(boost::any_cast<WrappedEnum<TYPE> >(v)) |
                                         static_cast<TYPE>(WrappedEnum(values[i]))
                                      )
                                  )
                             );
                    }
               }
          }

          //! Type-cast operator to the type of the wrapped enum
          operator TYPE() const {
               return val;
          }

          //! Overload output operator     
          friend std::ostream & operator<<(std::ostream &o, const WrappedEnum &we) {
               o << we.val;
               return o;
          }          
     };


     //! Wrap a double pointer.
     //! This type makes it possible to override the lexical_cast behaviour
     //! of a double, making it possible to read "inf", and "uninitialized" values
     //! The trick is to reinterpret_cast the double pointer to this type. 
     //! lexical cast will create an empty copy of this type (hence the default constructor)
     //! and will fill this value using the >> operator. Finally, this value is assigned
     //! to the original pointer destination using operator=. We therefore reinterpret_cast
     //! back to the original type before assigning the value
     class WrappedDoublePointer {
     public:

          //! Wrapped value
          double value;

          //! Constructor (default)
          WrappedDoublePointer(){}

          //! Overload assignment operator
          //! Here, it is assumed that the target object is in reality
          //! a double pointer, so we reinterpret and then assign.
          //!
          //! \param other Source object from which assignment is made.
          //! \return Current object (this)
          WrappedDoublePointer &operator=(const WrappedDoublePointer &other) {
               *(reinterpret_cast<double*>(this)) = other.value;
               return (*this);
          }

          //! Overload input operator
          friend std::istream &operator>>(std::istream &input, WrappedDoublePointer &ed) {
               std::string raw_string;
               std::getline(input, raw_string);
               boost::algorithm::to_lower(raw_string);
               if (raw_string.find("uninitialized") != std::string::npos) {
                    ed.value = UNINITIALIZED;
               } else if (raw_string.find("-inf") != std::string::npos) {
                    ed.value = -std::numeric_limits<double>::infinity();
               } else if (raw_string.find("inf") != std::string::npos) {
                    ed.value = std::numeric_limits<double>::infinity();
               } else {
                    ed.value = boost::lexical_cast<double>(raw_string);
               }
               return input;
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &output, const WrappedDoublePointer &ed) {
               output << ed.value;
               return output;
          }
     };


     //! Wrap a string pointer.
     //! This type was introduced to avoid the different behaviour of quoting
     //! in the program_options class of different versions of boost
     //! The trick is to reinterpret_cast the double pointer to this type. 
     //! lexical cast will create an empty copy of this type (hence the default constructor)
     //! and will fill this value using the >> operator. Finally, this value is assigned
     //! to the original pointer destination using operator=. We therefore reinterpret_cast
     //! back to the original type before assigning the value
     class WrappedStringPointer {
     public:

          //! Wrapped value
          std::string value;

          //! Constructor (default)
          WrappedStringPointer(){}

          //! Constructor
          WrappedStringPointer(std::string value)
               : value(value){}

          //! Overload assignment operator
          //! Here, it is assumed that the target object is in reality
          //! a double pointer, so we reinterpret and then assign.
          //!
          //! \param other Source object from which assignment is made.
          //! \return Current object (this)
          WrappedStringPointer &operator=(const WrappedStringPointer &other) {
               *(reinterpret_cast<std::string*>(this)) = other.value;
               return (*this);
          }

          //! Overload input operator
          friend std::istream &operator>>(std::istream &input, WrappedStringPointer &ed) {
               std::string raw_string;
               std::getline(input, raw_string);
               ed.value = program_option_parser_unquote_string(raw_string);
               return input;
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &output, const WrappedStringPointer &ed) {
               output << program_option_parser_unquote_string(ed.value);
               return output;
          }
     };


     //! Command line option type: WrappedEnumPointer - makes input of pointers to enums possible
     //! Note: the enum itself must define an operator<< and an operator>> (in namespace std)
     //!
     //! \tparam Wrapped type
     template <typename TYPE>
     class WrappedEnumPointer {
     public:

          //! Wrapped type
          TYPE value;

          //! Constructor
          WrappedEnumPointer() : value(TYPE(std::numeric_limits<int>::min())) {}

          //! Overload assignment operator
          //! Here, it is assumed that the target object is in reality
          //! a TYPE pointer, so we reinterpret and then assign.
          //!
          //! \param other Source object from which assignment is made.
          //! \return Current object (*this)
          WrappedEnumPointer &operator=(const WrappedEnumPointer &other) {
               *(reinterpret_cast<TYPE*>(this)) = other.value;
               return (*this);
          }

          //! Overload input operator
          friend std::istream &operator>>(std::istream &input, WrappedEnumPointer &ee) {
               // Get a copy of the line before passing it on.
               // Note that seekg and tellg cannot be used on a stream from lexical cast.
               std::string line;
               getline(input, line);

               // Pass the line as a stringstream
               std::istringstream strstream(line);
               strstream >> ee.value;

               // Check if ee.value has been changed
               if (ee.value==std::numeric_limits<int>::min()) {
                    throw boost::program_options::invalid_option_value(line);
               }
               return input;
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &output, const WrappedEnumPointer &ee) {
               output << ee.value;
               return output;
          }
     };


     //! Helper class allowing easy handling of C-style command lines
     class CommandLine {
     public:

          //! Raw command line options
          std::string raw_string;

          //! argv array (input from main())
          char **argv;

          //! argc value (input from main())
          int argc;

          //! Constructor (default)
          CommandLine():argv(NULL), argc(0){}

          //! Set command line
          //!
          //! \param command_line_str Input command line string
          void operator()(std::string command_line_str) {

               raw_string = command_line_str;

               // Delete old contents
               for (int i=0; i<argc; ++i) {
                    delete[] argv[i];
               }
               delete[] argv;

               std::vector<std::string> command_line_vector = split(command_line_str, " ");
               argv = new char*[command_line_vector.size()];
               for (unsigned int i=0; i<command_line_vector.size(); ++i) {
                    argv[i] = new char[command_line_vector[i].size()+1];
                    strcpy(argv[i], command_line_vector[i].c_str());
               }
               argc = command_line_vector.size();
          }

          //! Destructor
          ~CommandLine() {
               for (int i=0; i<argc; ++i) {
                    delete[] argv[i];
               }
               delete[] argv;
          }
     };


     //! Helper class making it possible to specify how many occurrences of each option
     //! are currently active
     class Filter {

          //! Main data container
          std::map<std::string, int> occurrences;

          //! Default value used when key is not found
          int default_value;

     public:

          //! Constructor
          //!
          //! \param default_value Value for options not explicitly present in the occurrences map
          Filter(int default_value = 0)
               : default_value(default_value) {}

          //! Constructor
          //!
          //! \param input_map Map of occurrences
          //! \param default_value Value for options not explicitly present in the occurrences map
          Filter(std::map<std::string, int> &input_map, 
                 int default_value = 0)
               : occurrences(input_map),
                 default_value(default_value) {}

          //! Overload [] operator for lookups
          //!
          //! \param key Name of option
          int operator[](std::string key) const {
               std::map<std::string,int>::const_iterator it = occurrences.find(key);
               if (it != occurrences.end()) {
                    return it->second;
               } else {
                    return default_value;
               }
          }
     };


     //! Wrapper class for boost::program_options::options_description,
     //! because this class doesn't allow us to access the title (caption)
     //! of a options_description object
     class OptionsDescription: public boost::program_options::options_description {
     public:
          //! Name of options description
          std::string title;

          //! Optional extra description
          std::string description;

          //! Whether this options_description should be shown in output
          bool hidden;

          //! Constructor
          //!
          //! \param title Title of options block
          //! \param description Description of options block
          //! \param hidden Whether these options are displayed in output
          OptionsDescription(std::string title, std::string description, bool hidden=false)
               : boost::program_options::options_description(title, 150),
                 title(title),
                 description(description),
                 hidden(hidden) {}

          //! Constructor
          //!
          //! \param title Title of options block
          //! \param hidden Whether these options are displayed in output
          OptionsDescription(std::string title, bool hidden=false)
               : boost::program_options::options_description(title, 150),
                 title(title),
                 description(""),
                 hidden(hidden) {}

          //! Copy constructor.
          //!
          //! \param other Source object from which copy is made.
          //! \param hidden Whether these options are displayed in output
          OptionsDescription(const OptionsDescription &other, bool hidden)
               : boost::program_options::options_description(other),
                 title(other.title),
                 description(other.description),
                 hidden(hidden) {}

          //! Overload output operator
          friend std::ostream & operator<<(std::ostream &o, OptionsDescription &od) {
               o << static_cast<boost::program_options::options_description&>(od);
               return o;
          }
     };


     //! Option type: Shorthand.
     //! Makes it possible to specify a collection of options as a single option
     //! Supported syntax examples: 
     //!    --main-option sub-option1[prm1:val1,prm2:val2,...] sub-option2[...]
     //!    --main-option sub-option1[prm1:val1,prm2:val2,...],sub-option2[...]
     //!    --main-option sub-option1[prm1[val1],prm2[val2],...] sub-option2[...]
     //!    --main-option sub-option1:value
     //!    --main-option value
     class OptionShorthand {

          //! Recursive definition: A shorthand consists of a (tag,shorthand) pairs
          std::vector<std::vector<std::pair<std::string, OptionShorthand> > > data;

     public:

          //! Raw command line string
          std::vector<std::string> raw_strings;

          //! Constructor
          //!
          //! \param command_line_str Input in command line format
          OptionShorthand(const std::string &command_line_str) {
               std::istringstream stream(command_line_str);
               stream >> (*this);
          }
               

          //! Constructor
          //!
          //! \param raw_strings Vector of raw input strings
          OptionShorthand(const std::vector<std::string> &raw_strings = std::vector<std::string>()) {
               for (unsigned int i=0; i<raw_strings.size(); ++i) {
                    std::istringstream stream(raw_strings[i]);
                    stream >> (*this);
               }
          }

          //! Determines whether Shorthand at index contains only a single value
          //!
          //! \param index Index into data vector
          //! \return True if data[i] is a singleton
          bool is_singleton(int index) const {
               if ((data[index].size()==1) && (data[index][0].second.data.size()==0)) {
                   int value;
                   std::istringstream stream(data[index][0].first);
                   return (stream >> value);
               } else {
                    return false;
               }
          }

          //! Interprets data at index as a singleton value
          //!
          //! \param index Index into data vector
          //! \return Value cast to integer
          int get_singleton(int index) const {
               return boost::lexical_cast<int>(data[index][0].first);
          }


          //! Return number of occurrences. 
          //! normally: the size of data vector.
          //! in case some of the data members are singletons, use their value
          //! \return Number of occurrences
          int occurrences() const {
               int value = 0;
               for (unsigned int i=0; i<data.size(); ++i) {
                    if (is_singleton(i))
                         value += get_singleton(i);
                    else
                         value += 1;
               }
               return value;
          }

          //! Return number of occurrences at specified index in data vector
          //!
          //! \param index Index into data vector
          int occurrences(int index) const {
               if (is_singleton(index))
                    return get_singleton(index);
               else
                    return data[index].size();
          }


          //! Return size of data vector
          int size() const {
               return data.size();
          }

          //! Output in command line format
          //! max_recursion_depth specifies how far in the recursion the command line output should go - the rest
          //! is outputted in unprocessed style, and can thus be parsed by a subsequent OptionShorthand parser.
          //! This makes it possible to parse an OptionShorthand in several stages, for instance:
          //! Input:                                        --move dbn[weight:1] cra[weight:2]
          //! Stage 1(using max_recursion 1):               --move-dbn weight[1] --move-cra weight[2]
          //! Stage 2(applied on each of the options):      --move-dbn-weight 1 --move-cra-weight 2
          //! This is the same output as we would have obtained by parsing it all in one go,
          //! but parsing in two stages gives slightly more flexibility, since the prefix and the 
          //! methods of replacement can be changed. For instance, we would like the output of phase 2
          //! to be:                     --move-dbn --move-dbn-weight 1 --move-cra --move-cra-weight 2
          //! so that the moves are actually turned on as well. 
          //! Note: assigning -1 to an unsigned in results in the maximum possible value
          //!
          //! \param prefix String to add when contructing command line
          //! \param label_multiple_occurrences Whether to add numeric tag if an option occurs multiple times
          //! \param omit_singletons Whether singletons whould be skipped
          //! \param active_indices_pointer Makes it possible to mark certain options as inactive, in which case they will be skipped
          //! \param max_recursion_depth Specifies how far into the recursion we should go in a single call.
          //! \param current_recursion_depth Current recursion level.
          std::string command_line_format(std::string prefix="-",
                                          bool label_multiple_occurrences = false,
                                          bool omit_singletons = false,
                                          std::vector<bool> *active_indices_pointer = NULL, 
                                          unsigned int max_recursion_depth=(unsigned)-1,
                                          unsigned int current_recursion_depth=0) const {

               // Use local active index vector if non if provided
               std::vector<bool> active_indices(data.size(), true);
               if (!active_indices_pointer)
                    active_indices_pointer = &active_indices;

               // Count number of inactive options - this will be used as an offset
               // since inactive options are assumed to have lowest indices
               // (since they are typically defined in a config file)
               int tag_offset = 0;
               for (unsigned int i=0; i<data.size(); ++i) {
                    if (!(*active_indices_pointer)[i])
                         tag_offset += is_singleton(i)?get_singleton(i):1;
               }

               std::string output = "";

               for (unsigned int i=0; i<data.size(); ++i) {

                    if (!(*active_indices_pointer)[i])
                         continue;

                    // If the content is just a single element, simply output this element
                    if (current_recursion_depth == 0 && is_singleton(i)) {
                         if (omit_singletons && current_recursion_depth == 0) {
                              continue;
                         } else {
                              output += prefix + " " + data[i][0].first + " ";
                              continue;
                         }
                    }

                    for (unsigned int j=0; j<data[i].size(); ++j) {

                         // Main recursion
                         if (current_recursion_depth<max_recursion_depth && data[i][j].second.data.size() > 0) {
                              std::string tag = data[i][j].first;
                              // This assumes that the inactive options are located at the end
                              // of the data vector. Currently, inactive options are only 
                              // created by config-files, which are parsed last, so this 
                              // should be ok
                              if (label_multiple_occurrences && (i>0 || tag_offset > 0))
                                   tag = boost::lexical_cast<std::string>(tag_offset+i)+"-"+tag;
                              output += data[i][j].second.command_line_format(prefix+"-"+tag, 
                                                                              label_multiple_occurrences,
                                                                              omit_singletons,
                                                                              NULL,
                                                                              max_recursion_depth,
                                                                              current_recursion_depth+1);

                         // In case recursion was interupted (due to max_recursion_depth parameter)
                         } else if (data[i][j].second.data.size() > 0) {

                              // At the top level of recursion, singular tags are treated as part of an option
                              // rather than an argument. That is --move dbn becomes --move-dbn. At any other
                              // level, a singular tag denotes a value (e.g. --moves dbn[value])
                              if (current_recursion_depth==0) {
                                   output += prefix + "-" + data[i][j].first + " ";
                              }

                              // The recursion was interrupted. The remaining data is printed out in unparsed form. 
                              // This makes it possible to do multi-state parsing
                              if (data[i][j].second.data.size() > 0) {
                                   output += boost::lexical_cast<std::string>(data[i][j].second)+" ";
                              } 

                         // Recursion end due to end of data
                         } else { 

                              // At the top level of recursion, singular tags are treated as part of an option
                              // rather than an argument. That is --moves dbn becomes --moves-dbn. At any other
                              // level, a singular tag denotes a value (e.g. --moves dbn[value])
                              if (current_recursion_depth==0) {
                                   output += prefix + "-" + data[i][j].first + " 1 ";
                              } else {
                                   output += prefix + " " + data[i][j].first + " ";
                              }
                         }
                    }
               }
               return output;
          }


          //! Overload + operator to combine two shorthand options
          friend OptionShorthand operator+(const OptionShorthand &os1, const OptionShorthand& os2) {

               // Concatenate raw strings
               std::vector<std::string> raw_strings = os1.raw_strings;
               raw_strings.insert(raw_strings.end(), os2.raw_strings.begin(), os2.raw_strings.end());
               
               return OptionShorthand(raw_strings);
          }


          //! Overload validate function. This function is necessary
          //! for program_options to use this type
          friend void validate(boost::any& v, 
                               const std::vector<std::string>& values,
                               OptionShorthand* target_type, int){

               std::string merged = boost::join(values, " ");

               if (!v.empty()) {

                    const OptionShorthand &previous = boost::any_cast<OptionShorthand>(v);
                    const OptionShorthand next(merged) ;
                    v = boost::any(previous+next);

               } else {
                    v = boost::any(OptionShorthand(merged));
               }
          }


          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &o, const OptionShorthand &os) {

               std::string output = "";

               for (unsigned int i=0; i<os.data.size(); ++i) {
                    for (unsigned int j=0; j<os.data[i].size(); ++j) {
                         output += os.data[i][j].first;
                         if (os.data[i][j].second.data.size() > 0)
                              output += "[" + boost::lexical_cast<std::string>(os.data[i][j].second) + "]";
                         if (j<(os.data[i].size()-1)) {
                              output += ",";
                         }
                    }
                    if (i<(os.data.size()-1)) {
                         output += ",";
                    }
               }

               // Newer versions of boost do not accept empty output from lexical case
               if (output.empty())
                    output = " ";

               o << output;

               return o;
          }

          //! Input functionality
          //!
          //! \param i Input stream
          //! \param open_parentheses The number of currently open parentheses
          void input(std::istream &i, int open_parentheses=0) {

               // Skip whitespace
               while(i.good() && i.peek() == ' ')
                    i.ignore(1);

               // Keep track of delimiter type
               bool colon_delimiter = false;
               bool square_bracket_begin = false;
               if (i.good()) {
                    if (i.peek() == '[' ) {
                         i.ignore(1);
                         square_bracket_begin = true;
                    } else if (i.peek() == ':') {
                         i.ignore(1);
                         colon_delimiter = true;
                    }
               }

               data.push_back(std::vector<std::pair<std::string, OptionShorthand> >());

               // Iterate over list of iterms
               while (i.good() && (i.peek() != ']')) {
               
                    // Read a tag
                    std::string tag = "";
                    while (i.good() && (i.peek() != '[' && i.peek() != ']' && i.peek() != ',' && i.peek() != ':' && i.peek() != ' ')) {
                         char tmp;
                         i.read(&tmp,1);
                         if (!i.fail()) {
                              tag.push_back(tmp);
                         }
                    }

                    // If no tag could be read using the normal approach
                    // we use a more flexible approach, allowing tags containing
                    // characters not normally allowed ('[',']',','). This is
                    // to allow parsing vectors as arguments
                    if (tag == "") {
                         // Do not allow ']' and ',' at outer level
                         int open_parentheses_inner = 0;
                         while (i.good() && (i.peek() != ':' && i.peek() != ' ' && 
                                             !(i.peek()==']' && open_parentheses_inner==0) && 
                                             !(i.peek()==',' && open_parentheses_inner==0))) {
                              char tmp;
                              i.read(&tmp,1);
                              if (!i.fail()) {
                                   tag.push_back(tmp);
                                   if (tmp == '[') {
                                        open_parentheses_inner++;
                                   } else if (tmp == ']') {
                                        open_parentheses_inner--;
                                   }
                              }
                         }
                    }

                    // Add element to vector
                    data.back().push_back(make_pair(tag,OptionShorthand()));
                    if (i.good() && (i.peek() == '[' || i.peek() == ':')) {
                         // Call recursively
                         data.back().back().second.input(i, open_parentheses+(square_bracket_begin?1:0));

                    }

                    if (i.good() && (i.peek() == ',' || i.peek() == ' ')) {
                         // Only allow more than one element if we are not using a colon delimiter
                         if (!colon_delimiter)
                              i.ignore(1);
                         else
                              break;
                    }
               }

               // Only allow closing bracket if beginning bracket was used
               if (square_bracket_begin && i.good() && i.peek() == ']') {
                    i.ignore(1);
               }

          }


          //! Input operator
          friend std::istream &operator>>(std::istream &i, OptionShorthand &os) {

               // Extract entire string
               std::string raw_string;
               std::getline(i, raw_string);
               os.raw_strings.push_back(raw_string);

               // Create new stream based on string
               std::istringstream istream(raw_string);

               os.input(istream);

               return i;
          }
     };


     //! Check whether all shorthand options of a specific OptionsDescription object are uninitialized
     //!
     //! \param options OptionsDescription input
     //! \return True if all options are unintitialized
     bool check_all_shorthand_uninitialized(const OptionsDescription &options);

     //! Check whether all shorthand options in a supergroup are uninitialized
     //!
     //! \param super_group_name Name of super group
     //! \return True if all options are unintitialized
     bool check_all_shorthand_uninitialized_in_super_group(std::string super_group_name);

     //! Uses default settings for super groups if all options in super group are uninitialized
     //!
     //! \return True if command line was modified
     bool apply_super_group_defaults();



     //////////////////////////////////////     
     // Boost::program_options variables //
     //////////////////////////////////////

     //! Map between options names and values
     boost::program_options::variables_map option_map;

     //! Container of option descriptions
     boost::program_options::options_description options_description_main;



     //! Vector of local options_description objects.
     //! All options descriptions are maintained both internally
     //! in program_options and in this vector, so we can iterate
     //! over them when creating config_file output.
     std::vector<boost::shared_ptr<OptionsDescription> > options_description_group_vector;

     //! Map allowing lookups in options_description_vector by string
     std::map<std::string, int> options_description_group_map;

     //! Vector of supergroups
     std::vector<std::pair<std::string,std::vector<int> > > options_description_super_group_vector;

     //! Map of supergroups. Maps string to index
     std::map<std::string, int> options_description_super_group_map;

     //! Map of supergroup default values
     std::map<std::string, std::string> super_group_defaults_map;

     //! Settings objects associated with various options (e.g. for move and energy objects)
     std::map<std::string,std::vector<boost::shared_ptr<Settings> > > settings_objects;

     //! Keep track of the number of copies of each option
     std::map<std::string, int> option_occurrences;

     //! Object keeping track of both textual and parsed version of command line
     CommandLine &command_line;

public:

     //! Constructor
     //!
     //! \param command_line Command line object input
     ProgramOptionParser(CommandLine &command_line);

     //! Desctructor
     ~ProgramOptionParser() {}

     //! Parse - with definition functor
     //!
     //! \param func Functor setting up options
     //! \param command_line Command line to parse
     //! \param allow_unregistered Whether to allow unregistered options
     //! \param config_file_name Name for optional config file
     template <typename FUNCTOR_TYPE>
     void parse(const FUNCTOR_TYPE &func, CommandLine &command_line, 
                bool allow_unregistered=false, const std::string config_file_name="") {

          // Execute definition function
          func(*this);
          
          parse(command_line, allow_unregistered, config_file_name);
     }

     //! Parse - with two definition functors
     //!
     //! \param func1 Functor setting up options
     //! \param func2 Functor setting up options
     //! \param command_line Command line to parse
     //! \param allow_unregistered Whether to allow unregistered options
     //! \param config_file_name Name for optional config file
     template <typename FUNCTOR_TYPE1, typename FUNCTOR_TYPE2>
     void parse(const FUNCTOR_TYPE1 &func1, const FUNCTOR_TYPE2 &func2, CommandLine &command_line, 
                bool allow_unregistered=false, const std::string config_file_name="") {

          // Execute definition function
          func1(*this);
          func2(*this);
          
          parse(command_line, allow_unregistered, config_file_name);
     }
          

     //! Main parse function
     //!     
     //! \param command_line Command line to parse
     //! \param allow_unregistered Whether to allow unregistered options
     //! \param config_file_name Name for optional config file
     void parse(CommandLine &command_line, bool allow_unregistered=false, const std::string config_file_name="");


     //! Parse selected super group options only - with definition functor
     //!     
     //! \param func Functor settings up options
     //! \param command_line Command line to parse
     //! \param super_group Super group name
     //! \param config_file_name Name for optional config file
     template <typename FUNCTOR_TYPE>
     void parse_selected_super_groups(const FUNCTOR_TYPE &func, const CommandLine &command_line,
                                      const std::string &super_group,
                                      const std::string config_file_name="") {
          
          // Execute definition function
          func(*this);
          
          parse_selected_super_groups(command_line, std::vector<std::string>(1,super_group), config_file_name);
     }
     
     //! Parse selected super group options only - with definition functor
     //!     
     //! \param func Functor settings up options
     //! \param command_line Command line to parse
     //! \param super_groups Vector of super group names
     //! \param config_file_name Name for optional config file
     template <typename FUNCTOR_TYPE>
     void parse_selected_super_groups(const FUNCTOR_TYPE &func, const CommandLine &command_line,
                                      const std::vector<std::string> &super_groups,
                                      const std::string config_file_name="") {
          
          // Execute definition function
          func(*this);
          
          parse_selected_super_groups(command_line, super_groups, config_file_name);
     }
     
     
     //! Parse selected super group options only - with two definition functors
     //!     
     //! \param func1 Functor settings up options
     //! \param func2 Functor settings up options
     //! \param command_line Command line to parse
     //! \param super_group Super group name
     //! \param config_file_name Name for optional config file
     template <typename FUNCTOR_TYPE1, typename FUNCTOR_TYPE2>
     void parse_selected_super_groups(const FUNCTOR_TYPE1 &func1, const FUNCTOR_TYPE2 &func2, const CommandLine &command_line,
                                      const std::string &super_group,
                                      const std::string config_file_name="") {
          
          // Execute definition function
          func1(*this);
          func2(*this);
          
          parse_selected_super_groups(command_line, std::vector<std::string>(1,super_group), config_file_name);
     }
     

     //! Parse selected super group options only - with two definition functors
     //!     
     //! \param func1 Functor settings up options
     //! \param func2 Functor settings up options
     //! \param command_line Command line to parse
     //! \param super_groups Vector of super group names
     //! \param config_file_name Name for optional config file
     template <typename FUNCTOR_TYPE1, typename FUNCTOR_TYPE2>
     void parse_selected_super_groups(const FUNCTOR_TYPE1 &func1, const FUNCTOR_TYPE2 &func2, const CommandLine &command_line,
                                      const std::vector<std::string> &super_groups,
                                      const std::string config_file_name="") {
          
          // Execute definition function
          func1(*this);
          func2(*this);
          
          parse_selected_super_groups(command_line, super_groups, config_file_name);
     }


     //! Parse selected super group options only - no definition functors
     //!
     //! \param command_line Command line to parse
     //! \param super_groups Vector of super group names
     //! \param config_file_name Name for optional config file
     void parse_selected_super_groups(const CommandLine &command_line, 
                                      const std::vector<std::string> &super_groups, 
                                      const std::string config_file_name="");


     //! Add a super group
     //!
     //! \param name Super group name
     void add_super_group(std::string name);

     //! Return default setting for super group
     //!
     //! \param name Super group name
     //! \return default value for the specified supergroup
     std::string &super_group_default(std::string name);

     //! Add option (using shared_ptr)
     //!
     //! \param options_description shared_ptr to an OptionsDescription object
     //! \param super_group Super group in which to add option
     //! \param register_final If option is present multiple times, this should only be true the last time
     void add(boost::shared_ptr<OptionsDescription> options_description, 
              std::string super_group = "", bool register_final=true);

     //! Add option 
     //!
     //! \param options_description OptionsDescription object
     //! \param hidden Whether option should be displayed in output
     //! \param super_group Super group in which to add option
     //! \param register_final If option is present multiple times, this should only be true the last time
     void add(OptionsDescription &options_description, bool hidden=false, std::string super_group="", 
              bool register_final=true);

     //! Overload [] operator
     //!
     //! \param key Option name
     //! \return Value assigned to option
     const OptionValue operator[](const std::string &key) const;

     //! Whether option key is present in parser
     //!
     //! \param key Option name
     //! \return True if option is present
     bool has_key(std::string key) const;


     //! Call-back functor - called when main option is activated
     struct SecondaryShorthandActionOnActivate;


     //! Creates a shorthand that can be parsed to provide new shorthands
     //! Examples are --moves and --energies
     //! The optional replacement_string specifies what the original
     //! tag should be replaced with as prefix in the output. For instance
     //! "energies" will have replacement string "energy", so that
     //! "--energies val1 val2" becomes "--energy-val1 --energy-val2
     //!
     //! \param group_name Name (displayed in output)
     //! \param option_name Name of option
     //! \param replacement_string Which string the option_name should be replaced with upon expansion
     //! \param multiple_occurrences_allowed Whether multiple occurrences of the same option name are allowed
     boost::shared_ptr<OptionsDescription> create_secondary_shorthand(std::string group_name,
                                                                      std::string option_name,
                                                                      std::string replacement_string = "",
                                                                      bool multiple_occurrences_allowed=true);


     //! Call-back functor - called when main option is activated
     struct CreateOptionsActionOnActivate {

          //! Overload () operator
          //!
          //! \param parser Parser object
          //! \param title Option title
          //! \param options_name Name of option
          //! \param option_map Map between options names and values
          //! \param command_line Command line
          //! \param val Parsed value
          void operator()(ProgramOptionParser *parser,
                          std::string title,
                          std::string options_name,
                          boost::program_options::variables_map &option_map,
                          CommandLine &command_line,
                          const OptionShorthand &val);
     };


     //! Create new option group
     //!
     //! \param default_settings  Functor settings up options
     //! \param group_name Name of group
     //! \param options_name Name of option
     //! \param settings Local settings object
     //! \param suboption_vector Object relating options to settings objects
     template <typename FUNCTOR_TYPE, typename SETTINGS_TYPE, typename SUBOPTION_VECTOR_TYPE>
     boost::shared_ptr<OptionsDescription> create_options(const FUNCTOR_TYPE &default_settings,
                                                          std::string group_name,
                                                          std::string options_name, 
                                                          SETTINGS_TYPE &settings,
                                                          const SUBOPTION_VECTOR_TYPE &suboption_vector) {

          if (!settings)
               return boost::shared_ptr<OptionsDescription>();

          bool previously_defined = settings_objects.count(options_name);

          if (!previously_defined) {
               settings_objects.insert(make_pair(options_name, std::vector<boost::shared_ptr<Settings> >()));
          } 

          settings_objects[options_name].push_back(settings);

          // Create functor
          CreateOptionsActionOnActivate action_on_activate;

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
                     ->implicit_value(OptionShorthand(""))
                     ->notifier(boost::bind<void>(action_on_activate,
                                                  this,
                                                  options_description->title,
                                                  options_name,
                                                  boost::ref(option_map), 
                                                  boost::ref(command_line), 
                                                  _1)),
                     ("Activate " + options_name + " [number of occurrences]").c_str());
          } else {
               options_description = options_description_group_vector[options_description_group_map[group_name]];
          }

          // Create functor
          AddOptions add_options(options_name, options_description);

          // Call default settings functor
          default_settings(add_options, settings);

          // Add specific options provided for this move type
          boost::fusion::for_each(suboption_vector, add_options);

          return options_description;
     }


     //! Create new option group - without Settings object
     //!
     //! \param default_settings  Functor settings up options
     //! \param group_name Name of group
     //! \param options_name Name of option
     //! \param suboption_vector Object relating options to settings objects
     template <typename FUNCTOR_TYPE, typename SUBOPTION_VECTOR_TYPE>
     boost::shared_ptr<OptionsDescription> create_options(const FUNCTOR_TYPE &default_settings,
                                                          std::string group_name,
                                                          std::string options_name,
                                                          const SUBOPTION_VECTOR_TYPE &suboption_vector) {


          bool previously_defined = options_description_group_map.count(group_name);

          // Create functor
          CreateOptionsActionOnActivate action_on_activate;

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
                     ->implicit_value(OptionShorthand(""))
                     ->notifier(boost::bind<void>(action_on_activate,
                                                  this,
                                                  options_description->title,
                                                  options_name,
                                                  boost::ref(option_map), 
                                                  boost::ref(command_line), 
                                                  _1)),
                     ("Activate " + options_name).c_str());
          } else {
               options_description = options_description_group_vector[options_description_group_map[group_name]];
          }


          // Create functor
          AddOptions add_options(options_name, options_description);

          // Call default settings functor
          default_settings(add_options);

          // Add specific options provided for this move type
          boost::fusion::for_each(suboption_vector, add_options);

          return options_description;
     }


     //! Create new option group - without Settings object - without default_settings
     //!
     //! \param group_name Name of group
     //! \param options_name Name of option
     //! \param suboption_vector Object relating options to settings objects
     template <typename SUBOPTION_VECTOR_TYPE>
     boost::shared_ptr<OptionsDescription> create_options(std::string group_name,
                                                          std::string options_name,
                                                          const SUBOPTION_VECTOR_TYPE &suboption_vector) {


          bool previously_defined = options_description_group_map.count(group_name);

          // Create functor
          CreateOptionsActionOnActivate action_on_activate;

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
                     ->implicit_value(OptionShorthand(""))
                     ->notifier(boost::bind<void>(action_on_activate,
                                                  this,
                                                  options_description->title,
                                                  options_name,
                                                  boost::ref(option_map), 
                                                  boost::ref(command_line), 
                                                  _1)),
                     ("Activate " + options_name).c_str());
          } else {
               options_description = options_description_group_vector[options_description_group_map[group_name]];
          }


          // Create functor
          AddOptions add_options(options_name, options_description);

          // Add specific options provided for this move type
          boost::fusion::for_each(suboption_vector, add_options);

          return options_description;
     }


     //! Call-back functor - called when main option is activated
     struct CreateOptionsReplacementActionOnActivate;


     //! Create a simple options expansion (e.g. search & replace in command line)
     //!
     //! \param group_name Name of group
     //! \param options_name Name of option
     //! \param replacement Replacement string
     boost::shared_ptr<OptionsDescription> create_option_expansion(std::string group_name,
                                                                   std::string options_name,
                                                                   std::string replacement);


     //! Remove substring from command line
     //!
     //! \param options_name Name of option
     //! \param value Optional value to be replaced as well
     int remove_option_from_command_line(std::string options_name, std::string value);


     //! Functor used to set a vector of options from fusion vectors to a specified OptionsDescription object
     class AddOptions {

          //! Main functionality - inner version
          //! Function taking a pair consisting of a string and a value of some type
          //! and assigning the corresponding command line option to it. Keeps track
          //! of whether the option has been found so far (called by accumulate)
          //!
          //! \param name Option Name
          //! \param description Option description
          //! \param value_pointer Where option value is stored
          //! \param default_value Optional default value
          //! \param implicit_value_pointer Optional implicit value
          template <typename TYPE>
          void execute_inner(std::string name, std::string description, 
                             TYPE *value_pointer, TYPE default_value=TYPE(), 
                             const TYPE *implicit_value_pointer=NULL) const {
               if (name=="")
                    return;
               
               std::string prefix_expanded = prefix;
               for (int counter=1; 
                    options_description->find_nothrow(prefix_expanded+"-"+name, false); 
                    ++counter) {
                    prefix_expanded = prefix+"-"+boost::lexical_cast<std::string>(counter);
               }

               std::stringstream default_value_string_stream;
               default_value_string_stream << default_value;

               if (implicit_value_pointer != NULL) {
                    if (value_pointer != NULL) {
                         options_description->add_options()
                              ((prefix_expanded+"-"+name).c_str(), 
                               boost::program_options::value(value_pointer)
                               ->default_value(default_value, default_value_string_stream.str())
                               ->implicit_value(*implicit_value_pointer),
                               description.c_str());
                    } else {
                         options_description->add_options()
                              ((prefix_expanded+"-"+name).c_str(), 
                               boost::program_options::value<TYPE>()
                               ->default_value(default_value, default_value_string_stream.str())
                               ->implicit_value(*implicit_value_pointer),
                               description.c_str());
                    }               
               } else {
                    if (value_pointer != NULL) {
                         options_description->add_options()
                              ((prefix_expanded+"-"+name).c_str(), 
                               boost::program_options::value(value_pointer)
                               ->default_value(default_value, default_value_string_stream.str()),
                               description.c_str());
                    } else {

                         options_description->add_options()
                              ((prefix_expanded+"-"+name).c_str(), 
                               boost::program_options::value<TYPE>()
                               ->default_value(default_value, default_value_string_stream.str()),
                               description.c_str());
                    }
               }
          }

          //! Main functionality
          //! Function taking a pair consisting of a string and a value of some type
          //! and assigning the corresponding command line option to it. Keeps track
          //! of whether the option has been found so far (called by accumulate)
          //!
          //! \param name Option Name
          //! \param description Option description
          //! \param value_pointer Where option value is stored
          //! \param default_value Optional default value
          //! \param implicit_value_pointer Optional implicit value
          template <typename TYPE>
          void execute(std::string name, std::string description, 
                       TYPE *value_pointer, TYPE default_value=TYPE(), 
                       const TYPE *implicit_value_pointer=NULL) const {
               if (value_pointer != NULL)
                    default_value = *value_pointer;
               execute_inner(name, description, value_pointer, default_value, implicit_value_pointer);
          }

          //! Main functionality
          //! Function taking a pair consisting of a string and a value of some type
          //! and assigning the corresponding command line option to it. Keeps track
          //! of whether the option has been found so far (called by accumulate)
          //!
          //! \param name Option Name
          //! \param description Option description
          //! \param value_pointer Where option value is stored
          //! \param default_value Optional default value
          //! \param implicit_value_pointer Optional implicit value
          void execute(std::string name, std::string description, 
                       std::string *value_pointer, std::string default_value="", 
                       const std::string *implicit_value_pointer=NULL) const {
               if (value_pointer != NULL && *value_pointer != "")
                    default_value = *value_pointer;
               execute_inner(name, description, 
                             reinterpret_cast<WrappedStringPointer*>(value_pointer), 
                             WrappedStringPointer(default_value),
                             reinterpret_cast<const WrappedStringPointer*>(implicit_value_pointer));
          }

     public:

          //! Option name prefix
          std::string prefix;

          //! Shared pointer to an OptionDescription object - into which options are added
          boost::shared_ptr<OptionsDescription> options_description;

          //! Constructor
          //!
          //! \param prefix Option name prefix
          //! \param options_description OptionsDescription to which options are added
          AddOptions(std::string prefix,
                     boost::shared_ptr<OptionsDescription> options_description)
               : prefix(prefix),
                 options_description(options_description){}

          //! Function taking a pair consisting of a string and a value of some type
          //! and assigning the corresponding command line option to it. Keeps track
          //! of whether the option has been found so far (called by accumulate)
          //! This version includes: (option-name, description, value_pointer, default_value, implicit_value)
          template <typename TYPE>
          void operator()(const boost::fusion::vector5<std::string,std::string,TYPE*,TYPE,TYPE> &v) const {
               execute(boost::fusion::at_c<0>(v),
                       boost::fusion::at_c<1>(v),
                       boost::fusion::at_c<2>(v),
                       boost::fusion::at_c<3>(v),
                       &boost::fusion::at_c<4>(v));
          }

          //! Function taking a pair consisting of a string and a value of some type
          //! and assigning the corresponding command line option to it. Keeps track
          //! of whether the option has been found so far (called by accumulate)
          //! This version includes: (option-name, description, value_pointer, default_value)
          template <typename TYPE>
          void operator()(const boost::fusion::vector4<std::string,std::string,TYPE*,TYPE> &v) const {
               execute(boost::fusion::at_c<0>(v),
                       boost::fusion::at_c<1>(v),
                       boost::fusion::at_c<2>(v),
                       boost::fusion::at_c<3>(v));
          }

          //! Function taking a pair consisting of a string and a value of some type
          //! and assigning the corresponding command line option to it. Keeps track
          //! of whether the option has been found so far (called by accumulate)
          //! This version includes: (option-name, description, value_pointer)
          template <typename TYPE>
          void operator()(const boost::fusion::vector3<std::string,std::string,TYPE*> &v) const {

               execute(boost::fusion::at_c<0>(v),
                       boost::fusion::at_c<1>(v),
                       boost::fusion::at_c<2>(v));
          }
     };


     //! Generate option help output
     std::string generate_help_output();


     //! Meta-function used to iterate over the typelist provided to generate_output
     //!
     //! \tparam PREV remaining list
     //! \tparam THIS current element
     template <typename PREV, typename THIS>
     struct generate_output_typelist_iteration {
          //! Output type of meta-function
          struct type {
               //! Run time functionality
               static std::string compute(const boost::any &a) {
                    if (a.type() == typeid(THIS)) {
                         return boost::lexical_cast<std::string>(boost::any_cast<THIS>(a));
                    } else if (a.type() == typeid(WrappedEnum<THIS>)) {
                         return boost::lexical_cast<std::string>(boost::any_cast<WrappedEnum<THIS> >(a));
                    } else if (a.type() == typeid(WrappedEnumPointer<THIS>)) {
                         return boost::lexical_cast<std::string>(boost::any_cast<WrappedEnumPointer<THIS> >(a));
                    } 
                    return PREV::compute(a);
               }
          };
     };

     //! Meta-function recursion end-point
     struct generate_output_typelist_base {
          //! Run time functionality
          static std::string compute(const boost::any &a) {
               return "";
          }
     };


     //! Output to string
     template <typename List>
     std::string generate_output() {

          std::string output = "";

          // Iterate over supergroups
          for (unsigned int sg = 0; sg < options_description_super_group_vector.size(); ++sg) {

               std::string super_group_name = options_description_super_group_vector[sg].first;
               std::vector<int> &super_group_indices = options_description_super_group_vector[sg].second;

               if (super_group_indices.size() > 0) {
                    bool all_hidden = true;
                    for (unsigned int sgi = 0; sgi < super_group_indices.size(); ++sgi) {
                         // get index
                         int i = super_group_indices[sgi];
                         
                         if (!options_description_group_vector[i]->hidden)
                              all_hidden = false;
                    }

                    if (!all_hidden) {
                         if (super_group_name == "") {
                              super_group_name = "main";
                         }
                         output += "\n";
                         output += "#####" + std::string(super_group_name.size(), '#') + "#############\n";
                         output += "#### " + boost::to_upper_copy(super_group_name)    + " OPTIONS ####\n";
                         output += "#####" + std::string(super_group_name.size(), '#') + "#############\n\n";
                    } 
               }

               // Iterate over groups
               for (unsigned int sgi = 0; sgi < super_group_indices.size(); ++sgi) {
          
                    // get index
                    int i = super_group_indices[sgi];
                    
                    if (options_description_group_vector[i]->hidden) {
                         continue;
                    }
               
                    output += "### " + options_description_group_vector[i]->title + " ###\n";
                    if (options_description_group_vector[i]->description != "") {
                         output += "# " + replace(options_description_group_vector[i]->description, "\n", "\n# ") + "\n";
                    }

                    const std::vector< boost::shared_ptr<boost::program_options::option_description> >& opt_vec = 
                         options_description_group_vector[i]->options();

                    for (unsigned int j=0; j< opt_vec.size(); j++) {
                         std::string out = "";

                         const boost::program_options::option_description& opt = *(opt_vec[j]);
                         out += opt.long_name() + " = ";

                         const boost::program_options::variable_value& value = option_map[opt.long_name()];

                         const std::type_info& type = value.value().type();

                         if (value.empty())
                              continue;

                         if (type == typeid(std::string) || type == typeid(WrappedStringPointer)) {
                              std::string str_value;
                              if (type == typeid(WrappedStringPointer))
                                   str_value = stringify(*reinterpret_cast<const std::string*>(&value.as<WrappedStringPointer>()));
                              else
                                   str_value = value.as<std::string>();
                              // Output string option WITH quotes
                              // This has been changed back and forth several times.
                              // If you experience problems with the current solution
                              // please let me know. /wb
                              if (str_value[0] != '\"')
                                   out += "\"" + str_value + "\""; 
                              else 
                                   out += str_value; 
                         } else if (type == typeid(WrappedStringPointer)) {
                              out += stringify(*reinterpret_cast<const double*>(&value.as<WrappedDoublePointer>()));                               
                         } else if (type == typeid(bool)) { 
                              out += stringify(value.as<bool>()?"true":"false"); 
                         } else if (type == typeid(int )) { 
                              out += stringify(value.as<int>()); 
                         } else if (type == typeid(unsigned int )) { 
                              out += stringify(value.as<unsigned int>()); 
                         } else if (type == typeid(long unsigned int )) { 
                              out += stringify(value.as<long unsigned int>()); 
                         } else if (type == typeid(PHAISTOS_LONG_LONG)) { 
                              out += stringify(value.as<PHAISTOS_LONG_LONG>()); 
                         } else if (type == typeid(double)) { 
                              out += stringify(value.as<double>()); 
                         } else if (type == typeid(WrappedDoublePointer)) { 
                              out += stringify(*reinterpret_cast<const double*>(&value.as<WrappedDoublePointer>())); 
                         } else if (type == typeid(PlainVector<std::pair<int,int> > ) ) {
                              out += stringify(value.as<PlainVector<std::pair<int,int> > > ());
                         } else if (type == typeid(OptionShorthand) ) {
                              out += stringify(value.as<OptionShorthand>());
                         } else if (type == typeid(std::vector<std::string> ) ) {
                              std::vector<std::string> vec=value.as< std::vector<std::string> >();
                              out += stringify(vec);
                         } else if (type == typeid(CallbackVariable<std::string>) ) {
                              out += stringify((value.as< CallbackVariable<std::string> >()).value());
                         } else {
                              // Iterate over provided type list
                              out += boost::mpl::fold<List, 
                                                      generate_output_typelist_base, 
                                                      generate_output_typelist_iteration<boost::mpl::_,boost::mpl::_> >::type::compute(value.value());
                         }


                         std::ostringstream o;
                         o.setf(std::ios::left);
                         o.width(60);
                         o << out;
                         output += o.str();
                         output += "\t# " + opt.description();
                         output += "\n";
                    }

                    output += "\n";
               }
          }
          return output;
     }
};

}

#endif
