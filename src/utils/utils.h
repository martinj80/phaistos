// utils.h --- Various utility functions
// Copyright (C) 2006-2008 Wouter Boomsma
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


#ifndef UTILS_H
#define UTILS_H

// Set maximum size for fusion vectors
#ifndef FUSION_MAX_VECTOR_SIZE
#define FUSION_MAX_VECTOR_SIZE 30
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#undef BOOST_MPL_LIMIT_VECTOR_SIZE
#define BOOST_MPL_LIMIT_VECTOR_SIZE 30
#endif

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "debug_levels.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <map>
#include <sys/time.h>
#include <time.h>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <cstdlib>
#include <string>

#define SGN(x) (((x)<0) ? -1 : 1)
#define UNINITIALIZED std::numeric_limits<double>::quiet_NaN()

// Always inline attribute - for gcc > 4.0
#if ((__GNUC__==4 && __GNUC_MINOR__>=0) || __GNUC__>4)
  #define PHAISTOS_ALWAYS_INLINE __attribute__((always_inline))
#else
  #define PHAISTOS_ALWAYS_INLINE
#endif

// Check for availability of long long
#ifdef LLONG_MAX
  #define PHAISTOS_LONG_LONG long long
#else
  #define PHAISTOS_LONG_LONG unsigned long
#endif

namespace phaistos {

// make_vector can create conflicts with other implementations. We
// therefore introduce a separate namespace for increased control
namespace vector_utils {

//! Create a vector containing a single element
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val) {
     return std::vector<TYPE>(1, val);
}

//! Create a vector containing two elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2) {
     std::vector<TYPE> vec(2);
     vec[0] = val1;
     vec[1] = val2;
     return vec;
}

//! Create a vector containing three elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3) {
     std::vector<TYPE> vec(3);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     return vec;
}

//! Create a vector containing four elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4) {
     std::vector<TYPE> vec(4);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     return vec;
}

//! Create a vector containing five elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5) {
     std::vector<TYPE> vec(5);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     return vec;
}

//! Create a vector containing six elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5, TYPE val6) {
     std::vector<TYPE> vec(6);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     return vec;
}

//! Create a vector containing seven elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5, TYPE val6, TYPE val7) {
     std::vector<TYPE> vec(7);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     return vec;
}

//! Create a vector containing eight elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5, TYPE val6, TYPE val7, TYPE val8) {
     std::vector<TYPE> vec(8);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     return vec;
}

// Create a vector containing nine elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9) { 
     std::vector<TYPE> vec(9);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     return vec;
}    

// Create a vector containing ten elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10) { 
     std::vector<TYPE> vec(10);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     return vec;
}    

// Create a vector containing eleven elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11) { 
     std::vector<TYPE> vec(11);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     return vec;
}    


// Create a vector containing twelve elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12) {
     std::vector<TYPE> vec(12);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     return vec;
}

// Create a vector containing thirteen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13) {
     std::vector<TYPE> vec(13);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     return vec;
}

// Create a vector containing fourteen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14) {
     std::vector<TYPE> vec(14);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     return vec;
}

// Create a vector containing fifteen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14, TYPE val15) {
     std::vector<TYPE> vec(15);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     vec[14] = val15;
     return vec;
}    

// Create a vector containing sixteen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14, TYPE val15,
                                     TYPE val16) {
     std::vector<TYPE> vec(16);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     vec[14] = val15;
     vec[15] = val16;
     return vec;
}    

// Create a vector containing seventeen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14, TYPE val15,
                                     TYPE val16, TYPE val17) {
     std::vector<TYPE> vec(17);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     vec[14] = val15;
     vec[15] = val16;
     vec[16] = val17;
     return vec;
}    

// Create a vector containing eighteen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14, TYPE val15,
                                     TYPE val16, TYPE val17, TYPE val18) {
     std::vector<TYPE> vec(18);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     vec[14] = val15;
     vec[15] = val16;
     vec[16] = val17;
     vec[17] = val18;
     return vec;
}    

// Create a vector containing nineteen elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5,
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14, TYPE val15,
                                     TYPE val16, TYPE val17, TYPE val18, TYPE val19) {
     std::vector<TYPE> vec(19);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     vec[14] = val15;
     vec[15] = val16;
     vec[16] = val17;
     vec[17] = val18;
     vec[18] = val19;
     return vec;
}    

//! Create a vector containing twenty elements
template <typename TYPE>
inline std::vector<TYPE> make_vector(TYPE val1, TYPE val2, TYPE val3, TYPE val4, TYPE val5, 
                                     TYPE val6, TYPE val7, TYPE val8, TYPE val9, TYPE val10,
                                     TYPE val11, TYPE val12, TYPE val13, TYPE val14, TYPE val15, 
                                     TYPE val16, TYPE val17, TYPE val18, TYPE val19, TYPE val20) {
     std::vector<TYPE> vec(20);
     vec[0] = val1;
     vec[1] = val2;
     vec[2] = val3;
     vec[3] = val4;
     vec[4] = val5;
     vec[5] = val6;
     vec[6] = val7;
     vec[7] = val8;
     vec[8] = val9;
     vec[9] = val10;
     vec[10] = val11;
     vec[11] = val12;
     vec[12] = val13;
     vec[13] = val14;
     vec[14] = val15;
     vec[15] = val16;
     vec[16] = val17;
     vec[17] = val18;
     vec[18] = val19;
     vec[19] = val20;
     return vec;
}

}

//! Check if value is initialized (double version)
//!
//! \param value Value to check
//! \return True if value is initialized
inline bool is_initialized(double value) {
     return !std::isnan(value);
}

//! Check if value is initialized (int version)
//!
//! \param value Value to check
//! \return True if value is initialized
inline bool is_initialized(int value) {
     return (value != std::numeric_limits<int>::max());
}

//! Set a variable to be uninitialized
//!
//! \param value Reference to value
template <typename TYPE>
inline void set_initial_value(TYPE &value) {
     value = TYPE();
}

//! Set a variable to be uninitialized (vector version)
//!
//! \param value Reference to value
template <typename TYPE>
inline void set_initial_value(std::vector<TYPE> &value) {
     for (unsigned int i=0; i<value.size(); ++i) {
          set_initial_value(value[i]);
     }
}

//! Set a variable to be uninitialized (int version)
//!
//! \param value Reference to value
template <>
inline void set_initial_value(int &value) {
     value = std::numeric_limits<int>::max();
}

//! Set a variable to be uninitialized (double version)
//!
//! \param value Reference to value
template <>
inline void set_initial_value(double &value) {
     value = UNINITIALIZED;
}

//! Return uninitialized value for given type
//!
//! \tparam TYPE Type of uninitialized value
//! \return Uninitialized value
template <typename TYPE>
inline TYPE uninitialized() {
     TYPE value;
     set_initial_value(value);
     return value;
}


//! Find value corresponding to key in map
//!
//! \param map Map
//! \param key Key value
//! \param default_val Value to return in case key is not found in map
//! \return value corresponding to key
template <typename KEYTYPE, typename VALTYPE>
inline const VALTYPE &map_lookup(const std::map<KEYTYPE,VALTYPE> &map, KEYTYPE key, const VALTYPE &default_val=VALTYPE()) {
     typename std::map<KEYTYPE,VALTYPE>::const_iterator it = map.find(key);
     if (it != map.end()) {
          return it->second;
     } else {
          return default_val;
     }
}

//! Return maximum of three values
//!
//! \param x Input value 1
//! \param y Input value 2
//! \param z Input value 3
//! \return maximum of x, y, and z
template <typename TYPE>
inline const TYPE &max(const TYPE &x, const TYPE &y, const TYPE &z) {
     return std::max(x, std::max(y, z));
}


//! Get time stamp
//!
//! \return Time in seconds
inline double get_time() {
     double time1;
     struct timeval t1;
     struct timezone tz1;
     gettimeofday(&t1, &tz1);
     time1 = t1.tv_sec+0.000001*t1.tv_usec;
     return time1;
}

//! Compare two doubles for equality
//!
//! \param v1 First value
//! \param v2 Second value
//! \param equality_cutoff Tolerance used when comparing
//! \return True if difference between values is smaller than cutoff
inline bool check_equality(double v1, double v2, double equality_cutoff = 0.01) {
     return(fabs(v1-v2) < equality_cutoff);
}

//! Compare two doubles for equality. Assert if not equal.
//!
//! \param v1 First value
//! \param v2 Second value
//! \param equality_cutoff Tolerance used when comparing
inline void check_equality_assert(double v1, double v2, double equality_cutoff = 0.01) {
     strong_assert(fabs(v1-v2) < equality_cutoff);
}

//! Compare two angles for equality
//!
//! \param v1 First value
//! \param v2 Second value
//! \param equality_cutoff Tolerance used when comparing
//! \return True if difference between values is smaller than cutoff
inline bool check_equality_angle(double v1, double v2, double equality_cutoff = 0.01) {
     double difference = fabs(v1-v2);
     difference = fmin(difference, 2*M_PI-difference);
     return(difference < equality_cutoff);
}

//! Compare two angles for equalityAssert if not equal.
//!
//! \param v1 First value
//! \param v2 Second value
//! \param equality_cutoff Tolerance used when comparing
inline void check_equality_angle_assert(double v1, double v2, double equality_cutoff = 0.01) {
     double difference = fabs(v1-v2);
     difference = fmin(difference, 2*M_PI-difference);
     if (!(difference < equality_cutoff))
          std::cout << v1 << " " << v2 << "\n";
     strong_assert(difference < equality_cutoff);
}

//! Execute Unix command and return output
//!
//! \param command String containing unix command
//! \return Output of command
inline std::string unix_command(std::string command) {
     const int maxLen = 2000;
     char buffer[maxLen];
     FILE *fp = popen(command.data(), "r");
     //char *dummy;
     //dummy = fgets(buffer, maxLen, fp);
     (void)fgets(buffer, maxLen, fp);
     pclose(fp);
     std::string res(buffer);
     return res.erase(res.size());
}



//! Read file into string
//!
//! \param filename Name of file
//! \return Contents of file as string
inline std::string file_to_string(std::string filename) {

     // If no filename is given, return empty string
     if (filename == "")
          return "";

     std::ifstream file(filename.c_str());
     std::stringstream input;
     if (file.is_open()) {     
          input << file.rdbuf();
          file.close();
     } else {
          std::cerr << "Couldn't open file " << filename << ". Exiting.\n";
          exit(1);
     }
     return input.str();     
}

//! Read file line-wise into vector of strings
//!
//! \param filename Name of file
//! \return Contents of file as a vector of strings
inline std::vector<std::string> file_to_string_vector(std::string filename) {
     std::vector<std::string> lines;
     std::string content = file_to_string(filename);
     boost::split(lines, content, boost::is_any_of("\n"), boost::token_compress_on);
     return lines;
}


//! Turn value into string
//!
//! \param value Input value
//! \param scientific Whether to output in scientific notation (for numerical values)
//! \return String version of value
template <typename TYPE>
std::string stringify(TYPE value, bool scientific=false) {
     std::ostringstream o;
     if (scientific)
          o << std::scientific;
     o << value;
     return o.str();
}

//! Turn string into value
//!
//! \tparam TYPE Target type
//! \param str Input string
//! \return Value of TYPE
template <typename TYPE>
TYPE string_to_val(std::string str) {
     std::istringstream buffer(str);
     TYPE val;
     buffer >> val;
     return val;
}


//! Return a string with leading/trailing whitespaced characters clipped 
//!
//! \param str Input string
//! \return Trimmed string
inline std::string strip(std::string str) {
     if(str.empty())
          return str;

     size_t start_index = str.find_first_not_of(" \n");
     size_t end_index = str.find_last_not_of(" \n");
     std::string tmp_string = str;
     str.erase();

     if (start_index!=std::string::npos && end_index!=std::string::npos )
          str = tmp_string.substr(start_index, (end_index-start_index+1));

     return str;
}

//! Replace all instances of a pattern in a string
//!
//! \param str Input string
//! \param pattern Pattern to search force
//! \param rep Replacement string
//! \return New string
inline std::string replace(std::string str, std::string pattern, std::string rep) {
     size_t start_index = 0;
     while (1) {

          // Find start position
          start_index = str.find(pattern, start_index);

          // Break if pattern not found
          if (start_index == std::string::npos) {
               break;
          }

          // Replace string and update startIndex
          str = str.replace(start_index, pattern.size(), rep);
          start_index = start_index + rep.size();

     }
     return str;
}

//! Extract a subvector from a vector
//!
//! \param v Input vector
//! \param start_index Start index
//! \param end_index End index (not included)
//! \return New vector containing specified range
template <typename TYPE>
std::vector<TYPE> slice(std::vector<TYPE> v, int start_index, int end_index=-1) {
     if ((start_index<0 && end_index<0) || start_index >= (int)v.size())
          return v;
     else if (start_index<0)
          start_index = 0;
     else if (end_index < 0)
          end_index = v.size();
     
     std::vector<TYPE> slice(end_index-start_index, v[start_index]);
     for (int i=start_index; i<end_index; i++) {
          slice[i-start_index] = v[i];
     }
     return slice;
}

//! Split string into a vector of strings
//!
//! \param str Input string
//! \param delimiter Delimiting string (separates fields in input)
//! \return Vector of substrings
inline std::vector<std::string> split(std::string str, std::string delimiter=" ") {

     std::vector<std::string> str_vec;
     
     // Find first non-delimiter in string (skip initial delimiters)
     size_t start = str.find_first_not_of(delimiter, 0);

     // Find first delimiter in string
     size_t end = str.find_first_of(delimiter, start);

     while (start != std::string::npos || end != std::string::npos) {

          // Add token
          str_vec.push_back(str.substr(start, end - start));

          // Move start to next non-delimiter
          start = str.find_first_not_of(delimiter, end);

          // Move end to next delimiter
          end = str.find_first_of(delimiter, start);
     }
     return str_vec;
}


//! Test whether file exists in the filesystem
//!
//! \param filename Name of file
//! \return True if file exists
inline bool file_exists (std::string filename) {
     std::ifstream fin;
     fin.open (filename.c_str());
     if (fin.fail()) return false;
     fin.close();
     return true;
}

//! Generate log-filename
//!
//! Expands formatting tags in string
//! \param filename_template Filename with %t and %p tags
//! \param thread_index Which thread this belongs to (optional)
//! \return expanded filename
inline std::string generate_log_filename(const std::string filename_template, 
                                         const int thread_index=0) {
     pid_t pid = getpid();
     std::string filename = filename_template;
     filename = replace(filename, "%p", boost::lexical_cast<std::string>(pid));
     filename = replace(filename, "%t", boost::lexical_cast<std::string>(thread_index));
     return filename;
}

//! Generate log-filename
//!
//! Expands formatting tags in string
//! \param filename_template Filename with %t,%p,%e and %i tags
//! \param energy Energy value to insert instead of %e
//! \param iteration_number Iteration numer to insert instead of %i
//! \param thread_index Which thread this belongs to (optional)
//! \return expanded filename
inline std::string generate_log_filename(const std::string filename_template, 
                                         const double energy,
                                         const int iteration_number,
                                         const int thread_index=0) {
     pid_t pid = getpid();

     char iteration_number_string[100];
     sprintf(iteration_number_string, "%012d", iteration_number);
     char energy_string[100];
     sprintf(energy_string, "%012.4f", energy);

     std::string filename = filename_template;
     filename = replace(filename, "%p", boost::lexical_cast<std::string>(pid));
     filename = replace(filename, "%e", energy_string);
     filename = replace(filename, "%i", iteration_number_string);
     filename = replace(filename, "%t", boost::lexical_cast<std::string>(thread_index));
     return filename;
}

}


// Output operators should be in the same namespace as their containers
namespace std {

//! Overload output operator for vectors of any type
template <typename TYPE>
inline std::ostream &operator<<(std::ostream &o, std::vector<TYPE> v) {
     o << "[";
     for (typename std::vector<TYPE>::iterator it = v.begin(); it != v.end(); ++it) {
          if (it != v.begin())
               o << ",";
          o<<*it;
     }
     o << "]";
     return o;
}


//! Overload input operator for vectors of any type
template <typename TYPE>
inline std::istream &operator>>(std::istream &i, std::vector<TYPE> &v) {

     // Read a parenthesis
     if (i.peek() == '[' || i.peek() == '(')
          i.ignore(1);

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     bool ok_flag = true;

     while (i.good() && (i.peek() != ']' && i.peek() != ')' && i.peek() != '\n')) {

          v.push_back(TYPE());

          // Parse inner value
          // bool entry_fully_parsed = false;
          while (!i.eof() && i.peek() != ' ' && i.peek() != ',' && i.peek() != ']' && i.peek() != ')' && i.peek() != '\n') {
               if (!(i >> v.back())) {

                    phaistos::set_initial_value(v.back());
                    if (!i.eof()) {
                         i.clear();

                         i.ignore(1);
                    }
               }
          }

          // These last peeks should not effect the flags of the stream
          ok_flag = (bool)i;

          // Skip whitespace
          while(i.peek() == ' ')
               i.ignore(1);

          // Skip comma
          if (i.peek() == ',')
               i.ignore(1);
     }

     // Read a parenthesis
     while (i.peek() == ']' || i.peek() == ')' || i.peek() == '\n')
          i.ignore(1);

     if (ok_flag)
          i.clear();

     return i;
}


//! Overload input operator for vectors of string
template <typename TYPE>
inline std::istream &operator>>(std::istream &i, std::vector<std::string> &v) {

     // Read a parenthesis
     if (i.peek() == '[' || i.peek() == '(')
          i.ignore(1);

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     while (i.good() && (i.peek() != ']' && i.peek() != ')')) {

          v.push_back(TYPE());

          // Parse inner value
          while (i.good() && (i.peek() != ',' &&  i.peek() != ']' && i.peek() != ')')) {
               v.back() += std::string(1,i.get());
          }

          // Skip whitespace
          while(i.peek() == ' ')
               i.ignore(1);

          // Skip comma
          if (i.peek() == ',')
               i.ignore(1);
     }

     // Read a parenthesis
     if (i.peek() == ']' || i.peek() == ')')
          i.ignore(1);
     return i;
}

//! Overload output operator for maps of any type
template <typename LTYPE, typename RTYPE>
inline std::ostream &operator<<(std::ostream &o, std::map<LTYPE,RTYPE> m) {
     o << "{";
     for (typename std::map<LTYPE,RTYPE>::iterator it = m.begin(); it != m.end(); ++it) {
          if (it != m.begin())
               o << ", ";
          o << "(" << it->first << ", " << it->second << ")";
     }
     o << "}";
     return o;
}

//! Overload output operator for pairs of any type
template <typename LTYPE, typename RTYPE>
inline std::ostream &operator<<(std::ostream &o, std::pair<LTYPE,RTYPE> p) {
     o << "(" << p.first << "," << p.second << ")";
     return o;
}

//! Overload input operator for pairs of any type
template <typename LTYPE, typename RTYPE>
inline std::istream &operator>>(std::istream &i, std::pair<LTYPE,RTYPE> &p) {

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     // Read a parenthesis
     if (i.peek() == '(' || i.peek() =='[')
          i.ignore(1);

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     // Read first entry
     i >> p.first;

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     // Skip comma
     while(i.peek() == ',')
          i.ignore(1);

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     // Read second entry
     i >> p.second;

     // Skip whitespace
     while(i.peek() == ' ')
          i.ignore(1);

     // Read a parenthesis
     if (i.peek() == ')' || i.peek()==']')
          i.ignore(1);

     return i;
}

}


#endif
