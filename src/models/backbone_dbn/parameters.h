// parameters.h --- DBN parameter parsing
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


#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include "utils/utils.h"

namespace phaistos {

//! Parameter class for BackboneDbn models
class Parameters {

     //! Read integer values
     std::map<std::string, std::map<std::string, std::map<std::string, int> > > int_val;

     //! Read double values
     std::map<std::string, std::map<std::string, std::map<std::string, double> > > double_val;

     //! Read double vector values
     std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > double_vec_val;

     //! Read double matrix values
     std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::vector<double> > > > > double_matrix_val;
     
public:
     //! filename
     std::string filename;
     
     //! Read single value from input stream
     //!
     //! \tparam TYPE type of the value to be read
     //! \param input Input stream
     //!
     //! \return value
     template <typename TYPE>
     TYPE read_single(std::istringstream &input);


     //! Read values from input stream into vector
     //!
     //! \tparam TYPE type of the value to be read
     //! \param input Input stream
     //! \param size Number of values to read
     //!
     //! \return value
     template <typename TYPE>
     std::vector<TYPE> read_vector(std::istringstream &input, int size);


     //! Read values from input stream into matrix
     //!
     //! \tparam TYPE type of the value to be read
     //! \param input Input stream
     //! \param m Number of rows to read
     //! \param n Number of columns to read
     //!
     //! \return value
     template <typename TYPE>
     std::vector<std::vector<TYPE> > read_matrix(std::istringstream &input, int m, int n);


     //! Add entry to parameter container
     //!
     //! \tparam TYPE type of the value
     //! \param map_val Map to insert value into
     //! \param val value to be inserted
     //! \param node_type Type of node (e.g. DISCRETE)
     //! \param node_name Name of node (e.g. AA)
     //! \param value_type Type of value (e.g. EMISSION)
     template <typename TYPE>
     void add_entry(std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > > *map_val,
                    TYPE &val, std::string node_type, std::string node_name, std::string value_type) {
	  if (map_val->find(node_type) == map_val->end()) {
	       (*map_val)[node_type] = std::map<std::string, std::map<std::string, TYPE> >();
	  }
	  if ((*map_val)[node_type].find(node_name) == (*map_val)[node_type].end()) {
	       (*map_val)[node_type][node_name] = std::map<std::string, TYPE>();
	  }

	  (*map_val)[node_type][node_name][value_type] = val;

     }

     //! Get the relevant map
     std::map<std::string, std::map<std::string, std::map<std::string, int> > > *get_map(std::map<std::string, std::map<std::string, std::map<std::string, int> > > *tmp) {
	  return &int_val;
     }

     //! Get the relevant map
     std::map<std::string, std::map<std::string, std::map<std::string, double> > > *get_map(std::map<std::string, std::map<std::string, std::map<std::string, double> > > *tmp) {
	  return &double_val;
     }

     //! Get the relevant map
     std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > *get_map(std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > *tmp) {
	  return &double_vec_val;
     }

     //! Get the relevant map
     std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::vector<double> > > > > *get_map(std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::vector<double> > > > > *tmp) {
	  return &double_matrix_val;
     }
     
     //! Get parameters for a specific Nodetype,Nodename,value_type. 
     //!
     //! \tparam TYPE value type
     //! \param node_type Type of node (e.g. DISCRETE)
     //! \param node_name Name of node (e.g. AA)
     //! \param value_type Type of value (e.g. EMISSION)
     //! \param raise_exception Indicates whether to raise an exception when a label is not found.
     //!
     //! \return value
     template <typename TYPE>     
     TYPE get_node_parameters(std::string node_type, std::string node_name, std::string value_type, bool raise_exception=true) {
	  std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > > *map_val = get_map((std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > > *)0);
	  TYPE res;

	  if (map_val->find(node_type) == map_val->end()) {
	       if (raise_exception) {
		    std::ostringstream s;	       
		    s << "Error (get_node_parameters): node_type " << node_type << " not found.";
		    throw Parameters::Error(s.str());
	       }
	       return res;
	  }
	  if ((*map_val)[node_type].find(node_name) == (*map_val)[node_type].end()) {
	       if (raise_exception) {
		    std::ostringstream s;	       
		    s << "Error (get_node_parameters): node_name " << node_name << " not found.";
		    throw Parameters::Error(s.str());
	       }
	       return res;
	  }
	  if ((*map_val)[node_type][node_name].find(value_type) == (*map_val)[node_type][node_name].end()) {
	       if (raise_exception) {
		    std::ostringstream s;	       
		    s << "Error (get_node_parameters): value_type " << value_type << " not found.";
		    throw Parameters::Error(s.str());
	       }
	       return res;
	  }

	  return (*map_val)[node_type][node_name][value_type];
     }

     //! Set parameters for a specific Nodetype,Nodename,value_type. 
     //!
     //! \tparam TYPE value type
     //! \param node_type Type of node (e.g. DISCRETE)
     //! \param node_name Name of node (e.g. AA)
     //! \param value_type Type of value (e.g. EMISSION)
     //! \param val value
     template <typename TYPE>     
     void set_node_parameters(std::string node_type, std::string node_name, std::string value_type, TYPE &val) {
	  std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > > *map_val = get_map((std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > > *)0);

	  add_entry(map_val, val, node_type, node_name, value_type);
     }

     
     //! Constructor. Attempts to read filename first. If unsuccesful, it reads the string data.
     //! \param filename Name of file
     //! \param data Data string. Read if reading from filename is unsuccesful.
     Parameters(const std::string filename, const std::string &data="")
          : filename(filename) {
          if (filename != "") {
               std::string data = file_to_string(filename);
               std::istringstream input(data);
               read_parameters(input);
          } else {
               std::istringstream input(data);
               read_parameters(input);
          }
     }

     //! Read parameters from input stream.
     //! \param input Input stream
     void read_parameters(std::istringstream &input) {
	  this->filename = filename;

          std::string line;
          while (getline(input, line)) {

               std::vector<std::string> words;
               boost::split(words, line, boost::is_any_of("\t "), boost::token_compress_on);
               
               if (words.size()==0) {continue;}
	       std::string label = words[0];

	       if (label[0] == '#') {
		    std::vector<std::string> strVec = split(label.substr(1,std::string::npos), "_");
		    if (strVec.size() > 4) {

			 std::string type = strVec[0];
			 std::string size = strVec[1];
			 std::string node_type = strVec[2];
			 std::string node_name = strVec[3];
			 std::string value_type = strVec[4];

			 std::vector<std::string> sizeVec = split(size, "x");

			 if (type=="INT") {
			      int val = read_single<int>(input);
			      add_entry(&int_val, val, node_type, node_name, value_type);
			 } else if (type=="DOUBLE") {
			      if (sizeVec.size() == 1) {
				   double val = read_single<double>(input);
				   add_entry(&double_val, val, node_type, node_name, value_type);
			      } else if (string_to_val<int>(sizeVec[0]) == 1) {
				   std::vector<double> val = read_vector<double>(input, string_to_val<int>(sizeVec[1]));
				   add_entry(&double_vec_val, val, node_type, node_name, value_type);
			      } else {
				   int m = string_to_val<int>(sizeVec[0]);
				   int n = string_to_val<int>(sizeVec[1]);
				   std::vector<std::vector<double> > val = read_matrix<double>(input,m,n);
				   add_entry(&double_matrix_val, val, node_type, node_name, value_type);
			      }
			 }		    
		    }	  
	       }
	  }     
     }


     //! Output style for integer values
     //! \param val value
     //! \param prefix pointer allowing extraction of prefix
     //! \return output string
     std::string output(int val, std::string *prefix) const {
	  *prefix = "INT_1_";
	  return stringify(val, true) + "\n";
     }

     //! Output style for double values
     //! \param val value
     //! \param prefix pointer allowing extraction of prefix
     //! \return output string
     std::string output(double val, std::string *prefix) const {
	  *prefix = "DOUBLE_1_";
	  return stringify(val, true) + "\n";
     }

     //! Output style for vector<double> values
     //! \param val value
     //! \param prefix pointer allowing extraction of prefix
     //! \return output string
     std::string output(std::vector<double> val, std::string *prefix) const {
	  *prefix = "DOUBLE_1x" + stringify(val.size()) + "_";
	  std::string out;
	  for (unsigned int i=0; i<val.size(); i++)
	       out += stringify(val[i], true) + " ";
	  out += "\n";
	  return out;
     }

     //! Output style for vector<vector<double> > values
     //! \param val value
     //! \param prefix pointer allowing extraction of prefix
     //! \return output string
     std::string output(std::vector<std::vector<double> > val, std::string *prefix) const {
	  *prefix = "DOUBLE_" + stringify(val.size()) + "x" + stringify(val[0].size()) + "_";
	  std::string out;	  
	  for (unsigned int i=0; i<val.size(); i++) {
	       for (unsigned int j=0; j<val[i].size(); j++) {
		    out += stringify(val[i][j], true) + " ";
	       }
	       out += "\n";
	  }
	  return out;
     }
     
     //! Output parameters in self-readable format
     //! \tparam TYPE value type
     //! \param out Output stream
     //! \param map_val selected map value
     template <typename TYPE>     
     void output(std::ostream &out,
		 std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > > map_val) const {

	  typename std::map<std::string, std::map<std::string, std::map<std::string, TYPE> > >::iterator it1;
	  for (it1 = map_val.begin(); it1 != map_val.end(); ++it1) {
	       std::string label1 = it1->first + "_";
	       typename std::map<std::string, std::map<std::string, TYPE> >::iterator it2;
	       for (it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
		    std::string label2 = label1 + it2->first + "_";
		    typename std::map<std::string, TYPE>::iterator it3;
		    for (it3 = it2->second.begin(); it3 != it2->second.end(); ++it3) {
			 std::string label3 = label2 + it3->first;
			 std::string prefix;
			 std::string value = output(it3->second, &prefix);
			 out << "#" + prefix + label3 << "\n";
			 out << value << "\n";
		    }
	       }
	  }
     }
     
     
     //! Output operator
     friend std::ostream& operator<<(std::ostream& out, const Parameters &parameters) {
	  parameters.output(out, parameters.int_val);
	  parameters.output(out, parameters.double_val);
	  parameters.output(out, parameters.double_vec_val);
	  parameters.output(out, parameters.double_matrix_val);
	  return out;
     }


     //! Save parameters to file
     //! \param filename
     void save_to_file(std::string filename="") {
	  if (filename == "")
	       filename = this->filename;

          if (filename == "")
               return;
	  std::ofstream f;
          f.open(filename.c_str());
	  f << *this << "\n";
          f.close();
     }

     
     //! Local exception class
     class Error: public std::exception {
     private:
          //! message string
	  std::string msg;
     public:
          //! Override what function
	  virtual const char* what() const throw() {
	       return(msg.c_str());
	  }
	  
          //! Constructor
          //! \param s Error string
	  Error(const std::string &s = "") : msg(s) {}
	  
          //! Destructor
	  ~Error() throw() {}
     };

};


template <typename TYPE>
TYPE Parameters::read_single(std::istringstream &input) {
     std::string line;
     TYPE val;
     if (getline(input, line)) {
          boost::trim(line);
          std::vector<std::string> words;
          boost::split(words, line, boost::is_any_of("\t "), boost::token_compress_on);
          assert(words.size() == 1);
          val = boost::lexical_cast<TYPE>(words[0]);
     } else {
          std::cerr << "Error in read_single(): cannot read value from file\n";
          exit(1);
     }
     return val;
}

template <typename TYPE>
std::vector<TYPE> Parameters::read_vector(std::istringstream &input, int size) {
     std::vector<TYPE> vec;

     std::string line;
     TYPE val;
     if (getline(input, line)) {
          boost::trim(line);
          std::vector<std::string> words;
          boost::split(words, line, boost::is_any_of("\t "), boost::token_compress_on);
     
          if (words.size() > 0 and !(words.size()==1 and (words[0]=="NULL"))) {
               if ((int)words.size() != size) {
                    fprintf(stderr, "Error: Specified size and length of input don't match.\n");
                    assert(false);
               }
               for (unsigned int i=0; i<words.size(); i++) {
                    val = boost::lexical_cast<TYPE>(words[i]);
                    vec.push_back(val);
               }
          }
     } else {
          std::cerr << "Error in read_vector(): cannot read value from file\n";
          assert(false);
     }
     return vec;
}

template <typename TYPE>
inline std::vector<std::vector<TYPE> > Parameters::read_matrix(std::istringstream &input, int m, int n) {
     std::vector<std::vector<TYPE> > matrix;

     std::string line;
     int lines = 0;
     std::vector<std::string> words;
     while (lines < m && getline(input, line)) {
          boost::trim(line);
          std::vector<std::string> words;
          boost::split(words, line, boost::is_any_of("\t "), boost::token_compress_on);

	  std::vector<TYPE> rowVec;
          if ((words[0]=="NULL") || (int)words.size() != n) {
               fprintf(stderr, "Error: matrix input format corrupt"); 
               assert(false);
          }
          for (int i=0; i<n; i++) {
	       TYPE val = boost::lexical_cast<TYPE>(words[i]);
	       rowVec.push_back(val);
          }
          lines++;
	  matrix.push_back(rowVec);
     }
     return matrix;
}

}

#endif
