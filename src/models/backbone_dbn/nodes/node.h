// node.h --- Base classes for node objects
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

#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/fusion/include/vector.hpp>

#include <boost/mpl/filter_view.hpp>
#include <boost/fusion/include/nview.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "utils/utils.h"

namespace phaistos {

//! Base class for node objects
template <typename DBN_TYPE, typename SEQUENCE_TYPE=void*>
class Node {

protected:

     //! Functor class for calculating emission probability (used by get_emission())
     struct _get_emission {

          //! Which component to use (typically hidden node value)
          int component_index; 

          //! position in chain
          int pos; 

          //! Whether backup up values should be used
          bool use_backup; 

          //! Reference to variable in which calculation should be written
          double &val;

          //! Constructor
          _get_emission(int component_index, 
                        int pos, 
                        double &val, 
                        bool use_backup=false)
               : component_index(component_index), pos(pos),
                 use_backup(use_backup), val(val){}

          //! Overload () operator. Calculated emission.
          //! \param node DBN node.
          template <typename NODE_TYPE>
          inline void operator()(const NODE_TYPE *node) const {
               val *= node->template emission<NODE_TYPE>(component_index, pos, use_backup);
          }
     };


     //! Functor class for calculating emission log-probability (used by get_emission_log())
     struct _get_emission_log {
          //! Define result type
          // typedef double result_type;
          
          //! Which component to use (typically hidden node value)
          int component_index; 

          //! position in chain
          int pos; 

          //! Whether backup up values should be used
          bool use_backup; 

          //! Reference to variable in which calculation should be written
          double &val;

          //! Constructor          
          _get_emission_log(int component_index, 
                            int pos, 
                            double &val, 
                            bool use_backup)
               : component_index(component_index), pos(pos),
                 use_backup(use_backup), val(val){}

          //! Overload () operator. Calculated log-emission.
          //! \param node DBN node.
          template <typename NODE_TYPE>
          inline void operator()(const NODE_TYPE *node) const {
               val += node->template emission_log<NODE_TYPE>(component_index, pos, use_backup);
          }
     };


public:

     //! Name of node
     std::string name;

     //! Pointer to dbn in which this node is contained
     DBN_TYPE *dbn;
     
     //! Hidden node size
     int h_size;
     
     //! Whether this node is to be resampled
     bool fixed;

     //! Number of components in node
     int size;     

     //! Length of sequence (number of slices)
     int sequence_length;

     //! Sampled/observed values
     SEQUENCE_TYPE sequence;

     //! Backup values (last accepted state)
     SEQUENCE_TYPE sequence_backup;

     //! Which of the sequence values are observed
     bool *sequence_observed;


     //! Constructor
     //!
     //! \param name Node name
     //! \param size Number of components in node
     //! \param dbn Pointer to dbn in which this node is contained
     Node(std::string name, int size, DBN_TYPE *dbn) {
	  this->name = name;
	  this->size = size;		// Set by derived classes
	  this->dbn = dbn;
     }

     //! Copy constructor helper function
     //! \param other Source object.
     void copy(const Node &other) {
	  this->name = other.name;
	  this->size = other.size;
	  this->h_size = other.h_size;
	  this->fixed = other.fixed;
     }

     //! Copy constructor
     //! \param other Source object.
     Node(const Node &other):dbn(other.dbn) {
          copy(other);
     }

     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     Node(const Node &other, DBN_TYPE *dbn):dbn(dbn) {
          copy(other);
     }

     //! Destructor
     ~Node() {
	  delete[] sequence;
	  delete[] sequence_backup;
	  delete[] sequence_observed;     
     }	  


     //! Calculate emission probability for node given component index and position
     //! \param component_index Index of component (typically hidden node index)
     //! \param pos Position in sequence
     //! \param use_backup Whether backup values should be used
     //! \return emission probability
     template <typename DERIVED_CLASS>
     double emission(const int component_index, const int pos, const bool use_backup) const {
          
          // Cast this pointer to derived class (avoid virtual functions for efficiency)
          DERIVED_CLASS *this_pointer = (DERIVED_CLASS *)this;

          double value = 1.0;
          if (fixed && sequence_observed[pos]) {

               if (use_backup) {
                    value = this_pointer->get_likelihood(component_index, sequence_backup[pos]);
               } else {
                    value = this_pointer->get_likelihood(component_index, sequence[pos]);
               }
          } 
          return value;
     }


     //! Calculate log-emission probability for node given component index and position
     //! \param component_index Index of component (typically hidden node index)
     //! \param pos Position in sequence
     //! \param use_backup Whether backup values should be used
     //! \return emission log-probability
     template <typename DERIVED_CLASS>
     double emission_log(const int component_index, const int pos, const bool use_backup) const {
          
          // Cast this pointer to derived class (avoid virtual functions for efficiency)
          DERIVED_CLASS *this_pointer = (DERIVED_CLASS *)this;

          double value = 0.0;
          if (fixed && sequence_observed[pos]) {
               if (use_backup)
                    value = this_pointer->get_log_likelihood(component_index, sequence_backup[pos]);
               else 
                    value = this_pointer->get_log_likelihood(component_index, sequence[pos]);
          } 
          return value;
     }


     //! Functor class returning emission state of a node.
     struct _get_emission_state {
          //! Reference to variable in which results should be written
          bool &value;

          //! Constructor
          _get_emission_state(bool &value)
               : value(value){}

          //! Overload () operator. Retrieve emission state from node.
          //! \param node DBN node.          
          template <typename TYPE>
          inline void operator()(TYPE *node) const {
               value = (value && node->fixed);
          }
     };

};


//! Node base class for 1-dimensional nodes
template <typename TYPE, typename DBN_TYPE, int PARENT_NODE=0>
class Node_1D: public Node<DBN_TYPE,TYPE *> {
public:

     //! The type of value modeled by this node 
     typedef TYPE Type;

     //! Node index of parent node
     typedef typename boost::mpl::int_<PARENT_NODE> ParentNode;

     //! Vector of node indices of children
     typedef boost::fusion::vector<> ChildNodes;

     //! Constructor
     //! \param name of node
     //! \param dbn Pointer to dbn in which this node is contained
     Node_1D(std::string name, DBN_TYPE *dbn)
          : Node<DBN_TYPE,TYPE*>(name,1,dbn) {

	  this->init(0);
     }

     //! Copy constructor helper function
     //! \param other Source object.
     void copy(const Node_1D &other) {
          init(other.sequence_length);
          Node<DBN_TYPE,TYPE*>::copy(other);
          for (int i=0; i<this->sequence_length; i++) {
               this->sequence[i] = other.sequence[i];
               this->sequence_backup[i] = other.sequence_backup[i];
               this->sequence_observed[i] = other.sequence_observed[i];
          }     
     }

     //! Copy constructor
     //! \param other Source object.
     Node_1D(const Node_1D &other)
          : Node<DBN_TYPE,TYPE*>(other) {
          copy(other);
     }

     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     Node_1D(const Node_1D &other, DBN_TYPE *dbn)
          : Node<DBN_TYPE,TYPE*>(other, dbn) {
          copy(other);
     }

     //! Initialize sequence data
     //! \param length Length of sequence
     //! \param complete_reinitialization Whether to reset all sequences
     void init(int length, bool complete_reinitialization=true) {

          if (!complete_reinitialization)
               return;

	  this->sequence_length = length;

	  this->fixed = false;

	  if (length == 0) {
	       this->sequence = NULL;
	       this->sequence_backup = NULL;
	       this->sequence_observed = NULL;
	  } else {
	       this->sequence = new TYPE[length];
	       this->sequence_backup = new TYPE[length];
	       this->sequence_observed = new bool[length];
	       for (int i=0; i<this->sequence_length; i++) {
                    set_initial_value(this->sequence[i]);
		    set_initial_value(this->sequence_backup[i]);
		    this->sequence_observed[i] = true;
	       }
	  }
     }

     
     //! Accept proposed sequence
     void accept(int start_index=0, int end_index=-1) {
          if (end_index<0)
               end_index = this->sequence_length;

	  for (int i=start_index; i<end_index; i++) {
	       this->sequence_backup[i] = this->sequence[i];
	  }
     }

     //! Reject proposed sequence
     void reject(int start_index=0, int end_index=-1) {
          if (end_index<0)
               end_index = this->sequence_length;

	  for (int i=start_index; i<end_index; i++) {
	       this->sequence[i] = this->sequence_backup[i];
	  }
     }


     //! Set value at specific position
     //! \param position Position in sequence
     //! \param value Value to set
     //! \param set_backup Whether value should be written to backup as well
     void set(int position, TYPE value, bool set_backup=false) {
	  this->sequence[position] = value;

	  if (set_backup) {
	       this->sequence_backup[position] = value;
	  }
     }
     
     //! Set input sequence.
     //!
     //! \param input_sequence Vector of input values
     //! \param input_sequence_observed Vector of booleans specifying which of the input values are observed
     //! \param start_index Which slice to start at
     //! \param fix_emission indicates whether the node should be made fixed after the sequence has been set
     //! this is often what you want, since the node will thereby become an input node. In some cases, it
     //! however convenient to initialize the values of an output sequence.
     //! \param set_observed Whether to update observed values
     void set_sequence(const std::vector<TYPE> &input_sequence,
		       const std::vector<bool> &input_sequence_observed=std::vector<bool>(),
                       int start_index=0,
		       bool fix_emission=true,
                       bool set_observed=true) {

	  if (input_sequence.size() > 0) {

               // Initialize DBN if this is is the first sequence that is set
               if (this->dbn->sequence_length == 0) {
                    this->dbn->init(input_sequence.size());
               }

	       for (unsigned int i=0; i<input_sequence.size(); i++) {
		    if (input_sequence_observed.size() > 0) {
			 this->sequence_observed[i+start_index] = input_sequence_observed[i];
		    } else {
                         if (set_observed) {
                              if (!is_initialized(input_sequence[i]))
                                   this->sequence_observed[i+start_index] = false;
                              else
                                   this->sequence_observed[i+start_index] = true;
                         }
		    }
		    this->set(i+start_index, input_sequence[i+start_index], true);
	       }
	       if (fix_emission) {
		    this->fixed = true;          // turn this node into an output node
	
                    // The hidden nodes should be resampled after this
                    this->dbn->inconsistent_regions = 
                         std::vector<typename DBN_TYPE::ModifiedRegion>(
                              1, typename DBN_TYPE::ModifiedRegion(-1,-1));
	       }
	  }
     }
     
     //! Get string representation of sequence
     //! \param value_to_str_functor Functor used to translate values into strings
     //! \return Output string
     template <typename VAL_TO_STR_FUNCTOR>
     std::string get_sequence(VAL_TO_STR_FUNCTOR *value_to_str_functor) {
	  // Make sure that hidden node sequence is initialized
          if (this->dbn->inconsistent_regions.size() > 0) {
	       // Reinitialize sequences
	       this->dbn->init_sequences();
	  }

          std::string output = "";
          for (int i=0; i<this->sequence_length; ++i) {
               output += value_to_str_functor(this->sequence[i]);
          }
	  return output;
     }

     //! Get sequence (directly)
     //! \return pointer to internal sequence representation
     TYPE *get_sequence() {
	  // Make sure that hidden node sequence is initialized
          if (this->dbn->inconsistent_regions.size() > 0) {
	       // Reinitialize sequences
	       this->dbn->init_sequences();
	  }
	  return this->sequence;
     }

     //! Get sequence as a vector (making a copy)
     //! \return Sequence as a vector of values
     std::vector<TYPE> get_sequence_vector(int start_index=-1, int end_index=-1) {
	  if(end_index < 0 || end_index > this->sequence_length)
	       end_index = this->sequence_length;
	  
	  if(start_index < 0 || start_index > this->sequence_length)
	       start_index = 0;

	  TYPE *sequence = this->get_sequence();
	  std::vector<TYPE> sequence_vec;
	  for (int i=start_index; i<end_index; i++) {
	       sequence_vec.push_back(sequence[i]);
	  }
	  return sequence_vec;
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, Node_1D &n) {
          o << "Node: " << n.name << "\n";

          o << "Sequence: [";
          for (int i=0; i<n.sequence_length; ++i) {
               if (i>0)
                    o << ",";
               o << (is_initialized(n.sequence[i])?boost::lexical_cast<std::string>(n.sequence[i]):"UNINITIALIZED");
          }
          o << "]\n";
          o << "Observed: [";
          for (int i=0; i<n.sequence_length; ++i) {
               if (i>0)
                    o << ",";
               o << n.sequence_observed[i];
          }
          o << "]\n";
          o << "Fixed: " << (n.fixed ? "true\n" : "false\n");
          return o;
     }


     //! Local iterator class
     class Iterator {

          //! Reference to node
          Node_1D &node;

     public:

          //! Current dimension index within node (always zero)
          const static int dimension_index=0;

          //! Current index in sequence (slice)
          int sequence_index;

          //! The end index
          int end_index;

          //! Constructor
          //! \param node Node to iterate over
          //! \param start_index Sequence start index
          //! \param end_index Sequence end index 
          Iterator(Node_1D &node, int start_index=0, int end_index=-1)
               : node(node),
                 sequence_index(start_index),
                 end_index(((end_index<0) ? node.sequence_length : end_index)) {}

          //! Overload * operator
          //! \return reference to value in node sequence
          typename Node_1D::Type &operator*() {
               return node.sequence[sequence_index];
          }

          //! Checks whether entry is observed
          bool &observed_entry() {
               return node.sequence_observed[sequence_index];
          }

          //! Increment
          //! \return reference to iterator
          Iterator &operator++() {
               sequence_index++;
               return *this;
          }

          //! Check whether we are at the end point
          bool end() {
               return sequence_index>=end_index;
          }
     };     
};



//! Node base class for 2 dimensional nodes
template <typename TYPE, typename DBN_TYPE, int PARENT_NODE=0>
class Node_2D: public Node<DBN_TYPE, TYPE **> {
public:

     //! The type of value modeled by this node 
     typedef TYPE Type;

     //! Node index of parent node     
     typedef typename boost::mpl::int_<PARENT_NODE> ParentNode;

     //! Vector of node indices of children     
     typedef boost::fusion::vector<> ChildNodes;

     //! Constructor
     //! \param name of node
     //! \param size Number of components in node
     //! \param dbn Pointer to dbn in which this node is contained
     Node_2D(std::string name, int size, DBN_TYPE *dbn)
          : Node<DBN_TYPE,TYPE**>(name,size,dbn) {

	  this->init(0);
     }

     //! Copy constructor helper function
     //! \param other Source object.
     void copy(const Node_2D &other) {
	  this->init(other.sequence_length);
          Node<DBN_TYPE,TYPE**>::copy(other);
	  for (int i=0; i<this->sequence_length; i++) {
	       for (int j=0; j<this->size; j++) {
		    this->sequence[i][j] = other.sequence[i][j];
		    this->sequence_backup[i][j] = other.sequence_backup[i][j];
	       }
	       this->sequence_observed[i] = other.sequence_observed[i];
	  }     
     }

     //! Copy constructor
     //! \param other Source object.
     Node_2D(const Node_2D &other)
          : Node<DBN_TYPE,TYPE**>(other) {
          copy(other);
     }

     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     Node_2D(const Node_2D &other, DBN_TYPE *dbn)
          : Node<DBN_TYPE,TYPE**>(other,dbn) {
          copy(other);
     }

     //! Destructor
     ~Node_2D() {
          for (int i=0; i<this->sequence_length; i++) {
               delete[] this->sequence[i];
               delete[] this->sequence_backup[i];
          }
     }

     //! Initialize sequence data
     //! \param length Length of sequence
     //! \param complete_reinitialization Whether to reset all sequences
     void init(int length, bool complete_reinitialization=true) {

          if (!complete_reinitialization)
               return;

	  this->sequence_length = length;

	  this->fixed = false;     

	  if (length == 0) {
	       this->sequence = NULL;
	       this->sequence_backup = NULL;
	       this->sequence_observed = NULL;
	  } else {
	       this->sequence = new TYPE*[length];
	       this->sequence_backup = new TYPE*[length];
	       for (int i=0; i<length; i++) {
		    this->sequence[i] = new TYPE[this->size];
		    this->sequence_backup[i] = new TYPE[this->size];
		    for (int j=0; j<this->size; j++) {
                         set_initial_value(this->sequence[i][j]);
                         set_initial_value(this->sequence_backup[i][j]);
		    }
	       }	  
	       this->sequence_observed = new bool[length];
	       for (int i=0; i<length; i++) {
		    this->sequence_observed[i] = true;
	       }
	  }
     }
	  
     //! Accept proposed sequence
     void accept(int start_index=0, int end_index=-1) {
          if (end_index<0)
               end_index = this->sequence_length;

	  for (int i=start_index; i<end_index; i++) {
	       for (int j=0; j<this->size; j++) {
		    this->sequence_backup[i][j] = this->sequence[i][j];
	       }
	  }
     }

     //! Reject proposed sequence
     void reject(int start_index=0, int end_index=-1) {
          if (end_index<0)
               end_index = this->sequence_length;

	  for (int i=start_index; i<end_index; i++) {
	       for (int j=0; j<this->size; j++) {
		    this->sequence[i][j] = this->sequence_backup[i][j];
	       }
	  }
     }

     //! Set values at specific position
     //! \param position Position in sequence
     //! \param input_values Vector of values to set
     //! \param set_backup Whether value should be written to backup as well
     void set(int position, std::vector<TYPE> input_values, bool set_backup=false) {
	  for (int j=0; j<this->size; j++) {
	       this->sequence[position][j] = input_values[j];
	       
	       if (set_backup)
		    this->sequence_backup[position][j] = input_values[j];
	  }
     }


     //! Set values at specific position - only useful for 2D nodes
     //! \param position Position in sequence
     //! \param value1 First value to set
     //! \param value2 Second value to set
     //! \param set_backup Whether value should be written to backup as well
     void set(int position, TYPE value1, TYPE value2, bool set_backup=false) {
	  this->sequence[position][0] = value1;
	  this->sequence[position][1] = value2;
	       
	  if (set_backup) {
	       this->sequence_backup[position][0] = value1;
	       this->sequence_backup[position][1] = value2;
	  }
     }
     
     
     //! Set input sequence.
     //!
     //! \param input_sequence Vector of input values
     //! \param input_sequence_observed Vector of booleans specifying which of the input values are observed
     //! \param start_index Which slice to start at
     //! \param fix_emission indicates whether the node should be made fixed after the sequence has been set
     //! this is often what you want, since the node will thereby become an input node. In some cases, it
     //! however convenient to initialize the values of an output sequence.
     //! \param set_observed Whether to update observed values
     void set_sequence(const std::vector<std::vector<TYPE> > &input_sequence,
		       const std::vector<bool> &input_sequence_observed=std::vector<bool>(),
                       int start_index=0,
		       bool fix_emission=true,
                       bool set_observed=true) {
	  if (input_sequence.size() > 0 && input_sequence[0].size() > 0) {

               // Initialize DBN if this is is the first sequence that is set
               if (this->dbn->sequence_length == 0) {
                    this->dbn->init(input_sequence.size());
               }

	       for (unsigned int i=0; i<input_sequence.size(); i++) {
                    if (set_observed) {
                         bool observed;
                         if (input_sequence_observed.size() > 0) {
                              observed = input_sequence_observed[i];
                         } else {
                              bool nanFound = false;
                              for (unsigned int j=0; j<input_sequence[i+start_index].size(); j++) {
                                   if (!is_initialized(input_sequence[i+start_index][j])) {
                                        nanFound = true;
                                   }
                              }
                              if (nanFound)
                                   observed = false;
                              else
                                   observed = true;			 
                         }
                         this->sequence_observed[i+start_index] = observed;     
                    }	       
		    this->set(i+start_index, input_sequence[i+start_index], true);
	       }
	       if (fix_emission) {
		    this->fixed = true;          // turn this node into an output node

                    // The hidden nodes should be resampled after this
                    this->dbn->inconsistent_regions = 
                         std::vector<typename DBN_TYPE::ModifiedRegion>(
                              1, typename DBN_TYPE::ModifiedRegion(-1,-1));
	       }
	  }
     }


     //! Set input sequence - only useful for 2D nodes
     //!
     //! \param input_sequence1 Vector of first entry input values
     //! \param input_sequence2 Vector of second entry input values
     //! \param input_sequence_observed Vector of booleans specifying which of the input values are observed
     //! \param start_index Which slice to start at
     //! \param fix_emission indicates whether the node should be made fixed after the sequence has been set
     //! this is often what you want, since the node will thereby become an input node. In some cases, it
     //! however convenient to initialize the values of an output sequence.
     //! \param set_observed Whether to update observed values
     void set_sequence(const std::vector<TYPE> &input_sequence1,
                       const std::vector<TYPE> &input_sequence2,
		       const std::vector<int> &input_sequence_observed=std::vector<int>(),
                       int start_index=0,
		       bool fix_emission=true,
                       bool set_observed=true) {
	  if (input_sequence1.size() > 0 && input_sequence2.size() > 0) {

               // Initialize DBN if this is is the first sequence that is set
               if (this->dbn->sequence_length == 0) {
                    this->dbn->init(std::max(input_sequence1.size(), input_sequence2.size()));
               }

	       for (unsigned int i=0; i<input_sequence1.size(); i++) {

                    if (set_observed) {
                         bool observed = true;
                         if (input_sequence_observed.size() > 0) {
                              observed = input_sequence_observed[i];
                         } else {
                              if (std::isnan(input_sequence1[i]) || std::isnan(input_sequence2[i])) {
                                   observed = false;
                              }
                         }
                         this->sequence_observed[i+start_index] = observed;     
                    }
                    std::vector<double> pair(2);
                    pair[0] = input_sequence1[i];
                    pair[1] = input_sequence2[i];

		    this->set(i+start_index, pair, true);
	       }
	       if (fix_emission) {
		    this->fixed = true;          // turn this node into an output node
		    this->dbn->initialized = false;    // The hidden nodes should be resampled after this
	       }
	  }
     }

     
     //! Get sequence (directly)
     //! \return pointer to internal sequence representation
     TYPE **get_sequence() {
	  // Make sure that hidden node sequence is initialized
          if (this->dbn->inconsistent_regions.size() > 0) {

	       // Reinitialize sequences
	       this->dbn->init_sequences();
	  } 
	  return this->sequence;
     }


     //! Get sequence as a vector (making a copy)
     //! \return Sequence as a vector of values
     std::vector<std::vector<TYPE> > get_sequence_vector(int start_index=-1, int end_index=-1) {
	  if(end_index < 0 || end_index > this->sequence_length)
	       end_index = this->sequence_length;
	  
	  if(start_index < 0 || start_index > this->sequence_length)
	       start_index = 0;
	  
	  TYPE **sequence = this->get_sequence();
	  std::vector<std::vector<TYPE> > sequence_vec;
	  for (int i=start_index; i<end_index; i++) {
	       sequence_vec.push_back(std::vector<TYPE>());
	       for (int j=0; j<this->size; j++) {
		    sequence_vec[i-start_index].push_back(sequence[i][j]);
	       }
	  }
	  return sequence_vec;
     }     


     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, Node_2D &n) {
          o << "Node: " << n.name << "\n";

          o << "Sequence: [";
          for (int i=0; i<n.sequence_length; ++i) {
               if (i>0)
                    o << ",";
               o << "[";
               for (int j=0; j<n.size; ++j) {
                    if (j>0)
                         o << ",";
                    o << (is_initialized(n.sequence[i][j])?boost::lexical_cast<std::string>(n.sequence[i][j]):"UNINITIALIZED");
               }
               o << "]";
          }
          o << "]\n";
          o << "Observed: [";
          for (int i=0; i<n.sequence_length; ++i) {
               if (i>0)
                    o << ",";
               o << n.sequence_observed[i];
          }
          o << "]\n";
          o << "Fixed: " << (n.fixed ? "true\n" : "false\n");
          return o;
     }

     
     //! Local iterator class
     class Iterator {

          //! Reference to node
          Node_2D &node;

     public:

          //! Current dimension index within node
          int dimension_index;

          //! Current index in sequence (slice)
          int sequence_index;

          //! The end index          
          int end_index;
          
          //! Constructor
          //! \param node Node to iterate over
          //! \param start_index Sequence start index
          //! \param end_index Sequence end index 
          Iterator(Node_2D &node, int start_index=0, int end_index=-1)
               : node(node),
                 dimension_index(0),
                 sequence_index(start_index),
                 end_index(((end_index<0) ? node.sequence_length : end_index)) {}


          //! Overload * operator
          //! \return reference to value in node sequence
          typename Node_2D::Type &operator*() {
               return node.sequence[sequence_index][dimension_index];
          }

          //! Checks whether entry is observed
          bool &observed_entry() {
               return node.sequence_observed[sequence_index];
          }

          //! Increment
          //! \return reference to iterator
          Iterator &operator++() {
               if (dimension_index == node.size-1) {
                    dimension_index = 0;
                    sequence_index++;
               } else {
                    dimension_index++;
               }
               return *this;
          }

          //! Check whether we are at the end point
          bool end() {
               return sequence_index>=end_index;
          }
     };
};

}

#endif     
