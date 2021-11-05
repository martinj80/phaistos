// backbone_dbn.h --- probabilistic model for the protein backbone
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


#ifndef BACKBONE_DBN_H
#define BACKBONE_DBN_H

#include <iostream>

#include "utils/utils.h"
#include "utils/settings.h"

#include <boost/mpl/vector.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/is_sequence.hpp>
#include <boost/fusion/include/transform.hpp>

#include "nodes/fb5.h"
#include "nodes/torus.h"
#include "nodes/discrete.h"
#include "nodes/gaussian.h"
#include "nodes/vonmises.h"
#include "nodes/hidden.h"

#include "default_parameters_torus_dbn.h"
#include "default_parameters_fb5_dbn.h"
#include "default_parameters_torus_cs_dbn.h"

#include "protein/definitions.h"
#include "protein/protein_data.h"

namespace phaistos {

//! Main class for models modeling backbone degrees of freedom (e.g. TorusDBN, Fb5Hmm)
template <typename DERIVED_CLASS, typename NODES>
class BackboneDBN {
public:

     //! Nodes template variable for external reference
     typedef NODES Nodes;

     //! Data structure to keep track of modifications made to the model
     class ModifiedRegion {
     public:
          //! Beginning of region
          int start_index;

          //! End of region
          int end_index;

          //! Constructor
          ModifiedRegion(int start_index,
                         int end_index)
               : start_index(start_index),
                 end_index(end_index) {}

          //! Comparison operator
          friend bool operator<(const ModifiedRegion &first, const ModifiedRegion &second) {
               if (first.start_index == second.start_index)
                    return first.end_index < second.end_index;
               else
                    return first.start_index < second.start_index;
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const ModifiedRegion &mr) {
               o << "(" << mr.start_index << " " << mr.end_index << ")";
               return o;
          }                                        
     };


     //! Index types for making node selections: size of node vector
     typedef typename boost::fusion::result_of::size<NODES>::type NODE_SIZE;
     //! Index types for making node selections: all nodes
     typedef boost::mpl::range_c<int, 0, NODE_SIZE::value > ALL_NODES;
     //! Index types for making node selections: all child nodes (hidden node disabled)
     typedef boost::mpl::range_c<int, 1, NODE_SIZE::value > ALL_CHILD_NODES;

     //! Node indices (must match the order in the node definition)
     typedef boost::mpl::int_<0> HIDDEN_NODE;

     //! Make sure that Hidden node is at first position
     BOOST_MPL_ASSERT((boost::is_same<typename boost::fusion::result_of::at_c<NODES,0>::type,
                                      HiddenNode<DERIVED_CLASS>*&>));

     //! Define child node type by removing hidden node
     typedef typename boost::fusion::result_of::as_vector<typename boost::mpl::pop_front<NODES>::type >::type CHILD_NODES;


     //! Random number engine
     RandomNumberEngine *random_number_engine;

     //! Fusion vector of nodes
     NODES nodes;

     //! Length of chain
     int sequence_length;                     

     //! number of hidden node states
     int hidden_node_size;

     //! Forward variables
     double **fw;

     //! Backward variables
     double **bw;

     //! Vector of regions that have been modified outside the DBN
     //! and for which the DBN needs to be reinitialized
     std::vector<ModifiedRegion>  inconsistent_regions;

     
     ///////////////////////////////////////////////////////////////////////////
     //// Internal duplication support (for example for threading purposes) ////
     ///////////////////////////////////////////////////////////////////////////

     //! Vector of copies of self.
     std::vector<DERIVED_CLASS*> copies;

     //! Index specifying where in the copies vector "this" is located.
     int this_index;


     //! Copy constructor helper function
     //! \param other Source object.
     void copy(const BackboneDBN &other) {
          // Copy forward and backward matrices
          fw = new double*[sequence_length];
          bw = new double*[sequence_length];     
          for (int i=0; i<sequence_length; i++) {
               fw[i] = new double[hidden_node_size];
               bw[i] = new double[hidden_node_size];
          }     
          for (int i=0; i<this->sequence_length; i++) {
               for (int j=0; j<hidden_node_size; j++) {
                    this->fw[i][j] = other.fw[i][j];
               }
          }

          for (int i=0; i<this->sequence_length; i++) {
               for (int j=0; j<hidden_node_size; j++) {
                    this->bw[i][j] = other.bw[i][j];
               }
          }
     }

     //! Internal functor for copying nodes (used by copy())
     struct _copy_node {

          //! Define result type
          template<class> struct result;

          //! Define result type
          template<class F, class T>
          struct result<F(T)> {
               //! Define result type
               typedef typename boost::remove_reference<T>::type type;
          };

          //! Internal copy of model
          DERIVED_CLASS *dbn;

          //! Constructor
          _copy_node(DERIVED_CLASS *dbn):dbn(dbn){}

          //! Overload () operator. Makes a copy of a node, passing along 
          //! the dbn as an argument.
          //! \param node DBN node of which copy should be made
          //! \return New node.
          template <typename TYPE>
          TYPE operator()(TYPE node) const {
               return new typename boost::remove_pointer<TYPE>::type(*node, dbn);
          }
     };

     //! Initialize nodes
     //! \param complete_reinitialization Whether sequence data in node should be reinitialized.
     void init_nodes(bool complete_reinitialization=true) {
          boost::fusion::for_each(nodes, _init_node(sequence_length, complete_reinitialization));
          this->hidden_node_size = boost::fusion::at<HIDDEN_NODE>(nodes)->size;
     }

     //! Internal functor used for initializing nodes (used by init_nodes())
     struct _init_node {
          //! Internal copy of sequence length
          int sequence_length;

          //! Whether to initialize sequences when initializing nodes
          bool complete_reinitialization;

          //! Constructor
          _init_node(int sequence_length, bool complete_reinitialization)
               : sequence_length(sequence_length), 
                 complete_reinitialization(complete_reinitialization){}

          //! Overload () operator. Initializes node.
          //! \param node DBN node to be initialized.
          template <typename TYPE>
          void operator()(TYPE *node) const {
               node->init(sequence_length, complete_reinitialization);               
          }
     };

     //! Internal functor used for deleting nodes (used by ~BackboneDBN())
     struct _delete_node {
          //! Overload () operator. Deletes node.
          //! \param node DBN node to be deleted
          template <typename TYPE>
          void operator()(TYPE *node) const {
               delete node;               
          }
     };


     //! Samples values for unobserved nodes.
     void init_sequences() {
          inconsistent_regions.clear();
          sample();
          accept();
     }


     //! Internal functor to sample from a node (used by sample())
     struct _sample {
          //! Internal copy of start index 
          int start_index;

          //! Internal copy of endindex 
          int end_index;

          //! Constructor
          _sample(int start_index, int end_index)
               : start_index(start_index),end_index(end_index){}

          //! Overload () operator. Samples values for a node.
          //! \param node DBN node to sample from.
          template <typename TYPE>
          void operator()(TYPE *node) const {
               node->sample(start_index, end_index);               
          }
     };


     //! Internal functor to get a sequence vector from a node
     struct _get_sequence_vector {
          //! Reference to vector in which result in stored
          std::vector<std::vector<double> > &target_vector;

          //! Internal copy of start index 
          int start_index;

          //! Internal copy of endindex 
          int end_index;

          //! Constructor
          _get_sequence_vector(int start_index, int end_index, 
                               std::vector<std::vector<double> > &target_vector)
               : target_vector(target_vector),
                 start_index(start_index), end_index(end_index) {}

          //! Overload () operator. Retrieves the sequence vector from the specified node.
          //! \param node DBN node from which to retrieve sequence.
          template <typename NODE_TYPE>
          inline void operator()(NODE_TYPE *node) const {

               for (typename NODE_TYPE::Iterator it(*node, start_index, end_index); !it.end(); ++it) {
                    target_vector[it.sequence_index-start_index].push_back(*it);
                    if (!is_initialized(*it)) {
                         set_initial_value(target_vector[it.sequence_index-start_index].back());
                    }
               }
          }
     };


     //! Internal functor to set a sequence in a node
     template <typename NODE_OUTPUT_TYPE>
     struct _set_sequence_vector {
          //! Internal reference to source vector
          const std::vector<std::vector<NODE_OUTPUT_TYPE> > &source_vector;

          //! Internal copy of start index 
          int start_index;

          //! Internal copy of endindex 
          int end_index;

          //! The index in the source vector at which data should be read
          int &dimension_offset;

          //! Whether the emission status of the node should be set to fixed after sequence is set.
          bool fix_emission;

          //! Whether to set the observed values based on the sequence values 
          bool set_observed;

          //! Constructor
          _set_sequence_vector(const std::vector<std::vector<NODE_OUTPUT_TYPE> > &source_vector,
                               int start_index, 
                               int &dimension_offset, 
                               bool fix_emission=true,
                               bool set_observed=true)
               : source_vector(source_vector),
                 start_index(start_index), end_index(start_index+source_vector.size()),
                 dimension_offset(dimension_offset),
                 fix_emission(fix_emission),
                 set_observed(set_observed){}

          
          //! Overload () operator. Sets the sequence for the specified node
          //! \param node DBN node.
          template <typename NODE_TYPE>
          inline void operator()(NODE_TYPE *node) const {

               int sequence_index=-1;
               for (typename NODE_TYPE::Iterator it(*node,start_index,end_index); !it.end(); ++it) {
                    if (((it.sequence_index-start_index) < (int)source_vector.size()) &&
                        (dimension_offset+it.dimension_index) < (int)source_vector[it.sequence_index-start_index].size()) 
                         *it = source_vector[it.sequence_index-start_index][dimension_offset+it.dimension_index];

                    // Check if we are starting with a new sequence index
                    if (sequence_index != it.sequence_index) {
                         if (set_observed)
                              it.observed_entry() = true;
                    }
                    sequence_index = it.sequence_index;

                    // if current value is uninitialized, we set the corresponding observed entry
                    if (!is_initialized(*it)) {
                         if (set_observed)
                              it.observed_entry() = false;
                    }
               }
               if (fix_emission)
                    node->fixed = true;

               dimension_offset += node->size;
          }
     };

     //! Internal functor to set the emission state of a node. See Node::fixed
     struct _add_offset {
          //! Internal copy of value to be set
          double offset;

          //! Constructor
          _add_offset(double offset):offset(offset){}

          //! Overload () operator. Sets the emission state of the specified node.
          //! \param node DBN node.          
          template <typename TYPE>
          inline void operator()(TYPE *node) const {
               for (int i=0; i<node->sequence_length; ++i) {
                    bool set_backup = true;
                    node->set(i, node->get_sequence()[i]+offset, set_backup);
               }
          }
     };

     //! Internal functor to set the emission state of a node. See Node::fixed
     struct _set_emission_state {
          //! Internal copy of value to be set
          bool value;

          //! Constructor
          _set_emission_state(bool value):value(value){}

          //! Overload () operator. Sets the emission state of the specified node.
          //! \param node DBN node.          
          template <typename TYPE>
          inline void operator()(TYPE *node) const {
               node->fixed = value;
          }
     };

     //! Internal functor to accept the new state of a node.
     struct _accept {
          //! Sequence index at which to start
          int start_index;

          //! Sequence index at which to end
          int end_index;

          //! Constructor
          _accept(int start_index=0, int end_index=-1)
               : start_index(start_index), end_index(end_index) {}

          //! Overload () operator. Calls accept in specified node.         
          //! \param node DBN node.          
          template <typename TYPE>
          void operator()(TYPE *node) const {
               node->accept(start_index, end_index);
          }
     };

     //! Internal functor to reject the new state of a node.
     struct _reject {
          //! Sequence index at which to start
          int start_index;

          //! Sequence index at which to end
          int end_index;

          //! Constructor
          _reject(int start_index=0, int end_index=-1)
               : start_index(start_index), end_index(end_index) {}

          //! Overload () operator. Calls reject in specified node.         
          //! \param node DBN node.          
          template <typename TYPE>
          void operator()(TYPE *node) const {
               node->reject(start_index, end_index);
          }
     };


     //! Functor for outputting node (used by operator<<())
     struct _output {
          //! Interal reference to ostream object
          std::ostream &o;

          //! Constructor
          _output(std::ostream &o):o(o){}

          //! Overload () operator. Calls << on the specified node.
          //! \param node DBN node.                    
          template <typename TYPE>
          inline void operator()(TYPE *node) const {
               o << *node << "\n";
          }
     };



     //! Viterbi algorithm (internal).
     //!
     //! \param start Start index
     //! \param end End index
     //!
     //! \return Vector with most probable state for each position.
     std::vector<int> _viterbi(int start, int end) {

          int length = (end-start);

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);

          // Get transition matrix for fast access
          // The alternative: hiddenNode->transition(i,j) is significantly slower
          double **cpd_transition = hidden_node->get_transition_matrix();

          int **backtraceMatrix = new int*[sequence_length];
          for (int i=0; i<sequence_length; i++) {
               backtraceMatrix[i] = new int[hidden_node_size];
          }

          if(start==0) {
               // Resample h at position 0        
               double sum = 0;
               for (int h=0; h<hidden_node_size; h++) {

                    fw[0][h] = hidden_node->start(h)*hidden_node->emission(h,0);
                    sum += fw[0][h];
               }

               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h] = fw[0][h]/sum;
               }
          } else {
          
               // Fill FW at position start+1        
               double sum = 0;
               int start_state = hidden_node->sequence[start-1];
               int l = start;

               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition[start_state][h];
                    fw[0][h] = tr*hidden_node->emission(h,l);
                    sum += fw[0][h];
               }
               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h]=fw[0][h]/sum;
               }
          }

          for(int l=1; l<length; l++) {
               double sum = 0;

               for (int h=0; h<hidden_node_size; h++) {

                    int maxInt = -1;
                    double max = 0.0;
                    for(int ph=0; ph<hidden_node_size; ph++) {
                         double val = fw[l-1][ph]*cpd_transition[ph][h];
                         if ((maxInt == -1) or (val > max)) {
                              max = val;
                              maxInt = ph;
                         }
                    }

                    fw[l][h] = max*hidden_node->emission(h, start+l);
                    sum+=fw[l][h];
                    backtraceMatrix[l][h] = maxInt;
               }

               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[l][h]=fw[l][h]/sum;
               }
          }

          if (end != sequence_length) {
               int end_state = hidden_node->sequence[end];
               double sum = 0.0;
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition[h][end_state];
                    fw[length-1][h] *= tr;
                    sum += fw[length-1][h];
               }
               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[length-1][h]=fw[length-1][h]/sum;
               }
          }

          // Find maximum h value at last index;
          int maxInt = -1;
          double max=0.0;
          for (int h=0; h<hidden_node_size; h++) {
               double val = fw[length-1][h];
               if ((maxInt == -1) || (val > max)) {
                    max = val;
                    maxInt = h;
               }
          }

          std::vector<int> path(sequence_length, -1);
          path[length-1] = maxInt;
          for (int i=length-2; i>=0; i--) {
               int prevH = path[i+1];
               int currentH = backtraceMatrix[i+1][prevH];
               path[i] = currentH;
          }
     
          for (int i=0; i<sequence_length; i++) {
               delete[] backtraceMatrix[i];
          }
          delete[] backtraceMatrix;

          return path;
     }


     //! Viterbi algorithm in log space (internal).
     //!
     //! \param start Start index
     //! \param end End index
     //!
     //! \return Vector with most probable state for each position.
     std::vector<int> _viterbi_log(int start, int end) {

          int length = (end-start);

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);

          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition_log = hidden_node->get_transition_log_matrix();


          int **backtraceMatrix = new int*[sequence_length];
          for (int i=0; i<sequence_length; i++) {
               backtraceMatrix[i] = new int[hidden_node_size];
          }
     
          if(start==0) {
               // Resample h at position 0
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h] = hidden_node->start_log(h) + hidden_node->emission_log(h,0);
               }
          } else {          
               // Fill FW at position start+1
               int start_state = hidden_node->sequence[start-1];
               int l = start;

               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition_log[start_state][h];
                    fw[0][h] = tr + hidden_node->emission_log(h,l);
               }
          }

          for(int l=1; l<length; l++) {
               for (int h=0; h<hidden_node_size; h++) {

                    int maxInt = -1;
                    double max=0.0;
                    for(int ph=0; ph<hidden_node_size; ph++) {
                         double val = fw[l-1][ph] + cpd_transition_log[ph][h];
                         if ((maxInt == -1) or (val > max)) {
                              max = val;
                              maxInt = ph;
                         }
                    }

                    fw[l][h] = max + hidden_node->emission_log(h,start+l);
                    backtraceMatrix[l][h] = maxInt;
               }
          }

          if (end != sequence_length) {
               int end_state = hidden_node->sequence[end];
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition_log[h][end_state];
                    fw[length-1][h] += tr;
               }
          }

          // Find maximum h value at last index;
          int maxInt = -1;
          double max=0.0;
          for (int h=0; h<hidden_node_size; h++) {
               double val = fw[length-1][h];
               if ((maxInt == -1) or (val > max)) {
                    max = val;
                    maxInt = h;
               }
          }

          std::vector<int> path(sequence_length, -1);
          path[length-1] = maxInt;
          for (int i=length-2; i>=0; i--) {
               int prevH = path[i+1];
               int currentH = backtraceMatrix[i+1][prevH];
               path[i] = currentH;
          }
     
          for (int i=0; i<sequence_length; i++) {
               delete[] backtraceMatrix[i];
          }
          delete[] backtraceMatrix;

          return path;
     }


     //! Posterior decoding.
     //!
     //! \param start Start index
     //! \param end End index
     //!
     //! \return Vector of discrete probability distributions over states for each position.
     std::vector<std::vector<double> > _posterior(int start, int end) {

          double ll1, ll2;
          bool include_end = false;

          // Calculate fw and bw tables
          ll1 = forward(start, end);
          ll2 = backward(start, end, include_end);
          assert(check_equality(ll1, ll2));

          std::vector<std::vector<double> > posterior_hidden;
     
          for (int i=start; i<end; i++) {
               std::vector<double> position_vec;
               double sum=0.0;
               for (int j=0; j<hidden_node_size; j++) {
                    double val = std::exp(std::log(fw[i][j])+std::log(bw[i][j]));
                    position_vec.push_back(val);
                    sum += val;
               }
               for (int j=0; j<hidden_node_size; j++) {
                    position_vec[j]/=sum;

               }
               posterior_hidden.push_back(position_vec);
          }
          return posterior_hidden;
     }

     
     //! Posterior decoding in log space.
     //!
     //! \param start Start index
     //! \param end End index
     //!
     //! \return Vector of discrete probability distributions over states for each position.
     std::vector<std::vector<double> > _posterior_log(int start, int end) {

          double ll1, ll2;
          bool include_end = false;

          // Calculate fw and bw tables
          ll1 = forward_log(start, end);
          ll2 = backward_log(start, end, include_end);
          assert(check_equality(ll1, ll2));

          std::vector<std::vector<double> > posterior_hidden;

          for (int i=start; i<end; i++) {
               std::vector<double> position_vec;
               for (int j=0; j<hidden_node_size; j++) {
                    position_vec.push_back(std::exp((fw[i][j]+bw[i][j])-ll1));
               }
               posterior_hidden.push_back(position_vec);
          }
          return posterior_hidden;     
     }




public:

     //! Local settings class.
     const class Settings: public ::Settings {
     public:

          //! Whether to do calculation in log space
          //! This is somewhat slower, but can in some rare cases be 
          //! useful to avoid overflow in the fw and bw matrices.
          bool log_space;		

          //! Distribution state of start distribution (N-terminus)
          TransitionEmissionState start_distribution;

          //! Distribution state of transition distribution
          TransitionEmissionState transition_distribution;

          //! Optional initialization data 
          //! input from PDB file
          std::string initial_pdb_file;

          //! Which polypeptide to use in pdb file
          int initial_pdb_file_index;

          //! Sequence length
          int sequence_length;

          //! Level of debug information
          int debug;          
     public:

          //! Constructor. Defines default values for settings object.
          Settings(bool log_space=false,
                   TransitionEmissionState start_distribution=NORMAL,
                   TransitionEmissionState transition_distribution=NORMAL,
                   std::string initial_pdb_file="",
                   int initial_pdb_file_index=0,
                   int sequence_length=-1,
                   int debug=0)
               : log_space(log_space),
                 start_distribution(start_distribution),
                 transition_distribution(transition_distribution),
                 initial_pdb_file(initial_pdb_file),
                 initial_pdb_file_index(initial_pdb_file_index),
                 sequence_length(sequence_length),
                 debug(debug) { 

          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "log-space:" << settings.log_space << "\n";
               o << "start-distribution:" << settings.start_distribution << "\n";
               o << "transition-distribution:" << settings.transition_distribution << "\n";
               o << "initial-pdb-file:" << settings.initial_pdb_file << "\n";
               o << "initial-pdb-file-index:" << settings.initial_pdb_file_index << "\n";
               o << "sequence-length:" << settings.sequence_length << "\n";
               o << "debug:" << settings.debug << "\n";
               return o;
          }
     } settings;    //!< Local settings object
     

     //! Constructor
     //!
     //! \param hidden_node Pointer to the hidden node object
     //! \param child_nodes Fusion vector of child nodes
     //! \param settings Local Settings object.
     //! \param random_number_engine Object from which random number generators can be created.
     BackboneDBN(HiddenNode<DERIVED_CLASS>* hidden_node,
                 const CHILD_NODES &child_nodes,
                 const Settings &settings=Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : random_number_engine(random_number_engine), 
            nodes(boost::fusion::insert(child_nodes, 
                                        boost::fusion::begin(child_nodes), 
                                        hidden_node)),
            sequence_length(0),
            hidden_node_size(0),
            fw(NULL), bw(NULL),
            inconsistent_regions(std::vector<ModifiedRegion>(1, ModifiedRegion(-1,-1))),
            this_index(0),
            settings(settings) {
     }


     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     BackboneDBN(const BackboneDBN &other)
          : random_number_engine(other.random_number_engine),
            nodes(boost::fusion::transform(other.nodes, _copy_node((DERIVED_CLASS*)this))),
            sequence_length(other.sequence_length),
            hidden_node_size(other.hidden_node_size),
            inconsistent_regions(other.inconsistent_regions),
            this_index(0),
            settings(other.settings) {

          // When copying, do not do a complete reinitializtion
          bool complete_reinitialization=false;
          init_nodes(complete_reinitialization);

          copy(other);
     }

     //! Copy constructor. Different random_number_generator specified.
     //!
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be created.
     BackboneDBN(const BackboneDBN &other,
                 RandomNumberEngine *random_number_engine)
          : random_number_engine(random_number_engine),
            nodes(boost::fusion::transform(other.nodes, _copy_node((DERIVED_CLASS*)this))),
            sequence_length(other.sequence_length),
            hidden_node_size(other.hidden_node_size),
            inconsistent_regions(other.inconsistent_regions),
            this_index(0),
            settings(other.settings) {

          // When copying, do not do a complete reinitializtion
          bool complete_reinitialization=false;
          init_nodes(complete_reinitialization);

          copy(other);
     }     

     //! Destructor
     ~BackboneDBN() {

          // Delete forward and backward matrices
          if (fw) {
               for (int i=0; i<sequence_length; i++) {
                    delete[] fw[i];
               }
               delete[] fw;
          }
          if (bw) {
               for (int i=0; i<sequence_length; i++) {
                    delete[] bw[i];
               }
               delete[] bw;
          }

          // Delete nodes
          boost::fusion::for_each(nodes, _delete_node());

          // Clean up threads
          if (this_index==0) {
               for (unsigned int i=1; i<copies.size(); i++) {
                    delete copies[i];
               }          
          }
     }

     //! Initialization
     //!
     //! \param sequence_length Number of slices in DBN.
     void init(int sequence_length){
          
          this->sequence_length = sequence_length;

          // Initialize nodes
          init_nodes();

          // The DBN is not initialized
          inconsistent_regions = std::vector<ModifiedRegion>(1, ModifiedRegion(-1,-1));

          // Forward and backward variables arrays
          fw = new double*[sequence_length];
          bw = new double*[sequence_length];     
          for (int i=0; i<sequence_length; i++) {
               fw[i] = new double[hidden_node_size];
               bw[i] = new double[hidden_node_size];
          }    
     }

     //! Accept proposed change
     void accept(int start_index=0, int end_index=-1) {
          if (end_index<0)
               end_index = this->sequence_length;

          boost::fusion::for_each(nodes, _accept(start_index, end_index));
     }


     //! Reject proposed change
     void reject(int start_index=0, int end_index=-1) {
          if (end_index<0)
               end_index = this->sequence_length;

          boost::fusion::for_each(nodes, _reject(start_index, end_index));
     }


     //! Forward algorith.
     //!
     //! Fills the fw 2D-array that is subsequently used by backtrack().
     //!
     //! \param start Start index
     //! \param end End index
     //! \param use_backup Whether to use backup values (values from before last update).
     //!
     //! \return log-likelihood value
     double forward(int start, int end, bool use_backup=false) {

          int length = (end-start);

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);

          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition = hidden_node->get_transition_matrix();

          double ll = 0.0;
          if(start==0) {
               // Resample h at position 0        
               double sum = 0;
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h] = hidden_node->start(h)*hidden_node->emission(h, 0, use_backup);
                    sum += fw[0][h];
               }

               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h] = fw[0][h]/sum;
               }
               ll += std::log(sum);

          } else {
          
               // Fill FW at position start+1        
               double sum = 0;
               int l = start;

               int start_state = hidden_node->sequence[start-1];
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h] = cpd_transition[start_state][h]*hidden_node->emission(h, l, use_backup);
                    sum += fw[0][h];
               }

               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h]=fw[0][h]/sum;
               }
               ll += std::log(sum);
          }

          for(int l=1; l<length; l++) {
               double sum = 0;

               for (int h=0; h<hidden_node_size; h++) {
                    double s = 0;

                    for(int ph=0; ph<hidden_node_size; ph++) {
                         s += fw[l-1][ph]*cpd_transition[ph][h];
                    }

                    fw[l][h] = s*hidden_node->emission(h, start+l, use_backup);
                    sum+=fw[l][h];
               }

               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[l][h]=fw[l][h]/sum;
               }
               ll += std::log(sum);          
          }

          if (end != sequence_length) {
               int end_state = hidden_node->sequence[end];
               double sum = 0.0;
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition[h][end_state];
                    fw[length-1][h] *= tr;
                    sum += fw[length-1][h];
               }
               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    fw[length-1][h]=fw[length-1][h]/sum;
               }
               ll += std::log(sum);          
          }

          return ll;
     }


     //! Forward algorith in log-space.
     //! 
     //! Fills the fw 2D-array that is subsequently used by backtrack_log().
     //!
     //! \param start Start index
     //! \param end End index
     //! \param use_backup Whether to use backup values (values from before last update).
     //!
     //! \return log-likelihood value
     double forward_log(int start, int end, bool use_backup=false) {

          int length = (end-start);

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);

          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition_log = hidden_node->get_transition_log_matrix();

          if(start==0) {
               // Resample h at position 0
               for (int h=0; h<hidden_node_size; h++) {
                    fw[0][h] = hidden_node->start_log(h) + hidden_node->emission_log(h, 0, use_backup);
               }
          } else {          
               // Fill FW at position start+1
               int start_state = hidden_node->sequence[start-1];
               int l = start;

               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition_log[start_state][h];
                    fw[0][h] = tr + hidden_node->emission_log(h, l, use_backup);
               }
          }

          for(int l=1; l<length; l++) {
               for (int h=0; h<hidden_node_size; h++) {
                    double s = std::log(0.0);

                    for(int ph=0; ph<hidden_node_size; ph++) {
                         s = ln_sum(s, fw[l-1][ph] + cpd_transition_log[ph][h]);
                    }

                    fw[l][h] = s + hidden_node->emission_log(h, start+l, use_backup);
               }
          }

          if (end != sequence_length) {
               int end_state = hidden_node->sequence[end];
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition_log[h][end_state];
                    fw[length-1][h] += tr;
               }
          }
     
          double ll = std::log(0.0);
          for (int h=0; h<hidden_node_size; h++) {
               ll = ln_sum(ll, fw[length-1][h]);
          }

          return ll;
     }



     //! Viterbi algorithm - Find optimal hidden node path
     //!
     //! \param start Start index
     //! \param end End index
     //!
     //! \return Vector with most probable state for each position.
     std::vector<int> viterbi(int start=0, int end=-1) {

          if(end < 0 || end > sequence_length)
               end = sequence_length;
     
          if (settings.log_space) {
               return _viterbi_log(start, end);
          } else{
               return _viterbi(start, end);
          }
     }


     //! Backward algorith.
     //!
     //! Fills the bw 2D-array and calculates the log-likelihood.
     //!
     //! \param start Start index
     //! \param end End index
     //! \param include_end Whether to include transition at end
     //!
     //! \return log-likelihood value
     double backward(int start, int end, bool include_end=true) {

          int length = (end-start);

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);
     
          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition = hidden_node->get_transition_matrix();
     
          double ll = 0.0;
     
          if (end == sequence_length) {
               double sum = 0.0;
               for (int h=0; h<hidden_node_size; h++) {
                    // Uniform initialization
                    bw[length-1][h] = 1.0;
                    sum += bw[length-1][h];
               }
	  
               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    bw[length-1][h] = bw[length-1][h]/sum;
               }
               ll += std::log(sum);
          } else {
               int end_state = hidden_node->sequence[end];
               double sum = 0.0;
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition[h][end_state];
                    bw[length-1][h] = tr;
                    sum += bw[length-1][h];
               }
	  
               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    bw[length-1][h] = bw[length-1][h]/sum;
               }
               ll += std::log(sum);
          }

          for(int l=length-2; l>=0; l--) {
               double sum = 0;

               // Calculate emissions once (cache)
               std::vector<double> emissions(hidden_node_size);
               for(int nh=0; nh<hidden_node_size; nh++) {
                    emissions[nh] = hidden_node->emission(nh, start+l+1);
               }

               for (int h=0; h<hidden_node_size; h++) {

                    double s = 0;
                    for(int nh=0; nh<hidden_node_size; nh++) {
                         s += cpd_transition[h][nh] * emissions[nh] * bw[l+1][nh];
                    }

                    bw[l][h] = s;
                    sum+=bw[l][h];
               }

               // Scale
               for (int h=0; h<hidden_node_size; h++) {
                    bw[l][h]=bw[l][h]/sum;
               }
               ll += std::log(sum);          
          }

          if (start==0) {
               double sum = 0.0;
               for (int h=0; h<hidden_node_size; h++) {
                    double val = hidden_node->start(h)*hidden_node->emission(h, start)*bw[0][h];
                    if (include_end)
                         bw[0][h] = val;
                    sum += val;
               }
               ll += std::log(sum);
          } else {
               int start_state = hidden_node->sequence[start-1];
               double sum = 0.0;
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition[start_state][h];
                    double val = tr*hidden_node->emission(h, start)*bw[0][h]; 
                    if (include_end)
                         bw[0][h] = val;
                    sum += val;
               }
               ll += std::log(sum);
          }
          return ll;
     }


     //! Backward algorith in log-space.
     //!
     //! Fills the bw 2D-array and calculates the log-likelihood.
     //!
     //! \param start Start index
     //! \param end End index
     //! \param include_end Whether to include transition at end
     //!
     //! \return log-likelihood value
     double backward_log(int start, int end, bool include_end=true) {

          int length = (end-start);

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);
     
          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition_log = hidden_node->get_transition_log_matrix();

          if (end == sequence_length) {
               for (int h=0; h<hidden_node_size; h++) {
                    // Uniform initialization
                    bw[length-1][h] = 0.0;
               }
          } else {
               int end_state = hidden_node->sequence[end];
               for (int h=0; h<hidden_node_size; h++) {
                    double tr = cpd_transition_log[h][end_state];
                    bw[length-1][h] = tr;
               }
          }

          for(int l=length-2; l>=0; l--) {
               for (int h=0; h<hidden_node_size; h++) {

                    // Calculate emissions once (cache)
                    std::vector<double> emissions_log(hidden_node_size);
                    for(int nh=0; nh<hidden_node_size; nh++) {
                         emissions_log[nh] = hidden_node->emission_log(nh, start+l+1);
                    }

                    double s = std::log(0.0);
                    for(int nh=0; nh<hidden_node_size; nh++) {
                         s = ln_sum(s, cpd_transition_log[h][nh] + emissions_log[nh] + bw[l+1][nh]);
                    }
                    bw[l][h] = s;
               }
          }

          double ll = std::log(0.0);
          if (include_end) {
               if (start==0) {
                    for (int h=0; h<hidden_node_size; h++) {
                         bw[0][h] += hidden_node->start(h) + hidden_node->emission_log(h, start);
                    }
               } else {
                    int start_state = hidden_node->sequence[start-1];
                    for (int h=0; h<hidden_node_size; h++) {
                         double tr = cpd_transition_log[start_state][h];
                         bw[0][h] += tr + hidden_node->emission_log(h, start);
                    }
               }

               for (int h=0; h<hidden_node_size; h++) {
                    ll = ln_sum(ll, bw[0][h]);
               }
          } else {
               if (start==0) {
                    for (int h=0; h<hidden_node_size; h++) {
                         ll = ln_sum(ll,
                                    hidden_node->start_log(h) + hidden_node->emission_log(h, start) + bw[0][h]);
                    }
               } else {
                    int start_state = hidden_node->sequence[start-1];
                    for (int h=0; h<hidden_node_size; h++) {
                         double tr = cpd_transition_log[start_state][h];
                         ll = ln_sum(ll,
                                    tr + hidden_node->emission_log(h, start) + bw[0][h]);
		    
                    }
               }
          }
     
          return ll;
     }

     //! Calculate likelihood by summing over hidden node sequences.
     //!
     //! Uses forward() or forward_log() to calculate the log-likelihood.
     //!
     //! \param start Start index
     //! \param end End index
     //! \param use_backup Whether to use backup values (values from before last update).
     //!
     //! \return log-likelihood value
     double get_log_likelihood(int start=-1, int end=-1, bool use_backup=false) {
          if(end < 0 || end > sequence_length)
               end = sequence_length;
     
          if(start < 0 || start > sequence_length)
               start = 0;

          if (settings.log_space)
               return forward_log(start, end, use_backup);
          else
               return forward(start, end, use_backup);
     }


     //! Calculate likelihood given the current hidden node sequence.
     //!
     //! \param start Start index
     //! \param end End index
     //! \param use_backup Whether to use backup values (values from before last update).
     //!
     //! \return log-likelihood value
     double get_log_likelihood_conditional(int start=-1, int end=-1, bool use_backup=false) {

          if(end < 0 || end > sequence_length)
               end = sequence_length;
     
          if(start < 0 || start > sequence_length)
               start = 0;

          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);
     
          // Get current hidden node sequence - initializes sequences if necessary
          int *hidden_sequence = hidden_node->get_sequence();
     
          if (use_backup) {
               // Use backup sequence
               hidden_sequence = hidden_node->sequence_backup;
          }
     
          double ll = 0.0;
          for (int i=start; i<end; i++) {
               ll += hidden_node->emission_log(hidden_sequence[i], i, use_backup);
          }
          return ll;
     }
     

     //! Posterior decoding.
     //!
     //! \param start Start index
     //! \param end End index
     //!
     //! \return Vector of discrete probability distributions over states for each position.
     std::vector<std::vector<double> > posterior(int start=-1, int end=-1) {

          if(end < 0 || end > sequence_length)
               end = sequence_length;
     
          if(start < 0 || start > sequence_length)
               start = 0;


          if (settings.log_space) {
               return _posterior_log(start, end);
          } else {
               return _posterior(start, end);
          }     
     }


     //! Backtrack algorithm. Samples values for a subsequence given the current forward matrix.
     //!
     //! Note: FB acts BETWEEN start and end.
     //!
     //! \param start Start index
     //! \param end End index
     //! \param sampling Whether sampling should be done, or whether we should just return the likelihood calculation
     //!
     //! \return log-likelihood value
     double backtrack(int start, int end, bool sampling=true) {
     
          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);

          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition = hidden_node->get_transition_matrix();
          
          int length = end-start;

          // Keep track of log-likelihood of sample
          double ll = 0.0;
     
          int end_state = 0;     
          double normalization = 0.0;
          if (sampling) {
               // Pick last state
               end_state = DiscreteSampler::sample(hidden_node_size, fw[length-1], &normalization, 
                                                   this->random_number_engine);

               // Set new end state in proposed seq
               hidden_node->sequence[end-1] = end_state;
          } else {
               for (int i=0; i<hidden_node_size; i++) {
                    normalization += fw[length-1][i];
               }
          }

          int h = hidden_node->sequence[end-1];     

          ll += std::log(fw[length-1][h]/normalization);
          
          for(int l=length-2; l>=0; l--) {

               double *p = new double[hidden_node_size];

               for(int ph=0; ph<hidden_node_size; ph++) {
                    p[ph] = cpd_transition[ph][h]*fw[l][ph];
               }

               double normalization = 0.0;
               if (sampling) {
                    for (int i=0; i<hidden_node_size; i++) {
                         normalization += p[i];
                    }

                    if (!std::isfinite(normalization)) {
                         std::cout << *this << "\n";
                         std::cout << "\n\nERROR: Nan value used by sampler - aborting\n";
                         assert(false);
                    }

                    // Update next hidden node value for next backtrack step
                    h = DiscreteSampler::sample(hidden_node_size, p, &normalization,
                                                this->random_number_engine);

                    // Assign new hidden node to temporary hidden node sequence
                    hidden_node->sequence[start+l] = h;
               } else {
                    for (int i=0; i<hidden_node_size; i++) {
                         normalization += p[i];
                    }
                    h = hidden_node->sequence[start+l];
               }
          
               ll += std::log(fw[l][hidden_node->sequence[start+l]]/normalization);

               delete[] p;
          }

          return ll;
     }


     //! Backtrack algorithm in log space. Samples values for a subsequence given the current forward matrix.
     //!
     //! Note: FB acts BETWEEN start and end.
     //!
     //! \param start Start index
     //! \param end End index
     //! \param sampling Whether sampling should be done, or whether we should just return the likelihood calculation
     //!
     //! \return log-likelihood value
     double backtrack_log(int start, int end, bool sampling) {
     
          // Extract hidden node from node vector
          HiddenNode<DERIVED_CLASS> *hidden_node = boost::fusion::at<HIDDEN_NODE>(nodes);

          // Get transition matrix for fast access
          // The alternative: hidden_node->transition(i,j) is significantly slower
          double **cpd_transition_log = hidden_node->get_transition_log_matrix();
          
          int length = end-start;

          // Keep track of log-likelihood of sample
          double ll = 0.0;
     
          int end_state = 0;     
          double normalization = 0;
          if (sampling) {     
               // Pick last state
               end_state = DiscreteSampler::sample_log(hidden_node_size, fw[length-1], &normalization,
                                                       this->random_number_engine);

               // Set new end state in proposed seq
               hidden_node->sequence[end-1] = end_state;
          } else {
               for (int i=0; i<hidden_node_size; i++) {
                    normalization += std::exp(fw[length-1][i]);
               }
          }

          int h = hidden_node->sequence[end-1];     

          ll += fw[length-1][h] - std::log(normalization);
               
          for(int l=length-2; l>=0; l--) {
               double *p = new double[hidden_node_size];

               bool allInf = true;
               for(int ph=0; ph<hidden_node_size; ph++) {
                    p[ph] = cpd_transition_log[ph][h] + fw[l][ph];
                    if (std::isinf(p[ph])==0)
                         allInf = false;
               }
               assert(!allInf);

               double normalization;
               if (sampling) {
                    // Update next hidden node value for next backtrack step
                    h = DiscreteSampler::sample_log(hidden_node_size, p, &normalization,
                                                    this->random_number_engine);

                    // Assign new hidden node to temporary hidden node sequence
                    hidden_node->sequence[start+l] = h;
               } else {
                    for (int i=0; i<hidden_node_size; i++) {
                         normalization += std::exp(p[i]);
                    }
                    h = hidden_node->sequence[start+l];
               }

               ll += fw[l][hidden_node->sequence[start+l]] - std::log(normalization);
                    
               delete[] p;
          }
     
          return ll;
     }


     //! Sample hidden node sequence and unobserved nodes.
     //!
     //! \param start Start index
     //! \param end End index
     void sample(int start=-1, int end=-1) {

          if (sequence_length==0)
               return;

          if(end < 0 || end > sequence_length)
               end = sequence_length;
     
          if(start < 0 || start > sequence_length)
               start = 0;
     
          // Make sure that the complete hidden node sequence has previously been sampled
          if (inconsistent_regions.size() > 0) {

               // Reinitialize sequences
               init_sequences();

               // Everything has been resampled. No need to continue.
               return;
          }

          // Resample nodes
          boost::fusion::for_each(nodes, _sample(start, end));
     }


     //! Retrieve sequence for specified nodes.
     //!
     //! \tparam SELECTION_NODE_INDICES boost::mpl::list of node indices (boost::mpl::int_). For instance ALL_ANGLE_NODES.
     //! \param start_index Start index
     //! \param end_index End index
     //!
     //! return A vector of values for each position.
     template <typename SELECTED_NODE_INDICES>
     std::vector<std::vector<double> > get_sequence_vector(int start_index=0, 
                                                           int end_index=-1) {

	  if(end_index < 0 || end_index > sequence_length)
	       end_index = sequence_length;
	  
          using namespace boost::fusion;
          typedef typename boost::mpl::if_<boost::fusion::traits::is_sequence<SELECTED_NODE_INDICES>, SELECTED_NODE_INDICES, boost::mpl::vector<SELECTED_NODE_INDICES> >::type SELECTED_NODE_INDICES_VECTOR;

          typedef nview<NODES,SELECTED_NODE_INDICES_VECTOR> SELECTED_NODES;

          std::vector<std::vector<double> > result(end_index-start_index, std::vector<double>());

          SELECTED_NODES selected_nodes(nodes);
          for_each(selected_nodes,
                   _get_sequence_vector(start_index, end_index, result));
          
          return result;
     }


     //! Set sequence for specified nodes - from string.
     //!
     //! \tparam SELECTION_NODE_INDICES boost::mpl::list of node indices (boost::mpl::int_). For instance ALL_ANGLE_NODES.
     //! \param input_string Input sequence.
     //! \param start_index Start index
     //! \param fix_emission Whether the emission status of the node should be set to fixed after sequence is set.
     template <typename SELECTED_NODE_INDICES>
     void set_sequence_vector(std::string input_string, int start_index=0, bool fix_emission=true) {
          using namespace boost::fusion;
          typedef typename boost::mpl::if_<boost::fusion::traits::is_sequence<SELECTED_NODE_INDICES>, SELECTED_NODE_INDICES, boost::mpl::vector<SELECTED_NODE_INDICES> >::type SELECTED_NODE_INDICES_VECTOR;

          typedef nview<NODES,SELECTED_NODE_INDICES_VECTOR> SELECTED_NODES;

          typedef typename boost::remove_pointer<typename boost::remove_reference<typename boost::mpl::at_c<SELECTED_NODES,0>::type >::type>::type::Type NodeOutputType;

          set_sequence_vector<SELECTED_NODE_INDICES_VECTOR>(input_string, &boost::lexical_cast<NodeOutputType,std::string>, 
                                                            start_index, fix_emission);
     }


     //! Set sequence for specified nodes - from string - using string-to-value-functor to translate strings to values.
     //!
     //! \tparam SELECTION_NODE_INDICES boost::mpl::list of node indices (boost::mpl::int_). For instance ALL_ANGLE_NODES.
     //! \param input_string Input sequence.
     //! \param string_to_value_functor Functor to translate from string to values.
     //! \param start_index Start index
     //! \param fix_emission Whether the emission status of the node should be set to fixed after sequence is set.
     //! \param set_observed Whether to set the observed values based on the sequence values 
     template <typename SELECTED_NODE_INDICES, typename STR_TO_VAL_FUNCTOR>
     void set_sequence_vector(std::string input_string, 
                              STR_TO_VAL_FUNCTOR *string_to_value_functor,
                              int start_index=0, 
                              bool fix_emission=true,
                              bool set_observed=true) {

          using namespace boost::fusion;
          typedef typename boost::mpl::if_<boost::fusion::traits::is_sequence<SELECTED_NODE_INDICES>, SELECTED_NODE_INDICES, boost::mpl::vector<SELECTED_NODE_INDICES> >::type SELECTED_NODE_INDICES_VECTOR;

          typedef nview<NODES,SELECTED_NODE_INDICES_VECTOR> SELECTED_NODES;

          typedef typename boost::remove_pointer<typename boost::remove_reference<typename boost::mpl::at_c<SELECTED_NODES,0>::type >::type>::type::Type NodeOutputType;

          if (input_string.size() == 0) {
               return;
          }

          std::vector<std::vector<NodeOutputType> > input_vector;

          // Remove white space at beginning and end
          boost::trim(input_string);

          // Check whether input string contains ',' or ' ' separators
          if (input_string.find(',') != std::string::npos || input_string.find(' ') != std::string::npos) {

               std::stringstream input_stream(input_string);

               // If there is only a single opening bracket, the input_string is parsed as 
               // a 1D vector
               if ((std::count(input_string.begin(), input_string.end(), '[')<=1) &&
                   (std::count(input_string.begin(), input_string.end(), '\n')==0)) {

                    std::vector<NodeOutputType> input_vector_tmp;
                    input_stream >> input_vector_tmp;
                    // Copy input to 2D array
                    for (unsigned int i=0; i<input_vector_tmp.size(); ++i) {
                         input_vector.push_back(std::vector<NodeOutputType>(1,input_vector_tmp[i]));
                    }

               } else {
                    std::vector<std::vector<std::string> > input_vector_strings;
                    input_stream >> input_vector_strings;
                    for (unsigned int i=0; i<input_vector_strings.size(); ++i) {
                         input_vector.push_back(std::vector<NodeOutputType>());
                         for (unsigned int j=0; j<input_vector_strings[i].size(); ++j) {
                              std::stringstream str_value(input_vector_strings[i][j]);
                              NodeOutputType value;
                              if (!(str_value >> value)) {
                                   set_initial_value(value);
                              }
                              input_vector.back().push_back(value);
                         }
                    }
               }

          // Otherwise parse string as vector of length-1 tokens 
          } else {
               for (unsigned int i=0; i<input_string.length(); i++) {
                    char ch = input_string[i];
                    if (ch==',' || ch==' ' || ch=='[' || ch==']' || ch=='(' || ch==')') {
                         continue;
                    } else {
                         NodeOutputType value;
                         value = (*string_to_value_functor)(std::string(1,ch));
                         input_vector.push_back(std::vector<NodeOutputType>(1,value));
                    }
               }
          }

          // Initialize DBN if this is is the first sequence that is set
          if (sequence_length == 0) {
               init(input_vector.size());
          }

          SELECTED_NODES selected_nodes(nodes);
          int dimension_offset = 0;
          for_each(selected_nodes,
                   _set_sequence_vector<NodeOutputType>(input_vector, start_index, dimension_offset, fix_emission, set_observed));
     }


     //! Set sequence for specified nodes - from 2D vector of values.
     //!
     //! \tparam SELECTION_NODE_INDICES boost::mpl::list of node indices (boost::mpl::int_). For instance ALL_ANGLE_NODES.
     //! \param input Input values.
     //! \param start_index Start index
     //! \param fix_emission Whether the emission status of the node should be set to fixed after sequence is set.
     //! \param set_observed Whether to set the observed values based on the sequence values 
     template <typename SELECTED_NODE_INDICES, typename NODE_OUTPUT_TYPE>
     void set_sequence_vector(const std::vector<std::vector<NODE_OUTPUT_TYPE> > &input, 
                              int start_index=0,
                              bool fix_emission=true,
                              bool set_observed=true) {

          using namespace boost::fusion;
          typedef typename boost::mpl::if_<boost::fusion::traits::is_sequence<SELECTED_NODE_INDICES>, SELECTED_NODE_INDICES, vector<SELECTED_NODE_INDICES> >::type SELECTED_NODE_INDICES_VECTOR;

          typedef nview<NODES,SELECTED_NODE_INDICES_VECTOR> SELECTED_NODES;

          // Initialize DBN if this is is the first sequence that is set
          if (sequence_length == 0) {
               init(input.size());
          }

          SELECTED_NODES selected_nodes(nodes);
          int dimension_offset = 0;
          for_each(selected_nodes,
                   _set_sequence_vector<NODE_OUTPUT_TYPE>(input, start_index, dimension_offset, fix_emission, set_observed));
     }


     //! Set sequence for specified nodes - from vector of values.
     //!
     //! \tparam SELECTION_NODE_INDEX node indicex (boost::mpl::int_). For instance ANGLE_NODE.
     //! \param input Input values.
     //! \param start_index Start index
     //! \param fix_emission Whether the emission status of the node should be set to fixed after sequence is set.
     //! \param set_observed Whether to set the observed values based on the sequence values 
     template <typename SELECTED_NODE_INDEX, typename NODE_OUTPUT_TYPE>
     void set_sequence_vector(const std::vector<NODE_OUTPUT_TYPE> &input, 
                              int start_index=0,
                              bool fix_emission=true,
                              bool set_observed=true) {

          // Initialize DBN if this is is the first sequence that is set
          if (sequence_length == 0) {
               init(input.size());
          }

          // Extract node from node vector
          boost::fusion::at<SELECTED_NODE_INDEX>(nodes)->set_sequence(input, 
                                                                      std::vector<bool>(),
                                                                      start_index,
                                                                      fix_emission,
                                                                      set_observed);
     }


     //! Add offset to sequence correspondint to node
     //!
     //! \tparam SELECTION_NODE_INDICES boost::mpl::list of node indices (boost::mpl::int_). For instance ALL_ANGLE_NODES.
     //! \param value Value (true=fixed)
     //! \param set_uninitialized Whether node should be flagged as being uninitialized after this operation.
     template <typename SELECTED_NODE_INDICES>
     void add_offset(double value, bool set_uninitialized=true) {
          using namespace boost::fusion;
          typedef typename boost::mpl::if_<boost::fusion::traits::is_sequence<SELECTED_NODE_INDICES>, SELECTED_NODE_INDICES, boost::mpl::vector<SELECTED_NODE_INDICES> >::type SELECTED_NODE_INDICES_VECTOR;

          typedef nview<NODES,SELECTED_NODE_INDICES_VECTOR> SELECTED_NODES;

          if (set_uninitialized) {
               this->inconsistent_regions = std::vector<ModifiedRegion>(1, ModifiedRegion(-1,-1));
          }

          SELECTED_NODES selected_nodes(nodes);
          for_each(selected_nodes,
                   _add_offset(value));
     }


     //! Set the emission (observed/non-observed) state of a selection of nodes. See Node::fixed
     //!
     //! \tparam SELECTION_NODE_INDICES boost::mpl::list of node indices (boost::mpl::int_). For instance ALL_ANGLE_NODES.
     //! \param value Value (true=fixed)
     //! \param set_uninitialized Whether node should be flagged as being uninitialized after this operation.
     template <typename SELECTED_NODE_INDICES>
     void set_emission_state(bool value, bool set_uninitialized=true) {
          using namespace boost::fusion;
          typedef typename boost::mpl::if_<boost::fusion::traits::is_sequence<SELECTED_NODE_INDICES>, SELECTED_NODE_INDICES, boost::mpl::vector<SELECTED_NODE_INDICES> >::type SELECTED_NODE_INDICES_VECTOR;

          typedef nview<NODES,SELECTED_NODE_INDICES_VECTOR> SELECTED_NODES;

          if (set_uninitialized) {
               this->inconsistent_regions = std::vector<ModifiedRegion>(1, ModifiedRegion(-1,-1));
          }

          SELECTED_NODES selected_nodes(nodes);
          for_each(selected_nodes,
                   _set_emission_state(value));
     }


     //! Return a specified node.
     //!
     //! \tparam NODE_INDEX node index (boost::mpl::int_). For instance ANGLE_NODE.
     //! return Pointer to node
     template <typename NODE_INDEX>
     typename boost::fusion::result_of::at<NODES, NODE_INDEX>::type get_node() {
          return boost::fusion::at<NODE_INDEX>(nodes);
     }


     //! Duplication - makes internal copies of self to support multithreading.
     //!
     //! \param ncopies Number of copies
     //! \param random_number_engines a vector of random number engines from which random number generators can be constructed.
     void duplicate(int ncopies, 
                    const std::vector<RandomNumberEngine *> &random_number_engines=std::vector<RandomNumberEngine *>()) {

          // Resize threads vector
          copies.resize(ncopies);

          // Cast this pointer to derived class (avoid virtual functions for efficiency)
          DERIVED_CLASS *this_pointer = (DERIVED_CLASS *)this;

          // Use the passed dbn as first thread
          copies[this_index] = this_pointer;

          if (this_index==0) {

               // Make nthreads-1 copies
               for (int i=1; i<ncopies; i++) {

                    if (random_number_engines.size() > 0) {
                         copies[i] = new DERIVED_CLASS(*this_pointer, random_number_engines[i]);
                    } else {
                         copies[i] = new DERIVED_CLASS(*this_pointer, this->random_number_engine);
                    }
                    copies[i]->this_index = i;
                    copies[i]->duplicate(ncopies);
               }

               // Make vectors in all copies consistent
               for (int i=1; i<ncopies; i++) {
                    for (int j=0; j<ncopies; j++) {
                         copies[i]->copies[j] = copies[j];
                    }
               }          

               // for (int i=0; i<ncopies; i++) {
               //      std::cout << i << ": \n";
               //      std::cout << *copies[i] << "\n";
               // }
          }     
     }

     //! Return nth copy of model
     //!
     //! \param index Index of copy
     //! \return nth duplicate
     DERIVED_CLASS &get_copy(unsigned int index) {
          if (copies.size() > index) {
               return *copies[index];
          } else if (index==0) {
               return *(DERIVED_CLASS *)this;
          } else {
               assert(false);

               // ONLY INCLUDED TO AVOID WARNING IN COMPILATION
               return *(DERIVED_CLASS *)this;
          }
     }


     //! Register a region of the sequence in the DBN as inconsistent.
     //!
     //! \param start_index Start index
     //! \param end_index End index
     void register_inconsistency(int start_index, int end_index) {
          if (inconsistent_regions.size()==0 ||
              (!(start_index == inconsistent_regions.back().start_index &&
                 end_index == inconsistent_regions.back().end_index))) {
               inconsistent_regions.push_back(ModifiedRegion(start_index, end_index));
          }
     }


     //! Bring the DBN back to a consistent state.
     //!
     //! It is not enough to merely resample the hidden node sequence
     //! for the region that has changed. The hidden node values outside
     //! this region will also be affected. The window_size option
     //! makes it possible to specify a window size to be considered. The
     //! default (-1) will resample the entire hidden node sequence.
     //!
     //! \param window_size Size of window around modified region to consider.
     void enforce_consistency(int window_size=-1) {

          // Do nothing if no regions are inconsistent
          if (inconsistent_regions.size() == 0) {
               return;
          }

          // If window size is negative, simply resample entire sequence
          if (window_size < 0) {
               this->init_sequences();
               return;
          }

          // Cast this pointer to derived class (avoid virtual functions for efficiency)
          DERIVED_CLASS *thisPointer = (DERIVED_CLASS *)this;

          // Range vector (modified regions extended with a window)
          std::vector<std::pair<int,int> > ranges;

          // Sort regions
          std::sort(inconsistent_regions.begin(), inconsistent_regions.end());

          for (unsigned int i=0; i<inconsistent_regions.size(); ++i) {

               // Create range
               std::pair<int,int> range = std::make_pair(std::max(0, 
                                                                  inconsistent_regions[i].start_index - window_size),
                                                         std::min(sequence_length, 
                                                                  inconsistent_regions[i].end_index   + window_size));
               
               // Merge ranges if possible
               if (ranges.size() > 0 && (range.first < ranges.back().second)) {
                    std::pair<int,int> &pair = ranges.back();
                    pair = std::make_pair(pair.first, range.second);
               } else {
                    ranges.push_back(std::make_pair(range.first, range.second));
               }
          }

          // Clear inconsistent regions vector (this must be done before resampling below)
          inconsistent_regions.clear();

          // Resample ranges
          for (unsigned int i=0; i<ranges.size(); ++i) {
               thisPointer->sample(ranges[i].first, ranges[i].second);
          }
          thisPointer->accept();

     }

     //! Output operator.
     friend std::ostream &operator<<(std::ostream &o, BackboneDBN<DERIVED_CLASS,NODES> &dbn) {
          o << "Sequence length: " << dbn.sequence_length << "\n";
          boost::fusion::for_each(dbn.nodes,
                                  _output(o));

          // o << "Forward matrix: \n";
          // for (int i=0; i<dbn.hidden_node_size; ++i) {
          //      for (int j=0; j<dbn.sequence_length; ++j) {
          //           o << dbn.fw[j][i] << " ";
          //      }
          //      o << "\n";
          // }
          return o;
     }
};


//! Settings class for amino acid node
class SettingsNodeAa {
public:

     //! amino acid input file
     std::string initial_aa_file;

     //! amino acid input sequence
     std::string initial_aa_sequence;

     //! Constructor
     SettingsNodeAa(std::string initial_aa_file="",
                    std::string initial_aa_sequence="")
          : initial_aa_file(initial_aa_file),                 
            initial_aa_sequence(initial_aa_sequence) {}

     //! Initialize a dbn with this Settings object
     //!
     //! \tparam NODE_INDEX index of node (boost::mpl::int_).
     //! \param dbn model object
     //! \param pdb_data Optional vector with data from a pdb file
     template <typename DBN_TYPE, typename NODE_INDEX>
     void initialize(DBN_TYPE &dbn, const std::vector<int> &pdb_data=std::vector<int>()) const {
          // Set AA input
          if (initial_aa_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX>(initial_aa_sequence, &definitions::str_to_aa);               
          else if (initial_aa_file != "") 
               dbn.template set_sequence_vector<NODE_INDEX>(file_to_string(initial_aa_file), &definitions::str_to_aa);
          else if (!pdb_data.empty()) 
               dbn.template get_node<NODE_INDEX>()->set_sequence(pdb_data);
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const SettingsNodeAa &settings) {
          o << "initial-aa-file:" << settings.initial_aa_file << "\n";
          o << "initial-aa-sequence:" << settings.initial_aa_sequence << "\n";
          return o;
     }
};


//! Settings class for secondary structure node
class SettingsNodeSs {
public:
     //! Secondary structure input file
     std::string initial_ss_file;

     //! Secondary structure input sequence
     std::string initial_ss_sequence;

     //! Constructor
     SettingsNodeSs(std::string initial_ss_file="",
                    std::string initial_ss_sequence="")
          : initial_ss_file(initial_ss_file),                 
            initial_ss_sequence(initial_ss_sequence) {}

     //! Initialize a dbn with this Settings object
     //!
     //! \tparam NODE_INDEX index of node (boost::mpl::int_).
     //! \param dbn model object
     template <typename DBN_TYPE,typename NODE_INDEX>
     void initialize(DBN_TYPE &dbn) const {
          // Set SS input
          if (initial_ss_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX>(initial_ss_sequence, &definitions::str_to_ss);               
          else if (initial_ss_file != "") 
               dbn.template set_sequence_vector<NODE_INDEX>(file_to_string(initial_ss_file), &definitions::str_to_ss);
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const SettingsNodeSs &settings) {
          o << "initial-ss-file:" << settings.initial_ss_file << "\n";
          o << "initial-ss-sequence:" << settings.initial_ss_sequence << "\n";
          return o;
     }
};


//! Settings class for angles node
class SettingsNodeAngles {
public:
     //! Angle input file
     std::string initial_angle_file;

     //! Angle input sequence
     std::string initial_angle_sequence;

     //! Constructor
     SettingsNodeAngles(std::string initial_angle_file="",
                        std::string initial_angle_sequence="")
          : initial_angle_file(initial_angle_file),                 
            initial_angle_sequence(initial_angle_sequence) {}


     //! Initialize a dbn with this Settings object
     //!
     //! \tparam NODE_INDEX index of node (boost::mpl::int_).
     //! \param dbn model object
     //! \param pdb_data Optional vector with data from a pdb file
     template <typename DBN_TYPE, typename NODE_INDEX>
     void initialize(DBN_TYPE &dbn, const std::vector<std::vector<double> > &pdb_data=std::vector<std::vector<double> >()) const {

          // Set angle input
          if (initial_angle_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX>(initial_angle_sequence);
          else if (initial_angle_file != "") {
               dbn.template set_sequence_vector<NODE_INDEX>(file_to_string(initial_angle_file));
          } else if (!pdb_data.empty()) {
               dbn.template set_sequence_vector<NODE_INDEX>(pdb_data);
          }
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const SettingsNodeAngles &settings) {
          o << "initial-angle-file:" << settings.initial_angle_file << "\n";
          o << "initial-angle-sequence:" << settings.initial_angle_sequence << "\n";
          return o;
     }
};


//! Settings class for cis node
class SettingsNodeCis {
public:
     //! Whether to initialize cis input from PDB input
     bool initialize_cis_from_pdb;

     //! Cis input file
     std::string initial_cis_file;

     //! Cis input sequence
     std::string initial_cis_sequence;
     
     //! Constructor
     SettingsNodeCis(bool initialize_cis_from_pdb=false,
                     std::string initial_cis_file="",
                     std::string initial_cis_sequence="")
          : initialize_cis_from_pdb(initialize_cis_from_pdb),
            initial_cis_file(initial_cis_file),                 
            initial_cis_sequence(initial_cis_sequence) {}

     //! Initialize a dbn with this Settings object
     //!
     //! \tparam NODE_INDEX index of node (boost::mpl::int_).
     //! \param dbn model object
     //! \param pdb_data Optional vector with data from a pdb file
     //! \param fix_emission Whether to fix the emission state of this node
     template <typename DBN_TYPE, typename NODE_INDEX>
     void initialize(DBN_TYPE &dbn, const std::vector<int> &pdb_data=std::vector<int>(),
                     bool fix_emission=true) const {
          // Set cis input
          if (initial_cis_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX>(initial_cis_sequence, 0, fix_emission);               
          else if (initial_cis_file != "") 
               dbn.template set_sequence_vector<NODE_INDEX>(file_to_string(initial_cis_file), 0, fix_emission);
          else if (!pdb_data.empty()) {
               dbn.template set_sequence_vector<NODE_INDEX>(pdb_data, 0, fix_emission);
          }
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const SettingsNodeCis &settings) {
          o << "initialialize-cis-from-pdb:" << settings.initialize_cis_from_pdb << "\n";
          o << "initial-cis-file:" << settings.initial_cis_file << "\n";
          o << "initial-cis-sequence:" << settings.initial_cis_sequence << "\n";
          return o;
     }
};


//! Settings class for omega node
class SettingsNodeOmega {
public:
     //! Omega input file
     std::string initial_omega_file;

     //! Omega input sequence
     std::string initial_omega_sequence;

     //! Constructor
     SettingsNodeOmega(std::string initial_omega_file="",
                       std::string initial_omega_sequence="")
          : initial_omega_file(initial_omega_file),                 
            initial_omega_sequence(initial_omega_sequence) {}

     //! Initialize a dbn with this Settings object
     //!
     //! \tparam NODE_INDEX index of node (boost::mpl::int_).
     //! \param dbn model object
     //! \param pdb_data Optional vector with data from a pdb file
     template <typename DBN_TYPE, typename NODE_INDEX>
     void initialize(DBN_TYPE &dbn, const std::vector<double> &pdb_data=std::vector<double>()) const {
          // Set omega input
          if (initial_omega_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX>(initial_omega_sequence);
          else if (initial_omega_file != "") {
               dbn.template set_sequence_vector<NODE_INDEX>(file_to_string(initial_omega_file));
          } else if (!pdb_data.empty()) {
               dbn.template set_sequence_vector<NODE_INDEX>(pdb_data);
          }
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const SettingsNodeOmega &settings) {
          o << "initial-omega-file:" << settings.initial_omega_file << "\n";
          o << "initial-omega-sequence:" << settings.initial_omega_sequence << "\n";
          return o;
     }
};


//! Settings class for chemical shift nodes
class SettingsNodeCs {
public:
     //! Chemical shift CA value input file
     std::string initial_cs_ca_file;
     //! Chemical shift CA value input sequence
     std::string initial_cs_ca_sequence;
     //! Chemical shift CB value input file
     std::string initial_cs_cb_file;
     //! Chemical shift CB value input sequence
     std::string initial_cs_cb_sequence;
     //! Chemical shift C value input file
     std::string initial_cs_c_file;
     //! Chemical shift C value input sequence
     std::string initial_cs_c_sequence;
     //! Chemical shift N value input file
     std::string initial_cs_n_file;
     //! Chemical shift N value input sequence
     std::string initial_cs_n_sequence;
     //! Chemical shift HA value input file
     std::string initial_cs_ha_file;
     //! Chemical shift HA value input sequence
     std::string initial_cs_ha_sequence;
     //! Chemical shift H value input file
     std::string initial_cs_h_file;
     //! Chemical shift H value input sequence
     std::string initial_cs_h_sequence;

     //! Chemical shift values input file (all values in one go, order: CA,CB,C,N,HA,H)
     std::string initial_cs_file;
     //! Chemical shift values input sequence (all values in one go, order: CA,CB,C,N,HA,H)
     std::string initial_cs_sequence;
     //! Chemical shift values input sequence in NMR STAR FORMAT (all values in one go)
     std::string initial_cs_nmr_star_file;

     //! Offset of CA chemical shift values (for re-referencing)
     double initial_cs_ca_offset;

     //! Offset of CB chemical shift values (for re-referencing)
     double initial_cs_cb_offset;

     //! Offset of C chemical shift values (for re-referencing)
     double initial_cs_c_offset;

     //! Offset of N chemical shift values (for re-referencing)
     double initial_cs_n_offset;

     //! Offset of HA chemical shift values (for re-referencing)
     double initial_cs_ha_offset;

     //! Offset of H chemical shift values (for re-referencing)
     double initial_cs_h_offset;

     //! Constructor
     SettingsNodeCs(std::string initial_cs_ca_file="",
                    std::string initial_cs_ca_sequence="",
                    std::string initial_cs_cb_file="",
                    std::string initial_cs_cb_sequence="",
                    std::string initial_cs_c_file="",
                    std::string initial_cs_c_sequence="",
                    std::string initial_cs_n_file="",
                    std::string initial_cs_n_sequence="",
                    std::string initial_cs_ha_file="",
                    std::string initial_cs_ha_sequence="",
                    std::string initial_cs_h_file="",
                    std::string initial_cs_h_sequence="",
                    std::string initial_cs_file="",
                    std::string initial_cs_sequence="",
                    std::string initial_cs_nmr_star_file="",
                    double initial_cs_ca_offset=0.,
                    double initial_cs_cb_offset=0.,
                    double initial_cs_c_offset=0.,
                    double initial_cs_n_offset=0.,
                    double initial_cs_ha_offset=0.,
                    double initial_cs_h_offset=0.)
          : initial_cs_ca_file(initial_cs_ca_file),                 
            initial_cs_ca_sequence(initial_cs_ca_sequence),
            initial_cs_cb_file(initial_cs_cb_file),                 
            initial_cs_cb_sequence(initial_cs_cb_sequence), 
            initial_cs_c_file(initial_cs_c_file),                 
            initial_cs_c_sequence(initial_cs_ca_sequence), 
            initial_cs_n_file(initial_cs_n_file),                 
            initial_cs_n_sequence(initial_cs_n_sequence), 
            initial_cs_ha_file(initial_cs_ha_file),                 
            initial_cs_ha_sequence(initial_cs_ha_sequence), 
            initial_cs_h_file(initial_cs_h_file),                 
            initial_cs_h_sequence(initial_cs_h_sequence),
            initial_cs_file(initial_cs_file),                 
            initial_cs_sequence(initial_cs_sequence),
            initial_cs_nmr_star_file(initial_cs_nmr_star_file),
            initial_cs_ca_offset(initial_cs_ca_offset),
            initial_cs_cb_offset(initial_cs_cb_offset),
            initial_cs_c_offset(initial_cs_c_offset),
            initial_cs_n_offset(initial_cs_n_offset),
            initial_cs_ha_offset(initial_cs_ha_offset),
            initial_cs_h_offset(initial_cs_h_offset) {}


     //! Initialize a dbn with this Settings object
     //!
     //! \tparam NODE_INDEX_CA index of CA node (boost::mpl::int_).
     //! \tparam NODE_INDEX_CB index of CB node (boost::mpl::int_).
     //! \tparam NODE_INDEX_C index of C node (boost::mpl::int_).
     //! \tparam NODE_INDEX_N index of N node (boost::mpl::int_).
     //! \tparam NODE_INDEX_HA index of HA node (boost::mpl::int_).
     //! \tparam NODE_INDEX_H index of H node (boost::mpl::int_).
     //! \tparam NODE_INDEX_ALL list of indices of all cs nodes.
     //! \param dbn model object
     template <typename DBN_TYPE, 
               typename NODE_INDEX_CA, typename NODE_INDEX_CB, typename NODE_INDEX_C, 
               typename NODE_INDEX_N, typename NODE_INDEX_HA, typename NODE_INDEX_H,
               typename NODE_INDEX_ALL>
     void initialize(DBN_TYPE &dbn) const {
          // Set Chemical shift input
          if (initial_cs_ca_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX_CA>(initial_cs_ca_sequence);
          else if (initial_cs_ca_file != "")
               dbn.template set_sequence_vector<NODE_INDEX_CA>(file_to_string(initial_cs_ca_file));
          if (initial_cs_cb_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX_CB>(initial_cs_cb_sequence);
          else if (initial_cs_cb_file != "")
               dbn.template set_sequence_vector<NODE_INDEX_CB>(file_to_string(initial_cs_cb_file));
          if (initial_cs_c_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX_C>(initial_cs_c_sequence);
          else if (initial_cs_c_file != "")
               dbn.template set_sequence_vector<NODE_INDEX_C>(file_to_string(initial_cs_c_file));
          if (initial_cs_n_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX_N>(initial_cs_n_sequence);
          else if (initial_cs_n_file != "")
               dbn.template set_sequence_vector<NODE_INDEX_N>(file_to_string(initial_cs_n_file));
          if (initial_cs_ha_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX_HA>(initial_cs_ha_sequence);
          else if (initial_cs_ha_file != "")
               dbn.template set_sequence_vector<NODE_INDEX_HA>(file_to_string(initial_cs_ha_file));
          if (initial_cs_h_sequence != "")
               dbn.template set_sequence_vector<NODE_INDEX_H>(initial_cs_h_sequence);
          else if (initial_cs_h_file != "")
               dbn.template set_sequence_vector<NODE_INDEX_H>(file_to_string(initial_cs_h_file));

          // all in one go.
          if (initial_cs_sequence != "") {
               dbn.template set_sequence_vector<NODE_INDEX_ALL>(initial_cs_sequence);
          } else if (initial_cs_file != "") {
               dbn.template set_sequence_vector<NODE_INDEX_ALL>(file_to_string(initial_cs_file));
          } else if (initial_cs_nmr_star_file != "") {

               // Extract amino acid sequence from model - for verification of consistency
               std::vector<int> aa_sequence;
               if (dbn.template get_node<typename DBN_TYPE::AA_NODE>()->fixed) {
                    aa_sequence = dbn.template get_node<typename DBN_TYPE::AA_NODE>()->get_sequence_vector();
               } else {
                    std::cout << "WARNING (BackboneDBN): Using chemical shift input in DBN without using amino acid information. Include both for best performance.\n";
               }

               std::vector<std::vector<double> > value_matrix;
               std::vector<std::string> lines = file_to_string_vector(initial_cs_nmr_star_file);
               int residue_index_current = uninitialized<int>();
               int node_index_offset = NODE_INDEX_CA::value;
               for (unsigned int i=0; i<lines.size(); ++i) {
                    std::string line = lines[i];
                    // Remove white space at beginning and end
                    boost::trim(line);
                    if (line.size() == 0)
                         continue;
                    std::vector<std::string> tokens;
                    boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);

                    // Detect whether this is a new residue
                    int residue_index = boost::lexical_cast<int>(tokens[1]);
                    if (residue_index != residue_index_current) {
                         int offset= (residue_index-residue_index_current);
                         if (!is_initialized(residue_index_current)) {
                              offset = 1;
                         }
                         value_matrix.resize(value_matrix.size()+offset, std::vector<double>(6, UNINITIALIZED));
                         residue_index_current = residue_index;
                    }
                    std::string atom_type = tokens[3];

                    // Attempt to convert CS value to a double
                    double value = uninitialized<double>();
                    if (tokens.size() > 5) {
                         std::stringstream ss(tokens[5]);
                         if ((ss >> value).fail() || !(ss >> std::ws).eof()) {
                              value = uninitialized<double>();
                         }
                    }
                              
                    if (atom_type == "CA") {
                         value_matrix.back()[NODE_INDEX_CA::value-node_index_offset] = value;
                    } else if (atom_type == "CB") {
                         value_matrix.back()[NODE_INDEX_CB::value-node_index_offset] = value;
                    } else if (atom_type == "C") {
                         value_matrix.back()[NODE_INDEX_C::value-node_index_offset] = value;
                    } else if (atom_type == "N") {
                         value_matrix.back()[NODE_INDEX_N::value-node_index_offset] = value;
                    } else if (atom_type == "HA") {
                         value_matrix.back()[NODE_INDEX_HA::value-node_index_offset] = value;
                    } else if (atom_type == "H") {
                         value_matrix.back()[NODE_INDEX_H::value-node_index_offset] = value;
                    } else {
                         // std::cout << "WARNING - NMR-STAR parser: Skipping unknown atom type " << atom_type << ".\n";
                    }

                    // Check for consistency with AA-sequence in model
                    if (!aa_sequence.empty()) {
                         std::string residue = tokens[2];
                         int aa = definitions::str_to_aa(residue);
                         if (aa != aa_sequence[value_matrix.size()-1]) {
                              std::cout << "WARNING (BackboneDBN): Amino acid mismatch in DBN at position " <<
                                   value_matrix.size()-1 << ": " << definitions::ResidueEnum(aa) << "!=" << 
                                   definitions::ResidueEnum(aa_sequence[value_matrix.size()-1]) << "\n";
                         }
                    }
               }

               // In case of missing values at the end, we explicitly add uninitialized values
               if ((int)value_matrix.size() < dbn.sequence_length) {
                    value_matrix.resize(dbn.sequence_length, std::vector<double>(6, UNINITIALIZED));
               }

               dbn.template set_sequence_vector<NODE_INDEX_ALL>(value_matrix);
          }

          dbn.template add_offset<NODE_INDEX_CA>(initial_cs_ca_offset);
          dbn.template add_offset<NODE_INDEX_CB>(initial_cs_cb_offset);
          dbn.template add_offset<NODE_INDEX_C>(initial_cs_c_offset);
          dbn.template add_offset<NODE_INDEX_N>(initial_cs_n_offset);
          dbn.template add_offset<NODE_INDEX_HA>(initial_cs_ha_offset);
          dbn.template add_offset<NODE_INDEX_H>(initial_cs_h_offset);
     }

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const SettingsNodeCs &settings) {
          o << "initial-cs-ca-file:" << settings.initial_cs_ca_file << "\n";
          o << "initial-cs-ca-sequence:" << settings.initial_cs_ca_sequence << "\n";
          o << "initial-cs-cb-file:" << settings.initial_cs_cb_file << "\n";
          o << "initial-cs-cb-sequence:" << settings.initial_cs_cb_sequence << "\n";
          o << "initial-cs-c-file:" << settings.initial_cs_c_file << "\n";
          o << "initial-cs-c-sequence:" << settings.initial_cs_c_sequence << "\n";
          o << "initial-cs-n-file:" << settings.initial_cs_n_file << "\n";
          o << "initial-cs-n-sequence:" << settings.initial_cs_n_sequence << "\n";
          o << "initial-cs-ha-file:" << settings.initial_cs_ha_file << "\n";
          o << "initial-cs-ha-sequence:" << settings.initial_cs_ha_sequence << "\n";
          o << "initial-cs-h-file:" << settings.initial_cs_h_file << "\n";
          o << "initial-cs-h-sequence:" << settings.initial_cs_h_sequence << "\n";
          o << "initial-cs-file:" << settings.initial_cs_file << "\n";
          o << "initial-cs-sequence:" << settings.initial_cs_sequence << "\n";
          o << "initial-cs-nmr-star-file:" << settings.initial_cs_nmr_star_file << "\n";
          o << "initial-cs-ca-offset:" << settings.initial_cs_ca_offset << "\n";
          o << "initial-cs-cb-offset:" << settings.initial_cs_cb_offset << "\n";
          o << "initial-cs-c-offset:" << settings.initial_cs_c_offset << "\n";
          o << "initial-cs-n-offset:" << settings.initial_cs_n_offset << "\n";
          o << "initial-cs-ha-offset:" << settings.initial_cs_ha_offset << "\n";
          o << "initial-cs-h-offset:" << settings.initial_cs_h_offset << "\n";
          return o;
     }
};


// Forward declaration
class TorusDBN;

//! Definition of node structure for TorusDBN
typedef boost::fusion::vector<HiddenNode  <TorusDBN, 0>*,
                              DiscreteNode<TorusDBN, 1>*,
                              DiscreteNode<TorusDBN, 2>*,
                              TorusNode   <TorusDBN, 3>*, 
                              DiscreteNode<TorusDBN, 4>* > TorusDBNNodes;

//! TorusDBN model class. Inherits from BackboneDBN passing the appropriate nodes as a template argument.
class TorusDBN: public BackboneDBN<TorusDBN,TorusDBNNodes> {
public:

     //! Model name
     static const char name[];

     //@{ 
     //! Node indices (must match the order in the node definition)
     typedef boost::mpl::int_<1> AA_NODE;
     typedef boost::mpl::int_<2> SS_NODE;
     typedef boost::mpl::int_<3> ANGLE_NODE;
     typedef boost::mpl::int_<4> CIS_NODE;
     //@}

     //! List of all angle node indices
     typedef boost::mpl::list<ANGLE_NODE> ALL_ANGLE_NODES;

     //! Local Settings object
     const class Settings: 
          public BackboneDBN<TorusDBN,TorusDBNNodes>::Settings,
          public SettingsNodeAa,
          public SettingsNodeSs,
          public SettingsNodeAngles,
          public SettingsNodeCis {

     public:
          
          //! Initialize a dbn with this Settings object
          //! \param dbn DBN model object
          void initialize(TorusDBN &dbn) const {

               // Initialize from sequence length option (if available)
               if (sequence_length > 0)
                    dbn.init(sequence_length);

               ProteinData pdb_data(this->initial_pdb_file);

               SettingsNodeAa::initialize<TorusDBN,AA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_sequence()[this->initial_pdb_file_index] : 
                     std::vector<int>()));
               SettingsNodeSs::initialize<TorusDBN,SS_NODE>(dbn);
               SettingsNodeAngles::initialize<TorusDBN,ALL_ANGLE_NODES>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_phi_psi()[this->initial_pdb_file_index] :
                     std::vector<std::vector<double> >()));

               SettingsNodeCis::initialize<TorusDBN, CIS_NODE>(
                    dbn,
                    (pdb_data && initialize_cis_from_pdb ? 
                     pdb_data.get_cis()[this->initial_pdb_file_index] : 
                     std::vector<int>()));

          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const BackboneDBN<TorusDBN,TorusDBNNodes>::Settings &>(settings);
               o << static_cast<const SettingsNodeAa &>(settings);
               o << static_cast<const SettingsNodeSs &>(settings);
               o << static_cast<const SettingsNodeAngles &>(settings);
               o << static_cast<const SettingsNodeCis &>(settings);
               return o;
          }
     } settings;


     //! Constructor
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be constructed
     //! \param parameters Model parameters for the TorusDBN
     TorusDBN(const Settings &settings = Settings(),
              RandomNumberEngine *random_number_engine = &random_global,
              Parameters parameters = default_parameters_torus_dbn):
          BackboneDBN<TorusDBN,TorusDBNNodes>(
               new HiddenNode<TorusDBN,0>("HIDDEN", parameters, this,
                                          settings.start_distribution, 
                                          settings.transition_distribution), 
               boost::fusion::make_vector(
                    new DiscreteNode<TorusDBN,1>("AA", parameters, this, random_number_engine), 
                    new DiscreteNode<TorusDBN,2>("SS", parameters, this, random_number_engine), 
                    new TorusNode<TorusDBN,3>("TORUS", parameters, this, random_number_engine),
                    new DiscreteNode<TorusDBN,4>("CIS", parameters, this, random_number_engine)),
               settings),
          settings(settings) {

          // Set initial values in DBN
          settings.initialize(*this);
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     TorusDBN(const TorusDBN &other)
          : BackboneDBN<TorusDBN,TorusDBNNodes>(other),
            settings(other.settings) {          
     }
     
     //! Copy constructor - using different random number engine.
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be constructed.
     TorusDBN(const TorusDBN &other,
              RandomNumberEngine *random_number_engine)
          : BackboneDBN<TorusDBN,TorusDBNNodes>(other, random_number_engine),
            settings(other.settings) {
     }
};
//! Name of TorusDBN model
const char TorusDBN::name[] = "torus";


// Forward declaration
class TorusOmegaDBN;

//! Definition of node structure for TorusOmegaDBN
typedef boost::fusion::vector<HiddenNode  <TorusOmegaDBN, 0>*,
                              DiscreteNode<TorusOmegaDBN, 1>*,
                              DiscreteNode<TorusOmegaDBN, 2>*,
                              TorusNode   <TorusOmegaDBN, 3>*, 
                              DiscreteNode<TorusOmegaDBN, 4>*,
                              VonMisesNode<TorusOmegaDBN, 5, 4>* > TorusOmegaDBNNodes;

//! TorusOmegaDBN model class. Inherits from BackboneDBN passing the appropriate nodes as a template argument.
class TorusOmegaDBN: public BackboneDBN<TorusOmegaDBN,TorusOmegaDBNNodes> {
public:

     //! Model name
     static const char name[];

     //@{ 
     //! Node indices (must match the order in the node definition)
     typedef boost::mpl::int_<1> AA_NODE;
     typedef boost::mpl::int_<2> SS_NODE;
     typedef boost::mpl::int_<3> ANGLE_NODE;
     typedef boost::mpl::int_<4> CIS_NODE;
     typedef boost::mpl::int_<5> OMEGA_NODE;
     //@}

     //! List of all angle node indices
     typedef boost::mpl::list<ANGLE_NODE,OMEGA_NODE> ALL_ANGLE_NODES;

     //! Local Settings object
     const class Settings: 
          public BackboneDBN<TorusOmegaDBN,TorusOmegaDBNNodes>::Settings,
          public SettingsNodeAa,
          public SettingsNodeSs,
          public SettingsNodeAngles,
          public SettingsNodeCis,
          public SettingsNodeOmega {

     public:
          
          //! Initialize a dbn with this Settings object
          //! \param dbn DBN model object
          void initialize(TorusOmegaDBN &dbn) const {

               // Initialize from sequence length option (if available)
               if (sequence_length > 0)
                    dbn.init(sequence_length);

               ProteinData pdb_data(this->initial_pdb_file);

               SettingsNodeAa::initialize<TorusOmegaDBN,AA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_sequence()[this->initial_pdb_file_index] : 
                     std::vector<int>()));
               SettingsNodeSs::initialize<TorusOmegaDBN,SS_NODE>(dbn);
               SettingsNodeAngles::initialize<TorusOmegaDBN,ALL_ANGLE_NODES>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_phi_psi()[this->initial_pdb_file_index] :
                     std::vector<std::vector<double> >()));

               // Do not fix emission for Cis node (omega is included instead)
               bool fix_emission = false;
               SettingsNodeCis::initialize<TorusOmegaDBN, CIS_NODE>(
                    dbn,
                    (pdb_data ? 
                     pdb_data.get_cis()[this->initial_pdb_file_index] : 
                     std::vector<int>()), fix_emission);
               SettingsNodeOmega::initialize<TorusOmegaDBN, OMEGA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_omega()[this->initial_pdb_file_index] : 
                     std::vector<double>()));
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const BackboneDBN<TorusOmegaDBN,TorusOmegaDBNNodes>::Settings &>(settings);
               o << static_cast<const SettingsNodeAa &>(settings);
               o << static_cast<const SettingsNodeSs &>(settings);
               o << static_cast<const SettingsNodeAngles &>(settings);
               o << static_cast<const SettingsNodeCis &>(settings);
               o << static_cast<const SettingsNodeOmega &>(settings);
               return o;
          }
     } settings;


     //! Constructor
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be constructed
     //! \param parameters Model parameters for the TorusDBN
     TorusOmegaDBN(const Settings &settings = Settings(),
              RandomNumberEngine *random_number_engine = &random_global,
              Parameters parameters = default_parameters_torus_dbn):
          BackboneDBN<TorusOmegaDBN,TorusOmegaDBNNodes>(
               new HiddenNode<TorusOmegaDBN,0>("HIDDEN", parameters, this,
                                          settings.start_distribution, 
                                          settings.transition_distribution), 
               boost::fusion::make_vector(
                    new DiscreteNode<TorusOmegaDBN,1>("AA", parameters, this, random_number_engine), 
                    new DiscreteNode<TorusOmegaDBN,2>("SS", parameters, this, random_number_engine), 
                    new TorusNode<TorusOmegaDBN,3>("TORUS", parameters, this, random_number_engine),
                    new DiscreteNode<TorusOmegaDBN,4>("CIS", parameters, this, random_number_engine), 
                    new VonMisesNode<TorusOmegaDBN,5,4>("OMEGA", 
                                                   // kappa = 1/sigma^2
                                                   vector_utils::make_vector<double>(M_PI, 0.0), 
                                                   vector_utils::make_vector<double>(1.0/Math<double>::sqr(5.0/180*M_PI),
                                                                                     1.0/Math<double>::sqr(5.0/180*M_PI)),
                                                   this, random_number_engine)),
               settings),
          settings(settings) {

          // Set initial values in DBN
          settings.initialize(*this);
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     TorusOmegaDBN(const TorusOmegaDBN &other)
          : BackboneDBN<TorusOmegaDBN,TorusOmegaDBNNodes>(other),
            settings(other.settings) {          
     }
     
     //! Copy constructor - using different random number engine.
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be constructed.
     TorusOmegaDBN(const TorusOmegaDBN &other,
              RandomNumberEngine *random_number_engine)
          : BackboneDBN<TorusOmegaDBN,TorusOmegaDBNNodes>(other, random_number_engine),
            settings(other.settings) {
     }
};
//! Name of TorusOmegaDBN model
const char TorusOmegaDBN::name[] = "torus-omega";




// Forward declaration
class FB5DBN;

//! Definition of node structure for FB5DBN
typedef boost::fusion::vector<HiddenNode  <FB5DBN,0>*,
                              DiscreteNode<FB5DBN,1>*,
                              DiscreteNode<FB5DBN,2>*,
                              Fb5Node     <FB5DBN,3>* > Fb5DBNNodes;

//! FB5DBN model class. Inherits from BackboneDBN passing the appropriate nodes as a template argument.
class FB5DBN:public BackboneDBN<FB5DBN,Fb5DBNNodes> {
public:

     //! Model name
     static const char name[];

     //@{ 
     //! Node indices (must match the order in the node definition)
     typedef boost::mpl::int_<1> AA_NODE;
     typedef boost::mpl::int_<2> SS_NODE;
     typedef boost::mpl::int_<3> ANGLE_NODE;
     //@}

     //! List of all angle node indices
     typedef boost::mpl::list<ANGLE_NODE> ALL_ANGLE_NODES;

     //! Local Settings object
     const class Settings: 
          public BackboneDBN<FB5DBN,Fb5DBNNodes>::Settings,
          public SettingsNodeAa,
          public SettingsNodeSs,
          public SettingsNodeAngles {
     public:

          //! Initialize a dbn with this Settings object
          //! \param dbn DBN model object
          void initialize(FB5DBN &dbn) const {

               // Initialize from sequence length option (if available)
               if (sequence_length > 0)
                    dbn.init(sequence_length);

               ProteinData pdb_data(this->initial_pdb_file);

               SettingsNodeAa::initialize<FB5DBN,AA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_sequence()[this->initial_pdb_file_index] : 
                     std::vector<int>()));
               SettingsNodeSs::initialize<FB5DBN,SS_NODE>(dbn);
               SettingsNodeAngles::initialize<FB5DBN,ALL_ANGLE_NODES>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_phi_psi()[this->initial_pdb_file_index] :
                     std::vector<std::vector<double> >()));
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const BackboneDBN<FB5DBN,Fb5DBNNodes>::Settings &>(settings);
               o << static_cast<const SettingsNodeAa &>(settings);
               o << static_cast<const SettingsNodeSs &>(settings);
               o << static_cast<const SettingsNodeAngles &>(settings);
               return o;
          }
     } settings;

     //! Constructor
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be constructed
     //! \param parameters Model parameters for the FB5DBN
     FB5DBN(const Settings &settings = Settings(),
            RandomNumberEngine *random_number_engine = &random_global,
            Parameters parameters = default_parameters_fb5_dbn):
          BackboneDBN<FB5DBN,Fb5DBNNodes>(
               new HiddenNode<FB5DBN,0>("HIDDEN", parameters, this,
                                        settings.start_distribution, 
                                        settings.transition_distribution), 
               boost::fusion::make_vector(
                    new DiscreteNode<FB5DBN,1>("AA", parameters, this, random_number_engine), 
                    new DiscreteNode<FB5DBN,2>("SS", parameters, this, random_number_engine), 
                    new Fb5Node<FB5DBN,3>("FB5", parameters, this, random_number_engine)),
               settings),
          settings(settings) {

          // Set initial values in DBN
          settings.initialize(*this);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     FB5DBN(const FB5DBN &other)
          : BackboneDBN<FB5DBN,Fb5DBNNodes>(other),
            settings(other.settings) {
     }
     
     //! Copy constructor - using different random number engine.
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be constructed.
     FB5DBN(const FB5DBN &other,
              RandomNumberEngine *random_number_engine)
          : BackboneDBN<FB5DBN,Fb5DBNNodes>(other, random_number_engine),
            settings(other.settings) {
     }
};
//! Name of TorusDBN model
const char FB5DBN::name[] = "fb5";




// Forward declaration
class TorusCsDBN;

//! Definition of node structure for TorusCsDBN
typedef boost::fusion::vector<HiddenNode  <TorusCsDBN, 0>*,
                              DiscreteNode<TorusCsDBN, 1>*,
                              DiscreteNode<TorusCsDBN, 2>*,
                              TorusNode   <TorusCsDBN, 3>*,
                              DiscreteNode<TorusCsDBN, 4>*,
                              GaussianNode<TorusCsDBN, 5>*,
                              GaussianNode<TorusCsDBN, 6>*,
                              GaussianNode<TorusCsDBN, 7>*,
                              GaussianNode<TorusCsDBN, 8>*,
                              GaussianNode<TorusCsDBN, 9>*,
                              GaussianNode<TorusCsDBN,10>* > TorusCsDBNNodes;


//! TorusCsDBN model class. Inherits from BackboneDBN passing the appropriate nodes as a template argument.
class TorusCsDBN:public BackboneDBN<TorusCsDBN,
                                    TorusCsDBNNodes> {
public:

     //! Model name
     static const char name[];

     //@{ 
     //! Node indices (must match the order in the node definition)
     typedef boost::mpl::int_< 1> AA_NODE;
     typedef boost::mpl::int_< 2> SS_NODE;
     typedef boost::mpl::int_< 3> ANGLE_NODE;
     typedef boost::mpl::int_< 4> CIS_NODE;
     typedef boost::mpl::int_< 5> CS_CA_NODE;
     typedef boost::mpl::int_< 6> CS_CB_NODE;
     typedef boost::mpl::int_< 7> CS_C_NODE;
     typedef boost::mpl::int_< 8> CS_N_NODE;
     typedef boost::mpl::int_< 9> CS_HA_NODE;
     typedef boost::mpl::int_<10> CS_H_NODE;
     //@}

     //! List of all angle node indices
     typedef boost::mpl::list<ANGLE_NODE> ALL_ANGLE_NODES;

     //! List of all chemical shift node indices
     typedef boost::mpl::list<CS_CA_NODE,CS_CB_NODE,CS_C_NODE,CS_N_NODE,CS_HA_NODE,CS_H_NODE> ALL_CS_NODES;

     //! Local Settings object
     const class Settings: 
          public BackboneDBN<TorusCsDBN,TorusCsDBNNodes>::Settings,
          public SettingsNodeAa,
          public SettingsNodeSs,
          public SettingsNodeAngles,
          public SettingsNodeCis,
          public SettingsNodeCs {
     public:

          //! Initialize a dbn with this Settings object
          //! \param dbn DBN model object
          void initialize(TorusCsDBN &dbn) const {

               // Initialize from sequence length option (if available)
               if (sequence_length > 0)
                    dbn.init(sequence_length);

               ProteinData pdb_data(this->initial_pdb_file);

               SettingsNodeAa::initialize<TorusCsDBN,AA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_sequence()[this->initial_pdb_file_index] : 
                     std::vector<int>()));
               SettingsNodeSs::initialize<TorusCsDBN,SS_NODE>(dbn);
               SettingsNodeAngles::initialize<TorusCsDBN,ALL_ANGLE_NODES>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_phi_psi()[this->initial_pdb_file_index] :
                     std::vector<std::vector<double> >()));
               SettingsNodeCis::initialize<TorusCsDBN, CIS_NODE>(
                    dbn,
                    (pdb_data ? 
                     pdb_data.get_cis()[this->initial_pdb_file_index] : 
                     std::vector<int>()));
               SettingsNodeCs::initialize<TorusCsDBN, 
                                          CS_CA_NODE, CS_CB_NODE, CS_C_NODE,
                                          CS_N_NODE, CS_HA_NODE, CS_H_NODE,
                                          ALL_CS_NODES>(dbn);
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const BackboneDBN<TorusCsDBN,TorusCsDBNNodes>::Settings &>(settings);
               o << static_cast<const SettingsNodeAa &>(settings);
               o << static_cast<const SettingsNodeSs &>(settings);
               o << static_cast<const SettingsNodeAngles &>(settings);
               o << static_cast<const SettingsNodeCis &>(settings);
               o << static_cast<const SettingsNodeCs &>(settings);
               return o;
          }
     } settings;


     //! Constructor
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be constructed
     //! \param parameters Model parameters for the TorusCsDBN
     TorusCsDBN(const Settings &settings = Settings(),
                RandomNumberEngine *random_number_engine = &random_global,
                Parameters parameters = default_parameters_torus_cs_dbn):
          BackboneDBN<TorusCsDBN,TorusCsDBNNodes>(
               new HiddenNode<TorusCsDBN,0>("HIDDEN", parameters, this,
                                            settings.start_distribution, 
                                            settings.transition_distribution), 
               boost::fusion::make_vector(
                    new DiscreteNode<TorusCsDBN,1>("AA", parameters, this, random_number_engine), 
                    new DiscreteNode<TorusCsDBN,2>("SS", parameters, this, random_number_engine), 
                    new TorusNode   <TorusCsDBN,3>("TORUS", parameters, this, random_number_engine),
                    new DiscreteNode<TorusCsDBN,4>("CIS", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsDBN,5>("CS-CA", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsDBN,6>("CS-CB", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsDBN,7>("CS-C", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsDBN,8>("CS-N", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsDBN,9>("CS-HA", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsDBN,10>("CS-H", parameters, this, random_number_engine)),
               settings),
          settings(settings) {

          // Set initial values in DBN
          settings.initialize(*this);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     TorusCsDBN(const TorusCsDBN &other)
          : BackboneDBN<TorusCsDBN,TorusCsDBNNodes>(other),
            settings(other.settings) {
     }
     
     //! Copy constructor - using different random number engine.
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be constructed.
     TorusCsDBN(const TorusCsDBN &other,
              RandomNumberEngine *random_number_engine)
          : BackboneDBN<TorusCsDBN,TorusCsDBNNodes>(other, random_number_engine),
            settings(other.settings) {
     }
};
//! Name of TorusCsDBN model
const char TorusCsDBN::name[] = "torus-cs";



// Forward declaration
class TorusCsOmegaDBN;

//! Definition of node structure for TorusCsOmegaDBN
typedef boost::fusion::vector<HiddenNode  <TorusCsOmegaDBN, 0>*,
                              DiscreteNode<TorusCsOmegaDBN, 1>*,
                              DiscreteNode<TorusCsOmegaDBN, 2>*,
                              TorusNode   <TorusCsOmegaDBN, 3>*,
                              DiscreteNode<TorusCsOmegaDBN, 4>*,
                              VonMisesNode<TorusCsOmegaDBN, 5, 4>*,
                              GaussianNode<TorusCsOmegaDBN, 6>*,
                              GaussianNode<TorusCsOmegaDBN, 7>*,
                              GaussianNode<TorusCsOmegaDBN, 8>*,
                              GaussianNode<TorusCsOmegaDBN, 9>*,
                              GaussianNode<TorusCsOmegaDBN,10>*,
                              GaussianNode<TorusCsOmegaDBN,11>* > TorusCsOmegaDBNNodes;


//! TorusCsOmegaDBN model class. Inherits from BackboneDBN passing the appropriate nodes as a template argument.
class TorusCsOmegaDBN:public BackboneDBN<TorusCsOmegaDBN,
                                    TorusCsOmegaDBNNodes> {
public:

     //! Model name
     static const char name[];

     //@{ 
     //! Node indices (must match the order in the node definition)
     typedef boost::mpl::int_< 1> AA_NODE;
     typedef boost::mpl::int_< 2> SS_NODE;
     typedef boost::mpl::int_< 3> ANGLE_NODE;
     typedef boost::mpl::int_< 4> CIS_NODE;
     typedef boost::mpl::int_< 5> OMEGA_NODE;
     typedef boost::mpl::int_< 6> CS_CA_NODE;
     typedef boost::mpl::int_< 7> CS_CB_NODE;
     typedef boost::mpl::int_< 8> CS_C_NODE;
     typedef boost::mpl::int_< 9> CS_N_NODE;
     typedef boost::mpl::int_<10> CS_HA_NODE;
     typedef boost::mpl::int_<11> CS_H_NODE;
     //@}

     //! List of all angle node indices
     typedef boost::mpl::list<ANGLE_NODE,OMEGA_NODE> ALL_ANGLE_NODES;

     //! List of all chemical shift node indices
     typedef boost::mpl::list<CS_CA_NODE,CS_CB_NODE,CS_C_NODE,CS_N_NODE,CS_HA_NODE,CS_H_NODE> ALL_CS_NODES;

     //! Local Settings object
     const class Settings: 
          public BackboneDBN<TorusCsOmegaDBN,TorusCsOmegaDBNNodes>::Settings,
          public SettingsNodeAa,
          public SettingsNodeSs,
          public SettingsNodeAngles,
          public SettingsNodeCis,
          public SettingsNodeOmega,
          public SettingsNodeCs {
     public:

          //! Initialize a dbn with this Settings object
          //! \param dbn DBN model object
          void initialize(TorusCsOmegaDBN &dbn) const {

               // Initialize from sequence length option (if available)
               if (sequence_length > 0)
                    dbn.init(sequence_length);

               ProteinData pdb_data(this->initial_pdb_file);

               SettingsNodeAa::initialize<TorusCsOmegaDBN,AA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_sequence()[this->initial_pdb_file_index] : 
                     std::vector<int>()));
               SettingsNodeSs::initialize<TorusCsOmegaDBN,SS_NODE>(dbn);
               SettingsNodeAngles::initialize<TorusCsOmegaDBN,ALL_ANGLE_NODES>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_phi_psi()[this->initial_pdb_file_index] :
                     std::vector<std::vector<double> >()));
               // Do not fix emission for Cis node (omega is included instead)
               bool fix_emission = false;
               SettingsNodeCis::initialize<TorusCsOmegaDBN, CIS_NODE>(
                    dbn,
                    (pdb_data ? 
                     pdb_data.get_cis()[this->initial_pdb_file_index] : 
                     std::vector<int>()), fix_emission);
               SettingsNodeOmega::initialize<TorusCsOmegaDBN, OMEGA_NODE>(
                    dbn, 
                    (pdb_data ? 
                     pdb_data.get_omega()[this->initial_pdb_file_index] : 
                     std::vector<double>()));
               SettingsNodeCs::initialize<TorusCsOmegaDBN, 
                                          CS_CA_NODE, CS_CB_NODE, CS_C_NODE,
                                          CS_N_NODE, CS_HA_NODE, CS_H_NODE,
                                          ALL_CS_NODES>(dbn);
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const BackboneDBN<TorusCsOmegaDBN,TorusCsOmegaDBNNodes>::Settings &>(settings);
               o << static_cast<const SettingsNodeAa &>(settings);
               o << static_cast<const SettingsNodeSs &>(settings);
               o << static_cast<const SettingsNodeAngles &>(settings);
               o << static_cast<const SettingsNodeCis &>(settings);
               o << static_cast<const SettingsNodeOmega &>(settings);
               o << static_cast<const SettingsNodeCs &>(settings);
               return o;
          }
     } settings;


     //! Constructor
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be constructed
     //! \param parameters Model parameters for the TorusCsDBN
     TorusCsOmegaDBN(const Settings &settings = Settings(),
                RandomNumberEngine *random_number_engine = &random_global,
                Parameters parameters = default_parameters_torus_cs_dbn):
          BackboneDBN<TorusCsOmegaDBN,TorusCsOmegaDBNNodes>(
               new HiddenNode<TorusCsOmegaDBN,0>("HIDDEN", parameters, this,
                                            settings.start_distribution, 
                                            settings.transition_distribution), 
               boost::fusion::make_vector(
                    new DiscreteNode<TorusCsOmegaDBN,1>("AA", parameters, this, random_number_engine), 
                    new DiscreteNode<TorusCsOmegaDBN,2>("SS", parameters, this, random_number_engine), 
                    new TorusNode   <TorusCsOmegaDBN,3>("TORUS", parameters, this, random_number_engine),
                    new DiscreteNode<TorusCsOmegaDBN,4>("CIS", parameters, this, random_number_engine),
                    new VonMisesNode<TorusCsOmegaDBN,5,4>("OMEGA", 
                                                     // kappa = 1/sigma^2
                                                     vector_utils::make_vector<double>(M_PI, 0.0), 
                                                     vector_utils::make_vector<double>(1.0/Math<double>::sqr(5.0/180*M_PI),
                                                                                       1.0/Math<double>::sqr(5.0/180*M_PI)),
                                                     this, random_number_engine),
                    new GaussianNode<TorusCsOmegaDBN,6>("CS-CA", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsOmegaDBN,7>("CS-CB", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsOmegaDBN,8>("CS-C", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsOmegaDBN,9>("CS-N", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsOmegaDBN,10>("CS-HA", parameters, this, random_number_engine),
                    new GaussianNode<TorusCsOmegaDBN,11>("CS-H", parameters, this, random_number_engine)),
               settings),
          settings(settings) {

          // Set initial values in DBN
          settings.initialize(*this);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     TorusCsOmegaDBN(const TorusCsOmegaDBN &other)
          : BackboneDBN<TorusCsOmegaDBN,TorusCsOmegaDBNNodes>(other),
            settings(other.settings) {
     }
     
     //! Copy constructor - using different random number engine.
     //! \param other Source object from which copy is made.
     //! \param random_number_engine Object from which random number generators can be constructed.
     TorusCsOmegaDBN(const TorusCsOmegaDBN &other,
              RandomNumberEngine *random_number_engine)
          : BackboneDBN<TorusCsOmegaDBN,TorusCsOmegaDBNNodes>(other, random_number_engine),
            settings(other.settings) {
     }
};
//! Name of TorusCsOmegaDBN model
const char TorusCsOmegaDBN::name[] = "torus-cs-omega";

}

#endif
