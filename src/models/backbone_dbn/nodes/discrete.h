// discrete.h --- Discrete node
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


#ifndef DISCRETE_H
#define DISCRETE_H

#include "utils/random.h"
#include "parameters.h"
#include "node.h"

namespace phaistos {

// Logarithm of sum of probabilities - calculated from log-probs
inline double ln_sum(double lnx, double lny) {
     double val;
     if (lnx > lny) {
	  val = lnx + std::log(1+std::exp(lny-lnx));
     } else {
	  if (std::isinf(lny))
	  // if (fabs(std::isinf(lny)) == 1)
	       val = lny;
	  else
	       val = lny + std::log(1+std::exp(lnx-lny));
     }

     // strong_assert(!std::isnan(val));
     return val;
}

//! Discrete Probability Distribution class
class DiscreteProbability {
private:

     //! Size of outcome space 
     int size;

public:
     //! Array of probabilities (weights)
     double *probs;

     //! Array of log-probabilities
     double *log_probs;

     //! Constructor 
     //! \param probs vector of probabilities
     DiscreteProbability(std::vector<double> &probs);     

     //! Copy constructor
     //! \param other Source object
     DiscreteProbability(const DiscreteProbability &other);

     //! Destructor
     ~DiscreteProbability();

     //! Calculate likelihood
     //! \param x observation value
     //! \return likelihood value
     inline double get_likelihood(int x) const {
          if (x>=0)
               return probs[x];
          else
               return std::numeric_limits<double>::quiet_NaN();
     }

     //! Calculate log-likelihood
     //! \param x observation value
     //! \return log-likelihood value
     inline double get_log_likelihood(int x) const {
          if (x>=0)
               return log_probs[x];
          else 
               return std::numeric_limits<double>::quiet_NaN();
     }
};


//! Discrete Probability Distribution - sampler class
class DiscreteSampler {
private:
     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > variate_generator_uniform;

     //! Size of outcome space 
     int size;

     //! Array of cumulative probabilities
     double *cumulative_probs;
public:

     //! Constructor
     //! \param probs Vector of probabilities
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     DiscreteSampler(std::vector<double> &probs,
                     RandomNumberEngine *random_number_engine = &random_global);     
     
     //! Copy constructor
     //! \param other Source object
     DiscreteSampler(const DiscreteSampler &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     DiscreteSampler(const DiscreteSampler &other,
                     RandomNumberEngine *random_number_engine);

     //! Destructor
     ~DiscreteSampler();

     //! Sample value given probability array
     //! \param size Size of input array
     //! \param probs Input probability array
     //! \param normalization Optional pointer with which normalization constant can be extracted
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     //! return sampled value
     static int sample(int size, double *probs, double *normalization=NULL,
                       RandomNumberEngine *random_number_engine= &random_global);

     //! Sample value given log-probability array
     //! \param size Size of input array
     //! \param log_probs Input probability array
     //! \param normalization Optional pointer with which normalization constant can be extracted
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     //! return sampled value
     static int sample_log(int size, double *log_probs, double *normalization=NULL,
                           RandomNumberEngine *random_number_engine= &random_global);

     //! Sample value vector of probabilities
     //! \param probs Input probability vector
     //! \param normalization Optional pointer with which normalization constant can be extracted
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     //! return sampled value
     static int sample(std::vector<double> &probs, double *normalization=NULL,
                       RandomNumberEngine *random_number_engine= &random_global);


     //! Sample value vector of log-probabilities
     //! \param log_probs Input probability vector
     //! \param normalization Optional pointer with which normalization constant can be extracted
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     //! return sampled value
     static int sample_log(std::vector<double> &log_probs, double *normalization=NULL,
                           RandomNumberEngine *random_number_engine= &random_global);

     //! Internal method for sampling given a cumulative probability array
     //! \param size Size of input array
     //! \param cumulative_probs Input cumulative probability array
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     //! return sampled value     
     static int _sample(int size, double *cumulative_probs, 
                        RandomNumberEngine *random_number_engine= &random_global);

     //! Internal method for sampling given a cumulative probability array and a specified random number generator.
     //! \param size Size of input array
     //! \param cumulative_probs Input cumulative probability array
     //! \param variate_generator_uniform Uniform random number generator.
     //! return sampled value     
     static int _sample(int size, double *cumulative_probs, 
                        boost::variate_generator<RandomNumberEngine&, boost::uniform_real<> >
                        *variate_generator_uniform);

     //! Calculate cumulative probabilities
     //! \param size Size of input array
     //! \param probs Input probability array
     //! \param cumulative_probs Output array for cumulative probabilities
     //! \param log_space Whether to return log-probs instead of probs
     //! \param normalization Optional pointer with which normalization constant can be extracted
     static void calc_cumulative_probs(int size, double *probs, double *cumulative_probs, bool log_space=false, double *normalization=NULL);

     //! Calculate cumulative probabilities
     //! \param probs Input probability vector
     //! \param cumulative_probs Output array for cumulative probabilities
     //! \param log_space Whether to return log-probs instead of probs
     //! \param normalization Optional pointer with which normalization constant can be extracted
     static void calc_cumulative_probs(std::vector<double> &probs, double *cumulative_probs, bool log_space=false, double *normalization=NULL);

     //! Sample value from sampler object
     //! \return sampled value
     int sample();
};

     
//! Discrete Distribution class - container of density and sampler objects
class DiscreteDistribution {
public:
     //! Probability object
     DiscreteProbability *probability;

     //! Sampler object
     DiscreteSampler *sampler;
     
     //! Constructor
     //! \param probs Input probability vector
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     DiscreteDistribution(std::vector<double> &probs,
                          RandomNumberEngine *random_number_engine = &random_global);

     //! Copy constructor
     //! \param other Source object
     DiscreteDistribution(const DiscreteDistribution &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     DiscreteDistribution(const DiscreteDistribution &other,
                          RandomNumberEngine *random_number_engine);

     //! Destructor
     ~DiscreteDistribution();

     //! Sample value
     //! \return value
     int sample();

     //! Calculate log likelihood
     //! \param x observed value
     //! \return log-likelihood value
     double get_log_likelihood(int x);     
};


//! Discrete node class
template <typename DBN_TYPE, int NODE_INDEX, int PARENT_NODE=0>
class DiscreteNode: public Node_1D<int,DBN_TYPE,PARENT_NODE> {
private:

     //! An array of distribution objects
     DiscreteDistribution **data;

     //! Local copy of emission matrix (for faster access)
     double **emission_matrix;

     //! Local copy of log-emissions (for faster access)
     double **emission_log_matrix; 


     //! Internal functor used for sampling (used by sample())
     struct _sample {
          
          //! Start index
          int start; 

          //! End index
          int end;

          //! Constructor
          _sample(int start, int end)
               :start(start),end(end){}

          //! Overload () operator. Sample values for node.
          //! \param node DBN node to be initialized.
          template <typename TYPE>
          void operator()(TYPE *node) const {
               node->sample(start, end);               
          }
     };
          
public:

     //! Meta-function: Retrieve ParentNode value from NODE_TYPE - phase2.
     //! NOTE: The two-phase approach was introduced because filter_view seems
     //! to return an unexpanded boost::mpl::_ value on certain compilers. Phase1
     //! explicitly filters these out
     template <typename NODE_TYPE>
     struct GetParentNodePhase2 {
          //! return value
          typedef typename NODE_TYPE::ParentNode type;
     };

     //! Meta-function: Retrieve ParentNode value from NODE_TYPE - phase1.
     //! NOTE: The two-phase approach was introduced because filter_view seems
     //! to return an unexpanded boost::mpl::_ value on certain compilers. Phase1
     //! explicitly filters these out
     template <typename NODE_TYPE, typename NUMBER>
     struct GetParentNodePhase1 {
          //! return value
          typedef typename boost::mpl::eval_if<boost::is_same<boost::mpl::_, NUMBER>, boost::mpl::identity<NODE_TYPE>, GetParentNodePhase2<typename boost::remove_pointer<NODE_TYPE>::type > >::type type;
     };

     //! Find indices of nodes which have the current node as parent
     typedef typename boost::mpl::filter_view<typename DBN_TYPE::ALL_CHILD_NODES, boost::is_same<typename boost::mpl::int_<NODE_INDEX>, GetParentNodePhase1<boost::mpl::at<typename DBN_TYPE::Nodes,boost::mpl::_>,boost::mpl::_> > >::type ChildNodeIndices;

     //! Child nodes type
     typedef typename boost::fusion::nview<typename DBN_TYPE::Nodes,ChildNodeIndices> ChildNodes;

     //! Child nodes
     ChildNodes *child_nodes;


     //! Constructor
     //! \param name Node name
     //! \param parameters Parameter object
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     DiscreteNode(std::string name, Parameters &parameters, 
                  DBN_TYPE *dbn,
                  RandomNumberEngine *random_number_engine = &random_global)
          : Node_1D<int,DBN_TYPE>(name, dbn),
            child_nodes(NULL) {

	  // Extract parameters
	  std::vector<std::vector<double> > emission = parameters.get_node_parameters<std::vector<std::vector<double> > > ("DISCRETE", name, "EMISSION");

	  this->h_size = emission.size();
	  if (this->h_size>0)
	       this->size = emission[0].size();
     
	  data = new DiscreteDistribution*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new DiscreteDistribution(emission[i], random_number_engine);
	  }

	  // For fast access
	  emission_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       emission_matrix[i] = data[i]->probability->probs;
	  }
	  emission_log_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       emission_log_matrix[i] = data[i]->probability->log_probs;
	  }	  
     }

     //! Copy constructor helper
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     void copy(const DiscreteNode &other,
               RandomNumberEngine *random_number_engine) {

	  data = new DiscreteDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new DiscreteDistribution(*other.data[i], 
                                                  random_number_engine);
	  }

	  // For fast access
	  emission_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       emission_matrix[i] = data[i]->probability->probs;
	  }     
	  emission_log_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       emission_log_matrix[i] = data[i]->probability->log_probs;
	  }     	  
     }

     //! Copy constructor
     //! \param other Source object
     DiscreteNode(const DiscreteNode &other)
          : Node_1D<int,DBN_TYPE>(other),
            child_nodes(NULL) {
          copy(other, other->dbn->random_number_engine);
     }

     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     DiscreteNode(const DiscreteNode &other, DBN_TYPE *dbn)
          : Node_1D<int,DBN_TYPE>(other,dbn),
            child_nodes(NULL) {
          copy(other, dbn->random_number_engine);
     }

     //! Destructor
     ~DiscreteNode() {
	  for (int i=0; i<this->h_size; i++) {
	       delete data[i];
	  }
	  delete[] data;

	  delete[] emission_matrix;
	  delete[] emission_log_matrix;	  

          delete child_nodes;
     }

     //! Initialize node
     //! \param length length of sequence
     //! \param complete_reinitialization Whether to reset all sequences
     void init(int length, bool complete_reinitialization=true) {
          Node_1D<int,DBN_TYPE>::init(length, complete_reinitialization);
          child_nodes = new ChildNodes(this->dbn->nodes);
     }     

     //! Overload indexing operator 
     //! \param index distribution component index
     //! \return distribution component
     DiscreteDistribution &operator[](const int index) const {
	  return *data[index];
     }

     //! Calculate emission probability
     //! \param h component index (typically hidden node value)
     //! \param pos Position in sequence
     //! \param use_backup Whether backup values (values from last accepted state) should be used
     //! \return likelihood
     template <typename DERIVED_CLASS>
     inline double emission(int h, int pos, bool use_backup) const {
          double result = 1.0;

          if (this->fixed && this->sequence_observed[pos]) {

               int node_value; 
               if (use_backup) {
                    node_value = this->sequence_backup[pos];
               } else {
                    node_value = this->sequence[pos];
               }

               result = data[h]->probability->get_likelihood(node_value);

               if (!boost::fusion::empty(*child_nodes)) {
                    for_each(*child_nodes,
                             typename Node_1D<int,DBN_TYPE,PARENT_NODE>::_get_emission(node_value, pos, result, use_backup));
               }

          } else if (!boost::fusion::empty(*child_nodes)) {

               bool any_children_fixed = true;
               for_each(*child_nodes,
                        typename Node<DBN_TYPE>::_get_emission_state(any_children_fixed));

               // If current node is not fixed, and unfixed children exist, sum out value
               if (any_children_fixed) {                        
                    result = 0.0;
                    for (int i=0; i<this->size; ++i) {
                         int node_value = i;
                         double partial_result = 1.0;
                         partial_result = data[h]->probability->get_likelihood(node_value);
                         for_each(*child_nodes,
                                  typename Node_1D<int,DBN_TYPE,PARENT_NODE>::_get_emission(i, pos, partial_result, use_backup));                    
                         result += partial_result;
                    }
               }
          }
          return result;
     }


     //! Calculate emission log-probability
     //! \param h Component index (typically hidden node value)
     //! \param pos Position in sequence
     //! \param use_backup Whether backup values (values from last accepted state) should be used
     //! \return log-likelihood
     template <typename DERIVED_CLASS>
     inline double emission_log(int h, int pos, bool use_backup) const {
          double result = 0.0;

          if (this->fixed && this->sequence_observed[pos]) {

               int node_value; 
               if (use_backup) {
                    node_value = this->sequence_backup[pos];
               } else {
                    node_value = this->sequence[pos];
               }

               result = data[h]->probability->get_log_likelihood(node_value);

               if (!boost::fusion::empty(*child_nodes)) {
                    for_each(*child_nodes,
                             typename Node_1D<int,DBN_TYPE,PARENT_NODE>::_get_emission_log(node_value, pos, result, use_backup));
               }

          } else if (!boost::fusion::empty(*child_nodes)) {

               bool any_children_fixed = true;
               for_each(*child_nodes,
                        typename Node<DBN_TYPE>::_get_emission_state(any_children_fixed));
               
               // If current node is not fixed, and unfixed children exist, sum out value
               if (any_children_fixed) {                        
                    result = std::log(0.0);
                    for (int i=0; i<this->size; ++i) {
                         int node_value = i;
                         double partial_result = data[h]->probability->get_log_likelihood(node_value);
                         for_each(*child_nodes,
                                  typename Node_1D<int,DBN_TYPE,PARENT_NODE>::_get_emission_log(i, pos, partial_result, use_backup));                    
                         result = ln_sum(result,partial_result);
                    }
               }
          }
          return result;
     }


     //! Get emission probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param i observation value
     //! \return likelihood
     inline double get_likelihood(int h, int i) const {
          return data[h]->probability->get_likelihood(i);
     }
     
     //! Get emission log-probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param i observation value
     //! \return log-likelihood
     inline double get_log_likelihood(int h, int i) const {
          return data[h]->probability->get_log_likelihood(i);
     }

     //! Resample values in sequence from start to end 
     //! \param start Sequence start index
     //! \param end Sequence end index
     void sample(int start=-1, int end=-1) {

	  if (!this->fixed) {
	  
	       // Make sure that hidden node sequence is initialized
               if (this->dbn->inconsistent_regions.size() > 0) {

		    // Reinitialize sequences
		    this->dbn->init_sequences();
		    return;
	       }

	       if(end < 0 || end > this->sequence_length)
		    end = this->sequence_length;
     
	       if(start < 0 || start > this->sequence_length)
		    start = 0;

               bool any_children_fixed = true;
               for_each(*child_nodes,
                        typename Node_1D<int,DBN_TYPE,PARENT_NODE>::_get_emission_state(any_children_fixed));

	       for(int l=start; l<end; l++) {

		    int h = this->dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->sequence[l];

                    // If there are no childnodes (or none are fixed), simply sample from node
                    if (boost::fusion::empty(*child_nodes) || !any_children_fixed) {                        
                         this->sequence[l] = data[h]->sampler->sample();
                    } else {

                         // Evaluate emission for each output value
                         std::vector<double> probs(this->size);
                         for (int i=0; i<this->size; ++i) {
                              int node_value = i;
                              probs[i] = data[h]->probability->get_likelihood(node_value);
                              boost::fusion::for_each(*child_nodes, 
                                                      typename Node_1D<int,DBN_TYPE,PARENT_NODE>::_get_emission(i,l,probs[i]));
                         }
                         
                         this->sequence[l] = DiscreteSampler::sample(probs, NULL, this->dbn->random_number_engine);
                    }
	       }
	  }

          // Sample child nodes
          boost::fusion::for_each(*child_nodes, _sample(start, end));
     }


     //! Return emission matrix
     //! \return pointer to matrix
     double **get_emission_matrix() {
     	  return emission_matrix;
     }
     
     //! Return log emission matrix
     //! \return pointer to matrix
     double **get_emission_log_matrix() {
     	  return emission_log_matrix;
     }

     //! Return the distribution corresponding to the viterbi hidden node path.
     //! \param start Sequence start index
     //! \param end Sequence end index
     //! \return Distribution over output values corresponding to the viterbi path over component values
     std::vector<std::vector<double> > get_viterbi_distribution(int start=-1, int end=-1) {

	  if(end < 0 || end > this->sequence_length)
	       end = this->sequence_length;
     
	  if(start < 0 || start > this->sequence_length)
	       start = 0;

	  std::vector<int> viterbi_path = this->dbn->viterbi(start, end);
     
	  std::vector<std::vector<double> > distribution;
	  for (unsigned int i=0; i<viterbi_path.size(); i++) {
	       std::vector<double> positionVec;
	       for (int j=0; j<this->size; j++) {
		    int h = viterbi_path[i];
		    positionVec.push_back(get_likelihood(h,j));
	       }
	       distribution.push_back(positionVec);
	  }
	  return distribution;     
     }


     //! Return the posterior distribution for each position.
     //! \param start Sequence start index
     //! \param end Sequence end index
     //! \return Posterior distribution
     std::vector<std::vector<double> > get_posterior_distribution(int start=-1, int end=-1) {

	  if(end < 0 || end > this->sequence_length)
	       end = this->sequence_length;
     
	  if(start < 0 || start > this->sequence_length)
	       start = 0;
     
	  const std::vector<std::vector<double> > &posterior_hidden = this->dbn->posterior(start, end);
	  std::vector<std::vector<double> > distribution;

	  for (unsigned int i=0; i<posterior_hidden.size(); i++) {
	       std::vector<double> positionVec;
	       for (int j=0; j<this->size; j++) {
		    double sum = 0.0;
		    for (int h=0; h<this->h_size; h++) {
			 sum += posterior_hidden[i][h] * get_likelihood(h, j);
		    }
		    positionVec.push_back(sum);
	       }
	       distribution.push_back(positionVec);
	  }

	  return distribution;     
     }


private:

     //! Hypothesis class used by one_best_decoding
     class Hypothesis {
     public:
          //! Internal data object
          std::vector<int> data;

          //! Keeps track of the probability that each hidden node assigns to this hypothesis
          std::vector<double> hidden_node_ll;

          //! Whether this hypothesis has been selected by any hidden node
          bool selected;

          //! Constructor
          Hypothesis(const int hidden_node_size)
               : hidden_node_ll(hidden_node_size),
                 selected(false) {}

          //! Constructor
          Hypothesis(const int hidden_node_size, 
                     const std::vector<int> &data)
               : data(data),
                 hidden_node_ll(hidden_node_size),
                 selected(false) {
          }
     };

public:

     //! Return the best sequence according to the "one best decoding" algorithm
     //! \return One best decoding of labels
     std::vector<int> one_best_decoding_log() {

          // This method has only been tested for direct child nodes to the hidden node
          BOOST_MPL_ASSERT((boost::is_same<boost::mpl::int_<PARENT_NODE>,
                                           typename DBN_TYPE::HIDDEN_NODE>));
          

          // The node for which prediction is done should be free
          // (this means the get_emission will not include this node)
          if (this->fixed) {
               std::cerr << "The node on which one-best decoding is done should be free (non-input)\n";
               return std::vector<int>();
          }

          // Get transition matrix for fast access
          // The alternative: hiddenNode->transition(i,j) is significantly slower
          // double **cpd_transition = this->dbn->hiddenNode->get_transition_matrix();
          double **cpd_transition_log = 
               this->dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->get_transition_log_matrix();


          // Vector of hypotheses
          std::vector<Hypothesis> hypotheses(this->size, this->dbn->hidden_node_size);
          for (unsigned int i=0; i<hypotheses.size(); ++i) {
               hypotheses[i].data = std::vector<int>(1,i);
          }
          std::vector<Hypothesis> hypotheses_old;

          for (unsigned int i=0; i<hypotheses.size(); ++i) {
          
               for (int h=0; h<this->dbn->hidden_node_size; ++h) {
                    hypotheses[i].hidden_node_ll[h] = 
                         this->dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->start_log(h) +
                         get_log_likelihood(h, i)+this->dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->emission_log(h, 0);          
               }
          }


          for (int pos=1; pos<this->dbn->sequence_length; ++pos) {

               // Store old hypotheses
               hypotheses_old = hypotheses;

               hypotheses.clear();

               // Create all new possible hypotheses
               for (unsigned int i=0; i<hypotheses_old.size(); ++i) {
                    for (int node_output=0; node_output<this->size; ++node_output) {

                         hypotheses.push_back(Hypothesis(this->dbn->hidden_node_size, 
                                                         hypotheses_old[i].data));
                         hypotheses.back().data.push_back(node_output);

                         for (int h=0; h<this->dbn->hidden_node_size; ++h) {

                              double transition_ll = std::log(0.0);
                              for (int h_prev=0; h_prev<this->dbn->hidden_node_size; ++h_prev) {
                                   transition_ll = ln_sum(transition_ll, cpd_transition_log[h_prev][h] + hypotheses_old[i].hidden_node_ll[h_prev]);
                              }
                              transition_ll += this->dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->emission_log(h, pos) + get_log_likelihood(h, node_output);

                              hypotheses.back().hidden_node_ll[h] = transition_ll;
                         }
                    }
               }

               for (int h=0; h<this->dbn->hidden_node_size; ++h) {
                    int max_index = 0;
                    double max_value = -std::numeric_limits<double>::infinity();
                    for (unsigned int i=0; i<hypotheses.size(); ++i) {
                         double ll = hypotheses[i].hidden_node_ll[h];
                         if (ll > max_value) {
                              max_value = ll;
                              max_index = i;
                         }
                    }
                    hypotheses[max_index].selected = true;
               }

               // Remove all hypotheses not chosen by any of the hidden nodes
               for (typename std::vector<Hypothesis>::iterator it = hypotheses.begin(); it !=hypotheses.end();) {
                    if (!((*it).selected)) {
                         it = hypotheses.erase(it);
                    } else {
                         ++it;
                    }
               }
          }

          int max_index = 0;
          double max_value = -std::numeric_limits<double>::infinity();
          for (unsigned int i=0; i<hypotheses.size(); ++i) {
               double tot_ll = 0.0;
               for (int h=0; h<this->dbn->hidden_node_size; ++h) {
                    tot_ll += hypotheses[i].hidden_node_ll[h];
               }

               if (tot_ll > max_value) {
                    max_index = i;
                    max_value = tot_ll;
               }
          }

          return hypotheses[max_index].data;
     }


};

}

#endif
