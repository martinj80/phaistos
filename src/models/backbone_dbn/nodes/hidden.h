// hidden.h --- Hidden node
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


#ifndef HIDDEN_H
#define HIDDEN_H

#include "discrete.h"
#include "node.h"

namespace phaistos {

// Type used for setting N-terminus and transition probabilities
enum TransitionEmissionState {NORMAL, UNIFORM, STATIONARY, TRANSITION_EMISSION_STATE_SIZE};
static const std::string TransitionEmissionStateNames[] = {"normal", "uniform", "stationary"};

//! Overload input operator for reading in TransitionEmissionState values.
inline std::istream &operator>>(std::istream &input, TransitionEmissionState &tes) {
     std::string raw_string;
     input >> raw_string;

     for (unsigned int i=0; i<TRANSITION_EMISSION_STATE_SIZE; ++i) {
          if (raw_string == TransitionEmissionStateNames[i]) {
               tes = TransitionEmissionState(i);
          }
     }
     return input;
}

//! Overload output operator for printing out TransitionEmissionState values.
inline std::ostream &operator<<(std::ostream &o, const TransitionEmissionState &tes) {
     o << TransitionEmissionStateNames[static_cast<unsigned int>(tes)];
     return o;
}


//! Hidden node class
template <typename DBN_TYPE, int NODE_INDEX=0>
class HiddenNode: public Node_1D<int,DBN_TYPE> {
private:

     //! Conditional probability distribution transition matrix
     DiscreteDistribution **cpd_transition;

     //! Conditional probability distribution initial distribution
     DiscreteDistribution *cpd_start;

public:

     //! Meta-function: Retrieve ParentNode value from NODE_TYPE - phase2.
     //! NOTE: The two-phase approach was introduced because filter_view seems
     //! to return an unexpanded boost::mpl::_ value on certain compilers. Phase1
     //! explicitly filters these out
     template <typename NODE_TYPE>
     struct GetParentNodePhase2 {
          // BOOST_MPL_ASSERT_MSG(boost::is_integral<double>::value, NON_INTEGRAL_TYPES_ARE_NOT_ALLOWED, (NODE_TYPE)); 
          //! return value (type)
          typedef typename NODE_TYPE::ParentNode type;
     };

     //! Meta-function: Retrieve ParentNode value from NODE_TYPE - phase1.
     //! NOTE: The two-phase approach was introduced because filter_view seems
     //! to return an unexpanded boost::mpl::_ value on certain compilers. Phase1
     //! explicitly filters these out
     template <typename NODE_TYPE, typename NUMBER>
     struct GetParentNodePhase1 {
          //! return value (type)
          typedef typename boost::mpl::eval_if<boost::is_same<boost::mpl::_, NUMBER>, boost::mpl::identity<NODE_TYPE>, GetParentNodePhase2<typename boost::remove_pointer<NODE_TYPE>::type > >::type type;          
     };

     //! Find indices of nodes which have the current node as parent
     typedef typename boost::mpl::filter_view<typename DBN_TYPE::ALL_CHILD_NODES, boost::is_same<typename boost::mpl::int_<NODE_INDEX>, GetParentNodePhase1<boost::mpl::at<typename DBN_TYPE::Nodes,boost::mpl::_>,boost::mpl::_> > >::type ChildNodeIndices;

     //! Child nodes type     
     typedef typename boost::fusion::nview<typename DBN_TYPE::Nodes,ChildNodeIndices> ChildNodes;

     //! Child nodes
     ChildNodes *child_nodes;

     //! Transition matrix (for faster access)
     double **transition_matrix; 
     
     //! Log-transition matrix (for faster access)     
     double **transition_log_matrix; 

     //! Constructor
     //! \param name Node name
     //! \param parameters Parameter object
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param initial_distribution Distribution at first slice (index 0)
     //! \param transition_distribution Type of transition distribution
     HiddenNode(std::string name, Parameters &parameters, DBN_TYPE *dbn,
		TransitionEmissionState initial_distribution=NORMAL,
		TransitionEmissionState transition_distribution=NORMAL)
          : Node_1D<int,DBN_TYPE>(name, dbn),
            child_nodes(NULL) {

	  // Extract transition parameters
	  std::vector<std::vector<double> > transition_parameters = parameters.get_node_parameters<std::vector<std::vector<double> > > ("DISCRETE", name, "TRANSITION");
	  
	  // Set uniform or stationary transitions if requested
	  if (transition_distribution == UNIFORM) {
	       for (unsigned int i=0; i<transition_parameters.size(); i++) {
		    for (unsigned int j=0; j<transition_parameters[i].size(); j++) {
			 transition_parameters[i][j] = 1/transition_parameters[i].size();
		    }
	       }
	  } else if (transition_distribution == STATIONARY) {
	       std::vector<double> stationary = parameters.get_node_parameters<std::vector<double> > ("DISCRETE", name, "DOMINANT-EIGENVECTOR");
	       for (unsigned int i=0; i<transition_parameters.size(); i++) {
		    for (unsigned int j=0; j<transition_parameters[i].size(); j++) {
			 transition_parameters[i][j] = stationary[j];		    
		    }
	       }
	  }

	  // Extract N-terminus probability parameters
	  std::vector<double> start_parameters = parameters.get_node_parameters<std::vector<double> > ("DISCRETE", name, "START");

	  // Set uniform or stationary N-terminus probability if requested
	  if (initial_distribution == UNIFORM) {
	       for (unsigned int i=0; i<start_parameters.size(); i++) {
		    start_parameters[i] = 1.0/start_parameters.size();
	       }
	  } else if (initial_distribution == STATIONARY) {
	       start_parameters = parameters.get_node_parameters<std::vector<double> > ("DISCRETE", name, "DOMINANT-EIGENVECTOR");
	  }

	  this->h_size = transition_parameters.size();
	  if (this->h_size>0)
	       this->size = transition_parameters[0].size();

	  
	  // Initialize components
	  cpd_transition = new DiscreteDistribution*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       cpd_transition[i] = new DiscreteDistribution(transition_parameters[i]);
	  }

	  cpd_start = new DiscreteDistribution(start_parameters);

	  // For fast access
	  transition_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       transition_matrix[i] = cpd_transition[i]->probability->probs;
	  }
	  transition_log_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       transition_log_matrix[i] = cpd_transition[i]->probability->log_probs;
	  }
     }
	  
     //! Copy constructor helper function
     //! \param other Source object
     void copy(const HiddenNode &other) {

	  cpd_transition = new DiscreteDistribution*[this->h_size];

	  for (int i=0; i<this->h_size; i++) {
	       cpd_transition[i] = new DiscreteDistribution(*other.cpd_transition[i]);
	  }

	  cpd_start = new DiscreteDistribution(*other.cpd_start);

	  // For fast access
	  transition_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       transition_matrix[i] = cpd_transition[i]->probability->probs;
	  }     
	  transition_log_matrix = new double*[this->h_size];
	  for (int i=0; i<this->h_size; i++) {
	       transition_log_matrix[i] = cpd_transition[i]->probability->log_probs;
	  }     
     }

     //! Copy constructor
     //! \param other Source object
     HiddenNode(const HiddenNode &other)
          : Node_1D<int,DBN_TYPE>(other),
            child_nodes(NULL) {
          copy(other);
     }

     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     HiddenNode(const HiddenNode &other, DBN_TYPE *dbn)
          : Node_1D<int,DBN_TYPE>(other,dbn),
            child_nodes(NULL) {
          copy(other);
     }

     //! Destructor
     ~HiddenNode() {
	  for (int i=0; i<this->h_size; i++) {
	       delete cpd_transition[i];
	  }
	  delete[] cpd_transition;

	  delete cpd_start;
	  delete[] transition_matrix;
	  delete[] transition_log_matrix;

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
	  return *cpd_transition[index];
     }

     //! Return N-terminus probability
     //! \param i hidden node value
     //! \return likelihood
     double start(int i) {
	  return cpd_start->probability->get_likelihood(i);
     }

     //! Return N-terminus log-probabilitiy
     //! \param i hidden node value
     //! \return log-likelihood
     double start_log(int i) {
	  return cpd_start->probability->get_log_likelihood(i);
     }

     //! Return transition probability from i to j
     //! \param i hidden node value
     //! \param j hidden node value
     //! \return likelihood     
     double transition(int i, int j) {
	  return cpd_transition[i]->probability->get_likelihood(j);

     }
     
     //! Return log transition probability from i to j
     //! \param i hidden node value
     //! \param j hidden node value
     //! \return log-likelihood     
     double transition_log(int i, int j) {
	  return cpd_transition[i]->probability->get_log_likelihood(j);

     }

     //! Return transition matrix
     //! \return pointer to transition matrix
     double **get_transition_matrix() {
	  return transition_matrix;
     }
     
     //! Return log transition matrix
     //! \return pointer to log-transition matrix
     double **get_transition_log_matrix() {
	  return transition_log_matrix;
     }


     //! Calculate emission probability
     //! \param h component index (typically hidden node value)
     //! \param pos Position in sequence
     //! \param use_backup Whether backup values (values from last accepted state) should be used
     //! \return likelihood
     inline double emission(int h, int pos, bool use_backup=false) {
          double val = 1.0;
          for_each(*child_nodes,
                   typename Node_1D<int,DBN_TYPE>::_get_emission(h, pos, val, use_backup));

          return val;
     }


     //! Calculate emission log-probability
     //! \param h Component index (typically hidden node value)
     //! \param pos Position in sequence
     //! \param use_backup Whether backup values (values from last accepted state) should be used
     //! \return log-likelihood
     inline double emission_log(int h, int pos, bool use_backup=false) {
          double val = 0.0;
          for_each(*child_nodes,
                   typename Node_1D<int,DBN_TYPE>::_get_emission_log(h, pos, val, use_backup));
          return val;
     }


     //! Resample hidden node sequence
     //! \param start Sequence start index
     //! \param end Sequence end index
     //! \param sampling Whether to do sampling or only calculate log-likelihood
     //! \return log-likelihood
     double sample(int start=-1, int end=-1, bool sampling=true) {

          double LL=0.0;
	  if (!this->fixed) {
	  
	       // Make sure that the complete hidden node sequence has previously been sampled
	       if (this->dbn->inconsistent_regions.size() > 0) {

		    // Reinitialize sequences
		    this->dbn->init_sequences();

		    // init_sequences causes everything to be resampled. No need to continue with this call
		    return 0.0; // Note: incorrect likelihood value
	       }

	       // Set up end points
	       if(end < 0 || end > this->sequence_length)
		    end = this->sequence_length;
     
	       if(start < 0 || start > this->sequence_length)
		    start = 0;

	       if (this->dbn->settings.log_space) {
		    this->dbn->forward_log(start, end);
		    LL = this->dbn->backtrack_log(start, end, sampling);
	       } else {
		    this->dbn->forward(start, end);
		    LL = this->dbn->backtrack(start, end, sampling);
               }
	  }
          return LL;
     }

};

}

#endif
