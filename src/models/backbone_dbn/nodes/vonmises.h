// vonmises.h --- von Mises node. Density and sampler for a univariate von Mises distritbution
// Copyright (C) 2006-2009 Wouter Boomsma
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

#ifndef VONMISES_H
#define VONMISES_H

#include "utils/random.h"
#include "utils/vector_matrix_3d.h"
#include "utils/vector_nd.h"
#include "parameters.h"
#include "node.h"

namespace phaistos {

//! Density class for univariate von Mises distribution
class VonMisesDensity {
private:
public:
     //! Kappa parameter
     double k;

     //! Mu parameter
     double mu;
     
     //! Logarithm of the normalization constant     
     double log_norm_const;

     //! Constructor
     //! \param k kappa parameter
     //! \param mu mu parameter
     VonMisesDensity(double k, double mu);

     //! Destructor
     ~VonMisesDensity();

     //! Initializer
     //! \param k kappa parameter
     //! \param mu mu parameter
     void init(double k, double mu);

     //! Calculate log-likelihood
     //! \param angle observed value
     //! \return log-likelihood value               
     double get_log_likelihood(double angle);

     //! Calculate likelihood
     //! \param angle observed value
     //! \return likelihood value          
     double get_likelihood(double angle);
};


//! Sampler class for univariate von Mises distribution
class VonMisesSampler {
private:

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > variate_generator_uniform;

     //! kappa parameter
     double k;

     //! mu parameter
     double mu;

     //@{
     //! Internals for sampling algorithm
     double a;
     double b;
     double r;
     //@}

public:
     //! Constructor
     //! \param k kappa parameter
     //! \param mu mu parameter
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     VonMisesSampler(double k, double mu,
                     RandomNumberEngine *random_number_engine = &random_global);     

     //! Copy constructor
     //! \param other Source object          
     VonMisesSampler(const VonMisesSampler &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     VonMisesSampler(const VonMisesSampler &other,
                     RandomNumberEngine *random_number_engine);

     //! Destructor     
     ~VonMisesSampler();

     //! Sample value
     //! \return sampled value
     double sample();
};


//! Distribution class for the von Mises distribution. Container of density and sampler objects
class VonMisesDistribution {
public:
     
     //! Density object
     VonMisesDensity *density;

     //! Sampler object     
     VonMisesSampler *sampler;
     
     //! Constructor
     //! \param k kappa parameter
     //! \param mu mu parameter
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     VonMisesDistribution(double k, double mu, 
                          RandomNumberEngine *random_number_engine = &random_global);

     //! Copy constructor
     //! \param other Source object     
     VonMisesDistribution(const VonMisesDistribution &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.          
     VonMisesDistribution(const VonMisesDistribution &other, 
                          RandomNumberEngine *random_number_engine);

     //! Destructor               
     ~VonMisesDistribution();

     //! Sample value
     //! \return value          
     double sample();

     //! Calculate log likelihood
     //! \param x observed value
     //! \return log-likelihood value               
     double get_log_likelihood(double x);
};



//! Node class for the von Mises distribution
template <typename DBN_TYPE, int NODE_INDEX, int PARENT_NODE=0 >
class VonMisesNode: public Node_1D<double,DBN_TYPE,PARENT_NODE> {
private:

     //! An array of distribution components
     VonMisesDistribution **data;
     
public:

     //! Constructor
     //! \param name Node name
     //! \param parameters Parameter object
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     VonMisesNode(std::string name, Parameters &parameters, 
                  DBN_TYPE *dbn,
                  RandomNumberEngine *random_number_engine = &random_global)
          : Node_1D<double,DBN_TYPE,PARENT_NODE>(name, dbn) {

	  // Extract parameters
	  
	  // kappa and mu are mandatory
	  bool raiseException = true;
	  std::vector<double> kappa_list = parameters.get_node_parameters<std::vector<double> > ("VONMISES", name, "KAPPA", raiseException);
	  std::vector<double> mu_list = parameters.get_node_parameters<std::vector<double> > ("VONMISES", name, "MU", raiseException);
     
	  this->h_size = kappa_list.size();
	  data = new VonMisesDistribution*[this->h_size];

	  for (unsigned int i=0; i<(unsigned int)this->h_size; i++) {

	       // Get parameters for position i. 
	       
	       double k = kappa_list[i];
	       double mu = mu_list[i];

	       // Initialize distribution component
	       data[i] = new VonMisesDistribution(k, mu,
                                                  random_number_engine);
          }
     }

     //! Constructor - without parameter file
     //! \param name Node name
     //! \param mu_list vector of mean values
     //! \param kappa_list vector of kappa values
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     VonMisesNode(std::string name, 
                  const std::vector<double> &mu_list, const std::vector<double> &kappa_list, 
                  DBN_TYPE *dbn,
                  RandomNumberEngine *random_number_engine = &random_global)
          : Node_1D<double,DBN_TYPE,PARENT_NODE>(name, dbn) {

	  this->h_size = kappa_list.size();
	  data = new VonMisesDistribution*[this->h_size];

	  for (unsigned int i=0; i<(unsigned int)this->h_size; i++) {

	       // Get parameters for position i. 
	       
	       double k = kappa_list[i];
	       double mu = mu_list[i];

	       // Initialize distribution component
	       data[i] = new VonMisesDistribution(k, mu,
                                                  random_number_engine);
          }
     }

     
     //! Copy constructor
     //! \param other Source object
     VonMisesNode(const VonMisesNode &other): Node_1D<double,DBN_TYPE,PARENT_NODE>(other) {
	  data = new VonMisesDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new VonMisesDistribution(*other.data[i],
                                                  other->dbn->random_number_engine);
	  }
     }

     
     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     VonMisesNode(const VonMisesNode &other, DBN_TYPE *dbn): Node_1D<double,DBN_TYPE,PARENT_NODE>(other,dbn) {
	  data = new VonMisesDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new VonMisesDistribution(*other.data[i],
                                                  dbn->random_number_engine);
	  }
     }
     

     //! Destructor
     ~VonMisesNode() {
	  for (int i=0; i<this->h_size; i++) {
	       delete data[i];
	  }
	  delete[] data;
     }
	  
     
     //! Overload indexing operator 
     //! \param index distribution component index
     //! \return distribution component
     VonMisesDistribution *operator[](const int index) const {
	  return data[index];
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

	       for(int l=start; l<end; l++) {
		    int h = this->dbn->template get_node<boost::mpl::int_<PARENT_NODE> >()->sequence[l];
                    if (h>=0 && h<this->h_size) {
                         this->sequence[l] = data[h]->sampler->sample();
                    } else {
                         this->sequence[l] = UNINITIALIZED;
                    }
	       }
	  }
     }


     //! Get emission probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param angle observation value
     //! \return likelihood
     inline double get_likelihood(int h, double angle) const {
	  return data[h]->density->get_likelihood(angle);     
     }

     //! Get emission log-probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param angle observation value
     //! \return log-likelihood
     inline double get_log_likelihood(int h, double angle) const {
	  return data[h]->density->get_log_likelihood(angle);     
     }
          
     //! Return mean values for all distribution components     
     //! \return mean value for each component
     std::vector<double> get_means() {
	  std::vector<double> means;
	  for (int i=0; i<this->h_size; i++) {
               means.push_back(data[i]->density->mu);
	  }
	  return means;
     }

     
     //! Return kappa parameters for all distributions components
     //! \return kappa value for each component
     std::vector<double> get_shape_parameters() {
	  std::vector<double> shape_parameters;
	  for (int i=0; i<this->h_size; i++) {
	       shape_parameters.push_back(data[i]->density->k);               
	  }
	  return shape_parameters;
     }
};

}

#endif
