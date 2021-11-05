// gaussian.h --- Gaussian node. Density and sampler for a normal distribution
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

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "utils/random.h"
#include "utils/vector_matrix_3d.h"
#include "utils/vector_nd.h"
#include "parameters.h"
#include "node.h"

namespace phaistos {

//! Density class for the Gaussian distribution
class GaussianDensity {
private:
     //! standard deviation
     double sigma;

     //! mean
     double mu;

     //! Logarithm of the normalization constant
     double log_norm_const;

public:

     //! Constructor
     //! \param sigma Standard deviation
     //! \param mu Mean
     GaussianDensity(double sigma, double mu);

     //! Destructor
     ~GaussianDensity();

     //! Initialize
     //! \param sigma Standard deviation
     //! \param mu Mean
     void init(double sigma, double mu);

     //! Calculate log-likelihood
     //! \param x observed value
     //! \return log-likelihood value          
     inline double get_log_likelihood(double x) const {
          double dens = -(x-mu)*(x-mu)/(2*sigma*sigma) - log_norm_const;
          // assert(std::isfinite(dens));
          return dens;
     }

     //! Calculate likelihood
     //! \param x observed value
     //! \return likelihood value     
     inline double get_likelihood(double x) const {
          return std::exp(get_log_likelihood(x));
     }
};


//! Sampler class for the Gaussian distribution
class GaussianSampler {
private:

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::normal_distribution<> > variate_generator_gaussian;

     //! standard deviation
     double sigma;

     //! mean
     double mu;

public:
     //! Constructor
     //! \param sigma standard deviation
     //! \param mu Mean
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     GaussianSampler(double sigma, double mu,
                     RandomNumberEngine *random_number_engine = &random_global);     

     //! Copy constructor
     //! \param other Source object     
     GaussianSampler(const GaussianSampler &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     GaussianSampler(const GaussianSampler &other,
                     RandomNumberEngine *random_number_engine);

     //! Destructor
     ~GaussianSampler();

     //! Sample value
     //! \return sampled value          
     double sample();
};


//! Distribution class for the Gaussian distribution. Container of density and sampler objects
class GaussianDistribution {
public:

     //! Density object
     GaussianDensity *density;

     //! Sampler object
     GaussianSampler *sampler;
     
     //! Constructor
     //! \param sigma standard deviation
     //! \param mu Mean
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     GaussianDistribution(double sigma, double mu,
                          RandomNumberEngine *random_number_engine = &random_global);

     
     //! Copy constructor
     //! \param other Source object
     GaussianDistribution(const GaussianDistribution &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     GaussianDistribution(const GaussianDistribution &other,
                          RandomNumberEngine *random_number_engine);

     //! Destructor          
     ~GaussianDistribution();

     //! Sample value
     //! \return value     
     double sample();

     //! Calculate log likelihood
     //! \param x observed value
     //! \return log-likelihood value          
     double get_log_likelihood(double x);
};



//! Node class for the Gaussian distribution
template <typename DBN_TYPE, int NODE_INDEX, int PARENT_NODE=0>
class GaussianNode: public Node_1D<double,DBN_TYPE,PARENT_NODE> {
private:

     //! An array of distribution components
     GaussianDistribution **data;
     
public:

     //! Constructor
     //! \param name Node name
     //! \param parameters Parameter object
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     GaussianNode(std::string name, Parameters &parameters, 
                  DBN_TYPE *dbn,
                  RandomNumberEngine *random_number_engine = &random_global)
          : Node_1D<double,DBN_TYPE,PARENT_NODE>(name, dbn) {

	  // Extract parameters
	  
	  // kappa and mu are mandatory
	  bool raiseException = true;
	  std::vector<std::vector<double> > sigma_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("GAUSSIAN", name, "SIGMA", raiseException);
	  std::vector<std::vector<double> > mu_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("GAUSSIAN", name, "MU", raiseException);
     
          // std::cout << sigma_list << " " << sigma_list.size() << " " << sigma_list[0].size() << "\n";

	  this->h_size = sigma_list.size();
	  data = new GaussianDistribution*[this->h_size];

	  for (unsigned int i=0; i<(unsigned int)this->h_size; i++) {

	       // Get parameters for position i. 
	       
               // It seems that Mocapy reports the variance rather that stddev
               double sigma = sqrt(sigma_list[i][0]);
               // double sigma = sigma_list[i][0];

               double mu = mu_list[i][0];

	       // Initialize distribution component
	       data[i] = new GaussianDistribution(sigma, mu,
                                                  random_number_engine);
          }
     }

     
     //! Copy constructor
     //! \param other Source object
     GaussianNode(const GaussianNode &other)
          : Node_1D<double,DBN_TYPE>(other) {

	  data = new GaussianDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new GaussianDistribution(*other.data[i],
                                                  other->dbn->random_number_engine);
	  }
     }
     
     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     GaussianNode(const GaussianNode &other, DBN_TYPE *dbn)
          : Node_1D<double,DBN_TYPE>(other,dbn) {

	  data = new GaussianDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new GaussianDistribution(*other.data[i],
                                                  dbn->random_number_engine);
	  }
     }
     
     //! Destructor
     ~GaussianNode() {
	  for (int i=0; i<this->h_size; i++) {
	       delete data[i];
	  }
	  delete[] data;
     }
	  
     
     //! Overload indexing operator 
     //! \param index distribution component index
     //! \return distribution component
     GaussianDistribution *operator[](const int index) const {
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
		    int h = this->dbn->template get_node<typename DBN_TYPE::HIDDEN_NODE>()->sequence[l];
		    this->sequence[l] = data[h]->sampler->sample();
	       }
	  }
     }


     //! Get emission probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param x observation value
     //! \return likelihood
     inline double get_likelihood(int h, double x) const {
	  return data[h]->density->get_likelihood(x);     
     }

     //! Get emission log-probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param x observation value
     //! \return log-likelihood
     inline double get_log_likelihood(int h, double x) const {
	  return data[h]->density->get_log_likelihood(x);     
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

     
     //! Return sigma parameters for all distributions components
     //! \return sigma value for each component
     std::vector<double> get_shape_parameters() {
	  std::vector<double> shape_parameters;
	  for (int i=0; i<this->h_size; i++) {
	       shape_parameters.push_back(data[i]->density->sigma);               
	  }
	  return shape_parameters;
     }
};

}

#endif
