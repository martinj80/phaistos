// torus.h --- Torus node. Density and sampler for a bivariate von Mises distritbution - cosine model
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

#ifndef TORUS_H
#define TORUS_H

#include "utils/random.h"
#include "utils/vector_matrix_3d.h"
#include "utils/vector_nd.h"
#include "parameters.h"
#include "node.h"

namespace phaistos {

class VM2CosMarginalSampler;

//! Density class for the bivariate von Mises cosine model distribution
class TorusDensity {
public:
     
     //! kappa parameters (k1,k2,k3)
     Vector_nD k;

     //! mu parameters (mu1,mu2)
     Vector_nD mu;

     //! logarithm of the normalization constant
     double log_norm_const;

     //! Constructor
     //! \param k kappa parameters (k1,k2,k3)
     //! \param mu mu parameters (mu1,mu2)
     TorusDensity(Vector_nD &k, Vector_nD &mu);

     //! Constructor - with log-norm-const provided.
     //! \param k Kappa parameters (k1,k2,k3)
     //! \param mu Mu parameters (mu1,mu2)     
     //! \param log_norm_const Logarithm of the normalization constant
     TorusDensity(Vector_nD &k, Vector_nD &mu, double log_norm_const);

     //! Destructor
     ~TorusDensity();

     //! Initialize node
     //! \param k kappa parameters (k1,k2,k3)
     //! \param mu mu parameters (mu1,mu2)
     void init(Vector_nD &k, Vector_nD &mu);

     //! Calculate log-likelihood
     //! \param angle_pair observed value
     //! \return log-likelihood value          
     inline double get_log_likelihood(double *angle_pair) const {
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          double mu1 = mu[0];
          double mu2 = mu[1];
          double phi = angle_pair[0];
          double psi = angle_pair[1];
          double dens;
          // if (logSpace) {
          dens = k1*cos(phi-mu1) + k2*cos(psi-mu2) - k3*cos(phi-mu1-psi+mu2) - log_norm_const;
          // } else {
          //      dens = std::exp(k1*cos(phi-mu1) + k2*cos(psi-mu2) - k3*cos(phi-mu1-psi+mu2))/std::exp(log_norm_const);
          // }
          return dens;
     }

     //! Calculate likelihood
     //! \param angle_pair observed value
     //! \return likelihood value          
     inline double get_likelihood(double *angle_pair) const {
          return std::exp(get_log_likelihood(angle_pair));
     }
};


//! Sampler class for the bivariate von Mises cosine model distribution
class TorusSampler {
private:

     //! Random number engine
     RandomNumberEngine *random_number_engine;

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > variate_generator_uniform;

     //! Kappa parameters (k1,k2,k3)
     Vector_nD k;

     //! Mu parameters (mu1,mu2)
     Vector_nD mu;

     //! Pointer to sampler of marginal distribution
     VM2CosMarginalSampler *marginal_phi_sampler;

public:
     //! Constructor
     //! \param k Kappa parameters (k1,k2,k3)
     //! \param mu Mu parameters (mu1,mu2)     
     //! \param log_norm_const Logarithm of the normalization constant
     //! \param envelope_scale The overall height of the envelope
     //! \param envelope_mu Envelope mu parameters
     //! \param envelope_k Envelope kappa parameters
     //! \param envelope_proportion Envelope proportion parameters
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     TorusSampler(Vector_nD &k, Vector_nD &mu,
                  double *log_norm_const, double *envelope_scale,
                  Vector_nD &envelope_k, Vector_nD &envelope_mu, Vector_nD &envelope_proportion,
                  RandomNumberEngine *random_number_engine = &random_global);

     //! Constructor
     //! \param k Kappa parameters (k1,k2,k3)
     //! \param mu Mu parameters (mu1,mu2)     
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     TorusSampler(Vector_nD &k, Vector_nD &mu,
                  RandomNumberEngine *random_number_engine = &random_global);

     //! Copy constructor
     //! \param other Source object          
     TorusSampler(const TorusSampler &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.          
     TorusSampler(const TorusSampler &other,
                  RandomNumberEngine *random_number_engine);

     //! Destructor
     ~TorusSampler();

     //! Sample value
     //! \return sampled value
     Vector_nD sample();
};


//! Distribution class for the bivariate von Mises cosine model distribution. Container of density and sampler objects
class TorusDistribution {
public:
     //! Density object     
     TorusDensity *density;

     //! Sampler object     
     TorusSampler *sampler;
     
     //! Constructor
     //! \param k Kappa parameters (k1,k2,k3)
     //! \param mu Mu parameters (mu1,mu2)     
     //! \param log_norm_const Logarithm of the normalization constant
     //! \param envelope_scale The overall height of the envelope
     //! \param envelope_mu Envelope mu parameters
     //! \param envelope_k Envelope kappa parameters
     //! \param envelope_proportion Envelope proportion parameters
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     TorusDistribution(Vector_nD &k, Vector_nD &mu,
                       double *log_norm_const, double *envelope_scale,
                       Vector_nD &envelope_k, Vector_nD &envelope_mu, Vector_nD &envelope_proportion,
                       RandomNumberEngine *random_number_engine = &random_global);

     //! Copy constructor
     //! \param other Source object     
     TorusDistribution(const TorusDistribution &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     TorusDistribution(const TorusDistribution &other,
                       RandomNumberEngine *random_number_engine);

     //! Destructor
     ~TorusDistribution();

     //! Sample value
     //! \return value     
     Vector_nD sample();

     //! Calculate log likelihood
     //! \param angle_pair observed value
     //! \return log-likelihood value          
     double get_log_likelihood(double *angle_pair);     
};



//! Node class for the bivariate von Mises cosine model distribution
template <typename DBN_TYPE, int NODE_INDEX, int PARENT_NODE=0>
class TorusNode: public Node_2D<double,DBN_TYPE,PARENT_NODE> {
private:

     //! An array of distribution components
     TorusDistribution **data;
     
public:

     //! \param name Node name
     //! \param parameters Parameter object
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     TorusNode(std::string name, Parameters &parameters, 
               DBN_TYPE *dbn,
               RandomNumberEngine *random_number_engine = &random_global)
          : Node_2D<double,DBN_TYPE,PARENT_NODE>(name, 2, dbn) {

	  // Extract parameters
	  
	  // kappa and mu are mandatory
	  bool raiseException = true;
	  std::vector<std::vector<double> > kappa_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "KAPPA", raiseException);
	  std::vector<std::vector<double> > mu_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "MU", raiseException);

	  // The rest are calculated if not already found in parameter file
	  raiseException = false;
	  std::vector<std::vector<double> > kappa_envelope_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "KAPPA-ENVELOPE", raiseException);
	  std::vector<std::vector<double> > envelope_mu_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "MU-ENVELOPE", raiseException);
	  std::vector<std::vector<double> > envelope_proportion_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "PROPORTION-ENVELOPE", raiseException);
	  std::vector<double> envelope_scale_list = parameters.get_node_parameters<std::vector<double> > ("TORUS", name, "SCALE-ENVELOPE", raiseException);
	  std::vector<double> log_norm_const_list = parameters.get_node_parameters<std::vector<double> > ("TORUS", name, "LOGNORMCONST", raiseException);
     

	  this->h_size = kappa_list.size();
	  data = new TorusDistribution*[this->h_size];

	  bool update_parameters = false;
	  for (unsigned int i=0; i<(unsigned int)this->h_size; i++) {

	       // Get parameters for position i. 
	       
	       Vector_nD k(3, kappa_list[i][0], kappa_list[i][1], kappa_list[i][2]);
	       Vector_nD mu(2, mu_list[i][0], mu_list[i][1]);


	       // If a vector was not initialized from the parameter file, its value is calculated by TorusDistribution
	       Vector_nD kappa_envelope;
	       if (kappa_envelope_list.size() > i) {
		    kappa_envelope = Vector_nD(kappa_envelope_list[i]);
	       }
	  
	       Vector_nD envelope_mu;
	       if (envelope_mu_list.size() > i) {
		    envelope_mu = Vector_nD(envelope_mu_list[i]);
	       }

	       Vector_nD envelope_proportion;
	       if (envelope_proportion_list.size() > i) {
		    envelope_proportion = Vector_nD(envelope_proportion_list[i]);
	       } 

	       double envelope_scale = 0.0;
	       if (envelope_scale_list.size() > i) {
		    envelope_scale = envelope_scale_list[i];
	       }
	  
	       double log_norm_const;
	       if (log_norm_const_list.size() > i) {
		    log_norm_const = log_norm_const_list[i];
	       }

	       // Initialize distribution component
	       data[i] = new TorusDistribution(k, mu,
                                               &log_norm_const, &envelope_scale,
                                               kappa_envelope, envelope_mu, envelope_proportion,
                                               random_number_engine);

	       // Update parameters if necessary
	       if (kappa_envelope_list.size() <= i) {
		    kappa_envelope_list.push_back(kappa_envelope.get_vector());
		    update_parameters = true;
	       }
	       if (envelope_mu_list.size() <= i) {
		    envelope_mu_list.push_back(envelope_mu.get_vector());
		    update_parameters = true;
	       }
	       if (envelope_proportion_list.size() <= i) {
		    envelope_proportion_list.push_back(envelope_proportion.get_vector());
		    update_parameters = true;
	       }
	       if (envelope_scale_list.size() <= i) {
		    envelope_scale_list.push_back(envelope_scale);
		    update_parameters = true;
	       }
	       if (log_norm_const_list.size() <= i) {
		    log_norm_const_list.push_back(log_norm_const);
		    update_parameters = true;
	       }
	  }

	  // Save updated parameters
	  if (update_parameters) {

	       parameters.set_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "KAPPA-ENVELOPE", kappa_envelope_list);
	       parameters.set_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "MU-ENVELOPE", envelope_mu_list);
	       parameters.set_node_parameters<std::vector<std::vector<double> > > ("TORUS", name, "PROPORTION-ENVELOPE", envelope_proportion_list);
	       parameters.set_node_parameters<std::vector<double> > ("TORUS", name, "SCALE-ENVELOPE", envelope_scale_list);
	       parameters.set_node_parameters<std::vector<double> > ("TORUS", name, "LOGNORMCONST", log_norm_const_list);

	       parameters.save_to_file();
	  }
     }

     
     //! Copy constructor
     //! \param other Source object
     TorusNode(const TorusNode &other)
          : Node_2D<double,DBN_TYPE>(other) {

	  data = new TorusDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new TorusDistribution(*other.data[i],
                                               other->dbn->random_number_engine);
	  }
     }
     
     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     TorusNode(const TorusNode &other, DBN_TYPE *dbn)
          : Node_2D<double,DBN_TYPE>(other, dbn) {
          
	  data = new TorusDistribution*[this->h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new TorusDistribution(*other.data[i],
                                               dbn->random_number_engine);
	  }
     }
     
     //! Destructor
     ~TorusNode() {
	  for (int i=0; i<this->h_size; i++) {
	       delete data[i];
	  }
	  delete[] data;
     }
	  
     
     //! Overload indexing operator 
     //! \param index distribution component index
     //! \return distribution component
     TorusDistribution *operator[](const int index) const {
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
		    Vector_nD p = data[h]->sampler->sample();

		    for (int i=0; i<p.size; i++) {
			 this->sequence[l][i] = p[i];
		    }
	       }
	  }
     }


     //! Get emission probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param angle_pair observation value
     //! \return likelihood
     inline double get_likelihood(int h, double *angle_pair) const {
	  return data[h]->density->get_likelihood(angle_pair);     
     }

     //! Get emission log-probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param angle_pair observation value
     //! \return log-likelihood
     inline double get_log_likelihood(int h, double *angle_pair) const {
	  return data[h]->density->get_log_likelihood(angle_pair);     
     }

     //! Get emission probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param angle1 first angle
     //! \param angle2 second angle
     //! \return likelihood
     // Return emission probability
     inline double get_likelihood(int h, double angle1, double angle2) {
	  double angle_pair[2];
	  angle_pair[0] = angle1;
	  angle_pair[1] = angle2;
	  return data[h]->density->get_likelihood(angle_pair);     
     }

     //! Get emission log-probability for a given observation value
     //! \param h Component index (typically hidden node value)
     //! \param angle1 first angle
     //! \param angle2 second angle
     //! \return log-likelihood
     inline double get_log_likelihood(int h, double angle1, double angle2) {
	  double angle_pair[2];
	  angle_pair[0] = angle1;
	  angle_pair[1] = angle2;
	  return data[h]->density->get_log_likelihood(angle_pair);     
     }
          
     //! Return mean values for all distribution components     
     //! \return vector of mean values for each component
     std::vector<std::vector<double> > get_means() {
	  std::vector<std::vector<double> > means;
	  for (int i=0; i<this->h_size; i++) {
	       std::vector<double> rowVec;
	       for (int j=0; j<this->size; j++) {
		    rowVec.push_back(data[i]->density->mu[j]);
	       }
	       means.push_back(rowVec);
	  }
	  return means;
     }

     
     //! Return kappa1, kappa2, kappa3 parameters for all distributions components
     //! \return kappa1,kappa2,kappa3 values for each component.
     std::vector<std::vector<double> > get_shape_parameters() {
	  int nParameters = 3;
	  std::vector<std::vector<double> > shape_parameters;
	  for (int i=0; i<this->h_size; i++) {
	       std::vector<double> rowVec;
	       for (int j=0; j<nParameters; j++) {
		    rowVec.push_back(data[i]->density->k[j]);
	       }
	       shape_parameters.push_back(rowVec);
	  }
	  return shape_parameters;
     }
};

}

#endif
