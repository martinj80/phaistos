// discrete.cpp --- Discrete node
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


#include "discrete.h"

namespace phaistos {

// Constructor
DiscreteProbability::DiscreteProbability(std::vector<double> &probs) {
     this->size = probs.size();
     this->probs = new double[size];
     this->log_probs = new double[size];
     for(int i=0; i<size; i++) {
	  this->probs[i] = probs[i];
	  this->log_probs[i] = log(probs[i]);
     }
}


// Copy constructor
DiscreteProbability::DiscreteProbability(const DiscreteProbability &other) {
     this->size = other.size;
     this->probs = new double[size];
     this->log_probs = new double[size];
     for(int i=0; i<size; i++) {
	  this->probs[i] = other.probs[i];
	  this->log_probs[i] = other.log_probs[i];
     }     
}


// Destructor
DiscreteProbability::~DiscreteProbability() {
     delete[] this->probs;
     delete[] this->log_probs;
}


/////////////
// Sampler //
/////////////

// Constructor
DiscreteSampler::DiscreteSampler(std::vector<double> &probs,
                                 RandomNumberEngine *random_number_engine)
     : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     this->size = probs.size();
     this->cumulative_probs = new double[size];
     calc_cumulative_probs(probs, this->cumulative_probs);
}


// Copy constructor
DiscreteSampler::DiscreteSampler(const DiscreteSampler &other): 
     variate_generator_uniform(other.variate_generator_uniform) {

     this->size = other.size;
     this->cumulative_probs = new double[size];
     for(int i=0; i<size; i++) {
	  this->cumulative_probs[i] = other.cumulative_probs[i];
     }     
}

// Copy constructor
DiscreteSampler::DiscreteSampler(const DiscreteSampler &other,
                                 RandomNumberEngine *random_number_engine)
     : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     this->size = other.size;
     this->cumulative_probs = new double[size];
     for(int i=0; i<size; i++) {
	  this->cumulative_probs[i] = other.cumulative_probs[i];
     }     
}


// Destructor
DiscreteSampler::~DiscreteSampler() {
     delete[] cumulative_probs;
}

// STATIC METHOD Internal sampling method
int DiscreteSampler::_sample(int size, double *cumulative_probs, 
                             RandomNumberEngine *random_number_engine) {
     boost::variate_generator<RandomNumberEngine&, boost::uniform_real<> > 
          variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1));
     double rnd = variate_generator_uniform();
     
     for(int i=0; i<size; i++) {
          if(rnd<cumulative_probs[i]) {
               return i;
          }
     }
     return size-1;
}


// STATIC METHOD Internal sampling method
int DiscreteSampler::_sample(int size, double *cumulative_probs, 
                             boost::variate_generator<RandomNumberEngine&, boost::uniform_real<> > 
                             *variate_generator_uniform) {
     double rnd = (*variate_generator_uniform)();
     
     for(int i=0; i<size; i++) {
          if(rnd<cumulative_probs[i]) {
               return i;
          }
     }
     return size-1;
}


// Return next sample from the discrete distribution 
int DiscreteSampler::sample() {
     return DiscreteSampler::_sample(this->size, this->cumulative_probs, &this->variate_generator_uniform);
}

// STATIC METHOD: Draw a sample from a given probs array
int DiscreteSampler::sample(int size, double *probs, double *normalization, 
                            RandomNumberEngine *random_number_engine) {
     bool log_space = false;
     double *cumulative_probs = new double[size];
     DiscreteSampler::calc_cumulative_probs(size, probs, cumulative_probs, log_space, normalization);     
     int val = DiscreteSampler::_sample(size, cumulative_probs, random_number_engine);
     delete[] cumulative_probs;
     return val;
}

// STATIC METHOD: Draw a sample from a given probs vector
int DiscreteSampler::sample(std::vector<double> &probs, double *normalization, 
                            RandomNumberEngine *random_number_engine) {
     bool log_space = false;
     double *cumulative_probs = new double[probs.size()];
     DiscreteSampler::calc_cumulative_probs(probs, cumulative_probs, log_space, normalization);     
     int val = DiscreteSampler::_sample(probs.size(), cumulative_probs, random_number_engine);
     delete[] cumulative_probs;
     return val;
}

// STATIC METHOD: Draw a sample from a given log-probs array
int DiscreteSampler::sample_log(int size, double *log_probs, double *normalization, 
                                RandomNumberEngine *random_number_engine) {
     bool log_space = true;
     double *cumulative_probs = new double[size];
     DiscreteSampler::calc_cumulative_probs(size, log_probs, cumulative_probs, log_space, normalization);     
     int val = DiscreteSampler::_sample(size, cumulative_probs, random_number_engine);
     delete[] cumulative_probs;
     return val;
}

// STATIC METHOD: Draw a sample from a given log-probs vector
int DiscreteSampler::sample_log(std::vector<double> &log_probs, double *normalization, 
                                RandomNumberEngine *random_number_engine) {
     bool log_space = true;
     double *cumulative_probs = new double[log_probs.size()];
     DiscreteSampler::calc_cumulative_probs(log_probs, cumulative_probs, log_space, normalization);     
     int val = DiscreteSampler::_sample(log_probs.size(), cumulative_probs, random_number_engine);
     delete[] cumulative_probs;
     return val;
}

// STATIC METHOD: Calculate cummulative probability mass
void DiscreteSampler::calc_cumulative_probs(int size, double *probs, double *cumulative_probs, bool log_space, double *normalization) {
     if (size==0)
	  return;

     if (!log_space) {
	  cumulative_probs[0] = probs[0];
	  for(int i=1; i<size; i++) {
	       cumulative_probs[i] = cumulative_probs[i-1] + probs[i];
	  }
     } else {
	  cumulative_probs[0] = exp(probs[0]);
	  for(int i=1; i<size; i++) {
	       cumulative_probs[i] = cumulative_probs[i-1] + exp(probs[i]);
	  }
     }

     double normalization_inner;
     if (!normalization)
          normalization = &normalization_inner;
     *normalization = cumulative_probs[size-1];

     // Fix possible rounding errors
     for(int i=0; i<size; i++) {
          cumulative_probs[i]/=*normalization;
     }     
}

// STATIC METHOD: Calculate cummulative probability mass
void DiscreteSampler::calc_cumulative_probs(std::vector<double> &probs, double *cumulative_probs, bool log_space, double *normalization) {
     return calc_cumulative_probs(probs.size(), &probs[0], cumulative_probs, log_space, normalization);
}




//////////////////
// Distribution //
//////////////////

// Constructor - Initialize sampler and density
DiscreteDistribution::DiscreteDistribution(std::vector<double> &probs,
                                           RandomNumberEngine *random_number_engine) {
     this->probability = new DiscreteProbability(probs);
     this->sampler = new DiscreteSampler(probs, random_number_engine);
}


// Copy constructor
DiscreteDistribution::DiscreteDistribution(const DiscreteDistribution &other) {
     this->probability = new DiscreteProbability(*other.probability);
     this->sampler = new DiscreteSampler(*other.sampler);
}


// Copy constructor
DiscreteDistribution::DiscreteDistribution(const DiscreteDistribution &other,
                                           RandomNumberEngine *random_number_engine) {
     this->probability = new DiscreteProbability(*other.probability);
     this->sampler = new DiscreteSampler(*other.sampler, random_number_engine);
}


// Destructor
DiscreteDistribution::~DiscreteDistribution() {
     delete probability;
     delete sampler;
}

// Draw sample
int DiscreteDistribution::sample() {
     return sampler->sample();
}

// Evaluate log-likelihood
double DiscreteDistribution::get_log_likelihood(int x) {
     return probability->get_log_likelihood(x);
}

}
