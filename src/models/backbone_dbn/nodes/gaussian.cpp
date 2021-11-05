// gaussian.cpp --- Gaussian node. Density and sampler for a normal distribution
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

#include "gaussian.h"

#include "utils/vector_nd.h"

#include <iostream>
#include <math.h>

namespace phaistos {

/*****************/
/**** DENSITY ****/
/*****************/

// Constructor
GaussianDensity::GaussianDensity(double sigma, double mu) {
     this->sigma = sigma;
     this->mu = mu;
     this->log_norm_const = 0.5*std::log(2*M_PI*sigma*sigma);
}

// Destructor
GaussianDensity::~GaussianDensity() {}



/*****************/
/**** SAMPLER ****/
/*****************/

// Sample from a Gaussian distribution
GaussianSampler::GaussianSampler(double sigma, double mu,
                                 RandomNumberEngine *random_number_engine)
     : variate_generator_gaussian(*random_number_engine, boost::normal_distribution<>(mu, sigma)),
       sigma(sigma), mu(mu) {}

// Copy constructor
GaussianSampler::GaussianSampler(const GaussianSampler &other)
     : variate_generator_gaussian(other.variate_generator_gaussian),
       sigma(other.sigma), mu(other.mu) {}

// Copy constructor
GaussianSampler::GaussianSampler(const GaussianSampler &other,
                                 RandomNumberEngine *random_number_engine)
     : variate_generator_gaussian(*random_number_engine, 
                                  boost::normal_distribution<>(other.mu, other.sigma)),
       sigma(other.sigma), mu(other.mu) {}

// Destructor
GaussianSampler::~GaussianSampler() {}

// Return a Gaussian distribution pseudo-random variate 
double GaussianSampler::sample() {
     return variate_generator_gaussian();
}



/*********************/
/**** DISTRIBUTION ***/
/*********************/

// Constructor - Initialize sampler and density
GaussianDistribution::GaussianDistribution(double sigma, double mu,
                                           RandomNumberEngine *random_number_engine) {

     this->density = new GaussianDensity(sigma, mu);
     this->sampler = new GaussianSampler(sigma, mu, random_number_engine);
}

// Copy constructor
GaussianDistribution::GaussianDistribution(const GaussianDistribution &other) {
     this->density = new GaussianDensity(*other.density);
     this->sampler = new GaussianSampler(*other.sampler);
}

// Copy constructor
GaussianDistribution::GaussianDistribution(const GaussianDistribution &other,
                                           RandomNumberEngine *random_number_engine) {
     this->density = new GaussianDensity(*other.density);
     this->sampler = new GaussianSampler(*other.sampler, random_number_engine);
}

// Destructor
GaussianDistribution::~GaussianDistribution() {
     delete density;
     delete sampler;
}

// Draw sample
double GaussianDistribution::sample() {
     return sampler->sample();
}

// Evaluate log-likelihood
double GaussianDistribution::get_log_likelihood(double x) {
     return density->get_log_likelihood(x);
}


}
