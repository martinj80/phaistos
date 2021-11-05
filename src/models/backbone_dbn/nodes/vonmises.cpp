// vonmises.cpp --- von Mises node. Density and sampler for a univariate von Mises distritbution
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

#include "vonmises.h"

#include "utils/vector_nd.h"

#include <iostream>
#include <math.h>

#include <boost/math/special_functions/bessel.hpp>

namespace phaistos {

/*****************/
/**** DENSITY ****/
/*****************/

// Constructor
VonMisesDensity::VonMisesDensity(double k, double mu) {
     this->k = k;
     this->mu = mu;
     // this->log_norm_const = std::log(2*M_PI*i0(k));
     this->log_norm_const = std::log(2*M_PI*boost::math::cyl_bessel_i(0, k));
}

// Destructor
VonMisesDensity::~VonMisesDensity() {}


// Calculate log likelihood
double VonMisesDensity::get_log_likelihood(double angle) {
     double dens = k*cos(angle-mu) - log_norm_const;
     assert(std::isfinite(dens));
     return dens;
}

// Calculate likelihood
double VonMisesDensity::get_likelihood(double angle) {
     return std::exp(get_log_likelihood(angle));
}



/*****************/
/**** SAMPLER ****/
/*****************/

// Sample from a von Mises distribution
VonMisesSampler::VonMisesSampler(double k, double mu,
                                 RandomNumberEngine *random_number_engine)
     : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     this->k = k;
     this->mu = mu;
     if (mu < -M_PI || mu > M_PI) {
          fprintf(stderr, "VonMises Error: mu must be in the interval (-pi,pi). Mu=%f\n", mu);
     }

     // Set sampling internals
     a = 1.0 + sqrt(1.0 + 4.0*k*k);
     b = (a - sqrt(2*a)) / (2*k);
     r = (1.0 + b*b)/(2*b);
}

// Copy constructor
VonMisesSampler::VonMisesSampler(const VonMisesSampler &other)
     : variate_generator_uniform(other.variate_generator_uniform) {

     k = other.k;
     mu = other.mu;

     a = other.a;
     b = other.b;
     r = other.r;
}

// Copy constructor
VonMisesSampler::VonMisesSampler(const VonMisesSampler &sampler,
                                 RandomNumberEngine *random_number_engine)
     : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     k = sampler.k;
     mu = sampler.mu;

     a = sampler.a;
     b = sampler.b;
     r = sampler.r;
}

// Destructor
VonMisesSampler::~VonMisesSampler() {}

// Return a von Mises distribution pseudo-random variate on [-pi, +pi].
// The implementation is similar to the algorithm by Best and Fisher,
// 1979; see N.I. Fisher, "Statistical Analysis of Circular Data",
// Cambridge University Press, 1993, p. 49.
double VonMisesSampler::sample() {

     double f,res;
     while(true) {
          double U1 = variate_generator_uniform();
          double z = cos(M_PI * U1);
          f = (1.0 + r*z)/(r + z);
          double c = k * (r - f);
          double U2 = variate_generator_uniform();
          if (((c*(2.0-c) - U2) > 0.0) || ((std::log(c/U2) + 1.0 - c >= 0.0))){
               break;           // accept
          }
     }
     double U3 = variate_generator_uniform();
     if (U3 > 0.5) {
          res = fmod(acos(f)+mu, 2*M_PI);
     } else {
          res = fmod(-acos(f)+mu, 2*M_PI);
     }
     return res;
}



/*********************/
/**** DISTRIBUTION ***/
/*********************/

// Constructor - Initialize sampler and density
VonMisesDistribution::VonMisesDistribution(double k, double mu,
                                           RandomNumberEngine *random_number_engine) {

     this->density = new VonMisesDensity(k, mu);
     this->sampler = new VonMisesSampler(k, mu, random_number_engine);
}

// Copy constructor
VonMisesDistribution::VonMisesDistribution(const VonMisesDistribution &other) {
     this->density = new VonMisesDensity(*other.density);
     this->sampler = new VonMisesSampler(*other.sampler);
}

// Copy constructor
VonMisesDistribution::VonMisesDistribution(const VonMisesDistribution &other,
                                           RandomNumberEngine *random_number_engine) {
     this->density = new VonMisesDensity(*other.density);
     this->sampler = new VonMisesSampler(*other.sampler, random_number_engine);
}

// Destructor
VonMisesDistribution::~VonMisesDistribution() {
     delete density;
     delete sampler;
}

// Draw sample
double VonMisesDistribution::sample() {
     return sampler->sample();
}

// Evaluate log-likelihood
double VonMisesDistribution::get_log_likelihood(double x) {
     return density->get_log_likelihood(x);
}

}
