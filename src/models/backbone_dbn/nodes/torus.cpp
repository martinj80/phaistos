// torus.cpp --- Torus node. Density and sampler for a bivariate von Mises distribution - cosine model
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

#include "torus.h"

#include "utils/netlib/integrate.h"
#include "utils/optimize.h"
#include "utils/vector_nd.h"

#include <iostream>
#include <math.h>

#include <boost/math/special_functions/bessel.hpp>

namespace phaistos {

// For values above this limit the calculations are scaled to avoid overflow
#define LIMITSCALING 300
// For values above this limit an asymptotic expression for the bessel function is used
#define LIMITAPPROXBESSEL 300

#define i0(X) boost::math::cyl_bessel_i(0, X)
#define i1(X) boost::math::cyl_bessel_i(1, X)

// Exact integrand for normalization
double log_norm_const_exact_integrand(double *phi, double *extra_arguments) {
     double k1 = extra_arguments[0];
     double k2 = extra_arguments[1];
     double k3 = extra_arguments[2];
     double k13 = sqrt((k1*k1) + (k3*k3) - 2*k1*k3*cos(*phi));
     double res = (2*M_PI*std::exp(std::log(i0(k13)) + k2*cos(*phi)));
     return res;
}

// Integrand using approximation for Bessel function (works for large values of k13)
double log_norm_const_approximate_integrand(double *phi, double *extra_arguments) {
     double k1 = extra_arguments[0];
     double k2 = extra_arguments[1];
     double k3 = extra_arguments[2];
     double scale_term = extra_arguments[3];
     double k13 = sqrt((k1*k1) + (k3*k3) - 2*k1*k3*cos(*phi));
     double e = std::exp(k13 + k2*cos(*phi) - scale_term);
     return sqrt((2*M_PI)/k13)*e;
}


// Calculate logarithm of normalize constant c(k1, k2, k3)
double compute_log_norm_const(Vector_nD &k, double lower=-M_PI, double upper=M_PI) {
     double k1 = k[0];
     double k2 = k[1];
     double k3 = k[2];
     double exact_limit = 400;
     if ((sqrt(k1*k1 + k2*k2 - 2*k1*k3) + k2) < exact_limit) {
          double extra_arguments[3] = {k1, k2, k3};
          double res_exact = integrate_quad(&log_norm_const_exact_integrand, extra_arguments, lower, upper);
          return (std::log(res_exact));
     } else {
          double scale_term = sqrt(k1*k1 + k3*k3 - 2*k1*k3) + k2;
          double extra_arguments[4] = {k1, k2, k3, scale_term};
          double res_approximate = integrate_quad(&log_norm_const_approximate_integrand, extra_arguments, lower, upper);
          return (std::log(res_approximate) + scale_term);
     }
     return 0;
}

// Density function for a univariate von Mises distribution
double vm1_cos_density(double phi, double k, double mu, double scale_term=0.0) {
     double res;
     // For k>limitApproxBessel an asymptotic expression for the bessel function is used
     if (k<LIMITAPPROXBESSEL) {
          // res = std::exp(k*cos(phi-mu) - k - std::log(i0e(k)) - std::log(2*M_PI));
          res = std::exp(k*cos(phi-mu) - std::log(i0(k)) - std::log(2*M_PI));
     } else {
          res = std::exp(k*cos(phi-mu) - k - 0.5*std::log(4*M_PI*M_PI*k) - scale_term);
     }
     return res;
}

//! Mixture of univariate von Mises. Used as comparison function in rejection sampling.
class VM1CosMixtureDensity {
public:
     
     //! kappa parameters
     Vector_nD k;
     
     //! mu parameters
     Vector_nD mu;
     
     //! Weight of the mixture components
     Vector_nD proportion;

     //! Cumulative weight of the mixture components
     Vector_nD cum_proportion;

     //! Constructor
     //! \param k kappa parameters (k1,k2,k3)
     //! \param mu mu parameters (mu1,mu2)
     //! \param proportion weight of the mixture components
     VM1CosMixtureDensity(Vector_nD k, Vector_nD mu, Vector_nD proportion) {
          this->k = k;
          this->mu = mu;
          this->proportion = proportion / sum(proportion);
          this->cum_proportion = cum_sum(proportion);
     }

     //! Calculate probability density value
     //! \param phi observed value
     //! \param scale_term Used to avoid overflow
     //! \return likelihood
     double density(double phi, double scale_term=0.0) {
          double res = 0.0;
          for (int i=0; i<k.size; i++) {
               res += proportion[i] * vm1_cos_density(phi, k[i], mu[i], scale_term);
          }
          return res;
     }


     //! Calculate log-probability density value
     //! \param phi observed value
     //! \return log-likelihood
     double log_density(double phi) {
          Vector_nD tmp = fmap(&cos, -(mu-phi))-1;
          double scale_term = -(-0.5*(max(tmp)+min(tmp)))*max(k);

          double res;
          // For k>limitScaling the calculations are scaled to avoid overflow
          if (max(k) < LIMITSCALING) {
               res = std::log(density(phi));
          } else {
               res = std::log(density(phi, scale_term)) + scale_term;
          }
          return res;
     }

     //! Calculate log-probability density value for a vector of values
     //! \param phi vector of observed values
     //! \return log-likelihood
     Vector_nD log_density(Vector_nD phi) {
          Vector_nD res(phi.size);
          for (int i=0; i<phi.size; i++) {
               res[i] = log_density(phi[i]);
          }
          return res;
     }
};


/*****************/
/**** DENSITY ****/
/*****************/

// Constructor
TorusDensity::TorusDensity(Vector_nD &k, Vector_nD &mu) {
     init(k, mu);
     this->log_norm_const = compute_log_norm_const(k);
}

// Constructor - logarithm of normalization constant specified
TorusDensity::TorusDensity(Vector_nD &k, Vector_nD &mu, double log_norm_const) {
     init(k, mu);
     this->log_norm_const = log_norm_const;     
}

// Destructor
TorusDensity::~TorusDensity() {}


// Initializer
void TorusDensity::init(Vector_nD &k, Vector_nD &mu) {
     this->k = k;
     this->mu = mu;
}


// Return a von Mises distribution pseudo-random variate on [-pi, +pi].
// The implementation is similar to the algorithm by Best and Fisher,
// 1979; see N.I. Fisher, "Statistical Analysis of Circular Data",
// Cambridge University Press, 1993, p. 49.
double von_mises_sampler(double k, double mean, 
                         boost::variate_generator<RandomNumberEngine&, 
                         boost::uniform_real<> > &variate_generator_uniform) {

     if (mean < -M_PI || mean > M_PI) {
          fprintf(stderr, "vonMises Error: mean must be in the interval (-pi,pi). Mean=%f\n", mean);
     }
     double res;
     double a = 1.0 + sqrt(1.0 + 4.0*k*k);
     double b = (a - sqrt(2*a)) / (2*k);
     double r = (1.0 + b*b)/(2*b);
     double f;
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
          res = fmod(acos(f)+mean, 2*M_PI);
     } else {
          res = fmod(-acos(f)+mean, 2*M_PI);
     }
     return res;
}



//! Mixture of univariate von Mises. Used as comparison function in rejection sampling.
class VM1CosMixtureSampler {
private:

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > variate_generator_uniform;

     //! kappa parameters
     Vector_nD k;

     //! mu parameters
     Vector_nD mu;

     //! Weight of the mixture components
     Vector_nD proportion;

     //! Cumulative weight of the mixture components
     Vector_nD cum_proportion;
     
public:
     //! Constructor
     //! \param k kappa parameters (k1,k2,k3)
     //! \param mu mu parameters (mu1,mu2)
     //! \param proportion weight of the mixture components
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     VM1CosMixtureSampler(Vector_nD &k, Vector_nD &mu, Vector_nD &proportion, 
                          RandomNumberEngine *random_number_engine)
          : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)),
            k(k), mu(mu) {
          
          this->proportion = proportion / sum(proportion);
          this->cum_proportion = cum_sum(proportion);
     }

     //! Copy constructor
     //! \param other Source object          
     VM1CosMixtureSampler(const VM1CosMixtureSampler &other)
          : variate_generator_uniform(other.variate_generator_uniform),
            k(other.k), mu(other.mu), 
            proportion(other.proportion), cum_proportion(other.cum_proportion) {}

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.          
     VM1CosMixtureSampler(const VM1CosMixtureSampler &other, 
                          RandomNumberEngine *random_number_engine)
          : variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)),
            k(other.k), mu(other.mu), 
            proportion(other.proportion), cum_proportion(other.cum_proportion) {}


     //! Return single sample
     //! \return angle
     double operator()() {
          double k=0.0;
          double mu=0.0;
          double u = variate_generator_uniform();
          for (int i=0; i<cum_proportion.size; i++) {
               if (cum_proportion[i] >= u) {
                    k = this->k[i];
                    mu = this->mu[i];
                    break;
               }
          }
          double res = von_mises_sampler(k, mu, variate_generator_uniform);
          return res;
     }

     //! Return n samples
     //! \return vector of angles
     Vector_nD operator()(int n) {
          Vector_nD res(n);
          for (int i=0; i<n; i++) {
               res[i] = (*this)();
          }
          return res;
     }
};

//! Bivariate von Mises - cosine model - Marginal of phi - Density function
class VM2CosMarginalDensity {
private:

     //! kappa parameters     
     Vector_nD k;

     //! logarithm of the normalization constant     
     double log_norm_const;
public:

     //! Constructor
     //! kappa parameters     
     VM2CosMarginalDensity(Vector_nD &k) {
          init(k);
          this->log_norm_const = compute_log_norm_const(k);
     }

     //! Constructor - logarithm of normalization contant specified
     //! \param k kappa parameters     
     //! \param log_norm_const Logarithm of the normalization constant
     VM2CosMarginalDensity(Vector_nD &k, double log_norm_const) {
          init(k);
          this->log_norm_const = log_norm_const;
     }

     //! Initializer
     //! kappa parameters     
     void init(Vector_nD &k) {
          this->k = k;
     }

     //! Return log density for single phi
     //! \param phi obverved value
     //! \return log-likelihood
     double log_density(double phi) {
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(phi));
          
          double res;
          // For values above limitApproxBessel an asymptotic expression for the bessel function is used
          if (k1*k1 + k2*k2 < LIMITAPPROXBESSEL) {
               res = std::log(2*M_PI) - log_norm_const + std::log(i0(k13)) + k2*cos(phi);
                    } else {
               res = -log_norm_const + std::log(sqrt((2*M_PI)/k13)) + k13 + k2*cos(phi);
          }
          return res;
     }

     //! Return log densities for vector of phis
     //! \param phi vector of obverved values
     //! \return log-likelihood
     Vector_nD log_density(Vector_nD &phi) {
          Vector_nD res(phi.size);
          for (int i=0; i<phi.size; i++) {
               res[i] = log_density(phi[i]);
          }
          return res;
     }

     //! Return density for single phi
     //! \param phi obverved value
     //! \return likelihood
     double density(double phi) {
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(phi));
          double res =  2*M_PI*i0(k13)*std::exp(k2*cos(phi));
          return res;
     }

     //! Overload () - returns the negative density function
     //! \param phi observed value
     //! \return negative log-likelihood
     double operator()(double phi) {
          return -log_density(phi); // Note that this is the NEGATIVE density function
     }
     
};


//! Bivariate von Mises - cosine model. Marginal of phi. Sampler
class VM2CosMarginalSampler {
private:

     //! Random number engine
     RandomNumberEngine *random_number_engine;

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > variate_generator_uniform;

     //! kappa parameters
     Vector_nD k;

     //! logarithm of the normalization constant     
     double log_norm_const;

     //! Root of derivate - maximum of density function
     double root;

     //! Brent optimizer for root finding
     Brent<double> brent_optimizer;
     
     //! Envelope sampler
     VM1CosMixtureSampler *envelope_sampler;

     //! Envelope density
     VM1CosMixtureDensity *envelope_density;

     //! Scale parameter
     double K;
     
public:

     //! Initializer
     //! \param k kappa values
     void init(Vector_nD &k) {
          this->k = k;
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          
          this->log_norm_const = compute_log_norm_const(k);

          // Determine whether marginal is bimodal and set roots accordingly
          // We can use the exponentially scaled Bessel functions
          // since the same scaling factor appears in the numerator
          // and the denominator
          root = 0.0;
          // if (((i1e(fabs(k1 - k3)) / i0e(fabs(k1-k3))) > ((fabs(k1-k3) * k2)/(k1*k3))) && // Multimodal criterium Theorem 4
          if (((i1(std::fabs(k1 - k3)) / i0(std::fabs(k1-k3))) > ((std::fabs(k1-k3) * k2)/(k1*k3))) && // Multimodal criterium Theorem 4
              ((derivativeFactorB(-M_PI) * derivativeFactorB(0.0)) < 0)) { // In case k1>k3>0 and k2>k3>0 is not fulfilled
               
               // Find root of b (theorem 4)
               VM2CosMarginalDensity density_object(k, this->log_norm_const);
               brent_optimizer.optimize(&density_object, Bracket<double>(&density_object, 0, -M_PI/2, -M_PI));

               root = brent_optimizer.xmin;
          }

          // find parameters for mixture of von Mises distributions that fits best
          envelope_sampler = NULL;
          find_optimal_comparison_function();
     }

     //! Constructor
     //! \param k kappa values
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     VM2CosMarginalSampler(Vector_nD &k, RandomNumberEngine *random_number_engine)
          : random_number_engine(random_number_engine),
            variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

          init(k);
     }

     //! Constructor
     //! \param k Kappa parameters (k1,k2,k3)
     //! \param log_norm_const Logarithm of the normalization constant
     //! \param envelope_scale The overall height of the envelope
     //! \param envelope_k Envelope kappa parameters
     //! \param envelope_mu Envelope mu parameters
     //! \param envelope_proportion Envelope proportion parameters
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     VM2CosMarginalSampler(Vector_nD &k, 
                           double *log_norm_const, double *envelope_scale, Vector_nD &envelope_k,
                           Vector_nD &envelope_mu, Vector_nD &envelope_proportion,
                           RandomNumberEngine *random_number_engine)
          : random_number_engine(random_number_engine),
            variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

          this->k = k;

          if (envelope_k.initialized() &&
              envelope_mu.initialized() &&
              envelope_proportion.initialized()) {

               this->K = *envelope_scale;
               this->log_norm_const = *log_norm_const;
               this->envelope_sampler = new VM1CosMixtureSampler(envelope_k, envelope_mu, envelope_proportion, 
                                                                random_number_engine);
               this->envelope_density = new VM1CosMixtureDensity(envelope_k, envelope_mu, envelope_proportion);
          } else {
               init(k);
               *envelope_scale = K;
               *log_norm_const = this->log_norm_const;
               envelope_k = envelope_density->k;
               envelope_mu = envelope_density->mu;
               envelope_proportion = envelope_density->proportion;
          }
     }

     //! Copy constructor
     //! \param other Source object          
     VM2CosMarginalSampler(const VM2CosMarginalSampler &other)
          : random_number_engine(other.random_number_engine),
            variate_generator_uniform(other.variate_generator_uniform) {

          k = other.k;
          log_norm_const = other.log_norm_const;
          root = other.root;

          envelope_sampler = new VM1CosMixtureSampler(*other.envelope_sampler, random_number_engine);
          envelope_density = new VM1CosMixtureDensity(*other.envelope_density);
          K = other.K;
     }

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.          
     VM2CosMarginalSampler(const VM2CosMarginalSampler &other, 
                           RandomNumberEngine *random_number_engine)
          : random_number_engine(random_number_engine),
            variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

          k = other.k;
          log_norm_const = other.log_norm_const;
          root = other.root;

          envelope_sampler = new VM1CosMixtureSampler(*other.envelope_sampler, random_number_engine);
          envelope_density = new VM1CosMixtureDensity(*other.envelope_density);
          K = other.K;
     }

     //! destructor
     ~VM2CosMarginalSampler() {
          if (envelope_sampler)
               delete envelope_sampler;
          if (envelope_density)
               delete envelope_density;
     }

     //! Generate a random sample from the marginal distribution.
     //! Rejection sampling: 1) Sample random value r from comparison distribution
     //!                     2) Sample random value u from uniform distribution
     //!                     3) Accept r if u < f_target(r)/f_comparison(r) """
     //! \param status pointer to retrieve status of sampling
     //! \return sampled value
     double operator()(bool *status) {

         VM2CosMarginalDensity density(k, log_norm_const);
         int count = 0;         
         while(1) {
              double r = (*this->envelope_sampler)();
              double u = variate_generator_uniform();
              double x = std::exp(density.log_density(r) - (K + envelope_density->log_density(r)));
              
              if (u<x) {
                   *status=true;
                   return(r);
              }
              
              if ((count > 0) && count%(10)==0) {
                   fprintf(stderr, "Large rejection count: %d\n", count);
                   if (count > 50) {
                        *status=false;
                        return 0.0;
                   }
              }
              count++;
         }
     }
     
     //! Generate n samples
     //! \param n Number of samples to generate
     //! \return vector of samples
     Vector_nD operator()(int n) {
          
          Vector_nD phi(n);
          int i=0;
          while (i<n) {
               bool success;
               phi[i] = (*this)(&success);
               
               if (success) {
                    i++;
               } else {
                    return Vector_nD();
               }
          }
          return phi;
     }
     
     //! Calculate factor b from derivative of density. Defined in proof of theorem 4
     //! \param phi angle
     double derivativeFactorB(double phi){
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];          
          double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(phi));
          // return (-k2 + k1*k3*(i1e(k13)/i0e(k13))/(k13));
          return (-k2 + k1*k3*(i1(k13)/i0(k13))/(k13));
     }


     //! Comparison function: Mixture of two univariate von Mises distributions.
     //! Find parameter k and scale factor K so the the target marginal density is at
     //! no point above the comparison function.
     class MaxRatio {
     private:
          //! kappa parameters
          Vector_nD k;

          //! Root of derivate - maximum of density function
          double root;

          //! logarithm of the normalization constant     
          double log_norm_const;

     public:

	  //! Constructor
          //! \param k kappa values
          //! \param root root value
          //! \param log_norm_const Logarithm of the normalization constant
          MaxRatio(Vector_nD &k, double root, double log_norm_const) {
               this->k = k;
               this->root = root;
               this->log_norm_const = log_norm_const;
          }

	  //! Reparameterize k parameters
          //! \param k kappa values
          void reparameterisation(double *k) {
               *k = std::fabs(*k);
          }

	  //! Overload () operator - calculate parameters. 
          //! \param x Input value
          //! \return distance between envelope and density
          double operator()(Vector_nD x) {
               return (*this)(x[0]);
          }

	  //! Overload () operator - calculate parameters. 
          //! \param x Input value
          //! \return distance between envelope and density
          double operator()(double x) {
               if (std::isnan(x)) {
                    fprintf(stderr, "k=nan proposed. Returning inf.");
                    assert(false);
                    return INFINITY;
               }

               reparameterisation(&x);

               Vector_nD phi = range(-M_PI, M_PI, 0.01);  // Hack: Only check function at 0.01 intervals
                                                         // instead of doing full minimization.

               VM1CosMixtureDensity envelope_density(Vector_nD(2,x,x), Vector_nD(2, root, -root), Vector_nD(2, 0.5, 0.5));
               VM2CosMarginalDensity density(k, log_norm_const);
               double res = max(density.log_density(phi) - envelope_density.log_density(phi));
               return(res);
          }
     };
          

     //! Comparison function: Mixture of four univariate von Mises distributions
     //! This method is used for large values of k."""
     //! Find parameter k and scale factor K so the the target marginal density is at
     //! no point above the comparison function.
     class MaxRatioLargerMixture {
     private:
          //! kappa parameters
          Vector_nD k;

          //! Root of derivate - maximum of density function
          double root;

          //! logarithm of the normalization constant     
          double log_norm_const;

     public:

	  //! Constructor
          //! \param k kappa values
          //! \param root root value
          //! \param log_norm_const Logarithm of the normalization constant
          MaxRatioLargerMixture(Vector_nD &k, double root, double log_norm_const) {
               this->k = k;
               this->root = root;
               this->log_norm_const = log_norm_const;
          }

          //! Reparameterisation of parameters.
          //! k values are forced to be positive
          //! proportion values are forced to be between 0 and 1"""
          void reparameterisation(Vector_nD *x, int k_size, bool inverse=false) {
               // first k_size elements of x are k values
               assert (k_size <= x->size);
               for (int i=0; i<k_size; i++) {
                    (*x)[i] = std::fabs((*x)[i]);
               }

               for (int i=k_size; i<x->size; i++) {
                    (*x)[i] = std::fabs((*x)[i]);
                    
                    if (!inverse) {
                         (*x)[i] = 1-(1/std::exp((*x)[i]));
                    } else {
                         (*x)[i] = -std::log(1-(*x)[i]);
                    }
               }
          }

	  //! Overload () operator - calculate parameters. 
          //! \param x Input value
          //! \return distance between envelope and density
          double operator()(Vector_nD x) {
               reparameterisation(&x, 2);

               Vector_nD phi = range(-M_PI, M_PI, 0.01);  // Hack: Only check function at 0.01 intervals
                                                          // instead of doing full minimization.
               VM1CosMixtureDensity envelope_density(Vector_nD(4, x[0], x[0], x[1], x[1]),
                                                    Vector_nD(4, root, -root, root, -root),
                                                    Vector_nD(4, 0.5, 0.5, x[2], x[2]));
               VM2CosMarginalDensity density(k, log_norm_const);
               double res = max(density.log_density(phi) - envelope_density.log_density(phi));
               return(res);          
          }
     };
     
     //! Minimize the necessary scale factor K
     void find_optimal_comparison_function() {

          MaxRatio function(k,root,log_norm_const);

          Powell powell;
          powell.optimize(&function, Vector_nD(1, min(k)));

          function.reparameterisation(&powell.xmin[0]);
          
          Vector_nD k(2, powell.xmin[0], powell.xmin[0]);
          Vector_nD mu(2, root, -root);
          Vector_nD proportion(2, 0.5, 0.5);

          // Scale factor
          this->K = powell.fmin;

          if (K > 5) {
               // Retry with mixture model with four components
               MaxRatioLargerMixture function2(k,root,log_norm_const);
               Vector_nD startValues(3, min(k), 0.1, 0.2);
               bool inverse = true;
               function2.reparameterisation(&startValues, 2, inverse);

               Powell powell2;
               powell2.optimize(&function2, startValues);

               function2.reparameterisation(&powell2.xmin, 2);
               k = Vector_nD(4, powell2.xmin[0], powell2.xmin[0], powell2.xmin[1], powell2.xmin[1]);
               mu = Vector_nD(4, root, -root, root, -root);
               proportion = Vector_nD(4, 0.5, 0.5, powell2.xmin[2], powell2.xmin[2]);

               this->K = powell2.fmin;
          }

          this->envelope_sampler = new VM1CosMixtureSampler(k, mu, proportion, random_number_engine);
          this->envelope_density = new VM1CosMixtureDensity(k, mu, proportion);
     }        
};
     



/*****************/
/**** SAMPLER ****/
/*****************/

// Sample from a bivariate von Mises distritbution cosine model
// The distribution and the techniques used in this code are described in:
// 'Bivariate Densities for Angular Data' (draft), Mardia, Subramaniam, Taylor(2005)
// Any references found in comments are to this paper.
TorusSampler::TorusSampler(Vector_nD &k, Vector_nD &mu,
                           double *log_norm_const, double *envelope_scale, Vector_nD &envelope_k,
                           Vector_nD &envelope_mu, Vector_nD &envelope_proportion,
                           RandomNumberEngine *random_number_engine)
     : random_number_engine(random_number_engine),
       variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     this->k = k;
     this->mu = mu;
     
     // Sampling of VM2Cos is done by first sampling from the marginal distribution.
     
     // Assign marginal sampler
     marginal_phi_sampler = new VM2CosMarginalSampler(k,
                                                      log_norm_const, envelope_scale,
                                                      envelope_k, envelope_mu, envelope_proportion,
                                                      random_number_engine);     
}

TorusSampler::TorusSampler(Vector_nD &k, Vector_nD &mu,
                           RandomNumberEngine *random_number_engine)
     : random_number_engine(random_number_engine),
       variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     this->k = k;
     this->mu = mu;
     
     // Sampling of VM2Cos is done by first sampling from the marginal distribution.
     
     // Assign marginal sampler
     marginal_phi_sampler = new VM2CosMarginalSampler(k, random_number_engine);     
}

// Copy constructor
TorusSampler::TorusSampler(const TorusSampler &other)
     : random_number_engine(other.random_number_engine),
       variate_generator_uniform(other.variate_generator_uniform) {

     k = other.k;
     mu = other.mu;

     if (other.marginal_phi_sampler) {
          marginal_phi_sampler = new VM2CosMarginalSampler(*other.marginal_phi_sampler, 
                                                           random_number_engine);
     } else {
          marginal_phi_sampler = NULL;
     }
}

// Copy constructor
TorusSampler::TorusSampler(const TorusSampler &sampler,
                           RandomNumberEngine *random_number_engine)
     : random_number_engine(random_number_engine),
       variate_generator_uniform(*random_number_engine, boost::uniform_real<>(0,1)) {

     k = sampler.k;
     mu = sampler.mu;

     if (sampler.marginal_phi_sampler) {
          marginal_phi_sampler = new VM2CosMarginalSampler(*sampler.marginal_phi_sampler, 
                                                           random_number_engine);
     } else {
          marginal_phi_sampler = NULL;
     }
}

// Destructor
TorusSampler::~TorusSampler() {
     delete marginal_phi_sampler;
}

// Generate a single sample
Vector_nD TorusSampler::sample() {
     // First sample the marginal of psi f(psi). Then use these
     // values to sample from the conditional f(phi|psi). See Section 2 in paper

     if (!marginal_phi_sampler) {
          marginal_phi_sampler = new VM2CosMarginalSampler(k, random_number_engine);
     }
     
     bool success;
     double psi = (*this->marginal_phi_sampler)(&success);

     if (!success)
          return Vector_nD();

     double k1 = k[0];
     double k3 = k[2];
     double mu1 = mu[0];
     double mu2 = mu[1];
     double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(psi));
     double psi_mu = atan((-k3*sin(psi))/(k1 - k3*cos(psi)));
     double phi = von_mises_sampler(k13, psi_mu, variate_generator_uniform);

     // Angles are in the interval [-pi, pi]. Add 3*pi to bring them
     // to: [2*pi, 4*pi] which with mu values in [-pi, pi] brings it
     // to [pi, 3pi], which can then readily be transformed back to [-pi, pi]
     psi = fmod(psi + (3*M_PI) + mu2, 2*M_PI) - M_PI;
     phi = fmod(phi + (3*M_PI) + mu1, 2*M_PI) - M_PI;

     return Vector_nD(2, phi, psi);
}



/*********************/
/**** DISTRIBUTION ***/
/*********************/

// Constructor - Initialize sampler and density
TorusDistribution::TorusDistribution(Vector_nD &k, Vector_nD &mu,
                                     double *log_norm_const, double *envelope_scale, Vector_nD &envelope_k,
                                     Vector_nD &envelope_mu, Vector_nD &envelope_proportion,
                                     RandomNumberEngine *random_number_engine) {
     this->density = new TorusDensity(k, mu);
     this->sampler = new TorusSampler(k, mu,
                                      log_norm_const, envelope_scale,
                                      envelope_k, envelope_mu, envelope_proportion,
                                      random_number_engine);
}

// Copy constructor
TorusDistribution::TorusDistribution(const TorusDistribution &other) {
     this->density = new TorusDensity(*other.density);
     this->sampler = new TorusSampler(*other.sampler);
}

// Copy constructor
TorusDistribution::TorusDistribution(const TorusDistribution &other, 
                                     RandomNumberEngine *random_number_engine) {
     this->density = new TorusDensity(*other.density);
     this->sampler = new TorusSampler(*other.sampler, random_number_engine);
}

// Destructor
TorusDistribution::~TorusDistribution() {
     delete density;
     delete sampler;
}

// Draw sample
Vector_nD TorusDistribution::sample() {
     return sampler->sample();
}

// Evaluate log-likelihood
double TorusDistribution::get_log_likelihood(double *angle_pair) {
     return density->get_log_likelihood(angle_pair);
}

}
