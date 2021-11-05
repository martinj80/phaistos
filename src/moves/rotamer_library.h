// rotamer_lib.h --- Rotamer library base class
// Copyright (C) 2008-2009 Jesper Ferkinghoff-Borg, Wouter Boomsma
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


#ifndef ROTAMER_LIBRARY_H
#define ROTAMER_LIBRARY_H

#include "nodes/discrete.h"
#include "nodes/gaussian.h"
#include "nodes/vonmises.h"
#include "protein/definitions.h"
#include "protein/pdb_input.h"

#include "nodes/vonmises.h"

namespace phaistos {

//! Rotamer Library class
//! Supports sampling of chi angles from a rotamer library of Gaussians or von Mises distributions
template <typename DISTRIBUTION_TYPE>
class RotamerLibrary {
protected:

     //! Random number engine from which random number generators can be created.
     RandomNumberEngine *random_number_engine;     

     //! Wrapper for probability distribution
     class Distribution: public DISTRIBUTION_TYPE {

          //! Random number engine from which random number generators can be created.
          RandomNumberEngine *random_number_engine;     

     public:

          //! Copy constructor
          //! \param distribution Base class object from which copy is made
          Distribution(const DISTRIBUTION_TYPE &distribution)
               : DISTRIBUTION_TYPE(distribution) {}

          //! Constructor - different random number generator
          //! \param distribution Base class object from which copy is made
          //! \param random_number_engine Object from which random number generators can be created.
          Distribution(const DISTRIBUTION_TYPE &distribution,
                       RandomNumberEngine *random_number_engine)
               : DISTRIBUTION_TYPE(distribution,
                                   random_number_engine) {
               assert(random_number_engine!=NULL);
          }

          //! Sample value
          double sample() {
               return DISTRIBUTION_TYPE::sample();
          }

          //! Get log likelihood
          //! Implementations found below
          //! \return log-probability
          double get_log_likelihood(double x);
     };
     

     //! Rotamer - basically a list of distribution objects
     class Rotamer {

          //! Random number engine from which random number generators can be created.
          RandomNumberEngine *random_number_engine;     

     public:
          //! Distribution objects
     	  std::vector<Distribution> distributions;

          //! Constructor
          //! \param distributions vector of distribution objects
          //! \param random_number_engine Object from which random number generators can be created.          
     	  Rotamer(std::vector<Distribution> distributions, 
                  RandomNumberEngine *random_number_engine) 
               : random_number_engine(random_number_engine),
                 distributions(distributions) {}

          // Copy constructor
          //! \param other Source object from which copy is made
          Rotamer(const Rotamer &other)
               : random_number_engine(other.random_number_engine) {
               for (unsigned int i=0; i<other.distributions.size(); i++) {
                    distributions.push_back(Distribution(other.distributions[i],
                                                         other.random_number_engine));
               }
          }

          //! Copy constructor. Different random_number_generator specified.
          //! \param other Source object from which copy is made
          //! \param random_number_engine Object from which random number generators can be created.
          Rotamer(const Rotamer &other,
                  RandomNumberEngine *random_number_engine)
               : random_number_engine(random_number_engine) {
               for (unsigned int i=0; i<other.distributions.size(); i++) {
                    distributions.push_back(Distribution(other.distributions[i],
                                                         random_number_engine));
               }
          }

          //! Sample a rotamer value
          //! \return vector of chi-angle values
	  std::vector<double> sample() {
	       std::vector<double> rValue(distributions.size());
	       for (uint i=0; i<distributions.size(); i++) {
	            rValue[i] = distributions[i].sample();
	       }
	       return rValue;
	  }

          //! Log-likelihood evaluation for vector of chi angles
          //! \return log-probability
          double get_log_likelihood(const std::vector<double> &chi_angles) {
               double LL = 0.0;
               for (unsigned int i=0; i<distributions.size(); i++) {
                    LL += distributions[i].get_log_likelihood(chi_angles[i]);
               }
               return LL;
          }
     };


     //! Mixture model of Rotamers - List of rotamers each associated with a weight
     class RotamerMixture {

          //! Random number engine from which random number generators can be created.
          RandomNumberEngine *random_number_engine;     

     public:

          //! Vector of rotamers
	  std::vector<Rotamer> rotamers;

          //! Vector of weights (one for each rotamer)
	  std::vector<double> weights;

          //! Cumulative probabilities (based on weight vector)
	  double *cumProbs;

          //! Constructor
          //! \param random_number_engine Object from which random number generators can be created.          
	  RotamerMixture(RandomNumberEngine *random_number_engine)
               : random_number_engine(random_number_engine),
                 cumProbs(NULL) {}

          // Copy constructor
          //! \param other Source object from which copy is made
	  RotamerMixture(const RotamerMixture &other)
               : random_number_engine(other.random_number_engine),
                 weights(other.weights) {

               for (unsigned int i=0; i<other.rotamers.size(); i++) {
                    rotamers.push_back(Rotamer(other.rotamers[i], other.random_number_engine));
               }

	       cumProbs = new double[rotamers.size()];
	       DiscreteSampler::calc_cumulative_probs(weights, cumProbs, false);
          }

          //! Copy constructor. Different random_number_generator specified.
          //! \param other Source object from which copy is made
          //! \param random_number_engine Object from which random number generators can be created.
	  RotamerMixture(const RotamerMixture &other, RandomNumberEngine *random_number_engine)
               : random_number_engine(random_number_engine), 
                 weights(other.weights) {

               for (unsigned int i=0; i<other.rotamers.size(); i++) {
                    rotamers.push_back(Rotamer(other.rotamers[i], random_number_engine));
               }

	       cumProbs = new double[rotamers.size()];
	       DiscreteSampler::calc_cumulative_probs(weights, cumProbs, false);
          }

          //! Destructor
          ~RotamerMixture() {
               delete[] cumProbs;
          }

          //! Add a rotamer
          //! \param rotamer Rotamer object
          //! \param weight Weight associated with rotamer
	  void add_rotamer(const Rotamer &rotamer, double weight) {
	       weights.push_back(weight);
	       rotamers.push_back(rotamer);

	       if (cumProbs)
		    delete[] cumProbs;
	       cumProbs = new double[rotamers.size()];
	       DiscreteSampler::calc_cumulative_probs(weights, cumProbs, false);
	  }

          //! Sample rotamer
          //! \param rotamer_index_pointer Optionally extract rotamer index to caller
          //! \return Rotamer object
	  Rotamer &sample(int *rotamer_index_pointer=NULL) {
	       int index = DiscreteSampler::_sample(rotamers.size(), cumProbs,
                                                    random_number_engine);
               if (rotamer_index_pointer)
                    *rotamer_index_pointer = index;

	       return rotamers[index];
	  } 

          //! Get log likelihood of a mixture of Distribution objects
          double get_log_likelihood(const std::vector<double> &chi_angles) {
               double L = 0.0;
               for (unsigned int i=0; i<rotamers.size(); i++) {
                    L += weights[i]*std::exp(rotamers[i].get_log_likelihood(chi_angles));
               }
               return std::log(L);
          }
     };
	  
     //! Parameter container (data[AA][ROTAMER][CHI]
     std::vector<RotamerMixture> aa_rotamers;

public:

     //! Constructor
     //! \param random_number_engine Object from which random number generators can be created.          
     RotamerLibrary(RandomNumberEngine *random_number_engine=&random_global)
          : random_number_engine(random_number_engine) {
          aa_rotamers.resize(20, random_number_engine);
     }

     //! Copy constructor. Different random_number_generator specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     RotamerLibrary(const RotamerLibrary &other,
                    RandomNumberEngine *random_number_engine=&random_global)
          : random_number_engine(random_number_engine) {
          for (unsigned int i=0; i<other.aa_rotamers.size(); i++) {
               aa_rotamers.push_back(RotamerMixture(other.aa_rotamers[i], random_number_engine));
          }
     }     

     //! Sample: residue ID -> chi vector
     //! \return vector of chi angles
     std::vector<double> sample(definitions::ResidueEnum res, int *rotamer_index_pointer=NULL) {

          RotamerMixture &rotamer_mixture = aa_rotamers[(int)res];
          if (rotamer_mixture.rotamers.size() > 0) {
               Rotamer &rotamer = rotamer_mixture.sample(rotamer_index_pointer);
               return rotamer.sample();
          } else {
               return std::vector<double>();
          }
     }


     //! Sample: residue ID -> chi vector - but from specified rotamer
     //! \return vector of chi angles
     std::vector<double> sample(definitions::ResidueEnum res, int rotamer_index) {

          RotamerMixture &rotamer_mixture = aa_rotamers[(int)res];
          if (rotamer_mixture.rotamers.size() > 0) {
               Rotamer &rotamer = rotamer_mixture.rotamers[rotamer_index];
               return rotamer.sample();
          } else {
               return std::vector<double>();
          }
     }

     //! Calculate the likelihood of a sample
     //! \param res Residue type
     //! \param chi_angles vector of chi angles values
     //! \return log-probability
     double get_log_likelihood(definitions::ResidueEnum res, const std::vector<double> &chi_angles) {
          
          RotamerMixture &rotamer_mixture = aa_rotamers[(int)res];
          if (rotamer_mixture.rotamers.size() > 0) {
               return rotamer_mixture.get_log_likelihood(chi_angles);
          } 
          return 0;
     }

     
     //! Add a rotamer
     //! \param res ResidueType
     //! \param rotamer Rotamer object
     //! \param weight Associated weight
     void add_rotamer(definitions::ResidueEnum res, const Rotamer &rotamer, double weight) {
          int aaIndex = (int)res;
          aa_rotamers[aaIndex].add_rotamer(rotamer, weight);
     }
};

//! get_log_likelihood: Default implementation
//! \return log-probability
template<typename DISTRIBUTION_TYPE>
double RotamerLibrary<DISTRIBUTION_TYPE>::Distribution::get_log_likelihood(double x) {
     return DISTRIBUTION_TYPE::get_log_likelihood(x);
}

//! get_log_likelihood: Specialization for Gaussian
//! Ensures that likelihoods of values outside legal range are evaluated correctly
//! \return log-probability
template<>
double RotamerLibrary<GaussianDistribution>::Distribution::get_log_likelihood(double x) {
     double LL1 = GaussianDistribution::get_log_likelihood(x);
     double LL2 = GaussianDistribution::get_log_likelihood(x-2*M_PI);
     double LL3 = GaussianDistribution::get_log_likelihood(x+2*M_PI);
     
     return max(LL1, LL2, LL3);
}

}

#endif
