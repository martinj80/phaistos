// rotamer_lib.h --- Dunbrack backbone independent rotamer library
// Copyright (C) 2008-2009 Wouter Boomsma
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

#ifndef ROTAMER_LIBRARY_DUNBRACK_H
#define ROTAMER_LIBRARY_DUNBRACK_H

#include "rotamer_library.h"
#include "rotamer_library_dunbrack_parameters.h"

namespace phaistos {

//! Interface to the Dunbrack backbone independent rotamer library - derives from RotamerLibrary class.
//! See: R. L. Dunbrack, Jr. and F. E. Cohen.
//! "Bayesian statistical analysis of protein sidechain rotamer preferences ."
//! Protein Science 6, 1661-1681 (1997).
template <typename DISTRIBUTION_TYPE=GaussianDistribution>
class RotamerLibraryDunbrack: public RotamerLibrary<DISTRIBUTION_TYPE> {
private:
public:

     //! Constructor
     //! \param sigma_scale_factor Can optionally be used to flatten/sharpen the distributions 
     //! \param random_number_engine Object from which random number generators can be created.
     RotamerLibraryDunbrack(double sigma_scale_factor=1.0,
                            RandomNumberEngine *random_number_engine=&random_global)
          : RotamerLibrary<DISTRIBUTION_TYPE>(random_number_engine) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Read in parameters
          std::vector<std::string> lines = split(rotamer_library_dunbrack, "\n");
          for (unsigned int i=0; i<lines.size(); i++) {

               // Split line
               std::vector<std::string> split_line = split(lines[i], " ");
               if (split_line.size() < 2) 
                    continue;

               // Extract residue type
               ResidueEnum res = (ResidueEnum)str_to_aa[split_line[0]];

               // Extract weight
               double weight = string_to_val<double>(split_line[1])*0.01;

               // Extract Gaussian (mode, variance)
               std::vector<typename RotamerLibrary<DISTRIBUTION_TYPE>::Distribution> distributions;
               for (unsigned int j=2; (j+1)<split_line.size(); j+=2) {
                    double mu = string_to_val<double>(split_line[j])/180.0*M_PI;
                    double sigma = sigma_scale_factor*string_to_val<double>(split_line[j+1])/180.0*M_PI;
                    distributions.push_back(init_distribution(sigma, mu));
               }

               // Create rotamer
               typename RotamerLibrary<DISTRIBUTION_TYPE>::Rotamer rotamer(distributions,
                                                                           random_number_engine);
 
               // Add rotamer
               this->add_rotamer(res, rotamer, weight);
          }
     }


     //! Copy constructor
     //! \param other Source object from which copy is made
     RotamerLibraryDunbrack(const RotamerLibraryDunbrack &other)
          : RotamerLibrary<DISTRIBUTION_TYPE>(other) {}


     //! Copy constructor. Different random_number_generator specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     RotamerLibraryDunbrack(const RotamerLibraryDunbrack &other,
                            RandomNumberEngine *random_number_engine)
          : RotamerLibrary<DISTRIBUTION_TYPE>(other, random_number_engine) {}

     
     //! Initialize the relevant distribution (Gaussian or von Mises) using a standard deviation and mean parameter
     //! \param sigma Standard deviation
     //! \param mu Mean
     //! \return Distribution object
     typename RotamerLibrary<DISTRIBUTION_TYPE>::Distribution init_distribution(double sigma, double mu);
};

//! init_distribution: Default implementation
//! For Gaussians, sigma and mu are used directly
//! \param sigma Standard deviation
//! \param mu Mean
//! \return Distribution object
template <typename DISTRIBUTION_TYPE>
typename RotamerLibrary<DISTRIBUTION_TYPE>::Distribution RotamerLibraryDunbrack<DISTRIBUTION_TYPE>::init_distribution(double sigma, double mu) {
     return typename RotamerLibrary<DISTRIBUTION_TYPE>::Distribution(DISTRIBUTION_TYPE(sigma,mu));
}

//! init_distribution: Specialization for VonMises
//! For Von Mises distribution, sigma is converted to the corresponding kappa value
//! \param sigma Standard deviation
//! \param mu Mean
//! \return Distribution object
template <>
RotamerLibrary<VonMisesDistribution>::Distribution RotamerLibraryDunbrack<VonMisesDistribution>::init_distribution(double sigma, double mu) {
     double kappa = 1.0/(sigma*sigma);
     return RotamerLibrary<VonMisesDistribution>::Distribution(VonMisesDistribution(kappa,mu));
}

}

#endif
