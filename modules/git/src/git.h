// git.h --- Interface to Gauss Integral Tuned code
// Copyright (C) 2011  Tim Harder, Thomas Hamelryck, Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef PHAISTOS_GIT_H
#define PHAISTOS_GIT_H

#include "../external/git/git.h"
#include "utils/random.h"
#include "utils/vector_matrix_3d.h"
#include "protein/definitions.h"

namespace phaistos {

//! Phaistos specific Gauss Integral Tuned (GIT) class 
//! Inherits from general GIT library
class Git: public git::Git {
private:

     //! Random number generator - uniform
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > random_generator_uniform_01;


     //! Draw a uniformly distributed random value.     
     //! Override default implementation to use the
     //! random number engine used in Phaistos
     virtual double sample_uniform_01() {
          return random_generator_uniform_01();
     }

     //! File to which output is written
     std::ofstream output_file;

public:

     //! Constructor (no output file)
     //!
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param parameter_filename Filename containing average gauss
     //! integrals parameters. If not specified, a meaningful default
     //! is used.
     //! \param smoothen_backbone_mode Whether to smoothen the backbone     
     Git(RandomNumberEngine *random_number_engine=&random_global, const std::string parameter_filename="", const bool smoothen_backbone_mode=1)
          : git::Git(parameter_filename, smoothen_backbone_mode),
            random_generator_uniform_01(*random_number_engine, boost::uniform_real<>(0,1)) {}


     //! Constructor (with output file)
     //!
     //! \param output_filename Filename to which output will be written
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param parameter_filename Filename containing average gauss
     //! integrals parameters. If not specified, a meaningful default
     //! is used.
     //! \param smoothen_backbone_mode Whether to smoothen the backbone     
     Git(const std::string output_filename, RandomNumberEngine *random_number_engine=&random_global, const std::string parameter_filename="", const bool smoothen_backbone_mode=1)
          : git::Git(parameter_filename, smoothen_backbone_mode),
            random_generator_uniform_01(*random_number_engine, boost::uniform_real<>(0,1)) {
          set_output_file(output_filename);
     }


     //! Destructor
     ~Git() {
          if (output_file.is_open()) {
               output_file.close();
          }
     }
     

     //! Associate output file with git object
     //!
     //! \param output_filename File name
     void set_output_file(std::string output_filename) {
          if (output_file.is_open()) {
               output_file.close();
          }
          output_file.open(output_filename.c_str());
     }


     //! Generate Gauss Integrals from molecule chain object
     //!
     //! \param chain Molecule chain object
     //! \param name Optional name (used when outputting to file)
     //! \return 30 dimensional GIT vector     
     template <typename CHAIN_TYPE>
     std::vector<double> generate_gauss_integrals(CHAIN_TYPE &chain, std::string name="") {

          std::vector<Vector_3D> positions = chain.template get_positions<definitions::CA_ONLY>();

          std::vector<git::Vec3> ca_trace(positions.size());
          for (unsigned int i=0; i<positions.size(); ++i) {
               ca_trace[i] = git::Vec3(positions[i].get_array());
          }

          std::vector<double> git_output_vector = git::Git::generate_gauss_integrals(ca_trace);

          // Output to file if filename was specified during construction
          if (output_file.is_open()) {
               if (name!="")
                    output_file << name << " ";
               output_file << positions.size() << " " << std::fixed << std::setprecision(4) << 19.11*pow(positions.size(),1.0/3.0) << " ";
               for (unsigned int i=0; i<30; ++i) {
                    output_file << std::fixed << std::setprecision(4) << git_output_vector[i] << " ";
               }
               output_file << "\n";
               output_file.flush();
          }
          

          return git_output_vector;
     }
};


}



#endif
