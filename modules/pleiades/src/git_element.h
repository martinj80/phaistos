// git_element.h --- Encapsulated GIT vector class for the clustering
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

#ifndef GIT_ELEMENT_H
#define GIT_ELEMENT_H

#include <vector>

#include "utils/vector_matrix_3d.h"
#include "protein/chain_ca.h"
#include "protein/chain_fb.h"

namespace phaistos {

//! Calculates Euclidean distance between two GIT vectors
inline double get_git_distance(std::vector<double> v1, std::vector<double> v2, double w1=1., double w2=1.) {

     if (v1.size() != v2.size()) {
          std::cerr << "WARNING: Vector size mismatch calculating the distance in cluster.cpp! \n\n";
          return 99999.;
     };
     double d = 0.;
     //
     // the distance is defined as the sqrt ( sum ( (v1[i]-v2[i])^2 ) )
     for (unsigned int i = 0; i < v1.size(); i++) {
          d += (v1[i] - v2[i]) * (v1[i] - v2[i]);
     }
     double tmp = sqrt(d);
     return tmp;
}

//! Calculates log-Euclidean distance between two GIT vectors
inline double get_git_log_distance(std::vector<double> v1, std::vector<double> v2) {
     return std::log(get_git_distance(v1, v2));
}



//! Representation of a GIT vector.
//! Allows to quickly access the GIT vector and accompanying attributes
//! such as the original filename, original chain length and so on.
//! Additionally it allows to store data relevant for clustering, such
//! as the weight or the distance to the mean structure.
class GitElement {

private:

     //! name description (ie the original filename)
     std::string name;

     //! Chain identifier
     int chain;
     
     //! Index to keep track of git elements when shuffling
     int index;

     //! Domain identifier
     int domain;

     //! Original protein chain length
     int length;

     //! The GIT vector itself
     std::vector<double> git;

     //! Distance to the mean structure of a cluster
     double distance_to_mean;

     //! The weight is used for the weighted clustering.
     double weight;


public:

     //! Standard constructor setting all the attributes available when
     //! reading in a standard GIT vector file.
     //!
     //! \param name allows to specify the name.
     //! \param chain specify the chain identifier.
     //! \param index Index to keep track of git elements when shuffling
     //! \param domain specify the domain identifier.
     //! \param length specify the original chain length.
     //! \param git the 30 dimensional vector itself.
     GitElement(std::string &name, int &chain, int &index, int &domain, int &length, std::vector<double> &git);

     //! Standard class Constructor
     GitElement() {};

     //! Standard class destructor
     ~GitElement() {};

     //! get name attribute
     std::string get_name();

     //! get chain attribute
     int get_chain();
     
     //! get index attribute
     int get_index();

     //! get domain attribute
     int get_domain();

     //! get length attribute
     int get_length();

     //! get entire GIT vector attribute
     std::vector<double> get_git();

     //! access GIT vector at a given index
     double get_git(int index);

     //! get weight attribute
     double get_weight();

     //! set weight attribute
     void set_weight(double weight) ;

     //! set name attribute (string)
     void set_name(std::string name) ;

     //! set name attribute (char*)
     void set_name(const char *name) ;

     //! set domain attribute
     void set_domain(int domain) ;

     //! set chain attribute
     void set_chain(int chain) ;
     
     //! get index attribute
     void set_index(int index);

     //! set length attribute
     void set_length(int length) ;

     //! set GIT vector (array)
     void set_git(double (&git)[31]);

     //! set GIT vector (vector)
     void set_git(std::vector<double> git);

     //! output vector the standard output stream
     void print();

     //! return distance to mean attribute
     double get_distance_to_mean() const;

     //! set distance to mean attribute
     void set_distance_to_mean(double distance);

     //! setting the distance to mean, calculated to the given mean structure
     void set_distance_to_mean(GitElement &mean);

     //! Comparison operator for two GitElements. The comparision
     //! is then made by comparing the distance to the mean
     //! for both contestants.
     //!
     //! \param a first GitElement
     //! \param b second GitElement
     friend inline bool operator<(const GitElement &a, const GitElement &b)  {
          return (a.get_distance_to_mean()<b.get_distance_to_mean());
     }
};

}

#endif
