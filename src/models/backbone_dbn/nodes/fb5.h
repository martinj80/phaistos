// fb5.h --- Fb5 Node. Density and sampler for the 5-parameter Fisher-Bingham distribution
// Copyright (C) 2006-2008 Wouter Boomsma, Thomas Hamelryck
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

#ifndef FB5_H
#define FB5_H

#include "utils/random.h"
#include "utils/vector_matrix_3d.h"
#include "utils/vector_nd.h"
#include "parameters.h"
#include "node.h"

namespace phaistos {

//! Density class for the 5-parameter Fisher-Bingham distribution
class Fb5Density {
private:
     //! Unit vector 1
     Vector_3D *e1;

     //! Unit vector 2
     Vector_3D *e2;

     //! Unit vector 3
     Vector_3D *e3;

     //! kappa parameter
     double kappa;

     //! beta parameter
     double beta;

     //! log-normalization constant
     double logc;

public:
     //! Constructor
     //! \param kappa Kappa parameter
     //! \param beta Beta parameter
     //! \param e1 e1 unit vector
     //! \param e2 e1 unit vector
     //! \param e3 e1 unit vector
     Fb5Density(double kappa, double beta,
		std::vector<double> &e1, 
                std::vector<double> &e2, 
                std::vector<double> &e3);

     //! Copy constructor
     //! \param other Source object
     Fb5Density(const Fb5Density &other);

     //! Destructor
     ~Fb5Density();

     //! Calculate likelihood
     //! \param angle_pair observed value
     //! \return likelihood value     
     inline double get_likelihood(double *angle_pair) const {
          return std::exp(get_log_likelihood(angle_pair));
     }

     //! Calculate log-likelihood
     //! \param angle_pair observed value
     //! \return log-likelihood value     
     inline double get_log_likelihood(double *angle_pair) const {
          double theta = angle_pair[0];
          double tau = angle_pair[1];

          Vector_3D v;
          Fb5Density::polar_to_xyz(theta, tau, &v);
     
          // Return log likelihood of the vector v
          double a, b, c, d;
     
          a=kappa*(v*(*e1));
          b=v*(*e2);
          c=v*(*e3);
     
          d=beta*(b*b-c*c);

          double dens = -logc+a+d;
          return dens;
     }

     //! Convert from cartesian to polar coordinates
     //! \param x cartesian coordinate
     //! \param tt target polar coordinates (two angles)
     static void xyz_to_polar(Vector_3D *x, double *tt);

     //! Convert from polar to cartesian coordinates
     //! \param theta theta angle
     //! \param tau tau angle
     //! \param v target cartesian coordinate
     static void polar_to_xyz(double theta, double tau, Vector_3D *v);
};

//! Sampler class for the 5-parameter Fisher-Bingham distribution
class Fb5Sampler {
private:

     //! Random number generator
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > variate_generator_uniform;

     //@{
     //! Internal parameters
     double a, b, c2, gamma, lam1, lam2;
     Matrix_3D* rho;
     //@}
public:
     
     //! Constructor
     //! \param kappa Kappa parameter
     //! \param beta Beta parameter
     //! \param e1 e1 unit vector
     //! \param e2 e1 unit vector
     //! \param e3 e1 unit vector
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     Fb5Sampler(double kappa, double beta,
		std::vector<double> &e1, 
                std::vector<double> &e2, 
                std::vector<double> &e3,
                RandomNumberEngine *random_number_engine);

     //! Copy constructor
     //! \param other Source object     
     Fb5Sampler(const Fb5Sampler &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     Fb5Sampler(const Fb5Sampler &other, 
                RandomNumberEngine *random_number_engine);

     //! Destructor     
     ~Fb5Sampler() { delete rho; };

     //! Sample value
     //! \return sampled value     
     Vector_nD sample();
};

//! Distribution class for the 5-parameter Fisher-Bingham distribution. Container of density and sampler objects
class Fb5Distribution {
public:
     //! Density object
     Fb5Density *density;

     //! Sampler object
     Fb5Sampler *sampler;
     
     //! Constructor
     //! \param kappa Kappa parameter
     //! \param beta Beta parameter
     //! \param e1 e1 unit vector
     //! \param e2 e1 unit vector
     //! \param e3 e1 unit vector
     //! \param random_number_engine Random number engine from which random number generators can be constructed.     
     Fb5Distribution(double kappa, double beta,
                     std::vector<double> &e1, 
                     std::vector<double> &e2, 
                     std::vector<double> &e3,
                     RandomNumberEngine *random_number_engine = &random_global);

     //! Copy constructor
     //! \param other Source object
     Fb5Distribution(const Fb5Distribution &other);

     //! Copy constructor - different random number engine
     //! \param other Source object
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     Fb5Distribution(const Fb5Distribution &other,
                     RandomNumberEngine *random_number_engine);

     //! Destructor     
     ~Fb5Distribution();

     //! Sample value
     //! \return value
     Vector_nD sample();

     //! Calculate log likelihood
     //! \param angle_pair observed value
     //! \return log-likelihood value     
     double get_log_likelihood(double *angle_pair);
};


//! Node class for the 5-parameter Fisher-Bingham distribution.
template <typename DBN_TYPE, int NODE_INDEX, int PARENT_NODE=0>
class Fb5Node: public Node_2D<double,DBN_TYPE,PARENT_NODE> {
private:

     //! An array of distribution components
     Fb5Distribution **data;

public:

     //! \param name Node name
     //! \param parameters Parameter object
     //! \param dbn Pointer to dbn in which this node is contained
     //! \param random_number_engine Random number engine from which random number generators can be constructed.
     Fb5Node(std::string name, Parameters &parameters, 
             DBN_TYPE *dbn,
             RandomNumberEngine *random_number_engine = &random_global)
          : Node_2D<double,DBN_TYPE,PARENT_NODE>(name, 2, dbn) {

	  // Extract parameters
	  std::vector<double> kappa_list = parameters.get_node_parameters<std::vector<double> > ("FB5", name, "KAPPA");
	  std::vector<double> beta_list = parameters.get_node_parameters<std::vector<double> > ("FB5", name, "BETA");
	  std::vector<std::vector<double> > e1_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("FB5", name, "E1");
	  std::vector<std::vector<double> > e2_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("FB5", name, "E2");
	  std::vector<std::vector<double> > e3_list = parameters.get_node_parameters<std::vector<std::vector<double> > > ("FB5", name, "E3");

	  this->h_size = kappa_list.size();
     
	  data = new Fb5Distribution*[this->h_size];

	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new Fb5Distribution(kappa_list[i], beta_list[i],
                                             e1_list[i], e2_list[i], e3_list[i],
                                             random_number_engine);
	  }
     }
	  

     //! Copy constructor
     //! \param other Source object
     Fb5Node(const Fb5Node &other)
          : Node_2D<double,DBN_TYPE>(other) {

	  data = new Fb5Distribution*[other.h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new Fb5Distribution(*other.data[i],
                                             other->dbn->random_number_engine);
	  }
     }


     //! Copy constructor - different dbn object.
     //! \param other Source object.
     //! \param dbn Pointer to a dbn object.
     Fb5Node(const Fb5Node &other, DBN_TYPE *dbn)
          : Node_2D<double,DBN_TYPE>(other, dbn) {

	  data = new Fb5Distribution*[other.h_size];
	  
	  for (int i=0; i<this->h_size; i++) {
	       data[i] = new Fb5Distribution(*other.data[i],
                                             dbn->random_number_engine);
	  }
     }

     //! Destructor
     ~Fb5Node() {
	  for (int i=0; i<this->h_size; i++) {
	       delete data[i];
	  }
	  delete[] data;
     }
     
     
     //! Overload indexing operator 
     //! \param index distribution component index
     //! \return distribution component
     Fb5Distribution *operator[](const int index) const {
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
	  for (int h=0; h<this->h_size; h++) {
	       std::vector<double> mu;
	       Vector_3D p((*data[h]->density->e1)[0],
			   (*data[h]->density->e1)[1],
			   (*data[h]->density->e1)[2]);
	       double tt[2];
	       Fb5Density::xyz_to_polar(&p, tt);
	       mu.push_back(tt[0]);
	       mu.push_back(tt[1]);
	       means.push_back(mu);
	  }
	  return means;
     }

     
     //! Return kappa,beta,e1,e2,e3 parameters for all distributions components
     //! \return vector of kappa,beta,e1,e2,e3 values for each component
     std::vector<std::vector<double> > get_shape_parameters() {
	  std::vector<std::vector<double> > shape_parameters_list;
	  for (int h=0; h<this->h_size; h++) {
	       std::vector<double> shape_parameters;
	       shape_parameters.push_back(data[h]->density->kappa);
	       shape_parameters.push_back(data[h]->density->beta);	  
	       shape_parameters.push_back((*data[h]->density->e1)[0]);	  
	       shape_parameters.push_back((*data[h]->density->e1)[1]);	  
	       shape_parameters.push_back((*data[h]->density->e1)[2]);	  
	       shape_parameters.push_back((*data[h]->density->e2)[0]);	  
	       shape_parameters.push_back((*data[h]->density->e2)[1]);	  
	       shape_parameters.push_back((*data[h]->density->e2)[2]);	  
	       shape_parameters.push_back((*data[h]->density->e3)[0]);	  
	       shape_parameters.push_back((*data[h]->density->e3)[1]);	  
	       shape_parameters.push_back((*data[h]->density->e3)[2]);	  
	       shape_parameters_list.push_back(shape_parameters);
	  }
	  return shape_parameters_list;
     }
};

}

#endif
