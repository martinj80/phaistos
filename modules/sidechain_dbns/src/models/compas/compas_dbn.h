// compas_dbn.h --- Encapsulates all dbn handling for Compas
// Copyright (C) 2008-2010 Tim Harder
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


#ifndef COMPAS_DBN_H
#define COMPAS_DBN_H

#include "utils/vector_matrix_3d.h"
#include "protein/chain_fb.h"
#include "mocapy.h"
#include "utils/random.h"

namespace phaistos {

//! CompasDBN class.
//! This class is occludes all the complexity of accessing
//! the trained DBN model from the user. The model is
//! properly initialized and can then be accessed in a
//! variety off different ways.
class CompasDBN {
private:

     //! Map from amino acid indices to the reduced amino acid indices used in the model
     static const int aa_to_model_aa_array[];

     //! Map from the amino acid indices in the model to the normal amino acid indices
     static const int model_aa_to_aa_array[];

     //! Side chain mean vectors
     static const double aa_mean_vectors[20 * 3];
     
     //! Dynamic Bayesian Network model object     
     mocapy::DBN dbn;

     //! Path to DBN model filename     
     std::string dbn_filename;

     //! Internal copies of model - used for multithreading purposes     
     std::vector<CompasDBN*> copies;

     //! Thread ID of current copy     
     int thread_id;

     //! Random number engine from which random number generators can be created.     
     RandomNumberEngine *random_number_engine;

     //! Keeps track of last resamples residue     
     int res_cache;

     //! Mocapy sampler object
     mocapy::SampleInfEngineMM *sampler;

     //! Mocapy likelihood evaluator
     mocapy::LikelihoodInfEngineMM ll_engine;

     //@{
     //! Maintain mismasks for different scenarios, like
     //! sampling, inference and backbone independent sampling.
     mocapy::MDArray<mocapy::eMISMASK> *mm;
     mocapy::MDArray<mocapy::eMISMASK> *mm_ll;
     mocapy::MDArray<mocapy::eMISMASK> *mm_no_bb;
     //@}

     //! Get data - initialize Mocapy sequence object
     //!
     //! This method initializes a mocapy sequence object in the
     //! correct fashion and also sets all the variables and the
     //! sequence correctly. This method expects a Residue pointer
     //! to set the angles correctly.
     //!
     //! \param res residue pointer, to extract the pseudo side chain angles and the amino acid type and so on
     //! \param ignore_bb whether or not to include the backbone angles
     mocapy::Sequence get_data(ResidueFB *res, bool ignore_bb = false);

     //! Get data
     //!
     //! This method initializes a mocapy sequence object in the
     //! correct fashipm and also sets the index variables and the
     //! sequence correctly, so that only the current angle
     //! values need to be added.
     //!
     //! \param res amino acid type
     mocapy::Sequence get_data(definitions::ResidueEnum res);

     //! Get mismask for sampling
     //!
     //! This method initializes a mocapy mismask object that
     //! can be used for sampling new angles.
     mocapy::MDArray<mocapy::eMISMASK> get_mismask();

     //! Get mismask for likelihood calculations
     //!
     //! This method initializes a mocapy mismask object that
     //! can be used when evaluating the likelihood of a
     //! given set of angles.
     mocapy::MDArray<mocapy::eMISMASK> get_mismask_ll();

     //! Get mismask for sampling without backbone dependencies
     //!
     //! This method initializes a mocapy mismask object that
     //! can be used for sampling new angles. In this case
     //! backbone independent sampling will be used.
     mocapy::MDArray<mocapy::eMISMASK> get_mismask_no_bb();

     //! Initializer
     //!
     //! Set the correct initial values for all
     //! member variables.
     void init();

     //! Initialize sampler
     //!
     //! This method intitializes the Mocapy sampler object.
     void init_sampler();

     //! Translate amino acid label to internal model index
     //!
     //! Translate the normal amino acid label to the indices used
     //! in the compas model (since we are not modeling some
     //! amino acids like glycine or alanin, the indices differ
     //! slightly).
     int aa_to_model_aa(int aa);

     //! Translate internal model index to amino acid index
     //!
     //! Translate the model amino acid indices to the regular
     //! labels used usually (since we are not modeling some
     //! amino acids like glycine or alanin, the indeces differ
     //! slightly)
     int model_aa_to_aa(int model_aa);

     //! Get the distance in Angstrom
     //!
     //! The compas model uses discrete bins to model the
     //! the distance from the alpha carbon to the pseudo
     //! atom. This method allows to relate the discrete
     //! bin number to a real distance.
     double get_distance(int aa, int dist_bin);

public:

     //! Name
     static const std::string name;

     //! Standard Constructor
     //!
     //! \param dbn_file filepath to the dbn pickle
     //! \param seed integer to seed the random number engine.
     CompasDBN(std::string dbn_file, uint seed);

     //! Destructor
     ~CompasDBN();

     //! Simple copy constructor
     //!
     //! \param other CompasDBN instance to create a copy from
     CompasDBN(const CompasDBN &other);

     //! Regular copy constructor
     //!
     //! This copy constructor should be used in a thread safe
     //! environment. It allows setting and using the correct
     //! random number engine.
     //!
     //! \param other CompasDBN instance to create a copy from
     //! \param random_number_engine instance of the RandomNumberEngine collection
     //! \param thread_id ID of the current running thread, to grab the correct  random number engine
     CompasDBN(const CompasDBN &other, RandomNumberEngine *random_number_engine, int thread_id = 0);

     //! Create internal copies
     //!
     //! When run in a multithreaded environment, this method makes it
     //! possible to create a number of internal copies. Each copy
     //! will then exclusively be used in a specific thread and only
     //! used a specific random number engine to ensure
     //! reproducability.
     //!
     //! \param ncopies number of internal copies to create
     //! \param random_number_engines random number engine collection to draw engines from
     void duplicate(int ncopies, const std::vector<RandomNumberEngine *> &random_number_engines);

     //! Retrieve one of the internal BasiliskDBN copies
     //!
     //! After duplication, each thread request a specific copy of the
     //! dbn object to continue working on, where index should be
     //! the thread_id to ensure only one thread is using a certain
     //! object.
     //!
     //! \param index Index of the internal copy to use (should be equal to the thread ID)
     CompasDBN &get_copy(unsigned int index);


     //! Generate a new sequence of angles for a specific residue in the chain
     //!
     //! Returns a newly sampled set of sidechain angles for a given
     //! amino acid flag. It is possible to condition sampling on the
     //! backbone angles of the residue. If left out backbone
     //! independent sampling will be used.
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     std::vector<double> get_angle_sequence_for_residue_with_bb(definitions::ResidueEnum res, 
                                                                double phi = UNINITIALIZED, 
                                                                double psi = UNINITIALIZED);

     //! Calculate the likelihood - for use immediately after sampling
     //!
     //! Returns the LL of the current network (only useful after sampling)
     //! This method returns the likelihood of the currently set sampler
     //! with the currently set values. This is useful directly after
     //! sampling, to evalute the likeihood.
     //!
     //! \param ignore_bb include the backbone angle information
     double calc_ll(bool ignore_bb = false);

     //! Calculates the likelihood
     //!
     //! Calculates the likelihood for a given set of angles, possibly
     //! including the backbone angles in the evaluation process.
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     //! \param dof dof in the side chain
     double calc_ll(int res, double phi, double psi, std::vector<double> dof);

     //! Calculates the likelihood
     //!
     //! Calculates the likelihood for a given set of angles, possibly
     //! including the backbone angles in the evaluation process.
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     //! \param dof dof in the side chain
     double calc_ll(definitions::ResidueEnum res, 
                    double phi, double psi, 
                    std::vector<double> dof);

     //! Get atom selection
     //!
     //! Returns the proper set of atoms required for the Compas model
     //! to work.
     definitions::AtomSelectionEnum get_atom_selection() {
          return definitions::PSEUDO_SIDECHAIN_ATOMS;
     }

     //! print
     //!
     //! A print function.
     //! Mostly useful in debugging scenarios.
     void print(int id = 0) {
          std::cout << "@@@@  ThreadID : " << this->thread_id << " \t debug id : " << id << "\t trying to access a node : \n" << *(this->dbn.getNodes1()[2])
                    << "\n";
          std::cout.flush();
     }

     //! printc
     //!
     //! A const version of print function.
     //! Mostly useful in debugging scenarios.
     void printc(int id = 0) const {
          std::cout << "@@@@  ThreadID : " << this->thread_id << " \t debug id : " << id << "\t trying to access a node : \n";
          std::cout.flush();
     }
};

}

#endif
