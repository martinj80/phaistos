// basilisk_dbn.h --- Encapsulates all dbn handling for Basilisk
// Copyright (C) 2010 Tim Harder
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

#ifndef BASILISK_DBN_H
#define BASILISK_DBN_H

#include "protein/chain_fb.h"
#include "utils/vector_matrix_3d.h"
#include "mocapy.h"
#include "utils/random.h"

namespace phaistos {

//! BasiliskDBN class
//! This class occludes all the complexity of accessing
//! the trained DBN model from the user. The model is
//! properly initialized and can then be accessed in a
//! variety off different ways.
class BasiliskDBN {
private:

     //! map from (aa,angle index) -> dbn index
     std::vector<std::vector<int> > aa_angle_to_index_matrix;

     //! Dynamic Bayesian Network model object
     mocapy::DBN dbn;

     //! Number of chi angles for all residue types
     int chi_lengths[20];

     //! Path to DBN model filename
     std::string dbn_filename;

     //! Internal copies of model - used for multithreading purposes
     std::vector<BasiliskDBN*> copies;

     //! Thread ID of current copy
     int thread_id;

     //! Random number engine from which random number generators can be created.
     RandomNumberEngine *random_number_engine;

     //! Keeps track of last resamples residue
     int res_cache;

     //! Keeps track of last phi angle
     double phi_cache;

     //! Keeps track of last psi angle
     double psi_cache;

     //! Keeps track of last chi angles     
     std::vector<double> chis_cache;

     //@{
     //! Maintain a SampleEngine for every possible
     //! chi length. This avoids expensive
     //! re-initializations when resetting.
     mocapy::SampleInfEngineHMM *sampler_len_1;
     mocapy::SampleInfEngineHMM *sampler_len_2;
     mocapy::SampleInfEngineHMM *sampler_len_3;
     mocapy::SampleInfEngineHMM *sampler_len_4;
     //@}

     //! Pointer to last SampleEngine used
     mocapy::SampleInfEngineHMM *sampler;

     //! Engine used for calculating likelihood values
     mocapy::LikelihoodInfEngineHMM ll_engine;

     //! Initialize the aa_angle_to_index_matrix map data structure
     void initialize_aa_angle_to_index_matrix();     

     //! Amino acid angle to index.
     //!
     //! Helper function, coverting the amino acid type
     //! plus the chi angle to be modeled into the correct
     //! index to be used in the dbn.
     int aa_angle_to_index(uint aa, uint chi_index);

     //! Get Chi length.
     //!
     //! Helper method returning the number of freely rotatable
     //! chi angles for a given amino acid label.
     int get_chi_len(int aa);

     //! Set sampler
     //!
     //! This method allows setting the current data and mismask
     //! into the sampler engine and also selects the correct
     //! sampler engine for the current amino acid (actually the
     //! chi length is the trigger here).
     //!
     //! \param data correctly formatted data object to be set in the sampler
     //! \param mismask correctly formatted mismask to be set in the sampler.
     //!
     void set_sampler(mocapy::Sequence &data, mocapy::MDArray<mocapy::eMISMASK> &mismask);

     //! Init function
     //!
     //! Sets the correct initial values for all
     //! member variables.
     void init();

     //! Init Chi Array
     //!
     //! Initializes an array to easily convert
     //! between the aminoacid type and the number of
     //! chi angles in the side chain.
     void init_chi_array();

     //! Get data method
     //!
     //! Initializes a mocapy sequence object in the correct length
     //! and also sets the index variables and the sequence length
     //! correctly, so that only the current angle values need to be
     //! added.
     //!
     //! \param res amino acid type flag
     mocapy::Sequence get_data(int res);

     //! Get mismask for sampling
     //!
     //! This method initializes a mocapy mismask object that
     //! can be used for sampling new angles.
     //!
     //! \param res amino acid type flag
     mocapy::MDArray<mocapy::eMISMASK> get_mismask(int res);

     //! Get mismask for likelihood evaluation
     //!
     //! This method initializes a mocapy mismask object that
     //! can be used when evaluating the likelihood of a
     //! given set of angles.
     //!
     //! \param res amino acid type flag
     mocapy::MDArray<mocapy::eMISMASK> get_mismask_ll(int res);

public:

     //! Name
     static const std::string name;

     //! Standard Constructor
     //!
     //! this is the standard constructor that should be used
     //!
     //! \param dbn_file filepath to the dbn pickle
     //! \param seed integer to seed the random number engine.
     BasiliskDBN(std::string dbn_file, uint seed = 666);


     //! Destructor
     ~BasiliskDBN();

     //! Simple copy constructor
     //!
     //! \param other BasiliskDBN instance to create a copy from
     BasiliskDBN(const BasiliskDBN &other);

     //! Regular copy constructor
     //!
     //! This copy constructor should be used in a thread safe
     //! environment. It allows setting and using the correct
     //! random number engine.
     //!
     //! \param other BasiliskDBN instance to create a copy from
     //! \param random_number_engine instance of the RandomNumberEngine collection
     //! \param thread_id ID of the current running thread, to grab the correct  random number engine
     BasiliskDBN(const BasiliskDBN &other, RandomNumberEngine *random_number_engine, int thread_id = 0);

     //! Create internal copies
     //!
     //! when run in a multithreaded environment, this method allows
     //! to create a number of internal copies. Each copy will then
     //! exclusively be used in a specific thread and only used a
     //! specific random number engine to ensure reproducability.
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
     BasiliskDBN &get_copy(unsigned int index);


     //! Generate a new sequence of angles for a specific residue in the chain
     //!
     //! Returns a newly sampled set of sidechain angles
     //! for a given amino acid flag. It is possible to condition
     //! sampling on the backbone angles of the residue. If left out
     //! backbone independent sampling will be used.
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     std::vector<double> get_angle_sequence_for_residue_with_bb(int res, 
                                                                double phi = UNINITIALIZED, 
                                                                double psi = UNINITIALIZED);

     //! Generate a new sequence of angles for a specific residue in the chain
     //!
     //! Returns a newly sampled set of sidechain angles
     //! for a given amino acid flag. It is possible to condition
     //! sampling on the backbone angles of the residue. If left out
     //! backbone independent sampling will be used.
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     std::vector<double> get_angle_sequence_for_residue_with_bb(definitions::ResidueEnum res, 
                                                                double phi = UNINITIALIZED, 
                                                                double psi = UNINITIALIZED);

     //! Generate a new sequence of angles for a specific residue in the chain - without resampling hidden nodes.
     //!
     //! Returns a newly sampled set of side chain angles
     //! for a given amino acid flag. It is possible to condition
     //! sampling on the backbone angles of the residue. If left out
     //! backbone independent sampling will be used.
     //!
     //! This method allows resampling the output nodes
     //! without resampling the hidden states. This will
     //! effectively lead to only modest changes in the
     //! angle.
     //!
     //! \param res amino acid type flag
     //! \param chis_old vector of current chi angles to condition sampling upon
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     std::vector<double> get_local_angle_sequence_for_residue(int res, std::
                                                              vector<double> chis_old, 
                                                              double phi = UNINITIALIZED, 
                                                              double psi = UNINITIALIZED) ;

     //! Generate a new sequence of angles for a specific residue in the chain - without resampling hidden nodes.
     //!
     //! Returns a newly sampled set of side chain angles
     //! for a given amino acid flag. It is possible to condition
     //! sampling on the backbone angles of the residue. If left out
     //! backbone independent sampling will be used.
     //!
     //! This method allows resampling the output nodes
     //! without resampling the hidden states. This will
     //! effectively lead to only modest changes in the
     //! angle.
     //!
     //! \param res amino acid type flag
     //! \param chis_old vector of current chi angles to condition sampling upon
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     std::vector<double> get_local_angle_sequence_for_residue(definitions::ResidueEnum res, 
                                                              std::vector<double> chis_old, 
                                                              double phi = UNINITIALIZED, 
                                                              double psi = UNINITIALIZED);

     //! Calculate the likelihood - for use immediately after sampling
     //!
     //! Calculates the likelihood for a given set of chi
     //! angles, possibly including the backbone angles in the
     //! evaluation process.
     //! respect_symmetry allows to include symmetry considerations
     //! in the side chain (ie in the phenylalanin ring).
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     //! \param chis chi side chain angle
     //! \param respect_symmetry whether Symmetry
     double calc_ll(int res, double phi, double psi, 
                    std::vector<double> chis, 
                    bool respect_symmetry = false);

     //! Calculates the likelihood
     //!
     //! Calculates the likelihood for a given set of chi
     //! angles, possibly including the backbone angles in the
     //! evaluation process.
     //! respect_symmetry allows to include symmetry considerations
     //! in the side chain (ie in the phenylalanin ring).
     //!
     //! \param res amino acid type flag
     //! \param phi phi backbone angle (if not set, the backbone angles will be ignored)
     //! \param psi psi backbone angle (if not set, the backbone angles will be ignored)
     //! \param chis chi side chain angle
     //! \param respect_symmetry whether Symmetry
     double calc_ll(definitions::ResidueEnum res, double phi, double psi, 
                    std::vector<double> chis, bool respect_symmetry = false);

     //! Calculates the likelihood
     //!
     //! Returns the LL of the current network (only useful after sampling)
     //! This method returns the likelihood of the currently set sampler
     //! with the currently set values. This is useful directly after
     //! sampling, to evalute the likeihood.
     //!
     //! \param ignore_bb include the backbone angle information
     //!
     double calc_ll(bool ignore_bb = true);

     //! Get atom selection
     //!
     //! Returns the proper set of atoms required for the Basilisk
     //! model to work.
     definitions::AtomSelectionEnum get_atom_selection() {
          return definitions::SIDECHAIN_ATOMS;
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
