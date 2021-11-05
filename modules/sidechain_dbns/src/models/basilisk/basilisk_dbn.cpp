// basilisk_dbn.cpp --- Encapsulates all dbn handling for Basilisk
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

#include "models/basilisk/basilisk_dbn.h"
#include "mocapy.h"

#include <stdio.h>
#include "utils/random.h"

#include <fstream>
#include <math.h>

namespace phaistos {

//! Name
const std::string BasiliskDBN::name = "basilisk";

//! Initialize the (aa,chi_index) -> dbn index data structure
void BasiliskDBN::initialize_aa_angle_to_index_matrix() {

     aa_angle_to_index_matrix[1][1] = 2;
     aa_angle_to_index_matrix[2][1] = 3;
     aa_angle_to_index_matrix[2][2] = 4;
     aa_angle_to_index_matrix[3][1] = 5;
     aa_angle_to_index_matrix[3][2] = 6;
     aa_angle_to_index_matrix[3][3] = 7;
     aa_angle_to_index_matrix[4][1] = 8;
     aa_angle_to_index_matrix[4][2] = 9;
     aa_angle_to_index_matrix[6][1] = 10;
     aa_angle_to_index_matrix[6][2] = 11;
     aa_angle_to_index_matrix[7][1] = 12;
     aa_angle_to_index_matrix[7][2] = 13;
     aa_angle_to_index_matrix[8][1] = 14;
     aa_angle_to_index_matrix[8][2] = 15;
     aa_angle_to_index_matrix[8][3] = 16;
     aa_angle_to_index_matrix[8][4] = 17;
     aa_angle_to_index_matrix[9][1] = 18;
     aa_angle_to_index_matrix[9][2] = 19;
     aa_angle_to_index_matrix[10][1] = 20;
     aa_angle_to_index_matrix[10][2] = 21;
     aa_angle_to_index_matrix[10][3] = 22;
     aa_angle_to_index_matrix[11][1] = 23;
     aa_angle_to_index_matrix[11][2] = 24;
     aa_angle_to_index_matrix[12][1] = 25;
     aa_angle_to_index_matrix[12][2] = 26;
     aa_angle_to_index_matrix[13][1] = 27;
     aa_angle_to_index_matrix[13][2] = 28;
     aa_angle_to_index_matrix[13][3] = 29;
     aa_angle_to_index_matrix[14][1] = 30;
     aa_angle_to_index_matrix[14][2] = 31;
     aa_angle_to_index_matrix[14][3] = 32;
     aa_angle_to_index_matrix[14][4] = 33;
     aa_angle_to_index_matrix[15][1] = 34;
     aa_angle_to_index_matrix[16][1] = 35;
     aa_angle_to_index_matrix[17][1] = 36;
     aa_angle_to_index_matrix[18][1] = 37;
     aa_angle_to_index_matrix[18][2] = 38;
     aa_angle_to_index_matrix[19][1] = 39;
     aa_angle_to_index_matrix[19][2] = 40;
}

//! Lookup dbn index based on amino acid and chi_index
int BasiliskDBN::aa_angle_to_index(uint aa, uint chi_index) {
     return aa_angle_to_index_matrix[aa][chi_index];
}

//! Return the number of chi angles in a given amino acid
int BasiliskDBN::get_chi_len(int aa) {
     return this->chi_lengths[aa];
}

//! Initializes a mocapy sequence object
mocapy::Sequence BasiliskDBN::get_data(int res) {
     mocapy::Sequence data;
     // get and check number of chi angles

     int len = get_chi_len(res);
     if (len == 0)
          return data;

     // create the array
     data.set_shape(2 + len, 3);
     data.set_wildcard(mocapy::vec(-1, 1), 0);
     data.set_wildcard(mocapy::vec(-1, 2), 0);

     /////////
     // setting all the indeces properly
     // Phi & Psi
     data.set(0, 0, 0.);
     data.set(1, 0, 1.);

     // And all the chi+aminoacid indeces
     data.set(2, 0, aa_angle_to_index(res, 1)); // x1

     if (len > 1) {
          data.set(3, 0, aa_angle_to_index(res, 2)); // x2

          if (len > 2) {
               data.set(4, 0, aa_angle_to_index(res, 3)); // x3

               if (len > 3) {
                    data.set(5, 0, aa_angle_to_index(res, 4)); // x4
               }
          }
     }

     return data;
}

//! Get mismask for sampling
mocapy::MDArray<mocapy::eMISMASK> BasiliskDBN::get_mismask(int res) {
     mocapy::MDArray<mocapy::eMISMASK> mism_sample;
     // get and check number of chi angle
     int len = get_chi_len(res);
     if (len == 0)
          return mism_sample;

     mism_sample.set_shape(2 + len, 3);
     mism_sample.set_wildcard(mocapy::vec(-1, 0), mocapy::MOCAPY_OBSERVED);
     mism_sample.set_wildcard(mocapy::vec(-1, 1), mocapy::MOCAPY_HIDDEN);
     mism_sample.set_wildcard(mocapy::vec(-1, 2), mocapy::MOCAPY_HIDDEN);

     return mism_sample;
}

//! Get mismask for likelihood evaluation
mocapy::MDArray<mocapy::eMISMASK> BasiliskDBN::get_mismask_ll(int res) {
     mocapy::MDArray<mocapy::eMISMASK> mism_sample;
     // get and check number of chi angle
     int len = get_chi_len(res);
     if (len == 0)
          return mism_sample;

     mism_sample.set_shape(2 + len, 3);
     mism_sample.set_wildcard(mocapy::vec(-1, 0), mocapy::MOCAPY_OBSERVED);
     mism_sample.set_wildcard(mocapy::vec(-1, 1), mocapy::MOCAPY_HIDDEN);
     mism_sample.set_wildcard(mocapy::vec(-1, 2), mocapy::MOCAPY_OBSERVED);

     return mism_sample;
}

//! Initializer
void BasiliskDBN::init() {

     this->sampler = NULL;
     this->sampler_len_1 = NULL;
     this->sampler_len_2 = NULL;
     this->sampler_len_3 = NULL;
     this->sampler_len_4 = NULL;

     this->res_cache = -1;
     this->phi_cache = UNINITIALIZED;
     this->psi_cache = UNINITIALIZED;

     this->random_number_engine = NULL;
     this->thread_id = 0;

     this->init_chi_array();
}

//! Initialized chi_array data structure
void BasiliskDBN::init_chi_array() {
     this->chi_lengths[0] = 0;
     this->chi_lengths[1] = 1;
     this->chi_lengths[2] = 2;
     this->chi_lengths[3] = 3;
     this->chi_lengths[4] = 2;
     this->chi_lengths[5] = 0;
     this->chi_lengths[6] = 2;
     this->chi_lengths[7] = 2;
     this->chi_lengths[8] = 4;
     this->chi_lengths[9] = 2;
     this->chi_lengths[10] = 3;
     this->chi_lengths[11] = 2;
     this->chi_lengths[12] = 2;
     this->chi_lengths[13] = 3;
     this->chi_lengths[14] = 4;
     this->chi_lengths[15] = 1;
     this->chi_lengths[16] = 1;
     this->chi_lengths[17] = 1;
     this->chi_lengths[18] = 2;
     this->chi_lengths[19] = 2;
}

//! Destructor
BasiliskDBN::~BasiliskDBN() {
     this->sampler = NULL;
     if (this->sampler_len_1 != NULL) {
          delete this->sampler_len_1;
          this->sampler_len_1 = NULL;
     }
     if (this->sampler_len_2 != NULL) {
          delete this->sampler_len_2;
          this->sampler_len_2 = NULL;
     }
     if (this->sampler_len_3 != NULL) {
          delete this->sampler_len_3;
          this->sampler_len_3 = NULL;
     }
     if (this->sampler_len_4 != NULL) {
          delete this->sampler_len_4;
          this->sampler_len_4 = NULL;
     }
     if (this->random_number_engine != NULL) {
          delete this->random_number_engine;
          this->random_number_engine = NULL;
     }
}

//! Constructor
BasiliskDBN::BasiliskDBN(std::string dbn_file, uint seed)
     : aa_angle_to_index_matrix(20, std::vector<int>(5, -1)) {

     initialize_aa_angle_to_index_matrix();

     // just to make sure that the old default does still work.
     if (dbn_file.substr((dbn_file.length()-27)) == "/src/corona/data/corona.dbn") {
          std::cerr << "################################################################################################# \n";
          std::cerr << "# DEPRECATION WARNING: You are trying to access the deprecated corona.dbn. \n";
          std::cerr << "# Please switch to the latest version and use basilisk.dbn to make this warning disappear. \n";
          std::cerr << "################################################################################################# \n";
          dbn_file = dbn_file.substr(0,(dbn_file.length()-27))+"/src/basilisk/data/basilisk.dbn";
     }

     assert(file_exists(dbn_file));
     this->init();
     this->init_chi_array();

     this->dbn.load(dbn_file.c_str());
     this->dbn.seed(seed);

     // store this in case we need to take copies
     this->dbn_filename = dbn_file;

     mocapy::LikelihoodInfEngineHMM infengine(&dbn, 1);
     this->ll_engine = infengine;
}

//! Simple copy constructor
BasiliskDBN::BasiliskDBN(const BasiliskDBN &other)
     : aa_angle_to_index_matrix(other.aa_angle_to_index_matrix) {

     // this->chain = other.chain;

     this->init_chi_array();

     // copy and reseed the DBN
     this->dbn = other.dbn;

     this->dbn_filename = other.dbn_filename;

     this->copies = other.copies;
     this->thread_id = other.thread_id;
     this->random_number_engine = other.random_number_engine;

     this->res_cache = other.res_cache;

     this->sampler_len_1 = other.sampler_len_1;
     this->sampler_len_2 = other.sampler_len_2;
     this->sampler_len_3 = other.sampler_len_3;
     this->sampler_len_4 = other.sampler_len_4;
     this->sampler = other.sampler;
     this->ll_engine = other.ll_engine;

}

//! Regular copy constructor
BasiliskDBN::BasiliskDBN(const BasiliskDBN &other, RandomNumberEngine *random_number_engine, int thread_id)
     : aa_angle_to_index_matrix(other.aa_angle_to_index_matrix) {

     int seed = (*random_number_engine)();

     // should really be Null ..
     // this->chain = other.chain;

     this->init();
     this->init_chi_array();

     // copy and reseed the DBN
     this->dbn = *(new mocapy::DBN);
     this->dbn.load(other.dbn_filename.c_str());
     this->dbn.seed(seed);

     mocapy::LikelihoodInfEngineHMM infengine(&(this->dbn), 1);
     this->ll_engine = infengine;

     this->res_cache = 0;
     this->dbn_filename = other.dbn_filename;
     this->thread_id = thread_id;
     this->copies = other.copies;
     this->random_number_engine = other.random_number_engine;
}

//! Duplicate DBN for multithreading
void BasiliskDBN::duplicate(int ncopies, const std::vector<RandomNumberEngine *> &random_number_engines) {
     // Resize threads vector
     assert(copies.size()==0);
     copies.resize(ncopies);

     if (thread_id == 0) {
          // Use the passed dbn as first thread
          copies[thread_id] = this;

          // Make nthreads-1 copies
          for (int i = 1; i < ncopies; i++) {

               if (random_number_engines.size() > 0) {
                    copies[i] = new BasiliskDBN(*this, random_number_engines[i], i);
               } else {
                    copies[i] = new BasiliskDBN(*this, this->random_number_engine, i);
               }
          }

          // Make vectors in all copies consistent
          for (int i = 1; i < ncopies; i++) {
               for (int j = 0; j < ncopies; j++) {
                    copies[i]->copies[j] = copies[j];
               }
          }
     }
}

//! Retrieve one of the internal BasiliskDBN copies
BasiliskDBN &BasiliskDBN::get_copy(unsigned int index) {
     if (copies.size() > index) {
          return *copies[index];
     } else if (index == 0) {
          return *this;
     } else {
          assert(false);
     }
}


//! Generate a new sequence of angles for a specific residue in the chain - without resampling hidden nodes.
std::vector<double> BasiliskDBN::get_local_angle_sequence_for_residue(int res, std::vector<double> chis_old, double phi, double psi) {
     return this->get_local_angle_sequence_for_residue((definitions::ResidueEnum)res,chis_old, phi, psi);
}


//! Generate a new sequence of angles for a specific residue in the chain - without resampling hidden nodes.
std::vector<double> BasiliskDBN::get_local_angle_sequence_for_residue(definitions::ResidueEnum res, std::vector<double> chis_old, 
                                                                      double phi, double psi) { 

     // For debugging purposes. Can be removed.
     assert((!is_initialized(phi) || phi > -4.) && (!is_initialized(psi) || psi > -4.));

     // Import protein definitions (such as residue names)
     using namespace definitions;

     std::vector<double> chis;

     this->res_cache = res;
     this->phi_cache = phi;
     this->psi_cache = psi;

     // Do nothing for amino acids without side chains
     if (res == ALA || res == GLY || chis_old.size() < 1) {
          this->res_cache = -1;
          this->chis_cache = chis;
          return chis;
     }

     // Correct for alternative definition of HIS chi2
     if (res==HIS) {
          chis_old[1] = fmod(chis_old[1] + 2*M_PI, 2*M_PI) - M_PI;
     }

     // setting up the data and mismask
     // array properly
     mocapy::Sequence data = this->get_data(res);
     mocapy::MDArray<mocapy::eMISMASK> mism_sample = this->get_mismask(res);
     //
     // lets see whether we got backbone angles as well
     double tmp = phi;
     if (is_initialized(tmp)) {
          // not at the terminus
          data.set(0, 2, tmp);
          mism_sample.set(0, 2, mocapy::MOCAPY_OBSERVED);
     }
     tmp = psi;
     if (is_initialized(tmp)) {
          // not at the terminus
          data.set(1, 2, tmp);
          mism_sample.set(1, 2, mocapy::MOCAPY_OBSERVED);
     }
     //
     // set the chi values properly and mark them observed
     // for (int i = 0; i < (int) chis_old.size(); i++) {
     for (int i = 0; i < (int)data.get_shape()[0]-2; i++) {
          data.set(i + 2, 2, chis_old[i]);
          mism_sample.set(i + 2, 2, mocapy::MOCAPY_OBSERVED);
     }

     // Setup the proper sampler
     this->set_sampler(data, mism_sample);
     //
     // calculate the hidden node sequence
     mocapy::MDArray<double> *sample;
     sample = &(this->sampler->sample_next());

     //
     mocapy::Node* node = this->dbn.getNodes1()[2];
     //
     // resampling the outputnodes only, but
     // keeping the hidden node states fixed.
     // for (int i = 0; i < (int) chis_old.size(); i++) {
     for (int i = 0; i < (int)data.get_shape()[0]-2; i++) {
          (*node).sample(i + 2);
     }
     //
     // storing the local chis.
     for (int i = 2; i < (int) (*sample).size() / 3; i++) {
          chis.push_back((*sample).get(i, 2));
     }

     // Correct for alternative definition of HIS chi2
     if (res==HIS) {
          chis[1] = fmod(chis[1] + 2*M_PI, 2*M_PI) - M_PI;
     }

     // Save chi values
     this->chis_cache = chis;

     return chis;
}


//! Generate a new sequence of angles for a specific residue in the chain
std::vector<double> BasiliskDBN::get_angle_sequence_for_residue_with_bb(int res, double phi, double psi) {
     return this->get_angle_sequence_for_residue_with_bb((definitions::ResidueEnum) res, phi, psi);
}

//! Generate a new sequence of angles for a specific residue in the chain
std::vector<double> BasiliskDBN::get_angle_sequence_for_residue_with_bb(definitions::ResidueEnum res, double phi, double psi) { //

     // For debugging purposes. Can be removed.
     assert((!is_initialized(phi) || phi > -4.) && (!is_initialized(psi) || psi > -4.));

     // Import protein definitions (such as residue names)
     using namespace definitions;

     std::vector<double> chis;

     this->res_cache = res;
     this->phi_cache = phi;
     this->psi_cache = psi;

     if (res == ALA || res == GLY) {
          this->res_cache = -1;
          this->chis_cache = chis;
          return chis;
     }
     // setting up the data and mismask
     // array properly
     mocapy::Sequence data = this->get_data(res);
     mocapy::MDArray<mocapy::eMISMASK> mism_sample = this->get_mismask(res);
     //
     // adding the backbone values .. if set
     double tmp = phi;
     if (is_initialized(tmp)) {
          // not at the terminus
          data.set(0, 2, tmp);
          mism_sample.set(0, 2, mocapy::MOCAPY_OBSERVED);
     }
     tmp = psi;
     if (is_initialized(tmp)) {
          // not at the terminus
          data.set(1, 2, tmp);
          mism_sample.set(1, 2, mocapy::MOCAPY_OBSERVED);
     }

     // Setup the proper sampler
     this->set_sampler(data, mism_sample);

     mocapy::MDArray<double> sample;
     sample = this->sampler->sample_next();

     // store the chi
     for (int i = 2; i < (int) sample.size() / 3; i++) {
          chis.push_back(sample.get(i, 2));
     }

     // Correct for alternative definition of HIS chi2
     if (res==HIS) {
          chis[1] = fmod(chis[1] + 2*M_PI, 2*M_PI) - M_PI;
     }

     this->chis_cache = chis;

     return chis;
}

//! Calculate the likelihood - based on cached values
double BasiliskDBN::calc_ll(bool ignore_bb) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     if (this->sampler == NULL || this->res_cache < 0) {
          return 0.;
     }
     // in case we have had a alanine or glycine last
     if (this->res_cache == ALA || this->res_cache == GLY || this->chis_cache.size() < 1) {
          return 0.;
     }

     return this->calc_ll(this->res_cache, this->phi_cache, this->psi_cache, this->chis_cache, false);
}


//! Calculate the likelihood
double BasiliskDBN::calc_ll(int aa, double phi, double psi, std::vector<double> chis, bool respect_symmetry) {
     definitions::ResidueEnum res = (definitions::ResidueEnum) aa;
     return this->calc_ll(res, phi, psi, chis, respect_symmetry);
}


//! Calculate the likelihood
double BasiliskDBN::calc_ll(definitions::ResidueEnum res, double phi, double psi, std::vector<double> chis, bool respect_symmetry) {

     // For debugging purposes. Can be removed.
     assert((!is_initialized(phi) || phi > -4.) && (!is_initialized(psi) || psi > -4.));

     // Import protein definitions (such as residue names)
     using namespace definitions;

     if (res == ALA || res == GLY)
          return 0.;
     // setting up the data and mismask
     // array properly
     mocapy::Sequence data = this->get_data((int) res);
     mocapy::MDArray<mocapy::eMISMASK> mism = this->get_mismask_ll((int) res);

     bool has_bb = false;

     // Correct for alternative definition of HIS chi2
     if (res==HIS) {
          chis[1] = fmod(chis[1] + 2*M_PI, 2*M_PI) - M_PI;
     }

     //
     // adding the backbone values .. if set
     if (is_initialized(phi)) {
          // not at the terminus
          data.set(0, 2, phi);
          mism.set(0, 2, mocapy::MOCAPY_OBSERVED);
          has_bb = true;
     } else {
          mism.set(0, 2, mocapy::MOCAPY_HIDDEN);
     }
     if (is_initialized(psi)) {
          // not at the terminus
          data.set(1, 2, psi);
          mism.set(1, 2, mocapy::MOCAPY_OBSERVED);
          has_bb = true;
     } else {
          mism.set(1, 2, mocapy::MOCAPY_HIDDEN);
     }
     //
     //
     double tmp = 0;
     for (int i = 2; i < (int) (data.size() / 3); i++) {
          tmp = chis[i - 2];

          // For debugging purposes. Can be removed.
          assert((!is_initialized(tmp) || tmp > -4.));

          if (is_initialized(tmp)) {
               data.set(i, 2, tmp);
               mism.set(i, 2, mocapy::MOCAPY_OBSERVED);
          } else {
               mism.set(i, 2, mocapy::MOCAPY_HIDDEN);
          }

     }

     double t1 = this->ll_engine.calc_ll(data, mism);
     /* Do some corrections for symetric sidechains */
     if (respect_symmetry && (res == TYR || res == ASP || res == PHE)) {
          double tmp = 0;
          for (int i = 2; i < (int) (data.size() / 3); i++) {
               tmp = chis[i - 2];

               // For debugging purposes. Can be removed.
               assert((!is_initialized(tmp) || tmp > -4.));

               if ((i - 2) == 1 && is_initialized(tmp)) {
                    // the angle in question (x2)
                    data.set(i, 2, tmp + M_PI);
                    mism.set(i, 2, mocapy::MOCAPY_OBSERVED);
               } else if (is_initialized(tmp)) {
                    // the other angles
                    data.set(i, 2, tmp);
                    mism.set(i, 2, mocapy::MOCAPY_OBSERVED);
               } else {
                    // we should never see that
                    mism.set(i, 2, mocapy::MOCAPY_HIDDEN);
               }
          }
          tmp = this->ll_engine.calc_ll(data, mism);
          //
          t1 = std::log(std::exp(t1) + std::exp(tmp));
     } else if (respect_symmetry && res == GLU) {
          double tmp = 0;
          for (int i = 2; i < (int) (data.size() / 3); i++) {
               tmp = chis[i - 2];

               // For debugging purposes. Can be removed.
               assert((!is_initialized(tmp) || tmp > -4.));

               if ((i - 2) == 2 && is_initialized(tmp)) {
                    data.set(i, 2, tmp + M_PI);
                    mism.set(i, 2, mocapy::MOCAPY_OBSERVED);
               } else if (is_initialized(tmp)) {
                    data.set(i, 2, tmp);
                    mism.set(i, 2, mocapy::MOCAPY_OBSERVED);
               } else {
                    mism.set(i, 2, mocapy::MOCAPY_HIDDEN);
               }
          }
          tmp = this->ll_engine.calc_ll(data, mism);
          //
          t1 = std::log(std::exp(t1) + std::exp(tmp));
     }

     // if we have a backbone dependent sample, we
     // have to subtract the backbone contribution
     if (has_bb) {
          for (int i = 2; i < (int) (data.size() / 3); i++) {
               mism.set(i, 2, mocapy::MOCAPY_HIDDEN);
          }
          double bb = this->ll_engine.calc_ll(data, mism);
          t1 = t1 - bb;
     }

     return t1;
}

//! Select the relevant sampler based on the sidechain length
void BasiliskDBN::set_sampler(mocapy::Sequence &data, mocapy::MDArray<mocapy::eMISMASK> &mismask) {

     int len = (data.size() / 3) - 2;

     if (len == 1) {
          if (sampler_len_1 == NULL) {
               this->sampler_len_1 = new mocapy::SampleInfEngineHMM(&dbn, data, mismask, 1);
          } else {
               this->sampler_len_1->set_seq_mismask(data, mismask);
          }
          this->sampler = this->sampler_len_1;
     } else if (len == 2) {
          if (sampler_len_2 == NULL) {
               this->sampler_len_2 = new mocapy::SampleInfEngineHMM(&dbn, data, mismask, 1);
          } else {
               this->sampler_len_2->set_seq_mismask(data, mismask);
          }
          this->sampler = this->sampler_len_2;
     } else if (len == 3) {
          if (sampler_len_3 == NULL) {
               this->sampler_len_3 = new mocapy::SampleInfEngineHMM(&dbn, data, mismask, 1);
          } else {
               this->sampler_len_3->set_seq_mismask(data, mismask);
          }
          this->sampler = this->sampler_len_3;
     } else if (len == 4) {
          if (sampler_len_4 == NULL) {
               this->sampler_len_4 = new mocapy::SampleInfEngineHMM(&dbn, data, mismask, 1);
          } else {
               this->sampler_len_4->set_seq_mismask(data, mismask);
          }
          this->sampler = this->sampler_len_4;
     }
}

}
