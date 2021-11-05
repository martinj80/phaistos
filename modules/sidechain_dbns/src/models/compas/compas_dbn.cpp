// compas_dbn.cpp --- Encapsulates all dbn handling for Compas
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


#include "models/compas/compas_dbn.h"
#include "mocapy.h"
#include "utils/random.h"

#include <stdio.h>
#include "utils/random.h"

namespace phaistos {

//! Name
const std::string CompasDBN::name = "compas";

//! Map from amino acid indices to the reduced amino acid indices used in the model
const int CompasDBN::aa_to_model_aa_array[] = { 99, 0, 1, 2, 3, 99, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 };

//! Map from the amino acid indices in the model to the normal amino acid indices
const int CompasDBN::model_aa_to_aa_array[] = { 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };

//! Side chain mean vectors
const double CompasDBN::aa_mean_vectors[20 * 3] = { 0., 0., 0.,          // ALA 0
                                                    2.80, 2.80, 2.80,    // CYS 1
                                                    2.92, 2.92, 2.92,    // ASP 2
                                                    3.125, 3.875, 3.125, // GLU 3
                                                    3.79, 3.79, 3.79,    // PHE 4
                                                    0., 0., 0.,          // GLY 5
                                                    3.57, 3.57, 3.57,    // HIS 6
                                                    2.25, 2.7, 2.7,      // ILE 7
                                                    3.8, 4.2, 4.6,       // LYS 8
                                                    3.05, 3.05, 3.05,    // LEU 9
                                                    3.185, 3.85, 3.185,  // MET 10
                                                    2.91, 2.91, 2.91,    // ASN 11
                                                    2.29, 2.29, 2.29,    // PRO 12
                                                    3.125, 3.875, 3.875, // GLN 13
                                                    4.31, 4.8, 5.31,     // ARG 14
                                                    2.43, 2.43, 2.43,    // SER 15
                                                    2.17, 2.17, 2.17,    // THR 16
                                                    2.19, 2.19, 2.19,    // VAL 17
                                                    4.2, 4.5, 4.2,       // TRP 18
                                                    4.27, 4.27, 4.27     // TYR 19
                                               };


//! Unit omega/phi angles to x,y,z vector
std::vector<double> polar_to_xyz(double omega, double phi) {
     std::vector<double> ar;
     ar.push_back(cos(omega)); // x
     ar.push_back(sin(omega) * cos(phi)); // y
     ar.push_back(sin(omega) * sin(phi)); // z
     return ar;
}

//! Unit vector to array of omega/phi
std::vector<double> xyz_to_polar(double x, double y, double z) {
     std::vector<double> ar;
     ar.push_back(acos(x));
     ar.push_back(atan2(z, y));
     return ar;
}

//! Convert distance to bin value
// This function could potentially be improved
// by saving a vector<vector<pair<double,int> > >,
// maintaining a vector of <distance_cutoff, bin>
// pairs for each amino acid.
int dist_to_bin(int aa, double dist) {
     int m3l = 0;

     //### Cystein , uni
     if (aa == 0) {
          m3l = 0;
     } else if (aa == 1) {
          //### Aspartate , uni
          m3l = 0;
     } else if (aa == 2) {
          //### Glutamate, two , 3.67
          if (dist < 3.67) {
               m3l = 0;
          } else {
               m3l = 1;
          }
     } else if (aa == 3) {
          //### Phenylalanin , uni
          m3l = 0;
     } else if (aa == 4) {
          //### Histidine, uni
          m3l = 0;
     } else if (aa == 5) {
          //### Isoleucine, two, 2.5
          if (dist < 2.5) {
               m3l = 0;
          } else {
               m3l = 1;
          }
     } else if (aa == 6) {
          //### Lysine , three, 4.04 4.375
          if (dist < 4.04) {
               m3l = 0;
          } else if (dist < 4.375) {
               m3l = 1;
          } else {
               m3l = 2;
          }
     } else if (aa == 7) {
          //### Leucine, uni
          m3l = 0;
     } else if (aa == 8) {
          //### Methionine, two , 3.625
          if (dist < 3.625) {
               m3l = 0;
          } else {
               m3l = 1;
          }
     } else if (aa == 9) {
          //### Asparagine, uni
          m3l = 0;
     } else if (aa == 10) {
          //### Proline, uni
          m3l = 0;
     } else if (aa == 11) {
          //### Glutamine, two 3.625
          if (dist < 3.625) {
               m3l = 0;
          } else {
               m3l = 1;
          }
     } else if (aa == 12) {
          //### Arginine, three, 4.5, 5.125
          if (dist < 4.5) {
               m3l = 0;
          } else if (dist < 5.125) {
               m3l = 0;
          } else {
               m3l = 1;
          }
     } else if (aa == 13) {
          //### Serine, uni
          m3l = 0;
     } else if (aa == 14) {
          //### Threonine, uni
          m3l = 0;
     } else if (aa == 15) {
          //### Valine, uni
          m3l = 0;
     } else if (aa == 16) {
          //### Tryptophan , two , 4.375
          if (dist < 4.375) {
               m3l = 0;
          } else {
               m3l = 1;
          }
     } else if (aa == 17) {
          //### Tyrosine, uni
          m3l = 0;
     }

     return m3l;
}

//! Translate amino acid label to internal model index
int CompasDBN::aa_to_model_aa(int aa) {
     return aa_to_model_aa_array[aa];
}

//! Translate internal model index to amino acid index
int CompasDBN::model_aa_to_aa(int model_aa) {
     return model_aa_to_aa_array[model_aa];
}

//! Get the distance in Angstrom
double CompasDBN::get_distance(int aa, int dist_bin) {
     return aa_mean_vectors[(aa * 3) + dist_bin];
}

//! Initializer
void CompasDBN::init() {
     this->mm = NULL;
     this->mm_ll = NULL;
     this->mm_no_bb = NULL;
}

//! Initialize Mocapy sampler object
void CompasDBN::init_sampler() {
     // inf engine
     this->ll_engine = *(new mocapy::LikelihoodInfEngineMM(&this->dbn, 1));
}

//! Constructor
CompasDBN::CompasDBN(std::string dbn_file, uint seed) {
     assert(file_exists(dbn_file));
     this->dbn.load(dbn_file.c_str());
     this->dbn.seed(seed);

     this->thread_id = 0;
     this->res_cache = 0;
     this->sampler = NULL;
     this->init();
     this->init_sampler();

     this->dbn_filename = dbn_file;
}

//! Destructor
CompasDBN::~CompasDBN() {
     if (this->sampler != NULL)
          delete this->sampler;
     if (this->mm_ll != NULL)
          delete this->mm_ll;
     if (this->mm_no_bb != NULL)
          delete this->mm_no_bb;
     if (this->mm != NULL)
          delete this->mm;

     this->sampler = NULL;
     this->mm = NULL;
     this->mm_ll = NULL;
     this->mm_no_bb = NULL;
}

//! Copy constructor
CompasDBN::CompasDBN(const CompasDBN &other) {

     this->res_cache = other.res_cache;

     this->sampler = other.sampler;
     this->ll_engine = other.ll_engine;

     this->mm = other.mm;
     this->mm_ll = other.mm_ll;
     this->mm_no_bb = other.mm_no_bb;

     // copy and reseed the DBN
     this->dbn = other.dbn;
     this->dbn_filename = other.dbn_filename;

     this->copies = other.copies;
     this->thread_id = other.thread_id;
     this->random_number_engine = other.random_number_engine;
}

//! Copy constructor
CompasDBN::CompasDBN(const CompasDBN &other, RandomNumberEngine *random_number_engine, int thread_id) {

     int seed = (*random_number_engine)();

     this->init();

     // copy and reseed the DBN
     this->dbn = *(new mocapy::DBN);
     this->dbn.load(other.dbn_filename.c_str());
     this->dbn.seed(seed);
     this->dbn_filename = other.dbn_filename;

     this->res_cache = other.res_cache;

     this->sampler = NULL;

     this->copies = other.copies;
     this->thread_id = thread_id;
     this->random_number_engine = random_number_engine;

     this->res_cache = 0;

     this->init_sampler();

}

//! Create internal copies of DBN (for multithreading)
void CompasDBN::duplicate(int ncopies, const std::vector<RandomNumberEngine *> &random_number_engines) {

     // Resize threads vector
     assert(copies.size()==0);
     copies.resize(ncopies);

     if (thread_id == 0) {
          // Use the passed dbn as first thread
          copies[thread_id] = this;

          // Make nthreads-1 copies
          for (int i = 1; i < ncopies; i++) {

               if (random_number_engines.size() > 0) {
                    copies[i] = new CompasDBN(*this, random_number_engines[i], i);
               } else {
                    copies[i] = new CompasDBN(*this, this->random_number_engine, i);
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
CompasDBN &CompasDBN::get_copy(unsigned int index) {
     if (copies.size() > index) {
          return *copies[index];
     } else if (index == 0) {
          return *this;
     } else {
          // die die die my darling .. (Misfits, 1984)
          assert(false);
     }
}


//! Get data - initialize Mocapy sequence object
mocapy::Sequence CompasDBN::get_data(ResidueFB *res, bool ignore_bb) {
     mocapy::Sequence data;

     std::vector<double> dof = (*res).get_sidechain_dof_values();
     std::vector<double> xyz = polar_to_xyz(dof[1], dof[2]);

     int aa = this->aa_to_model_aa((int) (*res).residue_type);

     // create the array
     data.set_shape(1, 9);
     // Phi
     data.set(0, 0, 0); // dummy
     data.set(0, 1, 0); // hidden
     if (ignore_bb) {
          data.set(0, 2, 0); // torus phi
          data.set(0, 3, 0); // torus psi
     } else {
          data.set(0, 2, (*res).get_phi()); // torus phi
          data.set(0, 3, (*res).get_psi()); // torus psi
     }
     data.set(0, 4, aa); // aa
     data.set(0, 5, xyz[0]); // kent x
     data.set(0, 6, xyz[1]); // kent y
     data.set(0, 7, xyz[2]); // kent z
     data.set(0, 8, dist_to_bin(aa, dof[0])); // dist bin

     return data;
}

//! Get data - initialize Mocapy sequence object
mocapy::Sequence CompasDBN::get_data(definitions::ResidueEnum res) {
     mocapy::Sequence data;

     int aa = this->aa_to_model_aa((int) res);

     // create the array
     data.set_shape(1, 9);
     // Phi
     data.set(0, 0, 0); // dummy
     data.set(0, 1, 0); // hidden
     // bb
     data.set(0, 2, 0); // torus phi
     data.set(0, 3, 0); // torus psi
     // aa
     data.set(0, 4, aa); // aa
     // unit vector
     data.set(0, 5, 1); // kent x
     data.set(0, 6, 0); // kent y
     data.set(0, 7, 0); // kent z
     //dist
     data.set(0, 8, 0); // dist bin

     return data;
}

//! Get mismask for sampling
mocapy::MDArray<mocapy::eMISMASK> CompasDBN::get_mismask() {
     if (this->mm == NULL) {
          mocapy::MDArray<mocapy::eMISMASK> m_seq;
          m_seq.set_shape(1, 6);
          m_seq.set(0, 0, mocapy::MOCAPY_HIDDEN); // dummy hidden
          m_seq.set(0, 1, mocapy::MOCAPY_HIDDEN); // hidden
          m_seq.set(0, 2, mocapy::MOCAPY_OBSERVED); // torus
          m_seq.set(0, 3, mocapy::MOCAPY_OBSERVED); // aa
          m_seq.set(0, 4, mocapy::MOCAPY_HIDDEN); // kent
          m_seq.set(0, 5, mocapy::MOCAPY_HIDDEN); // dist
          this->mm = new mocapy::MDArray<mocapy::eMISMASK>(m_seq);
     }
     return *(this->mm);
}


//! Get mismask for sampling without backbone dependencies
mocapy::MDArray<mocapy::eMISMASK> CompasDBN::get_mismask_no_bb() {
     if (this->mm_no_bb == NULL) {
          mocapy::MDArray<mocapy::eMISMASK> m_seq;
          m_seq.set_shape(1, 6);
          m_seq.set(0, 0, mocapy::MOCAPY_HIDDEN); // dummy hidden
          m_seq.set(0, 1, mocapy::MOCAPY_HIDDEN); // hidden
          m_seq.set(0, 2, mocapy::MOCAPY_HIDDEN); // torus
          m_seq.set(0, 3, mocapy::MOCAPY_OBSERVED); // aa
          m_seq.set(0, 4, mocapy::MOCAPY_HIDDEN); // kent
          m_seq.set(0, 5, mocapy::MOCAPY_HIDDEN); // dist
          this->mm_no_bb = new mocapy::MDArray<mocapy::eMISMASK>(m_seq);
     }
     return *(this->mm_no_bb);
}

//! Get mismask for likelihood calculations
mocapy::MDArray<mocapy::eMISMASK> CompasDBN::get_mismask_ll() {
     if (this->mm_ll == NULL) {
          mocapy::MDArray<mocapy::eMISMASK> m_seq;
          m_seq.set_shape(1, 6);
          m_seq.set(0, 0, mocapy::MOCAPY_HIDDEN);
          m_seq.set(0, 1, mocapy::MOCAPY_HIDDEN);
          m_seq.set(0, 2, mocapy::MOCAPY_HIDDEN);
          m_seq.set(0, 3, mocapy::MOCAPY_OBSERVED);
          m_seq.set(0, 4, mocapy::MOCAPY_OBSERVED);
          m_seq.set(0, 5, mocapy::MOCAPY_OBSERVED);
          this->mm_ll = new mocapy::MDArray<mocapy::eMISMASK>(m_seq);
     }
     return *(this->mm_ll);
}


//! Generate a new sequence of angles for a specific residue in the chain
std::vector<double> CompasDBN::get_angle_sequence_for_residue_with_bb(definitions::ResidueEnum res, 
                                                                      double phi, double psi) {

     // For debugging purposes. Can be removed.
     assert((!is_initialized(phi) || phi > -4.) && (!is_initialized(psi) || psi > -4.));

     // Import protein definitions (such as residue names)
     using namespace definitions;

     std::vector<double> angles;

     // Skip residues with no sidechain
     if (res == ALA || res == GLY) {
          angles.push_back(UNINITIALIZED);
          angles.push_back(UNINITIALIZED);
          angles.push_back(UNINITIALIZED);
          this->res_cache = -1;
          return angles;
     }

     // set up the data arrays
     mocapy::Sequence data;
     data = this->get_data(res);

     mocapy::MDArray<mocapy::eMISMASK> mism;

     if (is_initialized(phi) && is_initialized(psi)) {
          // only use the BB mismask if it really makes sense
          data.set(0, 2, phi); // torus phi
          data.set(0, 3, psi); // torus psi
          mism = this->get_mismask();
     } else {
          mism = this->get_mismask_no_bb();
     }
     //
     // setting up the sampler
     if (this->sampler != NULL) {
          this->sampler->set_seq_mismask(data, mism);
     } else {
          this->sampler = new mocapy::SampleInfEngineMM(&dbn, data, mism, 1);
     }

     // get a sample
     mocapy::MDArray<double> sample = this->sampler->sample_next();

     std::vector<double> tmp = xyz_to_polar(sample.get(0, 5), sample.get(0, 6), sample.get(0, 7));
     angles.push_back(this->get_distance((int) res, (int) sample.get(0, 8)));
     angles.push_back(tmp[0]); // omega
     angles.push_back(tmp[1]); // phi

     //! Save residue information
     this->res_cache = res;

     return angles;
}

//! Calculate the likelihood - for use immediately after sampling
double CompasDBN::calc_ll(bool ignore_bb) {
     if (this->sampler == NULL || this->res_cache < 0) {
          return 0.;
     }

     mocapy::MDArray<mocapy::eMISMASK> mism = get_mismask_ll();
     //
     // check whether or not to use the
     // backbone information.
     if (not ignore_bb) {
          mism.set(0, 2, mocapy::MOCAPY_OBSERVED);
     }
     return this->sampler->calc_ll(mism);
}

//! Calculates the likelihood (wrapper)
double CompasDBN::calc_ll(int aa, double phi, double psi, std::vector<double> dof) {
     definitions::ResidueEnum res = (definitions::ResidueEnum) aa;
     return this->calc_ll(res, phi, psi, dof);
}

//! Calculates the likelihood
double CompasDBN::calc_ll(definitions::ResidueEnum res, double phi, double psi, std::vector<double> dof) {

     // For debugging purposes. Can be removed.
     assert((!is_initialized(phi) || phi > -4.) && (!is_initialized(psi) || psi > -4.));

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Skip residues with no sidechain
     if (res == ALA || res == GLY) {
          return 0.;
     }

     assert(dof.size() == 3);

     // set up the data arrays
     mocapy::Sequence data;
     data = this->get_data(res);
     //
     // get the mismask
     mocapy::MDArray<mocapy::eMISMASK> mism;
     mism = this->get_mismask_ll();
     //
     // check whether or not to use the
     // backbone information.
     if (is_initialized(phi) && is_initialized(psi)) {
          // only use the BB mismask if it really makes sense
          data.set(0, 2, phi); // torus phi
          data.set(0, 3, psi); // torus psi
          mism.set(0, 2, mocapy::MOCAPY_OBSERVED);
     }
     //
     // turn the angles into coordinates.
     std::vector<double> xyz = polar_to_xyz(dof[1], dof[2]);
     data.set(0, 5, xyz[0]); // kent x
     data.set(0, 6, xyz[1]); // kent y
     data.set(0, 7, xyz[2]); // kent z
     // and finally the distance
     int aa = this->aa_to_model_aa((int) res);
     data.set(0, 8, dist_to_bin(aa, dof[0])); // dist bin

     return this->ll_engine.calc_ll(data, mism);;
}

}
