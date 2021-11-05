// opls_non_bonded.h --- OPLS non-bonded interaction (gbsa + vdw + charge) energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson, Sandro Bottaro, Wouter Boomsma
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

#ifndef OPLS_NONBONDED_H
#define OPLS_NONBONDED_H

#include <vector>
#include <iostream>
#include <cmath>

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "protein/iterators/pair_iterator_chaintree.h"
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "energy/parameter_data_structures.h"
#include "energy/opls_parameters.h"

namespace phaistos {

//! Parameter container class for pairwise non-bonded parameters used in NonbondedParameters class.
template <typename FLOAT_TYPE>
class PairParameter {

public:

     //! Geometric mean of the two atoms radii squared,
     //! note that sigma values are given in the OPLS parameter file
     FLOAT_TYPE vdw_radius_sq;

     //! Geometric mean of the two atoms epsilon value
     FLOAT_TYPE vdw_epsilon;

     //! GBSA radii for the two atoms added and squared
     FLOAT_TYPE gbsa_radius_sum_sq;

     //! Extra field to cache interatomic distance squared
     FLOAT_TYPE dist_sq;
     
     //! Boolean interpretation
     operator bool () {return (vdw_radius_sq > -0.0001); }; //radius=0.0 is used as off signal
     
     //! Default constructor
     PairParameter()
          : vdw_radius_sq(-1.0), 
            gbsa_radius_sum_sq(0.0), 
            dist_sq(0.0) {};
     
     //! Main constructor
     PairParameter(FLOAT_TYPE vdw_rad,FLOAT_TYPE vdw_eps,FLOAT_TYPE gbsa_rad)
          : vdw_radius_sq(vdw_rad), 
            vdw_epsilon(vdw_eps), 
            gbsa_radius_sum_sq(gbsa_rad), 
            dist_sq(0.0) {};
     
};

//! Container class for single atoms parameters used in NonbondedParameters class
template <typename FLOAT_TYPE>
class Parameter {

public:

     //! van der Waals radius
     FLOAT_TYPE radius;

     //! Partial charge
     FLOAT_TYPE charge;
     
     //! Boolean interpretation
     operator bool () {return (radius > -0.0001); }; //radius=0.0 is used as off signal
     
     //! Default constructor
     Parameter():radius(-1.0) {};
     
     //! Main constructor
     Parameter(FLOAT_TYPE rad,FLOAT_TYPE chg): radius(rad), charge(chg) {};
     
};

//! Parameter container class for OPLS nonbonded energy term
template <typename FLOAT_TYPE>
class NonbondedParameters: public TinkerParameterBase {

private:

     //! Vector to store parameters in
     std::vector<Parameter<FLOAT_TYPE> > param_cache;

     //! Vector to store pairwise parameters in
     std::vector<std::vector<PairParameter<FLOAT_TYPE> > > pair_cache;

     //! Binary map to store bond lentghs in
     ParameterData2D<FLOAT_TYPE> bond_map;

public:

     //! Number of single atom parameters
     int max_param_id;

     //! Number of pairwise parameters
     int max_pair_id;
     
     //! van der Waals scaling factor for third neighbour given in parameter file header
     FLOAT_TYPE vdw14scale;

     //! Charge scaling factor for third neighbour given in parameter file header
     FLOAT_TYPE chg14scale;

     //! Unit of the energy calculated from the dielectric constant given in parameter file header
     FLOAT_TYPE chargeUnit;

     //! Constructor
     NonbondedParameters() {

          // CHARGE PARAMETERS
          // get third neighbour scale factor from parameter file header
          chg14scale = 1.0;
          std::string chg14scaleStr = read_header(&parameters_opls,"chg-14-scale");
          chg14scale /= atof( chg14scaleStr.c_str() );
     
          // get charge unit from parameter file header
          chargeUnit = 332.05382; //1 e^2/Å
          std::string dielectricStr = read_header(&parameters_opls,"dielectric");
          chargeUnit /= atof( dielectricStr.c_str() );
     
          // get parameters from 'charge' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls,&raw_param,"charge",1);
     
          // GBSA PARAMETERS
          // get vdw radii parameters (with pair_param index!)
          std::vector<FLOAT_TYPE> gbsa_radii = get_gbsa_radii();
     
          // do not use index zero
          param_cache.push_back(Parameter<FLOAT_TYPE>());

          // store single-index parameters
          unsigned int id(0);
          unsigned int raw_param_size = raw_param.size();
          for (unsigned int i=0; i<raw_param_size; i++) {
               id = (raw_param[i].first)[0];
               FLOAT_TYPE radius = gbsa_radii[atom_id2param_id(id)];
               FLOAT_TYPE charge = (FLOAT_TYPE) raw_param[i].second[0];
               Parameter<FLOAT_TYPE> param(radius,charge);
               if (id != param_cache.size()) {
                    std::cerr<<"\nOPLS TERM NONBONDED - read error around charge parameter id"
                             <<id<<"\n\n";
               }
               param_cache.push_back(param);
          }
          max_param_id = id;
     
          // VDW PARAMETERS
          //Sixth Root of Two
          const FLOAT_TYPE srt = 1.122462048309372981;
     
          // get third neighbour scale factor from parameter file header
          vdw14scale = 1.0;
          std::string vdw14scaleStr = read_header(&parameters_opls,"vdw-14-scale");
          vdw14scale /= atof( vdw14scaleStr.c_str() );
     
          // get parameters from 'vdw' field
          raw_param.clear();
          read_param(&parameters_opls,&raw_param,"vdw",1);
     
          // convert vdw sigma values to radius
          raw_param_size = raw_param.size();
          for (unsigned int i=0; i<raw_param_size; i++) {
               /* radiussize              DIAMETER */
               raw_param[i].second[0] *= 0.5;
               /* radiustype              SIGMA    */
               raw_param[i].second[0] *= srt;
          }
     
          // store pair-index parameters
          unsigned int id_i(0);
          std::vector<PairParameter<FLOAT_TYPE> > buf;
          for (unsigned int i=0; i<raw_param_size; i++) {
               id_i = (raw_param[i].first)[0];
               std::vector<double> param_i = raw_param[i].second;
               for (unsigned int j=i; j<raw_param_size; j++) {
                    unsigned int id_j = (raw_param[j].first)[0];
                    std::vector<double> param_j = raw_param[j].second;
                    assert(id_j == buf.size()+1+i);
                    // if (id_j != buf.size()+1+i) {
                    //      std::cerr<<"\nOPLS TERM NONBONDED - read error around vdw parameter id "
                    //               <<id_j<<"\n\n";
                    // }
                    /* radiusrule  GEOMETRIC */
                    FLOAT_TYPE radius = (FLOAT_TYPE) 2*sqrt(param_i[0]*param_j[0]);
                    FLOAT_TYPE radius_sq = radius*radius;
                    /* epsilonrule GEOMETRIC */
                    FLOAT_TYPE epsilon = (FLOAT_TYPE) sqrt( param_i[1] * param_j[1]);
               
                    //GBSA params
                    FLOAT_TYPE sum = gbsa_radii[id_i] + gbsa_radii[id_j];
                    FLOAT_TYPE sum_sq = sum*sum;
               
                    buf.push_back( PairParameter<FLOAT_TYPE>(radius_sq,epsilon,sum_sq) );
               }
               assert(id_i == pair_cache.size()+1);
               // if (id_i != pair_cache.size()+1) {
               //      std::cerr<<"\nOPLS TERM NONBONDED - read error around vdw parameter id "
               //               <<id_i<<"\n\n";
               // }
               pair_cache.push_back(buf);
               buf.clear();
          }
          max_pair_id = id_i;
     
          // BOND PARAMETERS (eq length only)
          raw_param.clear();
          read_param(&parameters_opls,&raw_param,"bond",2);

          // put bond length params in map
          raw_param_size = raw_param.size();
          for (unsigned int i=0; i<raw_param_size; i++) {
               std::vector<int> id = raw_param[i].first;
               FLOAT_TYPE bondLen = (FLOAT_TYPE) raw_param[i].second[1];
               if (bond_map(id[0],id[1]))
                    std::cout<<"\nOPLS TERM NONBONDED - multiple definitions of bond parameter "
                             <<id[0]<<","<<id[1]<<"\n\n";
               else
                    bond_map(id[0],id[1]) = bondLen;
          }
     }

     //! Destructor
     ~NonbondedParameters() {};

     //! Get single-index parameter pointer
     const Parameter<FLOAT_TYPE> *get_param(Atom *atom) {
          int id = get_atom_id(atom);
          assert(id <= max_param_id);
          return ( &param_cache[id] );
     }

     //! Get pair-index parameter pointer
     const PairParameter<FLOAT_TYPE> *get_pair_param(Atom *atom1, Atom *atom2) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          assert(id1 <= max_pair_id);
          assert(id2 <= max_pair_id);
          // pair_cache is symmetric
          if (id1 > id2) {
               int buf=id1; id1=id2; id2=buf;
          }
          // pair_cache is zero-based
          id1--;
          id2--;
          // convert upper right to lower left triangle for easyer iteration
          id2 -= id1;
          return &pair_cache[id1][id2];
     }

     //! Get equilibrium bond length parameter as given in the parameter file
     FLOAT_TYPE get_bond(Atom *atom1, Atom *atom2) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          if (id1 > id2) {
               int buf=id1; id1=id2; id2=buf;
          }
          return bond_map(id1,id2);
     }

     //! Radii parameters from tinker GBSA implementation with tinker param id index
     std::vector<FLOAT_TYPE> get_gbsa_radii() {

          std::vector<FLOAT_TYPE> gbsa_radii;
          // Hydrogen
          gbsa_radii.push_back(0.0);     // index 0 - undef
          gbsa_radii.push_back(1.25);    // index 1
          gbsa_radii.push_back(1.15);
          gbsa_radii.push_back(1.05);
          // Carbon
          gbsa_radii.push_back(1.90);    // index 4
          gbsa_radii.push_back(1.875);
          gbsa_radii.push_back(1.825);
          // Nitrogen
          gbsa_radii.push_back(1.7063);  // index 7
          gbsa_radii.push_back(1.625);
          gbsa_radii.push_back(1.60);
          // Oxygen
          gbsa_radii.push_back(1.535);   // index 10
          gbsa_radii.push_back(1.48);
          // Sulfur
          gbsa_radii.push_back(1.775);   // index 12
     
          // Map to tinker param id
          std::vector<unsigned int> gbsa_radii_map;
          gbsa_radii_map.push_back(0);   // index zero not used in tinker

          gbsa_radii_map.push_back(4);   // index 1 
          gbsa_radii_map.push_back(1);   // index 2 
          gbsa_radii_map.push_back(0);   // index 3 
          gbsa_radii_map.push_back(0);   // index 4 
          gbsa_radii_map.push_back(5);   // index 5 
          gbsa_radii_map.push_back(1);   // index 6 
          gbsa_radii_map.push_back(10);  // index 7 
          gbsa_radii_map.push_back(3);   // index 8 
          gbsa_radii_map.push_back(0);   // index 9
          gbsa_radii_map.push_back(10);  // index 10
          gbsa_radii_map.resize(14+1,0); // index 11-14
          gbsa_radii_map.push_back(12);  // index 15
          gbsa_radii_map.push_back(0);  // index 16
          gbsa_radii_map.push_back(12);  // index 17
          gbsa_radii_map.push_back(1);   // index 18
          gbsa_radii_map.push_back(0);   // index 19
          gbsa_radii_map.push_back(0);   // index 20
          gbsa_radii_map.push_back(5);   // index 21
          gbsa_radii_map.push_back(11);  // index 22
          gbsa_radii_map.push_back(7);   // index 23
          gbsa_radii_map.push_back(2);   // index 24
          gbsa_radii_map.push_back(0);   // index 25
          gbsa_radii_map.push_back(0);   // index 26
          gbsa_radii_map.push_back(5);   // index 27
          gbsa_radii_map.push_back(11);  // index 28
          gbsa_radii_map.push_back(1);   // index 29
          gbsa_radii_map.push_back(8);   // index 30
          gbsa_radii_map.push_back(2);   // index 31 
          gbsa_radii_map.push_back(7);   // index 32
          gbsa_radii_map.resize(36+1,5); // index 33-36
          gbsa_radii_map.push_back(7);   // index 37
          gbsa_radii_map.resize(40+1,5); // index 38-40
          gbsa_radii_map.push_back(7);   // index 41
          gbsa_radii_map.push_back(4);   // index 42
          gbsa_radii_map.resize(64+1,0); // index 43-64
          gbsa_radii_map.resize(77+1,4); // index 65-77

          std::vector<FLOAT_TYPE> ret(78);
          for (int i=0; i<78; i++)
               ret[i] = gbsa_radii[gbsa_radii_map[i]];
          return ret;
     };
};


//! OPLS non-bonded interaction energy term - base class
template <typename DERIVED_CLASS, typename FLOAT_TYPE>
class TermOplsNonBondedBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS,ChainFB> EnergyTermCommon;     

protected:

     //! Parameters
     NonbondedParameters<FLOAT_TYPE> parameters;

     //! Number of interactions in the last evaluation
     int counter;

     //! Atom volumes. [res index][atoms index in res]
     std::vector<std::vector<FLOAT_TYPE> > volume;

     //! One divided by constant Born radius contribution from the first and second neighbour.
     std::vector<std::vector<FLOAT_TYPE> > inv_born_radius_neighbor;

     //! One divided by born radii of all atoms. [res index][atoms index in res]
     std::vector<std::vector<FLOAT_TYPE> > inv_born_radius;     

     //! Fill vector with "exclusive" volumes for every atom in the chain
     int volumes() {
          int size = this->chain->size();
          for (int r=0; r<size; r++) {
               Residue *res = &(*this->chain)[r];

               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom1 = res->atoms[a];
                    FLOAT_TYPE overlap = 0.0;
                    FLOAT_TYPE rvdw = parameters.get_param(atom1)->radius;
                    FLOAT_TYPE vol = ((4.0*M_PI*(rvdw*rvdw*rvdw))/3.0);
                    // Substract the spherical caps of the neighboring atoms
                    // for(deprecated::CovalentBondIterator it2(atom1, deprecated::CovalentBondIterator::DEPTH_1_ONLY); it2 != NULL; ++it2){
                    for(CovalentBondIterator<ChainFB> it2(atom1, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY); 
                        !it2.end(); ++it2){
                         FLOAT_TYPE rvdw_neigh = parameters.get_param(&*it2)->radius;
                         // FLOAT_TYPE mod_dist = 1.01*(parameters.get_bond(atom1,it2.getAtom()));
                         FLOAT_TYPE mod_dist = 1.01*(parameters.get_bond(atom1,&*it2));
                         FLOAT_TYPE h = rvdw*(1.0+(rvdw_neigh*rvdw_neigh - rvdw*rvdw - mod_dist*mod_dist)/(2.0*rvdw*mod_dist));
                         overlap += (M_PI*h*h*(3.0*rvdw-h))/3.0;
                    }

                    // Assign the volume
                    volume[atom1->residue->index][atom1->index] = vol-overlap;	
               }          
          }
          return 0;
     }

     //! Calculate the contributions to the born radii of the first and second neighbors of each atom.
     //! NOTE, the second neighbour term is not constant when using for instance the CRISP move.
     //! \param chain Molecule chain object
     //! \param inv_born_radius_neighbor Inverse constant Born radius contribution from the first and second neighbour.
     int born_radii_neighbors(ChainFB *chain,
                              std::vector<std::vector<FLOAT_TYPE> > *inv_born_radius_neighbor) {
          const FLOAT_TYPE p1 = 0.073;
          const FLOAT_TYPE p2 = 0.921;
          const FLOAT_TYPE dielectric = 332.05;
          const FLOAT_TYPE d_offset = -0.09;

          for (AtomIterator<ChainFB,definitions::ALL> it1(*chain); !it1.end(); ++it1) {
               Atom *atom1 = &*it1;
               FLOAT_TYPE rvdw_atom1 = parameters.get_param(atom1)->radius;

               // Add the first constant contribution
               (*inv_born_radius_neighbor)[atom1->residue->index][atom1->index] += -(0.5*dielectric)/(rvdw_atom1+d_offset+p1);

               for (AtomIterator<ChainFB,definitions::ALL> it2(it1+1); !it2.end(); ++it2) {
               
                    // Atom *atom2 = it2.getAtom();
                    Atom *atom2 = &*it2;
                    int distance = chain_distance<ChainFB>(atom1,atom2);
               
                    // Add the first neighbors contributions
                    if(distance == 1) {
                         FLOAT_TYPE bond_length = parameters.get_bond(atom1,atom2);
                    
                         FLOAT_TYPE bond4_first = bond_length*bond_length*bond_length*bond_length;
                         (*inv_born_radius_neighbor)[atom1->residue->index][atom1->index] += (p2*volume[atom2->residue->index][atom2->index])/bond4_first;
                         (*inv_born_radius_neighbor)[atom2->residue->index][atom2->index] += (p2*volume[atom1->residue->index][atom1->index])/bond4_first;

                    }
               
                    // // Add the second neighbors contributions
                    // else if(chain_distance == 2) {
                    //      FLOAT_TYPE dx = (FLOAT_TYPE) atom1->position[0] - atom2->position[0];
                    //      FLOAT_TYPE dy = (FLOAT_TYPE) atom1->position[1] - atom2->position[1];
                    //      FLOAT_TYPE dz = (FLOAT_TYPE) atom1->position[2] - atom2->position[2];
                    //      FLOAT_TYPE dist2 = dx*dx+dy*dy+dz*dz;
                    //      FLOAT_TYPE dist4 = dist2*dist2;
                    
                    //      (*inv_born_radius_neighbor)[atom1->residue->index][atom1->index] += (p3*volume[atom2->residue->index][atom2->index])/dist4;
                    //      (*inv_born_radius_neighbor)[atom2->residue->index][atom2->index] += (p3*volume[atom1->residue->index][atom1->index])/dist4;
                    // }
               }
          }
          return 0;
     }


     // Hacky bit-manipulating code that is killed by -O2 and -O3
     // Requires 32 bit IEEE 754 floats (?!)
     // #pragma requires gcc 4.4 or newer
     // from http://en.wikipedia.org/wiki/Fast_inverse_square_root
     /* #pragma optimize("",off) */
     /* float TermOplsNonBondedBase<FLOAT_TYPE>::fastInvSqrt(float x) { */
     /*           float xhalf = 0.5f*x; */
     /*           int i = *(int*)&x; */
     /*           i = 0x5f3759df - (i>>1); */
     /*           x = *(float*)&i ; */
     /*           // can be added for more accurate estimate */
     /*           /\* x = x*(1.5f - xhalf*x*x); *\/ */
     /*           return (x*(1.5f - xhalf*x*x)); */
     /* } */
     /* #pragma optimize("",on) */
          
     /* float fastInvSqrt(float x); */
     /* inline double invSqrt(double x) {return 1.0/sqrt(x);}; */
     /* inline float invSqrt(float x) {return 1.0/sqrtf(x);}; */
     
     // //! Exponential function with double presision
     // inline double expo(double x) {
     //      return std::exp(x);
     // };

     // //! Exponential function with single presision
     // inline float expo(float x) {
     //      return std::expf(x);
     // };


     //! Calculate non-polar solvation (solvent cavity) energy of an atom
     //! \param atom The atom
     //! \param inv_born_radius One divided by the Born radius of the atom
     //! \return energy in kcal/mol /(4*pi*sigma) 
     inline FLOAT_TYPE calc_npol_solv(Atom *atom, double inv_born_radius) {

          const FLOAT_TYPE r_solv = 1.4;   // A

          /* FLOAT_TYPE rvdw = vdw_radius[atom->biotype]; */
          FLOAT_TYPE rvdw = parameters.get_param(atom)->radius;
          /* FLOAT_TYPE ratio = rvdw/(*inv_born_radius)[atom->residue->index][atom->index]; */
          const FLOAT_TYPE ratio = rvdw*inv_born_radius;
          const FLOAT_TYPE pow2 = ratio*ratio;
          //sixth power for testing
          const FLOAT_TYPE pow6 = pow2*pow2*pow2;
          rvdw += r_solv;
          return rvdw*rvdw*pow6;
     };


     //! Allocate space in vectors
     void allocate_vectors() {
          volume.resize(this->chain->size());
          inv_born_radius.resize(this->chain->size());
          inv_born_radius_neighbor.resize(this->chain->size());
          for (int i=0; i<this->chain->size(); ++i) {
               volume[i].resize((*this->chain)[i].size());
               inv_born_radius[i].resize((*this->chain)[i].size());
               inv_born_radius_neighbor[i].resize((*this->chain)[i].size());
          }
     }     

public:

     //! Use same settings as base class
     typedef typename EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;     

     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsNonBondedBase(ChainFB *chain, std::string name, 
                           const Settings &settings=Settings(),
                           RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine) {

          // Annotate the chain with biotype
          parameters.annotate_chain(chain);

          // Allocate space in vectors
          allocate_vectors();

          // Cache the volume for each atom in the chain
          volumes();
     
          // Cache the first neighbour contributions for the born radii
          born_radii_neighbors(chain, &inv_born_radius_neighbor);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsNonBondedBase(const TermOplsNonBondedBase &other, 
                           RandomNumberEngine *random_number_engine,
                           int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters) {
          
          // Annotate the new chain with biotype
          parameters.annotate_chain(chain);

          // Allocate space in vectors
          allocate_vectors();

          // Cache the volume for each atom in the chain
          volumes();
     
          // Cache the first neighbour contributions for the born radii
          born_radii_neighbors(chain, &inv_born_radius_neighbor);          
     }

};


//! OPLS nonbonded interaction energy term
template <typename FLOAT_TYPE>
class TermOplsNonBonded: public TermOplsNonBondedBase<TermOplsNonBonded<FLOAT_TYPE>, 
                                                      FLOAT_TYPE> {

     //! Cache of all atom-atom distances (DOESN'T THIS BECOME EXCESSIVELY LARGE?)
     std::vector<std::vector<std::vector<std::vector<FLOAT_TYPE> > > > dist_sq_cache;

     //! For convenience, define local base class
     typedef phaistos::TermOplsNonBondedBase<TermOplsNonBonded<FLOAT_TYPE>, 
                                             FLOAT_TYPE> TermOplsNonBondedBase;

     //! Calculate pairwise polar solvation energy
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param chg1 partial charge of atom 1 in electron volt
     //! \param chg2 partial charge of atom 2 in electron volt
     //! \return energy in kcal/mol
     inline FLOAT_TYPE calc_pol_solv(Atom *atom1, Atom *atom2,
                                     FLOAT_TYPE born_radius1, FLOAT_TYPE born_radius2,
                                     const FLOAT_TYPE chg1, const FLOAT_TYPE chg2) {

          /* const FLOAT_TYPE dielectric = 332.05; */
          /* const FLOAT_TYPE epsilon = 78.3; */
          /* const FLOAT_TYPE fac = (-dielectric*(1.0 - 1.0/epsilon)); */
          const FLOAT_TYPE fac=-327.809;

     
          /* const FLOAT_TYPE dx = (FLOAT_TYPE) atom1->position[0] - atom2->position[0]; */
          /* const FLOAT_TYPE dy = (FLOAT_TYPE) atom1->position[1] - atom2->position[1]; */
          /* const FLOAT_TYPE dz = (FLOAT_TYPE) atom1->position[2] - atom2->position[2]; */
          /* const FLOAT_TYPE dist2 = dx*dx+dy*dy+dz*dz; */

          const FLOAT_TYPE dist2 = dist_sq_cache[atom1->residue->index][atom1->index][atom2->residue->index][atom2->index];
     
          FLOAT_TYPE c_ij = fac*chg1*chg2;
          if (atom1 == atom2)
               c_ij *= 0.5;

          const FLOAT_TYPE r_born12 = born_radius1*born_radius2;

          /* float x = dist2+expo(-(dist2*0.25*r_born12))/r_born12; */
          /* const FLOAT_TYPE f_ij = fastInvSqrt(x); */

          /* FLOAT_TYPE x = dist2+expo(-(dist2*0.25*r_born12))/r_born12; */
          /* const FLOAT_TYPE f_ij = invSqrt(x); */

          const FLOAT_TYPE f_ij = 1.0/sqrt(dist2+std::exp(-(dist2*0.25*r_born12))/r_born12);

          return c_ij*f_ij;
     }


     //! Evaluate pair interaction between 2 atoms in proximity
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param chg_e Pointer to charge energy sum
     //! \param vdw_e Pointer to vdw energy sum
     //! \param inv_born_radius One divided by born radii of all atoms. [res index][atoms index in res]
     //! \param pair_param Parameters for atom pair
     //! \param param1 Parameters for atom1
     //! \param param2 Parameters for atom2
     inline void calc_pair_interaction_prox(Atom *atom1, Atom *atom2,
                                            FLOAT_TYPE *chg_e, FLOAT_TYPE *vdw_e,
                                            std::vector<std::vector<FLOAT_TYPE> > *inv_born_radius,
                                            const PairParameter<FLOAT_TYPE> *pair_param,
                                            const Parameter<FLOAT_TYPE> *param1,
                                            const Parameter<FLOAT_TYPE> *param2) {

          const double p3 = 6.211;

          int d = this->parameters.template get_covalent_dist<ChainFB>(atom1,atom2);
          if (d > 3) {
               calc_pair_interaction(
                    atom1,atom2,chg_e,vdw_e,inv_born_radius,pair_param,param1,param2);
          } else if (d < 3) {
               FLOAT_TYPE dx = (FLOAT_TYPE) atom1->position[0] - atom2->position[0];
               FLOAT_TYPE dy = (FLOAT_TYPE) atom1->position[1] - atom2->position[1];
               FLOAT_TYPE dz = (FLOAT_TYPE) atom1->position[2] - atom2->position[2];
               FLOAT_TYPE dist2 = dx*dx + dy*dy + dz*dz;
               FLOAT_TYPE dist4 = dist2*dist2;
               dist_sq_cache[atom1->residue->index][atom1->index][atom2->residue->index][atom2->index] = dist2;

               if (d == 2) {
                    (*inv_born_radius)[atom1->residue->index][atom1->index] +=
                         (p3*this->volume[atom2->residue->index][atom2->index])/dist4;
                    (*inv_born_radius)[atom2->residue->index][atom2->index] +=
                         (p3*this->volume[atom1->residue->index][atom1->index])/dist4;
               }
               
               return;
          } else {//d==3
               FLOAT_TYPE chg_ep=0.0,vdw_ep=0.0;
               calc_pair_interaction(
                    atom1,atom2,&chg_ep,&vdw_ep,inv_born_radius,pair_param,param1,param2);
               *chg_e += chg_ep*this->parameters.chg14scale;
               *vdw_e += vdw_ep*this->parameters.vdw14scale;
          }
     }

     //! Evaluate pair interaction
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param chg_e Pointer to charge energy sum
     //! \param vdw_e Pointer to vdw energy sum
     //! \param inv_born_radius One divided by born radii of all atoms. [res index][atoms index in res]
     //! \param pair_param Parameters for atom pair
     //! \param param1 Parameters for atom1
     //! \param param2 Parameters for atom2
     inline void calc_pair_interaction(Atom *atom1, Atom *atom2,
                                       FLOAT_TYPE *chg_e, FLOAT_TYPE *vdw_e,
                                       std::vector<std::vector<FLOAT_TYPE> > *inv_born_radius,
                                       const PairParameter<FLOAT_TYPE> *pair_param,
                                       const Parameter<FLOAT_TYPE> *param1,
                                       const Parameter<FLOAT_TYPE> *param2) {

          this->counter++;

          const FLOAT_TYPE p4=15.236;
          const FLOAT_TYPE p5pi=3.939441;  // = 1.254*3.1415
          const FLOAT_TYPE p5inv=0.797448; // = 1.0/1.254

          /* FLOAT_TYPE dx = (FLOAT_TYPE) atom1->position[0] - atom2->position[0]; */
          /* FLOAT_TYPE dy = (FLOAT_TYPE) atom1->position[1] - atom2->position[1]; */
          /* FLOAT_TYPE dz = (FLOAT_TYPE) atom1->position[2] - atom2->position[2]; */
          /* FLOAT_TYPE r_sq = dx*dx + dy*dy + dz*dz; */

          /* *chg_e += (param1->charge)*(param2->charge)*invSqrt(r_sq); */

          /* *chg_e += (param1->charge)*(param2->charge)*fastInvSqrt(r_sq); */

          const FLOAT_TYPE r = (FLOAT_TYPE) (atom1->position - atom2->position).norm();
          *chg_e += (param1->charge)*(param2->charge)/r;
          const FLOAT_TYPE r_sq = r*r;

          dist_sq_cache[atom1->residue->index][atom1->index][atom2->residue->index][atom2->index] = r_sq;

          const FLOAT_TYPE ratio = (pair_param->vdw_radius_sq) / r_sq;
          const FLOAT_TYPE pow6 = ratio*ratio*ratio;
          *vdw_e += (pair_param->vdw_epsilon)*(pow6*pow6-2*pow6) ;


          const FLOAT_TYPE dist4 = r_sq*r_sq;

          // close contact function correction
          if (r_sq < p5inv*(pair_param->gbsa_radius_sum_sq)) {
               const FLOAT_TYPE born_ratio = r_sq / (pair_param->gbsa_radius_sum_sq);
               const FLOAT_TYPE ccf = 0.25*Math<double>::sqr(1.0 - cos(born_ratio*p5pi));
               (*inv_born_radius)[atom1->residue->index][atom1->index] +=
                    p4*this->volume[atom2->residue->index][atom2->index]*ccf/dist4;
               (*inv_born_radius)[atom2->residue->index][atom2->index] +=
                    p4*this->volume[atom1->residue->index][atom1->index]*ccf/dist4;
          } else {
               // add term to both atom1 and atom2 born radius
               (*inv_born_radius)[atom1->residue->index][atom1->index] +=
                    p4*this->volume[atom2->residue->index][atom2->index]/dist4;
               (*inv_born_radius)[atom2->residue->index][atom2->index] +=
                    p4*this->volume[atom1->residue->index][atom1->index]/dist4;
          }
     }


     //! Evaluate pair interaction between 2 atoms in proximity
     //! \param chain Molecule chain
     //! \return energy in kcal/mol
     FLOAT_TYPE evaluate_nonbonded(ChainFB *chain) {

          FLOAT_TYPE chg_e=0.0,vdw_e=0.0;
          this->counter=0;

          // Initialize Born radius vector with cached covalent neighbour terms
          this->inv_born_radius = this->inv_born_radius_neighbor;
          
          // Add all non-bonded interaction terms
          int size = chain->size();
          for (int i=0; i<size; i++) {
               Residue *res1 = &(*chain)[i];
               int res1_size = res1->size();
          
               // i == j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
                    for (int l=k+1; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         const Parameter<FLOAT_TYPE> *param2 = this->parameters.get_param(atom2);
                         const PairParameter<FLOAT_TYPE> *pair_param =
                              this->parameters.get_pair_param(atom1,atom2);
                    
                         calc_pair_interaction_prox(atom1,atom2,&chg_e,&vdw_e,&this->inv_born_radius,
                                                    pair_param,param1,param2);
                    }
               }
               // i == j+1 also in proximity of res1 (hence use calc_pair_interaction_prox)
               if (i+1 < size) {
                    Residue *res2 = &(*chain)[i+1];
                    int res2_size = res2->size();
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const Parameter<FLOAT_TYPE> *param2 = this->parameters.get_param(atom2);
                              const PairParameter<FLOAT_TYPE> *pair_param =
                                   this->parameters.get_pair_param(atom1,atom2);
                         
                              calc_pair_interaction_prox(
                                   atom1,atom2,&chg_e,&vdw_e,&this->inv_born_radius,pair_param,param1,param2);
                         }
                    }
               }
          
               // rest
               for (int j=i+2; j<size; j++) {
                    Residue *res2 = &(*chain)[j];
                    int res2_size = res2->size();
               
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const Parameter<FLOAT_TYPE> *param2 =
                                   this->parameters.get_param(atom2);
                              const PairParameter<FLOAT_TYPE> *pair_param =
                                   this->parameters.get_pair_param(atom1,atom2);
                         
                              calc_pair_interaction(atom1,atom2,&chg_e,&vdw_e,&this->inv_born_radius,
                                                    pair_param,param1,param2);
                         }
                    }
               }
          }
     
     
          chg_e *= this->parameters.chargeUnit;
     
          // GBSA doesn't necessarily use same unit as the parameter file
          /* const FLOAT_TYPE dielectric = -0.5*332.05; */
          const FLOAT_TYPE chargeUnitInv = -0.00602312; // 1 / (-0.5*332.05382)

          // Calc inverse born radius in inv_born_radius
          for (int r=0; r<size; r++) {
               Residue *res = &(*chain)[r];
               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom1 = res->atoms[a];
                    this->inv_born_radius[atom1->residue->index][atom1->index] *= chargeUnitInv;
                    /* inv_born_radius[atom1->residue->index][atom1->index] = */
                    /*      dielectric/inv_born_radius[atom1->residue->index][atom1->index]; */

                    assert(this->inv_born_radius[atom1->residue->index][atom1->index] > 0.0);
                    // if (this->inv_born_radius[atom1->residue->index][atom1->index] < 0.0) {
                    //      if (this->inv_born_radius[atom1->residue->index][atom1->index] < -1.0)
                    //           std::cerr<<"\nOPLS TERM NONBONDED - "<<atom1<<" has Born radius "
                    //                    <<this->inv_born_radius[atom1->residue->index][atom1->index]
                    //                    <<" which is very bad!"<<std::endl;
                    //      this->inv_born_radius[atom1->residue->index][atom1->index] = 1E-20;
                    // }
               }
          }

          //calculate solvation energy
          FLOAT_TYPE sol_e=this->solvation();
     
          /* std::cout<<"Nonbonded chg  "<<chg_e<<" kcal/mol  "<<counter<<" interactions\n"; */
          /* std::cout<<"Nonbonded vdw  "<<vdw_e<<" kcal/mol  "<<counter<<" interactions\n"; */
          /* std::cout<<"Nonbonded GBSA "<<sol_e<<" kcal/mol hep hep\n"; */
     
          return chg_e+vdw_e+sol_e;
     }


     //! Calculate the polar and non polar part of the solvation energy
     //! \return energy in kcal/mol
     FLOAT_TYPE solvation() {

          const FLOAT_TYPE sigma = 0.0049; // cal/(mol*A^2)
          FLOAT_TYPE e_solv_pol = 0.0;
          FLOAT_TYPE e_solv_npol = 0.0;
     
          // Looping over all pairs of charged atoms
          int size = this->chain->size();
          for (int i=0; i<size; i++) {
               Residue *res1 = &(*this->chain)[i];
               int res1_size = res1->size();
          
               // i == j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    e_solv_npol += this->calc_npol_solv(atom1, this->inv_born_radius[atom1->residue->index][atom1->index]);
                    const FLOAT_TYPE chg1 = this->parameters.get_param(atom1)->charge;
                    if (std::fabs(chg1) < 0.001)
                         continue;
                    for (int l=k; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         const FLOAT_TYPE chg2 = this->parameters.get_param(atom2)->charge;
                         if (std::fabs(chg2) < 0.001)
                              continue;
                         e_solv_pol += calc_pol_solv(atom1,atom2,
                                                     this->inv_born_radius[atom1->residue->index][atom1->index],
                                                     this->inv_born_radius[atom2->residue->index][atom2->index],
                                                     chg1,chg2);
                    }
               }
          
               // j != i
               for (int j=i+1; j<size; j++) {
                    Residue *res2 = &(*this->chain)[j];
                    int res2_size = res2->size();
               
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const FLOAT_TYPE chg1 = this->parameters.get_param(atom1)->charge;
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const FLOAT_TYPE chg2 = this->parameters.get_param(atom2)->charge;
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += calc_pol_solv(atom1,atom2,
                                                          this->inv_born_radius[atom1->residue->index][atom1->index],
                                                          this->inv_born_radius[atom2->residue->index][atom2->index],
                                                          chg1,chg2);
                         }
                    }
               }
          }
          e_solv_npol *= 4.0*M_PI*sigma;
          /* std::cout<<"GBSA: Polar "<<e_solv_pol<<", non-polar "<<e_solv_npol<<"\n"; */
          return e_solv_pol + e_solv_npol;
     }

     //! Allocate space in vectors
     void allocate_vectors() {
          // Allocate distance cache
          int size = (this->chain)->size();
          std::vector<FLOAT_TYPE> vec1;
          std::vector<std::vector<FLOAT_TYPE> > vec2(size,vec1);
          for (int i=0; i<size; i++) {
               vec2[i].resize((*this->chain)[i].size(),0.0);
          }
          std::vector<std::vector<std::vector<FLOAT_TYPE> > > vec3;
          dist_sq_cache.resize(size,vec3);
          for (int i=0; i<size; i++) {
               dist_sq_cache[i].resize((*this->chain)[i].size(),vec2);
          }
     }     

public:

     //! Use same settings as base class
     typedef typename EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;     

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsNonBonded(ChainFB *chain, 
                       const Settings &settings=Settings(),
                       RandomNumberEngine *random_number_engine = &random_global)
          : TermOplsNonBondedBase(chain, "opls-non-bonded", settings, random_number_engine) {
          allocate_vectors();
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsNonBonded(const TermOplsNonBonded &other, 
                       RandomNumberEngine *random_number_engine,
                       int thread_index, ChainFB *chain)
          : TermOplsNonBondedBase(other, random_number_engine, thread_index, chain) {
          allocate_vectors();
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return non-bonded potential energy of the chain in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {
          FLOAT_TYPE energy = evaluate_nonbonded( this->chain );
          return (double) energy;
     }

};


//! GBSA implicit solvation energy term - cached version
template <typename FLOAT_TYPE>
class TermOplsNonBondedCached: public TermOplsNonBondedBase<TermOplsNonBondedCached<FLOAT_TYPE>, 
                                                            FLOAT_TYPE> {
private:

     //! For convenience, define local base class
     typedef phaistos::TermOplsNonBondedBase<TermOplsNonBondedCached<FLOAT_TYPE>, 
                                             FLOAT_TYPE> TermOplsNonBondedBase;     

     //! Type of contribution (charge, vdw, born-radius)
     typedef boost::fusion::vector<FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE> Contribution;

     //! Cached iteration return type
     class NonBondedReturnType: public AtomVectorReturnType<ChainFB,FLOAT_TYPE> {

          //! Local definition of NodeType for ease of reference
          typedef ChainFB::ChainTree::NodeType NodeType;

     public:
          //! Specifies how a contribution is stored in this type - given
          //! the current node in the chaintree and the entity type
          void register_contribution(const NodeType *node, 
                                     const Atom *atom_self, 
                                     const Atom *atom_other, 
                                     const Contribution &contribution,
                                     bool potential_duplicate) {
               this->data[node->get_atom_index(atom_self->atom_type)] += boost::fusion::at_c<2>(contribution);
          }          
     };

     //! Datastructure used by cache: A vector of born radius contributions, a charge sum and a vdw sum
     class ExtendedVector {

     public:

          //! Inverse Born radius contributions
          std::vector<NonBondedReturnType> inv_born_radius_data;

          //! Charge contribution
          FLOAT_TYPE charge_sum;

          //! Van der Waal contribution
          FLOAT_TYPE vdw_sum;

          //! Main constructor
          ExtendedVector(int size)
               : inv_born_radius_data(size),
                 charge_sum(0),
                 vdw_sum(0) {}

          //! Default constructor
          ExtendedVector()
               : charge_sum(0),
                 vdw_sum(0) {}

          //! Overload [] operator
          NonBondedReturnType operator[](const unsigned int index) const {
               return inv_born_radius_data[index];
          }

          //! Overload [] operator
          NonBondedReturnType &operator[](const unsigned int index) {
               return inv_born_radius_data[index];
          }

          //! Return size
          unsigned int size() const {
               return inv_born_radius_data.size();
          }

          //! Resize vector
          void resize(unsigned int size) {
               inv_born_radius_data.resize(size);
          }

     };

     //! Specific CachedIterator for Born radii calculation
     //! This class is basically the same as CachedIteratorVectorBase
     //! but keeps track of the maximum percentwise change in the 
     //! vector of values
     class CachedIteratorNonBonded
          : public chaintree::CachedIteratorVectorBase<CachedIteratorNonBonded,
                                                       ChainFB, Atom, Atom,
                                                       Contribution,
                                                       ExtendedVector >  { 

          //! Define BvType and NodeType locally for ease of reference
          typedef ChainFB::ChainTree::BvType BvType;     
          typedef ChainFB::ChainTree::NodeType NodeType;

          //! Allowed deviation on born radius before node pairs are recalculated
          FLOAT_TYPE deviation_percentage_cutoff;

          //! Distance cutoff used for van der Waals contributions
          FLOAT_TYPE vdw_cutoff;

          //! Distance cutoff used for charge contributions
          FLOAT_TYPE charge_cutoff;

     public:

          //! Datastructure keeping track of at which iteration the change in born radii of the different
          //! nodes were last seen to exceed the cutoff
          std::vector<std::vector<long int> > node_status;

          //! Constructor
          CachedIteratorNonBonded(ChainFB &chain, 
                                  FLOAT_TYPE deviation_percentage_cutoff=0.0,
                                  FLOAT_TYPE vdw_cutoff = std::numeric_limits<FLOAT_TYPE>::infinity(),
                                  FLOAT_TYPE charge_cutoff = std::numeric_limits<FLOAT_TYPE>::infinity())
               : chaintree::CachedIteratorVectorBase<CachedIteratorNonBonded,
                                                     ChainFB, Atom, Atom,
                                                     Contribution,
                                                     ExtendedVector >(chain),
                 deviation_percentage_cutoff(deviation_percentage_cutoff),
                 vdw_cutoff(vdw_cutoff), charge_cutoff(charge_cutoff) {

               // Allocate space for node_status
               node_status.resize(this->ct->get_height());
               for (unsigned int i=0; i<node_status.size(); ++i) {
                    node_status[i].resize(this->ct->get_nodes_at_level(i), 0);
               }
          }

          //! Override the set_cache_entry functionality of the iterator
          //! Registers change in born radii for root node
          //!
          //! \param node1 First node
          //! \param node2 Second node
          //! \param value Value to set in cache
          //! \param fully_initialized_subtree Whether subtree with this node-pair as root is fully initialized          
          void set_cache_entry(NodeType *node1,
                               NodeType *node2,
                               const ExtendedVector &value,
                               bool fully_initialized_subtree=true) {
          
               NodeType *root = this->ct->nodes[this->ct->nodes.size()-1];
               if (node1->level == root->level) {

                    ExtendedVector value_old = this->get_cache_entry(node1, node2).first;
                    for (uint i=0; i<value.size(); ++i) {
                         FLOAT_TYPE current_max = 0;
                         for (uint j=0; j<value[i].size(); ++j) {
                              FLOAT_TYPE deviation_percentage = std::fabs((value[i][j] - value_old[i][j])/(value_old[i][j]));
                              if (deviation_percentage > current_max)
                                   current_max = deviation_percentage;
                         }
                         if (current_max > deviation_percentage_cutoff)
                              set_node_status(this->ct->nodes[i]);
                    }
               }

               // Call base class version
               chaintree::CachedIteratorVectorBase<CachedIteratorNonBonded,
                    ChainFB,Atom,Atom,
                    Contribution, ExtendedVector >::
                    set_cache_entry(node1, node2, value, fully_initialized_subtree);
          }

          //! Called by user during iteration to register new
          //! contribution to the cache
          //!
          //! \param contribution Value to register
          //! \return reference to cache entry in which value is stored
          ExtendedVector &register_contribution(const std::pair<Contribution,
                                                                Contribution> &contribution) {

               this->node_pair_sum.charge_sum += boost::fusion::at_c<0>(contribution.first);
               this->node_pair_sum.vdw_sum += boost::fusion::at_c<1>(contribution.first);
               
               // Call base class version
               return chaintree::CachedIteratorVectorBase<CachedIteratorNonBonded,
                    ChainFB,Atom,Atom,
                    Contribution, ExtendedVector >::
                    register_contribution(contribution);
          }

          //! Copy child value entries in vector to corresponding parent values.
          //! As inherited version, but registering charge and vdw contributions
          //!
          //! \param child1 Child node 1
          //! \param child2 Child node 2
          //! \param parent_pair_vector Destination vector
          void copy_child_pair_to_parent_pair(NodeType *child1, NodeType *child2,
                                              ExtendedVector &parent_pair_vector) {

               parent_pair_vector.charge_sum += this->get_cache_entry(child1, child2).first.charge_sum;
               parent_pair_vector.vdw_sum += this->get_cache_entry(child1, child2).first.vdw_sum;

               // Call base class version
               return chaintree::CachedIteratorVectorBase<CachedIteratorNonBonded,
                    ChainFB,Atom,Atom,
                    Contribution, ExtendedVector >::
                    copy_child_pair_to_parent_pair(child1, child2, parent_pair_vector);
          }

          //! Sets node status for current leaf node and propagates informatation
          //! to internal nodes
          //!
          //! \param node Node object
          void set_node_status(NodeType *node) {
               node_status[node->level][node->index] = this->ct->time;
               while (node->parent) {
                    node = node->parent;
                    if (node_status[node->level][node->index] != this->ct->time) {
                         node_status[node->level][node->index] = this->ct->time;
                    } else {
                         break;
                    }
               }
          }

          //! Overload assignment operator.
          //! This method is inherited automatically by GCC compilers, but needed explicitly by ICC.
          //!
          //! \param other Source object from which assignment is made.
          //! \return Current iterator (this)
          CachedIteratorNonBonded &operator=(const CachedIteratorNonBonded &other) {
               chaintree::CachedIteratorVectorBase<CachedIteratorNonBonded,
                                                   ChainFB, Atom, Atom,
                                                   Contribution,
                                                   ExtendedVector>::operator=(other);
               return *this;
          }
     };


     //! Specific CachedIterator for Born radii calculation
     //! This type is basically the same as CachedIteratorBase
     //! but keeps excludes atom-pairs for which the maximum
     //! percentwise bondradius change is small
     class CachedIteratorSolvent
          : public chaintree::CachedIteratorBase<CachedIteratorSolvent,
                                                 ChainFB, 
                                                 typename ChainFB::ChainTree::NodeType,
                                                 typename ChainFB::ChainTree::NodeType,
                                                 FLOAT_TYPE,FLOAT_TYPE> { 

          //! Define BvType and NodeType locally for ease of reference
          typedef ChainFB::ChainTree::BvType BvType;     
          typedef ChainFB::ChainTree::NodeType NodeType;

          //! This is a reference to the corresponding vector in CachedIteratorNonBonded
          std::vector<std::vector<long int> > &node_status;

     public:

          //! Constructor
          //!
          //! \param chain Molecule chain
          //! \param node_status Datastructure keeping track of at which iteration the change in born radii of the different
          //! nodes were last seen to exceed the cutoff
          CachedIteratorSolvent(ChainFB &chain, 
                          std::vector<std::vector<long int> > &node_status)
               : chaintree::CachedIteratorBase<CachedIteratorSolvent,
                                               ChainFB,
                                               NodeType,
                                               NodeType,
                                               FLOAT_TYPE,
                                               FLOAT_TYPE>(chain),
                 node_status(node_status) {
          }     

          //! Add additional criterion for
          //! excluding child nodes from iteration. In this case, we
          //! exclude a node pair if the maximum born radii change is below cutoff
          //! for both nodes
          //! NOTE that in addition to this criterion, the distance between
          //! the nodes must also be unaltered in order for exclusion to take place
          //!
          //! \param node1 First node
          //! \param node2 Second node
          //! return True if node pair is to be excluded
          bool extra_exclusion_criterion(NodeType *node1, NodeType *node2) {
               if ((node_status[node1->level][node1->index] != this->ct->time) &&
                   (node_status[node2->level][node2->index] != this->ct->time)) {
                    return true;
               }
               return false;
          }

          //! Overload assignment operator.
          //! This method is inherited automatically by GCC compilers, but needed explicitly by ICC.
          //!
          //! \param other Source object from which assignment is made.
          //! \return Current iterator (this)
          CachedIteratorSolvent &operator=(const CachedIteratorSolvent &other) {
               chaintree::CachedIteratorBase<CachedIteratorSolvent,
                                             ChainFB,
                                             NodeType,
                                             NodeType,
                                             FLOAT_TYPE,
                                             FLOAT_TYPE>::operator=(other);
               return *this;
          }
     };


     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;

     //! Return type used by cached iterator
     typedef NonBondedReturnType AtomVector;

     //! Cached iterator for born radii calculation
     CachedIteratorNonBonded cached_it_phase1;

     //! Cached iterator for solvent energi calculation
     CachedIteratorSolvent cached_it_phase2;

     //! Chaintree iterator settings - phase1
     typename chaintree::PairIterator<ChainFB,Atom,Atom>::Settings iterator_settings_phase1;

     //! Chaintree iterator settings - phase2
     typename chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings iterator_settings_phase2;

     //! Evaluate pair interaction between 2 atoms in proximity
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param pair_param Parameters for atom pair
     //! \param param1 Parameters for atom1
     //! \param param2 Parameters for atom2
     //! \param distance Distance between atom1 and atom2
     //! \param chg_e Pointer to charge energy sum
     //! \param vdw_e Pointer to vdw energy sum
     //! \param inv_born_radius_contribution1 Pointer to born radius contribution of atom1
     //! \param inv_born_radius_contribution2 Pointer to born radius contribution of atom2
     void calc_pair_interaction_prox(Atom *atom1, Atom *atom2,
                                     const PairParameter<FLOAT_TYPE> *pair_param,
                                     const Parameter<FLOAT_TYPE> *param1,
                                     const Parameter<FLOAT_TYPE> *param2,
                                     const FLOAT_TYPE distance,
                                     FLOAT_TYPE *chg_e, 
                                     FLOAT_TYPE *vdw_e,
                                     FLOAT_TYPE *inv_born_radius_contribution1,
                                     FLOAT_TYPE *inv_born_radius_contribution2) {

          const double p3 = 6.211;

          int d = this->parameters.template get_covalent_dist<ChainFB>(atom1,atom2);
          if (d > 3) {
               calc_pair_interaction(
                    atom1,atom2,
                    pair_param,param1,param2,
                    distance,
                    chg_e,vdw_e,
                    inv_born_radius_contribution1, inv_born_radius_contribution2);
          } else if (d == 2 ) {

               FLOAT_TYPE dist2 = distance*distance;
               FLOAT_TYPE dist4 = dist2*dist2;

               *inv_born_radius_contribution1 = (p3*this->volume[atom2->residue->index][atom2->index])/dist4;
               *inv_born_radius_contribution2 = (p3*this->volume[atom1->residue->index][atom1->index])/dist4;

          } else if (d==3) {//d==3
               FLOAT_TYPE chg_ep=0.0,vdw_ep=0.0;
               calc_pair_interaction(atom1, atom2, pair_param, param1, param2, 
                                     distance,
                                     &chg_ep, &vdw_ep, 
                                     inv_born_radius_contribution1, inv_born_radius_contribution2 );
               *chg_e += chg_ep*this->parameters.chg14scale;
               *vdw_e += vdw_ep*this->parameters.vdw14scale;
          }
     }

     //! Evaluate pair interaction
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param param1 Parameters for atom1
     //! \param param2 Parameters for atom2
     //! \param distance Distance between atom1 and atom2
     //! \param chg_e Pointer to charge energy sum
     //! \param vdw_e Pointer to vdw energy sum
     //! \param inv_born_radius_contribution1 Pointer to born radius contribution of atom1
     //! \param inv_born_radius_contribution2 Pointer to born radius contribution of atom2
     //! \return energy in kcal/mol
     void calc_pair_interaction(Atom *atom1, Atom *atom2,
                                const PairParameter<FLOAT_TYPE> *pair_param,
                                const Parameter<FLOAT_TYPE> *param1,
                                const Parameter<FLOAT_TYPE> *param2,
                                const FLOAT_TYPE distance,
                                FLOAT_TYPE *chg_e, 
                                FLOAT_TYPE *vdw_e,
                                FLOAT_TYPE *inv_born_radius_contribution1,
                                FLOAT_TYPE *inv_born_radius_contribution2) {

          const FLOAT_TYPE p4=15.236;
          const FLOAT_TYPE p5pi=3.939441;  // = 1.254*3.1415
          const FLOAT_TYPE p5inv=0.797448; // = 1.0/1.254

          // const FLOAT_TYPE r = (FLOAT_TYPE) (atom1->position - atom2->position).norm();
          if (distance < settings.charge_cutoff_distance) {
               *chg_e = (param1->charge)*(param2->charge)/distance;
          } else {
               *chg_e = 0.0;
          }

          const FLOAT_TYPE dist2 = distance*distance;

          if (distance < settings.vdw_cutoff_distance) {
               const FLOAT_TYPE ratio = (pair_param->vdw_radius_sq) / dist2;
               const FLOAT_TYPE pow6 = ratio*ratio*ratio;
               *vdw_e = (pair_param->vdw_epsilon)*(pow6*pow6-2*pow6) ;
          } else {
               *vdw_e = 0.0;
          }

          const FLOAT_TYPE dist4 = dist2*dist2;

          // close contact function correction
          if (dist2 < p5inv*(pair_param->gbsa_radius_sum_sq)) {
               const FLOAT_TYPE born_ratio = dist2 / (pair_param->gbsa_radius_sum_sq);
               const FLOAT_TYPE ccf = 0.25*Math<double>::sqr(1.0 - cos(born_ratio*p5pi));

               *inv_born_radius_contribution1 = (p4*this->volume[atom2->residue->index][atom2->index]*ccf)/dist4;
               *inv_born_radius_contribution2 = (p4*this->volume[atom1->residue->index][atom1->index]*ccf)/dist4;
          } else {
               // add term to both atom1 and atom2 born radius
               *inv_born_radius_contribution1 = (p4*this->volume[atom2->residue->index][atom2->index])/dist4;
               *inv_born_radius_contribution2 = (p4*this->volume[atom1->residue->index][atom1->index])/dist4;
          }
     }


     //! Evaluate pair interactions
     FLOAT_TYPE evaluate_nonbonded() {

          for (cached_it_phase1(*this->chain, iterator_settings_phase1); 
               !cached_it_phase1.end(); ++cached_it_phase1) {

               Atom *atom1 = cached_it_phase1->first;
               Atom *atom2 = cached_it_phase1->second;

               const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
               const Parameter<FLOAT_TYPE> *param2 = this->parameters.get_param(atom2);
               const PairParameter<FLOAT_TYPE> *pair_param =
                    this->parameters.get_pair_param(atom1,atom2);

               std::pair<Contribution,Contribution> contribution;
               
               FLOAT_TYPE &charge_energy = boost::fusion::at_c<0>(contribution.first);
               FLOAT_TYPE &vdw_energy = boost::fusion::at_c<1>(contribution.first);
               FLOAT_TYPE &inv_born_radius_contribution1 = boost::fusion::at_c<2>(contribution.first);
               FLOAT_TYPE &inv_born_radius_contribution2 = boost::fusion::at_c<2>(contribution.second);

               if (std::fabs(atom1->residue->index - atom2->residue->index) <= 1) {
                    calc_pair_interaction_prox(atom1,atom2,
                                               pair_param,
                                               param1,param2,
                                               cached_it_phase1->distance,
                                               &charge_energy,
                                               &vdw_energy,
                                               &inv_born_radius_contribution1,
                                               &inv_born_radius_contribution2);
               } else {
                    calc_pair_interaction(atom1,atom2,
                                          pair_param,
                                          param1,param2,
                                          cached_it_phase1->distance,
                                          &charge_energy,
                                          &vdw_energy,
                                          &inv_born_radius_contribution1,
                                          &inv_born_radius_contribution2);
               }

               // Register with cached iterator
               cached_it_phase1.register_contribution(contribution);
          }
          // Retrieve result from cached iterator
          ExtendedVector &sum_cached = cached_it_phase1.compute_total();



          // Post iteration calculations

          // GBSA doesn't necessarily use same unit as the parameter file
          const FLOAT_TYPE charge_unit_inv = -0.00602312; // 1 / (-0.5*332.05382)
          // const double dielectric = -0.5*332.05;     

          int size = this->chain->chain_tree->get_nodes_at_level(0);
          this->inv_born_radius.resize(size);
          for (int i=0; i<size; ++i) {
               this->inv_born_radius[i].resize(this->chain->chain_tree->nodes[i]->size());
               for (unsigned int k=0; k<this->chain->chain_tree->nodes[i]->size(); ++k) {
                    Atom *atom1 = this->chain->chain_tree->nodes[i]->atoms[k];
                    int atom1_index = this->chain->chain_tree->nodes[i]->get_atom_index(atom1->atom_type);
                    this->inv_born_radius[i][atom1_index] = 
                         (this->inv_born_radius_neighbor[atom1->residue->index][atom1->index] +
                          sum_cached[i][atom1_index])*charge_unit_inv;
                    if (this->inv_born_radius[i][atom1_index] < 0) {
                         std::cerr << "WARNING: GBSA: born radius<0: " << this->inv_born_radius[i][atom1_index] << ". Setting to ";
                         this->inv_born_radius[i][atom1_index] = std::numeric_limits<double>::infinity();
                         // this->inv_born_radius[i][atom1_index] = 1.0/1E-20;
                         std::cerr << this->inv_born_radius[i][atom1_index] << "\n";
                         std::cerr.flush();
                    }
               }
          }

          return sum_cached.vdw_sum + sum_cached.charge_sum*this->parameters.chargeUnit;
     }

     //! Evaluate pair interactions.
     //! \param Molecule chain object
     //! \param inv_born_radius Inverse born radii of all atoms. [res index][atoms index in res]
     FLOAT_TYPE evaluate_nonbonded_uncached(ChainFB *chain,
                                            std::vector<std::vector<FLOAT_TYPE> > &inv_born_radius) {

          FLOAT_TYPE charge_energy = 0.0;
          FLOAT_TYPE vdw_energy = 0.0;
          
          // Add all non-bonded interaction terms
          int size = chain->size();
          for (int i=0; i<size; i++) {
               Residue *res1 = &(*chain)[i];
               int res1_size = res1->size();
          
               // i == j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
                    for (int l=k+1; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         const Parameter<FLOAT_TYPE> *param2 = this->parameters.get_param(atom2);
                         const PairParameter<FLOAT_TYPE> *pair_param =
                              this->parameters.get_pair_param(atom1,atom2);
                    
                         calc_pair_interaction_prox(atom1,atom2,&charge_energy,&vdw_energy,&inv_born_radius,
                                                    pair_param,param1,param2);
                    }
               }
               // i == j+1 also in proximity of res1 (hence use calc_pair_interaction_prox)
               if (i+1 < size) {
                    Residue *res2 = &(*chain)[i+1];
                    int res2_size = res2->size();
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const Parameter<FLOAT_TYPE> *param2 = this->parameters.get_param(atom2);
                              const PairParameter<FLOAT_TYPE> *pair_param =
                                   this->parameters.get_pair_param(atom1,atom2);
                         
                              calc_pair_interaction_prox(
                                   atom1,atom2,&charge_energy,&vdw_energy,&inv_born_radius,pair_param,param1,param2);
                         }
                    }
               }
          
               // rest
               for (int j=i+2; j<size; j++) {
                    Residue *res2 = &(*chain)[j];
                    int res2_size = res2->size();
               
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const Parameter<FLOAT_TYPE> *param1 = this->parameters.get_param(atom1);
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const Parameter<FLOAT_TYPE> *param2 =
                                   this->parameters.get_param(atom2);
                              const PairParameter<FLOAT_TYPE> *pair_param =
                                   this->parameters.get_pair_param(atom1,atom2);
                         
                              calc_pair_interaction(atom1,atom2,&charge_energy,&vdw_energy,&inv_born_radius,
                                                    pair_param,param1,param2);
                         }
                    }
               }
          }
     
     
          charge_energy *= this->parameters.chargeUnit;
     
          // GBSA doesn't necessarily use same unit as the parameter file
          /* const FLOAT_TYPE dielectric = -0.5*332.05; */
          const FLOAT_TYPE chargeUnitInv = -0.00602312; // 1 / (-0.5*332.05382)
     
          // Calc inverse born radius in inv_born_radius
          for (int r=0; r<size; r++) {
               Residue *res = &(*chain)[r];
               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom1 = res->atoms[a];
                    inv_born_radius[atom1->residue->index][atom1->index] *= chargeUnitInv;
                    /* inv_born_radius[atom1->residue->index][atom1->index] = */
                    /*      dielectric/inv_born_radius[atom1->residue->index][atom1->index]; */
                    if (inv_born_radius[atom1->residue->index][atom1->index] < 0.0) {
                         if (inv_born_radius[atom1->residue->index][atom1->index] < -1.0)
                              std::cerr<<"\nOPLS TERM NONBONDED - "<<atom1<<" has Born radius "
                                       <<inv_born_radius[atom1->residue->index][atom1->index]
                                       <<" which is very bad!"<<std::endl;
                         inv_born_radius[atom1->residue->index][atom1->index] = 1E-20;
                    }
               }
          }
     
          return charge_energy + vdw_energy;
     }


     //! Calculate pairwise polar solvation energy
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param born_radius1 Born radius of atom1
     //! \param born_radius2 Born radius of atom2
     //! \param chg1 partial charge of atom 1 in electron volt
     //! \param chg2 partial charge of atom 2 in electron volt
     //! \return energy in kcal/mol
     inline FLOAT_TYPE calc_pol_solv(Atom *atom1, Atom *atom2,
                                     FLOAT_TYPE born_radius1, FLOAT_TYPE born_radius2,
                                     const FLOAT_TYPE chg1, const FLOAT_TYPE chg2) {

          /* const FLOAT_TYPE dielectric = 332.05; */
          /* const FLOAT_TYPE epsilon = 78.3; */
          /* const FLOAT_TYPE fac = (-dielectric*(1.0 - 1.0/epsilon)); */
          const FLOAT_TYPE fac=-327.809;

          const FLOAT_TYPE dx = (FLOAT_TYPE) atom1->position[0] - atom2->position[0];
          const FLOAT_TYPE dy = (FLOAT_TYPE) atom1->position[1] - atom2->position[1];
          const FLOAT_TYPE dz = (FLOAT_TYPE) atom1->position[2] - atom2->position[2];
          const FLOAT_TYPE dist2 = dx*dx+dy*dy+dz*dz;
          // const FLOAT_TYPE dist2 = dist_sq_cache[atom1->residue->index][atom1->index][atom2->residue->index][atom2->index];
     
          FLOAT_TYPE c_ij = fac*chg1*chg2;
          if (atom1 == atom2)
               c_ij *= 0.5;

          const FLOAT_TYPE r_born12 = born_radius1*born_radius2;

          /* float x = dist2+expo(-(dist2*0.25*r_born12))/r_born12; */
          /* const FLOAT_TYPE f_ij = fastInvSqrt(x); */

          /* FLOAT_TYPE x = dist2+expo(-(dist2*0.25*r_born12))/r_born12; */
          /* const FLOAT_TYPE f_ij = invSqrt(x); */

          const FLOAT_TYPE f_ij = 1.0/sqrt(dist2+std::exp(-(dist2*0.25*r_born12))/r_born12);

          return c_ij*f_ij;
     }


     //! Calculate the polar and non polar part of the solvation energy
     double solvation() {

          const double sigma = 0.0049; // cal/(mol*A^2)

          for (cached_it_phase2(*this->chain, iterator_settings_phase2); 
               !cached_it_phase2.end(); ++cached_it_phase2) {

               NodeType *node1 = cached_it_phase2->first;
               NodeType *node2 = cached_it_phase2->second;

               double e_solv_pol = 0.0;
               double e_solv_npol = 0.0;

               // Special case for identical nodes
               if (node1 == node2) {
                    for (unsigned int k=0; k<node1->size(); ++k) {
                         Atom *atom1 = node1->atoms[k];
                         int atom1_index = node1->get_atom_index(atom1->atom_type);

                         e_solv_npol += this->calc_npol_solv(atom1, this->inv_born_radius[node1->index][atom1_index]);
                         const FLOAT_TYPE chg1 = this->parameters.get_param(atom1)->charge;
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (unsigned int l=k; l<node1->size(); l++) {
                              Atom *atom2 = node1->atoms[l];
                              int atom2_index = node1->get_atom_index(atom2->atom_type);
                              const FLOAT_TYPE chg2 = this->parameters.get_param(atom2)->charge;
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += this->calc_pol_solv(atom1,atom2,
                                                                this->inv_born_radius[node1->index][atom1_index],
                                                                this->inv_born_radius[node1->index][atom2_index],
                                                                chg1,chg2);
                         }                    
                    }
               } else {
                    for (unsigned int k=0; k<node1->size(); ++k) {
                         Atom *atom1 = node1->atoms[k];
                         int atom1_index = node1->get_atom_index(atom1->atom_type);
                         const FLOAT_TYPE chg1 = this->parameters.get_param(atom1)->charge;
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (unsigned int l=0; l<node2->size(); ++l) {
                              Atom *atom2 = node2->atoms[l];
                              int atom2_index = node2->get_atom_index(atom2->atom_type);
                              const FLOAT_TYPE chg2 = this->parameters.get_param(atom2)->charge;
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += this->calc_pol_solv(atom1,atom2,
                                                                this->inv_born_radius[node1->index][atom1_index],
                                                                this->inv_born_radius[node2->index][atom2_index],
                                                                chg1,chg2);
                         }
                    }
               }

               e_solv_npol *= 4.0*M_PI*sigma;

               // Register contribution with cache
               cached_it_phase2.register_contribution(e_solv_pol + e_solv_npol);
          }

          // Retrieve total value from cache
          double energy = cached_it_phase2.compute_total();
          return energy;
     }

     // //! Calculate the polar and non polar part of the solvation energy
     // //! This method is provided for testing purposes only. It calculates the born radii in exactly 
     // //! the same way as the cached version, but without caching.
     // FLOAT_TYPE solvation_uncached() {

     //      const FLOAT_TYPE sigma = 0.0049; // cal/(mol*A^2)
     //      FLOAT_TYPE e_solv_pol = 0.0;
     //      FLOAT_TYPE e_solv_npol = 0.0;

     //      // Looping over all pairs of charged atoms
     //      int size = this->chain->chain_tree->get_nodes_at_level(0);
     //      for (int i=0; i<size; ++i) {

     //           NodeType *node1 = this->chain->chain_tree->nodes[i];

     //           // i == j only upper triangle:
     //           for (unsigned int k=0; k<node1->size(); ++k) {
     //                Atom *atom1 = node1->atoms[k];
     //                int atom1_index = node1->get_atom_index(atom1->atom_type);

     //                e_solv_npol += this->calc_npol_solv(atom1, this->inv_born_radius[node1->index][atom1_index]);
     //                const FLOAT_TYPE chg1 = this->parameters.get_param(atom1)->charge;
     //                if (fabs(chg1) < 0.001)
     //                     continue;
     //                for (unsigned int l=k; l<node1->size(); l++) {
     //                     Atom *atom2 = node1->atoms[l];
     //                     int atom2_index = node1->get_atom_index(atom2->atom_type);
     //                     const FLOAT_TYPE chg2 = this->parameters.get_param(atom2)->charge;
     //                     if (fabs(chg2) < 0.001)
     //                          continue;
     //                     e_solv_pol += this->calc_pol_solv(atom1,atom2,
     //                                                       this->inv_born_radius[node1->index][atom1_index],
     //                                                       this->inv_born_radius[node1->index][atom2_index],
     //                                                       chg1,chg2);
     //                }                    
     //           }

     //           // j != i
     //           for (int j=i+1; j<size; ++j) {

     //                NodeType *node2 = this->chain->chain_tree->nodes[j];

     //                for (unsigned int k=0; k<node1->size(); ++k) {
     //                     Atom *atom1 = node1->atoms[k];
     //                     int atom1_index = node1->get_atom_index(atom1->atom_type);
     //                     /* e_solv_npol += calc_npol_solv(atom1); */
     //                     const FLOAT_TYPE chg1 = this->parameters.get_param(atom1)->charge;
     //                     if (fabs(chg1) < 0.001)
     //                          continue;
     //                     for (unsigned int l=0; l<node2->size(); ++l) {
     //                          Atom *atom2 = node2->atoms[l];
     //                          int atom2_index = node2->get_atom_index(atom2->atom_type);
     //                          const FLOAT_TYPE chg2 = this->parameters.get_param(atom2)->charge;
     //                          if (fabs(chg2) < 0.001)
     //                               continue;
     //                          // std::cout << atom1->atom_type << " " << atom2->atom_type << " " << this->calc_pol_solv(atom1,atom2,
     //                          //                                   this->inv_born_radius[node1->index][atom1_index],
     //                          //                                   this->inv_born_radius[node2->index][atom2_index],
     //                          //                                  chg1,chg2) << "\n";
     //                          e_solv_pol += this->calc_pol_solv(atom1,atom2,
     //                                                            this->inv_born_radius[node1->index][atom1_index],
     //                                                            this->inv_born_radius[node2->index][atom2_index],
     //                                                            chg1,chg2);
     //                     }
     //                }
     //           }               
     //      }

     //      e_solv_npol *= 4.0*M_PI*sigma;
     //      return e_solv_pol + e_solv_npol;
     // }

public:

     //! Settings
     const class Settings: public TermOplsNonBondedBase::Settings {

     public:

          //! van der Waals interaction cutoff distance
          FLOAT_TYPE vdw_cutoff_distance;

          //! charge interaction cutoff distance
          FLOAT_TYPE charge_cutoff_distance;

          //! Maximum deviation allowed in born radii in two subtrees of the chaintree before it is recalculated
          FLOAT_TYPE gbsa_maximum_deviation_cutoff;

          //! Distance beyond which gbsa contributions are set to zero in phase1.
          FLOAT_TYPE gbsa_cutoff_distance_phase2;

          //! Distance beyond which gbsa contributions are set to zero in phase2.
          FLOAT_TYPE gbsa_cutoff_distance_phase1;

          //! Constructor
          Settings(FLOAT_TYPE vdw_cutoff_distance = std::numeric_limits<FLOAT_TYPE>::infinity(),
                   FLOAT_TYPE charge_cutoff_distance = std::numeric_limits<FLOAT_TYPE>::infinity(),
                   FLOAT_TYPE gbsa_maximum_deviation_cutoff = 0.1,
                   FLOAT_TYPE gbsa_cutoff_distance_phase2 = std::numeric_limits<FLOAT_TYPE>::infinity(),
                   FLOAT_TYPE gbsa_cutoff_distance_phase1 = std::numeric_limits<FLOAT_TYPE>::infinity())
               : vdw_cutoff_distance(vdw_cutoff_distance),
                 charge_cutoff_distance(charge_cutoff_distance),
                 gbsa_maximum_deviation_cutoff(gbsa_maximum_deviation_cutoff),
                 gbsa_cutoff_distance_phase2(gbsa_cutoff_distance_phase2),
                 gbsa_cutoff_distance_phase1(gbsa_cutoff_distance_phase1) {}

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "vdw-cutoff-distance:" << settings.vdw_cutoff_distance << "\n";
               o << "charge-cutoff-distance:" << settings.charge_cutoff_distance << "\n";
               o << "gbsa-maximum-deviation-cutoff:" << settings.gbsa_maximum_deviation_cutoff << "\n";
               o << "gbsa-cutoff-distance-phase1:" << settings.gbsa_cutoff_distance_phase1 << "\n";
               o << "gbsa-cutoff-distance-phase2:" << settings.gbsa_cutoff_distance_phase2 << "\n";
               o << (typename TermOplsNonBondedBase::Settings &)settings;
               return o;
          }                    
     } settings;     


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsNonBondedCached(ChainFB *chain, const Settings &settings=Settings(),
                             RandomNumberEngine *random_number_engine = &random_global)
          : TermOplsNonBondedBase(chain, "opls-non-bonded-cached", settings, random_number_engine),
            cached_it_phase1(*chain, settings.gbsa_maximum_deviation_cutoff, 
                             settings.vdw_cutoff_distance, settings.charge_cutoff_distance),
            cached_it_phase2(*chain, cached_it_phase1.node_status),
            settings(settings) {

          this->inv_born_radius.clear();

          // Only evaluate modified pairs
          bool only_modified_pairs = true;
          iterator_settings_phase1 = typename chaintree::PairIterator<ChainFB,Atom,Atom>::Settings(settings.gbsa_cutoff_distance_phase1,
                                                                                                      only_modified_pairs);
          // Evaluate all pairs because born radii changes occur non-locally
          only_modified_pairs = false;
          iterator_settings_phase2 = typename chaintree::PairIterator<ChainFB,NodeType,NodeType>::
               Settings(settings.gbsa_cutoff_distance_phase2,
                        only_modified_pairs);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsNonBondedCached(const TermOplsNonBondedCached &other, 
                             RandomNumberEngine *random_number_engine,
                             int thread_index, ChainFB *chain)
          : TermOplsNonBondedBase(other, random_number_engine, thread_index, chain),
            cached_it_phase1(*chain, other.settings.gbsa_maximum_deviation_cutoff, 
                             other.settings.vdw_cutoff_distance, other.settings.charge_cutoff_distance),
            cached_it_phase2(*chain, cached_it_phase1.node_status),
            iterator_settings_phase1(other.iterator_settings_phase1),
            iterator_settings_phase2(other.iterator_settings_phase2),
            settings(other.settings) {}


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          // Initialize Born radius vector with cached covalent neighbour terms
          std::vector<std::vector<FLOAT_TYPE> > inv_born_radius = this->inv_born_radius_neighbor;

          FLOAT_TYPE charge_vdw_energy = evaluate_nonbonded();

          //calculate solvation energy
          // FLOAT_TYPE solvation_energy = solvation_uncached();
          FLOAT_TYPE solvation_energy = solvation();

          // Cast to double
          return (double) charge_vdw_energy + solvation_energy;
     }

     //! Accept last energy evaluation
     void accept() {
          cached_it_phase1.accept();
          cached_it_phase2.accept();
     }

     //! Reject last energy evaluation
     void reject() {
          cached_it_phase1.reject();
          cached_it_phase2.reject();
     }
};

}
#endif
