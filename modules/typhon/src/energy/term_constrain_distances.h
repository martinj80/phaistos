// constrain_distances.h --- allows to define different distance constrains
//
// Copyright (C) 2011  Tim Harder, Martin Paluszewski, Thomas Hamelryck, Wouter Boomsma
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

#ifndef CONSTRAIN_DISTANCES_H
#define CONSTRAIN_DISTANCES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#include "protein/atom.h"
#include "energy/energy_term.h"
#include "protein/definitions.h"
#include "protein/iterators/residue_iterator.h"


namespace phaistos {
namespace constrain_distances {

//@{
//! DSSP constants (Kabsch W, Sander C, Biopolymers, 1983)
const static double dssp_qee = 332.0*0.42*0.20;
const static double dssp_min_HO_dssp = 1.7;
//@}

//! Define interaction types
enum InteractionType {NO_INTERACTION=0, BB_HBOND=1, SC_HBOND=2, BB_SC_HBOND=3, DISTANCE_GAUSSIAN=4, 
                      DISTANCE_GAUSSIAN_FIXED_POINT=5, SS_BOND=6, LOG_NORMAL=7, ANY_HBOND=8};

//! Gaussian parameter to model HBonds
const static double parameter_hbond_native_sd = 0.15;

//! top = estimated from the Top 500 dataset
const static double parameter_hbond_NO_top_mean = 3.023;
const static double parameter_hbond_NO_top_sd = 0.20;
const static double parameter_hbond_NC_top_mean = 4.056;
const static double parameter_hbond_NC_top_sd = 0.24;
const static double parameter_hbond_HO_top_mean = 2.1;
const static double parameter_hbond_HO_top_sd = 0.24;
const static double parameter_hbond_HC_top_mean = 3.2;
const static double parameter_hbond_HC_top_sd = 0.235;


//! class dummy atom
//!
//! This class provides a simple container to store
//! the essential parameters of an atom, such as the
//! atom name, coordinates and position in the chain.
class DummyAtom {
public:

     //! Name of atom
     std::string atom_name;
     
     //! Position (vector of coordinates)
     Vector_3D position;

     //! Position in the chain
     int index;

     //! An id
     int id;

     //! Contructor
     //!
     //! \param atom_name Atom Name
     //! \param position Vector of the coordinates
     //! \param index position in the chain
     //! \param id an ID
     DummyAtom (std::string atom_name, Vector_3D position, 
                int index, int id=0)
          : atom_name(atom_name), position(position), 
            index(index), id(id){}
};


//! Distance between two atoms
//! \param first First atom
//! \param second Second atom
inline double get_distance(Atom* first, Atom* second) {
     return (first->position - second->position).norm();
}

//! Distance between atom and dummy atom
//! \param first Dummy atom
//! \param second normal atom
inline double get_distance(DummyAtom* first, Atom* second) {
     return (first->position - second->position).norm();
}



//! Gaussian Normal Distribution
//!
//! Simple Gaussian normal distribution implementation.
class Gaussian {
public:

     //! Mean value
     double mean;

     //! Standard deviation
     double sd;

     //! Normalization constant
     double gauss_const;

     //! standard constructor
     //!
     //! \param mean set the mean
     //! \param sd set the standard deviation
     Gaussian(double mean = -1., double sd = -1.) {
          this->mean = mean;
          this->sd = sd;
          this->gauss_const = 0.5 * std::log(2 * M_PI * this->sd * this->sd);
     }

     //! get_energy
     //!
     //! returns the negative log likelihood for a given distance
     //! !param distance value to calculate the likelihood for
     double get_energy(double distance) {
          double dens = -(distance - this->mean) * (distance - this->mean) / (2 * this->sd * this->sd) - this->gauss_const;
          assert(std::isfinite(dens));
          return -dens;
     }
};

//! Atom pair class
//!
//! This class is the base class to the more specialized restraint classes
//! such as the Gaussian contact or the hydrogen or disulfide bond.
//! It stores all important information describing the atom
//! pair, plus it implements a simple caching functionality.
template<typename CHAIN_TYPE, typename ATOM_TYPE1=Atom, typename ATOM_TYPE2=Atom>
class AtomPair {
public:

     //! First atom
     ATOM_TYPE1 *first;

     //! Second atom
     ATOM_TYPE2 *second;

     //! Index of first atom (for convenience)
     int index1;

     //! Index of second atom (for convenience)
     int index2;

     //! How solvated the interaction is
     int desolvation_factor;

     //! Whether the interaction was updated
     bool was_updated;

     //! Energy contribution backup
     double contribution_cache;

     //! Energy contribution
     double contribution;

     //! Type of interaction
     InteractionType interaction_type;

     //! Whether the interaction has been initialized
     bool initialized;

     //! Constructor
     //!
     //! Creates an uninitialized, empty atom pair. Note, that all
     //! members need to be set manually afterwards.
     AtomPair(InteractionType interaction_type=NO_INTERACTION)
          : first(NULL),
            second(NULL),
            index1(0),
            index2(0),
            desolvation_factor(-1),
            was_updated(false),
            contribution_cache(0.0),
            contribution(0.0),
            interaction_type(interaction_type),
            initialized(false) {}

     //! Constructor - specific for Atom-Atom interaction
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param interaction_type specify the interaction type.
     AtomPair(Atom* first, Atom* second, 
              InteractionType interaction_type = BB_HBOND)
          : first(first),
            second(second),
            index1(first->residue->index),
            index2(second->residue->index),
            desolvation_factor(-1),
            was_updated(false),
            contribution_cache(0.0),
            contribution(0.0),
            interaction_type(interaction_type),
            initialized(true) {}

     //! Constructor - specific for DummyAtom-Atom interaction
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param interaction_type specify the interaction type.
     AtomPair(DummyAtom* first, Atom* second, 
              InteractionType interaction_type = BB_HBOND)
          : first(first),
            second(second),
            index1(-1),
            index2(second->residue->index),
            desolvation_factor(-1),
            was_updated(false),
            contribution_cache(0.0),
            contribution(0.0),
            interaction_type(interaction_type),
            initialized(true) {}
 
     //! Constructor
     //!
     //! \param chain pointer to the chain object, to extract the atoms from
     //! \param index1 residue index of the first atom in the chain
     //! \param index2 residue index of the second atom in the chain.
     //! \param interaction_type specify the interaction type
     AtomPair(CHAIN_TYPE *chain, int index1, int index2, 
              InteractionType interaction_type = BB_HBOND)
          : index1(index1),
            index2(index2),
            desolvation_factor(-1),
            was_updated(false),
            contribution_cache(0.0),
            contribution(0.0),
            interaction_type(interaction_type),
            initialized(false) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          Residue *rf = &(*chain)[index1];
          Residue *rs = &(*chain)[index2];

          if (this->interaction_type == BB_HBOND) {
               // first is pointing to the H atom, second to the O
               if (rf->has_atom(H) && rs->has_atom(O)) {
                    this->first = (*rf)[H];
                    this->second = (*rs)[O];
               } else {
                    std::cerr << "ERROR: missing atoms at H: " << this->index1 << " or O: " << this->index2 << ". Bailing out. \n";
                    assert(false);
               }
          } else if (this->interaction_type == DISTANCE_GAUSSIAN) {
               // first is pointing to the H atom, second to the O
               if (rf->has_atom(CA) && rs->has_atom(CA)) {
                    this->first = (*rf)[CA];
                    this->second = (*rs)[CA];
               } else {
                    std::cerr << "ERROR: missing atoms at CA: " << this->index1 << " or CA: " << this->index2 << ". Bailing out. \n";
                    assert(false);
               }
          } else if (this->interaction_type == SS_BOND) {
               // Assuming CYS residue
               if (rf->has_atom(SG) && rs->has_atom(SG)) {
                    this->first = (*rf)[SG];
                    this->second = (*rs)[SG];
               } else {
                    std::cerr << "ERROR: missing atoms at SG: " << this->index1 << " or SG: " << this->index2 << ". Bailing out. \n";
                    assert(false);
               }
          } else {
               std::cerr << "ERROR: Unknown interaction type detected. Bailing out.\n";
               assert(false);
          }

          this->initialized = true;
     }

     //! Copy constructor
     //!
     //! \param other AtomPair instance to create a copy from
     AtomPair(const AtomPair &other)
          : first(other.first),
            second(other.second),
            index1(other.index1),
            index2(other.index2),
            desolvation_factor(other.desolvation_factor),
            was_updated(other.was_updated),
            contribution_cache(other.contribution_cache),
            contribution(other.contribution),
            interaction_type(other.interaction_type),
            initialized(other.initialized) {}

     //! Check whether object is initialized
     bool is_initialized() {
          return this->initialized;
     }

     //! Determine whether one of the atoms partners
     //! was within the current move range.
     //!
     //! \param moveInfo moveInfo object from the MCMC class,
     //! containing start and end points of the current move
     bool in_range(MoveInfo *moveInfo = NULL) {

          if (moveInfo->move_type == definitions::NON_LOCAL)
               return true;
          if (moveInfo != NULL && ((this->index1 >= moveInfo->modified_angles_start && 
                                    this->index1 < moveInfo->modified_angles_end) || 
                                   (this->index2 >= moveInfo->modified_angles_start && 
                                    this->index2 < moveInfo->modified_angles_end))) {
               return true;
          }
          return false;
     }

     //! Accept the move by invalidating the cache
     void accept() {
          this->contribution_cache = 0.;
     };

     //! Rejecting the move and rolling back to the cached values.
     void reject() {
          if (this->contribution_cache != 0.) {
               this->contribution = this->contribution_cache;
          }
          this->contribution_cache = 0.;
     };

     //! Clear cache
     //!
     //! Invalidate the current cache content.
     void clear_cache() {
          this->contribution = 0.;
          this->contribution_cache = 0.;
     };
};


//! Atom pair - hydrogen bond
//!
//! More specific implementation of the atom pair class.
//! this pair contains the four atoms part of a hydrogen
//! bond, namely the C, the O, the N and the H.
//! the four mutual distances NC, NO, HC and HO are being
//! evaluated.
template<typename CHAIN_TYPE>
class AtomPairHBond: public AtomPair<CHAIN_TYPE,Atom> {
public:

     //! Pointer to Nitrogen atom
     Atom *atomN;

     //! Pointer to Hydrogen atom
     Atom *atomH;

     //! Pointer to Carbon atom
     Atom *atomC;

     //! Pointer to Oxygen atom
     Atom *atomO;

     //! Gaussian: N-O distance
     Gaussian *NO;

     //! Gaussian: N-C distance
     Gaussian *NC;

     //! Gaussian: H-O distance
     Gaussian *HO;

     //! Gaussian: H-C distance
     Gaussian *HC;

     //! Initialize object
     //!
     //! \param init_from_native Whether to initialize gaussians from original structure
     void init(bool init_from_native=false) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          this->atomH = this->first;
          this->atomO = this->second;
          this->atomN = NULL;
          this->atomC = NULL;

          if (this->atomH->residue->has_atom(N))
               this->atomN = (*(this->atomH->residue))[N];

          if (this->atomO->residue->has_atom(C))
               this->atomC = (*(this->atomO->residue))[C];

          if (atomN == NULL || atomC == NULL) {
               // couldn't find all the atoms ..
               this->initialized = false;
          }

          // Either Gaussians are initialized with the distances in the starting structure
          if (init_from_native) {
               double distanceNO = (this->atomN->position - this->atomO->position).norm();
               double distanceNC = (this->atomN->position - this->atomC->position).norm();
               double distanceHO = (this->atomH->position - this->atomO->position).norm();
               double distanceHC = (this->atomH->position - this->atomC->position).norm();

               this->NO = new Gaussian(distanceNO, parameter_hbond_native_sd);
               this->NC = new Gaussian(distanceNC, parameter_hbond_native_sd);
               this->HO = new Gaussian(distanceHO, parameter_hbond_native_sd);
               this->HC = new Gaussian(distanceHC, parameter_hbond_native_sd);

          // Or using values from the Top500 dataset
          } else {
               /* Values for the Gaussians derrived
                * from the top 500 dataset ..
                * They are not particularly well
                * fitted but should do OK for now!
                */
               NO = new Gaussian(parameter_hbond_NO_top_mean, parameter_hbond_NO_top_sd);
               NC = new Gaussian(parameter_hbond_NC_top_mean, parameter_hbond_NC_top_sd);
               HO = new Gaussian(parameter_hbond_HO_top_mean, parameter_hbond_HO_top_sd);
               HC = new Gaussian(parameter_hbond_HC_top_mean, parameter_hbond_HC_top_sd);
          }

          this->interaction_type = BB_HBOND;
     }

     //! Constructor - using only H and O atoms
     //!
     //! \param atomH Hydrogen atom
     //! \param atomO Oxygen atom
     //! \param init_from_native Whether to initialize gaussians from original structure
     AtomPairHBond(Atom *atomH, Atom *atomO, bool init_from_native=false) 
          : AtomPair<CHAIN_TYPE>(atomH, atomO) {
          
          this->init(init_from_native);
          this->get_energy();
     }

     //! Constructor - using N, H, C, and O atoms
     //!
     //! \param atomN Nitrogen atom
     //! \param atomH Hydrogen atom
     //! \param atomC Carbon atom
     //! \param atomO Oxygen atom
     //! \param interaction_type Type of hydrogen bond interaction
     //! \param init_from_native Whether to initialize gaussians from original structure
     AtomPairHBond(Atom *atomN, Atom *atomH, Atom *atomC, Atom *atomO, InteractionType interaction_type, bool init_from_native=false) 
          : AtomPair<CHAIN_TYPE>(atomH, atomO) {

          this->atomH = atomH;
          this->atomO = atomO;
          this->atomN = atomN;
          this->atomC = atomC;

          // Either Gaussians are initialized with the distances in the starting structure
          if (init_from_native) {
               double distanceNO = (this->atomN->position - this->atomO->position).norm();
               double distanceNC = (this->atomN->position - this->atomC->position).norm();
               double distanceHO = (this->atomH->position - this->atomO->position).norm();
               double distanceHC = (this->atomH->position - this->atomC->position).norm();

               this->NO = new Gaussian(distanceNO, parameter_hbond_native_sd);
               this->NC = new Gaussian(distanceNC, parameter_hbond_native_sd);
               this->HO = new Gaussian(distanceHO, parameter_hbond_native_sd);
               this->HC = new Gaussian(distanceHC, parameter_hbond_native_sd);

          // Or using values from the Top500 dataset
          } else {
               /* Values for the Gaussians derrived
                * from the top 500 dataset ..
                * They are not particularly well
                * fitted but should do OK for now!
                */
               NO = new Gaussian(parameter_hbond_NO_top_mean , parameter_hbond_NO_top_sd);
               NC = new Gaussian(parameter_hbond_NC_top_mean , parameter_hbond_NC_top_sd);
               HO = new Gaussian(parameter_hbond_HO_top_mean , parameter_hbond_HO_top_sd);
               HC = new Gaussian(parameter_hbond_HC_top_mean , parameter_hbond_HC_top_sd);
          }

          this->interaction_type = interaction_type;
          this->initialized = true;
          this->get_energy();
     }

     //! Constructor - using atom names, indices, and chain, and manually specifying ideal distances
     //!     
     //! \param atomNname Nitrogen atom name
     //! \param atomHname Hydrogen atom name
     //! \param atomOname Oxygen atom name
     //! \param atomCname Carbon atom name
     //! \param index1 Index of residue1 in chain
     //! \param index2 Index of residue2 in chain
     //! \param chain Molecule chain
     //! \param distanceNO Distance between N and O
     //! \param distanceNC Distance between N and C
     //! \param distanceHO Distance between H and O
     //! \param distanceHC Distance between H and C
     //! \param interaction_type Type of hydrogen bond interaction
     //! \param init_from_native Whether to initialize gaussians from original structure
     AtomPairHBond(CHAIN_TYPE *chain,
                   std::string atomNname, std::string atomHname, 
                   std::string atomOname, std::string atomCname, 
                   int index1, int index2, 
                   double distanceNO, double distanceNC, 
                   double distanceHO, double distanceHC, 
                   InteractionType interaction_type, bool init_from_native=false)
          : AtomPair<CHAIN_TYPE>() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          this->atomH = NULL;
          this->atomN = NULL;
          this->atomC = NULL;
          this->atomO = NULL;

          this->index1 = index1;
          this->index2 = index2;

          this->contribution = 0.;
          this->initialized = false;

          Residue *rf = &(*chain)[index1];
          Residue *rs = &(*chain)[index2];

          if (rf->has_atom(string_to_atom(atomHname)) && rs->has_atom(string_to_atom(atomOname)) ) {
               this->first = (*rf)[string_to_atom(atomHname)];
               this->second = (*rs)[string_to_atom(atomOname)];
          } else {
               std::cerr << "Missing atoms in the chains .. " << atomHname << " or " << atomOname << "\n";
               assert(false);
          }

          this->atomH = this->first;
          this->atomO = this->second;

          if (this->atomH->residue->has_atom(string_to_atom(atomNname)))
               this->atomN = (*(this->atomH->residue))[string_to_atom(atomNname)];

          if (this->atomO->residue->has_atom(string_to_atom(atomCname)))
               this->atomC = (*(this->atomO->residue))[string_to_atom(atomCname)];


          // Either Gaussians are initialized with default standard deviations 
          if (init_from_native) {
               this->NO = new Gaussian(distanceNO, parameter_hbond_native_sd);
               this->NC = new Gaussian(distanceNC, parameter_hbond_native_sd);
               this->HO = new Gaussian(distanceHO, parameter_hbond_native_sd);
               this->HC = new Gaussian(distanceHC, parameter_hbond_native_sd);

          // Or using values extracted from the Top500 dataset
          } else {
               /* Values for the Gaussians derrived
                * from the top 500 dataset ..
                * They are not particularly well
                * fitted but should do OK for now!
                */
               this->NO = new Gaussian(distanceNO, parameter_hbond_NO_top_sd);
               this->NC = new Gaussian(distanceNC, parameter_hbond_NC_top_sd);
               this->HO = new Gaussian(distanceHO, parameter_hbond_HO_top_sd);
               this->HC = new Gaussian(distanceHC, parameter_hbond_HC_top_sd);

          }

          this->interaction_type = interaction_type;
          if (this->atomH != NULL && this->atomO != NULL && this->atomN != NULL && this->atomC != NULL) {
               this->initialized = true;
               this->get_energy();
          } else {
               std::cerr << "Could not initialize Hbond properly .. \n" ;
          }
     }

     //! Constructor - from chain and indices
     //!
     //! \param chain Molecule chain
     //! \param first Index of residue1 in chain
     //! \param second Index of residue2 in chain
     //! \param init_from_native Whether to initialize gaussians from original structure
     AtomPairHBond(CHAIN_TYPE *chain, int first, int second, bool init_from_native=false)
          : AtomPair<CHAIN_TYPE> (chain, first, second, BB_HBOND) {

          this->init(init_from_native);
          this->get_energy();
     }

     //! Copy constructor.
     //! The gaussians are copied as pointers to save memory.
     //! \param other Source object from which copy is made
     AtomPairHBond(const AtomPairHBond &other)
          : AtomPair<CHAIN_TYPE>(other),
            atomN(other.atomN),
            atomH(other.atomH),
            atomC(other.atomC),
            atomO(other.atomO) {

          // Copy pointers
          this->NO = new Gaussian(*other.NO);
          this->NC = new Gaussian(*other.NC);
          this->HO = new Gaussian(*other.HO);
          this->HC = new Gaussian(*other.HC);
     }

     //! Assignment operator
     //!
     //! \param other Source object from which assignment is made
     const AtomPairHBond &operator=(const AtomPairHBond &other) {

          AtomPair<CHAIN_TYPE>::operator=(other);

          atomN = other.atomN;
          atomH = other.atomH;
          atomC = other.atomC;
          atomO = other.atomO;

          // Copy pointers
          this->NO = new Gaussian(*other.NO);
          this->NC = new Gaussian(*other.NC);
          this->HO = new Gaussian(*other.HO);
          this->HC = new Gaussian(*other.HC);
          
          return *this;
     }
     
     //! Destructor
     ~AtomPairHBond() {
          delete NO;
          delete NC;
          delete HO;
          delete HC;
     }

     //! DSSP energy factor
     //! (Kabsch W, Sander C, Biopolymers, 1983)
     //!
     //! \param move_info Object containing information about the last executed move
     double get_energy_dssp(MoveInfo *move_info = NULL) {
          if (not this->initialized )
               return std::numeric_limits<double>::infinity();

          double e = 0;

          double distanceNO = (this->atomN->position - this->atomO->position).norm();
          double distanceNC = (this->atomN->position - this->atomC->position).norm();
          double distanceHO = (this->atomH->position - this->atomO->position).norm();
          double distanceHC = (this->atomH->position - this->atomC->position).norm();

          if (distanceHO < dssp_min_HO_dssp)
               return 0;

          // Within H-bonding distance - continue
          e = dssp_qee * (1 / distanceNO + 1 / distanceHC - 1 / distanceHO - 1 / distanceNC);

          this->contribution = e;
          return e;
     }

     //! Retrieve energy
     //!
     //! \param move_info Object containing information about the last executed move
     double get_energy(MoveInfo *move_info = NULL) {
          if (not this->initialized )
               return std::numeric_limits<double>::infinity();
          this->contribution_cache = 0;
          if (move_info != NULL  && not this->in_range(move_info) && this->contribution != 0.) {
               return this->contribution;
          }
          double e = 0;

          //! Calculate energy from Gaussians
          e+= NO->get_energy((this->atomN->position - this->atomO->position).norm());
          e+= NC->get_energy((this->atomN->position - this->atomC->position).norm());
          e+= HO->get_energy((this->atomH->position - this->atomO->position).norm());
          e+= HC->get_energy((this->atomH->position - this->atomC->position).norm());

          this->contribution_cache = this->contribution;
          this->contribution = e;
          return e;
     }

};


//! Atom pair - Gaussian interaction
//!
//! More specific implementation of the atom pair class, 
//! for a Gaussian controlled distance between two atoms.
template<typename CHAIN_TYPE, typename ATOM_TYPE1=Atom, typename ATOM_TYPE2=Atom>
class AtomPairGaussian: public AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2> {
public:

     //! Gaussian
     Gaussian *gaussian;

     //! Initializer
     //! \param mean Mean value
     //! \param sd Standard deviation
     void init(double mean = -1., double sd = -1.) {

          if (mean < 0 && this->first != NULL)
               mean = get_distance(this->first, this->second);
          if (sd <= 0)
               sd = (mean / 6.);

          // this->mean = mean;
          // this->sd = sd;

          // this->gauss_const = 0.5 * std::log(2 * M_PI * this->sd * this->sd);

          gaussian = new Gaussian(mean, sd);

     }

     //! Constructor
     AtomPairGaussian()
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2>() {}

     //! Constructor - from atom pointers
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param interaction_type specify the interaction type.
     AtomPairGaussian(ATOM_TYPE1 *first, ATOM_TYPE2 *second, 
                      InteractionType interaction_type = DISTANCE_GAUSSIAN)
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2> (first, second, interaction_type) {
          this->init();
          this->get_energy();
     }

     //! Constructor - from atom pointers - mean and stddev provided
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     //! \param interaction_type specify the interaction type.
     AtomPairGaussian(ATOM_TYPE1 *first, ATOM_TYPE2 *second, 
                      double mean, double sd, 
                      InteractionType interaction_type = DISTANCE_GAUSSIAN)
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2> (first, second, interaction_type) {
          this->init(mean, sd);
          this->get_energy();
     }

     //! Constructor - from chain and residue indices
     //!
     //! \param chain Molecule chain.
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     //! \param interaction_type specify the interaction type.
     AtomPairGaussian(CHAIN_TYPE *chain, 
                      int first, int second, 
                      double mean, double sd, 
                      InteractionType interaction_type = DISTANCE_GAUSSIAN)
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2> (chain, first, second, interaction_type) {
          this->init(mean, sd);
          this->get_energy();
     }

     //! Constructor - from chain and residue indices
     //!
     //! \param chain Molecule chain.
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param interaction_type specify the interaction type.
     AtomPairGaussian(CHAIN_TYPE *chain, 
                      int first, int second, 
                      InteractionType interaction_type = DISTANCE_GAUSSIAN) 
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2> (chain, first, second, interaction_type) {}

     //! Constructor - using atom names, indices, and chain, and manually specifying ideal distances
     //!     
     //! \param chain Molecule chain.
     //! \param first_name Name of first atom
     //! \param second_name Name of second atom
     //! \param index1 Index of residue1 in chain
     //! \param index2 Index of residue2 in chain
     //! \param chain Molecule chain
     //! \param mean Mean distance for gaussian
     //! \param sd Standard deviation for gaussian
     //! \param interaction_type Type of interaction
     AtomPairGaussian(CHAIN_TYPE *chain,
                      std::string first_name, std::string second_name, 
                      int index1, int index2, 
                      double mean, double sd, 
                      InteractionType interaction_type=DISTANCE_GAUSSIAN)
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2>(interaction_type) {

          // Import protein definitions (such as residue names)
          using namespace definitions;
          
          this->index1 = index1;
          this->index2 = index2;
          
          Residue *rf = &(*chain)[index1];
          Residue *rs = &(*chain)[index2];

          if (rf->has_atom(string_to_atom(first_name)) && 
              rs->has_atom(string_to_atom(second_name))) {
               this->first = (*rf)[string_to_atom(first_name)];
               this->second = (*rs)[string_to_atom(second_name)];
          } else {
               std::cerr << "Cannot find atom " << string_to_atom(second_name) << 
                    " (" << index2 << ") or " << string_to_atom(first_name) << " (" << index1 << ") \n";
               assert(false);
          }
          
          init(mean, sd);
          this->interaction_type = interaction_type;
          this->initialized = true;
          this->get_energy();
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made
     AtomPairGaussian(const AtomPairGaussian &other)
          : AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2>(other) {

          // Copy pointer
          this->gaussian = new Gaussian(*other.gaussian);
     }

     //! Assignment operator
     //!
     //! \param other Source object from which assignment is made
     const AtomPairGaussian &operator=(const AtomPairGaussian &other) {

          AtomPair<CHAIN_TYPE,ATOM_TYPE1,ATOM_TYPE2>::operator=(other);

          // Copy pointer
          this->gaussian = new Gaussian(*other.gaussian);
          
          return *this;
     }

     //! Destructor
     ~AtomPairGaussian() {
          delete gaussian;
     }

     //! Retrieve energy
     //!
     //! \param move_info Object containing information about the last executed move
     double get_energy(MoveInfo *move_info = NULL) {
          if (not this->initialized )
               return std::numeric_limits<double>::infinity();

          this->contribution_cache = 0;
          if (move_info != NULL  && not this->in_range(move_info) && this->contribution != 0.) {
               return this->contribution;
          }

          double distance = get_distance(this->first, this->second);
          double energy = gaussian->get_energy(distance);

          this->contribution_cache = this->contribution;
          this->contribution = energy;
          return this->contribution;
     }
};



//! Atom pair - Log-normal interaction
//!
//! More specific implementation of the atom pair class, 
//! for a log-normalcontrolled distance between two atoms.
template<typename CHAIN_TYPE>
class AtomPairLogNormal: public AtomPair<CHAIN_TYPE> {
public:
     //! Mean value
     double mean;

     //! Standard deviation     
     double sd;

     //! Normalization constant          
     double gauss_const;

     //! Initializer
     //! \param mean Mean value
     //! \param sd Standard deviation
     void init(double mean = -1., double sd = -1.) {
          if (mean < 0)
               mean = get_distance(this->first, this->second);
          if (sd <= 0)
               sd = (mean / 6.);

          this->mean = mean;
          this->sd = sd;

          this->gauss_const = sqrt(2 * M_PI * this->sd * this->sd);

          this->interaction_type = LOG_NORMAL;
     }

     //! Constructor - from atom pointers
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom     
     AtomPairLogNormal(Atom *first, Atom *second) :
          AtomPair<CHAIN_TYPE> (first, second, LOG_NORMAL) {
          this->init();
          this->get_energy();
     }

     //! Constructor - from atom pointers - mean and stddev provided
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     AtomPairLogNormal(Atom *first, Atom *second, 
                       double mean, double sd) :
          AtomPair<CHAIN_TYPE> (first, second, LOG_NORMAL) {
          this->init(mean, sd);
          this->get_energy();
     }

     //! Constructor - from chain and residue indices
     //!
     //! \param chain Molecule chain.
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     AtomPairLogNormal(CHAIN_TYPE *chain, 
                       int first, int second, 
                       double mean, double sd) :
          AtomPair<CHAIN_TYPE> (chain, first, second, LOG_NORMAL) {
          this->init(mean, sd);
          this->get_energy();
     }

     //! Retrieve energy
     //!
     //! \param move_info Object containing information about the last executed move
     double get_energy(MoveInfo *move_info = NULL) {
          if (not this->initialized )
               return std::numeric_limits<double>::infinity();

          this->contribution_cache = 0;
          if (move_info != NULL  && not this->in_range(move_info) && this->contribution != 0.) {
               return this->contribution;
          }
          double distance = get_distance(this->first, this->second);
          double dens = -(std::log(distance) - this->mean) * (std::log(distance) - this->mean) / 
               (2 * this->sd * this->sd) - (std::log(distance * this->gauss_const));
          this->contribution_cache = this->contribution;
          this->contribution = -dens;
          return this->contribution;
     }
};


//! Atom pair - Disulfide interactions
//!
//! More specific implementation of the Gaussian atom pair class, 
//! for modeling disulfide interactions.
template<typename CHAIN_TYPE>
class AtomPairSSBond: public AtomPairGaussian<CHAIN_TYPE> {

     //! Default standard deviation
     static double default_sd(double sd=-1) {
          if (sd <= 0)
               sd = 0.2;
          return sd;
     }

     //! Default standard deviation
     static double default_mean(double mean=-1) {
          if (mean < 0)
               mean = 2.05;
          return mean;
     }

public:

     //! Constructor - from atom pointers - mean and stddev provided
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     AtomPairSSBond(Atom *first, Atom *second, 
                    double mean=-1, double sd=-1) :
          AtomPairGaussian<CHAIN_TYPE> (first, second, default_mean(mean), default_sd(sd), SS_BOND) {
     }

     //! Constructor - from chain and residue indices
     //!
     //! \param chain Molecule chain.
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     AtomPairSSBond(CHAIN_TYPE *chain, 
                    int first, int second, 
                    double mean, double sd) :
          AtomPairGaussian<CHAIN_TYPE> (chain, first, second, default_mean(mean), default_sd(sd), SS_BOND) {
     }

};


//! Atom pair - Between atom and a fixed point (implemented as a Dummy atom)
//!
//! More specific implementation of the Gaussian atom pair class, 
//! for modeling interactions between an atom and a fixed point
template<typename CHAIN_TYPE>
class AtomPairGaussianFixedPoint: public AtomPairGaussian<CHAIN_TYPE,DummyAtom,Atom> {

     //! Default standard deviation
     static double default_sd(double sd=-1) {
          if (sd <= 0)
               sd = 0.2;
          return sd;
     }

     //! Lookup atom pointer in chain
     static Atom* get_atom_pointer(const CHAIN_TYPE *chain, const int index, const std::string &atom_name) {
          Residue *res = &(*chain)[index];
          if (res->has_atom(definitions::string_to_atom(atom_name))) {
               return (*res)[definitions::string_to_atom(atom_name)];
          } else {
               // Can't find atoms .. bailing out
               std::cerr << "Cannot find atom with name " << atom_name << " in residue " << *(res) << "\n";
               assert(false);
          }          
     }

public:

     //! Constructor - from atom pointers
     //!
     //! \param first atom pointer to the first atom
     //! \param second atom pointer to the second atom
     //! \param mean Mean value
     //! \param sd Standard deviation
     AtomPairGaussianFixedPoint(DummyAtom *first, Atom *second, double mean=-1, double sd=-1) 
          : AtomPairGaussian<CHAIN_TYPE,DummyAtom,Atom> (first, second, mean, default_sd(sd), DISTANCE_GAUSSIAN_FIXED_POINT) {
     }


     //! Constructor - using atom names, indices, and chain, and manually specifying ideal distances
     //!     
     //! \param first Dummy atom
     //! \param atom_name2 Name of second atom
     //! \param index2 Index of residue2 in chain
     //! \param chain Molecule chain
     //! \param mean Mean distance for gaussian
     //! \param sd Standard deviation for gaussian
     AtomPairGaussianFixedPoint(DummyAtom *first, std::string atom_name2, int index2, CHAIN_TYPE *chain, double mean, double sd)
          : AtomPairGaussian<CHAIN_TYPE,DummyAtom,Atom> (first, 
                                                         get_atom_pointer(chain, index2, atom_name2),
                                                         mean, sd,
                                                         DISTANCE_GAUSSIAN_FIXED_POINT) {
     }

};

} // End of constrain_distances namespace




//! Energy term for constraining distances between atom pairs
template<typename CHAIN_TYPE>
class TermConstrainDistances: public EnergyTermCommon<TermConstrainDistances<CHAIN_TYPE>, CHAIN_TYPE> {

     // For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermConstrainDistances<CHAIN_TYPE> , CHAIN_TYPE> EnergyTermCommon;

protected:

     //! Hydrogen bond network
     std::vector<constrain_distances::AtomPairHBond<CHAIN_TYPE> > hb_network;

     //! Gaussian interaction network
     std::vector<constrain_distances::AtomPairGaussian<CHAIN_TYPE> > g_network;

     //! Disulfide bond network
     std::vector<constrain_distances::AtomPairSSBond<CHAIN_TYPE> > ss_network;

     //! Network between atoms and fixed points
     std::vector<constrain_distances::AtomPairGaussianFixedPoint<CHAIN_TYPE> > gfix_network;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! Path to network input file
          std::string network_filename;

          //! Path to PDB input file
          std::string pdb_file;

          //! Whether to include the hydrogen bond contacts
          bool include_bb_hbond;

          //! Whether to include the sidechain-sidechain hydrogen bond contacts
          bool include_sc_hbond;

          //! Whether to include the backbone-sidechain hydrogen bond contacts
          bool include_bb_sc_hbond;

          //! Whether to include the Calpha-Calpha contacts
          bool include_ca_contacts;

          //! Whether to include Atom-fixedpoint contacts
          bool include_fixed_point_contacts;

          //! Whether to include disulfide contacts
          bool include_ss_bond;

          //! Whether to initialize ideal hbond distances from PDB file (rather than averages from the Top500 database)
          bool init_hbond_from_native;

          //! Whether to prune/optimize the hbond network
          bool prune_network;

          //! Cutoff below which backbone-backbone hydrogen bonds are considered as dehydrated aka weak/broken
          int dehydron_bb_cutoff;

          //! Cutoff below which backbone-sidechain hydrogen bonds are considered as dehydrated aka weak/broken
          int dehydron_bb_sc_cutoff;

          //! Cutoff below which sidechain-sidechain hydrogen bonds are considered as dehydrated aka weak/broken
          int dehydron_sc_cutoff;

          //! Cutoff distance in Angstrom to be considered a Calpha contact
          double ca_distance;

          //! Cutoff distance in Angstrom to be considered a disulfide contact
          double ss_distance;

          //! How many residues to skip along the chain before considering a Calpha contact
          int ca_skip;

          //! How many residues to skip along the chain before considering a backbone-backbone hydrogen bond contact
          int bb_hbond_skip;

          //! How many residues to skip along the chain before considering a sidechain-sidechain hydrogen bond contact
          int sc_hbond_skip;

          //! How many residues to skip along the chain before considering a backbone-sidechain hydrogen bond contact
          int bb_sc_hbond_skip;

          //! How many residues to skip along the chain before considering a disulfide contact
          int ss_skip;

          //! Whether to generate a Python script that will visualize the network in PyMOL
          bool generate_pymol;

          //! Whether to cache interations
          bool use_caching;

          //! Whether to print out additional information
          bool verbose;

          //! Constructor
          Settings(std::string network_filename = "",
                   std::string pdb_file = "", 
                   bool include_bb_hbond = true, 
                   bool include_sc_hbond = true, 
                   bool include_bb_sc_hbond = true,
                   bool include_ca_contacts = true, 
                   bool include_fixed_point_contacts = false,  
                   bool include_ss_bond = true, 
                   bool init_hbond_from_native = false, 
                   bool prune_network=true, 
                   int dehydron_bb_cutoff = 14,
                   int dehydron_bb_sc_cutoff = 9, 
                   int dehydron_sc_cutoff = 7,
                   double ca_distance = 6., 
                   double ss_distance = 3., 
                   int ca_skip = 5, 
                   int bb_hbond_skip = 1, 
                   int sc_hbond_skip = 1, 
                   int bb_sc_hbond_skip = 1, 
                   int ss_skip = 1, 
                   bool generate_pymol=false, 
                   bool use_caching=true, 
                   bool verbose=false) 
               : network_filename(network_filename), 
                 pdb_file(pdb_file), 
                 include_bb_hbond(include_bb_hbond), 
                 include_sc_hbond(include_sc_hbond), 
                 include_bb_sc_hbond(include_sc_hbond), 
                 include_ca_contacts(include_ca_contacts),
                 include_fixed_point_contacts(include_fixed_point_contacts), 
                 include_ss_bond(include_ss_bond), 
                 init_hbond_from_native(init_hbond_from_native), 
                 prune_network(prune_network), 
                 dehydron_bb_cutoff(dehydron_bb_cutoff),
                 dehydron_bb_sc_cutoff(dehydron_bb_sc_cutoff), 
                 dehydron_sc_cutoff(dehydron_sc_cutoff),
                 ca_distance(ca_distance), 
                 ss_distance(ss_distance), 
                 ca_skip(ca_skip), 
                 bb_hbond_skip(bb_hbond_skip), 
                 sc_hbond_skip(sc_hbond_skip), 
                 bb_sc_hbond_skip(bb_sc_hbond_skip), 
                 ss_skip(ss_skip), 
                 generate_pymol(generate_pymol), 
                 use_caching(use_caching), 
                 verbose(verbose) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "network-filename:" << settings.network_filename << "\n";
               o << "pdb-file:" << settings.pdb_file << "\n";
               o << "include-bb-hbond:" << settings.include_bb_hbond << "\n";
               o << "include-sc-hbond:" << settings.include_sc_hbond << "\n";
               o << "include-bb-sc-hbond:" << settings.include_bb_sc_hbond << "\n";
               o << "include-ca-contacts:" << settings.include_ca_contacts << "\n";
               o << "include-fixed-point-contacts:" << settings.include_fixed_point_contacts << "\n";
               o << "include-ss-bond:" << settings.include_ss_bond << "\n";
               o << "init-hbond-from-native:" << settings.init_hbond_from_native << "\n";
               o << "prune-network:" << settings.prune_network << "\n";
               o << "dehydron-bb-bond:" << settings.dehydron_bb_cutoff << "\n";
               o << "dehydron-bb-sc-bond:" << settings.dehydron_bb_sc_cutoff << "\n";
               o << "dehydron-sc-bond:" << settings.dehydron_sc_cutoff << "\n";
               o << "ca-distance:" << settings.ca_distance << "\n";
               o << "ss-distance:" << settings.ss_distance << "\n";
               o << "ca-skip:" << settings.ca_skip << "\n";
               o << "bb-hbond-skip:" << settings.bb_hbond_skip << "\n";
               o << "sc-hbond-skip:" << settings.sc_hbond_skip << "\n";
               o << "bb-sc-hbond-skip:" << settings.bb_sc_hbond_skip << "\n";
               o << "ss-skip:" << settings.ss_skip << "\n";
               o << "generate-pymol:" << settings.generate_pymol << "\n";
               o << "use-caching:" << settings.use_caching << "\n";
               o << "verbose:" << settings.verbose << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings> (settings);
               return o;
          }
     } settings;    //!< Local settings object 


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object     
     //! \param random_number_engine Object from which random number generators can be created.
     TermConstrainDistances(CHAIN_TYPE *chain, 
                            const Settings &settings = Settings(),
                            RandomNumberEngine *random_number_engine = &random_global) 
          : EnergyTermCommon(chain, "constrain-distances", settings, random_number_engine), 
            settings(settings) {

          if (settings.use_caching) {
               this->clear_cache();
          }

          if (settings.verbose)
               std::cout << "# initializing constrain distances energy term ..  \n";

          // Initialize network
          if (settings.network_filename == "") {
               this->init_network();
          } else {
               this->init_network_from_file(settings.network_filename);
          }

          // Display network
          if (this->settings.verbose) {
               std::cout << "######################################################################################################\n";
               this->print_network();
               std::cout << "######################################################################################################\n";
          }

          // Display network in Pymol format
          if (this->settings.generate_pymol) {
               this->generate_pymol_network();
          }
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermConstrainDistances(const TermConstrainDistances &other, 
                            RandomNumberEngine *random_number_engine,
                            int thread_index, CHAIN_TYPE *chain = NULL) 
          : EnergyTermCommon(other, random_number_engine, thread_index, chain), 
            settings(other.settings) {

          for (unsigned int i = 0; i < other.hb_network.size(); i++) {
              this->hb_network.push_back(other.hb_network[i]);
          }
          for (unsigned int i = 0; i < other.g_network.size(); i++) {
               this->g_network.push_back(other.g_network[i]);
          }
          for (unsigned int i = 0; i < other.ss_network.size(); i++) {
               this->ss_network.push_back(other.ss_network[i]);
          }
          for (unsigned int i = 0; i < other.gfix_network.size(); i++) {
               this->gfix_network.push_back(other.gfix_network[i]);
          }
          this->chain = other.chain;
     }

     //! Check whether atom is in a carbonaceous group.
     //! This definition of carbonaceous groups follows the definition of
     //! Fernandez et al. (Biophysical Journal, Vol 85. Sept 2003).
     //! \param atom Atom pointer
     bool is_carbonaceous(Atom* atom) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // well it should be a carbon
          if (not is_atom_XC(atom->atom_type) )// and not is_atom_XH(atom->atom_type) )
               return false;
          // for some reason backbones are not included
          if (atom->atom_type == C || atom->atom_type == CA )
               return false;
          // and finally also exclude some sidechain carbons,
          // that have no hydrogens attached.
          if (atom->residue->residue_type == ARG && atom->atom_type == CZ)
               return false;
          if (atom->residue->residue_type == ASN && atom->atom_type == CG)
               return false;
          if (atom->residue->residue_type == ASP && atom->atom_type == CG)
               return false;
          if (atom->residue->residue_type == GLU && atom->atom_type == CD)
               return false;
          if (atom->residue->residue_type == GLN && atom->atom_type == CD)
               return false;

          return true;
     }

     //! Update the hydrogen bond network with a measure of the degree of desolvation
     void calc_desolvation() {

          // Import local definitions
          using namespace constrain_distances;

          // Create a list of all CHn groups (where n=1,2,3)
          std::vector<Atom*> carbonaceous;
          for (ResidueIterator<CHAIN_TYPE> res1(*(this->chain)); !res1.end(); ++res1) {
               for (AtomIterator<CHAIN_TYPE,definitions::ALL> at1(*res1); !at1.end(); ++at1) {
                    if (is_carbonaceous(&(*at1)))
                         carbonaceous.push_back(&(*at1));
               }
          }

          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               // for every hbond in the network we calculate
               // the number of carbogenious groups in the
               // desolvation shell .. (Fernandez, Sheraga (2003))
               Residue *rf = &(*this->chain)[this->hb_network[i].atomN->residue->index];
               Residue *rs = &(*this->chain)[this->hb_network[i].atomC->residue->index];
               Atom* ca_f = (*rf)[this->hb_network[i].atomN->atom_type];
               Atom* ca_s = (*rs)[this->hb_network[i].atomC->atom_type];

               int desolv = 0;

               for (int j=0; j<(int)carbonaceous.size(); j++) {

                    if (get_distance(carbonaceous[j], ca_f) < 6.5 )  {
                         desolv++;
                    } else if (get_distance(carbonaceous[j], ca_s) < 6.5 ) {
                         desolv++;
                    }
               }
               this->hb_network[i].desolvation_factor = desolv;
          }
     }

     //! Removes duplicate and weak bonds from the hydrogen bond network
     void prune_network() {

          // Import local definitions
          using namespace constrain_distances;

          if (!this->settings.prune_network)
               return;

          if (this->settings.verbose)
               std::cout << "\n# Removing duplicate hbonds : \n";

          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               int res11 = this->hb_network[i].atomN->residue->index;
               int res12 = this->hb_network[i].atomC->residue->index;

               for (unsigned int j = i+1; j < this->hb_network.size(); j++) {
                    int res21 = this->hb_network[j].atomN->residue->index;
                    int res22 = this->hb_network[j].atomC->residue->index;

                    // lets check more multiple HBonds to the same atom
                    if ( (res11 == res21 && this->hb_network[i].atomN->atom_type == this->hb_network[j].atomN->atom_type ) ||
                         (res12 == res22 && this->hb_network[i].atomC->atom_type == this->hb_network[j].atomC->atom_type ) ) {
                         // and which is the "stronger" Bond
                         if (this->hb_network[i].get_energy_dssp() <  this->hb_network[j].get_energy_dssp()) {
                              if (this->settings.verbose) {
                                   std::cout << "# Duplicate hbond will be removed (j): ";
                                   std::cout << this->hb_network[j].atomN->residue->index << " (" << this->hb_network[j].atomN->atom_type << " in " << this->hb_network[j].atomN->residue->residue_type << ")\t";
                                   std::cout << this->hb_network[j].atomC->residue->index << " (" << this->hb_network[j].atomC->atom_type << " in " << this->hb_network[j].atomC->residue->residue_type << ")\t";
                                   std::cout << "[ " << this->hb_network[j].get_energy_dssp() << "\t> " << this->hb_network[i].get_energy_dssp() << " ] \t";
                                   std::cout << this->hb_network[i].atomN->residue->index << " (" << this->hb_network[i].atomN->atom_type << " in " << this->hb_network[i].atomN->residue->residue_type << ")\t";
                                   std::cout << this->hb_network[i].atomC->residue->index << " (" << this->hb_network[i].atomC->atom_type << " in " << this->hb_network[i].atomC->residue->residue_type << ")\n";
                              }
                              this->hb_network.erase(this->hb_network.begin()+j);
                              j--;
                         } else {
                              if (this->settings.verbose) {
                                   std::cout << "# Duplicate hbond will be removed (i): ";
                                   std::cout << this->hb_network[i].atomN->residue->index << " (" << this->hb_network[i].atomN->atom_type << " in " << this->hb_network[i].atomN->residue->residue_type << ")\t";
                                   std::cout << this->hb_network[i].atomC->residue->index << " (" << this->hb_network[i].atomC->atom_type << " in " << this->hb_network[i].atomC->residue->residue_type << ")\t";
                                   std::cout << "[ " << this->hb_network[i].get_energy_dssp() << "\t> " << this->hb_network[j].get_energy_dssp() << " ]\t";
                                   std::cout << this->hb_network[j].atomN->residue->index << " (" << this->hb_network[j].atomN->atom_type << " in " << this->hb_network[j].atomN->residue->residue_type << ")\t";
                                   std::cout << this->hb_network[j].atomC->residue->index << " (" << this->hb_network[j].atomC->atom_type << " in " << this->hb_network[j].atomC->residue->residue_type << ")\n";
                              }
                              this->hb_network.erase(this->hb_network.begin()+i);
                              i--;
                              break;
                         }
                    }
               }
          }
          if (this->settings.verbose)
               std::cout << "\n\n# Removing weak hbonds : \n";
          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               if (this->hb_network[i].interaction_type == BB_HBOND && this->hb_network[i].desolvation_factor <= settings.dehydron_bb_cutoff) {
                    if (this->settings.verbose) {
                         std::cout << "# Weak BB-BB hbond will be removed : ";
                         std::cout << this->hb_network[i].atomN->residue->index << " (" << this->hb_network[i].atomN->residue->residue_type << ")\t";
                         std::cout << this->hb_network[i].atomC->residue->index << " (" << this->hb_network[i].atomC->residue->residue_type << ")\t";
                         std::cout << "[ " << this->hb_network[i].desolvation_factor << "<=" << settings.dehydron_bb_cutoff << " ]\n";
                    }
                    this->hb_network.erase(this->hb_network.begin()+i);
                    i--;
                    continue;
               }
               if (this->hb_network[i].interaction_type == BB_SC_HBOND && this->hb_network[i].desolvation_factor <= settings.dehydron_bb_sc_cutoff) {
                    if (this->settings.verbose) {
                         std::cout << "# Weak BB-SC hbond will be removed : ";
                         std::cout << this->hb_network[i].atomN->residue->index << " (" << this->hb_network[i].atomN->residue->residue_type << ")\t";
                         std::cout << this->hb_network[i].atomC->residue->index << " (" << this->hb_network[i].atomC->residue->residue_type << ")\t";
                         std::cout << "[ " << this->hb_network[i].desolvation_factor << "<=" << settings.dehydron_bb_sc_cutoff << " ]\n";
                    }
                    this->hb_network.erase(this->hb_network.begin()+i);
                    i--;
                    continue;
               }
               if (this->hb_network[i].interaction_type == SC_HBOND && this->hb_network[i].desolvation_factor <= settings.dehydron_sc_cutoff) {
                    if (this->settings.verbose) {
                         std::cout << "# Weak SC-SC hbond will be removed : ";
                         std::cout << this->hb_network[i].atomN->residue->index << " (" << this->hb_network[i].atomN->residue->residue_type << ")\t";
                         std::cout << this->hb_network[i].atomC->residue->index << " (" << this->hb_network[i].atomC->residue->residue_type << ")\t";
                         std::cout << "[ " << this->hb_network[i].desolvation_factor << "<=" << settings.dehydron_sc_cutoff << " ]\n";
                    }
                    this->hb_network.erase(this->hb_network.begin()+i);
                    i--;
                    continue;
               }

          }

          if (this->settings.verbose)
               std::cout << "\n\n# Removing redundant link from the ligands : \n";
          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               DummyAtom *now_dummy1 = this->gfix_network[i].first;
               int res1 = this->gfix_network[i].second->residue->index;
               double d1 = get_distance(now_dummy1, this->gfix_network[i].second);

               for (unsigned int j = i + 1; j < this->gfix_network.size(); j++) {
                    DummyAtom *now_dummy2 = this->gfix_network[j].first;
                    int res2 = this->gfix_network[j].second->residue->index;

                    if (res1 != res2)
                         continue;

                    // if (now_dummy1->position[0] != now_dummy2->position[0] || now_dummy1->position[1] != now_dummy2->position[1] || now_dummy1->position[2]
                    //          != now_dummy2->position[2])
                    if (now_dummy1->id != now_dummy2->id || now_dummy1->id == 0 || now_dummy2->id == 0)
                         continue;

                    // and which is the "stronger" Bond .. aka the shorter
                    double d2 = get_distance(now_dummy2, this->gfix_network[j].second);
                    //
                    if (d1 < d2) {
                         if (this->settings.verbose) {
                              std::cout << "# Duplicate ligand contact will be removed (j): ";
                              std::cout << this->gfix_network[j].first->atom_name << " (" << this->gfix_network[j].first->position << ")\t";
                              std::cout << this->gfix_network[j].second->residue->index << " (" << this->gfix_network[j].second->atom_type << " in "
                                        << this->gfix_network[j].second->residue->residue_type << ")\t";
                              std::cout << "[ " << d2 << "\t> " << d1 << " ] \t";
                              std::cout << this->gfix_network[i].first->atom_name << " (" << this->gfix_network[i].first->position << ")\t";
                              std::cout << this->gfix_network[i].second->residue->index << " (" << this->gfix_network[i].second->atom_type << " in "
                                        << this->gfix_network[i].second->residue->residue_type << ")\n";
                         }
                         this->gfix_network.erase(this->gfix_network.begin() + j);
                         j--;
                    } else {
                         if (this->settings.verbose) {
                              std::cout << "# Duplicate ligand contact will be removed (i): ";
                              std::cout << this->gfix_network[i].first->atom_name << " (" << this->gfix_network[i].first->position << ")\t";
                              std::cout << this->gfix_network[i].second->residue->index << " (" << this->gfix_network[i].second->atom_type << " in "
                                        << this->gfix_network[i].second->residue->residue_type << ")\t";
                              std::cout << "[ " << d2 << "\t> " << d1 << " ] \t";
                              std::cout << this->gfix_network[j].first->atom_name << " (" << this->gfix_network[j].first->position << ")\t";
                              std::cout << this->gfix_network[j].second->residue->index << " (" << this->gfix_network[j].second->atom_type << " in "
                                        << this->gfix_network[j].second->residue->residue_type << ")\n";
                         }
                         this->gfix_network.erase(this->gfix_network.begin() + i);
                         i--;
                         break;
                    }

               }
          }

     }

     //! Read in network information from file
     //!
     //! \param filename Path to network file
     void init_network_from_file(std::string filename) {

          // Import local definitions
          using namespace constrain_distances;

          // Clear networks
          this->hb_network.clear();
          this->g_network.clear();
          this->ss_network.clear();
          this->gfix_network.clear();

          std::ifstream ifs(filename.c_str());
          if (ifs.fail()) {
               std::cout << "Error reading " << filename << "\n";
               exit(1);
          }
          std::string temp;

          bool found_start = false;

          while (getline(ifs, temp)) {

               // filter out comments lines
               if (!found_start && temp.substr(0, 38) != "# [ Distance Constrain Network Start ]") {
                    continue;
               }

               if (!found_start) {
                    // first line .. we just need
                    // to remember that we found a
                    // the start token
                    found_start = true;
                    continue;
               }
               // still want to filter comments
               if (temp.substr(0, 1) == "#") {
                    continue;
               }

               int index1 = 0;
               int index2 = 0;
               InteractionType interaction_type = NO_INTERACTION;
               double e = 0.;

               double mean = 0, sd = 0;
               double x = 0, y = 0, z = 0;

               char atom1[10];
               char atom2[10];
               char atom3[10];
               char atom4[10];

               double NOmean =0.;
               double NCmean =0.;
               double HOmean =0.;
               double HCmean =0.;

               // Check for hydrogen bond entry: four atoms, including the
               // according mean values for the Gaussians.
               int interaction_type_int = 0;
               int r_value = sscanf(temp.c_str(), "%d %d %d %le %s %s %s %s %le %le %le %le", &interaction_type_int, &index1, &index2, &e, atom1, atom2, atom3, atom4, &NOmean, &NCmean, &HOmean, &HCmean);
               interaction_type = InteractionType(interaction_type_int);
               if (r_value < 12) {

                    // Check for gaussian-fixedpoint entry
                    interaction_type = DISTANCE_GAUSSIAN_FIXED_POINT;
                    r_value = sscanf(temp.c_str(), "5 %d %d %le %s %s %le %le %le %le %le", &index1, &index2, &e, atom1, atom2, &x, &y, &z, &mean, &sd);
                    if (r_value < 10) {

                         // Check for entry specifying arbitrary gaussian springs between any
                         // two atoms - including the mean and sd of the gaussian
                         r_value = sscanf(temp.c_str(), "%d %d %d %le %s %s %le %le", &interaction_type_int, &index1, &index2, &e, atom1, atom2, &mean, &sd);
                         interaction_type = InteractionType(interaction_type_int);
                         if (r_value < 8) {

                              //! Check for general format: interaction type, indices, mean, stddev, energy contribution
                              r_value = sscanf(temp.c_str(), "%d %d %d %le %le %le", &interaction_type_int, &index1, &index2, &e, &mean, &sd);
                              interaction_type = InteractionType(interaction_type_int);

                              if (r_value < 6) {
                                   mean = -1;
                                   sd = -1;
                                   //! Check for general format: interaction type, indices, energy contribution
                                   int r_value = sscanf(temp.c_str(), "%d %d %d %le", &interaction_type_int, &index1, &index2, &e);
                                   interaction_type = InteractionType(interaction_type_int);

                                   if (r_value < 3) {
                                        //! Check for general format: interaction type, indices
                                        int r_value = sscanf(temp.c_str(), "%d %d %d ", &interaction_type_int, &index1, &index2);
                                        interaction_type = InteractionType(interaction_type_int);
                                        if (r_value < 2)
                                             continue;
                                   }
                              }
                         }
                    }
               }

               // Add hydrogenbond interaction if applicable
               if ((settings.include_bb_hbond && NOmean > 0 && NCmean > 0 && HOmean > 0 && HCmean > 0 && interaction_type == BB_HBOND ) || 
                   (settings.include_sc_hbond && NOmean > 0 && NCmean > 0 && HOmean > 0 && HCmean > 0 && interaction_type == SC_HBOND ) ||
                   (settings.include_bb_sc_hbond && NOmean > 0 && NCmean > 0 && HOmean > 0 && HCmean > 0 && interaction_type == BB_SC_HBOND)) {
                    this->hb_network.push_back(*(new AtomPairHBond<CHAIN_TYPE> (this->chain, 
                                                                                std::string(atom1), std::string(atom2), 
                                                                                std::string(atom4), std::string(atom3),
                                                                                index1, index2, 
                                                                                NOmean, NCmean, HOmean, HCmean, 
                                                                                interaction_type, this->settings.init_hbond_from_native)));

               // Add Gaussian interaction if applicable
               } else if (settings.include_ca_contacts && mean != 0. && sd != 0. && interaction_type == DISTANCE_GAUSSIAN) {
                    this->g_network.push_back(*(new AtomPairGaussian<CHAIN_TYPE> (this->chain, std::string(atom1), std::string(atom2), index1, index2, mean, sd, interaction_type)));

               // Add Gaussian-fixedpoint interaction if applicable
               } else if (settings.include_fixed_point_contacts && mean != 0. && sd != 0.  && interaction_type == DISTANCE_GAUSSIAN_FIXED_POINT) {
                    this->gfix_network.push_back(*(new AtomPairGaussianFixedPoint<CHAIN_TYPE> (
                                                        new DummyAtom(std::string(atom1), Vector_3D(x,y,z), 
                                                                      index1, 0) ,std::string(atom2),  index2, this->chain, mean, sd ) ) );

               // Add hydrogen bond interaction if applicable - only residue indices specified
               } else if (settings.include_bb_hbond && interaction_type == BB_HBOND) {
                    this->hb_network.push_back(*(new AtomPairHBond<CHAIN_TYPE> (this->chain, index1, index2, this->settings.init_hbond_from_native)));

               // Add Gaussian interaction if applicable - only residue indices specified
               } else if (settings.include_ca_contacts && interaction_type == DISTANCE_GAUSSIAN) {
                    this->g_network.push_back(*(new AtomPairGaussian<CHAIN_TYPE> (this->chain, index1, index2, mean, sd)));

               // Add disulfide interaction if applicable     
               } else if (settings.include_ss_bond && interaction_type == SS_BOND) {
                    this->ss_network.push_back(*(new AtomPairSSBond<CHAIN_TYPE> (this->chain, index1, index2, mean, sd)));
               } else {
                    std::cerr << "Interaction Type unknown in \"" << filename << "\"  line: \"" << temp << "\" \n";
               }
          }

          // lets calculate the stability of the
          // the hbonds we found
          this->calc_desolvation();

          // Close the file
          ifs.close();
     }


     //! Initialize network without specific input.
     //! Make a guess at what a usefull network would look like
     void init_network() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Import local definitions
          using namespace constrain_distances;

          if (this->settings.verbose)
               std::cout << "# Initializing Network :";
          this->hb_network.clear();
          this->g_network.clear();
          this->ss_network.clear();

          this->chain->update_positions();

          for (ResidueIterator<CHAIN_TYPE> res1(*(this->chain)); !res1.end(); ++res1) {
               for (ResidueIterator<CHAIN_TYPE> res2(*(this->chain)); !res2.end(); ++res2) {
                    // look at each pair only once
                    if (res2->index <= res1->index)
                         continue;

                    // Backbone H-Bond
                    if (settings.include_bb_hbond && (res1)->has_atom(H) && (res2)->has_atom(O) && is_bb_hbond_dssp((*res1)[H], (*res2)[O])) {
                         // one way (H->O)
                         this->hb_network.push_back(AtomPairHBond<CHAIN_TYPE> ((*res1)[H], (*res2)[O], this->settings.init_hbond_from_native));
                    }

                    if (settings.include_bb_hbond && (res2)->has_atom(H) && (res1)->has_atom(O) && is_bb_hbond_dssp((*res2)[H], (*res1)[O])) {
                         // or the other (O->H)
                         this->hb_network.push_back(AtomPairHBond<CHAIN_TYPE> ((*res2)[H], (*res1)[O], this->settings.init_hbond_from_native));
                    }

                    if (settings.include_ca_contacts && (res1)->has_atom(CA) && (res2)->has_atom(CA) && is_ca_contact((*res1)[CA], (*res2)[CA])) {
                         this->g_network.push_back(AtomPairGaussian<CHAIN_TYPE> ((*res1)[CA], (*res2)[CA]));
                    }

                    // and look for SS bonds
                    if (settings.include_ss_bond && (res1)->has_atom(SG) && (res2)->has_atom(SG) && is_ss_bond((*res1)[SG], (*res2)[SG])) {
                         this->ss_network.push_back(AtomPairSSBond<CHAIN_TYPE> ((*res1)[SG], (*res2)[SG]));
                         // we need to additionally fix the geometry of the entire bond
                         this->g_network.push_back(AtomPairGaussian<CHAIN_TYPE> ((*res1)[CB], (*res2)[SG], 3.0, 0.3));
                         this->g_network.push_back(AtomPairGaussian<CHAIN_TYPE> ((*res1)[SG], (*res2)[CB], 3.0, 0.3));
                         this->g_network.push_back(AtomPairGaussian<CHAIN_TYPE> ((*res1)[CB], (*res2)[CB], 3.8, 0.38));
                    }

                    if (settings.include_bb_sc_hbond) {
                         // and look for main chain -- side chain Hbond
                         Atom *atomN = NULL;
                         Atom *atomH = NULL;
                         Atom *atomC = NULL;
                         Atom *atomO = NULL;
                         is_bb_sc_hbond_dssp(*res1, *res2, &atomN, &atomH, &atomC, &atomO);

                         if (atomN != NULL && atomH != NULL && atomC != NULL && atomO != NULL) {
                              this->hb_network.push_back(AtomPairHBond<CHAIN_TYPE> (atomN, atomH, atomC, atomO, BB_SC_HBOND, this->settings.init_hbond_from_native));
                         }
                    }
                    if (settings.include_sc_hbond) {
                         // and look for intra sidechain Hbond
                         Atom *atomN = NULL;
                         Atom *atomH = NULL;
                         Atom *atomC = NULL;
                         Atom *atomO = NULL;
                         is_sc_hbond_dssp(*res1, *res2, &atomN, &atomH, &atomC, &atomO);

                         if (atomN != NULL && atomH != NULL && atomC != NULL && atomO != NULL) {
                              this->hb_network.push_back(AtomPairHBond<CHAIN_TYPE> (atomN, atomH, atomC, atomO, SC_HBOND, this->settings.init_hbond_from_native));
                         }
                    }
               }
          }

          // finally lets check whether we know about
          // the pdb file and whether it contains
          // ligands 
          if (settings.include_fixed_point_contacts && settings.pdb_file != "") {
               if (settings.verbose)
                    std::cout << "# Checking for ligands in the pdb file (" << settings.pdb_file << ")\n";

               std::ifstream ifs(settings.pdb_file.c_str());
               if (ifs.fail()) {
                    std::cout << "# Error reading " << settings.pdb_file << " .. skipping this step.\n";
               } else {
                    std::string temp;
                    // int i=0;
                    while (getline(ifs, temp)) {
                         //  1        10        20        30        40        50        60
                         // HETATM 3443  O1E AP5 A 215      24.469  42.554  19.660  0.50 34.55           O
                         if (temp.substr(0, 6) != "HETATM")
                              continue;
                         double x = atof(temp.substr(30,8).c_str());
                         double y = atof(temp.substr(38,8).c_str());
                         double z = atof(temp.substr(46,8).c_str());
                         int id = atoi(temp.substr(6,5).c_str());
                         int index = atoi(temp.substr(22,4).c_str())-1;
                         std::string name = temp.substr(12,4);
                         DummyAtom try_me(name, Vector_3D(x,y,z), index,id);
                         for (ResidueIterator<CHAIN_TYPE> res1(*(this->chain)); !res1.end(); ++res1) {
                              for (AtomIterator<CHAIN_TYPE,ALL> at1(*res1); !at1.end(); ++at1) {
                                   if (is_atom_XH(at1->atom_type))
                                        continue;
                                   if (get_distance(&try_me, &(*at1)) < 4) {
                                        this->gfix_network.push_back(AtomPairGaussianFixedPoint<CHAIN_TYPE> (new DummyAtom(name, Vector_3D(x,y,z), index, id), &(*at1)));
                                   }
                              }
                         }
                    }
               }
          }

          // lets calculate the stability of the
          // the hbonds we found
          this->calc_desolvation();

          // Prune networks (for instance making sure
          // that we only have a single bond per atom.
          this->prune_network();
     }


     //! Evaluate energy term.
     //!
     //! Evaluates the energy of all the networks.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info = NULL) {
          double e = 0;

          // if we are not supposed to use caching
          // we just ignore the move_info
          if (not settings.use_caching)
               move_info = NULL;

          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               double tmp = this->hb_network[i].get_energy(move_info);
               e += tmp;
          }

          for (unsigned int i = 0; i < this->g_network.size(); i++) {
               double tmp = this->g_network[i].get_energy(move_info);
               e += tmp;
          }

          for (unsigned int i = 0; i < this->ss_network.size(); i++) {
               double tmp = this->ss_network[i].get_energy(move_info);
               e += tmp;
          }

          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               double tmp = this->gfix_network[i].get_energy(move_info);
               e += tmp;
          }

          return e;
     }

     //! Display the networks.
     //! This output can be used as input.
     void print_network() {

          std::cout << "\n# [ Distance Constrain Network Start ]" << std::endl;
          std::cout << "\n# "<< this->hb_network.size() << " hbonds found. " << std::endl;
          std::cout << "\n# HBonds : " << std::endl;
          double tmpE = this->evaluate();
          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               std::cout << "" << this->hb_network[i].interaction_type << "  " << this->hb_network[i].first->residue->index << "  "
                         << this->hb_network[i].second->residue->index << "  " << this->hb_network[i].contribution << "  "
                         << this->hb_network[i].atomN->atom_type << "  " << this->hb_network[i].atomH->atom_type << "  " << this->hb_network[i].atomC->atom_type
                         << "  " << this->hb_network[i].atomO->atom_type << "  "
                         << (this->hb_network[i].NO->mean) << "  "
                         << (this->hb_network[i].NC->mean) << "  "
                         << (this->hb_network[i].HO->mean) << "  "
                         << (this->hb_network[i].HC->mean) << "  "
                         << this->hb_network[i].desolvation_factor <<  std::endl;
          }

          std::cout << "\n# Gaussian distances : " << std::endl;
          for (unsigned int i = 0; i < this->g_network.size(); i++) {
               std::cout << this->g_network[i].interaction_type << "  " << this->g_network[i].first->residue->index << "  "
                         << this->g_network[i].second->residue->index << " " << this->g_network[i].contribution << "  " 
                         << this->g_network[i].first->atom_type << "  " << this->g_network[i].second->atom_type << "  " 
                         << this->g_network[i].gaussian->mean << "  " << this->g_network[i].gaussian->sd << "  "
                         << std::endl;
          }

          std::cout << "\n# Gaussian distances to fixed points: " << std::endl;
          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               std::cout << this->gfix_network[i].interaction_type << "  " << this->gfix_network[i].first->index << "  "<< this->gfix_network[i].second->residue->index << "  "
                         << this->gfix_network[i].contribution << "  " << this->gfix_network[i].first->atom_name << "  "
                         << this->gfix_network[i].second->atom_type << "  "
                         << this->gfix_network[i].first->position[0] << "  "
                         << this->gfix_network[i].first->position[1] << "  "
                         << this->gfix_network[i].first->position[2] << "  "
                         << this->gfix_network[i].gaussian->mean << "  " << this->gfix_network[i].gaussian->sd << "  "
                         << std::endl;
          }

          std::cout << "\n# Disulfid Bonds  : " << std::endl;
          for (unsigned int i = 0; i < this->ss_network.size(); i++) {
               std::cout << this->ss_network[i].interaction_type << "  " << this->ss_network[i].first->residue->index << "  "
                         << this->ss_network[i].second->residue->index << " " << this->ss_network[i].contribution << "  " 
                         << this->ss_network[i].gaussian->mean << "  " << this->ss_network[i].gaussian->sd << "  " 
                         << std::endl;
          }
          std::cout << "\n#\n";
          std::cout << "\n# Total energy : " << tmpE << "\n";
          std::cout << "#\n";
          std::cout << "# [ Distance Constrain Network End ]" << "\n";
          if (settings.debug)
        	  std::cout << settings << "\n\n";
     }


     //! Dump the network as a Pymol script
     void generate_pymol_network() {

          // lets begin by evaluating the energy once
          this->evaluate();

          std::string fname = "load_pymol_network.py";
          std::ofstream pymol(fname.c_str());

          ChainFB native(settings.pdb_file.c_str());
          char filename[200];
          sprintf(filename, "%s_network.pdb", (settings.pdb_file.substr(0,settings.pdb_file.size()-4)).c_str());
          this->chain->output_as_pdb_file(filename, &native);

          pymol << "# " << std::endl;
          pymol << "# PyMOL script generated by the Typhon Protein Monte Carlo Dynamics Tool" << std::endl;
          pymol << "# Part of the Phaistos structural bioinformatics framework." << std::endl;
          pymol << "# by Tim Harder, Wouter Boomsma, Thomas Hamelryck 2010, 2011" << std::endl;
          pymol << "# " << std::endl;
          pymol << "import sys " << std::endl;
          pymol << "try :" << std::endl;
          pymol << "    from pymol import cmd" << std::endl;
          pymol << "    from pymol.cgo import *" << std::endl;
          pymol << "except ImportError :" << std::endl;
          pymol << "    print 'Cannot load PyMOL libraries.\\nTry > $ pymol "<<fname<<"\\n\\n'  " << std::endl;
          pymol << "    sys.exit()" << std::endl;
          pymol << "cmd.load(\"" << filename <<"\") " << std::endl;
          pymol << "obj = [ " << std::endl;

          // [ CYLINDER, x1, y1, z1, x2, y2, z2, dr, rr, gg, bb, rr, gg, bb ]
          // define the colors
          float rr = 0.9;  //red
          float gg = 0.;   //green
          float bb = 0.;   //blue
          float dr = 0.2;  //stick radius
          // only print two digits after the decimal point
          //
          // origin
          //pymol << "CYLINDER, 0, 0, 0, 0, 0, 3, 0.3, 1, 1, 1, 1, 1, 1, " << std::endl;
          //pymol << "CYLINDER, 0, 0, 0, 0, 3, 0, 0.3, 1, 1, 1, 1, 1, 1, " << std::endl;
          //pymol << "CYLINDER, 0, 0, 0, 3, 0, 0, 0.3, 1, 1, 1, 1, 1, 1, " << std::endl;
          //
          // include all the Hbonds
          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               if (this->hb_network[i].interaction_type == 1) {
                    gg =0.; bb=0.;
               } else if (this->hb_network[i].interaction_type == 2) {
                    gg =0.9; bb=0.;
               } else if (this->hb_network[i].interaction_type == 3) {
                    gg =0.; bb=0.9;
               }
               pymol << "CYLINDER, " ;
               pymol << this->hb_network[i].atomH->position[0] << ", ";
               pymol << this->hb_network[i].atomH->position[1] << ", ";
               pymol << this->hb_network[i].atomH->position[2] << ", ";
               pymol << this->hb_network[i].atomO->position[0] << ", ";
               pymol << this->hb_network[i].atomO->position[1] << ", ";
               pymol << this->hb_network[i].atomO->position[2] << ", ";
               pymol << dr << ", " << rr << ", " << gg << ", " << bb << ", ";
               pymol << rr	<< ", " << gg << ", " << bb << "," << "    # a little test here :-)" << std::endl;
          }

          // all the gaussian contacts
          rr = 0.;  //red
          gg = 0.9;   //green
          bb = 0.;   //blue
          for (unsigned int i = 0; i < this->g_network.size(); i++) {
               pymol << "CYLINDER, " ;
               pymol << this->g_network[i].first->position[0] << ", ";
               pymol << this->g_network[i].first->position[1] << ", ";
               pymol << this->g_network[i].first->position[2] << ", ";
               pymol << this->g_network[i].second->position[0] << ", ";
               pymol << this->g_network[i].second->position[1] << ", ";
               pymol << this->g_network[i].second->position[2] << ", ";
               pymol << dr << ", " << rr << ", " << gg << ", " << bb << ", ";
               pymol << rr	<< ", " << gg << ", " << bb << "," << std::endl;

          }

          // all the fixed point contacts
          rr = 0.; //red
          gg = 0.9; //green
          bb = 0.9; //blue
          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               pymol << "CYLINDER, ";
               pymol << this->gfix_network[i].second->position[0] << ", ";
               pymol << this->gfix_network[i].second->position[1] << ", ";
               pymol << this->gfix_network[i].second->position[2] << ", ";
               pymol << this->gfix_network[i].first->position[0] << ", ";
               pymol << this->gfix_network[i].first->position[1] << ", ";
               pymol << this->gfix_network[i].first->position[2] << ", ";
               pymol << dr << ", " << rr << ", " << gg << ", " << bb << ", ";
               pymol << rr	<< ", " << gg << ", " << bb << "," << std::endl;
          }

          // all the disulfide bridges
          rr = 0.; //red
          gg = 0.; //green
          bb = 0.9; //blue
          for (unsigned int i = 0; i < this->ss_network.size(); i++) {
               pymol << "CYLINDER, ";
               pymol << this->ss_network[i].first->position[0] << ", ";
               pymol << this->ss_network[i].first->position[1] << ", ";
               pymol << this->ss_network[i].first->position[2] << ", ";
               pymol << this->ss_network[i].second->position[0] << ", ";
               pymol << this->ss_network[i].second->position[1] << ", ";
               pymol << this->ss_network[i].second->position[2] << ", ";
               pymol << dr << ", " << rr << ", " << gg << ", " << bb << ", ";
               pymol << rr	<< ", " << gg << ", " << bb << "," << std::endl;
          }
          pymol << "] " << std::endl;
          pymol << "cmd.load_cgo(obj, 'cgo') " << std::endl;
          pymol << std::endl;
          pymol.close();

          if (settings.verbose)
               std::cout << "#" << std::endl << "# generated PyMOL script." << std::endl << "# Try: $> pymol " << fname  << std::endl << "#"  << std::endl;

          return;
     }


     //! Clear the cache
     void clear_cache() {
          //std::cout << "# Clear Cache .. \n";
          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               this->hb_network[i].clear_cache();
          }
          for (unsigned int i = 0; i < this->g_network.size(); i++) {
               this->g_network[i].clear_cache();
          }
          for (unsigned int i = 0; i < this->ss_network.size(); i++) {
               this->ss_network[i].clear_cache();
          }
          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               this->gfix_network[i].clear_cache();
          }
     }

     //! Accept last energy evaluation
     void accept() {
          if (not settings.use_caching)
               return;

          //std::cout << "# Accept .. \n";
          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               this->hb_network[i].accept();
          }
          for (unsigned int i = 0; i < this->g_network.size(); i++) {
               this->g_network[i].accept();
          }
          for (unsigned int i = 0; i < this->ss_network.size(); i++) {
               this->ss_network[i].accept();
          }
          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               this->gfix_network[i].accept();
          }
     }

     //! Reject last energy evaluation
     void reject() {
          if (not settings.use_caching)
               return;

          //std::cout << "# Reject .. \n";
          for (unsigned int i = 0; i < this->hb_network.size(); i++) {
               this->hb_network[i].reject();
          }
          for (unsigned int i = 0; i < this->g_network.size(); i++) {
               this->g_network[i].reject();
          }
          for (unsigned int i = 0; i < this->ss_network.size(); i++) {
               this->ss_network[i].reject();
          }
          for (unsigned int i = 0; i < this->gfix_network.size(); i++) {
               this->gfix_network[i].reject();
          }
     }


     //! Find hydrogen bond donors in current residue.
     //!
     //! \param res Residue
     //! \param donor Output vector of donor atom pairs
     void get_hbond_donor(ResidueFB &res, std::vector<std::pair<Atom*, Atom*> > *donor) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          if (res.residue_type == ARG) {
               if (res.has_atom(NE) && res.has_atom(HE))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NE], res[HE]));
               if (res.has_atom(NH1) && res.has_atom(HH11))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NH1], res[HH11]));
               if (res.has_atom(NH1) && res.has_atom(HH12))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NH1], res[HH12]));
               if (res.has_atom(NH2) && res.has_atom(HH21))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NH2], res[HH21]));
               if (res.has_atom(NH2) && res.has_atom(HH22))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NH2], res[HH22]));

          } else if (res.residue_type == ASN) {
               if (res.has_atom(ND2) && res.has_atom(HD22))
                    donor->push_back(std::pair<Atom*, Atom*>(res[ND2], res[HD22]));
               if (res.has_atom(ND2) && res.has_atom(HD21))
                    donor->push_back(std::pair<Atom*, Atom*>(res[ND2], res[HD21]));

          } else if (res.residue_type == CYS) {
               if (res.has_atom(SG) && res.has_atom(HG))
                    donor->push_back(std::pair<Atom*, Atom*>(res[SG], res[HG]));

          } else if (res.residue_type == GLN) {
               if (res.has_atom(NE2) && res.has_atom(HE21))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NE2], res[HE21]));
               if (res.has_atom(NE2) && res.has_atom(HE22))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NE2], res[HE22]));

          } else if (res.residue_type == HIS) {
               if (res.has_atom(ND1) && res.has_atom(HD1))
                    donor->push_back(std::pair<Atom*, Atom*>(res[ND1], res[HD1]));
               if (res.has_atom(NE2) && res.has_atom(HE2))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NE2], res[HE2]));

          } else if (res.residue_type == LYS) {
               if (res.has_atom(NZ) && res.has_atom(HZ1))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NZ], res[HZ1]));
               if (res.has_atom(NZ) && res.has_atom(HZ2))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NZ], res[HZ2]));
               if (res.has_atom(NZ) && res.has_atom(HZ3))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NZ], res[HZ3]));

          } else if (res.residue_type == SER) {
               if (res.has_atom(OG) && res.has_atom(HG))
                    donor->push_back(std::pair<Atom*, Atom*>(res[OG], res[HG]));

          } else if (res.residue_type == THR) {
               if (res.has_atom(OG1) && res.has_atom(HG1))
                    donor->push_back(std::pair<Atom*, Atom*>(res[OG1], res[HG1]));

          } else if (res.residue_type == TRP) {
               if (res.has_atom(NE1) && res.has_atom(HE1))
                    donor->push_back(std::pair<Atom*, Atom*>(res[NE1], res[HE1]));

          } else if (res.residue_type == TYR) {
               if (res.has_atom(OH) && res.has_atom(HH))
                    donor->push_back(std::pair<Atom*, Atom*>(res[OH], res[HH]));

          }
     }

     
     //! Find hydrogen bond acceptors in current residue.
     //!
     //! \param res Residue
     //! \param acceptor Output vector of acceptor atom pairs
     void get_hbond_acceptor(ResidueFB &res, std::vector<std::pair<Atom*, Atom*> > *acceptor) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          if (res.residue_type == ASN) {
               if (res.has_atom(CG) && res.has_atom(OD1))
                    acceptor->push_back(std::pair<Atom*, Atom*>(res[CG], res[OD1]));

          } else if (res.residue_type == ASP) {
               if (res.has_atom(CG) && res.has_atom(OD1))
                    acceptor->push_back(std::pair<Atom*, Atom*>(res[CG], res[OD1]));
               if (res.has_atom(CG) && res.has_atom(OD2))
                    acceptor->push_back(std::pair<Atom*, Atom*>(res[CG], res[OD2]));

          } else if (res.residue_type == GLU) {
               if (res.has_atom(CD) && res.has_atom(OE1))
                    acceptor->push_back(std::pair<Atom*, Atom*>(res[CD], res[OE1]));
               if (res.has_atom(CD) && res.has_atom(OE2))
                    acceptor->push_back(std::pair<Atom*, Atom*>(res[CD], res[OE2]));

          } else if (res.residue_type == GLN) {
               if (res.has_atom(CD) && res.has_atom(OE1))
                    acceptor->push_back(std::pair<Atom*, Atom*>(res[CD], res[OE1]));

          }
     }

     //! Check whether interaction is a hydrogen bond according to DSSP criteria.
     //! (Kabsch W, Sander C, Biopolymers, 1983)
     //! \param first First atom
     //! \param second Second atom
     bool is_bb_hbond_dssp(Atom* first, Atom* second) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Import local definitions
          using namespace constrain_distances;

          if ((first->atom_type == H && second->atom_type == O)) {
               double d = get_distance(first, second);
               int skip = abs(first->residue->index - second->residue->index);
               if (skip > settings.bb_hbond_skip && d > 1 && d < 3) {
                    double e = 0;
                    // ok this is a quite promising candidate
                    AtomPairHBond<CHAIN_TYPE> tmp(first, second);
                    e = tmp.get_energy_dssp();
                    if (e < -0.5) {
                         return true;
                    }
               }
          }
          return false;
     }


     //! Check whether two residues have a hydrogen bond according to DSSP criteria.
     //! (Kabsch W, Sander C, Biopolymers, 1983)
     //! \param res1 First residue
     //! \param res2 Second residue
     //! \param donor_atom1 Output donor atom participating in hydrogen bond (e.g. nitrogen)
     //! \param donor_atom2 Output donor atom participating in hydrogen bond (e.g. hydrogen)
     //! \param acceptor_atom1 Output acceptor atom participating in hydrogen bond (e.g. carbon)
     //! \param acceptor_atom2 Output acceptor atom participating in hydrogen bond (e.g. oxygen)
     bool is_sc_hbond_dssp(ResidueFB &res1, ResidueFB &res2, Atom **donor_atom1, Atom **donor_atom2, Atom **acceptor_atom1, Atom **acceptor_atom2) {

          // Import local definitions
          using namespace constrain_distances;

          int skip = abs(res1.index - res2.index);
          if (skip <= settings.sc_hbond_skip) {
               return false;
          }
          std::vector<std::pair<Atom*, Atom*> > donor;
          std::vector<std::pair<Atom*, Atom*> > acceptor;
          donor.clear();
          acceptor.clear();

          // get all the possible donors
          get_hbond_donor(res1, &donor);
          get_hbond_donor(res2, &donor);
          // quick escape
          if (donor.size() < 1) {
               return false;
          }

          // get all the possible acceptors
          get_hbond_acceptor(res1, &acceptor);
          get_hbond_acceptor(res2, &acceptor);

          // Exit if no acceptors are found
          if (acceptor.size() < 1)
               return false;

          for (int i = 0; i < (int) donor.size(); i++) {
               for (int j = 0; j < (int) acceptor.size(); j++) {
                    // we need to redo this .. since we asked both sides
                    // for donor and acceptor atoms and otherwise will
                    // find contacts in the same residue.
                    int skip = abs(donor[i].second->residue->index - acceptor[j].second->residue->index);
                    if (skip <= settings.sc_hbond_skip) {
                         continue;
                    }
                    double d = get_distance(donor[i].second, acceptor[j].second);

                    if (d < 3.) {
                         AtomPairHBond<CHAIN_TYPE> tmp(donor[i].first, donor[i].second, acceptor[j].first, acceptor[j].second, SC_HBOND);
                         double e = tmp.get_energy_dssp();
                         //double e = tmp.get_energy();
                         if (e < -0.5) {
                              *donor_atom1 = donor[i].first;
                              *donor_atom2 = donor[i].second;
                              *acceptor_atom1 = acceptor[j].first;
                              *acceptor_atom2 = acceptor[j].second;
                              return true;
                         }
                    }
               }
          }
          // just to be save since nothing else worked anyways
          *donor_atom1 = NULL;
          *donor_atom2 = NULL;
          *acceptor_atom1 = NULL;
          *acceptor_atom2 = NULL;
          return false;
     }


     //! Check whether two residues have a hydrogen bond according to DSSP criteria.
     //! (Kabsch W, Sander C, Biopolymers, 1983)
     //! \param res1 First residue
     //! \param res2 Second residue
     //! \param donor_atom1 Output donor atom participating in hydrogen bond (e.g. nitrogen)
     //! \param donor_atom2 Output donor atom participating in hydrogen bond (e.g. hydrogen)
     //! \param acceptor_atom1 Output acceptor atom participating in hydrogen bond (e.g. carbon)
     //! \param acceptor_atom2 Output acceptor atom participating in hydrogen bond (e.g. oxygen)
     bool is_bb_sc_hbond_dssp(ResidueFB &res1, ResidueFB &res2, Atom **donor_atom1, Atom **donor_atom2, Atom **acceptor_atom1, Atom **acceptor_atom2) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Import local definitions
          using namespace constrain_distances;

          int skip = abs(res1.index - res2.index);
          if (skip <= settings.bb_sc_hbond_skip) {
               return false;
          }
          std::vector<std::pair<Atom*, Atom*> > donor;
          std::vector<std::pair<Atom*, Atom*> > acceptor;
          donor.clear();
          acceptor.clear();

          // get all the possible donors
          get_hbond_donor(res1, &donor);
          get_hbond_donor(res2, &donor);
          if (res1.has_atom(N) && res1.has_atom(H))
               donor.push_back(std::pair<Atom*, Atom*>(res1[N], res1[H]));
          if (res2.has_atom(N) && res2.has_atom(H))
               donor.push_back(std::pair<Atom*, Atom*>(res2[N], res2[H]));
          // quick escape
          if (donor.size() < 3) {
               return false;
          }

          // get all the possible acceptors
          get_hbond_acceptor(res1, &acceptor);
          get_hbond_acceptor(res2, &acceptor);
          if (res1.has_atom(C) && res1.has_atom(O))
               acceptor.push_back(std::pair<Atom*, Atom*>(res1[C], res1[O]));
          if (res2.has_atom(C) && res2.has_atom(O))
               acceptor.push_back(std::pair<Atom*, Atom*>(res2[C], res2[O]));
          // quick escape
          if (acceptor.size() < 3)
               return false;

          for (int i = 0; i < (int) donor.size(); i++) {
               for (int j = 0; j < (int) acceptor.size(); j++) {
                    if (donor[i].first->is_sidechain_atom && acceptor[j].second->is_sidechain_atom) {
                         // already covered those
                         continue;
                    }
                    if (donor[i].first->atom_type == N && acceptor[j].second->atom_type == O) {
                         // already covered those as well
                         continue;
                    }

                    // we need to redo this .. since we asked both sides
                    // for donor and acceptor atoms and otherwise will
                    // find contacts in the same residue.
                    int skip = abs(donor[i].second->residue->index - acceptor[j].second->residue->index);
                    if (skip <= settings.bb_sc_hbond_skip) {
                         continue;
                    }
                    double d = get_distance(donor[i].second, acceptor[j].second);

                    if (d < 3.) {
                         AtomPairHBond<CHAIN_TYPE> tmp(donor[i].first, donor[i].second, acceptor[j].first, acceptor[j].second, BB_SC_HBOND);
                         double e = tmp.get_energy_dssp();
                         //double e = tmp.get_energy();
                         if (e < -0.5) {
                              *donor_atom1 = donor[i].first;
                              *donor_atom2 = donor[i].second;
                              *acceptor_atom1 = acceptor[j].first;
                              *acceptor_atom2 = acceptor[j].second;
                              return true;
                         }
                    }
               }
          }
          // just to be save since nothing else worked anyways
          *donor_atom1 = NULL;
          *donor_atom2 = NULL;
          *acceptor_atom1 = NULL;
          *acceptor_atom2 = NULL;
          return false;
     }

     //! Check whether interaction is a CA-CA contact.
     //! \param first First atom
     //! \param second Second atom
     bool is_ca_contact(Atom* first, Atom* second) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Import local definitions
          using namespace constrain_distances;

          if ((first->atom_type == CA && second->atom_type == CA)) {
               double d = get_distance(first, second);
               int skip = abs(first->residue->index - second->residue->index);
               // lets keep them close in space, but far in the chain

               if (skip > settings.ca_skip && d < settings.ca_distance) {
                    return true;
               }
          }
          return false;
     }

     //! Check whether interaction is a disulfide bond.
     //! \param first First atom
     //! \param second Second atom
     bool is_ss_bond(Atom* first, Atom* second) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Import local definitions
          using namespace constrain_distances;

          if ((first->atom_type == SG && second->atom_type == SG)) {
               double d = get_distance(first, second);
               int skip = abs(first->residue->index - second->residue->index);
               // lets keep them close in space, but far in the chain
               if (skip > settings.ss_skip && d < settings.ss_distance) {
                    return true;
               }
          }
          return false;
     }

};

}

#endif

