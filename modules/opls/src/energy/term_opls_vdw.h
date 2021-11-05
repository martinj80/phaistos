// opls_vdw.h --- OPLS van der Waals interaction energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson, Wouter Boomsma
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

#ifndef OPLS_VDW_H
#define OPLS_VDW_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "energy/opls_parameters.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace phaistos {

//! Parameter container for OPLS van der Waals pair energy term
//! This class reads the OPLS parameter file, processes the parameters and stores them in a
//! biotype*biotype matrix for fast lookup
class VdwParameters: public TinkerParameterBase {

public:

     //! Small class to store a single pair parameter
     class Parameter {

     public:

          //! Geometric mean of the two atoms radii squared
          //! note that sigma values are given in the OPLS parameter file
          double radius_sq;

          //! Geometric mean of the two atoms epsilon value
          double epsilon;
          
          //! Boolean interpretation signaling if the parameter is non-zero
          //! radius_sq=0.0 is used as off signal in the parameter file
          operator bool () {return (radius_sq > -0.0001); };
          
          //! Default constructor
          //! radius_sq=-1 signals uninitialized
          Parameter():radius_sq(-1.0) {};
          
          //! Main constructor
          Parameter(double radius_sq,double eps): radius_sq(radius_sq), epsilon(eps) {};
     };

private:

     //! Vector to store pair parameters in
     std::vector<std::vector<Parameter> > param_cache;

public:

     //! Constructor
     VdwParameters() {

          // Sixth Root of Two - always useful
          const double srt = 1.122462048309372981;

          // get third neighbour scale factor from parameter file header
          vdw14scale = 1.0;
          std::string vdw14scale_str = read_header(&parameters_opls,"vdw-14-scale");
          vdw14scale /= atof( vdw14scale_str.c_str() );

          // get parameters from 'vdw' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls,&raw_param,"vdw",1);

          // convert sigma values to radius
          for (unsigned int i=0; i<raw_param.size(); i++) {
               /* radiustype              SIGMA    */
               raw_param[i].second[0] *= srt;
               /* radiussize              DIAMETER */
               raw_param[i].second[0] *= 0.5;
          }

          // store as pairwise parameters in a biotype*biotype matrix for fast lookup
          unsigned int id_i(0);
          std::vector<Parameter> buf;
          for (unsigned int i=0; i<raw_param.size(); i++) {
               id_i = (raw_param[i].first)[0];
               std::vector<double> param_i = raw_param[i].second;
               for (unsigned int j=i; j<raw_param.size(); j++) {
                    unsigned int id_j = (raw_param[j].first)[0];
                    std::vector<double> param_j = raw_param[j].second;
                    // check if vdw param_id is consequitive
                    assert(id_j == buf.size()+1+i);
                    /* radiusrule  GEOMETRIC */
                    double radius = 2*sqrt(param_i[0]*param_j[0]);
                    double radius_sq = radius*radius;
                    /* epsilonrule GEOMETRIC */
                    double epsilon = sqrt( param_i[1] * param_j[1]);
                    buf.push_back( Parameter(radius_sq,epsilon) );
               }
               // check if vdw param_id is consequitive
               assert(id_i == param_cache.size()+1);
               param_cache.push_back(buf);
               buf.clear();
          }
          max_param_id = id_i;
     };

     //! Destructor
     ~VdwParameters() {};

     //! Parameter getter
     VdwParameters::Parameter get(Atom *atom1, Atom *atom2) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          assert(id1 <= max_param_id);
          assert(id2 <= max_param_id);

          // param_cache is symmetric
          if (id1 > id2) {
               int buf=id1; id1=id2; id2=buf;
          }
          
          // param_cache is zero-based
          id1--;
          id2--;
          
          // convert upper right to lower left triangle for easyer iteration
          id2 -= id1;
          return param_cache[id1][id2];
     };
     
     //! Number of parameters
     int max_param_id;

     //! Scaling factor for third neighbour given in parameter file header
     double vdw14scale;
};



//! OPLS van der Waals interaction term
class TermOplsVdw: public EnergyTermCommon<TermOplsVdw, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsVdw, ChainFB> EnergyTermCommon;

     //! van der Waals pair parameters
     //! initialized on construction
     VdwParameters parameters;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsVdw(ChainFB *chain,
                 const Settings &settings = Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "opls-vdw", settings, random_number_engine) {

          // Annotate chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsVdw(const TermOplsVdw &other, 
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters),
            counter(other.counter) {

          // Annotate the new chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Evaluate van der Waal interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \return van der Waals energy in kcal/mol
     double calc_vdw_pair_energy(Atom *atom1, Atom *atom2) {
          int d = parameters.get_covalent_dist<ChainFB>(atom1,atom2); //6s @ 500 iterations

          const VdwParameters::Parameter param = parameters.get(atom1,atom2);
          //radius=0 for charged hydrogens is 'term off' signal
          if (param.radius_sq < 0.001) //1s @ 500 iterations
               return 0.0;
          if (d > 3) {
               return calc_lennard_jones_energy(atom1,atom2,param);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3 
               return calc_lennard_jones_energy(atom1,atom2,param) * parameters.vdw14scale;
          }
     }
     
     //! Evaluate Lennard-Jones potential
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param param van der Waals parameters
     //! \return lennard Jones contribution in kcal/mol
     double calc_lennard_jones_energy(Atom *atom1, Atom *atom2,
                                      const VdwParameters::Parameter param) {
          counter++;
          const double dx = atom1->position[0] - atom2->position[0];
          const double dy = atom1->position[1] - atom2->position[1];
          const double dz = atom1->position[2] - atom2->position[2];
          const double r_sq = dx*dx + dy*dy + dz*dz;
          const double ratio = param.radius_sq / r_sq;
          const double pow6 = ratio*ratio*ratio;
          return ( param.epsilon*(pow6*pow6-2*pow6) );
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_sum = 0.0;
          this->counter = 0;

          int size = (this->chain)->size();
          for (int i=0; i<size; i++) {
               Residue *res1 = &(*(this->chain))[i];
               int res1_size = res1->size();

               // i==j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    for (int l=k+1; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         energy_sum += calc_vdw_pair_energy(atom1,atom2);
                    }
               }

               // rest
               for (int j=i+1; j<size; j++) {
                    Residue *res2 = &(*(this->chain))[j];
                    int res2_size = res2->size();

                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              energy_sum += calc_vdw_pair_energy(atom1,atom2);
                         }
                    }
               }
          }

          return energy_sum;
     }
};       


//! OPLS van der Waals pair interaction term -- Cached version
class TermOplsVdwCached: public EnergyTermCommon<TermOplsVdwCached, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsVdwCached, ChainFB> EnergyTermCommon;            

     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;

     //! Cached iterator (used when iterating over node pairs)
     CachedIterator<chaintree::PairIterator<ChainFB,NodeType,NodeType> > cached_it;

     //! Chaintree iterator settings
     chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings iterator_settings;

     //! van der Waals pair parameters
     VdwParameters parameters;

public:
     
     //! Local settings class.     
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Cutoff distance beyond which energies are 0
          double cutoff_distance;

          //! Constructor
          Settings(double cutoff_distance = std::numeric_limits<double>::infinity())
               : cutoff_distance(cutoff_distance) {}

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "cutoff-distance:" << settings.cutoff_distance << "\n";
               o << static_cast<const EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }                    
     } settings;     //!< Local settings object 

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsVdwCached(ChainFB *chain, 
                       const Settings &settings = Settings(),
                       RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "opls-vdw-cached", settings, random_number_engine),
            cached_it(*chain),
            settings(settings) {

          // CACHE SETTINGS

          // annotate chain with biotype
          parameters.annotate_chain(chain);

          // Only evaluate modified pairs (as always when caching)
          bool only_modified_pairs = true;

          iterator_settings = chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings(settings.cutoff_distance,
                                                                                           only_modified_pairs);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsVdwCached(const TermOplsVdwCached &other, 
                       RandomNumberEngine *random_number_engine,
                       int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            iterator_settings(other.iterator_settings),
            parameters(other.parameters),
            settings(other.settings) {
          
          // annotate the new chain with biotype
          parameters.annotate_chain(chain);
     }

     //! Evaluate van der Waal interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param distance Distance between atom1 and atom2
     //! \return van der Waals energy in kcal/mol
     double calc_vdw_pair_energy(Atom *atom1, Atom *atom2, double distance) {
          int d = this->parameters.get_covalent_dist<ChainFB>(atom1,atom2); //6s @ 500 iterations

          const VdwParameters::Parameter param = this->parameters.get(atom1,atom2);
          //radius=0 for charged hydrogens is 'term off' signal
          if (param.radius_sq < 0.001) //1s @ 500 iterations
               return 0.0;
          if (d > 3) {
               return calc_lennard_jones_energy(atom1,atom2,distance,param);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3 
               return calc_lennard_jones_energy(atom1,atom2,distance,param) * this->parameters.vdw14scale;
          }
     }
     
     //! Evaluate Lennard-Jones potential
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param distance Distance between atom1 and atom2
     //! \param param van der Waals parameters
     //! \return lennard Jones contribution in kcal/mol
     double calc_lennard_jones_energy(Atom *atom1, Atom *atom2, double distance,
                                      const VdwParameters::Parameter param) {

          // const double dx = atom1->position[0] - atom2->position[0];
          // const double dy = atom1->position[1] - atom2->position[1];
          // const double dz = atom1->position[2] - atom2->position[2];
          // const double r_sq = dx*dx + dy*dy + dz*dz;
          const double r_sq = distance*distance;
          const double ratio = param.radius_sq / r_sq;
          const double pow6 = ratio*ratio*ratio;
          return ( param.epsilon*(pow6*pow6-2*pow6) );
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          // Cached iterator evaluation
          // We iterate over nodes - manually iterating over the atoms
          // this is slightly faster than using the chaintree iterator to directly
          // iterate over atom pairs (see example below)
          for (cached_it(*this->chain, iterator_settings); !cached_it.end(); ++cached_it) {

               double contribution = 0.0;
               if (cached_it->first == cached_it->second) {
                    for (unsigned int i=0; i<cached_it->first->size(); ++i) {
                         for (unsigned int j=i+1; j<cached_it->second->size(); ++j) {
                              double distance = ((*cached_it->first)[i]->position - 
                                                 (*cached_it->second)[j]->position).norm();
                              if (distance < settings.cutoff_distance) {
                                   contribution += calc_vdw_pair_energy((*cached_it->first)[i], 
                                                                        (*cached_it->second)[j],
                                                                        distance);
                              }
                         }
                    }
               } else {
                    for (unsigned int i=0; i<cached_it->first->size(); ++i) {
                         for (unsigned int j=0; j<cached_it->second->size(); ++j) {
                              double distance = ((*cached_it->first)[i]->position - 
                                                 (*cached_it->second)[j]->position).norm();
                              if (distance < settings.cutoff_distance) {
                                   contribution += calc_vdw_pair_energy((*cached_it->first)[i], 
                                                                        (*cached_it->second)[j],
                                                                        distance);
                              }
                         }
                    }
               }
               cached_it.register_contribution(contribution);
          }
          double energy = cached_it.compute_total();

          // Cached iterator evaluation
          // Here we iterate directly over atom-pairs. The code is simpler,
          // but slightly slower
          // for (cached_it(*this->chain, iterator_settings); !cached_it.end(); ++cached_it) {
          //      double contribution = calc_vdw_pair_energy(cached_it->first, 
          //                                                 cached_it->second,
          //                                                 cached_it->distance);
          //      cached_it.register_contribution(contribution);
          // }
          // double energy = cached_it.compute_total();

          // Check with noncached version
          // double uncached_energy = TermOplsVdw<ChainFB>::evaluate();
          // if (!(fabs(energy-uncached_energy)<(0.001*std::max(fabs(energy),fabs(uncached_energy))))) {
          //      std::cout << energy << " " << uncached_energy << "\n";
          //      assert(false);
          // }          
          
          return energy;
     }

     //! Accept last energy evaluation
     void accept() {
          cached_it.accept();
     }

     //! Reject last energy evaluation
     void reject() {
          cached_it.reject();
     }
};       

}

#endif
