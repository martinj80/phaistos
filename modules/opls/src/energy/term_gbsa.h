// gbsa.h --- GBSA implicit solvation energy term
// Copyright (C) 2009-2011 Sandro Bottaro, Kristoffer Enøe Johansson, Wouter Boomsma
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
// Literature references:
//
// The GB/SA Continuum Model for Solvation. A Fast Analytical Method for the Calculation of Approximate Born Radii 
// Di Qiu, Peter S. Shenkin, Frank P. Hollinger, and W. Clark Still
// J. Phys. Chem. A 1997, 101, 3005-3014 
//
// VdW radii are extracted from Tinker sourcecode
//

#ifndef GBSA_H
#define GBSA_H

#include <vector>
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "energy/parameter_data_structures.h"
#include "energy/opls_parameters.h"

#include "protein/iterators/pair_iterator_chaintree.h"

namespace phaistos {

//! Parameter reader and container class for GBSA
class GbsaParameters: public TinkerParameterBase {

private:

     //! Vectors to store partial charges in
     std::vector<double> charge;

     //! 2D binary map to store bond lengths in
     ParameterData2D<double> bond_length_map;
     
public:
     
     //! Main Constructor
     GbsaParameters() {
          
          // calculate charge unit
          charge_unit = 332.05382; //1 e^2/Å
          std::string dielectric_str = read_header(&parameters_opls,"dielectric");
          charge_unit /= atof( dielectric_str.c_str() );
          
          // get parameters from 'charge' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls,&raw_param,"charge",1);
          // do not use index zero
          charge.push_back(0.0);
          // put parameters in map
          unsigned int id(0);
          for (unsigned i=0; i<raw_param.size(); i++) {
               id = (raw_param[i].first)[0];
               double param = raw_param[i].second[0];
               assert(id == charge.size());
               // if (id != charge.size()) {
               //      std::cerr<<"\nOPLS TERM GBSA - read error near charge parameter id"
               //               <<id<<"\n\n";
               // }
               charge.push_back(param);
          }
          max_charge_id = id;

          // Get Bond length params
          raw_param.clear();
          read_param(&parameters_opls,&raw_param,"bond",2);
          // put bond length params in map
          for (unsigned int i=0; i<raw_param.size(); i++) {
               std::vector<int> id = raw_param[i].first;
               double bond_len = raw_param[i].second[1];
               if (bond_length_map(id[0],id[1])) {
                    std::cout<<"\nOPLS TERM GBSA - multiple definitions of bond parameter "<<id[0]<<","<<id[1]<<"\n\n";
                    assert(false);
               } else {
                    bond_length_map(id[0],id[1]) = bond_len;
               }
          }
     };

     //! Destructor
     ~GbsaParameters() {};

     //! Get bond length
     double get_bond_eq_length(Atom *atom1, Atom *atom2) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          if (id1 > id2) {
               int buf=id1; id1=id2; id2=buf;
          }
          return bond_length_map(id1,id2);
     }

     //! Get partial charge
     double get_charge(Atom *atom) {
          int id = get_atom_id(atom);
          assert(id <= max_charge_id);
          return ( charge[id] );
     };
     
     //! Number og partial charge parameters
     int max_charge_id;

     //! Overall scaling factor of the charge energy
     double charge_unit;
};


template <typename DERIVED_CLASS>
class TermGbsaBase: public EnergyTermCommon<DERIVED_CLASS,ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS,ChainFB> EnergyTermCommon;     

protected:

     //! read bond lengths and partial charges from the OPLS forcefield
     GbsaParameters parameters;

     //! Atom volumes. [res index][atoms index in res]
     std::vector<std::vector<double> > volume;

     //! Born radii of all atoms. [res index][atoms index in res]
     std::vector<std::vector<double> > born_radius;

     //! Constant Born radius contribution from the first and second neighbour. [res index][atoms index in res]
     std::vector<std::vector<double> > born_radius_neighbor;

     //! local biotype cache for all atoms. [res index][atoms index in res]
     std::vector<std::vector<int> > biotype;

     //! van der Waals radii from TINKER's implementation
     std::vector<double> vdw_radius;

     //! Calculate non-polar solvation (solvent cavity) energy of an atom
     //! \param atom The atom
     //! \param born_radius Born radius of the atom
     //! \return energy in kcal/mol /(4*pi*sigma) 
     double calc_npol_solv(Atom *atom, double born_radius) {
          static const double r_solv = 1.4;   // A
          double rvdw = vdw_radius[biotype[atom->residue->index][atom->index]];
          double ratio = rvdw/born_radius;
          double pow2 = ratio*ratio;
          //sixth power for testing
          double pow6 = pow2*pow2*pow2;
          return Math<double>::sqr(rvdw + r_solv)*pow6;
     }

     //! Calculate pairwise polar solvation energy
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param born_radius1 Born radius of first atom
     //! \param born_radius2 Born radius of second atom
     //! \param chg1 Partial charge of atom 1 in electron volt
     //! \param chg2 Partial charge of atom 2 in electron volt
     //! \return energy in kcal/mol
     double calc_pol_solv(Atom *atom1, Atom *atom2, 
                          double born_radius1, double born_radius2,
                          const double chg1, const double chg2) {
          /* const double dielectric = 332.05; */
          /* const double epsilon = 78.3; */
          /* const double fac = (-dielectric*(1.0 - 1.0/epsilon)); */
          /* const double off2 = 1E12; // square of cutoff */
          const double fac = -327.809;
     
          double dx = atom1->position[0] - atom2->position[0];
          double dy = atom1->position[1] - atom2->position[1];
          double dz = atom1->position[2] - atom2->position[2];
          double dist2 = dx*dx+dy*dy+dz*dz;
          double c_ij = fac*chg1*chg2;
          double f_ij = 0.0;
          double r_born12 = born_radius1*born_radius2;
          /* if (dist2 <= off2){ */
          f_ij = sqrt(1.0/(dist2+r_born12*std::exp(-(dist2/(4.0*r_born12)))));
          if (atom1 == atom2) c_ij = 0.5*c_ij;
          /* } */
          // assert(std::isfinite(c_ij*f_ij));
          return c_ij*f_ij;
     }


     //! Calculate pairwise (for optimization) contribution of two atoms to born radius
     //! \param atom1 first atom
     //! \param atom2 second atom
     //! \param rvdw1 vdw radius of atom 1 in A
     //! \param rvdw2 vdw radius of atom 2 in A
     //! \param born_radius_element1 Born radius of atom 1 to which terms are added
     //! \param born_radius_element2 Born radius of atom 2 to which terms are added
     //! \return energy in kcal/mol
     void calc_and_add_born_term(Atom *atom1, Atom *atom2,
                                 const double rvdw1, const double rvdw2,
                                 double &born_radius_element1,
                                 double &born_radius_element2) {

          const double p3 = 6.211;
          const double p4 = 15.236;
          const double p5 = 1.254;

          int chain_distance = parameters.get_covalent_dist<ChainFB>(atom1,atom2);
          if(chain_distance < 2) {
               return;
          }

          if (chain_distance < 3) {
               double dx = atom1->position[0] - atom2->position[0];
               double dy = atom1->position[1] - atom2->position[1];
               double dz = atom1->position[2] - atom2->position[2];
               double dist2 = dx*dx+dy*dy+dz*dz;
               double dist4 = dist2*dist2;
               
               born_radius_element1 += (p3*volume[atom2->residue->index][atom2->index])/dist4;
               born_radius_element2 += (p3*volume[atom1->residue->index][atom1->index])/dist4;

          } else {

               assert(rvdw1 > 1E-20);
               assert(rvdw2 > 1E-20);

               double dx = atom1->position[0] - atom2->position[0];
               double dy = atom1->position[1] - atom2->position[1];
               double dz = atom1->position[2] - atom2->position[2];
               double dist2 = dx*dx+dy*dy+dz*dz;
               double dist4 = dist2*dist2;
               double ccf = 1.0;
               double rvdw_sum = rvdw1 + rvdw2;
               double rvdw_sum2 = rvdw_sum*rvdw_sum;
               double ratio = dist2 / rvdw_sum2;
     
               // close contact function correction
               if (ratio < (1.0/p5))
                    ccf = 0.25*Math<double>::sqr(1.0 - cos(ratio*p5*M_PI));

               born_radius_element1 += (p4*volume[atom2->residue->index][atom2->index]*ccf)/dist4;
               born_radius_element2 += (p4*volume[atom1->residue->index][atom1->index]*ccf)/dist4;
          }
     }

     //! Allocate space in vectors
     void allocate_vectors() {
          int biotype_size = 648;

          vdw_radius.resize(biotype_size+1);

          biotype.resize(this->chain->size());
          volume.resize(this->chain->size());
          born_radius.resize(this->chain->size());
          born_radius_neighbor.resize(this->chain->size());
          for (int i=0; i<this->chain->size(); ++i) {
               biotype[i].resize((*this->chain)[i].size());
               volume[i].resize((*this->chain)[i].size());
               born_radius[i].resize((*this->chain)[i].size());
               born_radius_neighbor[i].resize((*this->chain)[i].size(), 0);
          }
     }

     
     //! Fill vector with biotypes for each atom
     void biotypes() {
          for (unsigned int i=0; i<biotype.size(); ++i) {
               for (unsigned int j=0; j<biotype[i].size(); ++j) {
                    Atom *atom = (*this->chain)[i][j];
                    biotype[i][atom->index] = parameters.get_biotype(atom);
                    // biotype[i][atom->index] = parameters.find_biotype(atom);
               }
          }
     }

     //! Fill vector with VdW Radii for each biotype
     int vdw_radii() {
          for (unsigned int index = 1; index < vdw_radius.size(); index++) {
               int atomic_number = parameters.get_atomic_number(index);
               switch(atomic_number){
               case 1:
                    vdw_radius[index] = 1.25;
                    if( parameters.get_tinker_atom_type(index) == "HN" || // backbone
                        index == 151 || //TRP
                        index == 171 || index == 177 || index == 188 || index == 209 || //HYS
                        index == 231 || // ASN
                        index == 257 || // GLN
                        index == 286 || // LYS
                        index == 300 || index == 303) // ARG
                         vdw_radius[index] = 1.15;
                    else if (index == 64 || index == 74 || index == 137)
                         vdw_radius[index] = 1.05;
                    break;
               case 6:
                    vdw_radius[index] = 1.90;
                    if (parameters.get_bond_n(index) == 3)
                         vdw_radius[index] = 1.875;
                    else if (parameters.get_bond_n(index) == 2)
                         vdw_radius[index] = 1.825;
                    break;
               case 7:
                    vdw_radius[index] = 1.7063;
                    if (parameters.get_bond_n(index) == 4)
                         vdw_radius[index] = 1.625;
                    else if (parameters.get_bond_n(index) == 1)
                         vdw_radius[index] = 1.60;
                    break;
               case 8:
                    vdw_radius[index] = 1.535;
                    if (parameters.get_bond_n(index) ==  1)
                         vdw_radius[index] = 1.48;
                    break;
               case 16:
                    vdw_radius[index] = 1.775;
                    break;
               }
          }

          return 0;
     }
    
     //! Fill vector with "exclusive " volumes for every atom in the chain
     int volumes() {
          int size = this->chain->size();
          for (int r=0; r<size; r++) {
               Residue *res = &(*this->chain)[r];

               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom1 = res->atoms[a];
                    double overlap = 0.0;
                    double rvdw = vdw_radius[biotype[atom1->residue->index][atom1->index]];
                    double vol = ((4.0*M_PI*(rvdw*rvdw*rvdw))/3.0);

                    for(CovalentBondIterator<ChainFB> it2(atom1, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY); 
                        !it2.end(); ++it2){
                         //std::cout << it2 << "\n";
                         double rvdw_neigh =  vdw_radius[biotype[it2->residue->index][it2->index]];
                         double mod_dist = 1.01*(parameters.get_bond_eq_length(atom1,&*it2));
                         double h = rvdw*(1.0+(rvdw_neigh*rvdw_neigh - rvdw*rvdw - mod_dist*mod_dist)/(2.0*rvdw*mod_dist));
                         overlap += (M_PI*h*h*(3.0*rvdw-h))/3.0;
                    }

                    assert(std::isfinite(vol-overlap));

                    // Assign the volume
                    volume[atom1->residue->index][atom1->index] = vol-overlap;	
               }
          }
          return 0;
     }

     //! Calculate the contributions to the born radii given by the first and second neighbors of each atom
     int born_radii_neighbors(ChainFB *chain, std::vector<std::vector<double> > &born_radius_neighbor){
            
          const double p1 = 0.073;
          const double p2 = 0.921;
          // const double p3 = 6.211;
          const double dielectric = 332.05;
          const double d_offset = -0.09;

          // for (deprecated::AtomIteratorAll it1 = chain->atomBeginAll(); it1 != chain->atomEndAll(); it1++) {
          for (AtomIterator<ChainFB,definitions::ALL> it1(*chain); !it1.end(); ++it1) {
               Atom *atom1 = &*it1;
               double rvdw_atom1 = vdw_radius[biotype[atom1->residue->index][atom1->index]];

               // Add the first constant contribution
               born_radius_neighbor[atom1->residue->index][atom1->index] += -(0.5*dielectric)/(rvdw_atom1+d_offset+p1);

               // for (deprecated::AtomIteratorAll it2 = it1+1 ; it2 != chain->atomEndAll(); it2++) {
               for (AtomIterator<ChainFB,definitions::ALL> it2(it1+1); !it2.end(); ++it2) {
	  
                    // Atom *atom2 = it2.getAtom();
                    Atom *atom2 = &*it2;
                    int distance = chain_distance<ChainFB>(atom1,atom2);
	  
                    // Add the first neighbors contributions
                    if(distance == 1) {
                         double bond_length = parameters.get_bond_eq_length(atom1,atom2);

                         double bond4_first = bond_length*bond_length*bond_length*bond_length;

                         born_radius_neighbor[atom1->residue->index][atom1->index] += 
                              (p2*volume[atom2->residue->index][atom2->index])/bond4_first;

                         born_radius_neighbor[atom2->residue->index][atom2->index] += 
                              (p2*volume[atom1->residue->index][atom1->index])/bond4_first;

                    }

                    // // Add the second neighbors contributions
                    // else if(chain_distance == 2) {
                    //      /* double dist = (atom1->position - atom2->position).norm(); */
                    //      /* double dist2 = dist*dist; */
                    //      double dx = atom1->position[0] - atom2->position[0];
                    //      double dy = atom1->position[1] - atom2->position[1];
                    //      double dz = atom1->position[2] - atom2->position[2];
                    //      double dist2 = dx*dx+dy*dy+dz*dz;
                    //      double dist4 = dist2*dist2;
	    
                    //      born_radius_neighbor[atom1->residue->index][atom1->index] += (p3*volume[atom2->residue->index][atom2->index])/dist4;
                    //      born_radius_neighbor[atom2->residue->index][atom2->index] += (p3*volume[atom1->residue->index][atom1->index])/dist4;
                    // }
               }
          }
          return 0;
     }
    
     //! Calculate the born radii
     void born_radii(){
      
          const double dielectric = -0.5*332.05;

          // Add the terms from the non-neighboring atoms
          int size = this->chain->size();
          for (int r=0; r<size; r++) {
               Residue *res = &(*this->chain)[r];
               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom1 = res->atoms[a];
                    born_radius[atom1->residue->index][atom1->index] =
                         born_radius_neighbor[atom1->residue->index][atom1->index];
               }
          }

          for (int i=0; i<size; i++) {
               Residue *res1 = &(*this->chain)[i];
               int res1_size = res1->size();

               // i == j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    const double rvdw1 = vdw_radius[biotype[atom1->residue->index][atom1->index]];
                    for (int l=k+1; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         const double rvdw2 = vdw_radius[biotype[atom2->residue->index][atom2->index]];
                         // calc_and_add_born_term(atom1,atom2,rvdw1,rvdw2);
                         calc_and_add_born_term(atom1,atom2,rvdw1,rvdw2,
                                                born_radius[atom1->residue->index][atom1->index],
                                                born_radius[atom2->residue->index][atom2->index]);
                    }
               }
               
               // j != i
               for (int j=i+1; j<size; j++) {
                    Residue *res2 = &(*this->chain)[j];
                    int res2_size = res2->size();
                    
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         const double rvdw1 = vdw_radius[biotype[atom1->residue->index][atom1->index]];
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const double rvdw2 = vdw_radius[biotype[atom2->residue->index][atom2->index]];
                              calc_and_add_born_term(atom1,atom2,rvdw1,rvdw2,
                                                     born_radius[atom1->residue->index][atom1->index],
                                                     born_radius[atom2->residue->index][atom2->index]);
                         }
                    }
               }
          }

          for (int r=0; r<size; r++) {
               Residue *res = &(*this->chain)[r];
               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom1 = res->atoms[a];
                    // std::cout << born_radius[atom1->residue->index][atom1->index] << "\n";
                    born_radius[atom1->residue->index][atom1->index] =
                         dielectric/born_radius[atom1->residue->index][atom1->index];
                    assert(born_radius[atom1->residue->index][atom1->index] >= 0.0);
                    // if (born_radius[atom1->residue->index][atom1->index] < 0) {
                    //      std::cerr << "WARNING: GBSA: born radius<0: " << born_radius[atom1->residue->index][atom1->index] << ". Setting to ";
                    //      born_radius[atom1->residue->index][atom1->index] = 1E-20;
                    //      std::cerr << born_radius[atom1->residue->index][atom1->index] << "\n";
                    // }
               }
          }
     }
    
   
     //! Calculate the polar and non polar part of the solvation energy
     double solvation(){

          const double sigma = 0.0049; // cal/(mol*A^2)
          double e_solv_pol = 0.0;
          double e_solv_npol = 0.0;

          // Looping over all pairs of charged atoms
          int size = this->chain->size();
          for (int i=0; i<size; i++) {
               Residue *res1 = &(*this->chain)[i];
               int res1_size = res1->size();

               // i == j only upper triangle:
               for (int k=0; k<res1_size; k++) {
                    Atom *atom1 = res1->atoms[k];
                    e_solv_npol += calc_npol_solv(atom1, born_radius[atom1->residue->index][atom1->index]);
                    const double chg1 = parameters.get_charge(atom1);
                    if (std::fabs(chg1) < 0.001)
                         continue;
                    for (int l=k; l<res1_size; l++) {
                         Atom *atom2 = res1->atoms[l];
                         const double chg2 = parameters.get_charge(atom2);
                         if (std::fabs(chg2) < 0.001)
                              continue;
                         e_solv_pol += calc_pol_solv(atom1,atom2,
                                                     born_radius[atom1->residue->index][atom1->index],
                                                     born_radius[atom2->residue->index][atom2->index],
                                                     chg1,chg2);
                    }
               }
               
               // j != i
               for (int j=i+1; j<size; j++) {
                    Residue *res2 = &(*this->chain)[j];
                    int res2_size = res2->size();
                    
                    for (int k=0; k<res1_size; k++) {
                         Atom *atom1 = res1->atoms[k];
                         /* e_solv_npol += calc_npol_solv(atom1); */
                         const double chg1 = parameters.get_charge(atom1);
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (int l=0; l<res2_size; l++) {
                              Atom *atom2 = res2->atoms[l];
                              const double chg2 = parameters.get_charge(atom2);
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += calc_pol_solv(atom1,atom2,
                                                          born_radius[atom1->residue->index][atom1->index],
                                                          born_radius[atom2->residue->index][atom2->index],
                                                          chg1,chg2);
                         }
                    }
               }
          }

          e_solv_npol *= 4.0*M_PI*sigma;

          return e_solv_pol + e_solv_npol;
     }

public:

     //! Use same settings as base class
     typedef typename EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;     

     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGbsaBase(ChainFB *chain, std::string name, const Settings &settings=Settings(),
                  RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine) {

          // Annotate chain with biotype
          parameters.annotate_chain(chain);

          // Allocate space in vectors
          allocate_vectors();

          // Create biotypes vector
          biotypes();

          // Get the vdW Radii for each biotype
          vdw_radii();
      
          // Get the volume for each atom in the chain
          volumes();
     
          // Calculate the first neighbour contributions for the born radii
          born_radii_neighbors(chain, born_radius_neighbor);                
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGbsaBase(const TermGbsaBase &other, 
                 RandomNumberEngine *random_number_engine,
                  int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters) {

          // Annotate new chain with biotype
          parameters.annotate_chain(chain);

          // Allocate space in vectors
          allocate_vectors();

          // Create biotypes vector
          biotypes();

          // Get the VdW Radii for each biotype
          vdw_radii();
      
          // Get the volume for each atom in the chain
          volumes();
     
          // Calculate the first neighbour contributions for the born radii
          born_radii_neighbors(chain, born_radius_neighbor);                
     }
     
};
     
//! GBSA implicit solvation energy term
class TermGbsa: public TermGbsaBase<TermGbsa> {
          
     //! For convenience, define local base class
     typedef phaistos::TermGbsaBase<TermGbsa> TermGbsaBase;     

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;     
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGbsa(ChainFB *chain, const Settings &settings=Settings(),
              RandomNumberEngine *random_number_engine = &random_global)
          : TermGbsaBase(chain, "gbsa", settings, random_number_engine) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGbsa(const TermGbsa &other, 
              RandomNumberEngine *random_number_engine,
              int thread_index, ChainFB *chain)
          : TermGbsaBase(other, random_number_engine, thread_index, chain) {
     }

     //! Evaluate
     //! \param move_info object containing information about last move
     //! \return solvation energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
      
          // Update the (global) born radii
          this->born_radii();

          //calculate solvation energy
	  double energy = this->solvation();

	  /* std::cout<<"Continuum Solvation "<<energy<<" kcal/mol hep hep\n"; */
          return energy;
     }

};




//! GBSA implicit solvation energy term - cached version
class TermGbsaCached: public TermGbsaBase<TermGbsaCached> {
public:
     //! For convenience, define local base class
     typedef phaistos::TermGbsaBase<TermGbsaCached> TermGbsaBase;     

private:

     //! Specific CachedIterator for Born radii calculation
     //! This class is basically the same as CachedIteratorVectorBase
     //! but keeps track of the maximum percentwise change in the 
     //! vector of values
     class CachedIteratorBornRadii
          : public chaintree::CachedIteratorVectorBase<CachedIteratorBornRadii,
                                                       ChainFB, Atom, Atom,
                                                       double,
                                                       std::vector<AtomVectorReturnType<ChainFB,double> > >  { 

          //! Define BvType and NodeType locally for ease of reference
          typedef ChainFB::ChainTree::BvType BvType;     
          typedef ChainFB::ChainTree::NodeType NodeType;

          //! 
          double deviation_percentage_cutoff;

     public:

          //! Datastructure keeping track of at which iteration the change in born radii of the different
          //! nodes were last seen to exceed the cutoff
          std::vector<std::vector<long int> > node_status;

          //! Constructor
          //! \param chain Molecule chain
          //! \param deviation_percentage_cutoff Allowed deviation within which cached values will be used
          CachedIteratorBornRadii(ChainFB &chain, 
                                  double deviation_percentage_cutoff=0.0)
               : chaintree::CachedIteratorVectorBase<CachedIteratorBornRadii,
                                                     ChainFB, Atom, Atom,
                                                     double,
                                                     std::vector<AtomVectorReturnType<ChainFB,double> > >(chain),
                 deviation_percentage_cutoff(deviation_percentage_cutoff) {

               // Allocate space for node_status
               node_status.resize(this->ct->get_height());
               for (unsigned int i=0; i<node_status.size(); ++i) {
                    node_status[i].resize(this->ct->get_nodes_at_level(i), 0);
               }
          }

          //! Override the set_cache_entry functionality of the iterator
          //! Registers change in born radii for root node
          void set_cache_entry(NodeType *node1,
                               NodeType *node2,
                               const std::vector<AtomVectorReturnType<ChainFB,double> > &value,
                               bool fully_initialized_subtree=true) {
          
               NodeType *root = this->ct->nodes[this->ct->nodes.size()-1];
               if (node1->level == root->level) {

                    std::vector<AtomVectorReturnType<ChainFB,double> > value_old = get_cache_entry(node1, node2).first;
                    for (uint i=0; i<value.size(); ++i) {
                         double current_max = 0;
                         for (uint j=0; j<value[i].size(); ++j) {
                              double deviation_percentage = std::fabs((value[i][j] - value_old[i][j])/(value_old[i][j]));
                              if (deviation_percentage > current_max)
                                   current_max = deviation_percentage;
                         }
                         if (current_max > deviation_percentage_cutoff)
                              set_node_status(this->ct->nodes[i]);
                    }
               }

               // Call base class version
               chaintree::CachedIteratorVectorBase<CachedIteratorBornRadii,
                    ChainFB,Atom,Atom,
                    double, std::vector<AtomVectorReturnType<ChainFB,double> > >::
                    set_cache_entry(node1, node2, value, fully_initialized_subtree);
          }

          //! Sets node status for current leaf node and propagates informatation
          //! to internal nodes
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
          CachedIteratorBornRadii &operator=(const CachedIteratorBornRadii &other) {
               chaintree::CachedIteratorVectorBase<CachedIteratorBornRadii,
                                                   ChainFB, Atom, Atom,
                                                   double,
                                                   std::vector<AtomVectorReturnType<ChainFB,double> > >::operator=(other);
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
                                                 ChainFB::ChainTree::NodeType,
                                                 ChainFB::ChainTree::NodeType,
                                                 double,double> { 

          //! Define BvType and NodeType locally for ease of reference
          typedef ChainFB::ChainTree::BvType BvType;     
          typedef ChainFB::ChainTree::NodeType NodeType;

          //! This is a reference to the corresponding vector in CachedIteratorBornRadii
          std::vector<std::vector<long int> > &node_status;

     public:

          //! Constructor
          //! \param chain Molecule chain
          //! \param node_status Indicates when nodes have been modified
          CachedIteratorSolvent(ChainFB &chain, 
                                std::vector<std::vector<long int> > &node_status)
               : chaintree::CachedIteratorBase<CachedIteratorSolvent,
                                               ChainFB,
                                               NodeType,
                                               NodeType,
                                               double,
                                               double>(chain),
                 node_status(node_status) {
          }     

          //! Add additional criterion for
          //! excluding child nodes from iteration. In this case, we
          //! exclude a node pair if the maximum born radii change is below cutoff
          //! for both nodes
          //! NOTE that in addition to this criterion, the distance between
          //! the nodes must also be unaltered in order for exclusion to take place
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
                                             double,
                                             double>::operator=(other);
               return *this;
          }
     };



     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;
     typedef AtomVectorReturnType<ChainFB,double> AtomVector;

     //! Cached iterator for born radii calculation
     CachedIteratorBornRadii cached_it_phase1;

     //! Cached iterator for solvent energi calculation
     CachedIteratorSolvent cached_it_phase2;

     //! Chaintree iterator settings 
     chaintree::PairIterator<ChainFB,Atom,Atom>::Settings iterator_settings_phase1;
     chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings iterator_settings_phase2;

     //! Calculate pairwise contribution of two atoms to born radius
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param rvdw1 First van der Waals radius
     //! \param rvdw2 Second van der Waals radius
     //! \param distance Measured distance between atoms
     std::pair<double,double> calc_born_term(Atom *atom1, Atom *atom2,
                                             const double rvdw1, const double rvdw2,
                                             double distance) {

          const double p3 = 6.211;
          const double p4 = 15.236;
          const double p5 = 1.254;

          int chain_distance = this->parameters.get_covalent_dist<ChainFB>(atom1,atom2);
          if(chain_distance < 2)
               return std::make_pair(0,0);

          double dist2 = distance*distance;
          double dist4 = dist2*dist2;

          if (chain_distance < 3) {
               
               return std::make_pair((p3*this->volume[atom2->residue->index][atom2->index])/dist4,
                                     (p3*this->volume[atom1->residue->index][atom1->index])/dist4);
          } else {

               double ccf = 1.0;
               double rvdw_sum = rvdw1 + rvdw2;
               double rvdw_sum2 = rvdw_sum*rvdw_sum;
               double ratio = dist2 / rvdw_sum2;
     
               // close contact function correction
               if (ratio < (1.0/p5))
                    ccf = 0.25*Math<double>::sqr(1.0 - cos(ratio*p5*M_PI));
               
               return std::make_pair((p4*this->volume[atom2->residue->index][atom2->index]*ccf)/dist4,
                                     (p4*this->volume[atom1->residue->index][atom1->index]*ccf)/dist4);
          }
     }


     //! Calculate the born radii
     void born_radii() {

          for (cached_it_phase1(*this->chain, iterator_settings_phase1); 
               !cached_it_phase1.end(); ++cached_it_phase1) {

               Atom *atom1 = cached_it_phase1->first;
               Atom *atom2 = cached_it_phase1->second;

               double distance = cached_it_phase1->distance;

               const double rvdw1 = this->vdw_radius[this->biotype[atom1->residue->index][atom1->index]];
               const double rvdw2 = this->vdw_radius[this->biotype[atom2->residue->index][atom2->index]];

               // Register with cached iterator
               cached_it_phase1.register_contribution(calc_born_term(atom1, atom2, rvdw1, rvdw2, distance));
          }
          // Retrieve result from cached iterator
          std::vector<AtomVector> &sum_cached = cached_it_phase1.compute_total();

          // Post iteration calculations
          const double dielectric = -0.5*332.05;     
          int size = this->chain->chain_tree->get_nodes_at_level(0);
          this->born_radius.resize(size);
          for (int i=0; i<size; ++i) {
               this->born_radius[i].resize(this->chain->chain_tree->nodes[i]->size());
               for (unsigned int k=0; k<this->chain->chain_tree->nodes[i]->size(); ++k) {
                    Atom *atom1 = this->chain->chain_tree->nodes[i]->atoms[k];
                    int atom1_index = this->chain->chain_tree->nodes[i]->get_atom_index(atom1->atom_type);
                    this->born_radius[i][atom1_index] = dielectric/
                         (this->born_radius_neighbor[atom1->residue->index][atom1->index] +
                          sum_cached[i][atom1_index]);
                    if (this->born_radius[i][atom1_index] < 0) {
                         std::cerr << "WARNING: GBSA: born radius<0: " << this->born_radius[i][atom1_index] << ". Setting to ";
                         this->born_radius[i][atom1_index] = 1E-20;
                         std::cerr << this->born_radius[i][atom1_index] << "\n";
                         std::cerr.flush();
                    }
               }
          }
     }


     //! Calculate the born radii - uncached version.
     //! This method is provided for testing purposes only. It calculates the born radii in exactly 
     //! the same way as the cached version, but without caching.
     void born_radii_uncached() {
      
          const double dielectric = -0.5*332.05;
     
          int size = this->chain->chain_tree->get_nodes_at_level(0);
          this->born_radius.resize(size);
          for (int i=0; i<size; ++i) {
               this->born_radius[i].resize(this->chain->chain_tree->nodes[i]->size());
               for (unsigned int k=0; k<this->chain->chain_tree->nodes[i]->size(); ++k) {
                    Atom *atom1 = this->chain->chain_tree->nodes[i]->atoms[k];
                    int atom1_index = this->chain->chain_tree->nodes[i]->get_atom_index(atom1->atom_type);
                    this->born_radius[i][atom1_index] = this->born_radius_neighbor[atom1->residue->index][atom1->index];
               }
          }

          for (int i=0; i<this->chain->chain_tree->get_nodes_at_level(0); ++i) {

               NodeType *node1 = this->chain->chain_tree->nodes[i];

               // i == j only upper triangle:
               for (unsigned int k=0; k<node1->atoms.size(); ++k) {
                    Atom *atom1 = node1->atoms[k];
                    int atom1_index = node1->get_atom_index(atom1->atom_type);
                    const double rvdw1 = this->vdw_radius[this->biotype[atom1->residue->index][atom1->index]];
                    for (unsigned int l=k+1; l<node1->atoms.size(); l++) {
                         Atom *atom2 = node1->atoms[l];
                         int atom2_index = node1->get_atom_index(atom2->atom_type);
                         const double rvdw2 = this->vdw_radius[this->biotype[atom2->residue->index][atom2->index]];
                         this->calc_and_add_born_term(atom1,atom2,rvdw1,rvdw2,
                                                      this->born_radius[i][atom1_index],
                                                      this->born_radius[i][atom2_index]);
                    }
               }

               // j != i
               for (int j=i+1; j<this->chain->chain_tree->get_nodes_at_level(0); ++j) {
                    NodeType *node2 = this->chain->chain_tree->nodes[j];
                    
                    for (unsigned int k=0; k<node1->atoms.size(); k++) {
                         Atom *atom1 = node1->atoms[k];
                         int atom1_index = node1->get_atom_index(atom1->atom_type);
                         const double rvdw1 = this->vdw_radius[this->biotype[atom1->residue->index][atom1->index]];
                         for (unsigned int l=0; l<node2->atoms.size(); l++) {
                              Atom *atom2 = node2->atoms[l];
                              int atom2_index = node2->get_atom_index(atom2->atom_type);
                              const double rvdw2 = this->vdw_radius[this->biotype[atom2->residue->index][atom2->index]];
                              this->calc_and_add_born_term(atom1,atom2,rvdw1,rvdw2,
                                                           this->born_radius[i][atom1_index],
                                                           this->born_radius[j][atom2_index]);
                         }
                    }
               }               
          }

          for (unsigned int i=0; i<this->born_radius.size(); ++i) {
               for (unsigned int k=0; k<this->born_radius[i].size(); ++k) {
                    this->born_radius[i][k] =
                         dielectric/this->born_radius[i][k];
                    assert(this->born_radius[i][k] >= 0.0);
                    // if (this->born_radius[i][k] < 0) {
                    //      std::cerr << "WARNING: GBSA: born radius<0: " << this->born_radius[i][k] << ". Setting to ";
                    //      this->born_radius[i][k] = 1E-20;
                    //      std::cerr << this->born_radius[i][k] << "\n";
                    //      std::cerr.flush();
                    // }
               }
          }
     }

     //! Calculate the polar and non polar part of the solvation energy
     //! \return Solvation energy
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

                         e_solv_npol += this->calc_npol_solv(atom1, this->born_radius[node1->index][atom1_index]);
                         const double chg1 = this->parameters.get_charge(atom1);
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (unsigned int l=k; l<node1->size(); l++) {
                              Atom *atom2 = node1->atoms[l];
                              int atom2_index = node1->get_atom_index(atom2->atom_type);
                              const double chg2 = this->parameters.get_charge(atom2);
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += this->calc_pol_solv(atom1,atom2,
                                                                this->born_radius[node1->index][atom1_index],
                                                                this->born_radius[node1->index][atom2_index],
                                                                chg1,chg2);
                         }                    
                    }
               } else {
                    for (unsigned int k=0; k<node1->size(); ++k) {
                         Atom *atom1 = node1->atoms[k];
                         int atom1_index = node1->get_atom_index(atom1->atom_type);
                         const double chg1 = this->parameters.get_charge(atom1);
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (unsigned int l=0; l<node2->size(); ++l) {
                              Atom *atom2 = node2->atoms[l];
                              int atom2_index = node2->get_atom_index(atom2->atom_type);
                              const double chg2 = this->parameters.get_charge(atom2);
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += this->calc_pol_solv(atom1,atom2,
                                                                this->born_radius[node1->index][atom1_index],
                                                                this->born_radius[node2->index][atom2_index],
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

     //! Calculate the polar and non polar part of the solvation energy
     //! This method is provided for testing purposes only. It calculates the born radii in exactly 
     //! the same way as the cached version, but without caching.
     //! \return Solvation energy
     double solvation_uncached() {

          const double sigma = 0.0049; // cal/(mol*A^2)
          double e_solv_pol = 0.0;
          double e_solv_npol = 0.0;

          // Looping over all pairs of charged atoms
          int size = this->chain->chain_tree->get_nodes_at_level(0);
          for (int i=0; i<size; ++i) {

               NodeType *node1 = this->chain->chain_tree->nodes[i];

               // i == j only upper triangle:
               for (unsigned int k=0; k<node1->size(); ++k) {
                    Atom *atom1 = node1->atoms[k];
                    int atom1_index = node1->get_atom_index(atom1->atom_type);

                    e_solv_npol += this->calc_npol_solv(atom1, this->born_radius[node1->index][atom1_index]);
                    const double chg1 = this->parameters.get_charge(atom1);
                    if (std::fabs(chg1) < 0.001)
                         continue;
                    for (unsigned int l=k; l<node1->size(); l++) {
                         Atom *atom2 = node1->atoms[l];
                         int atom2_index = node1->get_atom_index(atom2->atom_type);
                         const double chg2 = this->parameters.get_charge(atom2);
                         if (std::fabs(chg2) < 0.001)
                              continue;

                         e_solv_pol += this->calc_pol_solv(atom1,atom2,
                                                           this->born_radius[node1->index][atom1_index],
                                                           this->born_radius[node1->index][atom2_index],
                                                           chg1,chg2);
                    }                    
               }

               // j != i
               for (int j=i+1; j<size; ++j) {

                    NodeType *node2 = this->chain->chain_tree->nodes[j];

                    for (unsigned int k=0; k<node1->size(); ++k) {
                         Atom *atom1 = node1->atoms[k];
                         int atom1_index = node1->get_atom_index(atom1->atom_type);
                         /* e_solv_npol += calc_npol_solv(atom1); */
                         const double chg1 = this->parameters.get_charge(atom1);
                         if (std::fabs(chg1) < 0.001)
                              continue;
                         for (unsigned int l=0; l<node2->size(); ++l) {
                              Atom *atom2 = node2->atoms[l];
                              int atom2_index = node2->get_atom_index(atom2->atom_type);
                              const double chg2 = this->parameters.get_charge(atom2);
                              if (std::fabs(chg2) < 0.001)
                                   continue;
                              e_solv_pol += this->calc_pol_solv(atom1,atom2,
                                                                this->born_radius[node1->index][atom1_index],
                                                                this->born_radius[node2->index][atom2_index],
                                                                chg1,chg2);
                         }
                    }
               }               
          }

          e_solv_npol *= 4.0*M_PI*sigma;
          return e_solv_pol + e_solv_npol;
     }

public:

     //! Settings object
     const class Settings: public TermGbsaBase::Settings {

     public:

          //! Maximum deviation allowed in born radii in two subtrees of the chaintree before it is recalculated
          double maximum_deviation_cutoff;

          //! Distance beyond which gbsa contributions are set to zero in phase1.
          double cutoff_distance_phase2;

          //! Distance beyond which gbsa contributions are set to zero in phase2.
          double cutoff_distance_phase1;

          //! Constructor
          Settings(double maximum_deviation_cutoff = 0.1,
                   double cutoff_distance_phase2 = std::numeric_limits<double>::infinity(),
                   double cutoff_distance_phase1 = std::numeric_limits<double>::infinity())
               : maximum_deviation_cutoff(maximum_deviation_cutoff),
                 cutoff_distance_phase2(cutoff_distance_phase2),
                 cutoff_distance_phase1(cutoff_distance_phase1) {}

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "maximum-deviation-cutoff:" << settings.maximum_deviation_cutoff << "\n";
               o << "cutoff-distance-phase1:" << settings.cutoff_distance_phase1 << "\n";
               o << "cutoff-distance-phase2:" << settings.cutoff_distance_phase2 << "\n";
               o << (TermGbsaBase::Settings &)settings;
               return o;
          }                    
     } settings;   //!< Local settings object      


     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGbsaCached(ChainFB *chain,
                    const Settings &settings=Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : TermGbsaBase(chain, "gbsa-cached", settings, random_number_engine),
            cached_it_phase1(*chain, settings.maximum_deviation_cutoff),
            cached_it_phase2(*chain, cached_it_phase1.node_status),
            settings(settings) {

          this->born_radius.clear();

          // Only evaluate modified pairs
          bool only_modified_pairs = true;
          iterator_settings_phase1 = chaintree::PairIterator<ChainFB,Atom,Atom>::Settings(settings.cutoff_distance_phase1,
                                                                                          only_modified_pairs);
          // Evaluate all pairs because born radii changes occur non-locally
          only_modified_pairs = false;
          iterator_settings_phase2 = chaintree::PairIterator<ChainFB,NodeType,NodeType>::
               Settings(settings.cutoff_distance_phase2,
                        only_modified_pairs);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGbsaCached(const TermGbsaCached &other, 
                    RandomNumberEngine *random_number_engine,
                    int thread_index, ChainFB *chain)
          : TermGbsaBase(other,random_number_engine,thread_index,chain),
            cached_it_phase1(*chain, settings.maximum_deviation_cutoff),
            cached_it_phase2(*chain, cached_it_phase1.node_status),
            iterator_settings_phase1(other.iterator_settings_phase1),
            iterator_settings_phase2(other.iterator_settings_phase2),
            settings(other.settings) {}


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return solvation energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          // Update the (global) born radii
          this->born_radii();
          // born_radii_uncached();

          //calculate solvation energy
	  double energy = this->solvation();
	  // double energy = solvation_uncached();

          return energy;
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
    
