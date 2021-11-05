// opls_imptor.h --- OPLS energy: improper torsion or out-of-plane bending energy term
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

#ifndef OPLS_IMPTOR_H
#define OPLS_IMPTOR_H

#include "energy/energy_term.h"

namespace phaistos {

//! Parameter reader and container for OPLS torsion energy term
class ImptorParameters: public TinkerParameterBase {

public:

     //! Parameter container
     class Parameter {

     private:

          //! Which of the 3 components to write next          
          int write_next;
          
     public:
          
          //! The order of the 3 components
          int o[3];

          //! Amplitude of the 3 components
          double amp[3];

          //! Phase of the 3 components
          double ang[3];

          //! Boolean interpretation          
          operator bool () {return (o[0] >= 0);}; // return false only if default constructed

          //! Check if more components can be written to the parameter containre
          bool good() {
               return (write_next<3);
          }

          //! How many components in parameter container
          unsigned int size() {
               return write_next;
          }
          
          //! Fill in a component
          void assign(int order, double aamp, double aang) {
               if ( this->good() ) {
                    o[write_next] = order;
                    amp[write_next] = aamp;
                    ang[write_next] = aang;
                    write_next++;
               }
          }
          
          //! Reset container
          void reset() {
               write_next=0;
               for (int i=0; i<3; i++) {
                    o[i] = 0;
                    amp[i] = 0.0;
                    ang[i] = 0.0;
               }
          };
    
          //! Default constructor
          Parameter(): write_next(0) {
               for(int i=0; i<3; i++) {
                    o[i] = 0;
                    amp[i] = 0.0;
                    ang[i] = 0.0;
               }
               // mark as default constructed
               o[0] = -1;
          }
    
          //! Construcor to fill in a single component
          Parameter(int o1, double amp1, double ang1): write_next(1) {
               o[0] = o1; o[1] = 0; o[2] = 0;
               amp[0] = amp1; amp[1] = 0.0; amp[2] = 0.0;
               ang[0] = ang1; ang[1] = 0.0; ang[2] = 0.0;
          }
          
          //! Construcor to fill in two components
          Parameter(int o1, double amp1, double ang1,
                    int o2, double amp2, double ang2): write_next(2) {
               o[0] = o1; o[1] = o2; o[2] = 0;
               amp[0] = amp1; amp[1] = amp2; amp[2] = 0.0;
               ang[0] = ang1; ang[1] = ang2; ang[2] = 0.0;
          }
          
          //! Construcor to fill in all three components
          Parameter(int o1, double amp1, double ang1,
                            int o2, double amp2, double ang2,
                            int o3, double amp3, double ang3): write_next(3) {
               o[0] = o1; o[1] = o2; o[2] = 0;
               amp[0] = amp1; amp[1] = amp2; amp[2] = amp3;
               ang[0] = ang1; ang[1] = ang2; ang[2] = ang3;
          }
          
          //! Output torsion param as string
          friend std::ostream &operator<<(std::ostream &o, Parameter p) {
               o<<"order "<<p.o[0]<<": "<<p.amp[0]<<","<<p.ang[0]<<". order "<<
                    p.o[1]<<": "<<p.amp[1]<<","<<p.ang[1]<<". order "<<
                    p.o[2]<<": "<<p.amp[2]<<","<<p.ang[2]<<"\n";
               return o;
          }
          
     };
     
private:

     //! Map to store parameters in
     ParameterData4D<Parameter> parameter_map;     
     
public:

     //! Constructor
     ImptorParameters() {

          // get parameters from 'imptors' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > raw_param;
          read_param(&parameters_opls,&raw_param,"imptors",4);

          // store parameters
          for (unsigned int i=0; i<raw_param.size(); i++) {
               Parameter param;
               std::vector<int> id = (raw_param[i].first);
               std::vector<double> param_val = raw_param[i].second;
               unsigned int n_param = param_val.size();
               if (n_param == 0) {
                    // make empty parameter
                    param.assign(0,0.0,0.0);
               } else {
                    for (unsigned int k=0; k<n_param; k+=3) {
                         param_val[k+1] = deg2rad(param_val[k+1]);
                         param.assign( (int)param_val[k+2], param_val[k], param_val[k+1] );
                    }
               }
               // Check for previous definitions of this parameter
               assert( ! parameter_map(id[0],id[1],id[2],id[3]) );
               // if ( parameter_map(id[0],id[1],id[2],id[3]) ) {
               //      std::cerr<<"\nOPLS TERM IMPTOR - multiple definitions of parameter "
               //               <<id[0]<<","<<id[1]<<","<<id[2]<<","<<id[3]<<"\n\n";
               // }
               parameter_map(id[0],id[1],id[2],id[3]) = param;
          }
     }

     //! Destructor
     ~ImptorParameters() {};

     //! Getter
     std::pair<Parameter, std::vector<Atom*> > get_param_and_sorted_atoms(Atom *atom1, Atom *atom2,
                                                                          Atom *atom3, Atom *atom4) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          int id3 = get_param_id(atom3);
          int id4 = get_param_id(atom4);

          std::vector<Atom*> sorted_atoms;
          sorted_atoms.push_back(atom1);
          sorted_atoms.push_back(atom2);
          sorted_atoms.push_back(atom3);
          sorted_atoms.push_back(atom4);
          Parameter param = parameter_map(id1,id2,id3,id4);
          // Filter out parameters that are believed to be wrong in the parameter file
          if ((id1==28 && id2==28 && id3==27 && id4==73) ||
              (id1==76 && id2==37 && id3==40 && id4==40) ||
              (id1==67 && id2==22 && id3==21 && id4==23) ||
              (id1==72 && id2==5  && id3==34 && id4==35))
               return(std::pair<Parameter, std::vector<Atom*> > (param,sorted_atoms) );
          if (!param) {
               param = parameter_map(id2,id1,id3,id4);
               sorted_atoms[0] = atom2;
               sorted_atoms[1] = atom1;
          }
          if (!param) {
               param = parameter_map(id2,id4,id3,id1);
               sorted_atoms[1] = atom4;
               sorted_atoms[3] = atom1;
          }
          if (!param) {
               param = parameter_map(id4,id2,id3,id1);
               sorted_atoms[0] = atom4;
               sorted_atoms[1] = atom2;
          }
          if (!param) {
               param = parameter_map(id4,id1,id3,id2);
               sorted_atoms[1] = atom1;
               sorted_atoms[3] = atom2;
          }
          if (!param) {
               param = parameter_map(id1,id4,id3,id2);
               sorted_atoms[0] = atom1;
               sorted_atoms[1] = atom4;
          }
          assert( param );
          // if (!param) {
          //      std::cerr<<"\nOPLS TERM IMPTOR - missing parameter for param id "
          //               <<id1<<","<<id2<<","<<id3<<","<<id4<<"\n\n";
          // }
          return std::pair<Parameter,std::vector<Atom*> > (param,sorted_atoms);
     };

     // int get_param_id(Atom *atom) {
     //      int biotype = get_biotype(atom);
     //      unsigned int id = biotype_info->at(biotype).param_id;
     //      return id;
     // };
};


//! OPLS out-of-plane bending for trivalent atoms (improper torsion) energy term
class TermOplsImptor: public EnergyTermCommon<TermOplsImptor, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsImptor, ChainFB> EnergyTermCommon;

     //! Parameters
     ImptorParameters parameters;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsImptor(ChainFB *chain, const Settings &settings=Settings(),
                    RandomNumberEngine *random_number_engine = &random_global) 
          : EnergyTermCommon(chain, "opls-imptor", settings, random_number_engine) {

          // Annotate chain with biotype
          parameters.annotate_chain(chain);          
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsImptor(const TermOplsImptor &other,  
                    RandomNumberEngine *random_number_engine,
                    int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters),
            counter(other.counter) {

          // Annotate new chain with biotype
          parameters.annotate_chain(chain);          
     }     

     //! Evaluate imptor energy for a atom
     //! \param atom3 central atom
     //! \return improper torsional or out-of-plane bending energy
     inline double calc_imptor_energy(Atom *atom3) {
          //     std::cout << atom3->covalent_neighbours.size() << " " << parameters.get_bond_n(atom3) << "\n";
          //     assert(atom3->covalent_neighbours.size() == 3);

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // OPLS uses potonated histidines in which case there is a imptor at ND1 and
          // NE2. Phaistos has depotonated hisidines and no imptor at ND1
          // Should probably check all titratable hydrogens before continuing in case of native initialization
          // with hydrogens
          // if(atom3->atom_type == ND1 && atom3->residue->residue_type == HIS)
          if(atom3->atom_type == NE2 && atom3->residue->residue_type == HIS)
               return 0.0;
     
          CovalentBondIterator<ChainFB> it1(atom3, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
          int i=0;
          Atom *neighbour[3];
          for (; !it1.end(); ++it1) {
               neighbour[i] = &*it1;
               i += 1;
          }
          
          std::pair<ImptorParameters::Parameter, std::vector<Atom*> > ogly =
               parameters.get_param_and_sorted_atoms(neighbour[0],neighbour[1],atom3,neighbour[2]);
          ImptorParameters::Parameter param = ogly.first;
          std::vector<Atom*> sorted_atoms = ogly.second;
          double angle = calc_dihedral(sorted_atoms[0]->position,sorted_atoms[1]->position,
                                       atom3->position,sorted_atoms[3]->position);
          double e = calc_spring_energy(angle,param);
          int p1 = parameters.get_param_id(sorted_atoms[0]);
          int p2 = parameters.get_param_id(sorted_atoms[1]);
          /* int p3 = parameters.get_param_id(atom3); */
          int p4 = parameters.get_param_id(sorted_atoms[3]);
          
          if (p1 == p2 && p2 == p4) {
               // six-fold symmetri
               double energies[6];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/6.0); */
               angle = calc_dihedral(sorted_atoms[1]->position,sorted_atoms[0]->position,
                                     sorted_atoms[2]->position,sorted_atoms[3]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[1]/6.0); */
               angle = calc_dihedral(sorted_atoms[1]->position,sorted_atoms[3]->position,
                                     sorted_atoms[2]->position,sorted_atoms[0]->position);
               energies[2] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[2]/6.0); */
               angle = calc_dihedral(sorted_atoms[3]->position,sorted_atoms[1]->position,
                                     sorted_atoms[2]->position,sorted_atoms[0]->position);
               energies[3] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[3]/6.0); */
               angle = calc_dihedral(sorted_atoms[3]->position,sorted_atoms[0]->position,
                                     sorted_atoms[2]->position,sorted_atoms[1]->position);
               energies[4] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[4]/6.0); */
               angle = calc_dihedral(sorted_atoms[0]->position,sorted_atoms[3]->position,
                                     sorted_atoms[2]->position,sorted_atoms[1]->position);
               energies[5] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[5]/6.0); */
               for (i=1; i<6; i++)
                    e += energies[i];
               e /= 6.0;
          } else if (p1 == p2) {
               double energies[2];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/2.0); */
               angle = calc_dihedral(sorted_atoms[1]->position,sorted_atoms[0]->position,
                                     sorted_atoms[2]->position,sorted_atoms[3]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p2,p1,p3,p4,angle*180.0/3.14159265,energies[1]/2.0); */
               e = (energies[0]+energies[1])/2.0;
          } else if (p2 == p4) {
               double energies[2];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/2.0); */
               angle = calc_dihedral(sorted_atoms[0]->position,sorted_atoms[3]->position,
                                     sorted_atoms[2]->position,sorted_atoms[1]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p4,p3,p2,angle*180.0/3.14159265,energies[1]/2.0); */
               e = (energies[0]+energies[1])/2.0;
          } else if (p1 == p4) {
               double energies[2];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/2.0); */
               angle = calc_dihedral(sorted_atoms[3]->position,sorted_atoms[1]->position,
                                     sorted_atoms[2]->position,sorted_atoms[0]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p4,p2,p3,p1,angle*180.0/3.14159265,energies[1]/2.0); */
               e = (energies[0]+energies[1])/2.0;
          } else {
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e); */
          }          
          return e;
     };
     
     //! Evaluate a single imptor term
     //! \param angle Dihedral angle
     //! \param param Parameter container
     inline double calc_spring_energy(double angle, ImptorParameters::Parameter param) {
          double energy=0.0;
          counter++;
          for (unsigned int i=0; i<param.size(); i++) {
               if (param.o[i] <= 0)
                    break;
               //cos is symmetric around angle so we need not flip sign when param is mirrored
               energy += param.amp[i]*( 1 + cos(param.o[i]*angle - param.ang[i]) );
          }
          return energy;
     }

     //! Evaluate
     //! \param move_info object containing information about last move
     //! \return improper torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
          counter = 0;
          double energy_sum = 0.0;
          // all atoms
          int size = (this->chain)->size();
          for (int r = 0; r < size; r++) {
               Residue *res = &(*(this->chain))[r];
               int res_size = res->size();
               for (int a = 0; a < res_size; a++) {
                    Atom *atom = res->atoms[a];
                    // Atom must have excatly 3 neighbours (sp2 hybridized) to have
                    //  a improper torsion 
                    if(parameters.get_bond_n(atom)==3)
                         energy_sum += calc_imptor_energy(atom);
               }
          }
          /* std::cout<<"Improper torsion: "<<energy_sum<<" kcal/mol  "<<counter<<" interactions\n"; */
          return energy_sum;
     };
};
     
}

#endif
