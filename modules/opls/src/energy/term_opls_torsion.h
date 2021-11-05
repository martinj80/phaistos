// opls_torsion.h --- OPLS torsion angle energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson
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

#ifndef OPLS_TORSION_H
#define OPLS_TORSION_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "energy/tinker_parameters.h"
#include "energy/parameter_data_structures.h"
#include "energy/opls_parameters.h"

namespace phaistos {

//! Parameter container for OPLS torsion energy term
class DihedralParameters: public TinkerParameterBase {

public:

     //! Dihedral angle parameter container
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
          operator bool () {
               return (o[0] >= 0); // return false only if default constructed
          }
          
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
          }
    
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
     DihedralParameters() {

          // get torsion unit from parameter file header
          std::string torsionStr = read_header(&parameters_opls,"torsionunit");
          torsion_unit = atof( torsionStr.c_str() );

          // get parameters from 'angle' field
          std::vector<std::pair<std::vector<int>, std::vector<double> > > rawParam;
          read_param(&parameters_opls,&rawParam,"torsion",4);

          // store parameters
          for (unsigned int i=0; i<rawParam.size(); i++) {
               Parameter param;
               std::vector<int> id = (rawParam[i].first);
               std::vector<double> paramVal = rawParam[i].second;
               unsigned int n_param = paramVal.size();
               if (n_param == 0) {
                    // make empty parameter
                    param.assign(0,0.0,0.0);
               } else {
                    for (unsigned int k=0; k<n_param; k+=3) {
                         paramVal[k+1] = deg2rad(paramVal[k+1]);
                         param.assign( (int)paramVal[k+2], paramVal[k], paramVal[k+1] );
                    }
               }
               // Check for previous definitions of this parameter
               assert( ! parameter_map(id[0],id[1],id[2],id[3]) );
               // if ( parameter_map(id[0],id[1],id[2],id[3]) ) {
               //      std::cerr<<"\nOPLS TERM TORSION - multiple definitions of parameter "
               //               <<id[0]<<","<<id[1]<<","<<id[2]<<","<<id[3]<<"\n\n";
               // }
               parameter_map(id[0],id[1],id[2],id[3]) = param;
          }
     };

     //! Destructor
     ~DihedralParameters() {};

     //! Getter
     Parameter get(Atom *atom1, Atom *atom2, Atom *atom3, Atom *atom4) {
          int id1 = get_param_id(atom1);
          int id2 = get_param_id(atom2);
          int id3 = get_param_id(atom3);
          int id4 = get_param_id(atom4);
          if (id2 > id3) {
               int buf=id2; id2=id3; id3=buf;
               buf=id1; id1=id4; id4=buf;
          } else if (id2==id3 && id1>id4) {
               int buf=id1; id1=id4; id4=buf;
          }
          Parameter param = parameter_map(id1,id2,id3,id4);
          assert( param );
          // if (!param) {
          //      std::cerr<<"\nOPLS TERM TORSION - missing parameter for param id "
          //               <<id1<<","<<id2<<","<<id3<<","<<id4<<"\n\n";
          // }
          return param;
     };

     //! Overall scaling factor given in parameter file header
     double torsion_unit;
};


//! OPLS torion energy term
class TermOplsTorsion: public EnergyTermCommon<TermOplsTorsion, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsTorsion, ChainFB> EnergyTermCommon;

     //! Parameters
     DihedralParameters parameters;

     //! Number of interactions calculated
     int counter;
     
public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsTorsion(ChainFB *chain, 
                     const Settings &settings=Settings(),
                     RandomNumberEngine *random_number_engine = &random_global) 
          : EnergyTermCommon(chain, "opls-torsion", settings, random_number_engine) {

          // Annotate chain with biotype
          parameters.annotate_chain(chain);          
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsTorsion(const TermOplsTorsion &other, 
                     RandomNumberEngine *random_number_engine,
                     int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            parameters(other.parameters),
            counter(other.counter) {
          
          // Annotate new chain with biotype
          parameters.annotate_chain(chain);          
     }

     //! Evaluate torsion energy for a bond (atom2-atom3)
     //! \param atom2 First atom defining the bond
     //! \param atom3 Second atom defining the bond
     //! \return Torsional energy of that bond
     inline double calc_torsion_energy(Atom *atom2, Atom *atom3) {
          double angle,energy = 0.0;
          CovalentBondIterator<ChainFB> it1(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
          for (; !it1.end(); ++it1) {
               Atom *atom1 = &*it1;
               if (atom1 == atom3)
                    continue;
               CovalentBondIterator<ChainFB> it4(atom3, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
               for (; !it4.end(); ++it4) {
                    Atom *atom4 = &*it4;
                    if (atom4 == atom2)
                         continue;
                    DihedralParameters::Parameter param = parameters.get(atom1,atom2,atom3,atom4);
                    angle = calc_dihedral(atom1->position,atom2->position,
                                          atom3->position,atom4->position);
                    energy += calc_spring_energy(angle,param);
               }
          }
          return energy;
     }
     
     //! Evaluate a single torsional term
     //! \param angle Dihedral angle
     //! \param param Parameter container
     inline double calc_spring_energy(double angle, DihedralParameters::Parameter param) {
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

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
          
          double energy_sum=0.0;
          counter=0;

          // all torsions
          int size = (this->chain)->size();
          for (int r=0; r<size; r++) {
               Residue *res = &(*(this->chain))[r];
               int res_size = res->size();
               for (int a=0; a<res_size; a++) {
                    Atom *atom2 = res->atoms[a];  
                    CovalentBondIterator<ChainFB> it(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
                    for (; !it.end(); ++it) {
                         Atom *atom3 = &*it;
                         if(atom2->index < atom3->index)
                              energy_sum += calc_torsion_energy(atom2,atom3);
                    }
               }
          }
          
          /* // phaistos degrees of freedom only */
          /* DofIterator::angle_selection_enum dofs = DofIterator::DIHEDRAL_DOFS + */
          /*      DofIterator::CHI_ANGLES + DofIterator::N_DIHEDRAL; */
          /* DofIterator dofIt((*chain)(0,N),DofIterator::DIHEDRAL,dofs); */
          /* DofIterator end = chain->dofIteratorEnd(); */
          /* for (; dofIt!=end; ++dofIt) { */
          /*      Atom *atom3 = dofIt.getAtom(); */
          /*      Atom *atom1init,*atom2,*atom4init; */
          /*      atom3->get_dihedral_atoms(&atom1init,&atom2,&atom3,&atom4init); */
          /*      energy_sum += calc_torsion_energy(atom2,atom3); */
          /* } */

          energy_sum *= parameters.torsion_unit;

          /* std::cout<<"Torsional Angle "<<energy_sum<<" kcal/mol "<<counter<<" interactions\n"; */

          return energy_sum;
     }
};

}

#endif
