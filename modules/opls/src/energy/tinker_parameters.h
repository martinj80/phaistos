// tinkerParameters.h --- Tinker parameter reader base class
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

#ifndef TINKERPARAMETERS_H
#define TINKERPARAMETERS_H

#include <sstream>
#include <iostream>
#include <string>
#include <cmath>

#include "protein/chain_fb.h"

namespace phaistos {

//! Abstract Tinker Parameter Reader and Container Base Class
//! Derived classes gets functionality to read parameter files in the format
//! provided by the TINKER molecular modeling suite
//! J. W. Ponder et al "Current Status of the AMOEBA Polarizable Force Field",
//! J. Phys. Chem. B 114, page 2549–2564, 2010
class TinkerParameterBase {

private:

     //! forward declaration
     class BiotypeInfo;
     
public:

     //! Abscract Constructor
     TinkerParameterBase() {};

     //! Destructor
     virtual ~TinkerParameterBase() {};

     //! Annotate chain with biotype
     //! biotype=-1 is uninitialized and biotype=0 is biotype undefined (init attempted but failed)
     //! \param chain pointer to the chain that should be annotated
     template <typename CHAIN_TYPE>     
     void annotate_chain(CHAIN_TYPE *chain) {
          for (AtomIterator<CHAIN_TYPE, definitions::ALL> it(*chain); !it.end(); ++it) {
               annotate_atom(&*it);
          }
     }

     //! Annotate atom with biotyoe
     //! \param atom pointer to the atom that should be annotated
     void annotate_atom(Atom *atom) {
          // assign biotype if uninitialized
          if (atom->biotype < 0) {
               atom->biotype = find_biotype(atom);
          }
          // check if assignment was successful (>0)
          assert(atom->biotype > 0);
     }
     
     //! Method to get the biotype of an biotype-assigned phaists Atom instance
     //! \param atom pointer to phaistos Atom intance
     inline int get_biotype(Atom *atom) {
          // atoms should be successfully annotated with biotype by now
          assert(atom->biotype > 0);
          return atom->biotype;
     }
     
     //! Get the tinker param_id of a biotype
     //! Tinker param_id is the second number of an atom line used for all but charge parameters
     inline int get_param_id(int biotype) {
          return ( biotype_info[biotype].param_id );
     }

     //! Get the tinker param_id of an biotype assigned phaists Atom instance
     //! Tinker param_id is the second number of an atom line used for all but charge parameters
     inline int get_param_id(Atom *atom) {
          int biotype = get_biotype(atom);
          return ( get_param_id(biotype) );
     }
     
     //! Get the tinker atom_id of a biotype
     //! Tinker atom_id is the first number of an atom line used for charge parameters
     inline int get_atom_id(int biotype) {
          return ( biotype_info[biotype].atom_id );
     }

     //! Get the tinker atom_id of an biotype assigned phaists Atom instance
     //! Tinker atom_id is the first number of an atom line used for charge parameters
     inline int get_atom_id(Atom *atom) {
          int biotype = get_biotype(atom);
          return ( get_atom_id(biotype) );
     }
     
     //! Get number of chemical bonds (i.e. number of covalent neighbours) of a biotype
     //! db. and triple bonds count as one
     inline int get_bond_n(int biotype) {
          return ( biotype_info[biotype].n_bonds );
     }

     //! Get number of chemical bonds (i.e. number of covalent neighbours) of an atom
     //! db. and triple bonds count as one
     inline int get_bond_n(Atom *atom) {
          int biotype = get_biotype(atom);
          return ( get_bond_n(biotype) );
     }
     
     //! Get atomic number (i.e. number of protons) of a biotype
     inline int get_atomic_number(int biotype) {
          return ( biotype_info[biotype].atomic_number );
     }

     //! Get atomic number (i.e. number of protons) of an atom
     //! \param atom pointer to phaistos atom instance
     inline int get_atomic_number(Atom *atom) {
          int biotype = get_biotype(atom);
          return ( get_atomic_number(biotype) );
     }
     
     //! Translate phaistos atomtype to tinker types (thouse in biotype_info)
     inline std::string get_tinker_atom_type(int biotype) {
          return ( biotype_info[biotype].atom_type );
     }

     //! Translate phaistos atomtype to tinker types (thouse in biotype_info)
     //! \param atom pointer to phaistos atom instance
     inline std::string get_tinker_atom_type(Atom *atom) {
          int biotype = get_biotype(atom);
          return ( get_tinker_atom_type(biotype) );
     }

     //! Convert atom id to param id. Returns zero if atom id dosn't represent a biotype
     //! \param atom_id tinker atom id - first number in a atom line used for charge parameters
     //! \return param_id tinker param id - second number in a atom line used for all but charge param
     inline int atom_id2param_id(int atom_id) {
          unsigned int size=biotype_info.size();
          for (unsigned int i=1; i<size; i++) {
               if (biotype_info[i].atom_id == atom_id)
                    return (biotype_info[i].param_id);
          }
          return 0;
     }

     //! Covalent distance cache. A bit faster than the one inphaistos
     template <typename CHAIN_TYPE>
     int get_covalent_dist(Atom *atom1, Atom *atom2) {
          unsigned int ri1 = (unsigned int) atom1->residue->index;
          unsigned int ai1 = (unsigned int) atom1->index;
          unsigned int ri2 = (unsigned int) atom2->residue->index;
          unsigned int ai2 = (unsigned int) atom2->index;

          //expand cache if necessary
          if (covalent_dist_cache.size() < ri1+1) {
               std::vector<std::vector<std::vector<short int> > > vec3;
               covalent_dist_cache.resize(ri1+1,vec3);
          }
          if (covalent_dist_cache[ri1].size() < ai1+1) {
               std::vector<std::vector<short int> > vec2;
               covalent_dist_cache[ri1].resize(ai1+1,vec2);
          }
          if (covalent_dist_cache[ri1][ai1].size() < ri2+1) {
               std::vector<short int> vec1;
               covalent_dist_cache[ri1][ai1].resize(ri2+1,vec1);
          }
          if (covalent_dist_cache[ri1][ai1][ri2].size() < ai2+1)
               covalent_dist_cache[ri1][ai1][ri2].resize(ai2+1,-1);
     
          //get new chain distance and cache it
          if (covalent_dist_cache[ri1][ai1][ri2][ai2] < 0) {
               int d = chain_distance<CHAIN_TYPE>(atom1,atom2);
               covalent_dist_cache[ri1][ai1][ri2][ai2] = (short int) d;
               return d;
          } else {
               // ...or return the cached one
               return (int) covalent_dist_cache[ri1][ai1][ri2][ai2];
          }
     }

     //! Utill - tinker often use degrees
     //! \param deg the angle in degrees to be converted
     inline double deg2rad(double deg) {
          return ((deg/180.)*M_PI);
     }
     
protected:
     
     //! Vector containing the info to convert and interpret biotypes
     //! This is automaticly assigned when a field of parameters is read using read_param
     //! For different derived tinker parameter classes using same input file this vector will be identical
     //! Not much to do about this when static's are not allowed in phaistos
     //! Public access through public getters
     std::vector<BiotypeInfo> biotype_info;

     //! Method to read a field, i.e. a group of lines with the same keyword, in a tinker parameter file
     //! \param parameter_file_content Content of a tinker parameter file
     //! \param raw_param Pointer to empty vector of (id,parameter) pairs to be filled with the content of the field
     //! \param target_keyword Keyword defining the field
     //! \param number_of_ids Number of atoms and hence ids needed to define the parameters in the field
     void read_param(const std::string *parameter_file_content,
                     std::vector<std::pair<std::vector<int>, std::vector<double> > > *raw_param,
                     std::string target_keyword,
                     int number_of_ids) {

          // For convenience, extract biotype info automaticly when a field of parameters is read
          if ( biotype_info.size() == 0 )
               read_biotype_info(parameter_file_content);

          int i,id;
          double param_val;
          std::pair<std::vector<int>,std::vector<double> > param;
          std::string keyword,line;
          std::istringstream line_stream,content(*parameter_file_content);
          
          while (getline(content,line)) {
               line_stream.clear();     //reset stream error flags
               keyword="";             //reset keyword in case of blank line
               line_stream.str(line);   //use line as stringBuf for line_stream
               
               line_stream >> keyword;
               if (keyword == target_keyword) {
                    for (i=0; i<number_of_ids; i++) {
                         line_stream >> id;
                         param.first.push_back(id);
                    }
                    while (line_stream >> param_val) {
                         param.second.push_back(param_val);
                    }
                    raw_param->push_back(param);
                    param.first.clear();
                    param.second.clear();
               }
          }
     }
     
     //! Method to read a header field in a tinker parameter file
     //! Note that only the first max_lines are read
     //! \param parameter_file_content Content of a tinker parameter file
     //! \param target_keyword header line to look for
     //! \return First word after the keyword
     std::string read_header(const std::string *parameter_file_content, std::string target_keyword) {

          int max_lines=50;
          std::istringstream line_stream,content(*parameter_file_content);
          std::string keyword,line;
          std::string ret="";
          
          for (int i=0; i<max_lines; i++) {
               if ( !getline(content,line) )
                    break;
          
               line_stream.clear();     //reset stream error flags
               keyword="";             //reset keyword in case of blank line
               line_stream.str(line);   //use line as string buffer for line_stream
               
               line_stream >> keyword;
               if (keyword == target_keyword) {
                    line_stream >> ret;
                    break;
               }
          }
          return ret;
     }
     
private:

     //! Container for chain distances     
     std::vector<std::vector<std::vector<std::vector<short int> > > > covalent_dist_cache;

     //! Container Class for tinker biotype info
     //! One instance of this coresponds to a biotype line in a TINKER forcefield parameter file
     //! together with the info from the corresponding atom line
     class BiotypeInfo {

     public:
          
          //! first (consecutive) id number in an atom line. Used for charge parameters
          int atom_id;
          
          //! second id number in an atom line. Used for all but charge parameters
          int param_id;

          //! Special specifier for the residue type
          //! spec%3    : =1: c-term,  =2: n-term,  =0: internal
          //! HIS spec/3: =0: double protonated (charged), =1: ND protonated, =2: NE protonated
          //! CYS spec/3: =0: SH, =1: SS (S-brigde)
          int spec;

          //! Type of residue the instance belongs to
          definitions::ResidueEnum residue_type;

          //! Tinker atom type. convert to atomEnum with tniker2pdbAtomType
          std::string atom_type;

          //! String in the atom line
          std::string atom_description;
          
          //! String in the biotype line
          std::string residue_description;

          //! Atom valence
          int n_bonds;
          
          //! Number of atom protons
          int atomic_number;
          
          //! Default constructor
          BiotypeInfo()
               : atom_id(0),param_id(0),residue_type(definitions::AA_UNDEF),atom_type(""),
                 atom_description("unknown atom"),residue_description("unknown residue"),
                 n_bonds(0),atomic_number(0) {};
          
          //! Overload output oprator
          friend std::ostream &operator<<(std::ostream &o, BiotypeInfo bm) {
               o<<bm.atom_description<<" ("<<bm.atom_type<<") in "<<bm.residue_description
                <<" with "<<bm.n_bonds<<" neighbours";
               return o;
          }
     };

     //! Container Class for tinker atom info
     //! One instance of this coresponds to a atom line in a tinker forcefield parameter file
     struct AtomInfo {
          
          //! first (consecutive) id number in an atom line. Used for charge parameters
          int atom_id;
          
          //! second id number in an atom line. Used for all but charge parameters
          int param_id;

          //! Number of covalent neigubours
          int n_bonds;
          
          //! Number of atom protons
          int number;

          //! First and second string in the atom line
          std::string dscp1,dscp2;
          
          //! Atomic mass
          double mass;
     };

     //! Determine the biotype of an atom only based on info from the phaistos representation
     //! \param atom pointer to phaistos Atom instance to have it's biotype determined
     //! \return biotype of the atom
     int find_biotype(Atom *atom) {

           //! Import protein definitions (such as residue names)
          using namespace definitions;

          int i=1;
          std::string tinker_atom_type;
          int n_biotypes = biotype_info.size();
          pdb2tinker_atom_type(atom,&tinker_atom_type);
          
          //find residue
          while (biotype_info[i].residue_type != atom->residue->residue_type && i<=n_biotypes)
               i++;
          //find backbone atom
          if (!atom->is_sidechain_atom) { //only N-C-C BB in phaistos. This also includes HN, HA and O.
               //if terminal forward to terminal backbone atoms
               if (atom->residue->terminal_status == NTERM) {
                    while ((biotype_info[i].residue_type != atom->residue->residue_type ||
                            biotype_info[i].spec%3 != 2) && i<n_biotypes)
                         i++;
               }
               else if (atom->residue->terminal_status == CTERM) {
                    while ((biotype_info[i].residue_type != atom->residue->residue_type ||
                            biotype_info[i].spec%3 != 1) && i<n_biotypes)
                         i++;
               }
               //find backbone atom
               while (biotype_info[i].atom_type != tinker_atom_type && i<n_biotypes)
                    i++;
          }
          //find sidechain atom
          else { //atom->is_sidechain_atom == true
               //i += 5; //discard bb atoms - HA has isSC=1 so more robust without this
               while (biotype_info[i].atom_type != tinker_atom_type && i<n_biotypes)
                    i++;
          } //fi find atom
          
          //special cases:
          if (atom->residue->residue_type == PRO && atom->residue->terminal_status == NTERM && atom->atom_type == CD) {
               i = 410; //the special PRO CD (105) is not used in the N-terminal
          }
          
          //check again if the right biotype is found and set atom biotype
          if (biotype_info[i].atom_type == tinker_atom_type &&
              biotype_info[i].residue_type == atom->residue->residue_type) {
               return i;
          } else {
               // std::cout<<biotype_info[i].atom_type<<"  "<<tinker_atom_type<<"  "
               //          <<biotype_info[i].res<<"  "<<atom->residue->residue_type<<std::endl;
               return 0; //biotype undefined
          }
     }
     
     //! Fill biotype_info from input file
     //! \param parameter_file_content content of a tinker parameter file
     void read_biotype_info(const std::string *parameter_file_content) {

          AtomInfo atom_line;
          std::vector<AtomInfo> atom_line_vector;
          std::string keyword,line;
          std::istringstream line_stream,content(*parameter_file_content);
          
          if (biotype_info.size() != 0)
               biotype_info.clear();
          
          // biotype_info[0] is undef
          biotype_info.push_back(BiotypeInfo());
          
          // atom_line_vector[0] not used
          atom_line_vector.push_back(AtomInfo());
          
          // Read in atom and biotype fields
          while (getline(content,line)) {
               line_stream.clear();     //reset stream error flags
               keyword="";              //reset keyword in case of blank line
               line_stream.str(line);   //use line as stringBuf for line_stream
               
               line_stream >> keyword;
               // Read atom field
               if (keyword == "atom") {
                    unsigned int atom_n;
                    //atom_id (charge param only), param_id (others)
                    line_stream >> atom_n >> atom_line.param_id >> atom_line.dscp1;
                    
                    getline(line_stream,atom_line.dscp2,'"'); //discard whitespace and first quote (")
                    getline(line_stream,atom_line.dscp2,'"'); //extract description and discard quote
                    
                    //extract atomic number, mass and number of bonds
                    line_stream >> atom_line.number >> atom_line.mass >> atom_line.n_bonds;
                    
                    // check if atom_id numbering is consequitive
                    assert(atom_line_vector.size() == atom_n);
                    atom_line_vector.push_back(atom_line);
                    
               // Read biotype field
               } else if (keyword == "biotype") {
                    BiotypeInfo bf;
                    unsigned int bf_n;
                    line_stream >> bf_n >> bf.atom_type;
                    
                    getline(line_stream,bf.residue_description,'"');
                    getline(line_stream,bf.residue_description,'"');
                    biotype_description_interpret(bf.residue_description,&bf.residue_type,&bf.spec);
                    line_stream >> bf.atom_id;
                    
                    // check if biotype numbering is consequitive
                    assert(biotype_info.size() == bf_n);
                    biotype_info.push_back(bf);
               }
          }
          // put info from atom_line_vector in biotype_info
          for (unsigned i=0; i<biotype_info.size(); i++) {
               int atom_id = biotype_info[i].atom_id;
               biotype_info[i].param_id = atom_line_vector[atom_id].param_id;
               biotype_info[i].n_bonds = atom_line_vector[atom_id].n_bonds;
               biotype_info[i].atomic_number = atom_line_vector[atom_id].number;
               biotype_info[i].atom_description = atom_line_vector[atom_id].dscp2;
          }
     }

     //! Util functions to intrepret tinker parameter files - used in read_biotype_info
     //! The TinkerParameterBase class has the same atom type names (stored in biotype_info) 
     //! as the tinker parameter files. The names in tinker contain less information than 
     //! the pdb names (which e.g. have index on chemical identical hydrogen atoms) and 
     //! hence the obvious is to convert pdb to tinker atom types.
     //! Differences are commented in the following function.
     //! \param atom atom to be translated
     //! \param tinker_atom_type tinker atom name to be returned
     void pdb2tinker_atom_type(Atom *atom, std::string *tinker_atom_type) {

           //! Import protein definitions (such as residue names)
          using namespace definitions;
          
          std::string phaistos_atom_type = atom_name[atom->atom_type];
          ResidueEnum phaistos_res_type = atom->residue->residue_type;
          // fix hydrogen types
          if (phaistos_atom_type[0] == 'H') {
               if (phaistos_atom_type.length() == 1)
                    *tinker_atom_type = "HN"; //bb amide hydrogen
               else if (phaistos_atom_type[1] == '1' || phaistos_atom_type[1] == '2' || phaistos_atom_type[1] == '3')
                    *tinker_atom_type = "HN"; //N-term ammonium hydrogen
               else if (phaistos_atom_type.length() > 2) { //sc hydrogen
                    tinker_atom_type->assign(phaistos_atom_type,0,2);           //remove atom index
                    if (phaistos_res_type == VAL && phaistos_atom_type[1] == 'G') //and put back if nessesary
                         *tinker_atom_type += phaistos_atom_type[2];               //second index is never used
                    else if (phaistos_res_type == LEU && phaistos_atom_type[1] == 'D')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == ILE && phaistos_atom_type[1] == 'G')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == THR && phaistos_atom_type[1] == 'G')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == TRP && phaistos_atom_type[1] != 'B')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == HIS && phaistos_atom_type[1] != 'B')
                    *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == ASN && phaistos_atom_type[1] == 'D')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == GLN && phaistos_atom_type[1] == 'E')
                         *tinker_atom_type += phaistos_atom_type[2];
               }
               else //phaistos_atom_type.length() is 2 and should be ok
                    *tinker_atom_type = phaistos_atom_type;
          }
          //fix carbon
          else if (phaistos_atom_type[0] == 'C') {
               if (phaistos_atom_type.length() > 2) { //sc carbon with index
                    tinker_atom_type->assign(phaistos_atom_type,0,2);           //remove atom index
                    if (phaistos_res_type == VAL && phaistos_atom_type[1] == 'G') //and put back if nessesary
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == LEU && phaistos_atom_type[1] == 'D')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == ILE && phaistos_atom_type[1] == 'G')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == THR && phaistos_atom_type[1] == 'G')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == TRP && phaistos_atom_type[1] != 'B' && phaistos_atom_type[1] != 'G')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == HIS && (phaistos_atom_type[1] == 'D' || phaistos_atom_type[1] == 'E'))
                         *tinker_atom_type += phaistos_atom_type[2];
               }
               else //phaistos_atom_type.length() is <=2 and should be ok
                    *tinker_atom_type = phaistos_atom_type;
          }
          //fix nitrogen
          else if (phaistos_atom_type[0] == 'N') {
               if (phaistos_atom_type.length() > 2) { //sc oxygen with index
                    tinker_atom_type->assign(phaistos_atom_type,0,2);           //remove atom index
                    if (phaistos_res_type == ASN && phaistos_atom_type[1] == 'D') //and put back if nessesary
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == TRP)
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == HIS)
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == GLN && phaistos_atom_type[1] == 'E')
                         *tinker_atom_type += phaistos_atom_type[2];
               }
               else //phaistos_atom_type.length() is <=2 and should be ok
               *tinker_atom_type = phaistos_atom_type;
          }
          //fix oxygen
          else if (phaistos_atom_type[0] == 'O') {
               if (phaistos_atom_type == "OXT")
                    *tinker_atom_type = phaistos_atom_type;
               else if (phaistos_atom_type.length() == 1 && atom->residue->terminal_status  == CTERM)
                    *tinker_atom_type = "OXT"; //both C-term oxygens are called OXT in tinker
               else if (phaistos_atom_type.length() > 2) { //sc oxygen with index
                    tinker_atom_type->assign(phaistos_atom_type,0,2);           //remove atom index
                    if (phaistos_res_type == THR && phaistos_atom_type[1] == 'G') //and put back if nessesary
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == ASN && phaistos_atom_type[1] == 'D')
                         *tinker_atom_type += phaistos_atom_type[2];
                    else if (phaistos_res_type == GLN && phaistos_atom_type[1] == 'E')
                         *tinker_atom_type += phaistos_atom_type[2];
               }
               else //phaistos_atom_type.length() is <=2 and should be ok
                    *tinker_atom_type = phaistos_atom_type;
          }
          else { // sulfur need not be fixed
               *tinker_atom_type = phaistos_atom_type;
          }
     }

     //! Method to parse the residue description of a biotype line to a phaistos residue
     //! \param description string from a biotype line
     //! \param residue phaistos residue type to be returned
     //! \param spec special residue specifier to be returned
     void biotype_description_interpret(std::string description, definitions::ResidueEnum *residue, int *spec) {

           //! Import protein definitions (such as residue names)
          using namespace definitions;
          
          std::string line="";
          
          //make upcase
          std::string::iterator i=description.begin();
          for (; i != description.end(); ++i)
               line += ::toupper((unsigned char)*i);
          
          //find terminal
          if (line.find("C-TERMINAL") != std::string::npos)
               *spec = 1;
          else if (line.find("N-TERMINAL") != std::string::npos)
               *spec = 2;
          else 
               *spec=0;
          
          //find residue
          if ( line.find("HIS") != std::string::npos ) {
               *residue=HIS;
               if ( line.find("HD") != std::string::npos) //default is protonated (+)
                    *spec += 3;
               else if ( line.find("HE") != std::string::npos)
                    *spec += 6;
          }
          else if ( line.find("CYS") != std::string::npos ) {
               *residue=CYS;
               if (line.find("SS") != std::string::npos)
                    *spec += 3;
          }
          else if (line=="GLYCINE" || line.find(" GLY")!=std::string::npos)
               *residue=GLY;
          else if (line=="ALANINE"    || line.find(" ALA")!=std::string::npos)
               *residue=ALA;
          else if (line=="VALINE"     || line.find(" VAL")!=std::string::npos)
               *residue=VAL;
          else if (line=="LEUCINE"    || line.find(" LEU")!=std::string::npos)
               *residue=LEU;
          else if (line=="ISOLEUCINE" || line.find(" ILE")!=std::string::npos)
               *residue=ILE;
          else if (line=="SERINE"     || line.find(" SER")!=std::string::npos)
               *residue=SER;
          else if (line=="THREONINE"  || line.find(" THR")!=std::string::npos)
               *residue=THR;
          else if (line=="PROLINE"       || line.find(" PRO")!=std::string::npos)
               *residue=PRO;
          else if (line=="PHENYLALANINE" || line.find(" PHE")!=std::string::npos)
               *residue=PHE;
          else if (line=="TYROSINE"      || line.find(" TYR")!=std::string::npos)
               *residue=TYR;
          else if (line=="TRYPTOPHAN"    || line.find(" TRP")!=std::string::npos)
               *residue=TRP;
          else if (line=="ASPARTIC ACID" || line.find(" ASP")!=std::string::npos)
               *residue=ASP;
          else if (line=="ASPARAGINE"    || line.find(" ASN")!=std::string::npos)
               *residue=ASN;
          else if (line=="GLUTAMIC ACID" || line.find(" GLU")!=std::string::npos)
               *residue=GLU;
          else if (line=="GLUTAMINE"     || line.find(" GLN")!=std::string::npos)
               *residue=GLN;
          else if (line=="METHIONINE"    || line.find(" MET")!=std::string::npos)
               *residue=MET;
          else if (line=="LYSINE"        || line.find(" LYS")!=std::string::npos)
               *residue=LYS;
          else if (line=="ARGININE"      || line.find(" ARG")!=std::string::npos)
               *residue=ARG;
          else
               *residue = AA_UNDEF;
     }
     
};
     
}

#endif
