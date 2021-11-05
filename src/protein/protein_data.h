// atomData.h --- Protein input data
// Copyright (C) 2006-2008 Wouter Boomsma
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


#ifndef ATOMDATA_H
#define ATOMDATA_H

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include "utils/vector_matrix_3d.h"
#include "protein/definitions.h"
#include "utils/math.h"

namespace phaistos {

//! Input data for single atom
class AtomData {
public:

     //! Residue index (ResidueEnum)
     std::string *residue;

     //! Resseq index
     int *resseq;

     //! Atom string
     std::string *atom;

     //! Altloc string
     std::string *altloc;

     //! Secondary structure index (real)
     int *real_ss;

     //! Secondary structure index (predicted)
     int *pred_ss;

     //! Secondary structure propensity - C
     double *pred_prob_ss_C;

     //! Secondary structure propensity - H
     double *pred_prob_ss_H;

     //! Secondary structure propensity - E
     double *pred_prob_ss_E;
     
     //! 3D coordinates
     double *coordinates;

     //! Initialize (all member variables are optional)
     void init() {
          residue = NULL;
          resseq = NULL;
          atom = NULL;
          altloc = NULL;
          real_ss = NULL;
          pred_ss = NULL;
          pred_prob_ss_C = NULL;
          pred_prob_ss_H = NULL;
          pred_prob_ss_E = NULL;
          coordinates = NULL;
     }

     //! Constructor
     AtomData(){
          init();
     }

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     AtomData(const AtomData &other) {
          this->init();
          if (other.residue)
               add_residue(*other.residue);
          if (other.resseq)
               add_resseq(*other.resseq);
          if (other.atom)
               add_atom(*other.atom);
          if (other.altloc)
               add_altloc(*other.altloc);
          if (other.coordinates)
               add_coordinates(other.coordinates[0], other.coordinates[1], other.coordinates[2]);
          if (other.real_ss)
               add_real_ss(*other.real_ss);
          if (other.pred_ss)
               add_pred_ss(*other.pred_ss);
          if (other.pred_prob_ss_C)
               add_pred_prob_ss_C(*other.pred_prob_ss_C);
          if (other.pred_prob_ss_H)
               add_pred_prob_ss_H(*other.pred_prob_ss_H);
          if (other.pred_prob_ss_E)
               add_pred_prob_ss_E(*other.pred_prob_ss_E);
     }

     //! Destructor
     ~AtomData() {
          delete residue;
          delete resseq;
          delete atom;
          delete altloc;
          delete real_ss;
          delete pred_ss;
          delete pred_prob_ss_C;
          delete pred_prob_ss_H;
          delete pred_prob_ss_E;
          delete[] coordinates;
     }

     //! Set residue value
     //!
     //! \param residue Residue index (ResidueEnum)
     void add_residue(std::string residue) {
          this->residue = new std::string(residue);
     }

     //! Set resseq value
     //!
     //! \param resseq Resseq index
     void add_resseq(int resseq) {
          this->resseq = new int(resseq);
     }

     //! Set atom value
     //!
     //! \param atom Atom string
     void add_atom(std::string atom) {
          this->atom = new std::string(atom);
     }

     //! Set altloc value
     //!
     //! \param altloc Altloc identification string
     void add_altloc(std::string altloc) {
          std::string modified_altloc = altloc;
          if (altloc == std::string(" "))
               modified_altloc = "_";
          this->altloc = new std::string(modified_altloc);
     }

     //! Set 3D-coordinate value
     //!
     //! \param x x-coordinate     
     //! \param y y-coordinate     
     //! \param z z-coordinate     
     void add_coordinates(double x, double y, double z) {
          this->coordinates = new double[3];
          this->coordinates[0] = x;
          this->coordinates[1] = y;
          this->coordinates[2] = z;
     }

     //! Set secondary structure value (real)
     //!
     //! \param ss Secondary structure index (SecondaryStructureEnum)
     void add_real_ss(int ss) {
          this->real_ss = new int(ss);
     }

     //! Set secondary structure value (predicted)
     //!
     //! \param ss Secondary structure index (SecondaryStructureEnum)
     void add_pred_ss(int ss) {
          this->pred_ss = new int(ss);
     }

     //! Set secondary structure C propersity value
     //!
     //! \param ss_prob Secondary structure propensity
     void add_pred_prob_ss_C(double ss_prob) {
          this->pred_prob_ss_C = new double(ss_prob);
     }

     //! Set secondary structure H propensity value
     //!
     //! \param ss_prob Secondary structure propensity
     void add_pred_prob_ss_H(double ss_prob) {
          this->pred_prob_ss_H = new double(ss_prob);
     }

     //! Set secondary structure E propensity value
     //!
     //! \param ss_prob Secondary structure propensity
     void add_pred_prob_ss_E(double ss_prob) {
          this->pred_prob_ss_E = new double(ss_prob);
     }
     
     //! Overload input operator for AtomData
     friend std::ifstream &operator>>(std::ifstream &in, AtomData &data) {
          std::vector<std::string> words;
          if (in.good()) {
               std::string line;
               if (getline(in, line)) {
                    boost::trim(line);
                    std::vector<std::string> words;
                    boost::split(words, line, boost::is_any_of("\t "), boost::token_compress_on);
               
                    for (unsigned int i=0; i<words.size(); i++) {
                         std::string word = words[i];
                         if (word == "RES") {
                              data.add_residue(words[++i].data());
                         } else if (word == "RESSEQ") {
                              data.add_resseq(atoi(words[++i].c_str()));
                         } else if (word == "ATOM") {
                              data.add_atom(words[++i]);
                         } else if (word == "ALTLOC") {
                              data.add_altloc(words[++i]);
                         } else if (word == "COORD") {
                              double x = atof(words[++i].data());
                              double y = atof(words[++i].data());
                              double z = atof(words[++i].data());
                              data.add_coordinates(x, y, z);
                         } else if (word == "REALSS") {
                              data.add_real_ss(atoi(words[++i].data()));
                         } else if (word == "PREDSS") {
                              data.add_pred_ss(atoi(words[++i].data()));
                         } else if (word == "PREDPROBSS_C") {
                              data.add_pred_prob_ss_C(atof(words[++i].data()));
                         } else if (word == "PREDPROBSS_H") {
                              data.add_pred_prob_ss_H(atof(words[++i].data()));
                         } else if (word == "PREDPROBSS_E") {
                              data.add_pred_prob_ss_E(atof(words[++i].data()));
                         }
                    }
               }
          }
          return in;
     }

     //! Overload output operator for AtomData
     friend std::ostream & operator<<(std::ostream &o, AtomData &data) {
          char output[200];
          if (data.residue) {
               sprintf(output, "RES %4s ", data.residue->c_str());
               o << output;
          }

          if (data.resseq) {
               sprintf(output, "RESSEQ %4d ", *data.resseq);
               o << output;
          }

          if (data.atom) {
               sprintf(output, "ATOM %3s ", data.atom->data());
               o << output;
          }

          if (data.altloc) {
               sprintf(output, "ALTLOC %1s ", data.altloc->data());
               o << output;
          }

          if (data.coordinates) {
               sprintf(output, "COORD %20.15f %20.15f %20.15f ",
                       data.coordinates[0], data.coordinates[1], data.coordinates[2]);
               o << output;
          }

          if (data.real_ss) {
               sprintf(output, "REALSS %2d ", *data.real_ss);
               o << output;
          }

          if (data.pred_ss) {
               sprintf(output, "PREDSS %2d ", *data.pred_ss);
               o << output;
          }
          
          if (data.pred_prob_ss_C) {
               sprintf(output, "PREDPROBSS_C %20.15f ", *data.pred_prob_ss_C);
               o << output;
          }

          if (data.pred_prob_ss_H) {
               sprintf(output, "PREDPROBSS_H %20.15f ", *data.pred_prob_ss_H);
               o << output;
          }

          if (data.pred_prob_ss_E) {
               sprintf(output, "PREDPROBSS_E %20.15f ", *data.pred_prob_ss_E);
               o << output;
          }
          
          return o;
     }          
};


//! Input data for entire protein 
class ProteinData {
public:
     //! Internal data container.
     //! data: vector<polypeptide>
     //! polypeptide: vector<residue>
     //! residue: vector<atoms>
     std::vector<std::vector<std::vector<AtomData *> > > data;

     //! Constructor (default)
     ProteinData(){};

     //! Constructor
     //!
     //! \param filename Pdb input filename
     ProteinData(std::string filename){

          // Do nothing if filename is empty
          if (filename == "") {
               return;
          }

          // Open file stream
          std::ifstream ifs(filename.c_str());
          if (!ifs) {
               std::cerr<<"Cannot read pdb file " << filename << "\n";
               assert(false);
          }
     
          // Get contents of PDB file
          std::string pdbData = std::string(std::istreambuf_iterator<char>(ifs),
                                            std::istreambuf_iterator<char>());
          ifs.seekg(0);

          // Whether chain index is included in PDB file
          bool include_chain=0;

          std::string current_chain = "";
          std::string prev_chain;
          int chain_index = -1;
     
          std::string current_res = "-1000";
          std::string prev_res;
          int res_index = -1;

          bool only_CA=true;
          bool skip_to_next_chain = false;
          bool skip_to_next_res = false;

          std::string line;
          while (getline(ifs, line)) {

               boost::trim(line);

               // Check for termination 
               if (line.substr(0,3)=="TER") {
                    skip_to_next_chain = true;
               }

               // Check whether line contains atom coordinates
               if (line.substr(0,4)!="ATOM") {
                    continue;
               }
          
               // // Skip weird amino acids
               // if (map_lookup(strToRes,res,int(AA_UNDEF)) == AA_UNDEF) {
               //      continue;
               // }

               // std::string res = strip(line.substr(17,3));
               std::string res = boost::trim_copy(line.substr(17,3));

               prev_res = current_res;
               // current_res = strip(line.substr(22,5)).c_str();
               current_res = boost::trim_copy(line.substr(22,5)).c_str();

               // std::string atom = strip(line.substr(12, 4));
               std::string atom = boost::trim_copy(line.substr(12, 4));
               if (atom != "CA")
                    only_CA = false;
          
               // Check for chain index
               if (line[21]!=' ') {
                    include_chain=1;
                    prev_chain = current_chain;
                    current_chain = line[21];
               } else {
                    prev_chain = current_chain;
                    current_chain = " ";        
                    include_chain=0;
               }

               if ((prev_chain != current_chain)) {
                    this->new_polypeptide();
                    chain_index++;
                    res_index = -1;
                    prev_res = -1;
                    skip_to_next_chain = false;
               } else if (skip_to_next_chain) {
                    continue;
               }

               // Skip if first atom in a residue is not a N, CA, C, O, CB
               if (current_res!=prev_res) {
                    skip_to_next_res = false;
               
                    if (!(atom == "N")) {
                         skip_to_next_res = true;
                         continue;
                    } 
               } else if (skip_to_next_res) {
                    continue;
               }
                  
               if (current_res!=prev_res) {
                    this->new_residue(chain_index);
                    res_index++;
               }

               this->new_atom(chain_index, res_index);

               std::string altloc = line.substr(16, 1);
          

               double coordx = atof(line.substr(30, 8).c_str());
               double coordy = atof(line.substr(38, 8).c_str());
               double coordz = atof(line.substr(46, 8).c_str());

               this->current_atom()->add_atom(atom);
               this->current_atom()->add_altloc(altloc);
               this->current_atom()->add_residue(res);
               this->current_atom()->add_coordinates(coordx,coordy,coordz);
               this->current_atom()->add_resseq(atoi(current_res.c_str()));
               
          }

          // Scan data for chain breaks, and split data if necessary
          this->split_by_chain_breaks(only_CA);
     
     };

     //! Copy Constructor
     //!
     //! \param other Source object from which copy is made.
     ProteinData(const ProteinData &other) {
          for (int i=0; i<other.n_polypeptides(); i++) {
               data.push_back(std::vector<std::vector<AtomData *> >());
               for (int j=0; j<other.n_residues(i); j++) {
                    data[i].push_back(std::vector<AtomData *>());
                    for (int k=0; k<other.n_atoms(i,j); k++) {
                         data[i][j].push_back(new AtomData(*other.data[i][j][k]));
                    }
               }
          }
     }

     //! Destructor
     ~ProteinData() {
          for (int i=0; i<n_polypeptides(); i++) {
               for (int j=0; j<n_residues(i); j++) {
                    for (int k=0; k<n_atoms(i,j); k++) {
                         delete data[i][j][k];
                    }
               }
          }
     }

     //! Overload boolean evaluation operator
     operator bool() const {
          return !data.empty();
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current object (*this)
     ProteinData &operator=(const ProteinData &other) {
          if (this == &other)
               return *this;

          for (int i=0; i<other.n_polypeptides(); i++) {
               data.push_back(std::vector<std::vector<AtomData *> >());
               for (int j=0; j<other.n_residues(i); j++) {
                    data[i].push_back(std::vector<AtomData *>());
                    for (int k=0; k<other.n_atoms(i,j); k++) {
                         data[i][j].push_back(new AtomData(*other.data[i][j][k]));
                    }
               }
          }     
          return *this;
     }
     
     //! Return number of polypeptides
     //!
     //! \return Number of polypeptides in data object
     int n_polypeptides() const {
          return data.size();
     }

     //! Return number of residues in given polypeptide
     //!
     //! \param pp_index Polypeptide index
     //! \return Number of residues in polypeptide
     int n_residues(int pp_index) const {
          return data[pp_index].size();
     }

     //! Return number of atoms in given residue
     //!
     //! \param pp_index Polypeptide index
     //! \param res_index Residue index
     //! \return Number of atoms in residue
     int n_atoms(int pp_index, int res_index) const {
          return data[pp_index][res_index].size();
     }     

     //! Create new polypeptide
     void new_polypeptide() {
          data.push_back(std::vector<std::vector<AtomData *> >());
     }

     //! Create new residue
     //!
     //! \param pp_index Index of polypeptide in which residue should be placed
     void new_residue(int pp_index) {
          while (pp_index >= n_polypeptides())
               new_polypeptide();
          data[pp_index].push_back(std::vector<AtomData *>());
     }

     //! Create new atom
     //!
     //! \param pp_index Index of polypeptide in which atom should be placed
     //! \param res_index Index of residue in which atom should be placed
     void new_atom(int pp_index, int res_index) {
          while (pp_index >= n_polypeptides())
               new_polypeptide();
          while (res_index >= n_residues(pp_index))
               new_residue(pp_index);
          data[pp_index][res_index].push_back(new AtomData());
     }

     //! Return last added atom
     //!
     //! \return Pointer to atom
     AtomData *current_atom() {
          int lastPP = data.size()-1;
          int lastRes = data[lastPP].size()-1;
          int lastAtom = data[lastPP][lastRes].size()-1;
          return data[lastPP][lastRes][lastAtom];
     }

     //! Split polypeptides into smaller pieces if chainbreaks are detected - full atom version
     //!
     //! \param radius Maximum distance between consecutive atoms
     void split_by_chain_breaks_fb(double radius = 1.8) {
          for (int i=0; i<n_polypeptides(); i++) {
               std::vector<AtomData *> list_current_N;
               std::vector<AtomData *> list_prev_N;
               std::vector<AtomData *> list_current_C;
               std::vector<AtomData *> list_prev_C;

               for (int j=0; j<n_residues(i); j++) {
                    list_prev_C = list_current_C;
                    list_current_C.clear();
                    list_prev_N = list_current_N;
                    list_current_N.clear();
                    for (int k=0; k<n_atoms(i,j); k++) {
                         if (*data[i][j][k]->atom == std::string("C"))
                              list_current_C.push_back(data[i][j][k]);
                         else if (*data[i][j][k]->atom == std::string("N")) {
                              list_current_N.push_back(data[i][j][k]);
                         }
                    }

                    if (j>0) {
                         bool bond = false;
                         AtomData *current_N;
                         std::string altloc_N;
                         AtomData *prev_C;
                         std::string altloc_C;
                         for (unsigned int k1=0; k1<list_current_N.size() && !bond; k1++) {
                              current_N = list_current_N[k1];
                              altloc_N = *current_N->altloc;
                              for (unsigned int k2=0; k2<list_prev_C.size() && !bond; k2++) {
                                   prev_C = list_prev_C[k2];
                                   altloc_C = *prev_C->altloc;
                                   // To form a peptide bond, N and C must be 
                                   // within radius and have the same altloc
                                   // identifier or one altloc blanc
                                   if (altloc_N == altloc_C || altloc_N=="_" || altloc_C=="_") {
                                        double distance = sqrt(Math<double>::sqr((current_N->coordinates[0]) - (prev_C->coordinates[0])) +
                                                               Math<double>::sqr((current_N->coordinates[1]) - (prev_C->coordinates[1])) +
                                                               Math<double>::sqr((current_N->coordinates[2]) - (prev_C->coordinates[2])));
                                        if (distance <= radius) {
                                             // Even if several altloc path exist, we just pick the first found.
                                             bond = true;
                                        }
                                   }
                              }
                         }

                         if (bond) {
                              // Remove atoms in previous residue with wrong altLoc
                              if (altloc_C != "_") {
                                   int k1 = n_atoms(i,j-1)-1;
                                   std::string current_altloc = altloc_C;
                                   std::string prev_altloc;
                                   while (k1>=0 and current_altloc != std::string("_")) {
                                        AtomData *current_atom = data[i][j-1][k1];
                                        prev_altloc = current_altloc;
                                        current_altloc = *current_atom->altloc;
                                        if (!(current_altloc == "_" || (current_altloc == prev_altloc))) {
                                             std::vector<AtomData *>::iterator it = data[i][j-1].begin()+k1;
                                             delete data[i][j-1][k1];
                                             data[i][j-1].erase(it);
                                             current_altloc = prev_altloc;
                                        }
                                        k1--;
                                   }
                              }

                              // Remove N atoms in current residue with wrong altLoc
                              if (altloc_N != "_") {
                                   int k1 = 0;
                                   while (k1 < n_atoms(i,j)) {
                                        AtomData *current_atom = data[i][j][k1];
                                        if ((*current_atom->atom == std::string("N")) and current_atom != current_N) {
                                             std::vector<AtomData *>::iterator it = data[i][j].begin()+k1;
                                             delete data[i][j][k1];
                                             data[i][j].erase(it);
                                             list_current_N.clear();
                                             list_current_N.push_back(current_N);
                                        } else {
                                             k1++;
                                        }
                                   }
                              }
                         } else {
                              std::cerr << "WARNING: Chain break at residue " << j << " in chain " << i << ".\n";
                              std::vector<std::vector<std::vector<AtomData *> > >::iterator chain_iterator= data.begin() + i;
                              std::vector<std::vector<AtomData *> >::iterator res_iterator= data[i].begin() + j;
                         

                              std::vector<std::vector<AtomData *> > new_chain(res_iterator, data[i].end());
                              data[i].erase(res_iterator, data[i].end());

                              data.insert(chain_iterator+1, new_chain);
                         }
                    }
               }
          }
     }


     //! Split polypeptides into smaller pieces if chainbreaks are detected - CA-only version
     //!
     //! \param radius Maximum distance between consecutive atoms
     void split_by_chain_breaks_ca(double radius = 4.3) {
          for (int i=0; i<n_polypeptides(); i++) {
               std::vector<AtomData *> list_current_CA;
               std::vector<AtomData *> list_prev_CA;

               for (int j=0; j<n_residues(i); j++) {
                    list_prev_CA = list_current_CA;
                    list_current_CA.clear();
                    for (int k=0; k<n_atoms(i,j); k++) {
                         if (*data[i][j][k]->atom == std::string("CA"))
                              list_current_CA.push_back(data[i][j][k]);
                    }

                    if (j>0) {
                         bool bond = false;
                         AtomData *current_CA;
                         std::string current_altloc;
                         AtomData *prev_CA;
                         std::string prev_altloc;
                         for (unsigned int k1=0; k1<list_current_CA.size() && !bond; k1++) {
                              current_CA = list_current_CA[k1];
                              current_altloc = *current_CA->altloc;
                              for (unsigned int k2=0; k2<list_prev_CA.size() && !bond; k2++) {
                                   prev_CA = list_prev_CA[k2];
                                   prev_altloc = *prev_CA->altloc;
                                   // To form a peptide bond, currentCA and prevCA must be 
                                   // within radius and have the same altloc
                                   // identifier or one altloc blanc
                                   if (current_altloc == prev_altloc || current_altloc=="_" || prev_altloc=="_") {
                                        double distance = sqrt(Math<double>::sqr((current_CA->coordinates[0]) - (prev_CA->coordinates[0])) +
                                                               Math<double>::sqr((current_CA->coordinates[1]) - (prev_CA->coordinates[1])) +
                                                               Math<double>::sqr((current_CA->coordinates[2]) - (prev_CA->coordinates[2])));
                                        if (distance <= radius) {
                                             // Even if several altloc path exist, we just pick the first found.
                                             bond = true;
                                        }
                                   }
                              }
                         }

                         std::string chosen_altloc = current_altloc;
                         std::string chosen_altloc_prev = prev_altloc;

                         if (bond) {
                              // Remove atoms in previous residue with wrong altLoc
                              if (chosen_altloc_prev != "_") {
                                   int k1 = n_atoms(i,j-1)-1;
                                   current_altloc = chosen_altloc_prev;;
                                   while (k1>=0 and current_altloc != std::string("_")) {
                                        AtomData *current_atom = data[i][j-1][k1];
                                        prev_altloc = current_altloc;
                                        current_altloc = *current_atom->altloc;
                                        if (!(current_altloc == "_" || (current_altloc == prev_altloc))) {
                                             std::vector<AtomData *>::iterator it = data[i][j-1].begin()+k1;
                                             delete data[i][j-1][k1];
                                             data[i][j-1].erase(it);
                                             current_altloc = prev_altloc;
                                        }
                                        k1--;
                                   }
                              }

                              // Remove N atoms in current residue with wrong altLoc
                              if (chosen_altloc != "_") {
                                   int k1 = 0;
                                   while (k1 < n_atoms(i,j)) {
                                        AtomData *current_atom = data[i][j][k1];
                                        if ((*current_atom->atom == std::string("CA")) and current_atom != current_CA) {
                                             std::vector<AtomData *>::iterator it = data[i][j].begin()+k1;
                                             delete data[i][j][k1];
                                             data[i][j].erase(it);
                                             list_current_CA.clear();
                                             list_current_CA.push_back(current_CA);
                                        } else {
                                             k1++;
                                        }
                                   }
                              }
                         } else {
                              std::cerr << "WARNING: Chain break at residue " << j << " in chain " << i << ".\n";
                              std::vector<std::vector<std::vector<AtomData *> > >::iterator chain_iterator= data.begin() + i;
                              std::vector<std::vector<AtomData *> >::iterator res_iterator= data[i].begin() + j;
                         

                              std::vector<std::vector<AtomData *> > new_chain(res_iterator, data[i].end());
                              data[i].erase(res_iterator, data[i].end());

                              data.insert(chain_iterator+1, new_chain);
                         }
                    }
               }
          }
     }
     

     //! Split polypeptides into smaller pieces if chainbreaks are detected
     //!
     //! \param CA_only Determines whether to use full atom or CA_only version
     void split_by_chain_breaks(bool CA_only) {
          if (CA_only)
               split_by_chain_breaks_ca();
          else
               split_by_chain_breaks_fb();

     }

     //! Return chain size for each polypeptide
     //!
     //! \return Integer for each polypeptide
     std::vector<int> size(){
          std::vector<int> v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(n_residues(i));
          }
          return v;
     }
     
     //! Return AA sequences for each polypeptide
     //!
     //! \return sequence of amino acid indices for each polypeptide
     std::vector<std::vector<int> > get_sequence() const {
          bool status = false;
          std::vector<std::vector<int> > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<int>());
               for (int j=0; j<n_residues(i); j++) {
                    // Just use information from first atom
                    if (data[i][j][0]->residue) {
                         v[i].push_back(definitions::str_to_aa(*data[i][j][0]->residue));
                         status=true;
                    }
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No sequence data found\n");
          
          return v;
     }

     //! Return AA sequences for each polypeptide
     //!
     //! \return amino acid string for each polypeptide
     std::vector<std::string> get_sequence_string() const {
          bool status = false;
          std::vector<std::string> return_value;
          for (int i=0; i<n_polypeptides(); i++) {
               return_value.push_back(std::string());
               for (int j=0; j<n_residues(i); j++) {
                    // Just use information from first atom
                    if (data[i][j][0]->residue) {
                         return_value[i] += definitions::residue_name_short[definitions::str_to_aa[*data[i][j][0]->residue]];
                         status=true;
                    }
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No sequence data found\n");
          
          return return_value;
     }

     //! Return position sequences for each polypeptide
     //!
     //! \return vector of 3D coordinates for each polypeptide
     std::vector<std::vector<std::vector<Vector_3D> > > get_positions() const {
          bool status = false;
          std::vector<std::vector<std::vector<Vector_3D> > > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<std::vector<Vector_3D> >());
               for (int j=0; j<n_residues(i); j++) {
                    v[i].push_back(std::vector<Vector_3D>());
                    for (int k=0; k<n_atoms(i,j); k++) {             
                         v[i][j].push_back(Vector_3D(data[i][j][k]->coordinates));
                         status=true;
                    }
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No position data found\n");
          
          return v;
     }


     //! Return (phi, psi) sequences for each polypeptide
     //!
     //! \param include_omega Whether omega angle should be included
     //! \return vector of phi,psi pairs for each polypeptide
     std::vector<std::vector<std::vector<double> > > get_phi_psi(bool include_omega=false) const {
          bool status = false;

          std::vector<std::vector<std::vector<Vector_3D> > > positions = this->get_positions();

          std::vector<std::vector<std::vector <double> > > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<std::vector<double> >());
               for (int j=0; j<n_residues(i); j++) {
                    v[i].push_back(std::vector<double>());

                    double phi;
                    if (j==0) {
                         phi = std::numeric_limits<double>::quiet_NaN();

                         // Set first phi based on terminal Hydrogen position is available
                         for (unsigned int k=0; k<positions[i][j].size(); ++k) {
                              if ((*data[i][j][k]->atom) == "H" || (*data[i][j][k]->atom) == "H1" || (*data[i][j][k]->atom) == "1H") {
                                   phi = calc_dihedral(positions[i][j][k], positions[i][j][0], positions[i][j][1], positions[i][j][2]);
                                   break;
                              }
                         }
                    } else {
                         phi = calc_dihedral(positions[i][j-1][2], positions[i][j][0], positions[i][j][1], positions[i][j][2]);
                    }
                    
                    double psi;
                    if (j==n_residues(i)-1) {
                         psi = std::numeric_limits<double>::quiet_NaN();

                         // Set first phi based on terminal Hydrogen position is available
                         for (unsigned int k=0; k<positions[i][j].size(); ++k) {
                              if ((*data[i][j][k]->atom) == "O" || (*data[i][j][k]->atom) == "OXT") {
                                   psi = calc_dihedral(positions[i][j][0], positions[i][j][1], positions[i][j][2], positions[i][j][k]);
                                   break;
                              }
                         }

                    } else {
                         psi = calc_dihedral(positions[i][j][0], positions[i][j][1], positions[i][j][2], positions[i][j+1][0]);
                    }

                    double omega;
                    if (include_omega) {
                         if (j==0) {
                              omega = std::numeric_limits<double>::quiet_NaN();
                         } else {
                              omega = calc_dihedral(positions[i][j-1][1], positions[i][j-1][2], positions[i][j][0], positions[i][j][1]);
                         }
                    }
                    
                    v[i][j].push_back(phi);
                    v[i][j].push_back(psi);
                    if (include_omega) {
                         v[i][j].push_back(omega);
                    }
                    
                    status=true;
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No position data found\n");
          
          return v;
     }
     

     //! Return sequence of omega angles for each polypeptide
     //!
     //! \return vector of omega angles for each polypeptide
     std::vector<std::vector<double> > get_omega() const {
          bool status = false;
          
          std::vector<std::vector<std::vector <double> > > angles = get_phi_psi(true);

          std::vector<std::vector<double> > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<double>());
               for (int j=0; j<n_residues(i); j++) {
                    double omega = angles[i][j][2];
                    v[i].push_back(omega);
                    status = true;
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No cis data found\n");
          
          return v;
     }
     

     //! Return sequence of cis/trans states for each polypeptide
     //!
     //! \return vector of cis/trans states for each polypeptide
     std::vector<std::vector<int> > get_cis() const {
          bool status = false;
          
          std::vector<std::vector<std::vector <double> > > angles = get_phi_psi(true);

          std::vector<std::vector<int> > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<int>());
               for (int j=0; j<n_residues(i); j++) {
                    double omega = angles[i][j][2];
                    int cis = 0;
                    if (std::isnan(omega)) {
                         set_initial_value(cis);
                    } else if (std::fabs(omega) < M_PI/4) {
                         cis = 1;
                    }
                    v[i].push_back(cis);
                    status = true;
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No cis data found\n");
          
          return v;
     }
     
     //! Return (theta, tau) sequences for each polypeptide
     //!
     //! \return vector of theta,tau pairs for each polypeptide
     std::vector<std::vector<std::vector <double> > > get_theta_tau() const {
          bool status = false;

          std::vector<std::vector<std::vector<Vector_3D> > > positions = this->get_positions();

          std::vector<std::vector<std::vector <double> > > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<std::vector<double> >());
               for (int j=0; j<n_residues(i); j++) {
                    v[i].push_back(std::vector<double>());

                    double theta;
                    if ((j==0) || (j==n_residues(i)-1)) {
                         theta = std::numeric_limits<double>::quiet_NaN();
                    } else {
                         theta = calc_angle(positions[i][j-1][1], positions[i][j][1], positions[i][j+1][1]);
                    }
                    
                    double tau;
                    if ((j==0) || (j==1) || j==n_residues(i)-1) {
                         tau = std::numeric_limits<double>::quiet_NaN();
                    } else {
                         tau = calc_dihedral(positions[i][j-2][1], positions[i][j-1][1], positions[i][j][1], positions[i][j+1][1]);
                    }

                    v[i][j].push_back(theta);
                    v[i][j].push_back(tau);
                    
                    status=true;
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No position data found\n");
          
          return v;
     }

     //! Get resseq sequence for each polypeptide
     //!
     //! \return Vector of resseq vectors
     std::vector<std::vector<int> > get_resseq() const {
          std::vector<std::vector<int> > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<int>());
               for (int j=0; j<n_residues(i); j++) {
                    v.back().push_back(*data[i][j][0]->resseq);
               }
          }
          return v;
     }
     
     
     //! Get first resseq of each polypeptide
     //!
     //! \return Vector of resseq start indices
     std::vector<int> get_resseq_offset() const {
          std::vector<int> v;
          for (int i=0; i<n_polypeptides(); i++) {
               if (data[i][0][0]->resseq) {
                    v.push_back(*data[i][0][0]->resseq);
               }
          }
          return v;
     }

     
     
     //! Return real secondary structure sequences for each polypeptide
     //!
     //! \return Vector of secondary structure indices for each polypeptide
     std::vector<std::vector<int> > get_real_ss() const {
          bool status = false;
          std::vector<std::vector<int> > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<int>());
               for (int j=0; j<n_residues(i); j++) {
                    // Just use information from first atom
                    if (data[i][j][0]->real_ss) {
                         v[i].push_back(*data[i][j][0]->real_ss);
                         status=true;
                    }
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No realSS data found\n");
               
          return v;
     }

     //! Return predicted secondary sequences for each polypeptide
     //!
     //! \return Vector of secondary structure indices for each polypeptide
     std::vector<std::vector<int> > get_pred_ss() const {
          bool status = false;
          std::vector<std::vector<int> > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<int>());
               for (int j=0; j<n_residues(i); j++) {
                    // Just use information from first atom
                    if (data[i][j][0]->pred_ss) {
                         v[i].push_back(*data[i][j][0]->pred_ss);
                         status=true;
                    }
               }
          }
          if (!status)
               fprintf(stderr, "Warning: No predSS data found\n");

          return v;
     }

     //! Return predicted secondary structure propensity sequences for each polypeptide
     //!
     //! \return vector of H,E,C propensities for each polypeptide
     std::vector<std::vector<std::vector<double> > > get_pred_prob_ss() const {
          bool status = false;
          std::vector<std::vector<std::vector<double> > > v;
          for (int i=0; i<n_polypeptides(); i++) {
               v.push_back(std::vector<std::vector<double> >());
               for (int j=0; j<n_residues(i); j++) {
                    std::vector<double> probs;
                    // Just use information from first atom
                    if (data[i][j][0]->pred_prob_ss_C &&
                        data[i][j][0]->pred_prob_ss_H &&
                        data[i][j][0]->pred_prob_ss_E ) {

                         probs.push_back(*data[i][j][0]->pred_prob_ss_H);
                         probs.push_back(*data[i][j][0]->pred_prob_ss_E);
                         probs.push_back(*data[i][j][0]->pred_prob_ss_C);
                         status=true;
                    }
                    v[i].push_back(probs);
               }
          }

          if (!status)
               fprintf(stderr, "Warning: No predProbSS data found\n");
          return v;
     }

     //! Overload [] indexing operator for ProteinData
     //!
     //! \param pp_index polypeptide index
     //! \return Vector(residues) of vector of AtomData
     std::vector<std::vector<AtomData *> > operator[](int pp_index) {
          return data[pp_index];
     }

     //! Save data to file
     //!
     //! \param filename Destination filename
     void save(std::string filename) {
          std::ofstream out(filename.data());
          if (!out.good()) {
               std::cerr << "Couldn't open file " << filename.data() << "\n";
               return;
          }
          out << *this;
          out.close();
     }

     //! Load data from file
     //!
     //! \param filename Source filename
     void load(std::string filename) {
          std::ifstream in(filename.data());
          if (!in.good()) {
               std::cerr << "Couldn't open file " << filename.data() << "\n";
               return;
          }
          std::string pp_index_label;
          std::string pp_string;
          std::string res_index_label;
          std::string res_string;
          while (in >> pp_index_label >> pp_string >> res_index_label >> res_string) {
               int pp_index;
               std::stringstream(pp_string) >> pp_index;
               int res_index;
               std::stringstream(res_string) >> res_index;
               new_atom(pp_index, res_index);
               in >> *current_atom();
          }
     }

     //! Overload output operator for ProteinData
     friend std::ostream & operator<<(std::ostream &o, ProteinData &data) {
          for (int i=0; i<data.n_polypeptides(); i++) {
               for (int j=0; j<data.n_residues(i); j++) {
                    for (int k=0; k<data.n_atoms(i,j); k++) {
                         o << "PPINDEX " << std::setw(2) <<  i << " " << "RESINDEX " << std::setw(3) << j << " " << *data[i][j][k] << "\n";
                    }
               }
          }          
          return o;
     }

};

}

#endif
