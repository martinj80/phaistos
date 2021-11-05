// profasi_hydrophobicity.h --- PROFASI side chain hydrophobic attraction energy term
// Copyright (C) 2010 Pengfei Tian, Wouter Boomsma
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

// Implementation of the PROFASI force field
// An effective all-atom potential for proteins, Irback, Mitternacht, Mohanty, PMC Biophysics 2009


#ifndef PROFASI_HYDROPHOBICITY_H
#define PROFASI_HYDROPHOBICITY_H

#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"
#include "protein/iterators/pair_iterator_chaintree.h"
#include <cmath>

namespace phaistos {

//! PROFASI hydrophobicity energy term
template<typename DERIVED_CLASS>
class TermProfasiHydrophobicityBase: public EnergyTermCommon<DERIVED_CLASS,
                                                             ChainFB>,
                                     public EnergyTermCommonProfasi {
private:
     
     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;
     
protected:
     
     //! Sidechain atoms involved in hydrophobicity calculations
     std::vector<std::vector<definitions::AtomEnum> > side_chain_atoms;
     
public:
     
     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:
          
          //! m parameter
          std::vector<double> m;
          
          // There are two groups: determining the value (0.75 or 1) of gamma parameter (lower case)
          std::vector<int> group;
          
          //! Lower limit for linear interpolation
          double g_lower;

          //! Upper limit for linear interpolation
          double g_upper;
          
          //! Cutoff of distance for cached version
          double r_cut;
          
          //! Constructor. Defines default values for settings object
          Settings(std::vector<double> m = vector_utils::make_vector(0.0, // ALA
                                                                     0.0, // CYS
                                                                     0.0, // ASP
                                                                     0.0, // GLU
                                                                     1.6, // PHE
                                                                     0.0, // GLY
                                                                     0.0, // HIS
                                                                     0.8, // ILE
                                                                     0.4, // LYS
                                                                     0.8, // LEU
                                                                     0.4, // MET
                                                                     0.0, // ASN
                                                                     0.8, // PRO
                                                                     0.0, // GLN
                                                                     0.3, // ARG
                                                                     0.0, // SER
                                                                     0.0, // THR
                                                                     0.6, // VAL
                                                                     1.6, // TRP
                                                                     1.1),// TYR)
                   std::vector<int> group = vector_utils::make_vector(0, // ALA
                                                                      0, // CYS
                                                                      0, // ASP
                                                                      0, // GLU
                                                                      1, // PHE
                                                                      0, // GLY
                                                                      0, // HIS
                                                                      2, // ILE
                                                                      0, // LYS
                                                                      2, // LEU
                                                                      0, // MET
                                                                      0, // ASN
                                                                      1, // PRO
                                                                      0, // GLN
                                                                      0, // ARG
                                                                      0, // SER
                                                                      0, // THR
                                                                      2, // VAL
                                                                      1, // TRP
                                                                      1), // TYR
                   double g_lower = 3.7 * 3.7,
                   double g_upper = 4.5 * 4.5,
                   double r_cut = 4.5) 
               : m(m), group(group), g_lower(g_lower), 
                 g_upper(g_upper), r_cut(r_cut) {
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "m:" << settings.m << "\n";
               o << "group:" << settings.group << "\n";
               o << "g-lower:" << settings.g_lower << "\n";
               o << "g-upper:" << settings.g_upper << "\n";
               o << "r-cut:" << settings.r_cut << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }          
          
     } settings;     //!< Local settings object 
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrophobicityBase(ChainFB *chain, std::string name,
                                   const Settings &settings = Settings(),
                                   RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine) {
          
          //! Import protein definitions (such as residue names)
          using namespace definitions;
          
          //! Import phaistos::make_vector namespace
          using namespace vector_utils;
          
          //Atoms used in the calculation of the contact measure
          side_chain_atoms.resize(20);
          side_chain_atoms[PRO] = make_vector(CB, CG, CD);
          side_chain_atoms[TYR] = make_vector(CG, CD1, CE1, CZ, CE2, CD2);
          side_chain_atoms[VAL] = make_vector(CB, CG1, CG2);
          side_chain_atoms[ILE] = make_vector(CB, CG1, CG2, CD1);
          side_chain_atoms[LEU] = make_vector(CB, CG, CD1, CD2);
          side_chain_atoms[MET] = make_vector(CB, CG, SD, CE);
          side_chain_atoms[PHE] = make_vector(CG, CD1, CE1, CZ, CE2, CD2);
          side_chain_atoms[TRP] = make_vector(CG, CD1, CD2, CE3, CZ3, CH2);
          side_chain_atoms[ARG] = make_vector(CB, CG);
          side_chain_atoms[LYS] = make_vector(CB, CG, CD);
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrophobicityBase(const TermProfasiHydrophobicityBase &other,
                                   RandomNumberEngine *random_number_engine,
                                   int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            side_chain_atoms(other.side_chain_atoms) {}
     

     //! Calculate hydrophobicity protential of given two residues(normal version)
     //! \param res1 The first residue in the hydrophobicity interaction
     //! \param ni_atoms The atoms on the first residue which are involved in the contact measure
     //! \param res2 The second residue in the hydrophobicity interaction 
     //! \param nj_atoms The atoms on the second residue which are involved in the contact measure
     //! \return Total hydrophobicity protential between two residues
     double calculate_contribution(ResidueFB *res1,
                                   const std::vector<definitions::AtomEnum> &ni_atoms,
                                   ResidueFB *res2,
                                   const std::vector<definitions::AtomEnum> &nj_atoms) {
          
          double m_ij     = 0.0;
          double c_ij     = 0.0;
          double gamma_ij = 0.0;
          double E_hp     = 0.0;
          double GAMMA_ij = 0.0;
          unsigned int n;
          
          // Obtain the residue type
          int re1 = res1->residue_type;
          int re2 = res2->residue_type;
          
          // Calculate the parameter m_ij
          m_ij = this->settings.m[re1] + this->settings.m[re2];
          
          // If chain distance of the two residues is equal to  2
          if (res2->index - res1->index == 2)
               m_ij *= 0.5;
          
          unsigned int ni = ni_atoms.size();
          unsigned int nj = nj_atoms.size();
          double r2min[ni + nj];
          
          // If the two residues both have atoms that could have hydrophobicity contact
          if (ni != 0 && nj != 0) {
               
               if ( (this->settings.group[re1] != 0) && (this->settings.group[re1] == this->settings.group[re2]) )
                    gamma_ij = 0.75;
               else
                    gamma_ij = 1.0;
               
               for (n = 0; n < ni + nj; n++)
                    r2min[n] = 1e10;
               
               for (unsigned int i = 0; i < ni; i++) {
                    for (unsigned int j = 0; j < nj; j++) {
                         
                         Atom *atom1 = (*res1)[ni_atoms[i]];
                         Atom *atom2 = (*res2)[nj_atoms[j]];
                         
                         double r_ij =
                              (atom1->position - atom2->position).norm();
                         double r2 = r_ij * r_ij;
                         
                         if (r2 < r2min[i])
                              r2min[i] = r2;
                         
                         if (r2 < r2min[ni + j])
                              r2min[ni + j] = r2;
                         
                    }
               }
               
               for (n = 0; n < ni + nj; n++) {
                    if (r2min[n] > settings.g_upper)
                         continue;
                    if (r2min[n] < settings.g_lower)
                         GAMMA_ij++;
                    else
                         GAMMA_ij += ((settings.g_upper - r2min[n])
                                      / (settings.g_upper - settings.g_lower));
               }
               
               c_ij = std::min(1.0, GAMMA_ij / (gamma_ij * (ni + nj)));
               
               E_hp = -m_ij * c_ij;
          }
          
          return E_hp;
     }
     
     //! Calculate hydrophobicity protential of given two residues(improved version)
     //! \param res1 The first residue in the hydrophobicity interaction
     //! \param res2 The second residue in the hydrophobicity interaction 
     //! \return Total hydrophobicity protential between two residues
     double calculate_contribution_improved(ResidueFB *res1, ResidueFB *res2) {
          
          double m_ij     = 0.0;
          double c_ij     = 0.0;
          double gamma_ij = 0.0;
          double E_hp     = 0.0;
          double GAMMA_ij = 0.0;
          unsigned int n;
          
          int re1 = res1->residue_type;
          int re2 = res2->residue_type;
          
          const std::vector<definitions::AtomEnum> & ni_atoms = this->side_chain_atoms[res1->residue_type];
          const std::vector<definitions::AtomEnum> & nj_atoms = this->side_chain_atoms[res2->residue_type];
          
          m_ij = this->settings.m[re1] + this->settings.m[re2];
          
          if (res2->index - res1->index == 2)
               m_ij *= 0.5;
          
          unsigned int ni = ni_atoms.size();
          unsigned int nj = nj_atoms.size();
          double r2min[ni + nj];
          
          if ( (this->settings.group[re1] != 0) && (this->settings.group[re1] == this->settings.group[re2]) )
               gamma_ij = 0.75;
          else
               gamma_ij = 1.0;
          
          for (n = 0; n < ni + nj; n++)
               r2min[n] = 1e10;
          
          for (unsigned int i = 0; i < ni; i++) {
               for (unsigned int j = 0; j < nj; j++) {
                    
                    Atom *atom1 = (*res1)[ni_atoms[i]];
                    Atom *atom2 = (*res2)[nj_atoms[j]];
                    
                    double r_ij =
                         (atom1->position - atom2->position).norm();
                    double r2 = r_ij * r_ij;
                    
                    if (r2 < r2min[i])
                         r2min[i] = r2;
                    
                    if (r2 < r2min[ni + j])
                         r2min[ni + j] = r2;
                    
               }
          }
          
          for (n = 0; n < ni + nj; n++) {
               if (r2min[n] > settings.g_upper)
                    continue;
               if (r2min[n] < settings.g_lower)
                    GAMMA_ij++;
               else
                    GAMMA_ij += ((settings.g_upper - r2min[n])
                                 / (settings.g_upper - settings.g_lower));
          }
          
          c_ij = std::min(1.0, GAMMA_ij / (gamma_ij * (ni + nj)));
          
          E_hp = -m_ij * c_ij;
          
          return E_hp;
     }
     
};
     
//! PROFASI hydrophobicity energy term
class TermProfasiHydrophobicity: public TermProfasiHydrophobicityBase<TermProfasiHydrophobicity> {
     
     //! For convenience, define local base class
     typedef phaistos::TermProfasiHydrophobicityBase<TermProfasiHydrophobicity> TermProfasiHydrophobicityBase;
          
public:
          
     //! Use same settings as base class
     typedef TermProfasiHydrophobicityBase::Settings Settings;

     //! Local Settings object
     const Settings settings;
          
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrophobicity(ChainFB *chain,
                               const Settings &settings=Settings(),
                               RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiHydrophobicityBase(chain, "profasi-hydrophobicity", settings, random_number_engine),
            settings(settings) {}
          
     // Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrophobicity(const TermProfasiHydrophobicity &other, 
                               RandomNumberEngine *random_number_engine,
                               int thread_index, ChainFB *chain)
          : TermProfasiHydrophobicityBase(other, random_number_engine, thread_index, chain),
            settings(other.settings) {}
          
     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return hydrophobicity potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
               
          //! Import protein definitions (such as residue names)
          using namespace definitions;
               
          double energy = 0.0;
               
          // Go through all the residues pairs
          for (ResidueIterator<ChainFB> it1(*this->chain); !it1.end(); ++it1) {
               for (ResidueIterator<ChainFB> it2(it1+2); !it2.end(); ++it2) {
                         
                    ResidueFB *res1 = &*it1;
                    ResidueFB *res2 = &*it2;
                         
                    const std::vector<AtomEnum> & ni_atoms = this->side_chain_atoms[res1->residue_type];
                    const std::vector<AtomEnum> & nj_atoms = this->side_chain_atoms[res2->residue_type];
                         
                    energy += this->calculate_contribution(res1, ni_atoms,
                                                           res2, nj_atoms);
               }
          }
               
          return this->profasi_energy_in_kcal_per_mol * energy;
     }
          
};
     
//! PROFASI hydrophobicity energy term - cached version
class TermProfasiHydrophobicityCached: public TermProfasiHydrophobicityBase<TermProfasiHydrophobicityCached> {
          
     //! For convenience, define local base class
     typedef phaistos::TermProfasiHydrophobicityBase<TermProfasiHydrophobicityCached> TermProfasiHydrophobicityBase;
          
     //! Specifies, for each node, which atoms contribute to the energy term
     std::vector<std::pair<ResidueFB*,std::vector<definitions::AtomEnum> > > active_atoms;
          
     //! For all chaintree nodes at the lowest level, find the atoms that contribute
     void find_active_atoms() {
               
          //! Import protein definitions (such as residue names)
          using namespace definitions;
               
          //! The chaintree level (level 0 is omitted for efficiency reasons)
          int level = 1;
               
          active_atoms = std::vector<std::pair<ResidueFB*,std::vector<AtomEnum> > >(this->chain->get_chain_tree()->get_nodes_at_level(level),
                                                                                    std::make_pair((ResidueFB*)NULL,std::vector<AtomEnum>()));
          int offset = this->chain->get_chain_tree()->get_level_start_index(level);
               
          for (unsigned int i=0; i<active_atoms.size(); ++i) {

               NodeType *node = this->chain->get_chain_tree()->nodes[i+offset];
                    
               //! Merge atoms from the two child nodes
               std::vector<Atom*> merged_child_atoms = node->child1->atoms;
               if (node->child2) {
                    merged_child_atoms.insert(merged_child_atoms.end(),
                                              node->child2->atoms.begin(),
                                              node->child2->atoms.end());
               }
               
               //! Iterate over all child node atoms and detect active atoms
               for (unsigned int j=0; j<merged_child_atoms.size(); ++j) {
                    Atom *atom = merged_child_atoms[j];
                    ResidueFB *residue = &(*this->chain)[atom->residue->index];
                    AtomEnum atom_type = atom->atom_type;
                    ResidueEnum residue_type = residue->residue_type;
                    
                    std::vector<AtomEnum> residue_atoms = this->side_chain_atoms[residue_type];
                    
                    for (unsigned int k=0; k<residue_atoms.size(); ++k) {
                         if (residue_atoms[k] == atom_type) {
                              active_atoms[i].second.push_back(atom_type);
                              if (active_atoms[i].first == NULL) {
                                   active_atoms[i].first = residue;
                              } else {
                                   // Make sure that all atoms in a node point to the same residue
                                   assert(active_atoms[i].first == residue);
                              }
                              break;
                         }
                    }
               }
          }
     }
          
public:
          
     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;
          
     //! Cached iterator (used when iterating over node pairs)
     CachedIterator<chaintree::PairIterator<ChainFB,NodeType,NodeType> > cached_it;
          
     //! Use same settings as base class
     typedef TermProfasiHydrophobicityBase::Settings Settings;

     //! Local Settings object
     const Settings settings;
          
     //! Chaintree iterator settings - described in detail in chaintree.h
     chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings iterator_settings;
          
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrophobicityCached(ChainFB *chain,
                                     const Settings &settings=Settings(),
                                     RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiHydrophobicityBase(chain, "profasi-hydrophobicity-cached", settings, random_number_engine),
            cached_it(*chain),
            settings(settings) {
               
          // Only evaluate modified pairs (as always when caching)
          bool only_modified_pairs = true;
          int minimum_residue_distance = 2;
          bool omit_lowest_level = true;
          iterator_settings = chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings(settings.r_cut,
                                                                                           only_modified_pairs,
                                                                                           minimum_residue_distance,
                                                                                           omit_lowest_level);
               
          find_active_atoms();
     }
          
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrophobicityCached(const TermProfasiHydrophobicityCached &other, 
                                     RandomNumberEngine *random_number_engine,
                                     int thread_index, ChainFB *chain)
          : TermProfasiHydrophobicityBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings),
            iterator_settings(other.iterator_settings) {
               
          find_active_atoms();
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return hydrophobicity potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
               
          //! Import protein definitions (such as residue names)
          using namespace definitions;
               
          // Cached iterator evaluation
          for (cached_it(*this->chain, iterator_settings); !cached_it.end(); ++cached_it) {
                    
               double energy = 0.0;
                    
               NodeType *node1 = cached_it->first;
               NodeType *node2 = cached_it->second;
                    
               ResidueFB *res1 = active_atoms[node1->index].first;
               const std::vector<AtomEnum> &atoms_node1 = active_atoms[node1->index].second;
               ResidueFB *res2 = active_atoms[node2->index].first;
               const std::vector<AtomEnum> &atoms_node2 = active_atoms[node2->index].second;
                    
               if (res1 && res2) {
                    energy += this->calculate_contribution(res1,
                                                           atoms_node1,
                                                           res2,
                                                           atoms_node2);
                    cached_it.register_contribution(energy);
               }
          }
               
          double energy_ev = cached_it.compute_total();
               
          return this->profasi_energy_in_kcal_per_mol * energy_ev;
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
     
     
//! PROFASI hydrophobicity energy term - improved version
class TermProfasiHydrophobicityImproved: public TermProfasiHydrophobicityBase<TermProfasiHydrophobicityImproved> {
          
     //! For convenience, define local base class
     typedef phaistos::TermProfasiHydrophobicityBase<TermProfasiHydrophobicityImproved> TermProfasiHydrophobicityBase;
          
     //! Specifies, for each node, which atoms contribute to the energy term
     std::vector<std::pair<ResidueFB*,std::vector<definitions::AtomEnum> > > active_atoms;
          
     //! For all chaintree nodes at the lowest level, find the atoms that contribute
     std::vector< ResidueFB* > hp_residues;

     //! Put the residues which has hydrophobic atoms in vector          
     void find_hydrophobicity_residues() {
               
          //! Import protein definitions (such as residue names)
          using namespace definitions;
               
          int num_hp_residue = 0;
               
          for (ResidueIterator<ChainFB> it1(*this->chain); !it1.end(); ++it1) {
               switch (it1->residue_type) {
                         
               case PRO:
               case TYR:
               case VAL:
               case ILE:
               case LEU:
               case MET:
               case PHE:
               case TRP:
               case ARG:
               case LYS:
                    num_hp_residue += 1;
                    break;
               default:
                    break;
               }
          }
               
          hp_residues.resize(num_hp_residue);
               
          int i_hp_res = 0;
               
          for (ResidueIterator<ChainFB> it2(*this->chain); !it2.end(); ++it2) {
               ResidueFB* res = &(*it2);
                    
               switch (it2->residue_type) {
                         
               case PRO:
               case TYR:
               case VAL:
               case ILE:
               case LEU:
               case MET:
               case PHE:
               case TRP:
               case ARG:
               case LYS:
                    hp_residues[i_hp_res++] = res;
                    break;
               default:
                    break;
               }
          }
     }
          
public:
          
     //! Use same settings as base class
     typedef TermProfasiHydrophobicityBase::Settings Settings;
          
     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;
     
     //!Local Settings object
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrophobicityImproved(ChainFB *chain,
                                       const Settings &settings=Settings(),
                                       RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiHydrophobicityBase(chain, "profasi-hydrophobicity-improved", settings, random_number_engine),
            settings(settings) {
               find_hydrophobicity_residues();
     }
          
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrophobicityImproved(const TermProfasiHydrophobicityImproved &other, 
                                       RandomNumberEngine *random_number_engine,
                                       int thread_index, ChainFB *chain)
          : TermProfasiHydrophobicityBase(other, random_number_engine, thread_index, chain),
            settings(other.settings){

          find_hydrophobicity_residues();
     }
          
     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return hydrophobicity potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
               
          double energy = 0;

          // Iterate only the resiude pairs which both have atoms that can cause hydrophobic attraction
          for(unsigned int i=0; i< hp_residues.size(); i++ ){
               for(unsigned int j = i+1; j< hp_residues.size(); j++ ){
                         
                    ResidueFB *res1 = hp_residues[i];
                    ResidueFB *res2 = hp_residues[j];
                         
                    // If chain distance of two residues are longer than 1
                    if ((res2->index - res1->index) > 1)
                         energy += this->calculate_contribution_improved(res1, res2);
                         
               } 
          }
          return this->profasi_energy_in_kcal_per_mol * energy;      
     }
          
};
     
}

#endif
