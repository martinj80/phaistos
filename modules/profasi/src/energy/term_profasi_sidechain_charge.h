// profasi_sidechain_charge.h --- PROFASI side chain charged interaction energy term
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


#ifndef PROFASI_SIDECHAIN_CHARGE_H
#define PROFASI_SIDECHAIN_CHARGE_H

#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace phaistos {
     
//! PROFASI side chain charged interaction energy term - base class
template<typename DERIVED_CLASS>
class TermProfasiSidechainChargeBase: public EnergyTermCommon<DERIVED_CLASS,
                                                              ChainFB>,
                                      public EnergyTermCommonProfasi {
private:
     
     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;
     
public:
     
     //! Used to ave charged atoms according to residue type
     std::vector<std::vector<definitions::AtomEnum> > side_chain_charged_atoms;

     //! Used to save charged atoms(H1, H2, H3, O, OXT) which are on terminals
     std::vector<std::vector<definitions::AtomEnum> > terminal_charged_atoms;

     //! Used to save two kinds of charged atoms (on the internal and terminal of the chain) 
     std::vector<std::vector<definitions::AtomEnum> > all_sidechain_charged_atoms;

     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:
          
          //! m parameter
          std::vector<double> m;
          
          //! The group indicates whether a residue has a sidechain charge
          std::vector<int> group;
          
          //! Lower limit for linear interpolation
          double g_lower;

          //! upper limit for linear interpolation
          double g_upper;
          
          //! cutoff of distance for cached version
          double r_cut;
          
          //! Constructor. Defines default values for settings object
          Settings(std::vector<double> m = vector_utils::make_vector(0.0, // ALA
                                                                     0.0, // CYS
                                                                     -1.0, // ASP
                                                                     -1.0, // GLU
                                                                     0.0, // PHE
                                                                     0.0, // GLY
                                                                     0.0, // HIS
                                                                     0.0, // ILE
                                                                     1.0, // LYS
                                                                     0.0, // LEU
                                                                     0.0, // MET
                                                                     0.0, // ASN
                                                                     0.0, // PRO
                                                                     0.0, // GLN
                                                                     1.0, // ARG
                                                                     0.0, // SER
                                                                     0.0, // THR
                                                                     0.0, // VAL
                                                                     0.0, // TRP
                                                                     0.0), // TYR
                   std::vector<int> group = vector_utils::make_vector(0, // ALA
                                                                      0, // CYS
                                                                      1, // ASP
                                                                      1, // GLU
                                                                      0, // PHE
                                                                      0, // GLY
                                                                      0, // HIS
                                                                      0, // ILE
                                                                      1, // LYS
                                                                      0, // LEU
                                                                      0, // MET
                                                                      0, // ASN
                                                                      0, // PRO
                                                                      0, // GLN
                                                                      1, // ARG
                                                                      0, // SER
                                                                      0, // THR
                                                                      0, // VAL
                                                                      0, // TRP
                                                                      0), // TYR
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
          
     } settings;
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiSidechainChargeBase(ChainFB *chain, std::string name,
                                    const Settings &settings = Settings(),
                                    RandomNumberEngine *random_number_engine = &random_global) 
          : EnergyTermCommon(chain, name, settings, random_number_engine) {
          
          //! Import protein definitions (such as residue names)
          using namespace definitions;
          
          //! Import phaistos::make_vector namespace
          using namespace vector_utils;
          
          side_chain_charged_atoms.resize(20);
          side_chain_charged_atoms[ARG] = make_vector(NE, CZ, NH1, NH2);
          side_chain_charged_atoms[LYS] = make_vector(HZ1, HZ2, HZ3);
          side_chain_charged_atoms[ASP] = make_vector(OD1, OD2);
          side_chain_charged_atoms[GLU] = make_vector(OE1, OE2);
          
          terminal_charged_atoms.resize(2);
          terminal_charged_atoms[0] = make_vector(H1, H2, H3);
          terminal_charged_atoms[1] = make_vector(O, OXT);

          all_sidechain_charged_atoms.resize(22);
          all_sidechain_charged_atoms[ARG] = make_vector(NE, CZ, NH1, NH2);
          all_sidechain_charged_atoms[LYS] = make_vector(HZ1, HZ2, HZ3);
          all_sidechain_charged_atoms[ASP] = make_vector(OD1, OD2);
          all_sidechain_charged_atoms[GLU] = make_vector(OE1, OE2);
          all_sidechain_charged_atoms[20] = make_vector(H1, H2, H3);
          all_sidechain_charged_atoms[21] = make_vector(O, OXT);
     }
     
     // Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiSidechainChargeBase(const TermProfasiSidechainChargeBase &other, 
                                    RandomNumberEngine *random_number_engine,
                                    int thread_index,
                                    ChainFB *chain) 
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            side_chain_charged_atoms(other.side_chain_charged_atoms),
            terminal_charged_atoms(other.terminal_charged_atoms),
            all_sidechain_charged_atoms(other.all_sidechain_charged_atoms) {
     }
     
     // Calculate the gamma parameter   
     //! \param ni_atoms The charged atoms on the first residue which are involved in the contact measure
     //! \param nj_atoms The charged atoms on the second residue which are involved in the contact measure
     //! \param res1 The first residue which has charged sidechain
     //! \param res2 The second residue which has charged sidechain 
     //! \return the gamma parameter for each group of charged atoms
     double calculate_gamma(const std::vector<definitions::AtomEnum> &ni_atoms,
                            const std::vector<definitions::AtomEnum> &nj_atoms,
                            ResidueFB *res1,
                            ResidueFB *res2) {
          
          unsigned int n;
          unsigned int ni = ni_atoms.size();
          unsigned int nj = nj_atoms.size();
          double r2min[ni + nj];
          double gamma = 0.0;
          
          for (n = 0; n < ni + nj; n++)
               r2min[n] = 1e10;
          
          for (unsigned int i = 0; i < ni; i++) {
               for (unsigned int j = 0; j < nj; j++) {
                    
                    Atom *atom1 = (*res1)[ni_atoms[i]];
                    Atom *atom2 = (*res2)[nj_atoms[j]];
                    
                    double r_ij = (atom1->position - atom2->position).norm();
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
               
               if (r2min[n] < settings.g_lower) {
                    gamma++;
               } else {
                    gamma += ((settings.g_upper - r2min[n]) / (settings.g_upper
                                                               - settings.g_lower));
               };
          }
          return gamma;
     }
     
     // Calculate the sidechain charged interaction (normal version)   
     //! \param res1 The residue which has the first charged sidechain atoms group in the interaction
     //! \param res2 The residue which has the second charged sidechain atoms group in the interaction
     //! \return The charged sidechain potential between two groups
     double calculate_contribution(ResidueFB *res1, ResidueFB *res2) {
          
          // Import protein definitions (such as residue names)
          using namespace definitions;
          
          double m_ij = 0.0;
          double c_ij = 0.0;
          double E_sc = 0.0;
          
          int re1 = res1->residue_type;
          int re2 = res2->residue_type;
          
          m_ij = 1.5 * settings.m[re1] * settings.m[re2];
          
          const std::vector<AtomEnum> & n_i = side_chain_charged_atoms[re1];
          const std::vector<AtomEnum> & n_j = side_chain_charged_atoms[re2];
          const std::vector<AtomEnum> & n_N = terminal_charged_atoms[0];
          const std::vector<AtomEnum> & n_C = terminal_charged_atoms[1];
          
          int ni = n_i.size();
          int nj = n_j.size();
          int n_Nnum = 3; //number of charged atoms on the Nterminal
          int n_Cnum = 2; //number of charged atoms on the Cterminal

          double GAMMA_ij = calculate_gamma(n_i, n_j, res1, res2);
          
          if (ni != 0 && nj != 0) {
               c_ij = std::min(1.0, GAMMA_ij / (ni + nj));
               E_sc += m_ij * c_ij;
          }
          
          if (res1->terminal_status == NTERM && re1 != PRO) { //PRO dont have (H1 H2 H3)
               double GAMMA_Nj = calculate_gamma(n_N, n_j, res1, res2);
               
               m_ij = 1.5 * settings.m[re2];
               c_ij = std::min(1.0, GAMMA_Nj / (n_Nnum + nj));
               E_sc += m_ij * c_ij;
               
               if (res2->terminal_status == CTERM) {
                    double GAMMA_NC = calculate_gamma(n_N, n_C, res1, res2);
                    
                    m_ij = -1.5;
                    c_ij = std::min(1.0, GAMMA_NC / (n_Nnum + n_Cnum));
                    E_sc += m_ij * c_ij;
                    
               }
          }
          
          if (res2->terminal_status == CTERM) {
               if (ni != 0) {
                    double GAMMA_iC = calculate_gamma(n_i, n_C, res1, res2);
                    
                    m_ij = -1.5 * settings.m[re1];
                    c_ij = std::min(1.0, GAMMA_iC / (ni + n_Cnum));
               E_sc += m_ij * c_ij;
               }
          }
          return E_sc;
     }

     // Calculate the sidechain charged interaction (improved version)   
     //! \param res1 The residue which has the first charged sidechain atoms group in the interaction
     //! \param res2 The residue which has the second charged sidechain atoms group in the interaction
     //! \param m_i Total number of sidechain charged atoms on the first group
     //! \param m_j Total number of sidechain charged atoms on the second group
     //! \param sc_atom_index_i The index of the first charged sidechain atoms group
     //! \param sc_atom_index_j The index of the second charged sidechain atoms group
     //! \return The charged sidechain potential between two groups
     double calculate_contribution_improved(ResidueFB *res1, ResidueFB *res2, int m_i, int m_j, int sc_atom_index_i, int sc_atom_index_j) {
          
          // Import protein definitions (such as residue names)
          using namespace definitions;
          
          double m_ij = 0.0;
          double c_ij = 0.0;
          double E_sc = 0.0;
          
          m_ij = 1.5 * m_i * m_j;
          
          const std::vector<AtomEnum> & n_i = all_sidechain_charged_atoms[sc_atom_index_i];
          const std::vector<AtomEnum> & n_j = all_sidechain_charged_atoms[sc_atom_index_j];

          
          int ni = n_i.size();
          int nj = n_j.size();

          double GAMMA_ij = calculate_gamma(n_i, n_j, res1, res2);

          c_ij = std::min(1.0, GAMMA_ij / (ni + nj));
          E_sc += m_ij * c_ij;
          
          return E_sc;
     }
};

//! PROFASI side chain charged interaction energy term 
class TermProfasiSidechainCharge: public TermProfasiSidechainChargeBase<TermProfasiSidechainCharge> {
     
     //! For convenience, define local base class
     typedef phaistos::TermProfasiSidechainChargeBase<TermProfasiSidechainCharge> TermProfasiSidechainChargeBase;

public:
     
     //! Use same settings as base class
     typedef TermProfasiSidechainChargeBase::Settings Settings;

     //! Settings
     const Settings settings;
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiSidechainCharge(ChainFB *chain,
                                const Settings &settings=Settings(),
                                RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiSidechainChargeBase(chain, "profasi-sidechain-charge", settings, random_number_engine),
            settings(settings) {
          
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiSidechainCharge(const TermProfasiSidechainCharge &other, 
                                RandomNumberEngine *random_number_engine,
                                int thread_index, ChainFB *chain)
          : TermProfasiSidechainChargeBase(other, random_number_engine, thread_index, chain),
            settings(other.settings) {}
     
     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return charged sidechain potential energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {
          
          //! Import protein definitions (such as residue names)
          using namespace definitions;

          double energy = 0.0;
          int i = 0;
          int j = 0;
          
          // Go through all the residue pairs 
          for (ResidueIterator<ChainFB> it1(*this->chain); !it1.end(); ++it1) {
               for (ResidueIterator<ChainFB> it2(it1+1); !it2.end(); ++it2) {
                    
                    ResidueFB *residue1 = &*it1;
                    ResidueFB *residue2 = &*it2;
                    
                    // Get charged sidechain atoms group according to the residue type
                    const std::vector<AtomEnum> & atoms_sc1 = this->side_chain_charged_atoms[residue1->residue_type];
                    const std::vector<AtomEnum> & atoms_sc2 = this->side_chain_charged_atoms[residue2->residue_type];

                    // Get charged sidechain atoms group on the terminals
                    const std::vector<AtomEnum> & atoms_nterm = this->terminal_charged_atoms[0];
                    const std::vector<AtomEnum> & atoms_cterm = this->terminal_charged_atoms[1];
                    
                    if ((i == 0) &&
                        residue1->terminal_status == NTERM &&
                        (residue1->residue_type == ARG ||
                         residue1->residue_type == LYS ||
                         residue1->residue_type == ASP ||
                         residue1->residue_type == GLU)) {//residue with charged atoms are at the N-termainal
                         
                         double GAMMA_Nj = this->calculate_gamma(atoms_nterm, atoms_sc1, residue1, residue1);
                         
                         double m_ij = 1.5 * this->settings.m[residue1->residue_type];
                         double c_ij = std::min(1.0, GAMMA_Nj / (atoms_nterm.size() + atoms_sc1.size()));
                         energy += m_ij * c_ij;
                         i++;
               }
                    
                    if ((j == 0 ) &&
                        residue2->terminal_status == CTERM &&
                        (residue2->residue_type == ARG ||
                         residue2->residue_type == LYS ||
                         residue2->residue_type == ASP ||
                         residue2->residue_type == GLU)) {//residue with charged atoms are at the C-termainal
                         
                         double GAMMA_Nj = this->calculate_gamma(atoms_cterm, atoms_sc2, residue2, residue2);
                         
                         double m_ij = -1.5 * this->settings.m[residue2->residue_type];
                         double c_ij = std::min(1.0, GAMMA_Nj / (atoms_cterm.size() + atoms_sc2.size()));
                         energy += m_ij * c_ij;
                         j++;
                    }
                    
                    energy += this->calculate_contribution(residue1, residue2);
               }
          }
          
          return this->profasi_energy_in_kcal_per_mol * energy;
          
     }
     
};


//! PROFASI SidechainCharge energy term - cached version
class TermProfasiSidechainChargeCached: public TermProfasiSidechainChargeBase<TermProfasiSidechainChargeCached> {
     
     //! For convenience, define local base class
     typedef phaistos::TermProfasiSidechainChargeBase<TermProfasiSidechainChargeCached> TermProfasiSidechainChargeBase;
     
     //! Specifies, for each node, which atoms contribute to the energy term
     std::vector<std::pair<ResidueFB*,std::vector<definitions::AtomEnum> > > active_atoms;

     //! Specified the last iteration that an update was made for a given residue1 index
     std::vector<long int> observation_vector1;

     //! Specified the last iteration that an update was made for a given residue2 index
     std::vector<long int> observation_vector2;

     //! Specified the last iteration that an update was made for a given residue1,residue2 pair     
     std::vector<std::vector<long int> > observation_matrix;

     //! Keeps track of current iterations
     long int iterations_counter;

     //! For all chaintree nodes at the lowest level, find the atoms that contribute
     void find_active_atoms() {
          
          //! Import protein definitions (such as residue names)
          using namespace definitions;
          
          // The chaintree level (level 0 is omitted for efficiency reasons)
          int level = 1;
          active_atoms = std::vector<std::pair<ResidueFB*,std::vector<AtomEnum> > >(this->chain->get_chain_tree()->get_nodes_at_level(level),
                                                                                    std::make_pair((ResidueFB*)NULL,std::vector<AtomEnum>()));
          int offset = this->chain->get_chain_tree()->get_level_start_index(level);
          for (unsigned int i=0; i<active_atoms.size(); ++i) {
               
               NodeType *node = this->chain->get_chain_tree()->nodes[i+offset];
               
               // Merge atoms from the two child nodes
               std::vector<Atom*> merged_child_atoms = node->child1->atoms;
               if (node->child2) {
                    merged_child_atoms.insert(merged_child_atoms.end(),
                                              node->child2->atoms.begin(),
                                              node->child2->atoms.end());
               }
               
               // Iterate over all child node atoms and detect active atoms
               for (unsigned int j=0; j<merged_child_atoms.size(); ++j) {
                    Atom *atom = merged_child_atoms[j];
                    ResidueFB *residue = &(*this->chain)[atom->residue->index];
                    AtomEnum atom_type = atom->atom_type;
                    ResidueEnum residue_type = residue->residue_type;

                    std::vector<AtomEnum> residue_atoms = this->side_chain_charged_atoms[residue_type];
                    std::vector<AtomEnum> residue_atoms_nterm =
                         this->terminal_charged_atoms[0];
                    std::vector<AtomEnum> residue_atoms_cterm =
                         this->terminal_charged_atoms[1];
                    
                    
                    for (unsigned int k=0; k< residue_atoms.size(); ++k) {
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
                    
                    if (residue->terminal_status == NTERM) {
                         for (unsigned int k = 0; k < residue_atoms_nterm.size(); ++k) {
                              
                              if (residue_atoms_nterm[k] == atom_type) {
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
                    
                    else if (residue->terminal_status == CTERM) {
                         for (unsigned int k = 0; k < residue_atoms_cterm.size(); ++k) {
                              
                              if (residue_atoms_cterm[k] == atom_type) {
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
     }
     
public:
     
     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;
     
     //! Cached iterator (used when iterating over node pairs)
     CachedIterator<chaintree::PairIterator<ChainFB,NodeType,NodeType> > cached_it;
     
     //! Cached iterator - for local part of energy term
     CachedIterator<ResidueIterator<ChainFB> > cached_it_local;
     
     //! Use same settings as base class
     typedef TermProfasiSidechainChargeBase::Settings Settings;

     //! Local Settings object
     const Settings settings;
     
     //! Chaintree iterator settings - described in detail in chaintree.h
     chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings iterator_settings;
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiSidechainChargeCached(ChainFB *chain,
                                      const Settings &settings=Settings(),
                                      RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiSidechainChargeBase(chain, "profasi-sidechain-charge-cached", settings, random_number_engine),
            observation_vector1(std::vector<long int>(chain->size(), 0)),
            observation_vector2(std::vector<long int>(chain->size(), 0)),
            observation_matrix(std::vector<std::vector<long int> >(chain->size(), std::vector<long int>(chain->size(), 0))),
            iterations_counter(0),
            cached_it(*chain),
            cached_it_local(*chain),
            settings(settings) {
               
          // Only evaluate modified pairs (as always when caching)
          bool only_modified_pairs = true;
          int minimum_residue_distance = 1;
          bool omit_lowest_level = true;
          
          iterator_settings = chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings(settings.r_cut,
                                                                                           only_modified_pairs,
                                                                                           minimum_residue_distance,
                                                                                           omit_lowest_level);
          
          find_active_atoms();
          
          // Initialize local cache
          cached_it_local(*this->chain, 0, this->chain->size());
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiSidechainChargeCached(const TermProfasiSidechainChargeCached &other, 
                                      RandomNumberEngine *random_number_engine,
                                      int thread_index, ChainFB *chain)
          : TermProfasiSidechainChargeBase(other, random_number_engine, thread_index, chain),
            observation_vector1(other.observation_vector1),
            observation_vector2(other.observation_vector2),
            observation_matrix(other.observation_matrix),
            iterations_counter(other.iterations_counter),
            cached_it(*chain),
            cached_it_local(*chain),
            settings(other.settings),
            iterator_settings(other.iterator_settings) {
          
          find_active_atoms();
          
          // Initialize local cache
          cached_it_local(*this->chain, 0, this->chain->size());
     }
     
     //! Evaluate chain energy(cached version)
     //! \param move_info Object containing information about last move
     //! \return Charged sidechain potential energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {
          
          //! Import protein definitions (such as residue names)
          using namespace definitions;
          
          iterations_counter++;
          
          bool calculate_n_term = true;
          bool calculate_c_term = true;
          
          if (move_info) {
               if (move_info->modified_angles_start > 0) {
                    calculate_n_term = false;
               }
               if (move_info->modified_angles_end < this->chain->size()) {
                    calculate_c_term = false;
               }
          }
          
          if (calculate_n_term) {
               
               ResidueFB *res_n_term = &(*this->chain)[0];
               
               if((res_n_term->residue_type == ARG ||
                   res_n_term->residue_type == LYS ||
                   res_n_term->residue_type == ASP ||
                   res_n_term->residue_type == GLU)&&
                  (observation_vector1[res_n_term->index] < iterations_counter)){//residue with charged atoms are at the N-termainal
                    
                    const std::vector<AtomEnum> & sc1_atoms = this->side_chain_charged_atoms[res_n_term->residue_type];
                    const std::vector<AtomEnum> & nterm_atoms = this->terminal_charged_atoms[0];
                    
                    double GAMMA_NTERM = this->calculate_gamma(nterm_atoms, sc1_atoms, res_n_term, res_n_term);
                    double m_ij = 1.5 * this->settings.m[res_n_term->residue_type];
                    double c_ij = std::min(1.0, GAMMA_NTERM / (nterm_atoms.size() + sc1_atoms.size()));
                    
                    double energy_n_term = m_ij * c_ij;
                    observation_vector1[res_n_term->index] = iterations_counter;
                    
                    cached_it_local = ResidueIterator<ChainFB>(*this->chain, 0, 1);
                    cached_it_local.register_contribution(energy_n_term);
               }
               
          }
          
          if (calculate_c_term) {
               ResidueFB *res_c_term = &(*this->chain)[this->chain->size()-1];
               
               if((res_c_term->residue_type == ARG ||
                   res_c_term->residue_type == LYS ||
                   res_c_term->residue_type == ASP ||
                   res_c_term->residue_type == GLU) &&
                  (observation_vector2[res_c_term->index] < iterations_counter)){//residue with charged atoms are at the N-termainal
                    
                    const std::vector<AtomEnum> & sc2_atoms = this->side_chain_charged_atoms[res_c_term->residue_type];
                    const std::vector<AtomEnum> & cterm_atoms = this->terminal_charged_atoms[1];
                    
                    double GAMMA_CTERM = this->calculate_gamma(cterm_atoms, sc2_atoms, res_c_term, res_c_term);
                    double m_ij = -1.5 * this->settings.m[res_c_term->residue_type];
                    double c_ij = std::min(1.0, GAMMA_CTERM / (cterm_atoms.size() + sc2_atoms.size()));
                    
                    double energy_c_term = m_ij * c_ij;
                    
                    observation_vector2[res_c_term->index] = iterations_counter;
                    
                    cached_it_local = ResidueIterator<ChainFB>(*this->chain,
                                                                  this->chain->size()-1, this->chain->size());
                    cached_it_local.register_contribution(energy_c_term);
               }
          }
          
          // Cached iterator evaluation
          for (cached_it(*this->chain, iterator_settings); !cached_it.end(); ++cached_it) {
               
               double contribution = 0.0;
               NodeType *node1 = cached_it->first;
               NodeType *node2 = cached_it->second;
               
               ResidueFB *res_1 = active_atoms[node1->index].first;
               ResidueFB *res_2 = active_atoms[node2->index].first;
               
               if (res_1 && res_2) {
                    
                    if (observation_matrix[res_1->index][res_2->index] < iterations_counter) {
                         contribution += this->calculate_contribution(res_1, res_2);
                    }
                    
                    observation_vector1[res_1->index] = iterations_counter;
                    observation_vector2[res_2->index] = iterations_counter;
                    observation_matrix[res_1->index][res_2->index] = iterations_counter;
               }
               
               cached_it.register_contribution(contribution);
          }
          
          double energy_non_local = cached_it.compute_total();
          
          double energy_local = cached_it_local.compute_total();
          
          return this->profasi_energy_in_kcal_per_mol * (energy_non_local + energy_local);
     }
     
     //! Accept last energy evaluation
     void accept() {
          cached_it.accept();
          cached_it_local.accept();
     }
     
     //! Reject last energy evaluation
     void reject() {
          cached_it.reject();
          cached_it_local.reject();
     }
     
};

//! PROFASI sidechain charge energy term - improved version
class TermProfasiSidechainChargeImproved: public TermProfasiSidechainChargeBase<TermProfasiSidechainChargeImproved> {
     
     //! For convenience, define local base class
     typedef phaistos::TermProfasiSidechainChargeBase<TermProfasiSidechainChargeImproved> TermProfasiSidechainChargeBase;
      
     //! Specifies, for each node, which atoms contribute to the energy term
     std::vector<std::pair<ResidueFB*,std::vector<definitions::AtomEnum> > > active_atoms;
     
     //! For all chaintree nodes at the lowest level, find the atoms that contribute
     std::vector< ResidueFB* > sc_residues;

     //! Keeps track of charge sign for all residues
     std::vector< int > charge_sign;

     //! The index of the first charged sidechain atoms group for each residue in sc_residues
     std::vector< int > all_sidechain_charge_atoms_index;

     //! For all residue on the chain, find the atoms that contribute
     void find_sidechain_charge_residues() {
          
          using namespace definitions;
          
          int num_sc_residue = 0;
          
          // Go though all residues
          for (ResidueIterator<ChainFB> it1(*this->chain); !it1.end(); ++it1) {
               
               if((it1->terminal_status == NTERM && it1->residue_type!=PRO) || it1->terminal_status == CTERM )
                    num_sc_residue += 1;
               
               if(it1->residue_type == ARG ||
                  it1->residue_type == LYS ||
                  it1->residue_type == ASP ||
                  it1->residue_type == GLU  )
                    num_sc_residue += 1;
          }
          
          sc_residues.resize(num_sc_residue);
          charge_sign.resize(num_sc_residue);
          all_sidechain_charge_atoms_index.resize(num_sc_residue);
         
          int i_sc_res = 0;
          
          // Register all the charged sidechain atoms of the chain
          for (ResidueIterator<ChainFB> it2(*this->chain); !it2.end(); ++it2) {
               
               ResidueFB* res = &(*it2);
               
               if(it2->terminal_status == NTERM && it2->residue_type != PRO){
                    sc_residues[i_sc_res] = res;
                    all_sidechain_charge_atoms_index[i_sc_res] = 20;
                    charge_sign[i_sc_res] = 1;
                    i_sc_res += 1;
               }
               if(it2->terminal_status == CTERM){
                    sc_residues[i_sc_res] = res;
                    all_sidechain_charge_atoms_index[i_sc_res] = 21;
                    charge_sign[i_sc_res] = -1;
                    i_sc_res += 1;
               }
               if(it2->residue_type == ARG || it2->residue_type == LYS){
                    sc_residues[i_sc_res] = res;
                    all_sidechain_charge_atoms_index[i_sc_res] = it2->residue_type;
                    charge_sign[i_sc_res] = 1;
                    i_sc_res += 1;
               }
               if(it2->residue_type == ASP || it2->residue_type == GLU ){
                    sc_residues[i_sc_res] = res;
                    all_sidechain_charge_atoms_index[i_sc_res] = it2->residue_type;
                    charge_sign[i_sc_res] = -1;
                    i_sc_res += 1;
               }
          }
     }
     
public:
     
     //! Use same settings as base class
     typedef TermProfasiSidechainChargeBase::Settings Settings;
     
     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;
     
     //! Local Settings object
     const Settings settings;
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiSidechainChargeImproved(ChainFB *chain,
                                       const Settings &settings=Settings(),
                                        RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiSidechainChargeBase(chain, "profasi-sidechain-charge-improved", settings, random_number_engine),
            settings(settings){

               find_sidechain_charge_residues();
          }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiSidechainChargeImproved(const TermProfasiSidechainChargeImproved &other, 
                                        RandomNumberEngine *random_number_engine,
                                        int thread_index, ChainFB *chain)
          : TermProfasiSidechainChargeBase(other, random_number_engine, thread_index, chain),
            settings(other.settings){
          find_sidechain_charge_residues();
     }
     
     //! Evaluate chain energy
     //! \param move_info Object containing information about last move
     //! \return Charged sidechain potential energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {
          
          using namespace definitions;
          
          double energy = 0;
        
          // Go through all the charged sidechain group pairs
          for(unsigned int i=0; i< sc_residues.size(); i++ ){
               for(unsigned int j = i+1; j< sc_residues.size(); j++ ){
                    
                    ResidueFB *res1 = sc_residues[i];
                    ResidueFB *res2 = sc_residues[j];
      
                    int sign_i = charge_sign[i];
                    int sign_j = charge_sign[j];
                    
                    int sc_atom_index_i = all_sidechain_charge_atoms_index[i];
                    int sc_atom_index_j = all_sidechain_charge_atoms_index[j];
                    
                    energy += this->calculate_contribution_improved(res1, res2, sign_i, sign_j, sc_atom_index_i, sc_atom_index_j);
                    
               }
          } 
          
          return this->profasi_energy_in_kcal_per_mol * energy;      
          
     }
     
};

}

#endif
