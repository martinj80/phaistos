// profasi_hydrogen_bond --- PROFASI backbone-backbone and sidechain-backbone hydrogen bonding energy term
// Copyright (C) 2010-2011 Pengfei Tian, Wouter Boomsma, Jes Frellsen, Jesper Ferkinghoff-Borg
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


#ifndef PROFASI_HYDROGEN_BOND_H
#define PROFASI_HYDROGEN_BOND_H

#include "energy/energy_term.h"
#include "profasi_energy_term_common.h"
#include "protein/iterators/residue_iterator.h"
#include "protein/iterators/pair_iterator_chaintree.h"
#include "utils/vector_matrix_3d.h"
#include <cmath>

namespace phaistos {

//! PROFASI hydrogen bond energy term - base class
template<typename DERIVED_CLASS>
class TermProfasiHydrogenBondBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB>,
                                   public EnergyTermCommonProfasi {
private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

     //@{
     //! Internal variables
     double bhb;
     double ahb;
     //@}

protected:

     //! Vector of sidechain donors
     std::vector<std::vector<definitions::AtomEnum> > hb_sidechain_hatoms;

     //! Vector of sidechain acceptors
     std::vector<std::vector<definitions::AtomEnum> > hb_sidechain_oatoms;

public:

     //! Local settings class.     
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Strength parameter of backbone-backbone bonds
          double e_hb1;

          //! Strength parameter of sidechain-backbone bonds
          double e_hb2;

          //! Value of parameter sigma
          double sigma_hb;

          //! Cut off distanve
          double r_cut;

          //! Boolean indicating whether to use ideal length of C_H and C_O bond
          bool use_ideal_distances;

          //! Constructor          
          Settings(double e_hb1 = 3.0,
                   double e_hb2 = 2.3,
                   double sigma_hb = 2.0,
                   double r_cut = 4.5,
                   bool use_ideal_distances=true)
                 : e_hb1(e_hb1),
                   e_hb2(e_hb2),
                   sigma_hb(sigma_hb),
                   r_cut(r_cut),
                   use_ideal_distances(use_ideal_distances) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "e-hb1:" << settings.e_hb1 << "\n";
               o << "e-hb2:" << settings.e_hb2 << "\n";
               o << "sigma-hb:" << settings.sigma_hb << "\n";
               o << "r-cut:" << settings.r_cut << "\n";
               o << "use-ideal-distances:" << settings.use_ideal_distances << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }          

     } settings;     //!< Local settings object 


     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrogenBondBase(ChainFB *chain,
                                 std::string name,
                                 const Settings &settings=Settings(),
                                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          //! Import phaistos::make_vector namespace
          using namespace vector_utils;

	  bhb = -30*(pow(settings.sigma_hb / settings.r_cut, 10) -
		             pow(settings.sigma_hb / settings.r_cut, 12)) / pow(settings.r_cut, 2);
	  ahb = -(5*pow(settings.sigma_hb / settings.r_cut, 12) -
		             6*pow(settings.sigma_hb / settings.r_cut, 10)) - bhb*pow(settings.r_cut, 2);

          // Register sidechain donors
	  hb_sidechain_hatoms.resize(21);
	  hb_sidechain_hatoms[ARG] = make_vector( HE, HH11 , HH12 , HH22, HH21);
	  hb_sidechain_hatoms[LYS] = make_vector( HZ1 , HZ2 , HZ3);
	 
          // Index 20 used for N-terminal donor
	  hb_sidechain_hatoms[20]  = make_vector( H1 , H2 , H3);

          // Register sidechian acceptors
	  hb_sidechain_oatoms.resize(21);
	  hb_sidechain_oatoms[ASP] = make_vector( OD1, OD2 );
	  hb_sidechain_oatoms[GLU] = make_vector( OE1, OE2 );
	  
          // Index 20 used for C-terminal acceptor
	  hb_sidechain_oatoms[20]  = make_vector( O , OXT );
     }


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrogenBondBase(const TermProfasiHydrogenBondBase &other,
                                 RandomNumberEngine *random_number_engine,
                                 int thread_index, ChainFB *chain)
	  : EnergyTermCommon(other, random_number_engine, thread_index, chain),
	    bhb(other.bhb),
	    ahb(other.ahb),
	    hb_sidechain_hatoms(other.hb_sidechain_hatoms),
	    hb_sidechain_oatoms(other.hb_sidechain_oatoms)
	    {}


     //! Calculate hydrogen bond energy for given donor and acceptor
     //! \param atom_Hdon Dontor atom in the hydrogen bonding
     //! \param atom_Oacc Acceptor atom in the hydrogen bonding
     //! \return Hydrogen bond protential between the NH and CO groups
     inline double hbpotential( Atom * atom_Hdon, Atom * atom_Oacc) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          Atom *atom_Ndon = atom_Hdon->get_neighbour<-1>(BACKBONE);
          Atom *atom_Cacc = atom_Oacc->get_neighbour<-1>(BACKBONE);

          // Squared distances
          double r2_OH = (atom_Hdon->position - atom_Oacc->position).norm_squared();
          double r2_NH = (atom_Hdon->position - atom_Ndon->position).norm_squared();
          double r2_CO = (atom_Cacc->position - atom_Oacc->position).norm_squared();

          // Scalar product corresponding NHO angle (alpha)
          double dot_prod_nh_oh = ((atom_Hdon->position - atom_Ndon->position) *
                                   (atom_Hdon->position - atom_Oacc->position));

          // Scalar product corresponding HOC angle (beta)
          double dot_prod_oc_oh = ((atom_Cacc->position - atom_Oacc->position) *
                                   (atom_Hdon->position - atom_Oacc->position));

          // Variables represent function U and V 
          double u_r, v;

          if (settings.use_ideal_distances) {


               // In order to be consistent with the result with profasi, use the cutoff of the distance of NH and CO
               double r_CO_cut = 1.23;

               //double r_NH_cut = 1.0;

               if(atom_Oacc->atom_type == OE1 ||
                  atom_Oacc->atom_type == OE2 ||
                  atom_Oacc->atom_type == OD1 ||
                  atom_Oacc->atom_type == OD2) {

                    // Cut off of distance for OE1 OE2 in GLU, OD1 OD2 in ASP
                    r_CO_cut = 1.25; 
               }

               if (r2_OH > (settings.r_cut*settings.r_cut))
                    return 0;

               if (dot_prod_nh_oh > 0 || dot_prod_oc_oh > 0)
                    return 0;
               
               // Internal parameter, used to speed up the calculation
               double r2_internal = settings.sigma_hb * settings.sigma_hb / r2_OH;
               double r4_internal = r2_internal *r2_internal;
               double r8_internal = r4_internal *r4_internal;

               u_r = 5 * r4_internal * r8_internal - 6 * r2_internal * r8_internal ;
               v = pow(r2_NH * r2_CO / (r_CO_cut*r_CO_cut) * (dot_prod_nh_oh * dot_prod_nh_oh * dot_prod_oc_oh * dot_prod_oc_oh)/(r2_OH*r2_NH*r2_OH*r2_CO), 0.25);

          } else {

               if (r2_OH > (settings.r_cut*settings.r_cut))
                    return 0;

               if (dot_prod_nh_oh > 0 || dot_prod_oc_oh > 0)
                    return 0;

               double r2_internal = settings.sigma_hb * settings.sigma_hb / r2_OH;
               double r4_internal = r2_internal *r2_internal;
               double r8_internal = r4_internal *r4_internal;

               u_r = 5 * r4_internal * r8_internal -6 * r2_internal * r8_internal ;

               v = pow((dot_prod_nh_oh * dot_prod_nh_oh * dot_prod_oc_oh * dot_prod_oc_oh)/(r2_OH*r2_NH*r2_OH*r2_CO), 0.25);

          }

          return v * (u_r + ahb+bhb*r2_OH);
     }

};


//! PROFASI hydrogen bond energy term
class TermProfasiHydrogenBond: public TermProfasiHydrogenBondBase<TermProfasiHydrogenBond> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiHydrogenBondBase<TermProfasiHydrogenBond> TermProfasiHydrogenBondBase;

public:

     //! Use same settings as base class
     typedef TermProfasiHydrogenBondBase::Settings Settings;

     //! Settings
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrogenBond(ChainFB *chain,
                             const Settings &settings=Settings(),
                             RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiHydrogenBondBase(chain, "profasi-hydrogen-bond", settings, random_number_engine),
            settings(settings) {
          
     }
     

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrogenBond(const TermProfasiHydrogenBond &other, 
                             RandomNumberEngine *random_number_engine,
                             int thread_index, ChainFB *chain)
          : TermProfasiHydrogenBondBase(other, random_number_engine, thread_index, chain),
            settings(other.settings) {}

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     //! \return energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          // Hydrogen bond protential of backbone-backbone bonds 
          double HBMM = 0;

          // Hydrogen bond protential of backbone-sidechain bonds 
          double HBMS = 0;

          // Go through all the residues, it1 provides donors and it2 provides acceptors
          for (ResidueIterator<ChainFB> it1(*this->chain); !it1.end(); ++it1){
               for (ResidueIterator<ChainFB> it2(*this->chain); !it2.end(); ++it2){

                    if(it1->residue_type != PRO) {

                         if(it2->terminal_status != CTERM) {

                              const std::vector<AtomEnum> & res_H = this->hb_sidechain_hatoms[it1->residue_type];

                              // Choose sidechain donor and backbone acceptor
                              for(unsigned int i=0; i < res_H.size(); i++){

                                   Atom *atom_Hdonor = (*it1)[res_H[i]];
                                   Atom *atom_Oterm = (*it2)[O];

                                   HBMS += this->hbpotential(atom_Hdonor, atom_Oterm);
                              }

                              // Choose backbone acceptor and backbone donor
                              if(it1->terminal_status != NTERM){

                                   if( (it2->index - it1->index) >  0 || (it1->index - it2->index) >  2 ) {

                                        Atom *atom_Hdonor = (*it1)[H];
                                        Atom *atom_Oacceptor = (*it2)[O];

                                        HBMM += this->hbpotential(atom_Hdonor, atom_Oacceptor);
//                                        if(this->hbpotential(atom_Hdonor, atom_Oacceptor) !=0){
//                                             std::cout<<it1->index<<" "<<it1->residue_type<<" "<<atom_Hdonor->atom_type<<" "<<it2->index<<" "<<it2->residue_type<<" "<<atom_Oacceptor->atom_type<<"  "<<this->hbpotential(atom_Hdonor, atom_Oacceptor)<<"\n";
//                                        }
                                   }
                              }

                              //Choose N terminal donor and backbone acceptor
                              if(it1->terminal_status == NTERM){

                                   // Index 20 used for N-terminal donor
                                   const std::vector<AtomEnum> & res_H = this->hb_sidechain_hatoms[20];

                                   for(unsigned int i=0; i < res_H.size(); i++){

                                        Atom *atom_Hterm = (*it1)[res_H[i]];
                                        Atom *atom_Oacceptor = (*it2)[O];

                                        HBMS += this->hbpotential(atom_Hterm, atom_Oacceptor);
                                   }

                              }

                         }

                         if(it1->terminal_status != NTERM){

                              const std::vector<AtomEnum> & res_O = this->hb_sidechain_oatoms[it2->residue_type];

                              // Choose sidechain acceptor and backbone donor
                              for(unsigned int i = 0; i < res_O.size(); i++){

                                   Atom *atom_Hdonor = (*it1)[H];
                                   Atom *atom_Oterm = (*it2)[res_O[i]];

                                   HBMS += this->hbpotential(atom_Hdonor, atom_Oterm);
                              }

                              //Choose C terminal acceptor and backbone donor
                              if(it2->terminal_status == CTERM){

                                   const std::vector<AtomEnum> & res_O = this->hb_sidechain_oatoms[20];

                                   for(unsigned int i = 0; i < res_O.size(); i++){

                                        Atom *atom_Hdonor = (*it1)[H];
                                        Atom *atom_Oterm = (*it2)[res_O[i]];

                                        HBMS += this->hbpotential(atom_Hdonor, atom_Oterm);
                                   }
                              }
                         }
                    }
               }
          }

          HBMM *= this->settings.e_hb1;
          HBMS *= this->settings.e_hb2;

       return this->profasi_energy_in_kcal_per_mol * (HBMM + HBMS);

     }


};

//! PROFASI hydrogen bond energy term - cached version
class TermProfasiHydrogenBondCached: public TermProfasiHydrogenBondBase<TermProfasiHydrogenBondCached> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiHydrogenBondBase<TermProfasiHydrogenBondCached> TermProfasiHydrogenBondBase;

public:
     //! Define NodeType locally for ease of reference
     typedef ChainFB::ChainTree::NodeType NodeType;

     //! Cached iterator (used when iterating over node pairs)
     CachedIterator<chaintree::PairIterator<ChainFB,NodeType,NodeType> > cached_it;

     //! Use same settings as base class
     typedef TermProfasiHydrogenBondBase::Settings Settings;

     //! Settings object
     const Settings settings;

     //! Chaintree iterator settings - described in detail in chaintree.h
     chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings iterator_settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object     
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrogenBondCached(ChainFB *chain,
                                   const Settings &settings=Settings(),
                                   RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiHydrogenBondBase(chain, "profasi-hydrogen-bond-cached", settings, random_number_engine),
         cached_it(*chain),
         settings(settings) {

          //! Only evaluate modified pairs (as always when caching)
          bool only_modified_pairs = true;

          iterator_settings = chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings(settings.r_cut,
                                                                                           only_modified_pairs);
          //iterator_settings = typename chaintree::PairIterator<ChainFB,NodeType,NodeType>::Settings(std::numeric_limits<double>::infinity(),
                                          //                    only_modified_pairs);
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrogenBondCached(const TermProfasiHydrogenBondCached &other, 
                                   RandomNumberEngine *random_number_engine,
                                   int thread_index, ChainFB *chain)
          : TermProfasiHydrogenBondBase(other, random_number_engine, thread_index, chain),
            cached_it(*chain),
            settings(other.settings),
            iterator_settings(other.iterator_settings) {
     }

     //! Evaluate hydrogen bond potential produced by the acceptors and donors on the two nodes
     //! \param node1 The first node in the hydrogen bonding who can provide acceptors or donors  
     //! \param node2 The second node in the hydrogen bonding who can provide acceptors or donors
     //! \return hydrogen bond potential energy 
     double calculate_contribution(NodeType *node1, NodeType *node2) {

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          // Hydrogen bond protential of backbone-backbone bonds 
          double HBMM = 0;

          // Hydrogen bond protential of backbone-sidechain bonds 
          double HBMS = 0;

          ResidueFB *res1 = node1->frame->res;

          // Prolines have no hydrogen bond donor
          if(res1->residue_type != PRO) {

               // If the current node has an oxygen donor, but is not a C-terminal
               if(node2->has_atom(O) && !node2->has_atom(OXT)) {

                    Atom *atom_Oacceptor = (*node2)[O];

                    const std::vector<AtomEnum> & res_H = this->hb_sidechain_hatoms[res1->residue_type];

                    // Choose sidechain donor and backbone acceptor
                    for(unsigned int i=0; i < res_H.size(); i++){

                         if (node1->has_atom(res_H[i])) {
                              Atom *atom_Hdonor = (*node1)[res_H[i]];
                              HBMS += this->hbpotential(atom_Hdonor, atom_Oacceptor);
                         }
                    }

                    // Choose backbone acceptor and backbone donor
                    if(res1->terminal_status != NTERM){

                         if (node1->has_atom(H)) {
                              int hydrogen_res_index = res1->index;
                              int oxygen_res_index = node2->frame->res->index-1;

                              if( (oxygen_res_index - hydrogen_res_index) >  0 || (hydrogen_res_index - oxygen_res_index) >  2 ) {
                                   Atom *atom_Hdonor = (*node1)[H];
                                   HBMM += this->hbpotential(atom_Hdonor, atom_Oacceptor);

                              }
                         }

                         // Choose N terminal donor and backbone acceptor
                    } else {

                         // Index 20 used for N-terminal donor
                         const std::vector<AtomEnum> & res_H = this->hb_sidechain_hatoms[20];

                         for(unsigned int i=0; i < res_H.size(); i++){

                              if (node1->has_atom(res_H[i])) {
                                   Atom *atom_Hterm = (*node1)[res_H[i]];
                                   HBMS += this->hbpotential(atom_Hterm, atom_Oacceptor);
                              }
                         }

                    }

               }

               if(node1->has_atom(H) && res1->terminal_status != NTERM){

                    Atom *atom_Hdonor = (*node1)[H];

                    const std::vector<AtomEnum> & res_O = this->hb_sidechain_oatoms[node2->frame->res->residue_type];

                    // Choose sidechain acceptor and backbone donor
                    for(unsigned int i = 0; i < res_O.size(); i++){

                         if (node2->has_atom(res_O[i])) {
                              Atom *atom_Oterm = (*node2)[res_O[i]];
                              HBMS += this->hbpotential(atom_Hdonor, atom_Oterm);
                         }
                    }

                    // Choose C-terminal acceptor and backbone donor
                    if(node2->has_atom(OXT)) {

                         const std::vector<AtomEnum> & res_O = this->hb_sidechain_oatoms[20];

                         for(unsigned int i = 0; i < res_O.size(); i++){
                              Atom *atom_Oterm = (*node2)[res_O[i]];
                              HBMS += this->hbpotential(atom_Hdonor, atom_Oterm);
                         }
                    }
               }
          }

          HBMM *= this->settings.e_hb1;
          HBMS *= this->settings.e_hb2;

          return (HBMM + HBMS);

     }

     //! Evaluate chain energy
     //! \param move_info Object containing information about the last executed move
     //! \return energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          // Cached iterator evaluation
          for (cached_it(*this->chain, iterator_settings); !cached_it.end(); ++cached_it) {
               double energy = 0.0;

               // Node containing hydrogen bond donors
               NodeType *node1 = cached_it->first;

               // Node containing hydrogen bond acceptors
               NodeType *node2 = cached_it->second;

               energy += calculate_contribution(node1, node2);
               energy += calculate_contribution(node2, node1);

               cached_it.register_contribution(energy);
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

//! PROFASI hydrogen bond energy term - improved version
class TermProfasiHydrogenBondImproved: public TermProfasiHydrogenBondBase<TermProfasiHydrogenBondImproved> {

     //! For convenience, define local base class
     typedef phaistos::TermProfasiHydrogenBondBase<TermProfasiHydrogenBondImproved> TermProfasiHydrogenBondBase;

     //! Put backbone donor in the vector
     std::vector< Atom* > bb_don_atoms;

     //! Put backbone acceptor in the vector
     std::vector< Atom* > bb_acc_atoms;

     //! Put sidechain donor in the vector
     std::vector< Atom* > sc_don_atoms;

     //! Put sidechain acceptor in the vector
     std::vector< Atom* > sc_acc_atoms;

     //! Hydrogen bond - mainchain-mainchain contribution
     double HBMM;

     //! Hydrogen bond - mainchain-sidechain contribution
     double HBMS;

     //! Hydrogen bond - mainchain-mainchain contribution backup
     double HBMM_backup;

     //! Hydrogen bond - mainchain-sidechain contribution backup
     double HBMS_backup;

     //! Put acceptors and donors in vectors separately
     void find_hydrogenbond_atoms(){

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          int num_bbdon = 0;
          int num_bbacc = this->chain->size()-1;
          int num_scdon = 0;
          int num_scacc = 0;

          // Count the number of acceptors and donors
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){

              switch (it->residue_type) {
                  case ASP:
                  case GLU:
                      num_scacc+=2;
                      break;
                  case LYS:
                      num_scdon+=3;
                      break;
                  case ARG:
                      num_scdon+=5;
                      break;
                  default:
                      break;
              };

              if (it->terminal_status == NTERM && it->residue_type != PRO){
                   num_scdon += 3;
              }

              if(it->terminal_status == CTERM){
                   num_scacc += 2;
              }

              if(it->terminal_status != NTERM && it->residue_type != PRO ){
                   num_bbdon += 1;
              }

          }

          bb_don_atoms.resize(num_bbdon);
          bb_acc_atoms.resize(num_bbacc);
          sc_don_atoms.resize(num_scdon);
          sc_acc_atoms.resize(num_scacc);

          int ibbdon=0;
          int ibbacc=0;
          int iscdon=0;
          int iscacc=0;

          // Put the acceptros and donors in the vectors
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){

               if(it->residue_type!= PRO && it->terminal_status != NTERM){

                    bb_don_atoms[ibbdon++] = (*it)[H];
               }

               if(it->terminal_status != CTERM){

                    bb_acc_atoms[ibbacc++] = (*it)[O];
               }

               if(it->residue_type == ARG){
                    sc_don_atoms[iscdon++] = (*it)[HE];
                    sc_don_atoms[iscdon++] = (*it)[HH11];
                    sc_don_atoms[iscdon++] = (*it)[HH12];
                    sc_don_atoms[iscdon++] = (*it)[HH22];
                    sc_don_atoms[iscdon++] = (*it)[HH21];
               }
               else if(it->residue_type == LYS){
                    sc_don_atoms[iscdon++] = (*it)[HZ1];
                    sc_don_atoms[iscdon++] = (*it)[HZ2];
                    sc_don_atoms[iscdon++] = (*it)[HZ3];
               }
               else if(it->residue_type == ASP){
                    sc_acc_atoms[iscacc++] = (*it)[OD1];
                    sc_acc_atoms[iscacc++] = (*it)[OD2];
               }
               else if(it->residue_type == GLU){
                    sc_acc_atoms[iscacc++] = (*it)[OE1];
                    sc_acc_atoms[iscacc++] = (*it)[OE2];
               }

               if(it->residue_type != PRO && it->terminal_status == NTERM){
                    sc_don_atoms[iscdon++] = (*it)[H1];
                    sc_don_atoms[iscdon++] = (*it)[H2];
                    sc_don_atoms[iscdon++] = (*it)[H3];
               }

               if(it->terminal_status == CTERM){
                    sc_acc_atoms[iscacc++] = (*it)[O];
                    sc_acc_atoms[iscacc++] = (*it)[OXT];
               }

          }
     }

public:

     //! Use same settings as base class
     typedef TermProfasiHydrogenBondBase::Settings Settings;

     //! Local Settings object
     const Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiHydrogenBondImproved(ChainFB *chain,
                                     const Settings &settings=Settings(),
                                     RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiHydrogenBondBase(chain, "profasi-hydrogen-bond-improved", settings, random_number_engine),
            HBMM(UNINITIALIZED),
            HBMS(UNINITIALIZED),
            HBMM_backup(UNINITIALIZED),
            HBMS_backup(UNINITIALIZED),
            settings(settings) {

          find_hydrogenbond_atoms();

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiHydrogenBondImproved(const TermProfasiHydrogenBondImproved &other, 
                                     RandomNumberEngine *random_number_engine,
                                     int thread_index, ChainFB *chain)
          : TermProfasiHydrogenBondBase(other, random_number_engine, thread_index, chain),
            HBMM(other.HBMM),
            HBMS(other.HBMS),
            HBMM_backup(other.HBMM_backup),
            HBMS_backup(other.HBMS_backup),
            settings(other.settings) {

          find_hydrogenbond_atoms();
     }

     
     //! Check whether any of the moved residues involve potential hydrogen bonds
     bool hydrogen_bonding(MoveInfo *move_info){

          // Import protein definitions (such as residue names)
          using namespace definitions;

          bool sc_participate = false;
          
          for (int i=move_info->modified_angles_start; i< move_info->modified_angles_end; i++) {
               ResidueEnum re = (*this->chain)[i].residue_type;
               if (re == ARG ||
                   re == LYS||
                   re == ASP||
                   re == GLU) {
                    sc_participate = true;
               }
               }
          return sc_participate;
     }
     
     //! Evaluate chain energy
     //! \param move_info Object containing information about the last executed move
     //! \return energy in kcal/mol
     double evaluate(MoveInfo *move_info=NULL) {

          HBMM = 0;
          HBMS = 0;

          if (!is_initialized(HBMM_backup) || (!move_info) ||
              (move_info->move_type != definitions::SIDECHAIN)){
               //(move_info->modified_angles_start > 0 && move_info->modified_angles_end < (this->chain->size()-1)))) {
               for(unsigned int i=0; i< bb_don_atoms.size(); i++ ){
                    for(unsigned int j = 0; j< bb_acc_atoms.size(); j++){

                         Residue *res_don = bb_don_atoms[i]->residue;
                         Residue *res_acc = bb_acc_atoms[j]->residue;
                         
                         if((res_acc->index - res_don->index) > 0 || (res_don->index - res_acc->index) > 2 ){
                              HBMM += this->hbpotential(bb_don_atoms[i], bb_acc_atoms[j]);
                              
                         }
                         
                    }
               }
               HBMM *= this->settings.e_hb1;
          } else {
               HBMM = HBMM_backup;
          }
          

          if (!is_initialized(HBMS_backup) || (!move_info) || !(move_info->move_type == definitions::SIDECHAIN && !hydrogen_bonding(move_info))) {
               
                    for(unsigned int i=0; i< bb_don_atoms.size(); i++ ) {
                         for(unsigned int j = 0; j< sc_acc_atoms.size(); j++) {
                              HBMS += this->hbpotential(bb_don_atoms[i], sc_acc_atoms[j]);
                         }
                    }
                    
                    for(unsigned int i=0; i< sc_don_atoms.size(); i++ ) {
                         for(unsigned int j = 0; j< bb_acc_atoms.size(); j++) {
                              HBMS += this->hbpotential(sc_don_atoms[i], bb_acc_atoms[j]);
                         }
                    }
                    
                    HBMS *= this->settings.e_hb2;

          } else {

               HBMS = HBMS_backup;
          }

          return this->profasi_energy_in_kcal_per_mol * (HBMM + HBMS);
     }


     //! Accept last energy evaluation
     void accept() {
          HBMM_backup = HBMM;
          HBMS_backup = HBMS;
     }
};

}

#endif
