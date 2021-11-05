// Proline phi angle torsional energy(Gauss distribution) 
// Copyright (C) 2010-2011 Pengfei Tian, Wouter Boomsma, Jesper Ferkinghoff-Borg
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


#ifndef PROFASI_PROLINE_PHI_TORSION_H
#define PROFASI_PROLINE_PHI_TORSION_H

#include "energy/energy_term.h"
#include "protein/iterators/residue_iterator.h"
#include "utils/vector_matrix_3d.h"

namespace phaistos {
     
//! PROFASI proline phi torional energy term - base class
template<typename DERIVED_CLASS>
class TermProfasiProlinePhiTorsionBase: public EnergyTermCommon<DERIVED_CLASS, ChainFB>,
                                        public EnergyTermCommonProfasi {
private:
     
     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;
     
protected:
     
     //! Internal parameter equal to 1/2/beta/phi_dev/phi_dev
     double internal_param; 
     
public:
     //! Local setting class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
          
     public:
          
          //! Mean value for phi angle
          double mean_phi;

          //! Std deviation value for phi angle
          double dev_phi ;

          //! kT
          double inv_beta;

          //! Constructor. Defines default values for settings object
          Settings(double mean_phi = -1.1338,
                   double dev_phi = 0.1905,
                   double inv_beta = 0.596163)
               : mean_phi(mean_phi),
                 dev_phi(dev_phi),
                 inv_beta(inv_beta) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "mean-phi:" << settings.mean_phi << "\n";
               o << "dev-phi:" << settings.dev_phi << "\n";
               o << "inv-beta:" << settings.inv_beta << "\n";
               o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }          

     } settings;    //!< Local settings object 
     
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param name Energy term name
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiProlinePhiTorsionBase(ChainFB *chain,
                                      std::string name,
                                      const Settings &settings=Settings(),
                                      RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine) {

          internal_param = 0.5 / settings.dev_phi / settings.dev_phi * settings.inv_beta;
     }
     
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiProlinePhiTorsionBase(const TermProfasiProlinePhiTorsionBase &other,
                                      RandomNumberEngine *random_number_engine,
                                      int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            internal_param(other.internal_param) {}
     
     
     //! Calculate torsional energy
     //! \param x Proline phi angle
     //! \return Proline phi torsional angle potential energy in kcal/mol
     double calculate_torsion_energy(double x) {
          double eng = (x - settings.mean_phi)*(x - settings.mean_phi)*internal_param;;
          return eng; 
     }
     
};



//! PROFASI proline phi torsional energy term
class TermProfasiProlinePhiTorsion: public TermProfasiProlinePhiTorsionBase<TermProfasiProlinePhiTorsion> {
     
     //! For convenience, define local base class
     typedef phaistos::TermProfasiProlinePhiTorsionBase<TermProfasiProlinePhiTorsion> TermProfasiProlinePhiTorsionBase;
     
     //! CA atoms of Proline residues
     std::vector< Atom* > proline_ca;
     
     //! Proline torsion value
     double proline_torsion;

     //! Proline torsion value backup
     double proline_torsion_backup;
     
     //! Save all the CA atoms of Proline in a vector
     void find_proline_ca(){
          
          //Import protein definitions (such as residue names)
          using namespace definitions;
          
          int num_pro = 0;
          
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){
               if(it->residue_type == PRO)
                    num_pro += 1;
          }
          
          proline_ca.resize(num_pro);
          
          int ipro=0;
          
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){
               if(it->residue_type == PRO)
                    proline_ca[ipro++] = (*it)[CA];
          }
     }
     
     
public:
     
     //! Use same settings as base class
     typedef TermProfasiProlinePhiTorsionBase::Settings Settings;

     //! Local Settings object
     const Settings settings;
              
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermProfasiProlinePhiTorsion(ChainFB *chain,
                                  const Settings &settings=Settings(),
                                  RandomNumberEngine *random_number_engine = &random_global)
          : TermProfasiProlinePhiTorsionBase(chain, "profasi-proline-phi-torsion", settings, random_number_engine),
            proline_torsion(UNINITIALIZED),
            proline_torsion_backup(UNINITIALIZED),
            settings(settings) {
          
          find_proline_ca();
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermProfasiProlinePhiTorsion(const TermProfasiProlinePhiTorsion &other, 
                                  RandomNumberEngine *random_number_engine,
                                  int thread_index, ChainFB *chain)
          : TermProfasiProlinePhiTorsionBase(other, random_number_engine, thread_index, chain),
            proline_torsion(other.proline_torsion),
            proline_torsion_backup(other.proline_torsion_backup),
            settings(other.settings) {
          
          find_proline_ca();
     }
     
     //! Evaluate chain energy
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {
        
          //! Import protein definitions (such as residue names)   
          using namespace definitions;
          
          proline_torsion = 0;
          
          double phi = 0.0;
          
          if (!is_initialized(proline_torsion_backup) ||
              (move_info && move_info->move_type != definitions::SIDECHAIN)){
               //(move_info->modified_angles_start > 0 && move_info->modified_angles_end < (this->chain->size()-1)))) {
               
               for (unsigned int i=0; i< proline_ca.size(); ++i){
                    phi = proline_ca[i]->get_dihedral();  //get phi angle of Proline
                    proline_torsion += this->calculate_torsion_energy(phi);
               }
          } else{
               proline_torsion = proline_torsion_backup;
          }
          
          return proline_torsion;     
     }
     
     //! Accept last energy evaluation
     void accept() {
          proline_torsion_backup = proline_torsion;
     }
};
     
}

#endif
