// term_angle.h --- Generate angle distribution historgrams
// Copyright (C) 2008-2013 Jan B. Valentin, Wouter Boomsma
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


#ifndef ANGLE_H
#define ANGLE_H

#include "../energy/energy_term.h"
#include "../protein/iterators/dof_iterator.h"

namespace phaistos {

// Angle - angular statistics of chain
template <typename CHAIN_TYPE>
class TermAngle: public EnergyTermCommon<TermAngle<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermAngle<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;               

protected:

     std::vector<double> angles;

     DofIterator<ChainFB> begin;
     DofIterator<ChainFB> end;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          // Number of bins in histogram
          bool dihedral_only;

          Settings(bool dihedral_only=true)
               : dihedral_only(dihedral_only) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "dihedral-only:" << settings.dihedral_only << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }

     } settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermAngle(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "angle", settings, random_number_engine),
            begin((*chain)(0, definitions::CA),definitions::DIHEDRAL,
                  settings.dihedral_only?
                  (DofIterator<ChainFB>::DIHEDRAL_DOFS+DofIterator<ChainFB>::CHI_ANGLES)
                  :
                  (DofIterator<ChainFB>::STANDARD_DOFS+DofIterator<ChainFB>::CHI_ANGLES+
                   DofIterator<ChainFB>::N_DIHEDRAL+
                   DofIterator<ChainFB>::NTERM_CA_DIHEDRAL+DofIterator<ChainFB>::CTERM_N_DIHEDRAL)
            ),
            end((Atom*)NULL),
            settings(settings) {
            

          // Initialize angle bins
          for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
               angles.push_back(0.0);
          }
     }     

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermAngle(const TermAngle &other, 
                        RandomNumberEngine *random_number_engine,
                        int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            angles(other.angles),
            begin((*chain)(0, definitions::CA),definitions::DIHEDRAL,
                  settings.dihedral_only?
                  (DofIterator<ChainFB>::DIHEDRAL_DOFS+DofIterator<ChainFB>::CHI_ANGLES)
                  :
                  (DofIterator<ChainFB>::STANDARD_DOFS+DofIterator<ChainFB>::CHI_ANGLES+
                   DofIterator<ChainFB>::N_DIHEDRAL+
                   DofIterator<ChainFB>::NTERM_CA_DIHEDRAL+DofIterator<ChainFB>::CTERM_N_DIHEDRAL)
            ),
            end((Atom*)NULL),
            settings(other.settings) {
            
     }

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {
          int counter = 0;
          for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
               if (it.get_dof_type() == definitions::ANGLE) {
                    double value = *it;
                    while (value < -2*M_PI)
                         value += 2*M_PI;
                    while (value > 2*M_PI)
                         value -= 2*M_PI;
                    if (value < 0)
                         value *= -1;
                    if (value > M_PI)
                         value = 2*M_PI-value;
                    angles[counter] = value;
                    counter++;
               } else {
                    double value = *it+M_PI;
                    while (value < 0)
                         value += 2*M_PI;
                    while (value >= 2*M_PI)
                         value -= 2*M_PI;
                    angles[counter] = value;
                    counter++;
               }
          }

          // This energy term has no meaningful return value
          return 0;
     }
};

//! Observable specialization for TermAngle
template <typename CHAIN_TYPE>
class Observable<TermAngle<CHAIN_TYPE> >: public TermAngle<CHAIN_TYPE>, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermAngle<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}          

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<typename TermAngle<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     
     //! Constructor.
     //! \param energy_term Angle energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermAngle<CHAIN_TYPE> &energy_term, 
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermAngle<CHAIN_TYPE>(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
            
            std::string header = "(";
            for (DofIterator<ChainFB> it=this->begin; it!=this->end; ++it) {
                    std::string header_entry = ("(" + 
                                               boost::lexical_cast<std::string>(it.get_atom()->residue->index) + "_" + 
                                               boost::lexical_cast<std::string>(it.get_atom()->residue->residue_type) + "_" +
                                               boost::lexical_cast<std::string>(it.get_atom()->atom_type) + ")[" + 
                                               boost::lexical_cast<std::string>(it.get_dof_type()) + "]");
                    if (header != "(")
                         header += ",";
                    header += header_entry;
            }
            this->name = "angle"+header+")";
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index, typename TermAngle<CHAIN_TYPE>::ChainType *chain)
          : TermAngle<CHAIN_TYPE>(other, thread_index, chain),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermAngle<CHAIN_TYPE> *clone(int thread_index=0, typename TermAngle<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermAngle<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {

          this->evaluate();
          
          std::string output = "";
          // strip output-target for output-mode
          std::string ot = settings.output_target;
          ot = ot.substr(0, ot.find('#'));

          if (ot != "pdb-b-factor") {
              char buffer[128];
              for (unsigned int i = 0; i < this->angles.size(); i++) {
                  sprintf(buffer," %2.4f", this->angles[i]);
                  output+= buffer;
              }
          }
          
          return output;
     }
     
};

}

#endif
