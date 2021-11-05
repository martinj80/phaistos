// helix_content.h --- Measure the percentage of the chain that is in a helix conformation
// Copyright (C) 2008-2011 Wouter Boomsma
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

#ifndef HELIX_CONTENT_H
#define HELIX_CONTENT_H

#include "../energy/energy_term.h"
#include "observable.h"


namespace phaistos {

//! Measure the percentage of the chain that is in a helix conformation
template <typename CHAIN_TYPE>
class TermHelixContent: public EnergyTermCommon<TermHelixContent<CHAIN_TYPE>, CHAIN_TYPE> {
protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermHelixContent<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;               

protected:
     
     // Whether residues are in a helical state
     std::vector<bool> helical_state;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {

          //! Phi angle minimum boundary
          static double default_min_angle1(ChainFB *chain){return -90.0/180.0*M_PI;}
          //! Phi angle maximum boundary
          static double default_max_angle1(ChainFB *chain){return -30.0/180.0*M_PI;}
          //! Psi angle minimum boundary
          static double default_min_angle2(ChainFB *chain){return -77.0/180.0*M_PI;}
          //! Psi angle maximum boundary
          static double default_max_angle2(ChainFB *chain){return -17.0/180.0*M_PI;}

          //! Theta angle minimum boundary          
          static double default_min_angle1(ChainCA *chain){return UNINITIALIZED;}
          //! Theta angle maximum boundary          
          static double default_max_angle1(ChainCA *chain){return UNINITIALIZED;}
          //! Tau angle minimum boundary          
          static double default_min_angle2(ChainCA *chain){return UNINITIALIZED;}
          //! Tau angle maximum boundary          
          static double default_max_angle2(ChainCA *chain){return UNINITIALIZED;}

     public:

          //! Minimum angle1 used as boundary for helix region
          double min_angle1;

          //! Maximum angle1 used as boundary for helix region
          double max_angle1;

          //! Minimum angle2 used as boundary for helix region
          double min_angle2;

          //! Maximum angle2 used as boundary for helix region
          double max_angle2;

          //! Constructor
          Settings(double min_angle1 = default_min_angle1((CHAIN_TYPE*)NULL),
                   double max_angle1 = default_max_angle1((CHAIN_TYPE*)NULL),
                   double min_angle2 = default_min_angle2((CHAIN_TYPE*)NULL),
                   double max_angle2 = default_max_angle2((CHAIN_TYPE*)NULL))
               : min_angle1(min_angle1), max_angle1(max_angle1),
                 min_angle2(min_angle2), max_angle2(max_angle2) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "min-angle1:" << settings.min_angle1 << "\n";
               o << "max-angle1:" << settings.max_angle1 << "\n";
               o << "min-angle2:" << settings.min_angle2 << "\n";
               o << "max-angle2:" << settings.max_angle2 << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }

     } settings;    //!< Local settings object 

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermHelixContent(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                      RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "helix-content", settings, random_number_engine) {}


     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermHelixContent(const TermHelixContent &other, 
                      RandomNumberEngine *random_number_engine,
                      int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {}


     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {

          helical_state.resize(this->chain->size(), false);

          double helical_content = 0;
          int normalize_counter = 0;
          for (ResidueIterator<CHAIN_TYPE> it(*this->chain, 1, this->chain->size()-1); !it.end(); ++it) {

               double angle1 = it->get_angles()[0];
               double angle2 = it->get_angles()[1];

               if ((angle1 > settings.min_angle1) && (angle1 < settings.max_angle1) &&
                   (angle2 > settings.min_angle2) && (angle2 < settings.max_angle2)) {
                    helical_state[it->index] = true;
                    helical_content += 1;
               } 

               normalize_counter++;
          }
          helical_content /= (double)normalize_counter;

          return helical_content;          

     }

};


//! Observable specialization for TermHelixContent
template <typename CHAIN_TYPE>
class Observable<TermHelixContent<CHAIN_TYPE> >: public TermHelixContent<CHAIN_TYPE>, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermHelixContent<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Whether to split up observable into a per-residue vector
          bool per_residue;

          //! Constructor. Defines default values for settings object.
          Settings(bool per_residue=false)
               : per_residue(per_residue){}          

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "per-residue:" << settings.per_residue << "\n";
               o << static_cast<typename TermHelixContent<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     
     //! Constructor.
     //! \param energy_term HelixContent energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermHelixContent<CHAIN_TYPE> &energy_term, 
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermHelixContent<CHAIN_TYPE>(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index, typename TermHelixContent<CHAIN_TYPE>::ChainType *chain)
          : TermHelixContent<CHAIN_TYPE>(other, thread_index, chain),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermHelixContent<CHAIN_TYPE> *clone(int thread_index=0, typename TermHelixContent<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermHelixContent<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {
          this->evaluate_weighted(move_info);
          if (settings.per_residue) {

               std::string output = "";
               for (ResidueIterator<CHAIN_TYPE> it(*this->chain, 1, this->chain->size()-1); !it.end(); ++it) {
                    std::string output_entry = (ObservableBase::vector_output_tag(*it) + 
                                                boost::lexical_cast<std::string>(this->helical_state[it->index]));
                    if (output != "")
                         output += ",";
                    output += output_entry;
               }

               return output;
          } else {
               return boost::lexical_cast<std::string>(TermHelixContent<CHAIN_TYPE>::evaluate_weighted(move_info));
          }
     }

};


}

#endif


