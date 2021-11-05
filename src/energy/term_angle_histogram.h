// observable_angle_hist.h --- Generate angle distribution historgrams
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


#ifndef ANGLE_HISTOGRAM_H
#define ANGLE_HISTOGRAM_H

#include "../energy/energy_term.h"
#include "../protein/iterators/dof_iterator.h"

namespace phaistos {

// Angle Histogram - angular statistics of chain
template <typename CHAIN_TYPE>
class TermAngleHistogram: public EnergyTermCommon<TermAngleHistogram<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermAngleHistogram<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;               

protected:

     std::vector<std::vector<PHAISTOS_LONG_LONG> > angle_bins;

     DofIterator<ChainFB> begin;
     DofIterator<ChainFB> end;

public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          // Number of bins in histogram
          unsigned int bins;

          Settings(unsigned int bins=128)
               : bins(bins) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "bins:" << settings.bins << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }

     } settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermAngleHistogram(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "angle-histogram", settings, random_number_engine),
            begin((*chain)(0, definitions::CA),
                  definitions::DIHEDRAL,
                  (DofIterator<ChainFB>::STANDARD_DOFS+DofIterator<ChainFB>::CHI_ANGLES+
                   DofIterator<ChainFB>::N_DIHEDRAL+
                   DofIterator<ChainFB>::NTERM_CA_DIHEDRAL+DofIterator<ChainFB>::CTERM_N_DIHEDRAL)),
            end((Atom*)NULL),
            settings(settings) {

          // Initialize angle bins
          for (DofIterator<ChainFB> it=begin; it!=end; ++it) {
               angle_bins.push_back(std::vector<PHAISTOS_LONG_LONG>(settings.bins, 0));
          }
     }     

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermAngleHistogram(const TermAngleHistogram &other, 
                        RandomNumberEngine *random_number_engine,
                        int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            angle_bins(other.angle_bins),
            begin((*chain)(0, definitions::CA),
                  definitions::DIHEDRAL,
                  (DofIterator<ChainFB>::STANDARD_DOFS+DofIterator<ChainFB>::CHI_ANGLES+
                   DofIterator<ChainFB>::N_DIHEDRAL+
                   DofIterator<ChainFB>::NTERM_CA_DIHEDRAL+DofIterator<ChainFB>::CTERM_N_DIHEDRAL)),
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
                    int bin = settings.bins/2+(int)floor(((value)/(M_PI))*settings.bins/2 + 0.5);
                    angle_bins[counter][bin]++;
                    counter++;
               } else {
                    double value = *it+M_PI;
                    while (value < 0)
                         value += 2*M_PI;
                    while (value >= 2*M_PI)
                         value -= 2*M_PI;
                    int bin = static_cast<int>(floor(((value)/(2*M_PI))*settings.bins) + 0.5);
                    angle_bins[counter][bin]++;
                    counter++;
               }
          }

          // This energy term has no meaningful return value
          return 0;
     }
};

//! Observable specialization for TermAngleHistogram
template <typename CHAIN_TYPE>
class Observable<TermAngleHistogram<CHAIN_TYPE> >: public TermAngleHistogram<CHAIN_TYPE>, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermAngleHistogram<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}          

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<typename TermAngleHistogram<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     
     //! Constructor.
     //! \param energy_term AngleHistogram energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermAngleHistogram<CHAIN_TYPE> &energy_term, 
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermAngleHistogram<CHAIN_TYPE>(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index, typename TermAngleHistogram<CHAIN_TYPE>::ChainType *chain)
          : TermAngleHistogram<CHAIN_TYPE>(other, thread_index, chain),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermAngleHistogram<CHAIN_TYPE> *clone(int thread_index=0, typename TermAngleHistogram<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermAngleHistogram<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {

          this->evaluate_weighted(move_info);

          std::string output = "";
          if (!register_only) {
               int counter = 0;
               for (DofIterator<ChainFB> it=this->begin; it!=this->end; ++it) {
                    std::string output_entry = ("(" + 
                                                boost::lexical_cast<std::string>(it.get_atom()->residue->index) + "_" + 
                                                boost::lexical_cast<std::string>(it.get_atom()->residue->residue_type) + "_" +
                                                boost::lexical_cast<std::string>(it.get_atom()->atom_type) + ")[" + 
                                                boost::lexical_cast<std::string>(it.get_dof_type()) + "]" + 
                                                boost::lexical_cast<std::string>(this->angle_bins[counter]));
                    if (output != "")
                         output += ",";
                    output += output_entry;

                    ++counter;
               }
          }
          return output;
     }
     
};

}

#endif
