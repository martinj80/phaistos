// observable_xtc_trajectory.h --- Output Gromacs .xtc trajectory data
// Copyright (C) 2008-2012 Wouter Boomsma
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


#ifndef OBSERVABLE_XTC_TRAJECTORY_H
#define OBSERVABLE_XTC_TRAJECTORY_H

// #include <boost/filesystem.hpp>

// extern "C" {
// #include "xdrfile.h"
// #include "xdrfile_xtc.h"
// }

#include "protein/xtc_chain.h"

namespace phaistos {

//! Energy term for output trajectory information
//! This class only serves as an empty placeholder
//! for the specification of the corresponding
//! observable
template <typename CHAIN_TYPE>
class TermXtcTrajectory {};


//! Observable specialization for TermXtcTrajectory
//! Note: inheriting from IndependentObservable excludes this
//! observable from the normal I/O handling in the observable collection
template <typename CHAIN_TYPE>
class Observable<TermXtcTrajectory<CHAIN_TYPE> >: public EnergyTermCommon<Observable<TermXtcTrajectory<CHAIN_TYPE> >, CHAIN_TYPE>, 
                                                  public TermXtcTrajectory<CHAIN_TYPE>, 
                                                  public ObservableBase,
                                                  public IndependentObservable {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<Observable<TermXtcTrajectory<CHAIN_TYPE> >,CHAIN_TYPE> EnergyTermCommon;               

     //! Box matrix used by xtc
     matrix box;

     //! Log file name to which xtc content it dumped
     std::string log_filename;

     // std::string log_filename_backup;

     //! Initialization function
     void initialize(int thread_index) {

          log_filename = generate_log_filename(settings.output_target, thread_index);

          // size_t dot_position = log_filename.find_last_of(".");
          // log_filename_backup = log_filename.substr(0,dot_position) + "_backup" + log_filename.substr(dot_position);

          for(unsigned i=0; i<DIM; i++) {
               for(unsigned j=0; j<DIM; j++) {
                    box[i][j] = 0.0;
               }
          }
     }

public:

     //! Local settings class.
     // const class Settings: public TermXtcTrajectory<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     const class Settings: public EnergyTermCommon::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){
               // Override default value for output target
               this->output_target = "trajectory_%p_%t.xtc";
          }          

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<ObservableBase::Settings>(settings);
               return o;
          }

     } settings; //!< Local settings object 

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     Observable(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "xtc-trajectory", settings, random_number_engine),
            settings(settings) {}
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, 
                RandomNumberEngine *random_number_engine,
                int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {}

     //! Copy Constructor - assign reference function
     //! \param other Source object from which copy is made
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const Observable &other, 
                const ObservableBase::Settings &settings,
                Energy<CHAIN_TYPE> *reference_energy_function)
          : EnergyTermCommon(other),
            settings(dynamic_cast<const Settings&>(settings)) {
          initialize(other.thread_index);
     }     

     //! Make observation.
     //! \param move_info Object containing information about the last executed move     
     //! \param current_iteration Index of current iteration
     //! \param register_only Whether this call should only register data (not output it)
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {

          // if (boost::filesystem::exists(log_filename_backup))
          //     boost::filesystem::remove(log_filename_backup);
          // boost::filesystem::copy_file(log_filename, log_filename_backup);
          XDRFILE *xd = xdrfile_open(log_filename.c_str(),"a");

          int step = current_iteration;
          float time = current_iteration;
          float prec = 1000;
          xtc_write_chain(*this->chain, xd, step, time, prec);

          xdrfile_close(xd);          

          // boost::filesystem::remove(log_filename_backup);

          return "";
     }
     
};

}

#endif
