// observable_git.h --- Output GIT vector
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


#ifndef OBSERVABLE_GIT_H
#define OBSERVABLE_GIT_H

namespace phaistos {

//! Energy term for output trajectory information
//! This class only serves as an empty placeholder
//! for the specification of the corresponding
//! observable
template <typename CHAIN_TYPE>
class TermGit {};


//! Observable specialization for TermGit
//! Note: inheriting from IndependentObservable excludes this
//! observable from the normal I/O handling in the observable collection
template <typename CHAIN_TYPE>
class Observable<TermGit<CHAIN_TYPE> >: public EnergyTermCommon<Observable<TermGit<CHAIN_TYPE> >, CHAIN_TYPE>, 
                                        public TermGit<CHAIN_TYPE>, 
                                        public ObservableBase,
                                        public IndependentObservable {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<Observable<TermGit<CHAIN_TYPE> >,CHAIN_TYPE> EnergyTermCommon;               

     //! Reference to main energy function
     Energy<CHAIN_TYPE> *reference_energy_function;

     //! Log file name to which xtc content it dumped
     std::string log_filename_template;

     //! Whether to dump all PDB output to a single file
     bool single_file_output;

     // Random number generator
     RandomNumberEngine *random_number_engine;

     //! Git object
     Git git;
     

     //! Initialization function
     void initialize(int thread_index) {

          std::string output_target = settings.output_target;
          output_target = generate_log_filename(output_target, thread_index);

          boost::filesystem::path log_file_path(output_target);
          boost::filesystem::path directory_path = log_file_path.parent_path();

          if (directory_path.string() != "") {
               boost::filesystem::create_directory(directory_path.string());
          }

          log_filename_template = output_target;
          
          // If %i and %e do not occur, we dump a single file for the entire simulation
          single_file_output = true;
          if ((settings.output_target.find("%i") != std::string::npos) || 
              (settings.output_target.find("%e") != std::string::npos)) {
               single_file_output = false;
          } else {
               git.set_output_file(log_filename_template);
          }
          
     }

public:

     //! Local settings class.
     // const class Settings: public TermGit<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     const class Settings: public EnergyTermCommon::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){
               // Override default value for output target
               this->output_target = "samples_%p_%t/sample_%i_%e.gitvec";
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
          : EnergyTermCommon(chain, "git", settings, random_number_engine),
            random_number_engine(random_number_engine),
            git(random_number_engine),
            settings(settings) {}
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, 
                RandomNumberEngine *random_number_engine,
                int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            random_number_engine(random_number_engine),
            git(random_number_engine) {}

     //! Copy Constructor - assign reference function
     //! \param other Source object from which copy is made
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const Observable &other, 
                const ObservableBase::Settings &settings,
                Energy<CHAIN_TYPE> *reference_energy_function)
          : EnergyTermCommon(other),
            reference_energy_function(reference_energy_function),
            random_number_engine(other.random_number_engine),
            git(other.random_number_engine),
            settings(dynamic_cast<const Settings&>(settings)) {
          initialize(other.thread_index);
     }     

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {

          // Extract energy from reference function (has already been evaluated)
          double energy = settings.weight*this->reference_energy_function->sum;

          // Create output filename
          std::string log_filename = generate_log_filename(log_filename_template, energy, current_iteration);

          // Open file in case we are outputting to multiple files
          if (!single_file_output) {
               git.set_output_file(log_filename);               
          }

          git.generate_gauss_integrals(*this->chain, generate_log_filename("sample_%i_%t_%e",energy,current_iteration,this->thread_index));

          return "";
     }
     
};

}

#endif
