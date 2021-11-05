// observable_pdb.h --- Output PDB data
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


#ifndef OBSERVABLE_PDB_H
#define OBSERVABLE_PDB_H

#include <boost/filesystem.hpp>

namespace phaistos {

//! Energy term for output trajectory information
//! This class only serves as an empty placeholder
//! for the specification of the corresponding
//! observable
//! Note: inheriting from IndependentObservable
//! excludes this term from the energy collections.
template <typename CHAIN_TYPE>
class TermPdb{};


//! Observable specialization for TermPdb
//! Note: inheriting from IndependentObservable excludes this
//! observable from the normal I/O handling in the observable collection
template <typename CHAIN_TYPE>
class Observable<TermPdb<CHAIN_TYPE> >: public EnergyTermCommon<Observable<TermPdb<CHAIN_TYPE> >, CHAIN_TYPE>, 
                                        public TermPdb<CHAIN_TYPE>, 
                                        public ObservableBase,
                                        public IndependentObservable {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<Observable<TermPdb<CHAIN_TYPE> >,CHAIN_TYPE> EnergyTermCommon;               

     //! Reference to main energy function
     Energy<CHAIN_TYPE> *reference_energy_function;

     //! Log file name to which xtc content it dumped
     std::string log_filename_template;

     //! Whether to dump all PDB output to a single file
     bool single_file_output;

     //! Output stream 
     std::ofstream pdb_file;

     //! model index
     int model_index;

     //! Current iteration index
     PHAISTOS_LONG_LONG current_iteration;

     //! Previous index at which structure was dumped (in minimum energy mode)
     PHAISTOS_LONG_LONG dump_index_previous;

     //! Minimum energy observed (in minimum energy mode)
     double energy_min;

     //! Maximum energy observed (in minimum energy mode)
     double energy_max;


     //! Initialization function
     void initialize(int thread_index) {

          dump_index_previous = 0;
          energy_min = std::numeric_limits<double>::infinity();
          energy_max = -std::numeric_limits<double>::infinity();

          std::string output_target = settings.output_target;
          output_target = generate_log_filename(output_target, thread_index);

          boost::filesystem::path log_file_path(output_target);
          boost::filesystem::path directory_path = log_file_path.parent_path();

          if (directory_path.string() != "") {
               boost::filesystem::create_directory(directory_path.string());
          }
#if BOOST_VERSION / 100 % 1000 >= 46
          // boost version 1.46 or larger (boost::filesystem v3)
          std::string log_file = log_file_path.filename().string();
#else
          // for older versions of boost
          std::string log_file = log_file_path.filename();
#endif
          if (settings.minimum_energy_mode)
               log_file = "min_" + log_file;

          log_file_path = directory_path / boost::filesystem::path(log_file);
          log_filename_template = log_file_path.string();

          // If %i and %e do not occur, we dump a single file for the entire simulation
          single_file_output = true;
          if ((settings.output_target.find("%i") != std::string::npos) || 
              (settings.output_target.find("%e") != std::string::npos)) {
               single_file_output = false;
          } else {
               std::string log_filename = generate_log_filename(log_filename_template);
               pdb_file.open(log_filename.c_str());
          }

          model_index = 1;
     }

public:

     //! Local settings class.
     // const class Settings: public TermPdb<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     const class Settings: public EnergyTermCommon::Settings, public ObservableBase::Settings {
     public:

          //! Whether to keep track of minimum energy structures
          bool minimum_energy_mode;

          //! How far from the minimum energy a structure must be before being dumped
          double minimum_energy_fraction;

          //! Minimum interval between dumped structures
          int minimum_energy_interval;

          //! Constructor. Defines default values for settings object.
          Settings(bool minimum_energy_mode=false,
                   double minimum_energy_fraction=0.1,
                   int minimum_energy_interval=500)
               : minimum_energy_mode(minimum_energy_mode),
                 minimum_energy_fraction(minimum_energy_fraction),
                 minimum_energy_interval(minimum_energy_interval) {
               // Override default value for output target
               // this->output_target = "sample_%012lld_%d_%04f.pdb";
               this->output_target = "samples_%p_%t/sample_%i_%e.pdb";
          }          

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "minimum-energy-mode:" << settings.minimum_energy_mode << "\n";               
               o << "minimum-energy-fraction:" << settings.minimum_energy_fraction << "\n";               
               o << "minimum-energy-interval:" << settings.minimum_energy_interval << "\n";               
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
          : EnergyTermCommon(chain, "pdb", settings, random_number_engine),
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
            reference_energy_function(reference_energy_function),
            settings(dynamic_cast<const Settings&>(settings)) {
          initialize(other.thread_index);
     }     


     //! Destructor
     ~Observable() {
          if (single_file_output) {
               pdb_file.close();
          }
     }

     //! Make observation.
     //! Since the header in PDB output is potentially dependent on the output of other
     //! observables, all functionality is moved to the finalize function.
     //! \param move_info Object containing information about the last executed move     
     //! \param current_iteration Index of current iteration
     //! \param register_only Whether this call should only register data (not output it)
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {

          this->current_iteration = current_iteration;

          return "";
     }

     //! This method is called once all observables have been evaluated, allowing
     //! for communication between them.
     //! In the case of ObservablePdb, all functionality is located here
     virtual void observe_finalize(ObservableCollectionBase *observable_collection) {

          // Extract energy from reference function (has already been evaluated)
          double energy = settings.weight*this->reference_energy_function->sum;

          // Create output filename
          std::string log_filename = generate_log_filename(log_filename_template, energy, current_iteration);

          // In minimum energy mode, we only dump if the energy is below a certain cutoff
          bool dump_structure = true;
          if (settings.minimum_energy_mode) {

               // Default to not dumping
               dump_structure = false;

               // Update maximum energy value
               if (energy>energy_max) {
                    energy_max=energy;
               }

               // Update minimum energy value
               if (energy<energy_min) {
                    energy_min=energy;
               }

               double delta = (energy_max-energy_min)+1e-5;

               if (((current_iteration - dump_index_previous) >= settings.minimum_energy_interval) &&
                   ((std::fabs(energy-energy_min) / delta) <
                    settings.minimum_energy_fraction)) {

                    dump_structure = true;
                    dump_index_previous = current_iteration;
               }
          }

          // Abort if we are not dumpin structure
          if (!dump_structure)
               return;

          // Add energy output to header
          std::string header = "";          
          if (!this->reference_energy_function->terms.empty())
               header += "Energies:\n" + boost::lexical_cast<std::string>(this->reference_energy_function->get_data()) + "\n";

          // Add observable output to header
          header += observable_collection->
               gather_string_streams(ObservableCollectionBase::PDB_HEADER);

          // Add observable output to b factors
          std::string b_factor_string = observable_collection->
               gather_string_streams(ObservableCollectionBase::PDB_B_FACTOR);

          // Open file in case we are outputting to multiple files
          if (!single_file_output) {
               pdb_file.open(log_filename.c_str());               
          }

          // Output
          CHAIN_TYPE *superimpose_chain = NULL;
          int begin_offset = 0;
          int end_offset = 0;
          pdb_file << this->chain->output_as_pdb(superimpose_chain, 
                                                 begin_offset, 
                                                 end_offset, 
                                                 model_index, 
                                                 header, 
                                                 b_factor_string);
          pdb_file.flush();

          // Close file in case of multiple files. Otherwise, prepare
          // new section in existing pdb file
          if (!single_file_output) {
               pdb_file.close();
          } else {
               pdb_file << "END\nENDMDL\n";
               model_index++;
          }

     }
     
};

}

#endif
