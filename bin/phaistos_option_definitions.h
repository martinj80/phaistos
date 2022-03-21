// phaistos_option_definitions.h --- Definitions of options for front-end programs
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


#ifndef PHAISTOS_OPTION_DEFINITIONS_H
#define PHAISTOS_OPTION_DEFINITIONS_H

// By default - module support is enabled
#ifndef MODULE_SUPPORT
#define MODULE_SUPPORT 1
#endif

#include "utils/utils.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <cstdlib>
#include <cstring>
#include <time.h>

#include "protein/chain_ca.h"
#include "protein/chain_fb.h"

#include "backbone_dbn.h"
#include "utils/random.h"
#include "energy/energy.h"

#include "protein/pdb_input.h"
#include "utils/program_option_parser.h"

#include "utils/indented_stream.h"

#include "energy/term_clash_fast.h"
#include "energy/term_rmsd.h"
#include "energy/term_rg.h"
#include "energy/term_contact_map.h"
#include "energy/term_backbone_dbn.h"
#include "energy/term_helix_content.h"
#include "energy/term_q_factor.h"
#include "energy/term_angle_histogram.h"
#include "energy/term_angle.h"
#include "energy/term_energy_sum.h"
#include "energy/term_energy_term_wrapper.h"

#include "energy/observable.h"
#include "energy/observable_collection.h"
#include "energy/observable_pdb.h"

#include "monte_carlo/monte_carlo_metropolis_hastings.h"
#include "monte_carlo/monte_carlo_simulated_annealing.h"
#include "monte_carlo/monte_carlo_greedy_optimization.h"

#include "protein/pdb_input.h"
#include "utils/program_option_parser.h"

#include "utils/indented_stream.h"

#include "moves/move_backbone_dbn.h"
#include "moves/move_fixed_structure.h"
#include "moves/rotamer_library_dunbrack.h"
#include "moves/move_crankshaft.h"
#include "moves/move_sidechain_rotamer.h"
#include "moves/move_sidechain_uniform.h"
#include "moves/move_sidechain_local.h"
#include "moves/move_pivot_local.h"
#include "moves/move_pivot_uniform.h"
#include "moves/move_none.h"

namespace phaistos {
// Forward declaration
double temperature_to_one_over_k(const double temperature);
}

// Module includes
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/includes.cpp"
#endif

namespace phaistos {

//! Vector containing the different modes of execution
std::vector<std::string> phaistos_modes;

//! Mode definitions
//! All definitions associated with different modes are 
//! gathered here, for ease of reference
class ModeDefinitions {
public:

     // These are just examples of varibles that can be defined by a mode
     // Feel free to add others

     //! Atom types included in chain
     std::string atom_types_default;

     //! Whether implicit energies are allowed in this mode
     bool implicit_energies_allowed;

     //! Initializer (general)
     template <typename CHAIN_TYPE>
     void initialize(ProgramOptionParser &target, CHAIN_TYPE *chain=NULL) {

          // Default values
          atom_types_default = "ALL_PHYSICAL_ATOMS";
          implicit_energies_allowed = true;
          target.super_group_default("backbone-dbn")   = "--backbone-dbn torus";
          
          // Insert module-specific mode-definitions
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/mode_definition.cpp"
#endif
     }

     //! Constructor (general)
     ModeDefinitions(ProgramOptionParser &target) {

          // Call initialize with dummy chain (and base class chain type)
          initialize(target, (Chain<Residue>*)NULL);
     }

     //! Constructor (ChainFB specific)
     ModeDefinitions(ProgramOptionParser &target, ChainFB *chain)
          : implicit_energies_allowed(true) {

          // Insert if not already present
          if (std::find(phaistos_modes.begin(), phaistos_modes.end(), "casp") == phaistos_modes.end())
          phaistos_modes.push_back("casp");

          initialize(target, chain);

          if (target.has_key("mode") && target["mode"].as<std::string>() == "casp") {
               atom_types_default  = "BACKBONE_O_ATOMS BACKBONE_H_ATOMS CB_ATOMS PSEUDO_SIDECHAIN_ATOMS";
               target.super_group_default("move")   = "--move dbn";
               target.super_group_default("energy") = "--energy clash-fast";
          } 
     }

     //! Constructor (ChainCA specific)
     ModeDefinitions(ProgramOptionParser &target, ChainCA *chain) {

          // default_values();

          atom_types_default  = "ALL_ATOMS";
          target.super_group_default("move")   = "--move dbn";

          initialize(target, chain);
     }
};

//! Global container for threaded copies of energy_function
template <typename CHAIN_TYPE>
class EnergyFunctionThreaded {
public:
     static std::vector<Energy<CHAIN_TYPE>*> energy_functions;
};

//! Functor defining no options
struct DefineOptionsNone {
    //! Functor evaluate function
     //!
     //! \param add_options Functor used to set options from fusion vector
     //! \param settings Target settings object     
     template <typename SETTINGS_TYPE>
     void operator()(ProgramOptionParser::AddOptions &add_options, SETTINGS_TYPE &settings) const {
     }
};


//! Functor defining options common to all BackboneDbn objects
struct DefineBackboneDbnCommonOptions {

     //! Functor evaluate function
     //!
     //! \param add_options Functor used to set options from fusion vector
     //! \param settings Target settings object     
     template <typename SETTINGS_TYPE>
     void operator()(ProgramOptionParser::AddOptions &add_options, SETTINGS_TYPE &settings) const {

          using namespace boost::fusion;

          //! Add options common to all moves
          for_each(
               make_vector(
                    make_vector(std::string("debug"),
                                std::string("Debug level"),
                                &settings->debug),
                    make_vector(std::string("log-space"),
                                std::string("Whether DBN calculations should be done in log space."),
                                &settings->log_space),
                    make_vector(std::string("dbn-start-distribution"),
                                std::string("N-terminus probability: normal|uniform|stationary"),
                                reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<TransitionEmissionState>*>(&settings->start_distribution)),
                    make_vector(std::string("dbn-transition-distribution"),
                                std::string("Transition distribution: normal|uniform|stationary"),
                                reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<TransitionEmissionState>*>(&settings->transition_distribution)),
                    make_vector(std::string("parameter-file"),
                                std::string("parameter file name"),
                                (std::string*)NULL,
                                std::string("")),
                    make_vector(std::string("sequence-length"),
                                std::string("Set sequence length in model (only necessary when no input sequence is given)"),
                                &settings->sequence_length),
                    make_vector(std::string("initial-pdb-file"),
                                std::string("Input from PDB file"),
                                &settings->initial_pdb_file),
                    make_vector(std::string("initial-aa-sequence"),
                                std::string("Input from amino acid sequence"),
                                &settings->initial_aa_sequence),
                    make_vector(std::string("initial-aa-file"),
                                std::string("Input from amino acid sequence file"),
                                &settings->initial_aa_file),
                    make_vector(std::string("initial-ss-sequence"),
                                std::string("Input from secondary structure sequence"),
                                &settings->initial_ss_sequence),
                    make_vector(std::string("initial-ss-file"),
                                std::string("Input from secondary structure file"),
                                &settings->initial_ss_file)
               ),
               add_options);
     }
};

//! Functor defining options common to all Move objects
struct DefineMoveCommonOptions {

     //! Functor evaluate function
     //!
     //! \param add_options Functor used to set options from fusion vector
     //! \param settings Target settings object     
     template <typename SETTINGS_TYPE>
     void operator()(ProgramOptionParser::AddOptions &add_options, SETTINGS_TYPE &settings) const {

          using namespace boost::fusion;

          // Add options common to all moves
          for_each(
               make_vector(
                    make_vector(std::string("debug"),
                                std::string("Debug level"),
                                &settings->debug),
                    make_vector(std::string("weight"),
                                std::string("Weight used when selecting moves in move collection"),
                                &settings->weight,
                                1.0),
                    make_vector(std::string("move-length-min"),
                                std::string("Minimum move length"),
                                &settings->move_length_min),
                    make_vector(std::string("move-length-max"),
                                std::string("Maximum move length"),
                                &settings->move_length_max),
                    // Regions is of a vector type - program_options handles vectors
                    // in a special way, expecting multiple tokens. We thus define
                    // our own vector type which inherits from std::vector and cast the 
                    // pointer to the settings member to the real vector type
                    make_vector(std::string("regions"),
                                std::string("Regions of chain in which move will be applied"),
                                static_cast<ProgramOptionParser::PlainVector<std::pair<int,int> >*>
                                (&settings->regions))
                    ),
               add_options);
     }
};


//! Functor which makes it possible to change the default settings objekt used
struct SettingsModifier  {

     // Conversion functionality in case SETTINGS_PARENT is an observable
     template <typename SETTINGS_PARENT>
     typename SETTINGS_PARENT::Settings *convert_settings(const boost::true_type &t, 
                                                          typename SETTINGS_PARENT::Settings *settings) const {
          return settings;
     }

     // Conversion functionality in case SETTINGS_PARENT is not an observable
     template <typename SETTINGS_PARENT>
     typename SETTINGS_PARENT::Settings *convert_settings(const boost::false_type &t, 
                                                          typename SETTINGS_PARENT::Settings *settings) const {
          typedef typename SETTINGS_PARENT::Settings Settings;
          Settings *settings_new = static_cast<Settings*>(new typename Observable<SETTINGS_PARENT>::Settings());
          delete settings;
          return settings_new;
     }

     template <typename SETTINGS_PARENT>
     typename SETTINGS_PARENT::Settings *modify(typename SETTINGS_PARENT::Settings *settings, std::string prefix="") const {

          // Attempt to cast settings object to an energy-setting object
          typename EnergyTerm<typename SETTINGS_PARENT::ChainType>::Settings *energy_term_settings = 
               dynamic_cast<typename EnergyTerm<typename SETTINGS_PARENT::ChainType>::Settings*>(settings);

          // If the settings is an energy settings object, check if prefix is "observable"
          // and modify Settings object if necessary
          if (energy_term_settings && prefix=="observable") {
               return convert_settings<SETTINGS_PARENT>(boost::is_base_of<ObservableBase, SETTINGS_PARENT>(), settings);
          } else {
               // In case class inherits from ObservableBase and prefix != "observable", return
               // a null pointer, which signals that option will not be allowed (and will not appear
               // as an option in help output)
               if (boost::is_base_of<ObservableBase,SETTINGS_PARENT>::value) {
                    delete settings;
                    return NULL;
               }
               return settings;
          }
     }
};


//! Functor defining options common to all Energy objects
struct DefineEnergyCommonOptions {

     //! Functor evaluate function
     //!
     //! \param add_options Functor used to set options from fusion vector
     //! \param settings Target settings object     
     template <typename SETTINGS_TYPE>
     void operator()(ProgramOptionParser::AddOptions &add_options, SETTINGS_TYPE &settings) const {

          using namespace boost::fusion;

          // Add options common to all moves
          for_each(
               make_vector(
                    make_vector(std::string("debug"),
                                std::string("Debug level"),
                                &settings->debug),
                    // Weight is not saved in a settings object
                    make_vector(std::string("weight"),
                                std::string("Weight used when summing energy terms"),
                                &settings->weight)
               ),
               add_options);

          ObservableBase::Settings *observable_settings = 
               dynamic_cast<ObservableBase::Settings*>(&(*settings));
          // typename Observable<ENERGY_TERM>::Settings *observable_settings = 
          //      dynamic_cast<typename Observable<ENERGY_TERM>::Settings*>(&(*settings));
          if (observable_settings) {
               // Add options specific to observables
               for_each(
                    make_vector(
                         make_vector(std::string("register-interval"),
                                     std::string("How often to register/calculate observable."),
                                     &observable_settings->register_interval),
                         make_vector(std::string("register-burnin"),
                                     std::string("After how many iterations to start registering/calculating observable."),
                                     &observable_settings->register_burnin),
                         make_vector(std::string("output-target"),
                                     std::string("How/Where the observable should be reported (\"pdb-header\":output information to header in dumped pdb files; \"pdb-b-factor\":output information to b-factors in dumped pdb files; \"stdout|cout|stderr|cerr\": Output to stdout|stderr. Any other string is interpreted as a filename for a separate logfile."),
                                     &observable_settings->output_target),
                         make_vector(std::string("output-interval"),
                                     std::string("How often to register/calculate observable. This value will only be read whenever an observable is active (i.e. iteration number matches register-interval)"),
                                     &observable_settings->output_interval)
                         ),
                    add_options);
          }

     }
};


//! Functor defining options common to all Monte Carlo objects
struct DefineMonteCarloCommonOptions {

     //! Functor evaluate function
     //!
     //! \param add_options Functor used to set options from fusion vector
     //! \param settings Target settings object     
     template <typename SETTINGS_TYPE>
     void operator()(ProgramOptionParser::AddOptions &add_options, SETTINGS_TYPE &settings) const {

          using namespace boost::fusion;

          // Add options common to all moves
          for_each(
               make_vector(
                    make_vector(std::string("debug"),
                                std::string("Debug level"),
                                &settings->debug),
                    make_vector(std::string("declash-on-reinitialize"),
                                std::string("Whether to remove self-collisions from the chain when reinitializing"),
                                &settings->declash_on_reinitialize),
                    make_vector(std::string("maximum-declash-attempts"),
                                std::string("The number of times declashing is attempted before a complete reinitialization is done"),
                                &settings->maximum_declash_attempts),
                    make_vector(std::string("reinitialization-interval"),
                                std::string("How often reinitialization takes place"),
                                &settings->reinitialization_interval),
                    make_vector(std::string("consistency-check-interval"),
                                std::string("How often consistency of the chain is checked"),
                                &settings->consistency_check_interval)
               ),
               add_options);
     }
};

//! Functor defining options common to all Procedure objects
struct DefineProcedureCommonOptions {

     //! Functor evaluate function
     //!
     //! \param add_options Functor used to set options from fusion vector
     void operator()(ProgramOptionParser::AddOptions &add_options) const {

          using namespace boost::fusion;

          // Add options common to all moves
          for_each(
               make_vector(
                    make_vector(std::string("debug"),
                                std::string("Debug level"),
                                (int*)NULL,
                                0)
               ),
               add_options);
     }
};

//! Module global option definitions are inserted here
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/global_option_definitions.cpp"
#endif

//! Map temperature to 1/k
double temperature_to_one_over_k(const double temperature) {
     // 3.2976*10^-27 kcal/K * 6.022*10^23 mol^-1 * temperature
     double one_over_kt = 1.0/(3.2976E-27*6.022E23 * temperature);
     return one_over_kt;
}


// BackboneDbn initialization
struct BackboneDBNOptions {

     // Constructor
     BackboneDBNOptions(ProgramOptionParser &target,
                        const ProgramOptionParser::Filter &occurrences,
                        std::string super_group) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // TorusDBN options
          if (occurrences["backbone-dbn-torus"]) {

               // Create settings object
              typedef TorusDBN::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineBackboneDbnCommonOptions(),
                        "TorusDbn options",
                        "backbone-dbn-torus", settings,
                        make_vector(
                             make_vector(std::string("initial-cis-sequence"),
                                         std::string("Cis sequence"),
                                         &settings->initial_cis_sequence),
                             make_vector(std::string("initial-cis-file"),
                                         std::string("Cis filename"),
                                         &settings->initial_cis_file)
                             )), super_group);
          }

          // TorusDBN options
          if (occurrences["backbone-dbn-torus-omega"]) {

               // Create settings object
              typedef TorusOmegaDBN::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineBackboneDbnCommonOptions(),
                        "TorusOmegaDbn options",
                        "backbone-dbn-torus-omega", settings,
                        make_vector(
                             make_vector(std::string("initial-cis-sequence"),
                                         std::string("Cis sequence"),
                                         &settings->initial_cis_sequence),
                             make_vector(std::string("initial-cis-file"),
                                         std::string("Cis filename"),
                                         &settings->initial_cis_file),
                             make_vector(std::string("initial-omega-sequence"),
                                         std::string("Omega sequence"),
                                         &settings->initial_omega_sequence),
                             make_vector(std::string("initial-omega-file"),
                                         std::string("Omega filename"),
                                         &settings->initial_omega_file)
                             )), super_group);
          }


          // FB5DBN options
          if (occurrences["backbone-dbn-fb5"]) {

               // Create settings object
              typedef FB5DBN::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineBackboneDbnCommonOptions(),
                        "Fb5Dbn options",
                        "backbone-dbn-fb5", settings,
                        make_vector(
                             )), super_group);
          }


          // TorusDbnCs options
          if (occurrences["backbone-dbn-torus-cs"]) {

               // Create settings object
              typedef TorusCsDBN::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineBackboneDbnCommonOptions(),
                        "TorusDbnCs options",
                        "backbone-dbn-torus-cs", settings,
                        make_vector(
                             make_vector(std::string("initial-cis-sequence"),
                                         std::string("Cis sequence"),
                                         &settings->initial_cis_sequence),
                             make_vector(std::string("initial-cis-file"),
                                         std::string("Cis filename"),
                                         &settings->initial_cis_file),
                             make_vector(std::string("initial-cs-ca-sequence"),
                                         std::string("Chemical shift CA sequence"),
                                         &settings->initial_cs_ca_sequence),
                             make_vector(std::string("initial-cs-ca-file"),
                                         std::string("Chemical shift CA file"),
                                         &settings->initial_cs_ca_file),
                             make_vector(std::string("initial-cs-cb-sequence"),
                                         std::string("Chemical shift CB sequence"),
                                         &settings->initial_cs_cb_sequence),
                             make_vector(std::string("initial-cs-cb-file"),
                                         std::string("Chemical shift CB file"),
                                         &settings->initial_cs_cb_file),
                             make_vector(std::string("initial-cs-c-sequence"),
                                         std::string("Chemical shift C sequence"),
                                         &settings->initial_cs_c_sequence),
                             make_vector(std::string("initial-cs-c-file"),
                                         std::string("Chemical shift C file"),
                                         &settings->initial_cs_c_file),
                             make_vector(std::string("initial-cs-n-sequence"),
                                         std::string("Chemical shift N sequence"),
                                         &settings->initial_cs_n_sequence),
                             make_vector(std::string("initial-cs-n-file"),
                                         std::string("Chemical shift N file"),
                                         &settings->initial_cs_n_file),
                             make_vector(std::string("initial-cs-ha-sequence"),
                                         std::string("Chemical shift HA sequence"),
                                         &settings->initial_cs_ha_sequence),
                             make_vector(std::string("initial-cs-ha-file"),
                                         std::string("Chemical shift HA file"),
                                         &settings->initial_cs_ha_file),
                             make_vector(std::string("initial-cs-h-sequence"),
                                         std::string("Chemical shift H sequence"),
                                         &settings->initial_cs_h_sequence),
                             make_vector(std::string("initial-cs-h-file"),
                                         std::string("Chemical shift H file"),
                                         &settings->initial_cs_h_file),
                             make_vector(std::string("initial-cs-sequence"),
                                         std::string("Chemical shift sequence - all in one. Order: CA,CB,C,N,HA,H"),
                                         &settings->initial_cs_sequence),
                             make_vector(std::string("initial-cs-file"),
                                         std::string("Chemical shift file - all in one. Order: CA,CB,C,N,HA,H "),
                                         &settings->initial_cs_file),
                             make_vector(std::string("initial-cs-nmr-star-file"),
                                         std::string("Chemical shift file - all in one, using CS-section of nmr-star format (BMRB)"),
                                         &settings->initial_cs_nmr_star_file),
                             make_vector(std::string("initial-cs-ca-offset"),
                                         std::string("CA chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_ca_offset),
                             make_vector(std::string("initial-cs-cb-offset"),
                                         std::string("CB chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_cb_offset),
                             make_vector(std::string("initial-cs-c-offset"),
                                         std::string("C chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_c_offset),
                             make_vector(std::string("initial-cs-n-offset"),
                                         std::string("N chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_n_offset),
                             make_vector(std::string("initial-cs-ha-offset"),
                                         std::string("HA chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_ha_offset),
                             make_vector(std::string("initial-cs-h-offset"),
                                         std::string("H chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_h_offset)
                             )), super_group);
          }

          // TorusDbnCs options
          if (occurrences["backbone-dbn-torus-cs-omega"]) {

               // Create settings object
              typedef TorusCsOmegaDBN::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineBackboneDbnCommonOptions(),
                        "TorusCsOmegaDbn options",
                        "backbone-dbn-torus-cs-omega", settings,
                        make_vector(
                             make_vector(std::string("initial-cis-sequence"),
                                         std::string("Cis sequence"),
                                         &settings->initial_cis_sequence),
                             make_vector(std::string("initial-cis-file"),
                                         std::string("Cis filename"),
                                         &settings->initial_cis_file),
                             make_vector(std::string("initial-omega-sequence"),
                                         std::string("Omega sequence"),
                                         &settings->initial_omega_sequence),
                             make_vector(std::string("initial-omega-file"),
                                         std::string("Omega filename"),
                                         &settings->initial_omega_file),
                             make_vector(std::string("initial-cs-ca-sequence"),
                                         std::string("Chemical shift CA sequence"),
                                         &settings->initial_cs_ca_sequence),
                             make_vector(std::string("initial-cs-ca-file"),
                                         std::string("Chemical shift CA file"),
                                         &settings->initial_cs_ca_file),
                             make_vector(std::string("initial-cs-cb-sequence"),
                                         std::string("Chemical shift CB sequence"),
                                         &settings->initial_cs_cb_sequence),
                             make_vector(std::string("initial-cs-cb-file"),
                                         std::string("Chemical shift CB file"),
                                         &settings->initial_cs_cb_file),
                             make_vector(std::string("initial-cs-c-sequence"),
                                         std::string("Chemical shift C sequence"),
                                         &settings->initial_cs_c_sequence),
                             make_vector(std::string("initial-cs-c-file"),
                                         std::string("Chemical shift C file"),
                                         &settings->initial_cs_c_file),
                             make_vector(std::string("initial-cs-n-sequence"),
                                         std::string("Chemical shift N sequence"),
                                         &settings->initial_cs_n_sequence),
                             make_vector(std::string("initial-cs-n-file"),
                                         std::string("Chemical shift N file"),
                                         &settings->initial_cs_n_file),
                             make_vector(std::string("initial-cs-ha-sequence"),
                                         std::string("Chemical shift HA sequence"),
                                         &settings->initial_cs_ha_sequence),
                             make_vector(std::string("initial-cs-ha-file"),
                                         std::string("Chemical shift HA file"),
                                         &settings->initial_cs_ha_file),
                             make_vector(std::string("initial-cs-h-sequence"),
                                         std::string("Chemical shift H sequence"),
                                         &settings->initial_cs_h_sequence),
                             make_vector(std::string("initial-cs-h-file"),
                                         std::string("Chemical shift H file"),
                                         &settings->initial_cs_h_file),
                             make_vector(std::string("initial-cs-sequence"),
                                         std::string("Chemical shift sequence - all in one. Order: CA,CB,C,N,HA,H"),
                                         &settings->initial_cs_sequence),
                             make_vector(std::string("initial-cs-file"),
                                         std::string("Chemical shift file - all in one. Order: CA,CB,C,N,HA,H "),
                                         &settings->initial_cs_file),
                             make_vector(std::string("initial-cs-nmr-star-file"),
                                         std::string("Chemical shift file - all in one, using CS-section of nmr-star format (BMRB)"),
                                         &settings->initial_cs_nmr_star_file),
                             make_vector(std::string("initial-cs-ca-offset"),
                                         std::string("CA chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_ca_offset),
                             make_vector(std::string("initial-cs-cb-offset"),
                                         std::string("CB chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_cb_offset),
                             make_vector(std::string("initial-cs-c-offset"),
                                         std::string("C chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_c_offset),
                             make_vector(std::string("initial-cs-n-offset"),
                                         std::string("N chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_n_offset),
                             make_vector(std::string("initial-cs-ha-offset"),
                                         std::string("HA chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_ha_offset),
                             make_vector(std::string("initial-cs-h-offset"),
                                         std::string("H chemical shift offset (for re-referencing)"),
                                         &settings->initial_cs_h_offset)
                             )), super_group);
          }

     }
};


// Procedure initialization
struct ProcedureOptions {

     // Constructor
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     ProcedureOptions(ProgramOptionParser &target,
                      const ProgramOptionParser::Filter &occurrences,
                      std::string super_group,
                      CHAIN_TYPE *chain,
                      DBN_TYPE *dbn) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Specific options for fold-procedure
          if (occurrences["procedure-fold"]) {
               // Add options
               target.add(
                    target.create_options(
                         DefineProcedureCommonOptions(),
                         "Options for Phaistos fold procedure",
                         "procedure-fold",
                         make_vector(
                              make_vector(std::string("energy2-evaluation-interval"),
                                          std::string("How often to evaluate energy2 (0:never)"),
                                          (int*)NULL,
                                          0)
                              )
                         ), super_group);
          }

          // Specific options for comparison-procedure
          if (occurrences["procedure-comparison"]) {
               // Add options
               target.add(
                    target.create_options(
                         DefineProcedureCommonOptions(),
                         "Options for Phaistos comparison procedure",
                         "procedure-comparison",
                         make_vector(
                              )
                         ), super_group);
          }
     }
};

// Monte Carlo initialization
struct MonteCarloOptions {

     // Constructor
     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     MonteCarloOptions(ProgramOptionParser &target,
                       const ProgramOptionParser::Filter &occurrences,
                       std::string super_group,
                       std::vector<std::string> &mc_mode_description,
                       CHAIN_TYPE *chain,
                       DBN_TYPE *dbn) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Monte Carlo options
          ProgramOptionParser::OptionsDescription options_monte_carlo("General Monte Carlo options");
          options_monte_carlo.add_options()
               ("threads", po::value<int>()->default_value(1),
                "Number of threads (1: no multithreading)")
               ("identical-threads", po::value<bool>()->default_value(0)->implicit_value(1),
                "Make all threads identical")
               ("iterations", po::value<PHAISTOS_LONG_LONG>()->default_value(10000000),
                "Number of iterations pr. thread")
               ("steps-per-move", po::value<int>()->default_value(100),
                "Number of steps per each move in Monte Carlo run")
               ;
          bool hidden = false;
          target.add(options_monte_carlo, hidden, super_group);


          // Metropolis Hastings options
          if (occurrences["monte-carlo-metropolis-hastings"]) {

               // Create settings object
              typedef typename MonteCarloMetropolisHastings<CHAIN_TYPE>::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineMonteCarloCommonOptions(),
                        "Metropolis-Hastings options",
                        "monte-carlo-metropolis-hastings", settings,
                        make_vector()), super_group);
          }
          mc_mode_description.push_back("metropolis-hastings");


          // Simulated annealing options
          if (occurrences["monte-carlo-simulated-annealing"]) {

               // Create settings object
              typedef typename MonteCarloSimulatedAnnealing<CHAIN_TYPE>::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineMonteCarloCommonOptions(),
                        "Simulated-Annealing options",
                        "monte-carlo-simulated-annealing", settings,
                        make_vector(
                             make_vector(std::string("emin"),
                                         std::string("Lower bound on energy"),
                                         reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->energy_min)),
                             make_vector(std::string("emax"),
                                         std::string("Upper bound on energy"),
                                         reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->energy_max)),
                             make_vector(std::string("temperature-start"),
                                         std::string("Starting temperature (must be set to a finite value)"),
                                         &settings->temperature_start),
                             make_vector(std::string("temperature-end"),
                                         std::string("Ending temperature (must be set to a finite value)"),
                                         &settings->temperature_end)
                             )), super_group);
          }
          mc_mode_description.push_back("simulated-annealing");

          // Greedy optimization options
          if (occurrences["monte-carlo-greedy-optimization"]) {

               // Create settings object
              typedef typename MonteCarloGreedyOptimization<CHAIN_TYPE>::Settings Settings;
              boost::shared_ptr<Settings> settings(new Settings());

              // Add options
              target.add(
                   target.create_options(
                        DefineMonteCarloCommonOptions(),
                        "Greedy optimization options",
                        "monte-carlo-greedy-optimization", settings,
                        make_vector()), super_group);
          }
          mc_mode_description.push_back("greedy-optimization");
     }
};


// Move initialization
template <typename SETTINGS_MODIFIER>
struct MoveOptions {

     //! Initializer
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     void common_options(ProgramOptionParser &target,
                         const ProgramOptionParser::Filter &occurrences,
                         std::string super_group,
                         CHAIN_TYPE *chain,
                         DBN_TYPE *dbn) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // None move: no structural resampling
          for (int counter = occurrences["move-none"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveNone<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "None - no structural resampling",
                         "move-none", settings,
                         make_vector(
                              )), super_group, counter==1);
          }

//! Module move option definitions are inserted here
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/move_options.cpp"
#endif
     }

     // Constructor - general case
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 CHAIN_TYPE *chain,
                 DBN_TYPE *dbn) {
          common_options(target, occurrences, super_group, chain, dbn);
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 ChainFB *chain,
                 DBN_TYPE *dbn) {

          common_options(target, occurrences, super_group, chain, dbn);

          // Facilitate use of make_vector
          using namespace boost::fusion;
          
          // Define defaults for the different modes
          ModeDefinitions mode_definitions(target, chain);

          // MoveBackboneDBN
          for (int counter = occurrences["move-backbone-dbn"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveBackboneDBN<ChainFB,DBN_TYPE> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "DBN",
                         "move-backbone-dbn", settings,
                         make_vector(
                              make_vector(std::string("resample-mode"),
                                          std::string("Specifies what is resampled when the move is applied: resample-all|resample-hidden-only|resample-angles-only"),
                                          reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<MoveBackboneDBNResampleMode>*>(&settings->resample_mode)),
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the dihedral energy from the dbn should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy),
                              make_vector(std::string("dbn-consistency-window-size"),
                                          std::string("Size of window used (to each side) when bringing the dbn back to consistency. A good value for the window size is >7, and a negative window size means that the full hidden node sequence is resampled."),
                                          &settings->dbn_consistency_window_size),
                              make_vector(std::string("dbn-bias-window-size"),
                                          std::string("Size of window used when calculating bias. Approximates the move bias as (X[i-w,j+w])/P(X'[i-w,j+w]), where w in the window size and [i,j] is the interval where angles have been changed. A good value for the window size is >7, and a negative window size means that the full bias is used."),
                                          &settings->dbn_bias_window_size)
                         )), super_group, counter==1);
          }



          // Crankshaft
          for (int counter = occurrences["move-crankshaft"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveCrankshaft<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Crankshaft",
                         "move-crankshaft", settings,
                         make_vector(
                              make_vector(std::string("only-internal-moves"),
                                          std::string("only execute internal moves"),
                                          &settings->only_internal_moves),
                              make_vector(std::string("constrain-bond-angle"),
                                          std::string("constrain bond angle variation"),
                                          &settings->constrain_bond_angle),
                              make_vector(std::string("bond-angle-tolerance"),
                                          std::string("tolerance of  bond angle variation (deg)"),
                                          &settings->bond_angle_tolerance)
                              )), super_group, counter==1);
          }

          // Local pivot move - no priors
          for (int counter = occurrences["move-pivot-local"]; counter > 0; counter--) {

               // Create settings object
               typedef MovePivotLocal<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Local pivot move",
                         "move-pivot-local", settings,
                         make_vector(
                              make_vector(std::string("sample-bond-angle-dofs"),
                                          std::string("Whether bond angles should be sampled"),
                                          &settings->sample_bond_angle_dofs),
                              make_vector(std::string("sample-phi-psi-dofs"),
                                          std::string("Whether dihedral angles should be sampled"),
                                          &settings->sample_phi_psi_dofs),
                              make_vector(std::string("sample-omega-dofs"),
                                          std::string("Whether omega angles should be sampled"),
                                          &settings->sample_omega_dofs),
                              make_vector(std::string("std-dev-bond-angle"),
                                          std::string("Standard deviation in bond angle change (degr.) (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_bond_angle)),
                              make_vector(std::string("std-dev-phi-psi"),
                                          std::string("Standard deviation in dihedral angle change (degr.) (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_phi_psi)),
                              make_vector(std::string("std-dev-omega"),
                                          std::string("Standard deviation in omega angle change (degr.) (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_omega))
                              )), super_group, counter==1);
          }


          // Local pivot move - dihedral prior: DBN
          for (int counter = occurrences["move-pivot-local-dbn"]; counter > 0; counter--) {

               // Create settings object
               typedef MovePivotLocal<ChainFB,
                                      MovePriorDbn<BondAnglePriorEnghHuber,
                                                   DBN_TYPE> > Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Override default settings
               settings->sample_bond_angle_dofs = false;

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Local pivot move (DBN prior)",
                         "move-pivot-local-dbn", settings,
                         make_vector(
                              make_vector(std::string("sample-bond-angle-dofs"),
                                          std::string("Whether bond angles should be sampled"),
                                          &settings->sample_bond_angle_dofs),
                              make_vector(std::string("sample-phi-psi-dofs"),
                                          std::string("Whether dihedral angles should be sampled"),
                                          &settings->sample_phi_psi_dofs),
                              make_vector(std::string("sample-omega-dofs"),
                                          std::string("Whether omega angles should be sampled"),
                                          &settings->sample_omega_dofs),
                              make_vector(std::string("std-dev-bond-angle"),
                                          std::string("Standard deviation in bond angle change (degr.) (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_bond_angle)),
                              make_vector(std::string("std-dev-phi-psi"),
                                          std::string("Standard deviation in dihedral angle change (degr.) (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_phi_psi)),
                              make_vector(std::string("std-dev-omega"),
                                          std::string("Standard deviation in omega angle change (degr.) (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_omega))
                              )), super_group, counter==1);
          }


          // Pivot move using uniform phi,psi samples
          for (int counter = occurrences["move-pivot-uniform"]; counter > 0; counter--) {

               // Create settings object
               typedef MovePivotUniform<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Pivot - sampling changes to phi,psi angles from uniform distributions",
                         "move-pivot-uniform", settings,
                         make_vector(
                              make_vector(std::string("single-dof-only"),
                                          std::string("Whether only to resample a single dof in each iteration (selected randomly)"),
                                          &settings->single_dof_only),
                              make_vector(std::string("skip-proline-phi"),
                                          std::string("Whether to skip prolines phi angles (modifiying the proline phi angle introduces an improper torsion change)"),
                                          &settings->skip_proline_phi),
                              make_vector(std::string("max-delta"),
                                          std::string("Maximum change in angle value"),
                                          &settings->max_delta)
                              )), super_group, counter==1);
          }

          // Sidechain Rotamer
          for (int counter = occurrences["move-sidechain-rotamer"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveSidechainRotamer<ChainFB, RotamerLibraryDunbrack<> > Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Override default settings
               settings->implicit_energy = mode_definitions.implicit_energies_allowed;

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Rotamer sidechain move",
                         "move-sidechain-rotamer", settings,
                         make_vector(
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the rotamer bias (implicit energy) should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy),
                              make_vector(std::string("sigma-scale-factor"),
                                          std::string("Flatten/sharpen rotamer distributions by scaling sigma"),
                                          &settings->sigma_scale_factor),
                              make_vector(std::string("rotamer-state-resample-frequency"),
                                          std::string("Frequency with which move will resample rotamer state"),
                                          &settings->rotamer_state_resample_frequency),
                              make_vector(std::string("sample-hydrogen-chis"),
                                          std::string("Whether hydrogen chi angles should be resampled (uniformly)"),
                                          &settings->sample_hydrogen_chis),
                              make_vector(std::string("skip-proline"),
                                          std::string("Whether to skip prolines (prolines introduce a change in bond length which must be taken into account by the forcefield)"),
                                          &settings->skip_proline)
                              )), super_group, counter==1);
          }

          // Sidechain Uniform
          for (int counter = occurrences["move-sidechain-uniform"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveSidechainUniform<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Uniform sidechain move",
                         "move-sidechain-uniform", settings,
                         make_vector(
                              make_vector(std::string("single-dof-only"),
                                          std::string("Whether only to resample a single dof in each iteration (selected randomly)"),
                                          &settings->single_dof_only),
                              make_vector(std::string("skip-proline"),
                                          std::string("Whether to skip prolines (prolines introduce a change in bond length which must be taken into account by the forcefield)"),
                                          &settings->skip_proline)
                              )), super_group, counter==1);
          }


          // Sidechain local move
          for (int counter = occurrences["move-sidechain-local"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveSidechainLocal<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Local sidechain move",
                         "move-sidechain-local", settings,
                         make_vector(
                              make_vector(std::string("sigma-major-dofs"),
                                          std::string("Standard deviation of gaussian surrounding current sidechain for major degrees of freedom"),
                                          &settings->sigma_major_dofs),
                              make_vector(std::string("sigma-minor-dofs"),
                                          std::string("Standard deviation of gaussian surrounding current sidechain for minor degrees of freedom"),
                                          &settings->sigma_minor_dofs),
                              make_vector(std::string("sample-major-dofs"),
                                          std::string("Whether major dofs (chi angles) should be resampled"),
                                          &settings->sample_major_dofs),
                              make_vector(std::string("sample-minor-dofs"),
                                          std::string("Whether minor dofs (side chain bond angles, CB|O bond angle|dihedral, etc)"),
                                          &settings->sample_minor_dofs),
                              make_vector(std::string("mode"),
                                          std::string("Move mode (simple,scale-sigma,select-dofs,constrain-all-endpoints,constrain-one-endpoint)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<MoveSidechainLocalMode>*>(&settings->mode)),
                              make_vector(std::string("lagrange-multiplier"),
                                          std::string("When using constrain-endpoint mode, this controls the weight of the constraint."),
                                          &settings->lagrange_multiplier),
                              make_vector(std::string("skip-proline"),
                                          std::string("Whether to skip prolines (prolines introduce a change in bond length which must be taken into account by the forcefield)"),
                                          &settings->skip_proline)
                              )), super_group, counter==1);
          }

     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 ChainCA *chain,
                 DBN_TYPE *dbn) {

          common_options(target, occurrences, super_group, chain, dbn);

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // MoveBackboneDBN
          for (int counter = occurrences["move-backbone-dbn"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveBackboneDBN<ChainCA,FB5DBN> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "DBN",
                         "move-backbone-dbn", settings,
                         make_vector(
                              make_vector(std::string("resample-mode"),
                                          std::string("Specifies what is resampled when the move is applied: resample-all|resample-hidden-only|resample-angles-only"),
                                          reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<MoveBackboneDBNResampleMode>*>(&settings->resample_mode)),
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the dihedral energy from the dbn should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy),
                              make_vector(std::string("dbn-consistency-window-size"),
                                          std::string("Size of window used (to each side) when bringing the dbn back to consistency. A good value for the window size is >7, and a negative window size means that the full hidden node sequence is resampled."),
                                          &settings->dbn_consistency_window_size),
                              make_vector(std::string("dbn-bias-window-size"),
                                          std::string("Size of window used when calculating bias. Approximates the move bias as (X[i-w,j+w])/P(X'[i-w,j+w]), where w in the window size and [i,j] is the interval where angles have been changed. A good value for the window size is >7, and a negative window size means that the full bias is used."),
                                          &settings->dbn_bias_window_size)
                         )), super_group, counter==1);
          }

     }
};


// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     //! Initializer
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     void common_options(ProgramOptionParser &target,
                         const ProgramOptionParser::Filter &occurrences,
                         std::string super_group,
                         std::string prefix,
                         CHAIN_TYPE *chain,
                         DBN_TYPE *dbn) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Empty energy term (does nothing) - useful for detailed balance simulations
          if (occurrences[prefix+"-none"]) {
               target.add(
                    target.create_options(
                         "none - empty energy term (" + prefix + ")",
                         prefix+"-none",
                         make_vector()),
                    super_group);
          }

          
          if (prefix == "observable") {

               // Sum of all energy terms - for inclusion in observables
               for (int counter = occurrences[prefix+"-@energy-sum"]; counter > 0; counter--) {

                    // Create settings object
                    typedef TermEnergySum<CHAIN_TYPE> EnergyTerm;
                    typedef typename EnergyTerm::Settings Settings;
                    boost::shared_ptr<Settings> settings(
                         SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

                    // Add options
                    target.add(
                         target.create_options(
                              DefineEnergyCommonOptions(),
                              "Sum of all energy terms (for inclusion as observable)",
                              prefix+"-@energy-sum", settings,
                              make_vector(
                                   )), super_group, counter==1);
               }

               // include all energy terms in observable collection
               for (int counter = occurrences[prefix+"-@energy-terms"]; counter > 0; counter--) {

                    // Create settings object
                    typedef TermEnergyTermWrapper<CHAIN_TYPE> EnergyTerm;
                    typedef typename EnergyTerm::Settings Settings;
                    boost::shared_ptr<Settings> settings(
                         SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

                    // Add options
                    target.add(
                         target.create_options(
                              DefineEnergyCommonOptions(),
                              "Refer to energy terms in main energy",
                              prefix+"-@energy-terms", settings,
                              make_vector(
                                   )), super_group, counter==1);
               }
          }


          // Clash-fast
          for (int counter = occurrences[prefix+"-clash-fast"]; counter > 0; counter--) {

               // Create settings object
               typedef TermClashFast<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Clash-fast (" + prefix + ")",
                         prefix+"-clash-fast", settings,
                         make_vector(
                              make_vector(std::string("only-modified-pairs"),
                                          std::string("Specifies whether energy should consider only pairs modified by the last move. If a clashfree structure is maintained at all times, this is sufficient to detect all clashes"),
                                          &settings->only_modified_pairs),
                              make_vector(std::string("boolean-mode"),
                                          std::string("Specifies whether energy should work in clash/non-clash mode and return infinity/0 (true) or count all the clashes and return the number of clashes (false)"),
                                          &settings->boolean_mode),
                              make_vector(std::string("minimum-residue-distance"),
                                          std::string("Minimum distance along chain (measured in number of residues) before pair is taken into account by the energy"),
                                          &settings->minimum_residue_distance),
                              make_vector(std::string("clash-distance-h"),
                                          std::string("Distance within which pairs of hydrogens are considered to be clashing (in angstrom)"),
                                          &settings->clash_distance_H),
                              make_vector(std::string("clash-distance-no"),
                                          std::string("Distance within which pairs of Nitrogen-Oxygen pairs are considered to be clashing (in angstrom)"),
                                          &settings->clash_distance_NO),
                              make_vector(std::string("clash-distance-ps"),
                                          std::string("Distance within which pairs of pseudo-sidechain atom pairs are considered to be clashing (in angstrom)"),
                                          &settings->clash_distance_PS),
                              make_vector(std::string("clash-distance-sg"),
                                          std::string("Distance within which pairs of SG atom pairs are considered to be clashing (in angstrom)"),
                                          &settings->clash_distance_SG),
                              make_vector(std::string("clash-distance-any-pair"),
                                          std::string("Default distance within which pairs of atom pairs of all other types are considered to be clashing (in angstrom)"),
                                          &settings->clash_distance_any_pair)
                              )), super_group, counter==1);

          }


          // RG
          for (int counter = occurrences[prefix+"-rg"]; counter > 0; counter--) {

               // Create settings object
               typedef TermRg<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "rg: Radius of gyration (" + prefix + ")",
                         prefix+"-rg", settings,
                         make_vector(
                              )),
                    super_group, counter==1);
          }


          // RMSD
          for (int counter = occurrences[prefix+"-rmsd"]; counter > 0; counter--) {

               // Create settings object
               typedef TermRmsd<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "rmsd: Root mean square deviation (" + prefix + ")",
                         prefix+"-rmsd", settings,
                         make_vector(
                              make_vector(std::string("reference-pdb-file"),
                                          std::string("PDB containing reference structure."),
                                          &settings->reference_pdb_file),
                              make_vector(std::string("ca-only"),
                                          std::string("Whether to only include C-alphas in the RMSD calculation."),
                                          &settings->ca_only),
                              make_vector(std::string("residue-start"),
                                          std::string("Start residue."),
                                          &settings->residue_start),
                              make_vector(std::string("residue-end"),
                                          std::string("End residue."),
                                          &settings->residue_end)
                              )), super_group, counter==1);
          }


          // Contact map
          for (int counter = occurrences[prefix+"-contact-map"]; counter > 0; counter--) {

               // Create settings object
               typedef TermContactMap<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "contact-map (" + prefix + ")",
                         prefix+"-contact-map", settings,
                         make_vector(
                              make_vector(std::string("contact-map-file"),
                                          std::string("path to contact map file"),
                                          &settings->contact_map_filename),
                              make_vector(std::string("contact-map-string"),
                                          std::string("String containing contacts"),
                                          &settings->contact_map_string),
                              make_vector(std::string("pdb-file"),
                                          std::string("path to pdb file for constructing contact map"),
                                          &settings->pdb_filename),
                              make_vector(std::string("cutoff"),
                                          std::string("contact cutoff distance when constructing map from pdb"),
                                          &settings->cutoff),
                              make_vector(std::string("minimum-residue-distance"),
                                          std::string("Minimum separation between two residues before they are considered."),
                                          &settings->minimum_residue_distance),
                              make_vector(std::string("iteration-type"),
                                          std::string("Which atoms to include when constructing contactmap from chain"),
                                          &settings->iteration_type),
                              make_vector(std::string("dist-squared"),
                                          std::string("set contact energy ~ (r-r0)^2/width (default: exp(-(r-r0)^2/width)"),
                                          &settings->dist_squared),
                              make_vector(std::string("potential-width"),
                                          std::string("Width of energy potential"),
                                          &settings->potential_width),
                              ((prefix=="observable") ? 
                               make_vector(std::string("verbose"),
                                           std::string("Whether to output contact map in verbose mode"),
                                           &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->verbose) :
                               make_vector(std::string(""), std::string(""), (bool*)NULL)),
                              ((prefix=="observable") ? 
                               make_vector(std::string("average-mode"),
                                           std::string("Whether to report contact map averages"),
                                           &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->average_mode) :
                               make_vector(std::string(""), std::string(""), (bool*)NULL)),
                              ((prefix=="observable") ? 
                                    make_vector(std::string("persistency-cutoff"),
                                                std::string("In average mode: only contacts with a frequency greater than this value are reported "),
                                                &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->persistency_cutoff) :
                               make_vector(std::string(""), std::string(""), (double*)NULL)),
                              ((prefix=="observable") ? 
                                    make_vector(std::string("counts-as-weights"),
                                                std::string("In average mode: whether to use the occurence-counts as weights"),
                                                &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->counts_as_weights) :
                               make_vector(std::string(""), std::string(""), (bool*)NULL))
                              )),
                    super_group, counter==1);

               // if (prefix=="observable") {
               //      target.add(
               //           target.create_options(
               //                DefineOptionsNone(),
               //                "contact-map (" + prefix + ")",
               //                prefix+"-contact-map", settings,
               //                 make_vector(
               //                      make_vector(std::string("verbose"),
               //                                  std::string("Whether to output contact map in verbose mode"),
               //                                  &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->verbose),
               //                      make_vector(std::string("average-mode"),
               //                                  std::string("Whether to report contact map averages"),
               //                                  &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->average_mode),
               //                      make_vector(std::string("persistency-cutoff"),
               //                                  std::string("In average mode: only contacts with a frequency greater than this value are reported "),
               //                                  &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->persistency_cutoff),
               //                      make_vector(std::string("counts-as-weights"),
               //                                  std::string("In average mode: whether to use the occurence-counts as weights"),
               //                                  &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->counts_as_weights)
               //                      )),
               //           super_group, counter==1);
               // }

          }          

          // Q-factor (nativeness)
          for (int counter = occurrences[prefix+"-q-factor"]; counter > 0; counter--) {

               // Create settings object
               typedef TermQFactor<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "q-factor: percentage of native contacts in chain (" + prefix + ")",
                         prefix+"-q-factor", settings,
                         make_vector(
                              make_vector(std::string("contact-map-file"),
                                          std::string("path to contact map file"),
                                          &settings->contact_map_filename),
                              make_vector(std::string("contact-map-string"),
                                          std::string("String containing contacts"),
                                          &settings->contact_map_string),
                              make_vector(std::string("pdb-file"),
                                          std::string("path to pdb file for constructing contact map"),
                                          &settings->pdb_filename),
                              make_vector(std::string("cutoff"),
                                          std::string("contact cutoff distance when constructing map from pdb"),
                                          &settings->cutoff),
                              make_vector(std::string("minimum-residue-distance"),
                                          std::string("Minimum separation between two residues before they are considered."),
                                          &settings->minimum_residue_distance),
                              make_vector(std::string("iteration-type"),
                                          std::string("Which atoms to include when constructing contactmap from chain"),
                                          &settings->iteration_type),
                              make_vector(std::string("use-weights"),
                                          std::string("Whether to use contact weights to calculate a weighted version of the q factor"),
                                          &settings->use_weights)
                              )),
                    super_group, counter==1);
          }          

          // helix content
          for (int counter = occurrences[prefix+"-helix-content"]; counter > 0; counter--) {

               // Create settings object
               typedef TermHelixContent<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "helix-content: percentage of chain that is in a helical conformation (" + prefix + ")",
                         prefix+"-helix-content", settings,
                         make_vector(
                              make_vector(std::string("min-angle1"),
                                          std::string("Minimum boundary for angle1"),
                                          &settings->min_angle1),
                              make_vector(std::string("max-angle1"),
                                          std::string("Maximum boundary for angle1"),
                                          &settings->max_angle1),
                              make_vector(std::string("min-angle2"),
                                          std::string("Minimum boundary for angle2"),
                                          &settings->min_angle2),
                              make_vector(std::string("max-angle2"),
                                          std::string("Maximum boundary for angle2"),
                                          &settings->max_angle2),
                              ((prefix=="observable") ? 
                               make_vector(std::string("per-residue"),
                                           std::string("Whether to split up observable into a per-residue vector"),
                                           &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->per_residue) :
                               make_vector(std::string(""), std::string(""), (bool*)NULL))
                              )),
                    super_group, counter==1);

               // if (prefix=="observable") {
               //      target.add(
               //           target.create_options(
               //                DefineOptionsNone(),
               //                "helix-content: percentage of chain that is in a helical conformation (" + prefix + ")",
               //                prefix+"-helix-content", settings,
               //                 make_vector(
               //                      make_vector(std::string("per-residue"),
               //                                  std::string("Whether to split up observable into a per-residue vector"),
               //                                  &dynamic_cast<typename Observable<EnergyTerm>::Settings*>(&*settings)->per_residue)
               //                      )),
               //           super_group, counter==1);
               // }
          }          

          // angle histogram
          for (int counter = occurrences[prefix+"-angle-histogram"]; counter > 0; counter--) {

               // Create settings object
               typedef TermAngleHistogram<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "angle-histogram: statistics on angle distributions (" + prefix + ")",
                         prefix+"-angle-histogram", settings,
                         make_vector(
                              make_vector(std::string("bins"),
                                          std::string("Number of bins"),
                                          &settings->bins)
                              )),
                    super_group, counter==1);
          }
          
          // angle
          for (int counter = occurrences[prefix+"-angle"]; counter > 0; counter--) {

               // Create settings object
               typedef TermAngle<CHAIN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "angle: statistics on angle distributions (" + prefix + ")",
                         prefix+"-angle", settings,
                         make_vector(
                              make_vector(std::string("dihedral-only"),
                                          std::string("Dump only the dihedral angels."),
                                          &settings->dihedral_only)
                              )),
                    super_group, counter==1);
          } 

          // PDB output
          for (int counter = occurrences[prefix+"-pdb"]; counter > 0; counter--) {

               // Create settings object (note that this is an observable)
               typedef Observable<TermPdb<CHAIN_TYPE> > EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Since we are adding a term that can only be applied as observable, make sure 
               // that we do not add it to anything but the observable collection
               if (prefix!="observable")
                    continue;

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "pdb: dump pdb file(s) (" + prefix + ")",
                         prefix+"-pdb", settings,
                         make_vector(
                              make_vector(std::string("minimum-energy-mode"),
                                          std::string("Whether only to dump minimum energy structures"),
                                          &settings->minimum_energy_mode),
                              make_vector(std::string("minimum-energy-fraction"),
                                          std::string("How far from the minimum energy a structure must be before being dumped"),
                                          &settings->minimum_energy_fraction),
                              make_vector(std::string("minimum-energy-interval"),
                                          std::string("Minimum interval between dumped minimum energy structures"),
                                          &settings->minimum_energy_interval)
                              )),
                    super_group, counter==1);
          } 


//! Module energy option definitions are inserted here
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/energy_options.cpp"         
#endif

     }

     // Constructor - general case
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {
          common_options(target, occurrences, super_group, prefix, chain, dbn);
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          common_options(target, occurrences, super_group, prefix, chain, dbn);

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Local DBN energy term
          for (int counter = occurrences[prefix+"-backbone-dbn"]; counter > 0; counter--) {

               // Create settings object
               typedef TermBackboneDBN<ChainFB,DBN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "backbone-dbn: Backbone DBN energy term (" + prefix + ")",
                         prefix+"-backbone-dbn", settings,
                         make_vector(
                              make_vector(std::string("enable-dbn-update"),
                                          std::string("Whether to update the angles in the DBN from the chain when necessary."),
                                          &settings->enable_dbn_update),
                              make_vector(std::string("always-full-update"),
                                          std::string("Whether to always force update of all angles in the DBN from the chain."),
                                          &settings->always_full_update),
                              make_vector(std::string("window-size"),
                                          std::string("Size of window used when calculating the energy. A good value for the window size is >7, and a negative window size means that the full bias is used."),
                                          &settings->window_size),
                              make_vector(std::string("eliminate-move-bias"),
                                          std::string("Divide out the move-bias of the corresponding moves. Equivalent to (but faster than) setting implicit-energy to false."),
                                          &settings->eliminate_move_bias)
                              )), super_group, counter==1);
          }
     }

     // Constructor - ChainCA specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainCA *chain,
                   DBN_TYPE *dbn) {

          common_options(target, occurrences, super_group, prefix, chain, dbn);

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Local DBN energy term
          for (int counter = occurrences[prefix+"-backbone-dbn"]; counter > 0; counter--) {

               // Create settings object
               typedef TermBackboneDBN<ChainCA,DBN_TYPE> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "backbone-dbn: Local DBN energy term (" + prefix + ")",
                         prefix+"-backbone-dbn", settings,
                         make_vector(
                              make_vector(std::string("enable-dbn-update"),
                                          std::string("Whether DBN state should be copied from chain in each iteration."),
                                          &settings->enable_dbn_update),
                              make_vector(std::string("always-full-update"),
                                          std::string("Whether entire DBN state should be copied from chain in each iteration."),
                                          &settings->always_full_update),
                              make_vector(std::string("window-size"),
                                          std::string("Size of window used when calculating the energy. A good value for the window size is >7, and a negative window size means that the full bias is used."),
                                          &settings->window_size),
                              make_vector(std::string("eliminate-move-bias"),
                                          std::string("Divide out the move-bias of the corresponding moves. Equivalent to (but faster than) setting implicit-energy to false."),
                                          &settings->eliminate_move_bias)
                              )), super_group, counter==1);
          }
     }
};



//! Option class. Uses boost::program_options.
//! Interface to command line options and optionally
//! a container to Settings objects
class Options {
public:

     //! Local typedef for OptionShorthand
     typedef ProgramOptionParser::OptionShorthand OptionShorthand;

     //! Local typedef for OptionValue
     typedef ProgramOptionParser::OptionValue OptionValue;

     //@{
     //! Internal boost::program_options parser objects
     ProgramOptionParser *parser_0;
     ProgramOptionParser *parser_1;
     ProgramOptionParser *parser_2;
     ProgramOptionParser *parser_3;
     ProgramOptionParser *parser_4;
     ProgramOptionParser *parser_5;
     //@}

     //! Flag specifying whether full parse has completed
     bool fully_parsed;

     //! Flag specifying whether an empty command line was provided
     bool empty_command_line;

     //! Command line - a convenient representation of
     //! a command line both in string and argc,argv format
     ProgramOptionParser::CommandLine command_line;


     //! Initial parse function
     //!
     //! \tparam CHAIN_TYPE Molecule chain type
     //! \tparam DBN_TYPE Dynamic Bayesian Network type
     //! \tparam DEFINE_OPTIONS_PRIMARY Functor containing primary option definitions
     //! \tparam DEFINE_OPTIONS_SECONDARY Functor containing secondary option definitions
     //! \param force_reparsing Whether to force another round of parsing
     template <typename DEFINE_OPTIONS_PRIMARY,typename DEFINE_OPTIONS_SECONDARY>
     void parse_initial(bool force_reparsing=false) {
          bool allow_unregistered = true;
          for (unsigned int i=0; i<2; ++i) {
               // Start by a simple pass over the primary options
               // The rest is delayed until parse is called
               // Set mode-specific options
               parser_0->parse(DEFINE_OPTIONS_PRIMARY(), DEFINE_OPTIONS_SECONDARY(), command_line, allow_unregistered);

               // std::cout << "0. Expanded command line: " << parser_0->command_line.raw_string << "\n";

               std::string config_file = "";
               if (parser_0->has_key("config-file")) {
                    config_file = (*parser_0)["config-file"].as<std::string>();
               }
               parser_1->parse(DEFINE_OPTIONS_PRIMARY(), DEFINE_OPTIONS_SECONDARY(), command_line, allow_unregistered, config_file);

               // std::cout << "1. Expanded command line: " << parser_1->command_line.raw_string << "\n";

               if (i==0) {
                    // bool primary_options_only = true;
                    ModeDefinitions mode_definitions(*parser_1);
                    bool reparse_necessary = parser_1->apply_super_group_defaults();
                    if (reparse_necessary || force_reparsing) {
                         delete parser_0;
                         delete parser_1;
                         this->parser_0 = new ProgramOptionParser(command_line);
                         this->parser_1 = new ProgramOptionParser(command_line);
                         if (force_reparsing)
                              allow_unregistered = false;
                         continue;
                    } else {
                         break;
                    }
               }
          }

          // Set the random seed
          random_global.seed((*this)["seed"].as<unsigned int>());

     }

     //! Main parse function
     //!
     //! \tparam CHAIN_TYPE Molecule chain type
     //! \tparam DBN_TYPE Dynamic Bayesian Network type
     //! \tparam DEFINE_OPTIONS_PRIMARY Functor containing primary option definitions
     //! \tparam DEFINE_OPTIONS_SECONDARY Functor containing secondary option definitions
     //! \tparam DEFINE_OPTIONS_TERTIARY Functor containing tertiary option definitions
     //! \param include_first_primary Whether to include primary parse of first parser object
     //! \param reparsing Whether we are currently doing a re-parse.
     template <typename CHAIN_TYPE,typename DBN_TYPE,
               typename DEFINE_OPTIONS_PRIMARY,typename DEFINE_OPTIONS_SECONDARY,typename DEFINE_OPTIONS_TERTIARY>
     void parse(bool include_first_primary=false, bool reparsing=false) {

          std::string config_file = "";
          if (parser_0->has_key("config-file")) {
               config_file = (*parser_0)["config-file"].as<std::string>();
          }

          // Specifies whether unregistered options are allowed during parsing
          //  - for all preliminary parse phases, this is set to true
          bool allow_unregistered = true;

          // First, parse only purely shorthand options
          // Examples in this phase: --energy opls-vdw   =>   --energy-opls-vdw
          if (include_first_primary)
               parser_1->parse(DEFINE_OPTIONS_PRIMARY(), command_line, allow_unregistered);
          parser_1->parse_selected_super_groups(boost::bind<void>(DEFINE_OPTIONS_SECONDARY(),
                                                                  _1,
                                                                  ProgramOptionParser::Filter(1)),
                                                boost::bind<void>(DEFINE_OPTIONS_TERTIARY(),
                                                                  _1,
                                                                  ProgramOptionParser::Filter(1)),
                                                command_line,
                                                "shorthands");
          if ((*this)["debug"].as<int>())
               std::cout << "1. Expanded command line: " << parser_1->command_line.raw_string << "\n";

          // Since some of the shorthands have secondary expansions, another
          // pass is necessary
          // Examples in this phase: --energy-opls   =>   --energy-opls-clash-fast --energy-opls-non-bonded ...
          parser_2->parse(DEFINE_OPTIONS_PRIMARY(), command_line, allow_unregistered,
                          config_file);
          parser_2->parse_selected_super_groups(boost::bind<void>(DEFINE_OPTIONS_SECONDARY(),
                                                                  _1,
                                                                  ProgramOptionParser::Filter(1)),
                                                boost::bind<void>(DEFINE_OPTIONS_TERTIARY(),
                                                                  _1,
                                                                  ProgramOptionParser::Filter(1)),
                                                command_line,
                                                "shorthands",
                                                config_file);

          if ((*this)["debug"].as<int>())
               std::cout << "2. Expanded command line: " << parser_2->command_line.raw_string << "\n";

          // Now all options are parsed
          parser_3->parse(DEFINE_OPTIONS_PRIMARY(), command_line, allow_unregistered,
                          config_file);
          parser_3->parse(boost::bind<void>(DEFINE_OPTIONS_SECONDARY(),
                                            _1,
                                            ProgramOptionParser::Filter(1)),
                          boost::bind<void>(DEFINE_OPTIONS_TERTIARY(),
                                            _1,
                                            ProgramOptionParser::Filter(1)),
                          command_line,
                          allow_unregistered,
                          config_file);

          if ((*this)["debug"].as<int>())
               std::cout << "3. Expanded command line: " << parser_3->command_line.raw_string << "\n";

          // Check whether a complete reparse is necessary - this happens in case default values are used for
          // super group defaults
          if (!reparsing) {

               bool reparse_necessary = parser_3->apply_super_group_defaults();

               if (reparse_necessary) {

                    delete parser_1;
                    delete parser_2;
                    delete parser_3;
                    this->parser_1 = new ProgramOptionParser(command_line);
                    this->parser_2 = new ProgramOptionParser(command_line);
                    this->parser_3 = new ProgramOptionParser(command_line);
                    bool include_first_primary=true;
                    bool reparsing = true;
                    parse<CHAIN_TYPE,DBN_TYPE,
                          DEFINE_OPTIONS_PRIMARY,
                          DEFINE_OPTIONS_SECONDARY,
                          DEFINE_OPTIONS_TERTIARY>(include_first_primary, reparsing);
                    return;
               }
          }

          // All options parsed again
          parser_4->parse(DEFINE_OPTIONS_PRIMARY(), command_line, allow_unregistered,
                          config_file);
          parser_4->parse(boost::bind<void>(DEFINE_OPTIONS_SECONDARY(),
                                            _1,
                                            ProgramOptionParser::Filter(1)),
                          boost::bind<void>(DEFINE_OPTIONS_TERTIARY(),
                                            _1,
                                            ProgramOptionParser::Filter(1)),
                          command_line,
                          allow_unregistered,
                          config_file);

          if ((*this)["debug"].as<int>())
               std::cout << "4. Expanded command line: " << parser_4->command_line.raw_string << "\n";

          // The full command line has now been constructed, and parser_4.option_occurrences
          // contains the number of occurrences of each option. We construct the 5th parser
          // based on these values, and parse again. Since this is the final pass unregistered options 
          // no longer allowed
          parser_5->parse(DEFINE_OPTIONS_PRIMARY(), command_line, allow_unregistered,
                          config_file);
          parser_5->parse(boost::bind<void>(DEFINE_OPTIONS_SECONDARY(),
                                            _1,
                                            ProgramOptionParser::Filter(parser_4->option_occurrences, 1)),
                          boost::bind<void>(DEFINE_OPTIONS_TERTIARY(),
                                            _1,
                                            ProgramOptionParser::Filter(parser_4->option_occurrences, 1)),
                          command_line, !allow_unregistered, config_file);

          if ((*this)["debug"].as<int>())
               std::cout << "5. Expanded command line: " << parser_5->command_line.raw_string << "\n";


          // Set flag specifying that options are now fully parsed
          fully_parsed = true;

          // Check for unknown mode
          std::string mode = "";
          if (has_option("mode"))
               mode = (*this)["mode"].as<std::string>();
          if (mode != "" && 
              std::find(phaistos_modes.begin(), phaistos_modes.end(), mode) == phaistos_modes.end()) {
               std::cerr << "Unknown mode: " << mode << "\n";
               exit(1);
          }

          // Output help message if requested
          if (empty_command_line || this->has_option("help")) {
               std::cout << this->generate_help_output() << "\n";
               exit(1);
          }

     }

     //! Constructor
     //!
     //! \param argc Number of command line arguments
     //! \param argv Array of arguments
     Options(int argc, char* argv[])
          : fully_parsed(false),
            empty_command_line(false) {

          // Construct parsers
          parser_0 = new ProgramOptionParser(command_line);
          parser_1 = new ProgramOptionParser(command_line);
          parser_2 = new ProgramOptionParser(command_line);
          parser_3 = new ProgramOptionParser(command_line);
          parser_4 = new ProgramOptionParser(command_line);
          parser_5 = new ProgramOptionParser(command_line);

          // Set command line
          std::string command_line_str;
          for (int i=0; i<argc; ++i) {
               command_line_str += std::string(argv[i]) + " ";
          }
          command_line(command_line_str);

          // Check if there any command line options have been given
          if (argc == 1) {
               empty_command_line = true;
          }
     }


     //! Destructor
     ~Options() {
          delete parser_0;
          delete parser_1;
          delete parser_2;
          delete parser_3;
          delete parser_4;
          delete parser_5;
     }

     //! Overload [] indexing operator (using option name)
     //!
     //! \param option_name Name of option
     //! \return Option value object
     OptionValue operator[](const std::string &option_name) const {
          if (fully_parsed)
               return (*parser_5)[option_name];
          else
               return (*parser_1)[option_name];
     }

     //! Determine whether option is present
     //!
     //! \return True if option is present
     bool has_option(std::string option_name) const {
          if (fully_parsed)
               return parser_5->option_map.count(option_name);
          else
               return parser_1->option_map.count(option_name);
     }

     //! Check whether supergroup is initialized
     //! 
     //! \param name Name of supergroup
     //! \return True if supergroup is initialized
     bool check_super_group_uninitialized(std::string name) const {
          if (fully_parsed)
               return parser_5->check_all_shorthand_uninitialized_in_super_group(name);
          else
               return parser_1->check_all_shorthand_uninitialized_in_super_group(name);
     }     

     //! Retrieve settings object corresponding to option value and index
     //!
     //! \param option_value Option value object
     //! \param index Determines version of option - in case option is defined multiple times.
     //! \return Settings object
     template <typename SETTINGS_TYPE>
     SETTINGS_TYPE &get_settings(const OptionValue &option_value, unsigned int index=0) const {
          const std::vector<boost::shared_ptr<Settings> > &vec = (fully_parsed
                                                                  ? parser_5->settings_objects[option_value.name]
                                                                  : parser_1->settings_objects[option_value.name]);
          if (index >= vec.size() || !vec[index]) {
               std::cerr << "Error: Requested Settings object at index " << index << " not found. Exiting\n";
               assert(false);
          }
          return *dynamic_cast<SETTINGS_TYPE*>(&*vec[index]);
     }

     //! Retrieve settings object corresponding to option value and index
     //!
     //! \param option_value Option value object
     //! \param index Determines version of option - in case option is defined multiple times.
     //! \return Settings object
     template <typename SETTINGS_TYPE>
     SETTINGS_TYPE *get_settings_pointer(const OptionValue &option_value, unsigned int index=0) const {
          const std::vector<boost::shared_ptr<Settings> > &vec = parser_5->settings_objects[option_value.name];
          if (index >= vec.size() || !vec[index]) {
               std::cerr << "Error: Requested Settings object at index " << index << " not found. Exiting\n";
               assert(false);
          }
          return dynamic_cast<SETTINGS_TYPE*>(&*vec[index]);
     }


     //! Generate help output
     std::string generate_help_output() {
          std::string output = "";
          output += "####################### PHAISTOS HELP ######################\n";
          output += "#                                                          #\n";
          output += "#      All possible phaistos command line options are      #\n";
          output += "#                     displayed below.                     #\n";
          output += "#                                                          #\n";
          output += "#   Note that in some situations, some of these options    #\n";
          output += "#      might occur multiple times, using a suffix to       #\n";
          output += "#                   distinguish them.                      #\n";
          output += "#                                                          #\n";
          output += "#    Example:                                              #\n";
          output += "#       --move sidechain-rotamer[debug:10]                 #\n";
          output += "#              sidechain-rotamer[debug:0]                  #\n";
          output += "#                                                          #\n";
          output += "#    Will be expanded to:                                  #\n";
          output += "#       --move-sidechain-rotamer 2                         #\n";
          output += "#       --move-sidechain-rotamer-debug 10                  #\n";
          output += "#       --move-sidechain-rotamer-2-debug 0                 #\n";
          output += "#                                                          #\n";
          output += "############################################################\n";
          output += parser_1->generate_help_output();
          return output;
     }


     //! Generate output containing current option values
     std::string generate_current_options_output() {
          std::string output = "";
          output += "##################### PHAISTOS OPTIONS #####################\n";
          output += "#                                                          #\n";
          output += "#     This output can be directly copy&pasted into a       #\n";
          output += "#   configuration file, and applied using --config-file.   #\n";
          output += "#                                                          #\n";
          output += "############################################################\n";
          std::ostringstream o;

          ProgramOptionParser *output_parser = parser_5;
          if (!fully_parsed)
               output_parser = parser_1;

          // Call output of main parser
          // This parser will output options of all standard types - any
          // user-defined types (such as enums) should be passed along in the type list
          o << output_parser->generate_output<boost::mpl::list<
#if MODULE_SUPPORT
#ifdef HAVE_MUNINNLIB
                                              Muninn::GeEnum,
                                              Muninn::StatisticsLogger::Mode,
#endif
#endif
                                              definitions::AtomSelectionEnum,
                                              TransitionEmissionState,
                                              MoveSidechainLocalMode,
                                              MoveBackboneDBNResampleMode
                                              > >();
          output += o.str();
          output += "################### PHAISTOS OPTIONS END ###################";
          return output;
     }

     //! Output operator
     friend std::ostream & operator<<(std::ostream &o, Options &options) {
          o << options.generate_current_options_output();
          return o;
     }

};

#if MODULE_SUPPORT
#include "modules/phaistos_cpp/global_definitions.cpp"
#endif

//! A collection of random number engines
//! First engine is global engine, other engines are constructed from the global engine
class RandomNumberEngineCollection {
public:

     //! Vector of random number generators
     std::vector<RandomNumberEngine *> random_number_engines;

     //! Constructor
     //!
     //! \param ncopies Number of random number engines in collection
     //! \param identical_copies Whether all copies should be identical
     RandomNumberEngineCollection(int ncopies=1,
                                  bool identical_copies=false) {

          // Create vector of size #threads 
          random_number_engines.resize(ncopies);

          // First random number generator is global constant
          random_number_engines[0] = &random_global;

          // Generate remaining generators from first one
          for (int i=1; i<ncopies; i++) {
               random_number_engines[i] = new RandomNumberEngine();
          }

          // Reseed all copies
          reseed(identical_copies);
     }

     //! Destructor
     ~RandomNumberEngineCollection() {
          // Delete all random number engines except the first
          for (unsigned int i=1; i<random_number_engines.size(); i++) {
               delete random_number_engines[i];
          }
     }

     //! Reseed copies based on first copy
     void reseed(bool identical_copies=false) {

          // Generate remaining generators from first one
          for (unsigned int i=1; i<random_number_engines.size(); i++) {
               if (identical_copies) {
                    *(random_number_engines[i]) = *(random_number_engines[0]);
               } else {
                    random_number_engines[i]->seed(*random_number_engines[0]);
               }
          }
     }

     //! Dereference operator
     std::vector<RandomNumberEngine *> &operator*() {
          return random_number_engines;
     }

     //! Access a specific random number engine
     RandomNumberEngine *operator[](int index) {
          return random_number_engines[index];
     }
};




//! Initialize Dynamics Bayesian Network object
//!
//! \param options Command line options object
//! \param dbn Target dynamic bayesian network
//! \param pdb_filename PDB file used for initializing DBN (empty string means no initialization is done)
//! \param random_number_engines Vector of random number generators
template<typename DBN_TYPE>
void initialize_dbn(Options &options, 
                    DBN_TYPE **dbn,
                    std::string pdb_filename="",
                    std::vector<RandomNumberEngine *> *random_number_engines=NULL) {

     Options::OptionValue option;
     
     // Dbn
     option = options[std::string("backbone-dbn-")+DBN_TYPE::name];
     if (option.occurrences()) {
     
          // Settings typedef
          typedef typename DBN_TYPE::Settings Settings;
          Settings settings = options.get_settings<Settings>(option);

          // Insert global settings if not already set
          if (settings.initial_pdb_file == "" && pdb_filename != "") {
               settings.initial_pdb_file = pdb_filename;
          }
          if (settings.initial_aa_file == "" && options.has_option("aa-file")) {
               settings.initial_aa_file = options["aa-file"].as<std::string>();
          }
          if (settings.initial_ss_file == "" && options.has_option("ss-file")) {
               settings.initial_ss_file = options["ss-file"].as<std::string>();
          }

          RandomNumberEngine *random_number_engine;
          if (random_number_engines)
               random_number_engine = (*random_number_engines)[0];
          else
               random_number_engine = &random_global;

          std::string parameter_file_option_key = 
               std::string("backbone-dbn-")+DBN_TYPE::name+"-parameter-file";
          std::string parameter_file_path;
          if (options[parameter_file_option_key].occurrences() && 
	      options[parameter_file_option_key].as<std::string>() != "") {
               parameter_file_path = 
                    options["data-dir"].as<std::string>() + std::string("/backbone_dbn_parameters/") + 
                    options[parameter_file_option_key].as<std::string>();
               // If file doesn't exist, try to use parameter file name as a full path
               if (!file_exists(parameter_file_path)) {
                    parameter_file_path = options[parameter_file_option_key].as<std::string>();
               }
               // If file still doesn't exist, print error message
               if (!file_exists(parameter_file_path)) {
                    std::cerr << "Error: parameter file not found. Tried: \n\t" << 
                         options["data-dir"].as<std::string>() + std::string("/") + 
                         options[parameter_file_option_key].as<std::string>() << "\nand\n\t" << 
                         options[parameter_file_option_key].as<std::string>() << "\n\n";
                    assert(false);
               }

               *dbn = new DBN_TYPE(settings,
                                   random_number_engine,
                                   Parameters(parameter_file_path));
          } else {
               *dbn = new DBN_TYPE(settings,
                                   random_number_engine);
          }     

          // Initialize unobserved sequences
          (*dbn)->sample();

          // std::cout << **dbn << "\n";

          // if using threading, turn the DBN into a BackboneDBN_Parallel object
          if (options.has_option("threads") && options["threads"].as<int>() > 1) {

               (*dbn)->duplicate(options["threads"].as<int>(), *random_number_engines);
          }
     }
}





//! Initialize Molecule chain object
//!
//! \param options Command line options object
//! \param chain Target molecule chain
//! \param dbn Dynamic bayesian network used for initialization
//! \param pdb_filename PDB filename
template<typename CHAIN_TYPE, typename DBN_TYPE>
void initialize_chain(Options &options, CHAIN_TYPE **chain, DBN_TYPE *dbn, std::string pdb_filename) {

     if (*chain==NULL) {
          if ((!options.has_option("init-from-pdb") || options["init-from-pdb"].as<bool>()) &&
              pdb_filename!="") {
               *chain = new CHAIN_TYPE(pdb_filename,
                                       options["atom-types"].as<definitions::AtomSelectionEnum>());
          } else if (dbn->sequence_length > 0) {
               // Read resseq information from pdb file
               std::vector<int> res_seq;
               if (pdb_filename != "") {
                    ProteinData protein_data = read_pdb_input((char *)pdb_filename.data());
                    int pp_index = 0;
                    res_seq = protein_data.get_resseq()[pp_index];
               }
               *chain = new CHAIN_TYPE(dbn->template get_node<typename DBN_TYPE::ANGLE_NODE>()->get_sequence_vector(), 
                                       dbn->template get_node<typename DBN_TYPE::AA_NODE>()->get_sequence_vector(),
                                       res_seq);
          } else {

               std::cerr << "Error: no initialization information given for molecule chain. Aborting\n";
               exit(1);

          }

          
          // When reading in from PDB, report which atoms were added
          if (options["verbose"].as<bool>() && 
              ((!options.has_option("init-from-pdb") || options["init-from-pdb"].as<bool>()) &&
               pdb_filename!="")) {

               // Make backup of chain object
               CHAIN_TYPE chain_backup(**chain);

               (*chain)->add_atoms(options["atom-types"].as<definitions::AtomSelectionEnum>());

               // Report atoms that were added
               for (int i=0; i<(*chain)->size(); ++i) {

                    Residue &residue = (**chain)[i];

                    std::string header = "Atoms added (not found in PDB):\n";
                    std::string output = "";
                    for (AtomIterator<CHAIN_TYPE,definitions::ALL> it(residue); !it.end(); ++it) {
               
                         definitions::AtomEnum atom_type = it->atom_type;

                         if (!chain_backup[i].has_atom(atom_type)) {
                              output += boost::lexical_cast<std::string>(it->atom_type) + "  ";
                         } 
                    }

                    if (output != "") {
                         std::cerr << header
                                   << boost::lexical_cast<std::string>(residue.index) << " "
                                   << boost::lexical_cast<std::string>(residue.residue_type) << ": " 
                                   << output << "\n";
                    }
               }
          } else {

               (*chain)->add_atoms(options["atom-types"].as<definitions::AtomSelectionEnum>());

          }
     }
}



//! Initialize MonteCarlo object
//!
//! \param options Command line options object
//! \param monte_carlo Target MonteCarlo object
//! \param energy Energy object
//! \param energy_secondary Secondary energy object
//! \param move_collection Collection of moves
//! \param random_number_engines Vector of random number generators
template<typename CHAIN_TYPE>
void initialize_monte_carlo(Options &options,
                            MonteCarloBase<CHAIN_TYPE> **monte_carlo,
                            Energy<CHAIN_TYPE> *energy,
                            Energy<CHAIN_TYPE> *energy_secondary,
                            MoveCollection<CHAIN_TYPE> *move_collection,
                            std::vector<RandomNumberEngine *> *random_number_engines) {

     if (*monte_carlo==NULL) {

          Options::OptionValue option;

          if ((option = options["monte-carlo-metropolis-hastings"]).occurrences()) {

               typedef typename MonteCarloMetropolisHastings<CHAIN_TYPE>::Settings Settings;

               MonteCarloMetropolisHastings<CHAIN_TYPE> *monte_carlo_tmp =
                    new MonteCarloMetropolisHastings<CHAIN_TYPE>(move_collection->chain,
                                                                 energy,
                                                                 move_collection,
                                                                 options.get_settings<Settings>(option));

               if (options["threads"].as<int>() > 1) {
                    *monte_carlo =
                         new MonteCarloMultiThread<MonteCarloMetropolisHastings<CHAIN_TYPE> >(monte_carlo_tmp,
                                                                                              options["threads"].as<int>(),
                                                                                              options["steps-per-move"].as<int>(),
                                                                                              *random_number_engines);
               } else {
                    *monte_carlo = monte_carlo_tmp;
               }

          } else if ((option = options["monte-carlo-simulated-annealing"]).occurrences()) {
               // Parse the options and set the McmcSimulatedAnnealing settings object
               typedef typename MonteCarloSimulatedAnnealing<CHAIN_TYPE>::Settings Settings;

               // Setup the McmcSimulatedAnnealing class
               MonteCarloSimulatedAnnealing<CHAIN_TYPE> *monte_carlo_tmp =
                    new MonteCarloSimulatedAnnealing<CHAIN_TYPE>(move_collection->chain,
                                                                 energy,
                                                                 move_collection,
                                                                 options["iterations"].as<PHAISTOS_LONG_LONG>(),
                                                                 options["threads"].as<int>(),
                                                                 options.get_settings<Settings>(option));

               if (options["threads"].as<int>() > 1) {
                    *monte_carlo =
                         new MonteCarloMultiThread<MonteCarloSimulatedAnnealing<CHAIN_TYPE> >(monte_carlo_tmp,
                                                                                              options["threads"].as<int>(),
                                                                                              options["steps-per-move"].as<int>(),
                                                                                              *random_number_engines);
               } else {
                    *monte_carlo = monte_carlo_tmp;
               }
          } else if ((option = options["monte-carlo-greedy-optimization"]).occurrences()) {
               // Parse the options and set the OptimizeGreedy settings object
               typedef typename MonteCarloGreedyOptimization<CHAIN_TYPE>::Settings Settings;

               // Setup the OptimizeGreedy class
               MonteCarloGreedyOptimization<CHAIN_TYPE> *monte_carlo_tmp =
                         new MonteCarloGreedyOptimization<CHAIN_TYPE>(move_collection->chain,
                                                                      energy,
                                                                      move_collection,
                                                                      options["iterations"].as<PHAISTOS_LONG_LONG>(),
                                                                      options.get_settings<Settings>(option));

               if (options["threads"].as<int>() > 1) {
                    *monte_carlo =
                         new MonteCarloMultiThread<MonteCarloGreedyOptimization<CHAIN_TYPE> >(monte_carlo_tmp,
                                                                                              options["threads"].as<int>(),
                                                                                              options["steps-per-move"].as<int>(),
                                                                                              *random_number_engines);
               } else {
                    *monte_carlo = monte_carlo_tmp;
               }

//! Module monte carlo definition code is inserted here
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/monte_carlo.cpp"
#endif

          } else {
               std::cerr << "# ERROR mc-mode '" << options["mc-mode"].as<std::string>() << "' is not supported\n";
               exit(1);
          }
     }
}

//! Functor for initializing Move set
struct InitializeMoves {

     //! Initialization specific for ChainFB chains
     //!
     //! \param options Command line options object
     //! \param chain Molecule chain object
     //! \param dbn Dynamic Bayesian Network object
     //! \param move_collection Target move collection object
     //! \param random_number_engines Vector of random number generators     
     template <typename DBN_TYPE>
     void template_specific_initialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                           MoveCollection<ChainFB> *move_collection,
                                           std::vector<RandomNumberEngine *> *random_number_engines) {
          
          Options::OptionValue option;



          // Crankshaft move
          option = options["move-crankshaft"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveCrankshaft<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveCrankshaft<ChainFB>(chain, 
                                                                     options.get_settings<Settings>(option, i)));
          }


          // Pivot move - no priors
          option = options["move-pivot-local"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MovePivotLocal<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MovePivotLocal<ChainFB>(chain, 
                                                                     options.get_settings<Settings>(option, i)));
          }

          // Pivot move (TorusDBN + Engh-Huber priors) 
          option = options["move-pivot-local-dbn"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename MovePivotLocal<ChainFB,
                                               MovePriorDbn<BondAnglePriorEnghHuber,
                                                            DBN_TYPE> >::Settings Settings;

               // Add move
               move_collection->add_move(new MovePivotLocal<ChainFB,
                                                            MovePriorDbn<BondAnglePriorEnghHuber,
                                                                         DBN_TYPE> >(chain, 
                                                                                              dbn,
                                                                                              options.get_settings<Settings>(option, i)));
          }

          // Pivot move - uniform samples
          option = options["move-pivot-uniform"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MovePivotUniform<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MovePivotUniform<ChainFB>(chain, 
                                                                       options.get_settings<Settings>(option, i)));
          }

          // None move
          option = options["move-none"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveNone<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveNone<ChainFB>(chain,
                                                               options.get_settings<Settings>(option, i)));
          }



          // Sidechain move: rotamer
          option = options["move-sidechain-rotamer"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveSidechainRotamer<ChainFB, RotamerLibraryDunbrack<> >::Settings Settings;

               // Add move
               move_collection->add_move(new MoveSidechainRotamer<ChainFB, RotamerLibraryDunbrack<> >(chain, 
                                                                                              options.get_settings<Settings>(option, i)));
          }

          // Sidechain move: uniform
          option = options["move-sidechain-uniform"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveSidechainUniform<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveSidechainUniform<ChainFB>(chain, 
                                                                    options.get_settings<Settings>(option, i)));
          }

          // Sidechain move: local
          option = options["move-sidechain-local"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveSidechainLocal<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveSidechainLocal<ChainFB>(chain, 
                                                                  options.get_settings<Settings>(option, i)));
          }

     }


     //! Initialization specific for ChainCA chains
     //!
     //! \param options Command line options object
     //! \param chain Molecule chain object
     //! \param dbn Dynamic Bayesian Network object
     //! \param move_collection Target move collection object
     //! \param random_number_engines Vector of random number generators     
     void template_specific_initialization(const Options &options, ChainCA *chain, FB5DBN *dbn,
                                           MoveCollection<ChainCA> *move_collection,
                                           std::vector<RandomNumberEngine *> *random_number_engines) {


     }


     //! Functor evaluate function
     //!
     //! \param options Command line options object
     //! \param chain Molecule chain object
     //! \param dbn Dynamic Bayesian Network object
     //! \param move_collection Target move collection object
     //! \param random_number_engines Vector of random number generators     
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     void operator()(Options &options,
                     CHAIN_TYPE *chain,
                     DBN_TYPE *dbn,
                     MoveCollection<CHAIN_TYPE> *move_collection,
                     std::vector<RandomNumberEngine *> *random_number_engines) {

          /////////////////////
          // General options //
          /////////////////////

          Options::OptionValue option;

          // TorusDBN move
          option = options["move-backbone-dbn"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings for MoveBackboneDBN
               typedef typename MoveBackboneDBN<CHAIN_TYPE,DBN_TYPE>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveBackboneDBN<CHAIN_TYPE, 
                                         DBN_TYPE>(chain, dbn, options.get_settings<Settings>(option, i)));
          }

//! Module move definition code is inserted here
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/move.cpp"
#endif

          // Move - initializer. This move resamples the entire chain
          // it is used by the MCMC scheme at certain intervals
          if (options["init-from-pdb"].as<bool>() &&
              options["pdb-file"].as<std::string>() != "") {

               // Starting in native
               definitions::AtomSelectionEnum atom_types = options["atom-types"].as<definitions::AtomSelectionEnum>();
               CHAIN_TYPE native_chain(options["pdb-file"].as<std::string>(),
                                       atom_types);
               native_chain.add_atoms(atom_types);
               MoveFixedStructure<CHAIN_TYPE> *move_initializer_native =
                    new MoveFixedStructure<CHAIN_TYPE>(chain, native_chain);
               move_collection->set_initializer(move_initializer_native);          
          } else {

               // Settings for MoveBackboneDBN initializer
               typename MoveBackboneDBN<CHAIN_TYPE,DBN_TYPE>::Settings settings;

               // Always resample entire chain
               settings.move_length_min = chain->size();
               settings.move_length_max = chain->size();

               MoveBackboneDBN<CHAIN_TYPE, DBN_TYPE> *move_DBN_initializer = new MoveBackboneDBN<CHAIN_TYPE, DBN_TYPE>(chain, dbn,
                                                                                                         settings);

               move_collection->set_initializer(move_DBN_initializer);          
          }


          // Call template specific initialization code
          template_specific_initialization(options, chain, dbn, move_collection, random_number_engines);

     }
} initialize_moves;


//! Functor for initializing Energies
struct InitializeEnergy {

     //! Initialization specific for ChainFB chains
     //!
     //! \param options Command line options object
     //! \param chain Molecule chain object
     //! \param dbn Dynamic Bayesian Network object
     //! \param energy Target energy object
     //! \param random_number_engines Vector of random number generators     
     //! \param prefix Optional option prefix (so that both energy and secondary energy can be handled here)
     template <typename DBN_TYPE>
     void template_specific_initialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                           Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_engines, 
                                           std::string prefix="") {

          Options::OptionValue option;
     }


     //! Initialization specific for ChainCA chains
     //!
     //! \param options Command line options object
     //! \param chain Molecule chain object
     //! \param dbn Dynamic Bayesian Network object
     //! \param energy Target energy object
     //! \param random_number_engines Vector of random number generators     
     //! \param prefix Optional option prefix (so that both energy and secondary energy can be handled here)
     void template_specific_initialization(const Options &options, ChainCA *chain, FB5DBN *dbn,
                                           Energy<ChainCA> *energy, std::vector<RandomNumberEngine *> *random_number_engines, std::string prefix) {
     }


     //! Functor evaluate function
     //!
     //! \param options Command line options object
     //! \param chain Molecule chain object
     //! \param dbn Dynamic Bayesian Network object
     //! \param energy Target energy object
     //! \param random_number_engines Vector of random number generators     
     //! \param prefix Optional option prefix (so that both energy and secondary energy can be handled here)
     //! \param energy_main Main energy object - for optional reference when initializing other energy terms
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     void operator()(Options &options, CHAIN_TYPE *chain,
                     DBN_TYPE *dbn, 
                     Energy<CHAIN_TYPE> *energy,
                     std::vector<RandomNumberEngine *> *random_number_engines,
                     std::string prefix="energy",
                     Energy<CHAIN_TYPE> *energy_main=NULL) {

          /////////////////////
          // General options //
          /////////////////////

          Options::OptionValue option;

          std::map<std::string,bool> duplicate_setting_option_map;

          // Check whether energy is an observable_collection
          // if so: optionally add @energy sum term
          ObservableCollection<CHAIN_TYPE> *observables = dynamic_cast<ObservableCollection<CHAIN_TYPE>* >(energy);
          if (observables){
               if (options.has_option("observable-@energy-sum") && energy_main) {

                    option = options["observable-@energy-sum"];
                    for (int i=0; i<option.occurrences(); ++i) {
                         // Settings typedef
                         typedef typename TermEnergySum<CHAIN_TYPE>::Settings Settings;

                         // Add energy term
                         energy->add_term(new TermEnergySum<CHAIN_TYPE>(chain));
                    }
               }

               if (options.has_option("observable-@energy-terms") && energy_main) {

                    option = options["observable-@energy-terms"];

                    duplicate_setting_option_map.insert(std::make_pair("observable-@energy-terms", true));

                    for (unsigned int j=0; j<energy_main->terms.size(); ++j) {

                         for (int i=0; i<option.occurrences(); ++i) {
                              // Settings typedef
                              typedef typename TermEnergyTermWrapper<CHAIN_TYPE>::Settings Settings;

                              // Add energy term
                              energy->add_term(new TermEnergyTermWrapper<CHAIN_TYPE>(chain, j));
                         }
                    }
               }
          }

     
          // Clash-fast
          option = options[prefix+"-clash-fast"];
          for (int i=0; i<option.occurrences(); ++i) {
     
               // Settings typedef
               typedef typename TermClashFast<CHAIN_TYPE>::Settings Settings;
     
               // Add energy
               energy->add_term(new TermClashFast<CHAIN_TYPE>(chain,
                                                              options.get_settings<Settings>(option, i)));

          }


//! Module energy definition code is inserted here
#if MODULE_SUPPORT
#include "modules/phaistos_cpp/energy.cpp"
#endif
     
          // Radius of gyration
          option = options[prefix+"-rg"];
          for (int i=0; i<option.occurrences(); ++i) {
     
               // Settings typedef
               typedef typename TermRg<CHAIN_TYPE>::Settings Settings;
     
               // Add energy
               energy->add_term(new TermRg<CHAIN_TYPE>(chain,
                                                       options.get_settings<Settings>(option, i)));
          }
     
     
          // Root-mean-square-deviation
          option = options[prefix+"-rmsd"];
          for (int i=0; i<option.occurrences(); ++i) {
     
               // Settings typedef
               typedef typename TermRmsd<CHAIN_TYPE>::Settings Settings;
     
               // Add energy
               energy->add_term(new TermRmsd<CHAIN_TYPE>(chain,
                                                        options.get_settings<Settings>(option, i)));
          }
     
     
          // Backbone DBN as an energy term
          option = options[prefix+"-backbone-dbn"];
          for (int i=0; i<option.occurrences(); ++i) {
     
               // Settings typedef
               typedef typename TermBackboneDBN<CHAIN_TYPE,DBN_TYPE>::Settings Settings;
     
               // Add energy
               energy->add_term(new TermBackboneDBN<CHAIN_TYPE,DBN_TYPE>(chain, dbn,
                                                                      options.get_settings<Settings>(option, i)));
          }
     
          // contactmap
          option = options[prefix+"-contact-map"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermContactMap<CHAIN_TYPE>::Settings Settings;

               // Add energy term
               energy->add_term(new TermContactMap<CHAIN_TYPE>(chain,
                                                               options.get_settings<Settings>(option,i)));
          }

          // helix content
          option = options[prefix+"-helix-content"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermHelixContent<CHAIN_TYPE>::Settings Settings;

               // Add energy term
               energy->add_term(new TermHelixContent<CHAIN_TYPE>(chain,
                                                                 options.get_settings<Settings>(option,i)));
          }

          // q-factor
          option = options[prefix+"-q-factor"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermQFactor<CHAIN_TYPE>::Settings Settings;

               // Add energy term
               energy->add_term(new TermQFactor<CHAIN_TYPE>(chain,
                                                            options.get_settings<Settings>(option,i)));
          }

          // angle histogram
          option = options[prefix+"-angle-histogram"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermAngleHistogram<CHAIN_TYPE>::Settings Settings;

               // Add energy term
               energy->add_term(new TermAngleHistogram<CHAIN_TYPE>(chain,
                                                                   options.get_settings<Settings>(option,i)));
          }
          
          // angle
          option = options[prefix+"-angle"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermAngle<CHAIN_TYPE>::Settings Settings;

               // Add energy term
               energy->add_term(new TermAngle<CHAIN_TYPE>(chain,
                                                          options.get_settings<Settings>(option,i)));
          }

          // PDB output
          option = options[prefix+"-pdb"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef (note that this is an observable)
               typedef typename Observable<TermPdb<CHAIN_TYPE> >::Settings Settings;

               // Add energy term (note that this is an observable)
               energy->add_term(new Observable<TermPdb<CHAIN_TYPE> >(chain,
                                                                     options.get_settings<Settings>(option,i)));
          }

          // Call template specific initialization code
          template_specific_initialization(options, chain, dbn, energy, random_number_engines,  prefix);          


          // In case of observables: iterate over all energy terms and 
          // convert them to clas Observable (which inherits from EnergyTerm)
          if (observables) {

               for (unsigned int i=0; i<observables->terms.size(); ++i) {
                    EnergyTerm<CHAIN_TYPE> *energy_term = (*observables)[i];
                    std::string option_name = prefix+"-"+energy_term->name;
                    int occurrence_index = 0;
                    // Test if option name exists - if not, try to remove suffix
                    // (numerical suffix is added in case of multiple
                    if (!options.has_option(option_name)) {
                         int last_dash_index = option_name.find_last_of("-");
                         occurrence_index = boost::lexical_cast<int>(option_name.substr(last_dash_index+1, std::string::npos));
                         option_name = option_name.substr(0, last_dash_index);

                         // If this type of option allows duplicate settings objects, reset occurrence_index to 0
                         if (duplicate_setting_option_map.count(option_name) > 0) {
                              occurrence_index = 0;
                         }
                    }
                    if (options.has_option(option_name)) {               
                         option = options[option_name];

                         ObservableBase::Settings *settings = 
                              options.get_settings_pointer<ObservableBase::Settings>(option, occurrence_index);

                         if (settings) {
                              // observables->terms[i] = observables->terms[i]->clone_to_observable(*settings);
                              int thread_index = 0;
                              observables->convert_term_to_observable(i, *settings, thread_index, energy_main);
                         } else {
                              observables->convert_term_to_observable(i);
                         }

                    } else {
                         std::cerr << "Energy term name - option name mismatch: " << option_name << " not an option name\n";
                         std::exit(1);
                    }
               }
          }

     }

} initialize_energy;


}

#endif
