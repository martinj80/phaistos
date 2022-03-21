// typhon.cpp --- A program for the phaistos package. Runs a protein Monte Carlo dynamics simulation
// Copyright (C) 2011  Tim Harder, Martin Paluszewski, Thomas Hamelryck, Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#include <cstdlib>
#include <cstring>
#include <time.h>
#include <limits>
#include <exception>

#ifdef PHAISTOS_VERSION
#ifdef SVN_REVISION
#include "revision.h"
#endif
#endif

#include "typhon_options.h"
#include "utils/random.h"
#include "energy/energy.h"
#include "energy/term_clash_fast.h"
#include "energy/term_constrain_distances.h"
#include "energy/term_backbone_dbn.h"
#include "energy/term_opls_torsion.h"
#include "energy/term_opls_angle_bend.h"
#include "energy/term_opls_charge.h"
#include "energy/term_opls_vdw.h"
#include "energy/term_gbsa.h"

#include "energy/observable.h"
#include "energy/observable_collection.h"

#include "moves/move_backbone_dbn.h"
#include "moves/move_crisp.h"
#include "moves/move_cra.h"
#include "moves/move_bgs.h"
#include "moves/move_fixed_structure.h"

#include "monte_carlo/monte_carlo_metropolis_hastings.h"

#include "protein/pdb_input.h"
#include "git.h"

// sidechain moves
#include "moves/move_sidechain_dbn.h"
#include "moves/move_sidechain_dbn_local.h"


// Import Phaistos namespace
using namespace phaistos;

//! Map temperature to 1/k
double temperature_to_one_over_k(const double temperature) {
     // 3.2976*10^-27 kcal/K * 6.022*10^23 mol^-1 * temperature
     double one_over_kt = 1.0/(3.2976E-27*6.022E23 * temperature);
     return one_over_kt;
}

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
          for (int i = 1; i < ncopies; i++) {
               random_number_engines[i] = new RandomNumberEngine();
          }

          // Reseed all copies
          reseed(identical_copies);
     }

     // Destructor
     ~RandomNumberEngineCollection() {
          // Delete all random number engines except the first
          for (unsigned int i=1; i<random_number_engines.size(); i++) {
               delete random_number_engines[i];
          }
     }

     //! Reseed copies based on first copy
     void reseed(bool identical_copies = false) {

          // Generate remaining generators from first one
          for (unsigned int i = 1; i < random_number_engines.size(); i++) {
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

     //! Access a specific copy
     RandomNumberEngine *operator[](int index) {
          return random_number_engines[index];
     }
};

//! Main function: Conduct simulation
//!
//! \param options Command line options object
//! \param dbn Dynamic Bayesian Network object
//! \param chain Molecule chain object
//! \param monte_carlo MonteCarlo object
//! \param energy_function Energy object
//! \param move_collection Move set
void shake(TyphonOptions &options, TorusDBN *dbn, ChainFB *chain, 
           MonteCarlo<ChainFB> *monte_carlo, Energy<ChainFB> *energy_function, 
           MoveCollection<ChainFB> *move_collection) {

     ChainFB native(options.pdb_file, definitions::ALL_PHYSICAL_ATOMS - definitions::NON_BACKBONE_H_ATOMS);

     // dump the initial structure
     ChainFB *dchain = monte_carlo->chain;
     char filename[200];
     sprintf(filename, "%s/sample_%012lld_%d_%04f.pdb", options.output_directory.c_str(), monte_carlo->iteration_counter, monte_carlo->thread_index, 0.);
     dchain->output_as_pdb_file(filename, &native);

     int logoutput_interval = options.dump_interval;
     if (options.verbose && options.dump_interval > 10000)
      	 logoutput_interval=10000;
     if (options.debug)
       	 logoutput_interval=500;

     // dump the git as well if requested
     if (options.dump_git) {
          char filenameGIT[200];
          sprintf(filenameGIT, "%s/sample_%012lld_%d_%04f.gitvec", options.output_directory.c_str(), monte_carlo->iteration_counter,
                    monte_carlo->thread_index, 0.);
          Git git(filenameGIT);
          git.generate_gauss_integrals(*dchain, std::string(filename));
     }

     while (monte_carlo->iteration_counter < (unsigned int) options.iterations) {
          // Main loop
          // Make MCMC move using move collection
          monte_carlo->move();

          if (options.dump_interval != 0 && monte_carlo->iteration_counter % options.dump_interval == 0 && 
              ((options.reinitialize_interval == 0 &&  monte_carlo->iteration_counter > (uint) abs(options.burnin)) ||
               (options.reinitialize_interval > 0 && (monte_carlo->iteration_counter % options.reinitialize_interval) > (uint) abs(options.burnin)))) {
               // lets get our chain
               ChainFB *chain = monte_carlo->chain;
               double energy_value = monte_carlo->energy_current;

               char filename[200];
               sprintf(filename, "%s/sample_%012lld_%d_%04f.pdb", options.output_directory.c_str(), monte_carlo->iteration_counter, monte_carlo->thread_index,
                         energy_value);

               chain->output_as_pdb_file(filename, &native);

               if (options.dump_git) {
                    char filenameGIT[200];
                    sprintf(filenameGIT, "%s/sample_%012lld_%d_%04f.gitvec", options.output_directory.c_str(), monte_carlo->iteration_counter,
                              monte_carlo->thread_index, energy_value);

                    Git git(filenameGIT);
                    git.generate_gauss_integrals(*chain, std::string(filename));
               }
          }
          //if (options.dump_interval != 0 && monte_carlo->iteration_counter % options.dump_interval == 0)
          //        monte_carlo->reinitialize();

          // Output energy and move statistics
          if (options.verbose && (monte_carlo->iteration_counter % logoutput_interval == 0)) {
               std::cout << "\n# Iteration: " << monte_carlo->iteration_counter << "\n";
               std::cout << "    Monte Carlo:\n";
               std::cout.flush();
               std::cout << monte_carlo->get_statistics() << "\n";
               std::cout.flush();
               std::cout << "    Energy:\n";
               std::cout.flush();
               std::cout << monte_carlo->get_energy_data() << "\n";
               std::cout.flush();
               std::cout << "    Move statistics:\n";
               std::cout.flush();
               std::cout << monte_carlo->get_move_statistics() << "\n";
               std::cout.flush();
          }
     }
}

//! Initialize MonteCarlo object
//!
//! \param options Command line options object
//! \param monte_carlo Target MonteCarlo object
//! \param energy Energy object
//! \param move_collection Collection of moves
//! \param random_number_engines Vector of random number generators
void initialize_monte_carlo(TyphonOptions &options, MonteCarlo<ChainFB> **monte_carlo, 
                            Energy<ChainFB> *energy, MoveCollection<ChainFB> *move_collection,
                            std::vector<RandomNumberEngine *> *random_number_engines) {

     if (options.mode_mcmc == "metropolis") {
          MonteCarloMetropolisHastings<ChainFB>::Settings settings;
          settings.declash_on_reinitialize = false;
          settings.reinitialization_interval = options.reinitialize_interval;
          MonteCarloMetropolisHastings<ChainFB> *monte_carlo_tmp = new MonteCarloMetropolisHastings<ChainFB> (move_collection->chain, energy, move_collection, settings);
          *monte_carlo = monte_carlo_tmp;

     } else {
          std::cerr << "Unknown MCMC mode selected! Bailing out.\n";
          assert(false);
     }
}


//! Find loop, turn and coil regions
//!
//! \param options Command line options object
//! \param regions Detected regions 
void init_regions(TyphonOptions &options, std::vector<std::pair<int, int> > &regions) {

     std::vector<int> ss_seq;

     regions.clear();

     if (options.ss_file != "") {
          std::string input = file_to_string(options.ss_file);
          ss_seq = ss_str_to_vec(input);
          // lets make sure our sequences match
     } else {
          // ok we didn't get a ss file, so lets try to
          // read it from the pdb file then
          std::cerr << "# WARNING: Trying to read secondary structure from DSSP. Results may vary!\n";
          ss_seq = get_dssp(options.pdb_file)[0];
          //std::cerr << "ss seq from dssp :\n" << ss_seq << "\n\n";
     }

     int start = -1;
     int end = -1;

     for (int i = 0; i < (int) ss_seq.size(); i++) {
          if (ss_seq[i] == 2 && start < 0) {
               start = i;
          }
          if (ss_seq[i] == 2) {
               end = i;
          }
          if (ss_seq[i] != 2 && abs(end - start) >= 4) {
               regions.push_back(std::pair<int, int>(start, end));

          } else if (ss_seq[i] != 2 && abs(end - start) >= 2) {
               regions.push_back(std::pair<int, int>(start - 1, end - 1));
          }

          if (ss_seq[i] != 2) {
               start = -1;
               end = -1;
          }
          std::cout << i << " " << ss_seq[i] << "\n";
     }

     //regions.push_back (std::pair<int, int>(37, 45));
     if (regions.size() <= 0) {
          std::cerr << "ERROR: Trying to apply semi local move with regions, but could not identify any suitable region. \n";
     }
}


//! Initialize Moves (ChainFB specific code)
//!
//! \param options Command line options object
//! \param chain Molecule chain object
//! \param dbn Dynamic Bayesian Network object
//! \param move_collection Target move collection object
//! \param random_number_engines Vector of random number generators     
void initialize_moves(TyphonOptions &options, ChainFB *chain, TorusDBN *dbn, 
                      MoveCollection<ChainFB> *move_collection, 
                      std::vector<RandomNumberEngine *> *random_number_engines) {

     Move<ChainFB> *local_move;
     ChainFB *refchain;

     if(options.reinitialize_interval > 0){
             refchain = new ChainFB(*chain); // make a copy for reinitialization move
             MoveFixedStructure<ChainFB> *initmove = new MoveFixedStructure<ChainFB>(chain, *refchain);
             move_collection->set_initializer(initmove);
     }
     /////////////////////////////////////////////////////////////////////////////
     // always move the sidechains around please !
     int seed = (*(*random_number_engines)[0])();
     BasiliskDBN *scdbn = new BasiliskDBN(options.sc_dbn_file, seed);

     MoveSidechainDBN<BasiliskDBN>::Settings settings;

     // Default settings
     settings.implicit_energy = true;
     // We need to sample backbone independent here in order
     // to avoid problems when moving the backbone around.
     // In which case the side chain bias needs to be taken
     // into account when moving the BB. By sampling independent
     // we can avoid that whole ordeal.
     settings.ignore_bb = true;
     settings.sample_hydrogen_chis = true;
     settings.move_length_min = 1;
     settings.move_length_max = 1;
     settings.weight = options.sc_move_weight;
     MoveSidechainDBN<BasiliskDBN> *move_sc = new MoveSidechainDBN<BasiliskDBN> (chain, scdbn, settings);
     move_collection->add_move(move_sc);

     bool use_regions = false;

     if (options.use_semi_local_move) {

          //! Create settings object
          MoveBGS<ChainFB,bgs::GaussianStep<MovePrior<BondAnglePriorUninformative, DihedralPriorUninformative> > >::Settings settings;

          if (use_regions) {
               // calculate coil regions, big enought to do a
               // semi local move
               std::vector<std::pair<int, int> > tmp;
               init_regions(options, tmp);

               std::cout << "# Using Semi Local Moves with regions : \n";
               for (int i = 0; i < (int) tmp.size(); i++) {
                    std::cout << i << " : " << tmp[i].first << "  " << tmp[i].second << "\n";
               }
               std::cout << "# end regions \n";
               if (tmp.size() > 0) {
                    settings.regions = tmp;

               } else {
                    std::cerr << "ERROR: Trying to apply semi local move with regions, but could not identify any suitable region .. skipped adding the move\n";
               }
          }

          Move<ChainFB> *move_semi;
          settings.weight = options.semi_local_move_weight;
          move_semi = new MoveBGS<ChainFB, bgs::GaussianStep<MovePrior<BondAnglePriorUninformative,DihedralPriorUninformative> > > (chain,settings);
                                
          move_collection->add_move(move_semi);
     }

     if (options.mode_move == "crisp") {
          // CRISP - dihedral prior: uninformative, bondangle prior: Engh-Huber
          MoveCRISP<ChainFB, crisp::PreRotation<MovePrior<BondAnglePriorUninformative, DihedralPriorUninformative> >,
                    crisp::PostRotation>::Settings settings;

          settings.weight = 1.0;
          local_move = new MoveCRISP<ChainFB, crisp::PreRotation<MovePrior<BondAnglePriorUninformative, DihedralPriorUninformative> >,
                    crisp::PostRotation> (chain, settings);
          move_collection->add_move(local_move);

     } else if (options.mode_move == "crisp-with-prior") {
          MoveCRISP<ChainFB, crisp::PreRotation<MovePriorDbn<BondAnglePriorEnghHuber, TorusDBN> >,
                    crisp::PostRotation>::Settings settings;

          settings.weight = 1.0;
          local_move = new MoveCRISP<ChainFB, crisp::PreRotation<MovePriorDbn<BondAnglePriorEnghHuber, TorusDBN> >,
                    crisp::PostRotation> (chain, dbn, settings);
          move_collection->add_move(local_move);

     } else {
          std::cerr << "Unknown move mode selected! Bailing out.\n";
          assert(false);
     }

}

//! Initialize Dynamics Bayesian Network object
//!
//! \param options Command line options object
//! \param dbn Target dynamic bayesian network
void initialize_dbn(TyphonOptions &options, TorusDBN **dbn) {

     // Settings typedef
     TorusDBN::Settings settings;

     if (options.pdb_file != "") {
          settings.initial_pdb_file = options.pdb_file;
     }

     if (not options.ignore_ss) {
          if (options.ss_file != "") {
               settings.initial_ss_file = options.ss_file;
          } else {
               // ok we didn't get a ss file, so lets try to
               // read it from the pdb file then
               std::cerr << "# WARNING: Trying to read secondary structure from DSSP. Results may vary!\n";
               std::string ss_seq = split(get_dssp_string(options.pdb_file), " ")[0];
               settings.initial_ss_sequence = ss_seq;
          }
     }

     *dbn = new TorusDBN(settings, &random_global, Parameters(options.dbn_parameter_dir + "/" + options.dbn_parameter_file));

     (*dbn)->sample();
}

//! Initialization energies
//!
//! \param options Command line options object
//! \param chain Molecule chain object
//! \param dbn Dynamic Bayesian Network object
//! \param energy Target energy object
void initialize_energies(TyphonOptions &options, ChainFB *chain, TorusDBN *dbn, Energy<ChainFB> *energy) {
     if (options.mode_energy == "all") {
          //
          // we need some clash energy in both cases
          TermClashFast<ChainFB>::Settings settings_clash_fast;
          settings_clash_fast.only_modified_pairs = true;
          settings_clash_fast.boolean_mode = true;
          energy->add_term(new TermClashFast<ChainFB> (chain, settings_clash_fast), 1.);
          //
          // and make up for the proposal
          energy->add_term(new TermBackboneDBN<ChainFB, TorusDBN> (chain, dbn), 1);

          double one_over_kt = temperature_to_one_over_k(300);
          TermOplsAngleBendCached::Settings term_opls_angle_bend_settings;
          term_opls_angle_bend_settings.omit_sidechains = true;
          energy->add_term(new TermOplsAngleBendCached(chain, term_opls_angle_bend_settings), one_over_kt);
          //
          // constrain distances
          TermConstrainDistances<ChainFB>::Settings settings_constrain_distances;

          if (options.restore_network != "") {
               settings_constrain_distances.network_filename = options.restore_network;
          }

          if (options.ignore_hbonds)
               settings_constrain_distances.include_bb_hbond = false;
          if (options.ignore_bbsc_hbonds)
               settings_constrain_distances.include_bb_sc_hbond = false;
          if (options.ignore_sc_hbonds)
               settings_constrain_distances.include_sc_hbond = false;
          if (options.ignore_gaussian)
               settings_constrain_distances.include_ca_contacts = false;
          if (options.ignore_ssbond)
               settings_constrain_distances.include_ss_bond = false;
          if (options.ignore_fixed)
               settings_constrain_distances.include_fixed_point_contacts = false;
          if (!options.ignore_fixed)
               settings_constrain_distances.include_fixed_point_contacts = true;
          if (options.init_hbond_from_native)
               settings_constrain_distances.init_hbond_from_native = true;
          if (options.no_prune)
              settings_constrain_distances.prune_network = false;

          settings_constrain_distances.verbose = options.verbose;
          settings_constrain_distances.debug = options.debug;
          settings_constrain_distances.generate_pymol = options.generate_pymol;
          settings_constrain_distances.ca_distance = options.ca_cutoff;
          settings_constrain_distances.ca_skip = options.ca_skip;
          settings_constrain_distances.dehydron_bb_cutoff   = options.dehydron_bb_cutoff;
          settings_constrain_distances.dehydron_bb_sc_cutoff = options.dehydron_bbsc_cutoff;
          settings_constrain_distances.dehydron_sc_cutoff   = options.dehydron_sc_cutoff;
          settings_constrain_distances.pdb_file   = options.pdb_file;

          energy->add_term(new TermConstrainDistances<ChainFB> (chain, settings_constrain_distances));

     } else if (options.mode_energy == "opls") {
          double one_over_kt = temperature_to_one_over_k(300);
          energy->add_term(new TermGbsaCached(chain), one_over_kt);
          energy->add_term(new TermOplsChargeCached(chain), one_over_kt);
          energy->add_term(new TermOplsVdwCached(chain), one_over_kt);

          TermOplsAngleBendCached::Settings term_opls_angle_bend_settings;
          term_opls_angle_bend_settings.omit_sidechains = true;
          energy->add_term(new TermOplsAngleBendCached(chain, term_opls_angle_bend_settings), one_over_kt);

          energy->add_term(new TermOplsTorsion(chain), one_over_kt);
     } else {
          std::cerr << "Unknown energy mode selected! Bailing out.\n";
          assert(false);
     }

}

//! Main function
int main(int argc, char *argv[]) {

     /*****************************************************************************************
      *                     Initialization                                                    *
      *****************************************************************************************/

     // Read the command line
     TyphonOptions options(argc, argv);

     if (options.pdb_file == "") {
    	 options.print_usage();
		 return 0 ;
     }

     if (options.verbose) {
          printf("\n#############################################################################################################\n");
          std::cout<< "#    ___________             .__                    " << std::endl;
          std::cout<< "#    \\__    ___/__.__.______ |  |__   ____   ____   " << std::endl;
          std::cout<< "#      |    | <   |  |\\____ \\|  |  \\ /  _ \\ /    \\  " << std::endl;
          std::cout<< "#      |    |  \\___  ||  |_> >   Y  (  <_> )   |  \\ " << std::endl;
          std::cout<< "#      |____|  / ____||   __/|___|  /\\____/|___|  / " << std::endl;
          std::cout<< "#              \\/     |__|        \\/            \\/  " << std::endl << "# " << std::endl;
#ifdef PHAISTOS_VERSION
#ifdef SVN_REVISION          
          printf(     "#                       Version: %5s  Build: %5s \n" , PHAISTOS_VERSION, SVN_REVISION);
#endif
#endif
     }

     // create a place to dump out samples
     int system_result = system((std::string("mkdir -p ") + options.output_directory).c_str());
     if (system_result > 0) {
          std::cerr << "ERROR: Could not create output directory .. bailing out.\n";
          return 256;
     }

     ChainFB *chain = NULL;
     TorusDBN *dbn = NULL;
     MonteCarlo<ChainFB> *monte_carlo = NULL;

     /*******************************************************
      *                     Random                          *
      *******************************************************/
     // Initialize random number generators
     if (options.seed <= 0) {
          // Choose random seed
          options.seed = static_cast<unsigned int> (time(0));
     }
     random_global.seed(options.seed);

     RandomNumberEngineCollection random_number_engines(1, false);

     if (options.verbose) {

          options.print_options();
     }

     /*******************************************************
      *                     DBN                             *
      *******************************************************/
     // get Torus DBN up and running
     initialize_dbn(options, &dbn);

     /*******************************************************
      *                     Chain                           *
      *******************************************************/
     // read in the PDB file and make sure all the atoms are present
     chain = new ChainFB(options.pdb_file, definitions::ALL_PHYSICAL_ATOMS);// - NON_BACKBONE_H_ATOMS);
     // we do this slightly silly stunt here to make sure
     // all the heavy atoms are present ..
     //chain->add_atoms(ALL_PHYSICAL_ATOMS - NON_BACKBONE_H_ATOMS);
     // .. and then add the hydrogens in the phaistos
     // nomenclature, just to be save we have the correct names
     chain->add_atoms(definitions::ALL_PHYSICAL_ATOMS); // - NON_BACKBONE_H_ATOMS);
     chain->check_consistency();
     // if wanted we can dump a structure file
     if (options.create_phaistos_input) {
          char filename[200];
          sprintf(filename, "%s_phaistos.pdb", (options.pdb_file.substr(0,options.pdb_file.size()-4)).c_str());
          chain->output_as_pdb_file(filename);
     }


     /*******************************************************
      *                     Moves                           *
      *******************************************************/
     // initialize the moves
     MoveCollection<ChainFB> move_collection(chain);
     if (not options.mode_analyze)
          initialize_moves(options, chain, dbn, &move_collection, &(*random_number_engines));

     /*******************************************************
      *                     Energies                        *
      *******************************************************/
     // and we would also like to have some energies here.
     //
     // set up an energy container
     Energy<ChainFB> *energy = new Energy<ChainFB> (chain);
     if (options.mode_analyze && options.verbose)
               std::cout << "# ANALYZE \n";
     initialize_energies(options, chain, dbn, energy);

     /*******************************************************
      *                     Monte Carlo                     *
      *******************************************************/
     // Initialize Monte Carlo object
     if (not options.mode_analyze)
          initialize_monte_carlo(options, &monte_carlo, energy, &move_collection, &*random_number_engines);

     /*****************************************************************************************
      *                     Run                                                               *
      *****************************************************************************************/

     if (not options.mode_analyze) {
          // Calculate the native energy
          if (options.verbose) {
               ChainFB native(options.pdb_file, definitions::ALL_PHYSICAL_ATOMS);
               native.add_atoms(definitions::ALL_PHYSICAL_ATOMS);
               native.check_consistency();
               Energy<ChainFB> energy_native(*energy, &native);
               energy_native.evaluate();
               std::cout << "NATIVE ENERGY\n" << energy_native << "\n\n";
          }

          /*                     The Big Dance                   */
          shake(options, dbn, chain, monte_carlo, energy, &move_collection);
     }

     /*****************************************************************************************
      *                     Cleanup                                                           *
      *****************************************************************************************/
     delete chain;
     delete dbn;
     delete monte_carlo;
     delete energy;

     return 0;
}
;

