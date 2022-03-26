namespace module_sidechain_dbns {

// Module: energy term initialization
struct MoveInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     MoveInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                                   MoveCollection<CHAIN_TYPE> *move_collection,
                                   std::vector<RandomNumberEngine *> *random_number_generators) {
     }

     // Constructor - template specific case
     template <typename DBN_TYPE>
     MoveInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                   MoveCollection<ChainFB> *move_collection,
                                   std::vector<RandomNumberEngine *> *random_number_generators) {

          Options::OptionValue option;

          // SC move: basilisk
          option = options["move-sidechain-basilisk"];
          for (int i=0; i<option.occurrences(); ++i) {

               /////////////////////////////////////////////////////////////////////////
               //   Basilisk Full Atom Sidechain Move
               //

               // Settings typedef
               typedef MoveSidechainDBN<BasiliskDBN>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               int seed = (*(*random_number_generators)[0])();
               BasiliskDBN *dbn = new BasiliskDBN(settings.mocapy_dbn_dir + "/" + settings.model_filename, seed);
               // Allocate number of DBNs to match number of threads
               if (options["threads"].as<int> () > 1) {
                    dbn->duplicate(options["threads"].as<int> (), *random_number_generators);
               }

               // Add move
               move_collection->add_move(new MoveSidechainDBN<BasiliskDBN>(chain, dbn,
                                                                           settings));
          }


          // SC move: basilisk-local
          option = options["move-sidechain-basilisk-local"];
          for (int i=0; i<option.occurrences(); ++i) {

               /////////////////////////////////////////////////////////////////////////
               //   Basilisk Full Atom Sidechain Move
               //

               // Settings typedef
               typedef MoveSidechainDBNLocal<BasiliskDBN>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               int seed = (*(*random_number_generators)[0])();
               BasiliskDBN *dbn = new BasiliskDBN(settings.mocapy_dbn_dir + "/" + settings.model_filename, seed);
               // Allocate number of DBNs to match number of threads
               if (options["threads"].as<int> () > 1) {
                    dbn->duplicate(options["threads"].as<int> (), *random_number_generators);
               }


               // Add move
               move_collection->add_move(new MoveSidechainDBNLocal<BasiliskDBN>(chain, dbn,
                                                                                settings));
          }


          // SC move: basilisk-multi
          option = options["move-sidechain-basilisk-multi"];
          for (int i=0; i<option.occurrences(); ++i) {

               /////////////////////////////////////////////////////////////////////////
               //   Basilisk Full Atom Sidechain Move
               //

               // Settings typedef
               typedef MoveSidechainDBNMulti<BasiliskDBN>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               int seed = (*(*random_number_generators)[0])();
               BasiliskDBN *dbn = new BasiliskDBN(settings.mocapy_dbn_dir + "/" + settings.model_filename, seed);
               // Allocate number of DBNs to match number of threads
               if (options["threads"].as<int> () > 1) {
                    dbn->duplicate(options["threads"].as<int> (), *random_number_generators);
               }

               // Add move
               move_collection->add_move(new MoveSidechainDBNMulti<BasiliskDBN>(chain, dbn,
                                                                                settings));
          }


          // SC move: Compas (pseudo-sidechains)
          option = options["move-sidechain-compas"];
          for (int i=0; i<option.occurrences(); ++i) {

               /////////////////////////////////////////////////////////////////////////
               //   Compas Pseudo Sidechain Move
               //

               // Settings typedef
               typedef MoveSidechainDBN<CompasDBN>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               int seed = (*(*random_number_generators)[0])();
               CompasDBN *dbn = new CompasDBN(settings.mocapy_dbn_dir + "/" + settings.model_filename, seed);
               // Allocate number of DBNs to match number of threads
               if (options["threads"].as<int> () > 1) {
                    dbn->duplicate(options["threads"].as<int> (), *random_number_generators);
               }

               // Add move
               move_collection->add_move(new MoveSidechainDBN<CompasDBN>(chain, dbn,
                                                                         settings));
          }
     }
};

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                                       Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {
     }

     // Constructor - template specific case
     template <typename DBN_TYPE>
     EnergyInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                       Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {

          Options::OptionValue option;

          // Basilisk explicit energy term
          option = options[prefix + "-basilisk"];
          for (int i = 0; i < option.occurrences(); ++i) {

               /////////////////////////////////////////////////////////////////////////
               //   Basilisk Full Atom Sidechain Energy
               //

               if (not options.has_option(prefix + "-basilisk-model-filename")) {
                    std::cerr << "Cannot locate Basilisk DBN file while trying to initialize the energy function.\n";
                    break;
               }

               // Settings typedef
               typedef typename TermSidechainDBN<ChainFB, BasiliskDBN>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               int seed = (*(*random_number_generators)[0])();
               BasiliskDBN *scDBN = new BasiliskDBN(settings.mocapy_dbn_dir + "/" + settings.model_filename, seed);
               // Allocate number of DBNs to match number of threads
               if (options["threads"].as<int> () > 1) {
                    scDBN->duplicate(options["threads"].as<int> (), *random_number_generators);
               }

               // Add energy term
               energy->add_term(new TermSidechainDBN<ChainFB, BasiliskDBN> (chain, scDBN, settings));
          }


          // Compas explicit energy term
          option = options[prefix + "-compas"];
          for (int i = 0; i < option.occurrences(); ++i) {

               /////////////////////////////////////////////////////////////////////////
               //   Compas Pseudo Sidechain Energy
               //

               if (not options.has_option(prefix + "-compas-model-filename")) {
                    std::cerr << "Cannot locate Compas DBN file while trying to initialize the energy function.\n";
                    break;
               }

               // Settings typedef
               typedef typename TermSidechainDBN<ChainFB, CompasDBN>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               int seed = (*(*random_number_generators)[0])();
               CompasDBN *scDBN = new CompasDBN(settings.mocapy_dbn_dir + "/" + settings.model_filename, seed);
               // Allocate number of DBNs to match number of threads
               if (options["threads"].as<int> () > 1) {
                    scDBN->duplicate(options["threads"].as<int> (), *random_number_generators);
               }

               // Add energy term
               energy->add_term(new TermSidechainDBN<ChainFB, CompasDBN> (chain, scDBN, settings));
          }

     }

};


}
