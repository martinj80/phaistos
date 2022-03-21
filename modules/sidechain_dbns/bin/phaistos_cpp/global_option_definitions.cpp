namespace module_sidechain_dbns {

// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct MoveOptions {

     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 CHAIN_TYPE *chain,
                 DBN_TYPE *dbn) {
     }

     //! Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 ChainFB *chain,
                 DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;

          // Define defaults for the different modes
          ModeDefinitions mode_definitions(target, chain);

          // SC Basilisk
          for (int counter = occurrences["move-sidechain-basilisk"]; counter > 0; counter--) {
          
               // Create settings object
               typedef MoveSidechainDBN<BasiliskDBN> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));
          
               // Override Default settings
               settings->implicit_energy = mode_definitions.implicit_energies_allowed;

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";
          
               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Sidechain move Basilisk",
                         "move-sidechain-basilisk", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the dbn bias (implicit energy) should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy),
                              make_vector(std::string("ignore-bb"),
                                          std::string("Use the backbone independent version instead"),
                                          &settings->ignore_bb),
                              make_vector(std::string("reject-broken-prolines"),
                                          std::string("Whether to reject moves that produce prolines with broken rings"),
                                          &settings->reject_broken_prolines),
                              make_vector(std::string("sample-hydrogen-chis"),
                                          std::string("Sample sidechain residual dofs"),
                                          &settings->sample_hydrogen_chis),
                              make_vector(std::string("sample-hydrogen-chis-normal"),
                                          std::string("Resampling constrained to a gaussian around the current state, false=uniform"),
                                          &settings->sample_hydrogen_chis_normal),
                              make_vector(std::string("sample-hydrogen-chis-sigma"),
                                          std::string("Std. of gaussian for resampling residual dofs"),
                                          &settings->sample_hydrogen_chis_sigma)
                              )), super_group, counter==1);
          }
          
          
          // SC Basilisk-local
          for (int counter = occurrences["move-sidechain-basilisk-local"]; counter > 0; counter--) {
          
               // Create settings object
               typedef MoveSidechainDBNLocal<BasiliskDBN> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));
          
               // Override Default settings
               settings->implicit_energy = mode_definitions.implicit_energies_allowed;

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";
          
               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Sidechain move - Basilisk local",
                         "move-sidechain-basilisk-local", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the dbn bias (implicit energy) should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy)
                              )), super_group, counter==1);
          }
          
          // SC Basilisk-multi
          for (int counter = occurrences["move-sidechain-basilisk-multi"]; counter > 0; counter--) {
          
               // Create settings object
               typedef MoveSidechainDBNMulti<BasiliskDBN> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));
          
               // Override Default settings
               settings->implicit_energy = mode_definitions.implicit_energies_allowed;

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";
          
               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Sidechain move - Basilisk Multi",
                         "move-sidechain-basilisk-multi", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the dbn bias (implicit energy) should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy),
                              make_vector(std::string("multi-count"),
                                          std::string("How many side chains to resample in each move"),
                                          &settings->move_count)
                              )), super_group, counter==1);
          }
          
          
          // SC Compas
          for (int counter = occurrences["move-sidechain-compas"]; counter > 0; counter--) {

               typedef MoveSidechainDBN<CompasDBN> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));
          
               // Override Default settings
               settings->implicit_energy = mode_definitions.implicit_energies_allowed;

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";
          
               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "Sidechain move - Compas",
                         "move-sidechain-compas", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the dbn bias (implicit energy) should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy)
                              )), super_group, counter==1);
          }
          
          
     }
};


// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {                    
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;

          // Basilisk explicit energy function
          for (int counter = occurrences[prefix+"-basilisk"]; counter > 0; counter--) {

               // Create settings object
               typedef TermSidechainDBN<ChainFB,BasiliskDBN> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Basilisk explicit energy function (" + prefix + ")",
                         prefix+"-basilisk", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("use-cache"),
                                          std::string("Use the cached version of the term (should be faster)"),
                                          &settings->use_cache),
                              make_vector(std::string("ignore-bb"),
                                          std::string("Use the backbone independent version instead"),
                                          &settings->ignore_bb),
                              make_vector(std::string("eliminate-move-bias"),
                                          std::string("Divide out the move-bias of the corresponding moves. Equivalent to (but faster than) setting implicit-energy to false."),
                                          &settings->eliminate_move_bias)
                              )), super_group, counter==1);
          }

          // Basilisk explicit energy function
          for (int counter = occurrences[prefix+"-compas"]; counter > 0; counter--) {

               typedef TermSidechainDBN<ChainFB,CompasDBN> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Compas explicit energy function (" + prefix + ")",
                         prefix+"-compas", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("use-cache"),
                                          std::string("Use the cached version of the term (should be faster)"),
                                          &settings->use_cache),
                              make_vector(std::string("ignore-bb"),
                                          std::string("Use the backbone independent version instead"),
                                          &settings->ignore_bb),
                              make_vector(std::string("eliminate-move-bias"),
                                          std::string("Divide out the move-bias of the corresponding moves. Equivalent to (but faster than) setting implicit-energy to false."),
                                          &settings->eliminate_move_bias)
                              )), super_group, counter==1);
          }
     }
};

}
