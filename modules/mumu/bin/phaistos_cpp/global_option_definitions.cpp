namespace module_mumu {

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

          // Mumu energy function
          for (int counter = occurrences[prefix+"-mumu"]; counter > 0; counter--) {

               // Create settings object
               typedef TermMumu<ChainFB> EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
               // typedef typename TermMumu<ChainFB>::Settings Settings;
               // boost::shared_ptr<Settings> settings(new Settings());

               // Set mocapy directory value base on phaistos-wide setting
               settings->mocapy_dbn_dir = target["data-dir"].as<std::string>()+"/mocapy_dbns";

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Multibody multinomial (MuMu) model (" + prefix + ")",
                         prefix+"-mumu", settings,
                         make_vector(
                              make_vector(std::string("mocapy-dbn-dir"),
                                          std::string("Path to mocapy model file directory."),
                                          &settings->mocapy_dbn_dir),
                              make_vector(std::string("model-filename"),
                                          std::string("Model filename"),
                                          &settings->model_filename),
                              make_vector(std::string("reference-model-filename"),
                                          std::string("Reference model filename"),
                                          &settings->reference_model_filename),
                              make_vector(std::string("use-ratio"),
                                          std::string("Use the ration method"),
                                          &settings->use_ratio)
                              )), super_group, counter==1);
          }

     }
};

}
