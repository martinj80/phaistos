namespace module_saxs {

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

          // SAXS Debye
          for (int counter = occurrences[prefix+"-saxs-debye"]; counter > 0; counter--) {

               typedef TermSaxsDebye<ChainFB> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "saxs-debye (" + prefix + ")",
                         prefix+"-saxs-debye", settings,
                         make_vector(
                              make_vector(std::string("saxs-intensities-filename"),
                                          std::string("path to file containing SAXS intensities"),
                                          &settings->saxs_intensities_filename),
                              make_vector(std::string("saxs-form-factors-filename"),
                                          std::string("path to file containing SAXS form factors"),
                                          &settings->saxs_form_factors_filename),
                              make_vector(std::string("exp-errors-alpha"),
                                          std::string("alpha parameter for experimental error term"),
                                          &settings->exp_errors_alpha),
                              make_vector(std::string("exp-errors-beta"),
                                          std::string("beta parameter for experimental error term"),
                                          &settings->exp_errors_beta),
                              make_vector(std::string("q-bins"),
                                          std::string("number of q-bins used in the curve calculation"),
                                          &settings->q_bins),
                              make_vector(std::string("q-bins-first"),
                                          std::string("energy evaluation starts at this q-bin (default: 0). the first bin is used for normalization"),
                                          &settings->q_bins_first),
                              make_vector(std::string("one-body-model"),
                                          std::string("use one-body instead of two-body form factor model"),
                                          &settings->one_body_model),
                              make_vector(std::string("use-sine-lookup-table"),
                                          std::string("use sin() lookup table instead of the full trigonometric evaluation"),
                                          &settings->sine_lookup_table)
                              )),
                    super_group, counter==1);
          }
    }
};

}
