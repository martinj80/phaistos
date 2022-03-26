namespace module_opls {

//! Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {
     }

     //! Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {
          
          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Define defaults for the different modes
          ModeDefinitions mode_definitions(target, chain);
          
          // OPLS: charge term
          for (int counter = occurrences[prefix+"-opls-charge"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsCharge EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-charge: OPLS charge term (" + prefix + ")",
                         prefix+"-opls-charge", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // OPLS: charge term - cached
          for (int counter = occurrences[prefix+"-opls-charge-cached"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsChargeCached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-charge-cached: OPLS charge term - cached version (" + prefix + ")",
                         prefix+"-opls-charge-cached", settings,
                         make_vector(
                              make_vector(std::string("cutoff-distance"),
                                          std::string("Distance beyond which contributions are set to zero."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->cutoff_distance))
                              )),
                    super_group, counter==1);
          }
          
          
          // OPLS: vdw term
          for (int counter = occurrences[prefix+"-opls-vdw"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsVdw EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-vdw: OPLS van der Waals term (" + prefix + ")",
                         prefix+"-opls-vdw", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // OPLS: vdw term - cached
          for (int counter = occurrences[prefix+"-opls-vdw-cached"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsVdwCached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-vdw-cached: OPLS van der Waals term - cached version (" + prefix + ")",
                         prefix+"-opls-vdw-cached", settings,
                         make_vector(
                              make_vector(std::string("cutoff-distance"),
                                          std::string("Distance beyond which contributions are set to zero."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->cutoff_distance))
                              )),
                    super_group, counter==1);
          }
          
          
          // OPLS: angle-bend term
          for (int counter = occurrences[prefix+"-opls-angle-bend"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsAngleBend EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-angle-bend: OPLS angle bend term (" + prefix + ")",
                         prefix+"-opls-angle-bend", settings,
                         make_vector(
                              make_vector(std::string("omit-sidechains"),
                                          std::string("Exclude sidechain interactions"),
                                          &settings->omit_sidechains)
                              )),
                    super_group, counter==1);
          }
          
          
          // OPLS: angle-bend term - cached
          for (int counter = occurrences[prefix+"-opls-angle-bend-cached"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsAngleBendCached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-angle-bend-cached: OPLS angle bend term - cached version (" + prefix + ")",
                         prefix+"-opls-angle-bend-cached", settings,
                         make_vector(
                              make_vector(std::string("omit-sidechains"),
                                          std::string("Exclude sidechain interactions"),
                                          &settings->omit_sidechains)
                              )),
                    super_group, counter==1);
          }
          
          
          // OPLS: torsion term
          for (int counter = occurrences[prefix+"-opls-torsion"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsTorsion EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-torsion: OPLS torsion term (" + prefix + ")",
                         prefix+"-opls-torsion", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // OPLS: improper torsion term
          for (int counter = occurrences[prefix+"-opls-improper-torsion"]; counter > 0; counter--) {
                                       
               // Create settings object
               typedef TermOplsImptor EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
                                        
               // If temperature has been specified, set weight based on that value   
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
                                        
               // Add options    
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-improper-torsion: OPLS improper-torsion term (" + prefix + ")",
                         prefix+"-opls-improper-torsion", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // OPLS: bond-stretch term
          for (int counter = occurrences[prefix+"-opls-bond-stretch"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsBondStretch EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-bond-stretch: OPLS bond-stretch term (" + prefix + ")",
                         prefix+"-opls-bond-stretch", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // OPLS: all non-bonded terms (including GBSA)
          for (int counter = occurrences[prefix+"-opls-non-bonded"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsNonBonded<double> EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-non-bonded: gbsa, vdw and charge terms (" + prefix + ")",
                         prefix+"-opls-non-bonded", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // OPLS: all non-bonded terms (including GBSA) - cached version
          for (int counter = occurrences[prefix+"-opls-non-bonded-cached"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermOplsNonBondedCached<double> EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "opls-non-bonded-cached: gbsa, vdw and charge terms - cached version (" + prefix + ")",
                         prefix+"-opls-non-bonded-cached", settings,
                         make_vector(
                              make_vector(std::string("vdw-cutoff-distance"),
                                          std::string("Distance beyond which vdw contributions are set to zero."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->vdw_cutoff_distance)),
                              make_vector(std::string("charge-cutoff-distance"),
                                          std::string("Distance beyond which charge contributions are set to zero."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->charge_cutoff_distance)),
                              make_vector(std::string("gbsa-maximum-deviation-cutoff"),
                                          std::string("Maximum deviation allowed in born radii in two subtrees of the chaintree before it is recalculated"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->gbsa_maximum_deviation_cutoff)),
                              make_vector(std::string("gbsa-cutoff-distance-phase1"),
                                          std::string("Distance beyond which gbsa contributions are set to zero in phase1."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->gbsa_cutoff_distance_phase1)),
                              make_vector(std::string("gbsa-cutoff-distance-phase2"),
                                          std::string("Distance beyond which gbsa contributions are set to zero in phase2."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->gbsa_cutoff_distance_phase2))
                              )),
                    super_group, counter==1);
          }
          
          
          // GBSA implicit solvent term
          for (int counter = occurrences[prefix+"-gbsa"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermGbsa EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "gbsa: GBSA implicit solvent term (" + prefix + ")",
                         prefix+"-gbsa", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          
          
          // GBSA implicit solvent term - cached version
          for (int counter = occurrences[prefix+"-gbsa-cached"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermGbsaCached  EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
          
               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }
          
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "gbsa-cached: GBSA implicit solvent term - cached version (" + prefix + ")",
                         prefix+"-gbsa-cached", settings,
                         make_vector(
                              make_vector(std::string("maximum-deviation-cutoff"),
                                          std::string("Maximum deviation allowed in born radii in two subtrees of the chaintree before it is recalculated"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->maximum_deviation_cutoff)),
                              make_vector(std::string("cutoff-distance-phase1"),
                                          std::string("Distance beyond which contributions are set to zero in phase1."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->cutoff_distance_phase1)),
                              make_vector(std::string("cutoff-distance-phase2"),
                                          std::string("Distance beyond which contributions are set to zero in phase2."),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->cutoff_distance_phase2))
                              )),
                    super_group, counter==1);
          }
          
          
          // OPLS: automatically expands to all terms
          for (int counter = occurrences[prefix+"-opls"]; counter > 0; counter--) {
          
               // Add options
               target.add(
                    target.create_option_expansion(
                         "opls: All opls terms - incl gbsa (" + prefix + ")",
                         prefix+"-opls",
                         prefix+"-clash-fast "+prefix+"-opls-non-bonded "+prefix+"-opls-angle-bend "+prefix+"-opls-torsion "+prefix+"-opls-bond-stretch"),
                    "shorthands", counter==1);
          }
          
          // OPLS - cached: automatically expands to all terms (cached version)
          for (int counter = occurrences[prefix+"-opls"]; counter > 0; counter--) {
          
               // Add options
               target.add(
                    target.create_option_expansion(
                         "opls-cached: All opls terms - incl gbsa - cached version (" + prefix + ")",
                         prefix+"-opls-cached",
                         prefix+"-clash-fast "+prefix+"-opls-non-bonded-cached "+prefix+"-opls-angle-bend-cached "+prefix+"-opls-torsion "+prefix+"-opls-bond-stretch"),
                    "shorthands", counter==1);
          }
          
                    

     }

          
};

}
