namespace module_profasi {

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

          // Profasi: Local term
          for (int counter = occurrences[prefix+"-profasi-local"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiLocal EnergyTerm;
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
                         "profasi-local: Profasi local energy term (" + prefix + ")",
                         prefix+"-profasi-local", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Local term - cached
          for (int counter = occurrences[prefix+"-profasi-local-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiLocalCached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
                    //std::cout<<"dtu-setting->weight-profasi-local"<< settings->weight;
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "profasi-local-cached: Profasi local energy term - cached version(" + prefix + ")",
                         prefix+"-profasi-local-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }


          // Profasi: Local sidechain term
          for (int counter = occurrences[prefix+"-profasi-local-sidechain"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiLocalSidechain EnergyTerm;
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
                         "profasi-local-sidechain: Profasi local sidechain energy term (" + prefix + ")",
                         prefix+"-profasi-local-sidechain", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Local sidechain term - cached
          for (int counter = occurrences[prefix+"-profasi-local-sidechain-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiLocalSidechainCached EnergyTerm;
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
                         "profasi-local-sidechain-cached: Profasi local sidechain energy term - cached version(" + prefix + ")",
                         prefix+"-profasi-local-sidechain-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Excluded volume
          for (int counter = occurrences[prefix+"-profasi-excluded-volume"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiExcludedVolume EnergyTerm;
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
                         "profasi-excluded-volume: Profasi excluded volume energy term (" + prefix + ")",
                         prefix+"-profasi-excluded-volume", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Excluded volume - cached
          for (int counter = occurrences[prefix+"-profasi-excluded-volume-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiExcludedVolumeCached EnergyTerm;
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
                         "profasi-excluded-volume-cached: Profasi excluded volume energy term - cached version(" + prefix + ")",
                         prefix+"-profasi-excluded-volume-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }


          // Profasi: Excluded volume (local)
          for (int counter = occurrences[prefix+"-profasi-excluded-volume-local"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiExcludedVolumeLocal EnergyTerm;
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
                         "profasi-excluded-volume-local: Profasi local excluded volume energy term (" + prefix + ")",
                         prefix+"-profasi-excluded-volume-local", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Excluded volume (local) - cached
          for (int counter = occurrences[prefix+"-profasi-excluded-volume-local-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiExcludedVolumeLocalCached EnergyTerm;
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
                         "profasi-excluded-volume-local-cached: Profasi local excluded volume energy term - cached version(" + prefix + ")",
                         prefix+"-profasi-excluded-volume-local-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }


          // Profasi: Hydrogen bond
          for (int counter = occurrences[prefix+"-profasi-hydrogen-bond"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiHydrogenBond EnergyTerm;
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
                         "profasi-hydrogen-bond: Profasi hydrogen bond term (" + prefix + ")",
                         prefix+"-profasi-hydrogen-bond", settings,
                         make_vector(
			      make_vector(std::string("use-ideal-distances"),
					  std::string("Whether to use ideal distance for C-O and H-N."),
					  &settings->use_ideal_distances)
			      )),
                    super_group, counter==1);
          }

          // Profasi: Hydrogen bond - cached
          for (int counter = occurrences[prefix+"-profasi-hydrogen-bond-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiHydrogenBondCached EnergyTerm;
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
                         "profasi-hydrogen-bond-cached: Profasi hydrogen bond term - cached version(" + prefix + ")",
                         prefix+"-profasi-hydrogen-bond-cached", settings,
                         make_vector(
                    make_vector(std::string("use-ideal-distances"),
                          std::string("Whether to use ideal distance for C-O and H-N."),
                          &settings->use_ideal_distances)
                     )),
                    super_group, counter==1);
          }

          // Profasi: Hydrogen bond - improved
          for (int counter = occurrences[prefix+"-profasi-hydrogen-bond-improved"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiHydrogenBondImproved EnergyTerm;
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
                         "profasi-hydrogen-bond-improved: Profasi hydrogen bond term - improved version(" + prefix + ")",
                         prefix+"-profasi-hydrogen-bond-improved", settings,
                         make_vector(
                    make_vector(std::string("use-ideal-distances"),
                          std::string("Whether to use ideal distance for C-O and H-N."),
                          &settings->use_ideal_distances)
                     )),
                    super_group, counter==1);
          }


          // Profasi: Hydrophobicity
          for (int counter = occurrences[prefix+"-profasi-hydrophobicity"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiHydrophobicity EnergyTerm;
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
                         "profasi-hydrophobicity: Profasi hydrophobicity term (" + prefix + ")",
                         prefix+"-profasi-hydrophobicity", settings,
                         make_vector()),
                    super_group, counter==1);
          }


          // Profasi: Hydrophobicity cached
          for (int counter = occurrences[prefix+"-profasi-hydrophobicity-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiHydrophobicityCached EnergyTerm;
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
                         "profasi-hydrophobicity: Profasi hydrophobicity term - cached version(" + prefix + ")",
                         prefix+"-profasi-hydrophobicity-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Hydrophobicity improved
          for (int counter = occurrences[prefix+"-profasi-hydrophobicity-improved"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiHydrophobicityImproved EnergyTerm;
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
                         "profasi-hydrophobicity: Profasi hydrophobicity term - improved version(" + prefix + ")",
                         prefix+"-profasi-hydrophobicity-improved", settings,
                         make_vector()),
                    super_group, counter==1);
          }


          // Profasi: Sidechain charge
          for (int counter = occurrences[prefix+"-profasi-sidechain-charge"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiSidechainCharge EnergyTerm;
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
                         "profasi-sidechain-charge: Profasi sidechain charge term (" + prefix + ")",
                         prefix+"-profasi-sidechain-charge", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Profasi: Sidechain charge - cached
          for (int counter = occurrences[prefix+"-profasi-sidechain-charge-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiSidechainChargeCached EnergyTerm;
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
                         "profasi-sidechain-charge: Profasi sidechain charge term - cached version(" + prefix + ")",
                         prefix+"-profasi-sidechain-charge-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }


          // Profasi: Sidechain charge - improved
          for (int counter = occurrences[prefix+"-profasi-sidechain-charge-improved"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiSidechainChargeImproved EnergyTerm;
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
                         "profasi-sidechain-charge: Profasi sidechain charge term - improved version(" + prefix + ")",
                         prefix+"-profasi-sidechain-charge-improved", settings,
                         make_vector()),
                    super_group, counter==1);
          }




          // Profasi: proline phi torsion term
          for (int counter = occurrences[prefix+"-profasi-proline-phi-torsion"]; counter > 0; counter--) {

               // Create settings object
               typedef TermProfasiProlinePhiTorsion EnergyTerm;
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
                         "profasi-proline-phi-torsion: Profasi proline phi torsion term(" + prefix + ")",
                         prefix+"-profasi-proline-phi-torsion", settings,
                         make_vector()),
                    super_group, counter==1);
          }



          // Profasi: automatically expands to all terms
          for (int counter = occurrences[prefix+"-profasi"]; counter > 0; counter--) {

               // Add options
               target.add(
                    target.create_option_expansion(
                         "profasi: All profasi terms (" + prefix + ")",
                         prefix+"-profasi",
                         // prefix+"-clash-fast "+prefix+"-profasi-local "+prefix+"-profasi-local-sidechain "+prefix+"-profasi-excluded-volume "+prefix+"-profasi-excluded-volume-local "+prefix+"-profasi-hydrogen-bond "+prefix+"-profasi-hydrophobicity "+prefix+"-profasi-sidechain-charge"),
                         prefix+"-profasi-local "+prefix+"-profasi-local-sidechain "+prefix+"-profasi-excluded-volume "+prefix+"-profasi-excluded-volume-local "+prefix+"-profasi-hydrogen-bond "+prefix+"-profasi-hydrophobicity "+prefix+"-profasi-sidechain-charge"),
                    "shorthands", counter==1);
          }

          // Profasi-cached version: automatically expands to all cached terms
          for (int counter = occurrences[prefix+"-profasi-cached"]; counter > 0; counter--) {

               // Add options
               target.add(
                    target.create_option_expansion(
                         "profasi: All profasi terms - cached version(" + prefix + ")",
                         prefix+"-profasi-cached",
                         //prefix+"-clash-fast "+prefix+"-profasi-local-cached "+prefix+"-profasi-local-sidechain-cached "+prefix+"-profasi-excluded-volume-cached "+prefix+"-profasi-excluded-volume-local-cached "+prefix+"-profasi-hydrogen-bond-cached "+prefix+"-profasi-hydrophobicity-cached "+prefix+"-profasi-sidechain-charge-cached"),
                         prefix+"-profasi-local-cached "+prefix+"-profasi-local-sidechain-cached "+prefix+"-profasi-excluded-volume-cached "+prefix+"-profasi-excluded-volume-local-cached "+prefix+"-profasi-hydrogen-bond-improved "+prefix+"-profasi-hydrophobicity-cached "+prefix+"-profasi-sidechain-charge-cached"),
                    "shorthands", counter==1);
          }


     }

};

}
