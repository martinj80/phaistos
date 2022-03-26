namespace module_profasi {

//! Module: energy term initialization
struct EnergyInitialization {


     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
			  Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators,
			  std::string prefix="") {
     }

     //! Constructor - template specific case
     template <typename DBN_TYPE>
     EnergyInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
			  Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators,
			  std::string prefix="") {

          Options::OptionValue option;

          // Profasi - local term
          option = options[prefix+"-profasi-local"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiLocal::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiLocal(chain,
                                                     options.get_settings<Settings>(option,i)));
          }

          // Profasi - local term - cached version
          option = options[prefix+"-profasi-local-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiLocalCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiLocalCached(chain,
                                        options.get_settings<Settings>(option,i)));
          }


          // Profasi - local sidechain term
          option = options[prefix+"-profasi-local-sidechain"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiLocalSidechain::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiLocalSidechain(chain,
                                                              options.get_settings<Settings>(option,i)));
          }

          // Profasi - local sidechain term - cached version
          option = options[prefix+"-profasi-local-sidechain-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiLocalSidechainCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiLocalSidechainCached(chain,
                                                                    options.get_settings<Settings>(option,i)));
          }


          // Profasi - excluded volume term
          option = options[prefix+"-profasi-excluded-volume"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiExcludedVolume::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiExcludedVolume(chain,
                                                              options.get_settings<Settings>(option,i)));
          }

          // Profasi - excluded volume term - cached version
          option = options[prefix+"-profasi-excluded-volume-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiExcludedVolumeCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiExcludedVolumeCached(chain,
                                                                    options.get_settings<Settings>(option,i)));
          }

          // Profasi - local excluded volume term
          option = options[prefix+"-profasi-excluded-volume-local"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiExcludedVolumeLocal::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiExcludedVolumeLocal(chain,
                                                                   options.get_settings<Settings>(option,i)));
          }

          // Profasi - local excluded volume term - cached version
          option = options[prefix+"-profasi-excluded-volume-local-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiExcludedVolumeLocalCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiExcludedVolumeLocalCached(
                                     chain,
                                     options.get_settings<Settings>(option,i)));
          }

          // Profasi - hydrogen bond term
          option = options[prefix+"-profasi-hydrogen-bond"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiHydrogenBond::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiHydrogenBond(chain,
                                                            options.get_settings<Settings>(option,i)));
          }

          // Profasi - hydrogen bond term - cached version
          option = options[prefix+"-profasi-hydrogen-bond-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiHydrogenBondCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiHydrogenBondCached(chain,
                                                                  options.get_settings<Settings>(option,i)));
          }

          // Profasi - hydrogen bond term - improved version
          option = options[prefix+"-profasi-hydrogen-bond-improved"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiHydrogenBondImproved::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiHydrogenBondImproved(chain,
                                                                    options.get_settings<Settings>(option,i)));
          }


          // Profasi - hydrophobicity term
          option = options[prefix+"-profasi-hydrophobicity"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiHydrophobicity::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiHydrophobicity(chain,
                                                              options.get_settings<Settings>(option,i)));
          }

          // Profasi - hydrophobicity term - cached version
          option = options[prefix+"-profasi-hydrophobicity-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiHydrophobicityCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiHydrophobicityCached(chain,
                                                                    options.get_settings<Settings>(option,i)));
          }

          // Profasi - hydrophobicity term - cached version
          option = options[prefix+"-profasi-hydrophobicity-improved"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiHydrophobicityImproved::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiHydrophobicityImproved(chain,
                                                                      options.get_settings<Settings>(option,i)));
          }


          // Profasi - sidechain charge term
          option = options[prefix+"-profasi-sidechain-charge"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiSidechainCharge::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiSidechainCharge(chain,
                                                               options.get_settings<Settings>(option,i)));
          }

          // Profasi - sidechain charge term - cached version
          option = options[prefix+"-profasi-sidechain-charge-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiSidechainChargeCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiSidechainChargeCached(chain,
                                                                     options.get_settings<Settings>(option,i)));
          }

          // Profasi - sidechain charge term - improved version
          option = options[prefix+"-profasi-sidechain-charge-improved"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiSidechainChargeImproved::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiSidechainChargeImproved(chain,
                                                                       options.get_settings<Settings>(option,i)));
          }


          // Profasi - proline phi torsion term 
          option = options[prefix+"-profasi-proline-phi-torsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermProfasiProlinePhiTorsion::Settings Settings;

               // Add energy term
               energy->add_term(new TermProfasiProlinePhiTorsion(chain,
                                                                 options.get_settings<Settings>(option,i)));
          }

     }

};

}
