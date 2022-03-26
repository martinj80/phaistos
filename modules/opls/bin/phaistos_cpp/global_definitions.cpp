namespace module_opls {

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

          // OPLS - charge term
          option = options[prefix+"-opls-charge"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsCharge::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsCharge(chain,
                                                   options.get_settings<Settings>(option,i)));
          }


          // OPLS - charge term - cached
          option = options[prefix+"-opls-charge-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsChargeCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsChargeCached(chain,
                                                         options.get_settings<Settings>(option,i)));
          }


          // OPLS - vdw term
          option = options[prefix+"-opls-vdw"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsVdw::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsVdw(chain,
                                                options.get_settings<Settings>(option,i)));
          }


          // OPLS - vdw term - cached
          option = options[prefix+"-opls-vdw-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsVdwCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsVdwCached(chain,
                                                      options.get_settings<Settings>(option,i)));
          }


          // OPLS - angle-bend term
          option = options[prefix+"-opls-angle-bend"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsAngleBend::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsAngleBend(chain,
                                                      options.get_settings<Settings>(option,i)));
          }


          // OPLS - angle-bend term - cached
          option = options[prefix+"-opls-angle-bend-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsAngleBendCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsAngleBendCached(chain,
                                                            options.get_settings<Settings>(option,i)));
          }


          // OPLS - torsion term
          option = options[prefix+"-opls-torsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsTorsion::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsTorsion(chain,
                                                    options.get_settings<Settings>(option,i)));
          }


          // OPLS - improper-torsion term
          option = options[prefix+"-opls-improper-torsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsImptor::Settings Settings;
                    
               // Add energy term
               energy->add_term(new TermOplsImptor(chain,
                                                   options.get_settings<Settings>(option,i)));
          }


          // OPLS - bondstretch term
          option = options[prefix+"-opls-bond-stretch"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsBondStretch::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsBondStretch(chain,
                                                        options.get_settings<Settings>(option,i)));
          }


          // OPLS - nonbonded terms
          option = options[prefix+"-opls-non-bonded"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsNonBonded<double>::Settings Settings;

               // Add energy term
               energy->add_term(new TermOplsNonBonded<double>(chain,
                                                              options.get_settings<Settings>(option,i)));
          }


          // OPLS - nonbonded terms (cached version)
          option = options[prefix+"-opls-non-bonded-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermOplsNonBondedCached<double>::Settings Settings;
               
               // Add energy term
               energy->add_term(new TermOplsNonBondedCached<double>(chain,
                                                                    options.get_settings<Settings>(option,i)));
          }
          
          
          // GBSA solvent term
          option = options[prefix+"-gbsa"];
          for (int i=0; i<option.occurrences(); ++i) {
               
               // Settings typedef
               typedef typename TermGbsa::Settings Settings;
               
               // Add energy term
               energy->add_term(new TermGbsa(chain,
                                             options.get_settings<Settings>(option,i)));
          }


          // GBSA solvent term (cached version)
          option = options[prefix+"-gbsa-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermGbsaCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermGbsaCached(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
     }

};

}
