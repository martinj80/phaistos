namespace module_typhon {


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

          // Constrain distances Term (typhon energy function)
          option = options[prefix + "-constrain-distances"];
          for (int i = 0; i < option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermConstrainDistances<ChainFB>::Settings Settings;
               Settings settings = options.get_settings<Settings> (option, i);

               // Add energy term
               energy->add_term(new TermConstrainDistances<ChainFB> (chain, settings));
          }

          // Fix disulfide bridges
          option = options[prefix + "-constrain-disulfide-bridges"];
          for (int i = 0; i < option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermConstrainDistances<ChainFB>::Settings Settings;
               Settings ss_settings (options.get_settings<Settings> (option, i));
               ss_settings.include_bb_hbond = false;
               ss_settings.include_ca_contacts = false;
               ss_settings.include_ss_bond = true;

               // Add energy term
               energy->add_term(new TermConstrainDistances<ChainFB> (chain, ss_settings));
          }



     }

};


}
