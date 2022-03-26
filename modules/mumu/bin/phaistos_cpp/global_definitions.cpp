namespace module_mumu {

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

          // Mumu energy term
          option = options[prefix + "-mumu"];
          for (int i = 0; i < option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermMumu<ChainFB>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               // Add energy term
               energy->add_term(new TermMumu<ChainFB> (chain, settings));
          }

     }

};


}
