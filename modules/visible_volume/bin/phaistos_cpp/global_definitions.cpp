namespace module_visible_volume {

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain,
                          Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                          std::string prefix="") {
     }

     // Constructor - template specific case
     EnergyInitialization(const Options &options, ChainFB *chain,
                          Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                          std::string prefix="") {

          Options::OptionValue option;

          // Basilisk explicit energy term
          option = options[prefix + "-visible-volume"];
          for (int i = 0; i < option.occurrences(); ++i) {

               // Settings typedef
               typedef TermVisibleVolume<ChainFB>::Settings Settings;
               Settings settings = options.get_settings<Settings>(option, i);

               // Add energy term
               energy->add_term(new TermVisibleVolume<ChainFB> (chain, settings));
          }

     }

};


}
