namespace module_git {

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain,
                          Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                          std::string prefix="") {
          Options::OptionValue option;

          // Trajectory observable
          option = options[prefix+"-git"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef (note that this is an observable)
               typedef typename Observable<TermGit<CHAIN_TYPE> >::Settings Settings;

               // Add energy term (note that this is an observable)
               energy->add_term(new Observable<TermGit<CHAIN_TYPE> >(chain,
                                                                     options.get_settings<Settings>(option,i)));
          }
     }

};

}
