namespace module_git {

//! Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Trajectory observable
          for (int counter = occurrences[prefix+"-git"]; counter > 0; counter--) {
          
               // Create settings object
               typedef Observable<TermGit<CHAIN_TYPE> > EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Since we are adding a term that can only be applied as observable, make sure 
               // that we do not add it to anything but the observable collection
               if (prefix!="observable")
                    continue;

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "git: Dump git vectors (" + prefix + ")",
                         prefix+"-git", settings,
                         make_vector()),
                    super_group, counter==1);
          }          
     }
};

}
