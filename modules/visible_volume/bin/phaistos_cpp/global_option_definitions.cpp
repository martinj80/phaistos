namespace module_visible_volume {

// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain) {                    
     }

     // Constructor - ChainFB specific case
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain) {

          // Import namespace for make_vector
          using namespace boost::fusion;

          // Basilisk explicit energy function
          for (int counter = occurrences[prefix+"-visible-volume"]; counter > 0; counter--) {

               // Create settings object
               typedef TermVisibleVolume<ChainFB> EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
               // typedef TermVisibleVolume<ChainFB>::Settings Settings;
               // boost::shared_ptr<Settings> settings(new Settings());

               std::string title_str;
               if (prefix == "energy")
                    title_str  = "Visible Volume based solvation energy";
               else
                    title_str  = "Visible Volume AA (" + prefix + ")";
               
               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         title_str,
                         prefix+"-visible-volume", settings,
                         make_vector(
                              make_vector(std::string("sphere-radius"),
                                          std::string("Max distance to neighbors considered (angstroms)"),
                                          &settings->sphere_radius),
                              make_vector(std::string("min-sphere-points"),
                                          std::string("Minimum number of discrete angle points in model"),
                                          &settings->min_sphere_points),
                              // make_vector(std::string("design-mode"),
                              //             std::string("Assume that the amino acid type is unknown during contacts center positioning"),
                              //             &settings->design_mode),
                              // make_vector(std::string("neighbor-mode"),
                              //             std::string("The representation of neighbors"),
                              //             reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<NeighborModeEnum>*>(&settings->neighbor_mode)),
                              // make_vector(std::string("neighbor-list-every"),
                              //             std::string("How often should the neighbor list be updated"),
                              //             &settings->neighbor_list_every),
                              make_vector(std::string("atom-radius"),
                                          std::string("Disk radius of atom neighbors (angstroms)"),
                                          &settings->atom_radius),
                              make_vector(std::string("volume"),
                                          std::string("Calculate visible volume (observable only)"),
                                          &settings->volume),
                              make_vector(std::string("angexp"),
                                          std::string("Calculate angular exposure (relative residue SASA analog, observable only)"),
                                          &settings->angexp),
                              make_vector(std::string("fscn"),
                                          std::string("Calculate first shell coordination number (observable only)"),
                                          &settings->fscn),
                              make_vector(std::string("fscn-threshold"),
                                          std::string("Minimum number of angle points to define contact (observable only)"),
                                          &settings->fscn_threshold)
                              )), super_group, counter==1);
          }

     }
};

}
