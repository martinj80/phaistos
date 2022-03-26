namespace module_typhon {

// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {                    
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;


           // Constrain distances
          for (int counter = occurrences[prefix+"-constrain-distances"]; counter > 0; counter--) {

               // Create settings object
               typedef TermConstrainDistances<ChainFB> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Constrain distances in a protein like HBond networks, C alpha contacts or disulfide bonds (" + prefix + ")",
                         prefix+"-constrain-distances", settings,
                         make_vector(
                              make_vector(std::string("network-filename"),
                                          std::string("Restore network from a file."),
                                          &settings->network_filename),
                              make_vector(std::string("pdb-file"),
                                          std::string("Path to PDB input file."),
                                          &settings->pdb_file),
                              make_vector(std::string("include-bb-hbond"),
                                          std::string("Include backbone-backbone hydrogen bonds in the network"),
                                          &settings->include_bb_hbond),
                              make_vector(std::string("include-sc-hbond"),
                                          std::string("Include sidechain-sidechain hydrogen bonds in the network"),
                                          &settings->include_sc_hbond),
                              make_vector(std::string("include-bb-sc-hbond"),
                                          std::string("Include backbone-sidechain hydrogen bonds in the network"),
                                          &settings->include_bb_sc_hbond),
                              make_vector(std::string("include-ca-contacts"),
                                          std::string("Include CA-CA contacts in the network"),
                                          &settings->include_ca_contacts),
                              make_vector(std::string("include-fixed-point-contacts"),
                                          std::string("Include contacts between atoms and fixed points in space"),
                                          &settings->include_fixed_point_contacts),
                              make_vector(std::string("include-ss-bond"),
                                          std::string("Include disulfide bonds (SS bonds) in the network"),
                                          &settings->include_ss_bond),
                              make_vector(std::string("init-hbond-from-native"),
                                          std::string("Whether to initialize ideal hbond distances from PDB file (rather than averages from the Top500 database)"),
                                          &settings->init_hbond_from_native),
                              make_vector(std::string("prune-network"),
                                          std::string("Whether to prune/optimize the hbond network"),
                                          &settings->ca_distance),
                              make_vector(std::string("dehydron-bb-cutoff"),
                                          std::string("Cutoff below which backbone-backbone hydrogen bonds are considered as dehydrated aka weak/broken"),
                                          &settings->dehydron_bb_cutoff),
                              make_vector(std::string("dehydron-bb-sc-cutoff"),
                                          std::string("Cutoff below which backbone-sidechain hydrogen bonds are considered as dehydrated aka weak/broken"),
                                          &settings->dehydron_bb_sc_cutoff),
                              make_vector(std::string("dehydron-sc-cutoff"),
                                          std::string("Cutoff below which sidechain-sidechain hydrogen bonds are considered as dehydrated aka weak/broken"),
                                          &settings->dehydron_sc_cutoff),
                              make_vector(std::string("ca-distance"),
                                          std::string("Cutoff distance in Angstrom to be considered a Calpha contact"),
                                          &settings->ca_distance),
                              make_vector(std::string("ss-distance"),
                                          std::string("Cutoff distance in Angstrom to be considered a disulfide contact"),
                                          &settings->ss_distance),
                              make_vector(std::string("ca-skip"),
                                          std::string("How many residues to skip along the chain before considering a Calpha contact"),
                                          &settings->ca_skip),
                              make_vector(std::string("bb-hbond-skip"),
                                          std::string("How many residues to skip along the chain before considering a backbone-backbone hydrogen bond contact"),
                                          &settings->bb_hbond_skip),
                              make_vector(std::string("sc-hbond-skip"),
                                          std::string("How many residues to skip along the chain before considering a sidechain-sidechain hydrogen bond contact"),
                                          &settings->sc_hbond_skip),
                              make_vector(std::string("bb-sc-hbond-skip"),
                                          std::string("How many residues to skip along the chain before considering a backbone-sidechain hydrogen bond contact"),
                                          &settings->bb_sc_hbond_skip),
                              make_vector(std::string("ss-skip"),
                                          std::string("How many residues to skip along the chain before considering a disulfide contact"),
                                          &settings->ss_skip),
                              make_vector(std::string("generate-pymol"),
                                          std::string("Whether to generate a Python script that will visualize the network in PyMOL"),
                                          &settings->generate_pymol),
                              make_vector(std::string("use-caching (recommended)"),
                                          std::string("Whether to cache interations"),
                                          &settings->use_caching),
                              make_vector(std::string("verbose"),
                                          std::string("Whether to print out additional information (recommended)"),
                                          &settings->verbose)
                              )), super_group, counter==1);
          }

           // fix ss bonds in proteins
          for (int counter = occurrences[prefix+"-constrain-disulfide-bonds"]; counter > 0; counter--) {

               // Create settings object
               typedef TermConstrainDistances<ChainFB> EnergyTerm;
               typedef typename EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Lock disulfide bonds found in proteins in place (" + prefix + ")",
                         prefix+"-constrain-disulfide-bonds", settings,
                         make_vector(
                              make_vector(std::string("network-filename"),
                                          std::string("Restore network from a file."),
                                          &settings->network_filename),
                              make_vector(std::string("include-ss-bond"),
                                          std::string("Include disulfide bonds (SS bonds) in the network"),
                                          &settings->include_ss_bond),
                              make_vector(std::string("ss-distance"),
                                          std::string("Default CYS(SG)-CYS(SG) cutoff distance."),
                                          &settings->ss_distance),
                              make_vector(std::string("ss-skip"),
                                          std::string("Default minimum CYS-CYS distance in the chain."),
                                          &settings->ss_skip),
                              make_vector(std::string("generate-pymol"),
                                          std::string("Whether to generate a Python script that will visualize the network in PyMOL"),
                                          &settings->generate_pymol),
                              make_vector(std::string("use-caching (recommended)"),
                                          std::string("Whether to cache interations"),
                                          &settings->use_caching),
                              make_vector(std::string("verbose"),
                                          std::string("Whether to print out additional information (recommended)"),
                                          &settings->verbose)
                                   )), super_group, counter==1);
          }

     }
};

}
