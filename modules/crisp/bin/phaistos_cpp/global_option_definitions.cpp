namespace module_crisp {

//! Module move initialization
template <typename SETTINGS_MODIFIER>
struct MoveOptions {

     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 CHAIN_TYPE *chain,
                 DBN_TYPE *dbn) {
     }

     //! Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     MoveOptions(ProgramOptionParser &target,
                 const ProgramOptionParser::Filter &occurrences,
                 std::string super_group,
                 ChainFB *chain,
                 DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;

          // Define defaults for the different modes
          ModeDefinitions mode_definitions(target, chain);

          // CRISP - No priors
          for (int counter = occurrences["move-crisp"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveCRISP<ChainFB,
                    crisp::PreRotation<MovePrior<BondAnglePriorUninformative,
                    DihedralPriorUninformative> >,
                    crisp::PostRotation> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "CRISP (no priors)",
                         "move-crisp", settings,
                         make_vector(
                              make_vector(std::string("std-dev-bond-angle"),
                                          std::string("Standard deviation in bond angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_bond_angle)),
                              make_vector(std::string("std-dev-phi-psi"),
                                          std::string("Standard deviation in dihedral angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_phi_psi)),
                              make_vector(std::string("std-dev-omega"),
                                          std::string("Standard deviation in omega angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_omega)),
                              make_vector(std::string("only-internal-moves"),
                                          std::string("only execute internal moves"),
                                          &settings->only_internal_moves),
                              make_vector(std::string("sample-omega"),
                                          std::string("sample omega angles during prerotation"),
                                          &settings->sample_omega)
//       make_vector(std::string("prerotation-active-dofs"),
//                   std::string("How many dofs of the prerotation to activate (-1:all)"),
//                   &settings->prerotation_active_dofs)
                              )), super_group, counter==1);
          }


          // CRISP - bondangle prior: Engh-Huber
          for (int counter = occurrences["move-crisp-eh"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveCRISP<ChainFB,
                    crisp::PreRotation<MovePrior<BondAnglePriorEnghHuber,
                    DihedralPriorUninformative> >,
                    crisp::PostRotation> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "CRISP (Engh-Huber bond angle prior)",
                         "move-crisp-eh", settings,
                         make_vector(
                              make_vector(std::string("std-dev-bond-angle"),
                                          std::string("Standard deviation in bond angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_bond_angle)),
                              make_vector(std::string("std-dev-phi-psi"),
                                          std::string("Standard deviation in dihedral angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_phi_psi)),
                              make_vector(std::string("std-dev-omega"),
                                          std::string("Standard deviation in omega angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_omega)),
                              make_vector(std::string("only-internal-moves"),
                                          std::string("only execute internal moves"),
                                          &settings->only_internal_moves),
                              make_vector(std::string("sample-omega"),
                                          std::string("sample omega angles during prerotation"),
                                          &settings->sample_omega)
//       make_vector(std::string("prerotation-active-dofs"),
//                   std::string("How many dofs of the prerotation to activate (-1:all)"),
//                   &settings->prerotation_active_dofs)
                              )), super_group, counter==1);
          }



          // CRISP - dihedral prior: DBN, bondangle prior: Engh-Huber
          for (int counter = occurrences["move-crisp-dbn-eh"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveCRISP<ChainFB,
                    crisp::PreRotation<MovePriorDbn<BondAnglePriorEnghHuber,
                    DBN_TYPE> >,
                    crisp::PostRotation> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "CRISP (DBN + Engh-Huber prior)",
                         "move-crisp-dbn-eh", settings,
                         make_vector(
                              make_vector(std::string("std-dev-bond-angle"),
                                          std::string("Standard deviation in bondangle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_bond_angle)),
                              make_vector(std::string("std-dev-phi-psi"),
                                          std::string("Standard deviation in dihedral angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_phi_psi)),
                              make_vector(std::string("std-dev-omega"),
                                          std::string("Standard deviation in omega angle change (UNINITIALIZED => constraint off)"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->std_dev_omega)),
                              make_vector(std::string("dbn-consistency-window-size"),
                                          std::string("Size of window used (to each side) when bringing the dbn back to consistency. A good value for the window size is >7, and a negative window size means that the full hidden node sequence is resampled."),
                                          &settings->dbn_consistency_window_size),
                              make_vector(std::string("dbn-bias-window-size"),
                                          std::string("Size of window used when calculating bias. Approximates the move bias as P(X[i-w,j+w])/P(X'[i-w,j+w]), where w in the window size and [i,j] is the interval where angles have been changed. A good value for the window size is >7, and a negative window size means that the full bias is used."),
                                          &settings->dbn_bias_window_size),
                              make_vector(std::string("sample-omega"),
                                          std::string("sample omega angles during prerotation"),
                                          &settings->sample_omega),
                              make_vector(std::string("only-internal-moves"),
                                          std::string("only execute internal moves"),
                                          &settings->only_internal_moves)
                              )), super_group, counter==1);
          }

          // Semilocal-BGS
          for (int counter = occurrences["move-semilocal"]; counter > 0; counter--) {
               
               // Create settings object
               typedef MoveBGS<ChainFB,
                    bgs::GaussianStep<MovePrior<BondAnglePriorUninformative,
                    DihedralPriorUninformative> > > Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "semilocal",
                         "move-semilocal", settings,
                         make_vector(
                              make_vector(std::string("only-internal-moves"),
                                          std::string("only execute internal moves"),
                                          &settings->only_internal_moves),
                              make_vector(std::string("sample-omega"),
                                          std::string("sample omega angles during prerotation"),
                                          &settings->sample_omega),
                              make_vector(std::string("sample-bond-angle"),
                                          std::string("sample bond angle angles during prerotation"),
                                          &settings->sample_bond_angle),
                              make_vector(std::string("constraint-a"),
                                          std::string("value of parameter a (global constraint)"),
                                          &settings->constraint_a),
                              make_vector(std::string("constraint-b"),
                                          std::string("value of parameter a (locality constraint)"),
                                          &settings->constraint_b),
                              make_vector(std::string("omega-scaling"),
                                          std::string("scaling factor for omega angles"),
                                          &settings->omega_scaling),
                              make_vector(std::string("bond-angle-scaling"),
                                          std::string("scaling factor for bond angles"),
                                          &settings->bond_angle_scaling),
                              make_vector(std::string("skip-proline-phi"),
                                          std::string("Whether to skip prolines phi angles (modifiying the proline phi angle introduces an improper torsion change)"),
                                          &settings->skip_proline_phi)
                              )), super_group, counter==1);
          }

          for (int counter = occurrences["move-semilocal-dbn-eh"]; counter > 0; counter--) {
             

               // Create settings object
               typedef MoveBGS<ChainFB,
                    bgs::GaussianStep<MovePriorDbn<BondAnglePriorEnghHuber,
                    DBN_TYPE> > > Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "semilocal-dbn-eh",
                         "move-semilocal-dbn-eh", settings,
                         make_vector(
                              make_vector(std::string("only-internal-moves"),
                                          std::string("only execute internal moves"),
                                          &settings->only_internal_moves),
                              make_vector(std::string("sample-omega"),
                                          std::string("sample omega angles during prerotation"),
                                          &settings->sample_omega),
                              make_vector(std::string("sample-bond-angle"),
                                          std::string("sample bond angle angles during prerotation"),
                                          &settings->sample_bond_angle),
                              make_vector(std::string("constraint-a"),
                                          std::string("value of parameter a (global constraint)"),
                                          &settings->constraint_a),
                              make_vector(std::string("constraint-b"),
                                          std::string("value of parameter a (locality constraint)"),
                                          &settings->constraint_b),
                              make_vector(std::string("omega-scaling"),
                                          std::string("scaling factor for omega angles"),
                                          &settings->omega_scaling),
                              make_vector(std::string("bond-angle-scaling"),
                                          std::string("scaling factor for bond angles"),
                                          &settings->bond_angle_scaling),
                              make_vector(std::string("skip-proline-phi"),
                                          std::string("Whether to skip prolines phi angles (modifiying the proline phi angle introduces an improper torsion change)"),
                                          &settings->skip_proline_phi)
                              )), super_group, counter==1);
          }



          // CRA
          for (int counter = occurrences["move-cra"]; counter > 0; counter--) {

               // Create settings object
               typedef MoveCRA<ChainFB> Move;
               typedef typename Move::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<Move>(new Settings()));

               // Override default settings
               settings->implicit_energy = mode_definitions.implicit_energies_allowed;

               // Add options
               target.add(
                    target.create_options(
                         DefineMoveCommonOptions(),
                         "CRA",
                         "move-cra", settings,
                         make_vector(
                              make_vector(std::string("implicit-energy"),
                                          std::string("Whether the bond angle energy should be divided out (=false) or not (=true)"),
                                          &settings->implicit_energy)
                              )), super_group, counter==1);
          }

          
     }
};

}
