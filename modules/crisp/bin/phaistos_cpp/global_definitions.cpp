namespace module_crisp {

//! Module energy term initialization
struct MoveInitialization {


     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     MoveInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                        MoveCollection<CHAIN_TYPE> *move_collection,
                        std::vector<RandomNumberEngine *> *random_number_generators) {
     }

     //! Constructor - template specific case
     template <typename DBN_TYPE>
     MoveInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                        MoveCollection<ChainFB> *move_collection,
                        std::vector<RandomNumberEngine *> *random_number_generators) {

          Options::OptionValue option;


          // CRISP move (no prior)
          option = options["move-crisp"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveCRISP<ChainFB,
                                  crisp::PreRotation<MovePrior<BondAnglePriorUninformative,
                                                                   DihedralPriorUninformative> >,
                                  crisp::PostRotation>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveCRISP<ChainFB,
                                                       crisp::PreRotation<MovePrior<BondAnglePriorUninformative,
                                                                                        DihedralPriorUninformative> >,
                                                       crisp::PostRotation>(chain,options.get_settings<Settings>(option, i)));

          }

          // CRISP move (Engh-Huber prior)
          option = options["move-crisp-eh"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveCRISP<ChainFB,
                                  crisp::PreRotation<MovePrior<BondAnglePriorEnghHuber,
                                                                   DihedralPriorUninformative> >,
                                  crisp::PostRotation>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveCRISP<ChainFB,
                                                       crisp::PreRotation<MovePrior<BondAnglePriorEnghHuber,
                                                                                        DihedralPriorUninformative> >,
                                                       crisp::PostRotation>(chain,options.get_settings<Settings>(option, i)));

          }




          // CRISP move (TorusDBN + Engh-Huber priors) 
          option = options["move-crisp-dbn-eh"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename MoveCRISP<ChainFB,
                                          crisp::PreRotation<MovePriorDbn<BondAnglePriorEnghHuber,
                                                                                         DBN_TYPE> >,
                                          crisp::PostRotation>::Settings Settings;

               // Add move
               move_collection->add_move(
                    new MoveCRISP<ChainFB,
                                  crisp::PreRotation<MovePriorDbn<BondAnglePriorEnghHuber,
                                                                                 DBN_TYPE> >,
                                  crisp::PostRotation>(chain,dbn,options.get_settings<Settings>(option, i)));
          }

          // semilocal BGS move (no prior)
          option = options["move-semilocal"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveBGS<ChainFB,
                    bgs::GaussianStep<MovePrior<BondAnglePriorUninformative,
                    DihedralPriorUninformative> > > ::Settings Settings;
                                  
               // Add move
               move_collection->add_move(new MoveBGS<ChainFB,
                                         bgs::GaussianStep<MovePrior<BondAnglePriorUninformative,
                                         DihedralPriorUninformative> > > (chain,options.get_settings<Settings>(option, i)));

          }

          // semilocal BGS move (no prior)
          option = options["move-semilocal-dbn-eh"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename MoveBGS<ChainFB,
                    bgs::GaussianStep<MovePriorDbn<BondAnglePriorEnghHuber,
                    DBN_TYPE> > > ::Settings Settings;
                                  
               // Add move
               move_collection->add_move(new MoveBGS<ChainFB,
                                         bgs::GaussianStep<MovePriorDbn<BondAnglePriorEnghHuber,
                                         DBN_TYPE> > > (chain,dbn,options.get_settings<Settings>(option, i)));

          }


//          // CRISP semi-local move (no priors)
//          option = options["move-crisp-semi-local-np"];
//          for (int i=0; i<option.occurrences(); ++i) {
//
//               // Settings typedef
//               typedef MoveCRISP<ChainFB,
//                                  crisp::PreRotationSemiLocal<MovePrior<BondAnglePriorUninformative,
//                                                               DihedralPriorUninformative> >,
//                                  crisp::PostRotationSemiLocal>::Settings Settings;
//
//               // Add move
//               move_collection->add_move(new MoveCRISP<ChainFB,
//                                                       crisp::PreRotationSemiLocal<MovePrior<BondAnglePriorUninformative,
//                                                                                    DihedralPriorUninformative> >,
//                                                       crisp::PostRotationSemiLocal>(chain,options.get_settings<Settings>(option, i)),
//                                        options[option.name+"-weight"].as<double>());
//          }
//
//
//          // CRISP semi-local move (TorusDBN + Engh-Huber priors)
//          option = options["move-crisp-semi-local-dbn-eh"];
//          for (int i=0; i<option.occurrences(); ++i) {
//
//               // Settings typedef
//               typedef typename MoveCRISP<ChainFB,
//                                           crisp::PreRotationSemiLocal<MovePriorDbn<BondAnglePriorEnghHuber,
//                                                                           DBN_TYPE> >,
//                                           crisp::PostRotationSemiLocal>::Settings Settings;
//
//               // Add move
//               move_collection->add_move(new MoveCRISP<ChainFB,
//                                                       crisp::PreRotationSemiLocal<MovePriorDbn<BondAnglePriorEnghHuber,
//                                                                                       DBN_TYPE> >,
//                                                       crisp::PostRotationSemiLocal>(chain, dbn, options.get_settings<Settings>(option, i)),
//                                        options[option.name+"-weight"].as<double>());
//          }
//

          // CRA move
          option = options["move-cra"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef MoveCRA<ChainFB>::Settings Settings;

               // Add move
               move_collection->add_move(new MoveCRA<ChainFB>(chain,
                                                              options.get_settings<Settings>(option, i)));
          }

          
     }

};


}
