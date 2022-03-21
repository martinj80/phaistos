#ifdef HAVE_MUNINNLIB
          } else if ((option = options["monte-carlo-muninn"]).occurrences()) {

               // Parse the options and set the MonteCarloMuninn settings object
               typedef typename MonteCarloMuninn<CHAIN_TYPE>::Settings Settings;

               // Setup the MonteCarloMuninn class
               MonteCarloMuninn<CHAIN_TYPE> *monte_carlo_tmp =
                    new MonteCarloMuninn<CHAIN_TYPE>(move_collection->chain,
                                                     energy,
                                                     move_collection,
                                                     options.get_settings<Settings>(option),
                                                     energy_secondary
                         );

               if (options["threads"].as<int>() > 1) {
                    *monte_carlo =
                         new MonteCarloMultiThread<MonteCarloMuninn<CHAIN_TYPE> >(monte_carlo_tmp,
                                                                                  options["threads"].as<int>(),
                                                                                  options["steps-per-move"].as<int>(),
                                                                                  *random_number_engines);
               } else {
                    *monte_carlo = monte_carlo_tmp;
               }
#endif
