#ifdef HAVE_MUNINNLIB

                    // Muninn options
                    if (occurrences["monte-carlo-muninn"]) {

                         // Create settings object
                        typedef typename MonteCarloMuninn<CHAIN_TYPE>::Settings Settings;
                        boost::shared_ptr<Settings> settings(new Settings());

                        // // If enums are used within a settings object, we wrap them in an WrappedEnumPointer
                        // ProgramOptionParser::WrappedEnumPointer<GeEnum> *ge_enum_pointer = target.generate_wrapped_enum_pointer(&(settings->ensemble_type));

                        // Add options
                        target.add(
                             target.create_options(
                                  DefineMonteCarloCommonOptions(),
                                  "Muninn options",
                                  "monte-carlo-muninn", settings,
                                  make_vector(
                                       make_vector(std::string("energy-min"),
                                                   std::string("Lower bound on energy"),
                                                   reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->energy_min)),
                                       make_vector(std::string("energy-max"),
                                                   std::string("Upper bound on energy"),
                                                   reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->energy_max)),
                                       make_vector(std::string("histograms-per-reinit"),
                                                   std::string("Number of histogram updates between every reinitialization (zero or negative means no reinitialization will be done)"),
                                                   &settings->histograms_per_reinit),
                                       make_vector(std::string("burnin-fraction"),
                                                   std::string("The fraction of initmax used for burn-in (in each thread)"),
                                                   &settings->burnin_fraction),
                                       make_vector(std::string("use-energy2"),
                                                   std::string("Use the secondary energy as an additional energy in muninn - not used to estimate histograms"),
                                                   &settings->use_energy_secondary),
                                       make_vector(std::string("weight-scheme"),
                                                   std::string("Weight-scheme to use: invk|multicanonical|invkp"),
                                                   reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<Muninn::GeEnum>*>(&settings->weight_scheme)),
                                                   // ge_enum_pointer),
                                       make_vector(std::string("slope-factor-up"),
                                                   std::string("Slope factor use for the linear extrapolation of the weights, when the weights are increasing in the direction away from the main area of support"),
                                                   &settings->slope_factor_up),
                                       make_vector(std::string("slope-factor-down"),
                                                   std::string("Slope factor use for the linear extrapolation of the weights, when the weights are decreasing in the direction away from the main area of support"),
                                                   &settings->slope_factor_down),
                                       make_vector(std::string("min-beta"),
                                                   std::string("The minimal beta value to be used based on thermodynamics and in the extrapolation"),
                                                   reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->min_beta)),
                                       make_vector(std::string("max-beta"),
                                                   std::string("The maximal beta value to be used based on thermodynamics"),
                                                   reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->max_beta)),
                                       make_vector(std::string("initial-beta"),
                                                   std::string("{The initial beta value to be used (if uninitialized it takes the same value as min-beta)."),
                                                   reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->initial_beta)),
                                       make_vector(std::string("p"),
                                                   std::string("The p value for the invkp scheme, which should be between 0 and 1. If p=0 the invkp scheme resembles the multicanonical. If p=1 the invkp scheme resembles the invk."),
                                                   &settings->p),
                                       make_vector(std::string("resolution"),
                                                   std::string("The resolution for the non-uniform binner"),
                                                   &settings->resolution),
                                       make_vector(std::string("initial-width-is-max-left"),
                                                   std::string("Use the initial bin with as maximal bin width, when expanding to the left"),
                                                   &settings->initial_width_is_max_left),
                                       make_vector(std::string("initial-width-is-max-right"),
                                                   std::string("Use the initial bin with as maximal bin width, when expanding to the right"),
                                                   &settings->initial_width_is_max_right),
                                       make_vector(std::string("statistics-log-filename"),
                                                   std::string("The filename (including full path) for the Muninn statistics logfile (default is Muninn_[PID].txt and is obtained by removing this option from the config file or setting it to \"\")"),
                                                   &settings->statistics_log_filename),
                                       make_vector(std::string("log-mode"),
                                                   std::string("Muninn log mode (current|all)"),
                                                   reinterpret_cast<ProgramOptionParser::WrappedEnumPointer<Muninn::StatisticsLogger::Mode>*>(&settings->log_mode)),
                                       make_vector(std::string("read-statistics-log-filename"),
                                                   std::string("The filename for reading the Muninn statistics logfile. If the value difference from the empty string (\"\"), this log file is read and the history is set based on the content."),
                                                   &settings->read_statistics_log_filename),
                                       make_vector(std::string("read-fixed-weights-filename"),
                                                   std::string("The filename for reading a set of fixed weights for a given region. If the value difference from the empty string (\"\"), this file is read and the weights are fixed in the given region."),
                                                   &settings->read_fixed_weights_filename),
                                       make_vector(std::string("initial-max"),
                                                   std::string("Number of iterations used in first round of sampling"),
                                                   &settings->initial_max),
                                       make_vector(std::string("increase-factor"),
                                                   std::string("Scaling of iterations used in the subsequent round of sampling"),
                                                   &settings->increase_factor),
                                       make_vector(std::string("max-iterations-between-rounds"),
                                                   std::string("The maximum number of iterations between consecutive histograms"),
                                                   &settings->max_iterations_per_histogram),
                                       make_vector(std::string("memory"),
                                                   std::string("The number of consecutive histograms to keep in memory"),
                                                   &settings->memory),
                                       make_vector(std::string("min-count"),
                                                   std::string("The minimal number of counts in a bin in order to have support in that bin"),
                                                   &settings->min_count),
                                       make_vector(std::string("restricted-individual-support"),
                                                   std::string("Restrict the support of the individual histograms to only cover the support for the given histogram"),
                                                   &settings->restricted_individual_support),
                                       make_vector(std::string("dynamic-binning"),
                                                   std::string("Use dynamic binning"),
                                                   &settings->use_dynamic_binning),
                                       make_vector(std::string("max-number-of-bins"),
                                                   std::string("The maximal number of bins allowed to be used by the binner"),
                                                   &settings->max_number_of_bins),
                                       make_vector(std::string("bin-width"),
                                                   std::string("Bin width used for non-dynamic binning"),
                                                   &settings->bin_width),
                                       make_vector(std::string("verbose"),
                                                   std::string("The verbose level of Muninn"),
                                                   &settings->verbose)
                                       )), super_group);
                    }
                    mc_mode_description.push_back("muninn");

#endif
