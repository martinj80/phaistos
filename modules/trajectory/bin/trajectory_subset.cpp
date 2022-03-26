// trajectory_subset.cpp --- Extract a subset from a trajectory
// Copyright (C) 2008-2013 Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#include <numeric>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "utils/random.h"

#include "protein/chain_ca.h"
#include "protein/chain_fb.h"

#include "protein/xtc_chain.h"

#ifdef HAVE_MUNINNLIB     
#include "muninn/tools/CanonicalAveragerFromStatisticsLog.h"
#endif

// Import Phaistos namespace
using namespace phaistos;


//! Templated main function
template<typename CHAIN_TYPE>
void trajectory_subset_main(const po::variables_map &options) {

     // Check for pdb filename 
     std::string pdb_filename = "";
     if (options.count("pdb-file")) {
          pdb_filename = options["pdb-file"].as<std::string>();
     } else {
          std::cerr << "Please specify PDB file name (necessary for decoding trajectory)\n";
          exit(1);
     }

     // Initialize chain
     CHAIN_TYPE *chain = new CHAIN_TYPE(pdb_filename,
                                        definitions::ALL_PHYSICAL_ATOMS);
     chain->add_atoms(definitions::ALL_PHYSICAL_ATOMS);


     // Read observables matrix
     std::vector<std::vector<std::string> > observables;
     if (options.count("observable-files")) {
          std::vector<std::string> observable_filenames = options["observable-files"].as<std::vector<std::string> >();

          for (unsigned int i=0; i<observable_filenames.size(); ++i) {
               std::string observable_filename = observable_filenames[i];
               std::vector<std::string> observable_line_vector = file_to_string_vector(observable_filename);
               for (unsigned int j=0; j<observable_line_vector.size(); ++j) {
                    if (observable_line_vector[j] != "" && observable_line_vector[j][0] != '#' ) {
                         std::vector<std::string> split_line;
                         boost::split(split_line, observable_line_vector[j], boost::is_any_of(" \t"));
                         observables.push_back(split_line);
                    }
               }
          }
     }


     // Read trajectory filenames
     std::vector<std::string> trajectory_filenames;
     if (options.count("trajectory-files")) {
          trajectory_filenames = options["trajectory-files"].as<std::vector<std::string> >();
     }


     // Read iteration range
     PHAISTOS_LONG_LONG iteration_start = options["iteration-start"].as<PHAISTOS_LONG_LONG>();
     PHAISTOS_LONG_LONG iteration_end = options["iteration-end"].as<PHAISTOS_LONG_LONG>();


     // Container for all frames in range -- (trajectory_file_index, frame_index) pairs 
     std::vector<std::pair<int,int> > frames;
     // Map between index in frames index and corresponding entry in observables file
     std::vector<int> frame_observable_index;


     // Variables used by XTC reader/writer
     int n_atoms = 0;
     rvec *coordinates = NULL;
     bool dry_run=true;

     // Iterate over trajectory files, and trajectory frames. Register
     // all frames within iteration range
     unsigned int frame_counter = 0;
     for (unsigned int i=0; i<trajectory_filenames.size(); ++i) {
          std::string trajectory_filename = trajectory_filenames[i];

          if (trajectory_filename != "") {
               // Check consistency between chain and trajectory
               xtc_chain_check_consistency(*chain, trajectory_filename);

               int step;
               float time;
               float prec;
               XDRFILE *trajectory_file = xdrfile_open(trajectory_filename.c_str(),"r");
               while (xtc_read_chain(*chain, trajectory_file, step, time, prec, &n_atoms, &coordinates, dry_run)) {
                    if (step >= iteration_start && (iteration_end < 0 || step < iteration_end)) {
                         frames.push_back(std::make_pair(i, step));
                         frame_observable_index.push_back(frame_counter);
                    }
                    frame_counter++;
               }
          }
     }

     // Check whether number of frames equals number of observables
     if (options.count("observable-files")) {
          if (frame_counter != observables.size()) {
               std::cerr << "Mismatch between framecounts in trajectory (" << frame_counter << ") and length of observable vector (" << observables.size() << ")\n";
               exit(1);
          }
     }

#ifdef HAVE_MUNINNLIB     

     // Column index specifying which column in observable file to use as muninn reaction coordinate
     int muninn_column_index = 0;

     // Check if muninn log file has been specified
     std::string muninn_log_filename = "";
     if (options.count("muninn-log-file")) {
          muninn_log_filename = options["muninn-log-file"].as<std::string>();

          if (options.count("muninn-column")) {
               muninn_column_index = options["muninn-column"].as<int>();
          } else {
               std::cerr << "Please use muninn-column option to specify which column in observable contains the muninn observable.\n";
               exit(1);
          }
     }
     
     std::vector<double> muninn_observables;

     std::vector<double> frame_weights;
     if (muninn_log_filename != "") {

          for (unsigned int i=0; i<frames.size(); ++i) {
          
               int observable_index = frame_observable_index[i];

               double muninn_observable = boost::lexical_cast<double>(observables[observable_index][muninn_column_index]);
               muninn_observables.push_back(muninn_observable);
          }

          Muninn::CanonicalAveragerFromStatisticsLog canonical_averager(muninn_log_filename);
          frame_weights = canonical_averager.calc_weights(muninn_observables);

     } else {

#endif

          if (options.count("observable-files")) {
               
               // Column index specifying which column in observable file to use as weights
               int weight_column_index = options["weight-column"].as<int>();
               
               // Read in weights from observable file
               for (unsigned int i=0; i<frames.size(); ++i) {
                    
                    int observable_index = frame_observable_index[i];
                    
                    double weight = boost::lexical_cast<double>(observables[observable_index][weight_column_index]);
                    frame_weights.push_back(weight);
               }
          }


#ifdef HAVE_MUNINNLIB     
     }
#endif

     // Indices of randomly selected frames
     std::vector<int> sampled_frames;

     if (frame_weights.size() > 0) {

          // Calculate cumulative weights
          std::vector<double> cumulative_weights;
          std::partial_sum(frame_weights.begin(), frame_weights.end(),
                           std::back_inserter(cumulative_weights));

          // Construct random number generator
          boost::uniform_real<> dist(0, cumulative_weights.back());
          boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random_number_generator(random_global, dist);

          // Generate requested number of samples
          unsigned int samples = options["samples"].as<PHAISTOS_LONG_LONG>();
          for (unsigned int i=0; i<samples; ++i) {
               sampled_frames.push_back((std::lower_bound(cumulative_weights.begin(), cumulative_weights.end(), 
                                                          random_number_generator()) - cumulative_weights.begin()) + 1);
          }

          // Sort samples - single indices and (trajectory_file_index, frame_index) have same order
          std::sort(sampled_frames.begin(), sampled_frames.end());
     } else {
          
          // Copy all frame indices to sampled_frames
          for (unsigned int i=0; i<frames.size(); ++i) {
               sampled_frames.push_back(i);
          }

     }

     // Extract output filename
     std::string output_filename = options["output-file"].as<std::string>();

     // Output as pdb if extension is .pdb
     bool output_pdb = (output_filename.substr(output_filename.length()-3)=="pdb");

     // Open either PDB or XTC file
     std::ofstream pdb_file;
     XDRFILE *xd;
     if (output_pdb) {
          pdb_file.open(output_filename.c_str());
     } else {
          xd = xdrfile_open(output_filename.c_str(),"w");
     }

     // Iterate over all frames - selecting only those that match previously sampled indices
     unsigned int sample_index = 0;
     unsigned int trajectory_file_index = frames[sampled_frames[sample_index]].first;
     unsigned int frame_index = frames[sampled_frames[sample_index]].second;
     for (unsigned int i=0; i<trajectory_filenames.size() && sample_index < sampled_frames.size(); ++i) {
          std::string trajectory_filename = trajectory_filenames[i];

          if (trajectory_filename != "") {

               int step;
               float time;
               float prec;
               XDRFILE *trajectory_file = xdrfile_open(trajectory_filename.c_str(),"r");
               while (xtc_read_chain(*chain, trajectory_file, step, time, prec, &n_atoms, &coordinates, dry_run)) {

                    // Inner loop incase the same index appears multiple times
                    while (i==trajectory_file_index && step==(int)frame_index) {

                         xtc_coordinates_to_chain(*chain, coordinates);

                         // Output either in PDB or XTC format
                         if (output_pdb) {
                              CHAIN_TYPE *superimpose_chain = NULL;
                              int begin_offset = 0;
                              int end_offset = 0;

                              pdb_file << chain->output_as_pdb(superimpose_chain, 
                                                               begin_offset, 
                                                               end_offset, 
                                                               sample_index+1);
                              pdb_file << "END\nENDMDL\n";
                         } else {
                              xtc_write_chain(*chain, xd, step, time, prec);
                         }

                         // Output message to screen
                         if (observables.size() > 0)
                              std::cout << "Dumped: " << observables[frame_observable_index[sampled_frames[sample_index]]] << "\n";
                         else
                              std::cout << "Dumped: " << trajectory_file_index << " " << frame_index << "\n";
                              

                         // Increment sample counter
                         ++sample_index;

                         if (sample_index >= sampled_frames.size())
                              break;

                         trajectory_file_index = frames[sampled_frames[sample_index]].first;
                         frame_index = frames[sampled_frames[sample_index]].second;
                    }

                    if (sample_index >= sampled_frames.size())
                         break;

               }
          }
     }

     if (output_pdb) {
          pdb_file.close();
     } else {
          xdrfile_close(xd);          
     }

     if (coordinates)
          delete coordinates;

}


//! Main function
int main(int argc, char *argv[]) {

     // Declare the supported options.
     po::options_description desc("Options");
     desc.add_options()
          ("help", "produce help message")
          ("pdb-file", po::value<std::string>(), "PDB file")
          ("chain-type", po::value<std::string>()->default_value("chainfb"),
           "Chain type (chainfb|chainca)")
          ("trajectory-files", po::value<std::vector<std::string> >()->multitoken(), 
           "XTC files")
          ("iteration-start", po::value<PHAISTOS_LONG_LONG>()->default_value(0), 
           "Iteration range: start index")
          ("iteration-end", po::value<PHAISTOS_LONG_LONG>()->default_value(-1), 
           "Iteration range: end index")
          ("seed", po::value<unsigned int>()->default_value(time(NULL)),
           "Seed for random number generator. The default value is the current time. ")
          ("samples", po::value<PHAISTOS_LONG_LONG>(), 
           "In case weight file is specified, this specifies how many samples to draw.")
          ("observable-files", po::value<std::vector<std::string> >()->multitoken(), 
           "Observables files containing for instance weights or energies")
          ("weight-column", po::value<int>()->default_value(0),
           "Column index in observables file containing weights")
          ("output-file", po::value<std::string>()->default_value("trajectory_subset.xtc"),
           "Output file")
#ifdef HAVE_MUNINNLIB     
          ("muninn-log-file", po::value<std::string>(), 
           "Path to muninn log file")
          ("muninn-column", po::value<int>(), 
           "Column index in observables file containing reaction coordinate used with muninn.")
#endif
          ;

     // Parse options
     po::variables_map options;
     po::store(po::parse_command_line(argc, argv, desc), options);
     po::notify(options);    

     // Set the random seed
     random_global.seed(options["seed"].as<unsigned int>());

     // Help output
     if (options.count("help")) {
          std::cerr << desc << "\n";
          return 1;
     }

     // Call main function either in Torus or FB5 mode
     if (options["chain-type"].as<std::string>()=="chainfb") {
          trajectory_subset_main<ChainFB>(options);
     } else {
          trajectory_subset_main<ChainCA>(options);
     }

}

