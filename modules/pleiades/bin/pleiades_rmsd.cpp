// pleiades_rmsd.cpp --- rmsd based k-means clustering
// Copyright (C) 2011  Tim Harder, Thomas Hamelryck, Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#include <stdio.h>
#include <string>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <vector>
#include <iostream>	
#include <fstream>
#include <algorithm>
#include <cstring>

#include "utils/vector_matrix_3d.h"
#include "protein/chain_ca.h"
#include "protein/chain_fb.h"
#include "protein/chain.h"
#include "energy/term_rmsd.h"

#ifdef PHAISTOS_VERSION
#ifdef SVN_REVISION
#include "revision.h"
#endif
#endif

#include "models/backbone_dbn/nodes/discrete.h"

#ifdef HAVE_MUNINNLIB
#include "muninn/tools/CanonicalAveragerFromStatisticsLog.h"
#endif

#include "cluster_rmsd.h"
#include "pleiades_rmsd_options.h"

using namespace phaistos;

//! Load .pdb file as a ChainCA object.
//!
//! \param filename Filename of .pdb file.
//! \param chains vector pointer to return the chains found.
//! \param chain_names vector pointer to return all the filenames.
//! \param options runtime options.
//! \param energies Vector in which energies can be written.
//! \param extract_energies Whether to attempt to extract energy from pdb filename
bool get_pdb(std::string filename, std::vector<ChainCA> *chains, std::vector<std::string> *chain_names, PleiadesRmsdOptions options, std::vector<double> &energies, bool extract_energies) {

     if (options.verbose)
          std::cout << "Found file: " << filename << " ";

     // filename ends with .pdb
     if ((filename[filename.size() - 4] == '.' && 
          filename[filename.size() - 3] == 'p' && 
          filename[filename.size() - 2] == 'd' &&
          filename[filename.size() - 1] == 'b')) {
          if (options.verbose)
               std::cout << " .. reading. ";

          // If we haven't read energies yet, attempt to extract them from filename
          if (extract_energies) {
               std::string basename = filename.substr(0, filename.size() - 4);
               std::vector<std::string> tokens;
               boost::split(tokens, basename , boost::is_any_of("_"));
               energies.push_back(atof(tokens.back().c_str()));
               if (options.verbose)
                    std::cout << "Extracted energy=" << energies.back() << " ";
          }

          if (!file_exists(filename)) {
               std::cout << ".. not found.\n";
               return false;
          }

          ProteinData data = read_pdb_input(filename);
          if (data.n_polypeptides() > 1) {
               printf(" ... FAILED (multiple fragments)\n");
               return false;
          }

          ChainCA *tmp_chain = new ChainCA(data);
          tmp_chain->set_name(filename);
          chains->push_back(*tmp_chain);
          chain_names->push_back(std::string(filename));
          if (options.debug)
               std::cout << " length: " << chains->size();
          if (options.verbose)
               std::cout << std::endl;
     } else {
          if (options.verbose)
               std::cout << " .. skipped.\n";
          return false;
     }
     return true;
}

//! Read pdb file names from a file and load all .pdb files into a vector as
//! ChainCA objects. Further return all the filenames found
//! as a vector.
//!
//! \param pdb_list_filename Filename of file containing list of pdb filenames.
//! \param chains vector pointer to return the chains found.
//! \param chain_names vector pointer to return all the filenames.
//! \param options runtime options.
//! \param energies Vector in which energies can be written.
//! \param weights Vector in which weights can be written.
void get_pdb_files_from_file(std::string pdb_list_filename, std::vector<ChainCA> *chains, std::vector<std::string> *chain_names, PleiadesRmsdOptions options, std::vector<double> &energies, std::vector<double> &weights) {

     std::vector<std::string> pdb_filenames = file_to_string_vector(pdb_list_filename);

     for (unsigned int i=0; i<pdb_filenames.size(); i++) {
          std::string &pdb_filename = pdb_filenames[i];
          
          std::vector<std::string> tokens;
          boost::split(tokens, pdb_filename, boost::is_any_of(" \t"), boost::token_compress_on);
          if (tokens.size() > 1) {
               pdb_filename = tokens[0];               
          }

          bool extract_energies = true;
          if (tokens.size() > 2)
               extract_energies = false;

          bool success = get_pdb(pdb_filename, chains, chain_names, options, energies, extract_energies);

          if (success) {
               if (tokens.size() > 1) {
                    weights.push_back(atof(tokens[1].c_str()));
                    assert(weights.back() != 0);
               }
               if (tokens.size() > 2) {
                    energies.push_back(atof(tokens[2].c_str()));
               }
          }
     }

     assert(energies.size()==0 || chains->size() == energies.size());
     assert(weights.size()==0 || chains->size() == weights.size());
}

//! Parse a directory and load all .pdb files into a vector as
//! ChainCA objects. Further return all the filenames found
//! as a vector.
//!
//! \param dirname directory name to search for pdb files.
//! \param chains vector pointer to return the chains found.
//! \param chain_names vector pointer to return all the filenames.
//! \param options runtime options.
//! \param energies Vector in which energies can be written.
//! \param weights Vector in which weights can be written.
void get_pdb_files_from_directory(const char *dirname, std::vector<ChainCA> *chains, std::vector<std::string> *chain_names, PleiadesRmsdOptions options, std::vector<double> &energies, std::vector<double> &weights) {

     struct dirent *dp;
     DIR *dfd;
     char filename[200];
     std::string name;

     if ((dfd = opendir(dirname)) == NULL) {
          printf("Sorry - could not open directory <%s>\n", dirname);
          return;
     }
     int i = 0;
     std::vector<std::string> pdb_filenames;
     while ((dp = readdir(dfd)) != NULL) {
          if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
               continue;

          sprintf(filename, "%s/%s", dirname, dp->d_name);

          pdb_filenames.push_back(filename);

          i++;
     }

     bool extract_energies = true;
     for (unsigned int i=0; i<pdb_filenames.size(); ++i) {
          get_pdb(pdb_filenames[i], chains, chain_names, options, energies, extract_energies);
     }

     closedir(dfd);
     return;
}


//! parse a directory and return all the filenames found
//! as a vector.
//!
//! \param dirname directory name to search for pdb files.
//! \param filenames vector reference to return the found PDB files
void get_all_filenames(char *dirname, std::vector<std::string> &filenames) {

     struct dirent *dp;
     DIR *dfd;
     char filename[200];
     std::string name;

     if ((dfd = opendir(dirname)) == NULL) {
          printf("Sorry - could not open directory <%s>\n", dirname);
          return;
     }

     while ((dp = readdir(dfd)) != NULL) {
          if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
               continue;

          sprintf(filename, "%s/%s", dirname, dp->d_name);
          // filename ends with .pdb
          if ((dp->d_name[strlen(dp->d_name) - 4] == '.' && 
               dp->d_name[strlen(dp->d_name) - 3] == 'p' && 
               dp->d_name[strlen(dp->d_name) - 2] == 'd' && 
               dp->d_name[strlen(dp->d_name) - 1] == 'b')) {
               filenames.push_back(filename);
          }
     }

     closedir(dfd);
     return;
}

//! Output a header block containing the most important runtime variables.
//!
//! \param options Pleiades runtime options.
//! \param o output stream to write to.
void output_header(PleiadesRmsdOptions options, std::ofstream &o) {
     o << "############################################################" << std::endl;
     o << "#\t" << std::endl;
     o << "#\toutput generated by PLEIADES_RMSD	" << std::endl;
     o << "#\t" << std::endl;
     o << "#\tParameters used for the run" << std::endl;
     o << "#\t" << std::endl;
     o << "#\tinput dir :\t" << options.decoy_dir.c_str() << std::endl;
     o << "#\titerations :\t" << options.iterations << std::endl;
     o << "#\tK :\t\t" << options.k << std::endl;
     o << "#\t" << std::endl;
     o << "############################################################\n";
}

//! Output a vector of clusters. If the long_output runtime option has been
//! set, all member of the cluster will be written to the given stream,
//! otherwise, only the five members closest to the centroid are written.
//!
//! \param options Pleiades runtime options.
//! \param clusters vector of clusters.
//! \param o output stream to write to.
//! \param native Native molecule chain
//! \param energies Optional vector of energies - used for calculation of median
void output_cluster(PleiadesRmsdOptions options, std::vector<ClusterRmsd> &clusters, std::ofstream &o, ChainCA *native,
                    std::vector<double> energies=std::vector<double>()) {
     std::sort(clusters.begin(), clusters.end(), std::greater<ClusterRmsd>());
     output_header(options, o);
     // tmp string for formatting
     char tmp[1000];

     for (unsigned int i = 0; i < clusters.size(); i++) {
          sprintf(tmp, "[Cluster %4d (size: %10.4f, sd: %10.4f, dens: %10.4f ) ]", (clusters[i]).get_id(), (clusters[i]).get_weight_sum(), (clusters[i]).get_std_dev(),
                    (clusters[i]).get_density());
          o << tmp << std::endl;
          printf("%010.4f -> Cluster %4d (size: %10.4f, sd: %10.4f, dens: %10.4f \t%s)\n", (clusters[i]).get_density(), (clusters[i]).get_id(),
                 (clusters[i]).get_weight_sum(), (clusters[i]).get_std_dev(), (clusters[i]).get_density(), clusters[i].size() > 0 ? (((clusters[i]).get_median(energies))->get_name()).c_str() : "");

          if (options.long_output) {
               (clusters[i]).dump_member_list(o, native, -1, options.dump_mean, options.dump_median, energies);
          } else {
               (clusters[i]).dump_member_list(o, native, 5, options.dump_mean, options.dump_median, energies);
          }
     }
     printf("\n");
}


//! Weighted version of std::random_shuffle. Defaults to
//! std::random_shuffle if weights are not given.
//!
//! \param chains_access Vector of chain indices
//! \param weights Vector of weights
void random_shuffle(std::vector<int> &chains_access, std::vector<double> *weights) {
     
     if (!weights) {
          std::random_shuffle(chains_access.begin(), chains_access.end());
     } else {
          double *weights_local = new double[weights->size()];
          for (unsigned int i=0; i<weights->size(); i++) {
               chains_access[i] = i;
               weights_local[i] = (*weights)[i];
          }

          for (unsigned int i=0; i<weights->size(); i++) {
               double *weights_local_subset = weights_local+i;
               int index = DiscreteSampler::sample(weights->size()-i, weights_local_subset)+i;
               std::swap(chains_access[i], chains_access[index]);
               std::swap(weights_local[i], weights_local[index]);
          }
     }
}

//! Performes the actual k-means clustering on a set of GIT vectors.
//!
//! \param my_k number of cluster centers k.
//! \param options Pleiades runtime options.
//! \param chains vector of molecule chain objects
//! \param chain_names vector of chain names
//! \param chains_access Vector of chain indices
//! \param cluster Vector of clusters
//! \param weights Vector of weights
void k_means(int my_k, PleiadesRmsdOptions &options, std::vector<ChainCA> &chains, std::vector<std::string> &chain_names, 
             std::vector<int> &chains_access, std::vector<ClusterRmsd> &cluster, std::vector<double> *weights = NULL) {

     if (options.debug && weights == NULL) {
          std::cout << "D> K-mean setup. \n";
          std::cout << std::flush;
     } else if (options.debug) {
          std::cout << "D> weighted K-mean setup. \n";
          std::cout << std::flush;
     }

     double dist;
     int len = (int) chains.size();
     int cluster_id = 0;

     //ClusterRmsd *tmp = NULL;
     std::vector<Vector_3D> pos1;
     std::vector<Vector_3D> pos2;

     if (chains.size() != chains_access.size()) {
          std::cerr << "The chain vector and the access vector need to be equally long.\n";
          std::cerr.flush();
          assert(false);
     }

     if (options.debug) {
          std::cout << "D> Setting up initial K cluster ... \n";
          std::cout << std::flush;
     }

     if (my_k != (int) cluster.size()) {
          // not enough cluster center yet. Add additional random ones
          // std::random_shuffle(chains_access.begin(), chains_access.end());
          random_shuffle(chains_access, weights);
          cluster.clear();


          // Fill in k
          for (int i = 0; i < my_k; i++) {

               cluster_id++;

               ClusterRmsd *tmp2 = NULL;
               tmp2 = new ClusterRmsd(cluster_id, &chains, &chain_names, weights, options.residue_start, options.residue_end);

               if (options.debug) {
                    std::cout << "D>     cluster ID : " << cluster_id << " created (length:" << cluster.size() << " ) - " << chain_names[chains_access[i]] << "\n";
                    std::cout << std::flush;
               }

               tmp2->add_member(chains_access[i], false);
               cluster.push_back(*tmp2);
               delete tmp2;
          }
     }

     for (int k = 0; k < options.iterations; k++) {
          if (options.verbose) {
               printf("Calculating Kmeans cluster (with k = %d) [%2d / %2d]\n", (int) cluster.size(), k + 1, options.iterations);
          }
          // shuffle the order
          random_shuffle(chains_access, weights);
          // std::random_shuffle(chains_access.begin(), chains_access.end());

          // loop over all clusters and free their members
          // only interesting in the later rounds.
          for (unsigned int i = 0; i < cluster.size(); i++) {
               cluster[i].free_members();
          }

          // Assign all the points to the clusters again.
          for (int i = 0; i < len; i++) {

               double min_d = std::numeric_limits<double>::infinity();
               int best_cluster = -1;

               for (unsigned int j = 0; j < cluster.size(); j++) {
                    ChainCA *tmp_mean = cluster[j].get_mean();
                    dist = calc_rmsd<definitions::CA_ONLY> ((chains[chains_access[i]]), *tmp_mean, options.residue_start, options.residue_end);

                    if (dist < min_d) {
                         best_cluster = j;
                         min_d = dist;
                    }
               }

               if (best_cluster >= 0 && best_cluster < (int) cluster.size()) {
                    cluster[best_cluster].add_member(chains_access[i], false);
               }
          }

          // Recalculate mean centeroids
          for (unsigned int i = 0; i < cluster.size(); i++) {
               cluster[i].calc_mean();
          }
     }
}

void weight_clusters(PleiadesRmsdOptions &options, std::vector<ChainCA> &chains, std::vector<std::string> &chain_names, 
                     std::vector<int> &chains_access, std::vector<ClusterRmsd> &cluster, std::vector<double> *weights = NULL) {
     
     std::vector<ClusterRmsd> cluster_tmp;
     int len = (int) chains.size();
     double dist;
     
     for (int i = 0; i < (int)cluster.size(); i++) {

          ClusterRmsd *tmp2 = NULL;
          tmp2 = new ClusterRmsd(cluster[i].get_id(), &chains, &chain_names, weights, options.residue_start, options.residue_end);

          cluster_tmp.push_back(*tmp2);
          delete tmp2;
     }
     
     // Assign all the points to the new clusters.
     for (int i = 0; i < len; i++) {

          double min_d = std::numeric_limits<double>::infinity();
          int best_cluster = -1;

          for (unsigned int j = 0; j < cluster_tmp.size(); j++) {
               ChainCA *tmp_mean = cluster[j].get_mean();
               dist = calc_rmsd<definitions::CA_ONLY> ((chains[chains_access[i]]), *tmp_mean, options.residue_start, options.residue_end);

               if (dist < min_d) {
                    best_cluster = j;
                    min_d = dist;
               }
          }
          
          if (best_cluster >= 0 && best_cluster < (int) cluster_tmp.size()) {
               cluster_tmp[best_cluster].add_member(chains_access[i], false);
          }
     }
     
     cluster = cluster_tmp;
     
     // Recalculate mean centeroids
     for (unsigned int i = 0; i < cluster.size(); i++) {
          cluster[i].calc_mean();
     }
}


//! Procedure to guess a good number of k clusters given a distance threshold
//! \param options Pleiades runtime options.
//! \param chains vector of molecule chain objects
//! \param chain_names vector of chain names
//! \param chains_access Vector of chain indices
//! \param cluster Vector of clusters
void guess_k(PleiadesRmsdOptions &options, std::vector<ChainCA> &chains, std::vector<std::string> &chain_names, 
             std::vector<int> &chains_access, std::vector<ClusterRmsd> &cluster) {

     double dist;
     int len = (int) chains.size();
     int cluster_id = 0;
     bool is_stored;
     int threshold = 10;
     int opt_cluster = 10;

     ClusterRmsd *tmp = NULL;
     std::vector<Vector_3D> pos1;
     std::vector<Vector_3D> pos2;

     if (chains.size() != chains_access.size()) {
          std::cerr << "The chain vector and the access vector need to be equally long.\n";
          std::cerr.flush();
          assert(false);
     }

     for (int k = 0; k < opt_cluster; k++) {
          if (options.verbose) {
               printf("Calculating near cluster [%2d / %2d]\n", k + 1, opt_cluster);
          }

          // shuffle the order
          std::random_shuffle(chains_access.begin(), chains_access.end());

          // loop over all clusters and free their members
          // only interesting in the later rounds.
          for (unsigned int i = 0; i < cluster.size(); i++) {
               cluster[i].free_members();
          }

          for (int i = 0; i < len; i++) {
               is_stored = false;

               for (unsigned int j = 0; j < cluster.size(); j++) {
                    ChainCA *tmp_mean = cluster[j].get_mean();
                    dist = calc_rmsd<definitions::CA_ONLY> ((chains[chains_access[i]]), (*tmp_mean), options.residue_start, options.residue_end);

                    if (dist < options.guess_k_threshold) {
                         is_stored = true;
                         cluster[j].add_member(chains_access[i], false);
                         break;
                    }
               }
               if (!is_stored) {
                    if (options.debug)
                         printf("D> Created new cluster %d\n", cluster_id);

                    tmp = new ClusterRmsd(cluster_id++, &chains, &chain_names, NULL, options.residue_start, options.residue_end);
                    tmp->add_member(chains_access[i]);
                    cluster.push_back(*tmp);
               }
          }

          // Get rid of the mini or empty clusters.
          for (int i = 0; i < (int) cluster.size(); i++) {
               // keep all clusters in the last iteration
               if (cluster[i].get_density() < threshold) {
                    cluster.erase(cluster.begin() + i);
                    i--;
               }
          }

          // Recalculate mean centeroids
          for (unsigned int i = 0; i < cluster.size(); i++) {
               cluster[i].calc_mean();
          }
     }
     delete tmp;
}


//! Pleiades main function
//! \param argc commandline argument number
//! \param argv commandline argument value
int main(int argc, char *argv[]) {

     PleiadesRmsdOptions options(argc, argv);

#ifdef PHAISTOS_VERSION
#ifdef SVN_REVISION     
     if (options.verbose) {
          options.print_options();
          printf("\n######################################################################################################\n");
          printf("# Version: %5s  Build: %5s \n", PHAISTOS_VERSION, SVN_REVISION);
          printf("######################################################################################################\n");
     }
#endif
#endif

     uint now = (uint) time(NULL);
     if (options.verbose)
          std::cout << "Seed: " << now << "\n";
     srand(now);
     random_global.seed(now);

     if (options.decoy_dir == "" && options.pdb_list_file=="") {
          options.print_usage();
          return 0;
     }

     if (options.debug)
          options.verbose = true;

     std::vector<ChainCA> chains;
     std::vector<std::string> chain_names;
     std::vector<int> chains_access;
     std::vector<int> chains_access_copy;
     std::vector<ClusterRmsd> cluster;
     std::vector<double> energies;
     std::vector<double> weights;

     ChainCA *native = NULL;
     if (options.native_file != "") {
          native = new ChainCA(options.native_file);
     }

     // Parse input file
     if (options.debug)
          std::cout << "D> Reading the input file ... ";

     // Use pdb-list-file if specified
     if (options.pdb_list_file != "")
          get_pdb_files_from_file(options.pdb_list_file, &chains, &chain_names, options, energies, weights);
     else 
          get_pdb_files_from_directory(options.decoy_dir.c_str(), &chains, &chain_names, options, energies, weights);

     if (options.debug)
          std::cout << "done (read " << chains.size() << " vectors) \n";

     // Check whether parsing worked
     if (not chains.size()) {
          std::cerr << "ERROR: The input directory appears to be empty or only contain none-readable pdb files.\n\n";
          return 256;
     }

     // output file
     std::ofstream out_file(options.out_file.c_str());

     // Create a vector of ids for use when accessing the chains vector
     // i.e. we need randomized access, but the chains vector
     // needs to stay in the correct order.
     for (int i = 0; i < (int) chains.size(); i++) {
          chains_access.push_back(i);
     }
     chains_access_copy = chains_access;

     int my_k = options.k;

     // are we asked to guesstimate K?
     if (options.guess_k) {
          if (options.verbose)
               std::cout << "Trying to estimate K.\n";
          guess_k(options, chains, chain_names, chains_access, cluster);
          my_k = cluster.size();
          std::cout << "Using a threshold of " << options.guess_k_threshold << " - k was estimated to be " << my_k << "\n";
     }

     // do we want to run a real k-means
     if (options.k_means) {
          if (options.verbose)
               std::cout << "Starting K-Means (K=" << my_k << ")\n";

          k_means(my_k, options, chains, chain_names, chains_access, cluster);
     }

     if (options.w_k_means || options.evaluate_clusters_weighted) {

          if (weights.empty()) {

#ifdef HAVE_MUNINNLIB

               if (options.muninn_log != "" && energies.size() == chains.size()) {

                    if (options.verbose)
                         std::cout << "Starting weighted K-Means (K=" << my_k << ", beta=" << options.beta << ")\n";
                    // now we need to read in the histograms
                    // and have muninn calculate the correct weights from the histograms
                    if (options.debug)
                         std::cout << "Reading Muninn Weights from: " << std::string(options.muninn_log).c_str() << "\n";

                    Muninn::CanonicalAveragerFromStatisticsLog averager = Muninn::CanonicalAveragerFromStatisticsLog(options.muninn_log);
                    weights = averager.calc_weights(energies, options.beta);

                    double max = -std::numeric_limits<double>::infinity();
                    double min = std::numeric_limits<double>::infinity();
                    int zeros = 0;
                    for (int k = 0; k < (int) weights.size(); k++) {
                         if (weights[k] > max)
                              max = weights[k];
                         if (weights[k] < min && weights[k] != 0)
                              min = weights[k];
                         if (weights[k] == 0)
                              zeros++;
                    }

                    if (options.scale_weights) {
                         double factor = 1 / max;
                         for (int k = 0; k < (int) weights.size(); k++) {
                              if (weights[k] == 0)
                                   weights[k] = min;
                              weights[k] = weights[k] * factor;
                         }
                    } else if (zeros) {
                         std::cout << std::endl << std::endl << "WARNING: Performing weighted clustering and found " << zeros
                                   << " structures with weight zero.\nWARNING: This usually leads to weird results .. proceed with caution or enable scaling (--scale-weights). "
                                   << std::endl << std::endl << std::endl;
                    }

               } else 
#endif
               if (energies.size() == chains.size()) {
                    if (options.verbose)
                         std::cout << "Starting weighted K-Means (ab)using energies as weights (K=" << my_k << ", beta=" << options.beta << " [scaling factor])"
                                   << std::endl;
                    for (int k = 0; k < (int) energies.size(); k++) {
                         weights.push_back(energies[k] * options.beta);
                    }
               } else {
                    std::cerr << "ERROR: Requested weighted clustering, but no weights were set!" << std::endl;
                    return 256;
               }
          }
          if (options.w_k_means) {
               k_means(my_k, options, chains, chain_names, chains_access, cluster, &weights);
          }
     }
    
     //Evaluate clusters weighted:
     if (options.evaluate_clusters_weighted && !options.w_k_means) {
          weight_clusters(options, chains, chain_names, chains_access_copy, cluster, &weights);
     }

    
     output_cluster(options, cluster, out_file, native, energies);
     
     std::cout << "\n\n\n";

     return 0;
}

