// pdb2git.cpp --- convert pdb files into git vector
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
#include <dirent.h>
#include <vector>
#include <iostream>	
#include <cstring>

#ifdef PHAISTOS_VERSION
#ifdef SVN_REVISION
#include "revision.h"
#endif
#endif

#include "protein/chain_ca.h"
#include "protein/chain.h"
#include "git.h"

using namespace phaistos;

//! Simple commandline argument parser for the pdb2git program
class Options {
private:

     // For interal use
     std::string command_line;

public:

     //! output file in which to store the converted vector(s)
     std::string git_file;

     //! Specifies a single filename to be parsed
     std::string input_file;

     //! Specifies an entire directory of input files
     std::string input_dir;

     //! Whether to be verbose with the output
     bool verbose;

     //! Write debug information to standard out
     bool debug;


     //! Constructor
     //!
     //! \param argc commandline argument count
     //! \param argv commandline arguments
     Options(int argc, char *argv[]) {
          get_command_line(argc, argv);
     }


     //! Destructor
     ~Options() {
     }


     //! Specify default values for all
     //! arguments.
     void init() {
          git_file = "";
          input_file = "";
          input_dir = "";
          verbose = false;
          debug = false;
          command_line = "";
     }


     //! Parse the command line
     //!
     //! \param argc commandline argument count
     //! \param argv commandline arguments
     void get_command_line(int argc, char *argv[]) {
          init();

          for (int i = 0; i < argc; i++) {
               std::string arg = argv[i];
               std::string arg_name = "";

               // keyword or argument ?
               if (arg.substr(0, 2) == std::string("--")) {
                    arg_name = arg.substr(2, arg.length() - 2);
               } else if (arg.substr(0, 1) == std::string("-")) {
                    arg_name = arg.substr(1, arg.length() - 1);
               } else if (i != 0) {
                    // we want to actually store the progname
                    continue;
               }

               if (i == 0) {
                    command_line = arg;
               } else if (arg_name == std::string("help")) {
                    print_usage();
                    exit(0);
               } else if (arg_name == std::string("git-file")) {
                    i++;
                    git_file = argv[i];
               } else if (arg_name == std::string("input-file")) {
                    i++;
                    input_file = argv[i];
               } else if (arg_name == std::string("input-dir")) {
                    i++;
                    input_dir = argv[i];
               } else if (arg_name == std::string("verbose")) {
                    verbose = true;
               } else if (arg_name == std::string("debug")) {
                    debug = true;
               } else {
                    std::cerr << "Unknown commandline option " << arg << " will be ignored! \n";
               }
          }

     }

     //! Write the help message to the standard output.
     //! All the class attributes will be reset to their
     //! default settings in the process.
     void print_usage() {
          // reset those values
          init();
          // generate usage informations
          std::cout << "Usage " << command_line << " [options] \n";
          std::cout << "\toptions: \n";
          std::cout << "\t--git-file\t\tgit vector file [required] \n";
          std::cout << "\t--input-file\t\tconvert a single file [optinal] \n";
          std::cout << "\t--input-dir\t\tconvert all files in a directory[optional] \n";
          std::cout << "\n";
          std::cout << "\t--verbose\t\tbe talky? [default=" << verbose << "] \n";
          std::cout << "\t--debug\t\t\tprint debugging information [default=" << debug << "] \n";
          std::cout << "\n";
     }
};


//! Find all the .pdb files in a given directory and return them as
//! a vector of filenames.
//!
//! \param dirname name of the input directory to be searched
//! \param chains vector reference to return all the filenames found
void get_all_pdb_files(const char *dirname, std::vector<std::string> &chains) {

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
          if (dp->d_name[strlen(dp->d_name) - 4] == '.' && dp->d_name[strlen(dp->d_name) - 3] == 'p' && dp->d_name[strlen(dp->d_name) - 2] == 'd'
                    && dp->d_name[strlen(dp->d_name) - 1] == 'b') {
               chains.push_back(std::string(filename));
          }
     }

     closedir(dfd);
     return;
}

//! Main function
//!
//! \param argc commandline argument count
//! \param argv commandline arguments
int main(int argc, char *argv[]) {

     Options options = Options(argc, argv);

     if ((options.input_dir == "" && options.input_file == "") || options.git_file == "" ) {
          // otherwise print the usage
          options.print_usage();
          return 0;
     }
     
#ifdef PHAISTOS_VERSION
#ifdef SVN_REVISION
	if (options.verbose) {
    	 printf("\n######################################################################################################\n");
         printf("# Version: %5s  Build: %5s \n" , PHAISTOS_VERSION, SVN_REVISION);
         printf("######################################################################################################\n");
     }
#endif
#endif

     ///////////
     // Variables
     std::vector<std::string> chains;

     if (options.input_dir != "") {
          if (options.verbose)
               printf("Reading PDBfile(s) ..\n");
          get_all_pdb_files(options.input_dir.c_str(), chains);
     } else {
          if (options.verbose)
               printf("Reading single PDBfile ..\n");
          chains.push_back(options.input_file);
     }

     Git git(options.git_file.c_str());

     int len = chains.size();

     /* */
     for (int i = 0; i < len; i++) {
          if (options.verbose)
               printf("[%8d / %8d]%s", i, len, chains[i].c_str());

          ProteinData data = read_pdb_input(chains[i].c_str());
          if (data.n_polypeptides() > 1) {
        	  printf(" ... FAILED (multiple fragments)\n");
        	  continue;
          }
          printf("\n");

          ChainCA tmp_chain(data);
          git.generate_gauss_integrals(tmp_chain, chains[i].c_str());
     }

     if (options.verbose)
          printf("done\n");

     return 0;
}

