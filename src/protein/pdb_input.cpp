// input.cpp --- Simple PDB parser
// Copyright (C) 2006-2008 Wouter Boomsma
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

#include "pdb_input.h"
#include <boost/algorithm/string/trim.hpp>

namespace phaistos {

// Translate sequence string to vector
std::vector<int> seq_str_to_vec(std::string seq) {
     std::vector<int> seq_vec;

     for (unsigned int i=0; i<seq.length(); i++) {
          if (!(seq[i] == '\n')) {
               seq_vec.push_back(definitions::str_to_aa[seq.substr(i, 1)]);
          }
     }
     return seq_vec;
}


// Read input data
ProteinData read_pdb_input(std::string filename) {
     // Open file stream
     std::ifstream ifs(filename.c_str());
     if (!ifs) {
          std::cerr<<"Cannot read pdb file " << filename << "\n";
          assert(false);
     }
     
     // Get contents of PDB file
     std::string pdb_data = std::string(std::istreambuf_iterator<char>(ifs),
                                        std::istreambuf_iterator<char>());
     ifs.seekg(0);

     // Atom data
     ProteinData data;
     
     // Whether chain index is included in PDB file
     bool include_chain=0;

     std::string current_chain = "";
     std::string prev_chain;
     int chain_index = -1;
     
     std::string current_res = "-1000";
     std::string prev_res;
     int res_index = -1;

     bool only_CA=true;
     bool skip_to_next_chain = false;
     bool skip_to_next_res = false;

     std::string line;
     while (getline(ifs, line)) {

          boost::trim(line);

          // Check for termination 
          if (line.substr(0,3)=="TER") {
               skip_to_next_chain = true;
          }

          // Check whether line contains atom coordinates
          if (line.substr(0,4)!="ATOM") {
               continue;
          }
          
          // // Skip weird amino acids
          // if (map_lookup(strToRes,res,int(AA_UNDEF)) == AA_UNDEF) {
          //      continue;
          // }

          std::string res = strip(line.substr(17,3));

          prev_res = current_res;
          current_res = strip(line.substr(22,5)).c_str();

          std::string atom = strip(line.substr(12, 4));
          if (atom != "CA")
               only_CA = false;
          
          // Check for chain index
          if (line[21]!=' ') {
               include_chain=1;
               prev_chain = current_chain;
               current_chain = line[21];
          } else {
               prev_chain = current_chain;
               current_chain = " ";            
               include_chain=0;
          }

          if ((prev_chain != current_chain)) {
               data.new_polypeptide();
               chain_index++;
               res_index = -1;
               prev_res = -1;
               skip_to_next_chain = false;
          } else if (skip_to_next_chain) {
               continue;
          }

          // Skip if first atom in a residue is not a N, CA, C, O, CB
          if (current_res!=prev_res) {
               skip_to_next_res = false;
               
               if (!(atom == "N")) {
                    skip_to_next_res = true;
                    continue;
               } 
          } else if (skip_to_next_res) {
               continue;
          }
                  
          if (current_res!=prev_res) {
               data.new_residue(chain_index);
               res_index++;
          }

          data.new_atom(chain_index, res_index);

          std::string altloc = line.substr(16, 1);
          

          double coordx = atof(line.substr(30, 8).c_str());
          double coordy = atof(line.substr(38, 8).c_str());
          double coordz = atof(line.substr(46, 8).c_str());

          data.current_atom()->add_atom(atom);
          data.current_atom()->add_altloc(altloc);
          data.current_atom()->add_residue(res);
          data.current_atom()->add_coordinates(coordx,coordy,coordz);
          data.current_atom()->add_resseq(atoi(current_res.c_str()));
               
     }

     // Scan data for chain breaks, and split data if necessary
     data.split_by_chain_breaks(only_CA);
     
     // return data
     return data;
}


// Translate SS character to index
int ss_num(char ss) {
     if (ss == 'H' || ss == 'G' || ss == 'I')
          return 0;
     else if (ss=='E' || ss=='B')
          return 1;
     else
          return 2;
}

// Translate sequence string to vector
std::vector<int> ss_str_to_vec(std::string seq) {
     std::vector<int> seq_vec;

     for (unsigned int i=0; i<seq.length(); i++) {
          if (!(seq[i] == '\n'))
               seq_vec.push_back(ss_num(seq[i]));
     }
     return seq_vec;
}

// Call DSSP and parse output - returns string
std::string get_dssp_string(std::string pdb_filename, std::string dssp_command) {

     // Get DSSP SS-assignment 
     std::string commandline_ss = dssp_command + " " + pdb_filename + std::string(" 2>/dev/null | grep -Ev \"\\.$\" | grep -v \"  #\" | cut -b 17 | sed -e \"s/[ ]/C/g\" | tr -d \"\\n\"");
     std::string output_ss = unix_command(commandline_ss);

     if (output_ss.size() == 0) {
          std::cerr << "Warning: DSSP command failed. No secondary structure input\n";
          return "";
     }

     // Find chains in DSSP output
     std::string commandline_chain = dssp_command + " " + pdb_filename + std::string(" 2>/dev/null | grep -Ev \"\\.$\" | grep -v \"  #\" | cut -b 12 | tr -d \"\\n\"");
     std::string output_chain = unix_command(commandline_chain);
     std::string empty_string = std::string(" "); empty_string.resize(output_chain.size(), ' ');
     if (output_chain != empty_string) { // check for empty string: single chain
          size_t position = output_chain.size()-1;
          while (true) {
               position = output_chain.rfind(" ", position);
               
               if (position == std::string::npos)
                    break;
               else {
                    output_ss.replace(position, 1, " ");
                    position -= 1;
               }
          }
     }

     // Detect chain breaks in DSSP output
     // std::string commandline_chain_breaks = dsspCommand + " " + pdbFilename + std::string(" 2>/dev/null | grep -Ev \"\\.$\" | grep -v \"  #\" | cut -b 6-10 | tr \"\\n\" \"#\"");
     // std::string output_chain_breaks = unix_command(commandline_chain_breaks);
     // std::vector<std::string> output_chain_breaks_list = split(output_chain_breaks,"#");
     // for (unsigned int i=0; i<output_chain_breaks_list.size(); i++) {
     //           if (strip(output_chain_breaks_list[i]) == "") {
     //                output_SS.replace(i, 1, " ");
     //           }
     // }

     return output_ss;
}

// Call DSSP and parse output - splits according to chain
std::vector<std::vector<int> > get_dssp(std::string pdb_filename, std::string dssp_command) {

     std::vector<std::vector<int> > output_vec;

     std::string output_ss = get_dssp_string(pdb_filename, dssp_command);

     if (output_ss == "")
          return output_vec;

     // divide output according to chain and chain-breaks
     std::vector<std::string> split_output = split(output_ss, " ");

     for (unsigned int i=0; i<split_output.size(); i++) {
          output_vec.push_back(std::vector<int>());

          for (unsigned int j=0; j<split_output[i].size(); j++) {
               output_vec[i].push_back(ss_num(split_output[i][j]));
          }
     }
     return output_vec;
}


}
