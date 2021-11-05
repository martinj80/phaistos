// input.h --- Simple PDB parser
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

#ifndef INPUT_H
#define INPUT_H

#include "protein_data.h"

namespace phaistos {

//! Read input from a pdb file
//!
//! \param filename PDB filename
//! \return ProteinData structure
ProteinData read_pdb_input(std::string filename);

//! Convert amino acid sequence string to a vector
//!
//! \param seq Amino acid input string
//! \return Vector of amino acid indices (ResidueEnum)
std::vector<int> seq_str_to_vec(std::string seq);

//! Convert secondary structure sequence string to a vector
//!
//! \param seq Secondary structure input string
//! \return Vector of secondary structure indices (SecondaryStructureEnum)
std::vector<int> ss_str_to_vec(std::string seq);


//! Calculate DSSP secondary structure classification (calls external problem)
//!
//! \param pdb_filename PDB filename
//! \param dssp_command Location of dssp program executable
//! \return dssp string
std::string get_dssp_string(std::string pdb_filename, std::string dssp_command="dssp");

//! Calculate DSSP secondary structure classification (calls external problem)
//!
//! \param pdb_filename PDB filename
//! \param dssp_command Location of dssp program executable
//! \return Vector of secondary structure indices (SecondaryStructureEnum)
std::vector<std::vector<int> > get_dssp(std::string pdb_filename, std::string dssp_command="dssp");

}

#endif
