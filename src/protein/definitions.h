// definitions.h --- Values and enums used by proteins
// Copyright (C) 2006-2011 Wouter Boomsma
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


#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "utils/utils.h"
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string.hpp>

namespace phaistos {
namespace definitions {

////////////////////////////////////////////////////////////
////////////////////////// Atom ////////////////////////////
////////////////////////////////////////////////////////////


//@{
//! Atom weights
const static double atom_h_weight = 1.0079;
const static double atom_c_weight = 12.0107;
const static double atom_n_weight = 14.0067;
const static double atom_o_weight = 15.9994;
const static double atom_s_weight = 32.0655;
const static double atom_p_weight = 30.9738;
//@}

//Added by MJ: options for P, O1P, O2P, O3P
//! Atom type enumeration
//! PS: pseudo atom placed at the sidechain centre of mass
//! H1, H2, H3: 	N-terminal
//! OXT:         C-terminal
//! XX, XX_S, XX_O, XX_N, XX_C, XX_H: Wild card atom types
enum AtomEnum                                      {N=0,   CA,   C,   O,   CB,   SG,   OG,   OG1,   CG,   CG1,   CG2,   SD,   OD1,   OD2,   ND1,   ND2,   CD,   CD1,   CD2,   OE1,   OE2,   NE,   NE1,   NE2,   CE,   CE1,   CE2,   CE3,   NZ,   CZ,   CZ2,   CZ3,   OH,   NH1,   NH2,   CH2,   PS, P, O1P, O2P, O3P, H,   HA,   HA1,   HA2,   HA3,   HB,   HB1,   HB2,   HB3,   HG,   HG1,   HG2,   HG3,   HG11,   HG12,   HG13,   HG21,   HG22,   HG23,   HD1,   HD2,   HD3,   HD11,   HD12,   HD13,   HD21,   HD22,   HD23,   HE,   HE1,   HE2,   HE3,   HE21,   HE22,   HH,   HH2,   HH11,   HH12,   HH21,   HH22,   HZ,   HZ1,   HZ2,   HZ3,   H1,   H2,   H3,   OXT,   XX,   XS,   XO,   XN,   XC,   XH, ATOM_ENUM_SIZE};

//! Atom names
const static char *atom_name[ATOM_ENUM_SIZE]  =    {"N", "CA",  "C", "O", "CB", "SG", "OG", "OG1", "CG", "CG1", "CG2", "SD", "OD1", "OD2", "ND1", "ND2", "CD", "CD1", "CD2","OE1", "OE2", "NE", "NE1", "NE2", "CE", "CE1", "CE2", "CE3", "NZ", "CZ", "CZ2", "CZ3", "OH", "NH1", "NH2", "CH2",  "CM", "P",  "O1P",  "O2P",  "O3P", "H", "HA", "HA1", "HA2", "HA3", "HB", "HB1", "HB2", "HB3", "HG", "HG1", "HG2", "HG3", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23", "HD1", "HD2", "HD3", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "HE", "HE1", "HE2", "HE3", "HE21", "HE22", "HH", "HH2", "HH11", "HH12", "HH21", "HH22", "HZ", "HZ1", "HZ2", "HZ3", "H1", "H2", "H3", "OXT", "XX", "XS", "XO", "XN", "XC", "XH"};

//! Atom type wildcards - Sulfur
const static bool atom_type_XS[ATOM_ENUM_SIZE] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};

//! Atom type wildcards - Oxygen
const static bool atom_type_XO[ATOM_ENUM_SIZE] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0};

//! Atom type wildcards - Nitrogen
const static bool atom_type_XN[ATOM_ENUM_SIZE] =  {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};

//! Atom type wildcards - Carbon
const static bool atom_type_XC[ATOM_ENUM_SIZE] =  {0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};

//! Atom type wildcards - Hydrogen
const static bool atom_type_XH[ATOM_ENUM_SIZE] =  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1};

//! Atom type wildcards - Phosphorus
const static bool atom_type_XP[ATOM_ENUM_SIZE] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//! Atom type wildcards - Any
const static bool atom_type_X[ATOM_ENUM_SIZE]  =  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};

//! Size of atom_type_wildcards
const static int atom_type_wildcards_size = 46;

//! Atoms covered by each type of wildcard
const static AtomEnum atom_type_wildcards[6][atom_type_wildcards_size] = {
    { SG,  SD,  XS,  XS,  XS,  XS,  XS,  XS,  XS,  XS,  XS,  XS,   XS,   XS,   XS,   XS,   XS,   XS,  XS,  XS,  XS,   XS,   XS,   XS,   XS,   XS,   XS, XS,  XS,  XS,  XS,   XS,   XS, XS,  XS,   XS,   XS,   XS,   XS, XS,  XS,  XS,  XS, XS, XS, XS},
    {  O,  OG, OG1, OD1, OD2, OE1, OE2,  OH, OXT,  XO,  XO,  XO,   XO,   XO,   XO,   XO,   XO,   XO,  XO,  XO,  XO,   XO,   XO,   XO,   XO,   XO,   XO, XO,  XO,  XO,  XO,   XO,   XO, XO,  XO,   XO,   XO,   XO,   XO, XO,  XO,  XO,  XO, XO, XO, XO},
    {  N, ND1, ND2,  NE, NE1, NE2,  NZ, NH1, NH2,  XN,  XN,  XN,   XN,   XN,   XN,   XN,   XN,   XN,  XN,  XN,  XN,   XN,   XN,   XN,   XN,   XN,   XN, XN,  XN,  XN,  XN,   XN,   XN, XN,  XN,   XN,   XN,   XN,   XN, XN,  XN,  XN,  XN, XN, XN, XN},
    { CA,   C,  CB,  CG, CG1, CG2,  CD, CD1, CD2,  CE, CE1, CE2,  CE3,   CZ,  CZ2,  CZ3,  CH2,   XC,  XC,  XC,  XC,   XC,   XC,   XC,   XC,   XC,   XC, XC,  XC,  XC,  XC,   XC,   XC, XC,  XC,   XC,   XC,   XC,   XC, XC,  XC,  XC,  XC, XC, XC, XC},
    {  H,  HA, HA2, HA3,  HB, HB1, HB2, HB3,  HG, HG1, HG2, HG3, HG11, HG12, HG13, HG21, HG22, HG23, HD1, HD2, HD3, HD11, HD12, HD13, HD21, HD22, HD23, HE, HE1, HE2, HE3, HE21, HE22, HH, HH2, HH11, HH12, HH21, HH22, HZ, HZ1, HZ2, HZ3, H1, H2, H3}};

// Added by MJ: map P to something else?
//! Translation table to translate from different naming conventions
const std::map<std::string, std::map<std::string,std::string> > atom_name_translation_map = 
     boost::assign::map_list_of<std::string, std::map<std::string, std::string> >
     // General substitutions
     ("*", boost::assign::map_list_of
      ("HN" , "H")
      ("HT1", "H1")
      ("HT2", "H2")
      ("HT3", "H3")
      ("OT1", "O")
      ("OT2", "OXT"))
     // Residue specific substitutions
     ("ASP", boost::assign::map_list_of
      ("HB1", "HB3"))
     ("CYS", boost::assign::map_list_of
      ("HG1", "HG"))
     ("GLU", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HG1", "HG3"))
     ("PHE", boost::assign::map_list_of
      ("HB1", "HB3"))
     ("GLY", boost::assign::map_list_of
      ("HA1", "HA3"))
     ("HIS", boost::assign::map_list_of
      ("HB1", "HB3"))
     ("ILE", boost::assign::map_list_of
      ("CD", "CD1")
      ("HD1", "HD11")
      ("HD2", "HD12")
      ("HD3", "HD13")
      ("HG11", "HG13"))
     ("LYS", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HD1", "HD3")
      ("HE1", "HE3")
      ("HG1", "HG3"))
     ("LEU", boost::assign::map_list_of
      ("HB1", "HB3"))
     ("MET", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HG1", "HG3"))
     ("ASN", boost::assign::map_list_of
      ("HB1", "HB3"))
     ("PRO", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HD1", "HD3")
      ("HG1", "HG3"))
     ("GLN", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HG1", "HG3"))
     ("ARG", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HD1", "HD3")
      ("HG1", "HG3"))
     ("SER", boost::assign::map_list_of
      ("HB1", "HB3")
      ("HG1", "HG"))
     ("TRP", boost::assign::map_list_of
      ("HB1", "HB3"))
     ("TYR", boost::assign::map_list_of
      ("HB1", "HB3"))
     ;

//! Check if atom is a sulphur
//!
//! \param atom_type Type of atom
//! \return True if atom is of type sulfur
inline bool is_atom_XS(AtomEnum atom_type) {
     return atom_type_XS[atom_type];
}

//! Check if atom is a oxygen
//!
//! \param atom_type Type of atom
//! \return True if atom is of type oxygen
inline bool is_atom_XO(AtomEnum atom_type) {
     return atom_type_XO[atom_type];
}

//! Check if atom is a nitrogen
//!
//! \param atom_type Type of atom
//! \return True if atom is of type nitrogen
inline bool is_atom_XN(AtomEnum atom_type) {
     return atom_type_XN[atom_type];
}

//! Check if atom is a carbon
//!
//! \param atom_type Type of atom
//! \return True if atom is of type carbon
inline bool is_atom_XC(AtomEnum atom_type) {
     return atom_type_XC[atom_type];
}

//! Check if atom is a hydrogen
//!
//! \param atom_type Type of atom
//! \return True if atom is of type hydrogen
inline bool is_atom_XH(AtomEnum atom_type) {
     return atom_type_XH[atom_type];
}

//Added by MJ
//! Check if atom is a phosphorus
//!
//! \param atom_type Type of atom
//! \return True if atom is of type phosphorus
inline bool is_atom_XP(AtomEnum atom_type) {
    return atom_type_XP[atom_type];
}

//! Check if atom is a wildcard atom
//!
//! \param atom_type Type of atom
//! \return True if atom is a wildcard type
inline bool is_atom_wildcard(AtomEnum atom_type) {
     return atom_type_X[atom_type];
}

//! Determine whether atom_type corresponds to a backbone atom
//!
//! \param atom_type Type of atom
//! \return True if atom is a backbone type
inline bool is_backbone_atom_type(AtomEnum atom_type) {
     return (atom_type == N || atom_type == CA || atom_type == C);
}

//! Determine whether atom_type corresponds to a sidechain atom
//!
//! \param atom_type Type of atom
//! \return True if atom is a sidechain type
inline bool is_sidechain_atom_type(AtomEnum atom_type) {
     return (atom_type >= CB && atom_type <= HZ3 && atom_type != H);
}

//! Overload output operator for AtomEnum
inline std::ostream &operator<<(std::ostream &o, AtomEnum a) {
     o<<atom_name[a];
     return o;
}

//! Construct (string, enum) tree for quick atom name lookup
//!
//! \return map with (name,AtomEnum) pairs
inline std::map<std::string, AtomEnum> atom_str_to_name_construct() {
     std::map<std::string, AtomEnum> atom_str_to_name;
     for (int i=0; i< ATOM_ENUM_SIZE; i++) {
          atom_str_to_name.insert(std::make_pair(atom_name[i], AtomEnum(i)));
     }
     return atom_str_to_name;
}
const static std::map<std::string, AtomEnum> atom_str_to_name = atom_str_to_name_construct();

//! Convert a string to the corresponding AtomEnum
//!
//! \param str_name Name of atom
//! \param residue_name Optionally specify residue name
//! \return atom enum
inline AtomEnum string_to_atom(std::string str_name, std::string residue_name="") {

     // Translate non-standard into standard names (e.g. 1HG2 into HG21)
     int lastIndex = str_name.size()-1;
     if ((str_name[1] == 'H') && (std::isdigit(str_name[0])))  {
          str_name = str_name.substr(1, lastIndex) + str_name[0];
     }

     // Apply translation - any residue (*)
     str_name = map_lookup(map_lookup(atom_name_translation_map, std::string("*"), std::map<std::string, std::string>()),
                           str_name, str_name);
     if (residue_name != "") {
          // Apply translation - specific residue
          str_name = map_lookup(map_lookup(atom_name_translation_map, residue_name, std::map<std::string, std::string>()),
                                str_name, str_name);
     }
                           
     // Lookup value
     AtomEnum atom_type = map_lookup(atom_str_to_name, str_name, XX);

     if (atom_type==XX) {
          std::cerr << "Unknown atom type: " << str_name << ". Using XX.\n";
     }
     
     return atom_type;
}

//! Overload input operator for AtomEnum
inline std::istream &operator>>(std::istream &input, AtomEnum &a) {
     std::string input_str = "";
     while (isalnum(input.peek())) {
          input_str += input.peek();
          input.ignore(1);
     }
     if (input_str.size() == 0) {
          std::cout << "???: " << boost::lexical_cast<std::string>((char)input.peek()) << "\n";
          input >> input_str;
          std::cout << input_str << "\n";
     }
     a = string_to_atom(input_str);
     return input;
}


////////////////////////////////////////////////////////////
//////////////////////// Residue ///////////////////////////
////////////////////////////////////////////////////////////

//Added by MJ: SEP: Ser-PO3, TPO: Thr-PO3, PTR: Tyr-PO3 residues
//! Residue type enumeration
enum ResidueEnum {ALA=0,CYS,ASP,GLU,PHE,GLY,
		          HIS,ILE,LYS,LEU,MET,ASN,
		          PRO,GLN,ARG,SER,THR,VAL,
		          TRP,TYR,SEP,TPO,PTR,AA_UNDEF,

		          // Special residue names from ASTRAL RAF
		          _2AS,  _3AH,  _5HP,  ACL,  AGM,  AIB,  ALM,  ALO,  ALY,  ARM,  ASA,  ASB,  ASK,  ASL,  ASQ,  ASX,  AYA,  BCS,  BHD,  BMT,  BNN,  BUC,  BUG,  C5C,  C6C,  CCS,  CEA,  CGU,  CHG,  CLE,  CME,  CSD,  CSO,  CSP,  CSS,  CSW,  CSX,  CXM,  CY1,  CY3,  CYG,  CYM,  CYQ,  DAH,  DAL,  DAR,  DAS,  DCY,  DGL,  DGN,  DHA,  DHI,  DIL,  DIV,  DLE,  DLY,  DNP,  DPN,  DPR,  DSN,  DSP,  DTH,  DTR,  DTY,  DVA,  EFC,  FLA,  FME,  GGL,  GL3,  GLX,  GLZ,  GMA,  GSC,  HAC,  HAR,  HIC,  HIP,  HMR,  HPQ,  HTR,  HYP,  IIL,  IYR,  KCX,  LLP,  LLY,  LTR,  LYM,  LYZ,  MAA,  MEN,  MHS,  MIS,  MLE,  MPQ,  MSA,  MSE,  MVA,  NEM,  NEP,  NLE,  NLN,  NLP,  NMC,  OAS,  OCS,  OMT,  PAQ,  PCA,  PEC,  PHI,  PHL,  PR3,  PRR,    SAC,  SAR,  SCH,  SCS,  SCY,  SEL,  SET,  SHC,  SHR,  SMC,  SOC,  STY,  SVA,  TIH,  TPL,  TPQ,  TRG,  TRO,  TYB,  TYQ,  TYS,  TYY,  UNK, RESIDUE_ENUM_SIZE};


//! Three-letter residue names
const static char *residue_name[RESIDUE_ENUM_SIZE]={"ALA","CYS","ASP","GLU","PHE","GLY",
                                                    "HIS","ILE","LYS","LEU","MET","ASN",
                                                    "PRO","GLN","ARG","SER","THR","VAL",
                                                    "TRP","TYR","SEP","TPO","PTR","???",

            // Special residue names from ASTRAL RAF
            "2AS",  "3AH",  "5HP",  "ACL",  "AGM",  "AIB",  "ALM",  "ALO",  "ALY",  "ARM",  "ASA",  "ASB",  "ASK",  "ASL",  "ASQ",  "ASX",  "AYA",  "BCS",  "BHD",  "BMT",  "BNN",  "BUC",  "BUG",  "C5C",  "C6C",  "CCS",  "CEA",  "CGU",  "CHG",  "CLE",  "CME",  "CSD",  "CSO",  "CSP",  "CSS",  "CSW",  "CSX",  "CXM",  "CY1",  "CY3",  "CYG",  "CYM",  "CYQ",  "DAH",  "DAL",  "DAR",  "DAS",  "DCY",  "DGL",  "DGN",  "DHA",  "DHI",  "DIL",  "DIV",  "DLE",  "DLY",  "DNP",  "DPN",  "DPR",  "DSN",  "DSP",  "DTH",  "DTR",  "DTY",  "DVA",  "EFC",  "FLA",  "FME",  "GGL",  "GL3",  "GLX",  "GLZ",  "GMA",  "GSC",  "HAC",  "HAR",  "HIC",  "HIP",  "HMR",  "HPQ",  "HTR",  "HYP",  "IIL",  "IYR",  "KCX",  "LLP",  "LLY",  "LTR",  "LYM",  "LYZ",  "MAA",  "MEN",  "MHS",  "MIS",  "MLE",  "MPQ",  "MSA",  "MSE",  "MVA",  "NEM",  "NEP",  "NLE",  "NLN",  "NLP",  "NMC",  "OAS",  "OCS",  "OMT",  "PAQ",  "PCA",  "PEC",  "PHI",  "PHL",  "PR3",  "PRR",  "SAC",  "SAR",  "SCH",  "SCS",  "SCY",  "SEL",  "SET",  "SHC",  "SHR",  "SMC",  "SOC",  "STY",  "SVA",  "TIH",  "TPL",   "TPQ",  "TRG",  "TRO",  "TYB",  "TYQ",  "TYS",  "TYY",  "UNK"};

//Added by MJ: S corresponding to SEP enum - maybe causes problems???
//! One-letter residue names
const static char *residue_name_short[RESIDUE_ENUM_SIZE] = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","S","T","Y","X",

            // Special residue names from ASTRAL RAF
            "D",  "H",  "E",  "R",  "R",  "A",  "A",  "T",  "K",  "R",  "D",  "D",  "D",  "D",  "D",  "B",  "A",  "C",  "D",  "T",  "A",  "C",  "L",  "C",  "C",  "C",  "C",  "E",  "A",  "L",  "C",  "A",  "C",  "C",  "C",  "C",  "C",  "M",  "C",  "C",  "C",  "C",  "C",  "F",  "A",  "R",  "D",  "C",  "E",  "Q",  "A",  "H",  "I",  "V",  "L",  "K",  "A",  "F",  "P",  "S",  "D",  "T",  "W",  "Y",  "V",  "C",  "A",  "M",  "E",  "G",  "Z",  "G",  "E",  "G",  "A",  "R",  "H",  "H",  "R",  "F",  "W",  "P",  "I",  "Y",  "K",  "K",  "K",  "W",  "K",  "K",  "A",  "N",  "H",  "S",  "L",  "G",  "G",  "M",  "V",  "H",  "H",  "L",  "L",  "L",  "G",  "S",  "C",  "M",  "Y",  "E",  "C",  "F",  "F",  "C",  "A", "S",  "G",  "C",  "C",  "C",  "S",  "S",    "C",  "K",  "C",  "C",  "Y",  "S",  "A",  "W",  "A",  "K",  "W",  "Y",  "Y",  "Y",  "Y",  "X"};


//! Translate index to short residue names
//!
//! \param index Residue index (ResidueEnum)
//! \return Residue as char
inline char index_to_residue_short(int index) {
	return *(residue_name_short[index]);
}


//! Overload output operator for ResidueEnum (using three-letter code)
inline std::ostream &operator<<(std::ostream &o, ResidueEnum a) {
     o<<residue_name[a];
     return o;
}

//Added by MJ: SEP map
//! String to three-letter-code residue conversion
static const class StrToAa {

     //! Internal map
     std::map<std::string, int> str_to_aa_map;
public:

     //! Constructor
     StrToAa() {
          str_to_aa_map.insert(std::make_pair("ALA",ALA));
          str_to_aa_map.insert(std::make_pair("CYS",CYS));
          str_to_aa_map.insert(std::make_pair("ASP",ASP));
          str_to_aa_map.insert(std::make_pair("GLU",GLU));
          str_to_aa_map.insert(std::make_pair("PHE",PHE));
          str_to_aa_map.insert(std::make_pair("GLY",GLY));
          str_to_aa_map.insert(std::make_pair("HIS",HIS));
          str_to_aa_map.insert(std::make_pair("ILE",ILE));
          str_to_aa_map.insert(std::make_pair("LYS",LYS));
          str_to_aa_map.insert(std::make_pair("LEU",LEU));
          str_to_aa_map.insert(std::make_pair("MET",MET));
          str_to_aa_map.insert(std::make_pair("ASN",ASN));
          str_to_aa_map.insert(std::make_pair("PRO",PRO));
          str_to_aa_map.insert(std::make_pair("GLN",GLN));
          str_to_aa_map.insert(std::make_pair("ARG",ARG));
          str_to_aa_map.insert(std::make_pair("SER",SER));
          str_to_aa_map.insert(std::make_pair("THR",THR));
          str_to_aa_map.insert(std::make_pair("VAL",VAL));
          str_to_aa_map.insert(std::make_pair("TRP",TRP));
          str_to_aa_map.insert(std::make_pair("TYR",TYR));
          str_to_aa_map.insert(std::make_pair("A",ALA));
          str_to_aa_map.insert(std::make_pair("C",CYS));
          str_to_aa_map.insert(std::make_pair("D",ASP));
          str_to_aa_map.insert(std::make_pair("E",GLU));
          str_to_aa_map.insert(std::make_pair("F",PHE));
          str_to_aa_map.insert(std::make_pair("G",GLY));
          str_to_aa_map.insert(std::make_pair("H",HIS));
          str_to_aa_map.insert(std::make_pair("I",ILE));
          str_to_aa_map.insert(std::make_pair("K",LYS));
          str_to_aa_map.insert(std::make_pair("L",LEU));
          str_to_aa_map.insert(std::make_pair("M",MET));
          str_to_aa_map.insert(std::make_pair("N",ASN));
          str_to_aa_map.insert(std::make_pair("P",PRO));
          str_to_aa_map.insert(std::make_pair("Q",GLN));
          str_to_aa_map.insert(std::make_pair("R",ARG));
          str_to_aa_map.insert(std::make_pair("S",SER));
          str_to_aa_map.insert(std::make_pair("T",THR));
          str_to_aa_map.insert(std::make_pair("V",VAL));
          str_to_aa_map.insert(std::make_pair("W",TRP));
          str_to_aa_map.insert(std::make_pair("Y",TYR));

          // Protonation states of HIS in CHARMM format - read in as HIS
          str_to_aa_map.insert(std::make_pair("HSD",HIS));
          str_to_aa_map.insert(std::make_pair("HSE",HIS));

          // Special cases taken from ASTRAL RAF
          str_to_aa_map.insert(std::make_pair("2AS",_2AS));
          str_to_aa_map.insert(std::make_pair("3AH",_3AH));
          str_to_aa_map.insert(std::make_pair("5HP",_5HP));
          str_to_aa_map.insert(std::make_pair("ACL",ACL));
          str_to_aa_map.insert(std::make_pair("AGM",AGM));
          str_to_aa_map.insert(std::make_pair("AIB",AIB));
          str_to_aa_map.insert(std::make_pair("ALM",ALM));
          str_to_aa_map.insert(std::make_pair("ALO",ALO));
          str_to_aa_map.insert(std::make_pair("ALY",ALY));
          str_to_aa_map.insert(std::make_pair("ARM",ARM));
          str_to_aa_map.insert(std::make_pair("ASA",ASA));
          str_to_aa_map.insert(std::make_pair("ASB",ASB));
          str_to_aa_map.insert(std::make_pair("ASK",ASK));
          str_to_aa_map.insert(std::make_pair("ASL",ASL));
          str_to_aa_map.insert(std::make_pair("ASQ",ASQ));
          str_to_aa_map.insert(std::make_pair("ASX",ASX));
          str_to_aa_map.insert(std::make_pair("AYA",AYA));
          str_to_aa_map.insert(std::make_pair("BCS",BCS));
          str_to_aa_map.insert(std::make_pair("BHD",BHD));
          str_to_aa_map.insert(std::make_pair("BMT",BMT));
          str_to_aa_map.insert(std::make_pair("BNN",BNN));
          str_to_aa_map.insert(std::make_pair("BUC",BUC));
          str_to_aa_map.insert(std::make_pair("BUG",BUG));
          str_to_aa_map.insert(std::make_pair("C5C",C5C));
          str_to_aa_map.insert(std::make_pair("C6C",C6C));
          str_to_aa_map.insert(std::make_pair("CCS",CCS));
          str_to_aa_map.insert(std::make_pair("CEA",CEA));
          str_to_aa_map.insert(std::make_pair("CGU",CGU));
          str_to_aa_map.insert(std::make_pair("CHG",CHG));
          str_to_aa_map.insert(std::make_pair("CLE",CLE));
          str_to_aa_map.insert(std::make_pair("CME",CME));
          str_to_aa_map.insert(std::make_pair("CSD",CSD));
          str_to_aa_map.insert(std::make_pair("CSO",CSO));
          str_to_aa_map.insert(std::make_pair("CSP",CSP));
          str_to_aa_map.insert(std::make_pair("CSS",CSS));
          str_to_aa_map.insert(std::make_pair("CSW",CSW));
          str_to_aa_map.insert(std::make_pair("CSX",CSX));
          str_to_aa_map.insert(std::make_pair("CXM",CXM));
          str_to_aa_map.insert(std::make_pair("CY1",CY1));
          str_to_aa_map.insert(std::make_pair("CY3",CY3));
          str_to_aa_map.insert(std::make_pair("CYG",CYG));
          str_to_aa_map.insert(std::make_pair("CYM",CYM));
          str_to_aa_map.insert(std::make_pair("CYQ",CYQ));
          str_to_aa_map.insert(std::make_pair("DAH",DAH));
          str_to_aa_map.insert(std::make_pair("DAL",DAL));
          str_to_aa_map.insert(std::make_pair("DAR",DAR));
          str_to_aa_map.insert(std::make_pair("DAS",DAS));
          str_to_aa_map.insert(std::make_pair("DCY",DCY));
          str_to_aa_map.insert(std::make_pair("DGL",DGL));
          str_to_aa_map.insert(std::make_pair("DGN",DGN));
          str_to_aa_map.insert(std::make_pair("DHA",DHA));
          str_to_aa_map.insert(std::make_pair("DHI",DHI));
          str_to_aa_map.insert(std::make_pair("DIL",DIL));
          str_to_aa_map.insert(std::make_pair("DIV",DIV));
          str_to_aa_map.insert(std::make_pair("DLE",DLE));
          str_to_aa_map.insert(std::make_pair("DLY",DLY));
          str_to_aa_map.insert(std::make_pair("DNP",DNP));
          str_to_aa_map.insert(std::make_pair("DPN",DPN));
          str_to_aa_map.insert(std::make_pair("DPR",DPR));
          str_to_aa_map.insert(std::make_pair("DSN",DSN));
          str_to_aa_map.insert(std::make_pair("DSP",DSP));
          str_to_aa_map.insert(std::make_pair("DTH",DTH));
          str_to_aa_map.insert(std::make_pair("DTR",DTR));
          str_to_aa_map.insert(std::make_pair("DTY",DTY));
          str_to_aa_map.insert(std::make_pair("DVA",DVA));
          str_to_aa_map.insert(std::make_pair("EFC",EFC));
          str_to_aa_map.insert(std::make_pair("FLA",FLA));
          str_to_aa_map.insert(std::make_pair("FME",FME));
          str_to_aa_map.insert(std::make_pair("GGL",GGL));
          str_to_aa_map.insert(std::make_pair("GL3",GL3));
          str_to_aa_map.insert(std::make_pair("GLX",GLX));
          str_to_aa_map.insert(std::make_pair("GLZ",GLZ));
          str_to_aa_map.insert(std::make_pair("GMA",GMA));
          str_to_aa_map.insert(std::make_pair("GSC",GSC));
          str_to_aa_map.insert(std::make_pair("HAC",HAC));
          str_to_aa_map.insert(std::make_pair("HAR",HAR));
          str_to_aa_map.insert(std::make_pair("HIC",HIC));
          str_to_aa_map.insert(std::make_pair("HIP",HIP));
          str_to_aa_map.insert(std::make_pair("HMR",HMR));
          str_to_aa_map.insert(std::make_pair("HPQ",HPQ));
          str_to_aa_map.insert(std::make_pair("HTR",HTR));
          str_to_aa_map.insert(std::make_pair("HYP",HYP));
          str_to_aa_map.insert(std::make_pair("IIL",IIL));
          str_to_aa_map.insert(std::make_pair("IYR",IYR));
          str_to_aa_map.insert(std::make_pair("KCX",KCX));
          str_to_aa_map.insert(std::make_pair("LLP",LLP));
          str_to_aa_map.insert(std::make_pair("LLY",LLY));
          str_to_aa_map.insert(std::make_pair("LTR",LTR));
          str_to_aa_map.insert(std::make_pair("LYM",LYM));
          str_to_aa_map.insert(std::make_pair("LYZ",LYZ));
          str_to_aa_map.insert(std::make_pair("MAA",MAA));
          str_to_aa_map.insert(std::make_pair("MEN",MEN));
          str_to_aa_map.insert(std::make_pair("MHS",MHS));
          str_to_aa_map.insert(std::make_pair("MIS",MIS));
          str_to_aa_map.insert(std::make_pair("MLE",MLE));
          str_to_aa_map.insert(std::make_pair("MPQ",MPQ));
          str_to_aa_map.insert(std::make_pair("MSA",MSA));
          str_to_aa_map.insert(std::make_pair("MSE",MSE));
          str_to_aa_map.insert(std::make_pair("MVA",MVA));
          str_to_aa_map.insert(std::make_pair("NEM",NEM));
          str_to_aa_map.insert(std::make_pair("NEP",NEP));
          str_to_aa_map.insert(std::make_pair("NLE",NLE));
          str_to_aa_map.insert(std::make_pair("NLN",NLN));
          str_to_aa_map.insert(std::make_pair("NLP",NLP));
          str_to_aa_map.insert(std::make_pair("NMC",NMC));
          str_to_aa_map.insert(std::make_pair("OAS",OAS));
          str_to_aa_map.insert(std::make_pair("OCS",OCS));
          str_to_aa_map.insert(std::make_pair("OMT",OMT));
          str_to_aa_map.insert(std::make_pair("PAQ",PAQ));
          str_to_aa_map.insert(std::make_pair("PCA",PCA));
          str_to_aa_map.insert(std::make_pair("PEC",PEC));
          str_to_aa_map.insert(std::make_pair("PHI",PHI));
          str_to_aa_map.insert(std::make_pair("PHL",PHL));
          str_to_aa_map.insert(std::make_pair("PR3",PR3));
          str_to_aa_map.insert(std::make_pair("PRR",PRR));
          str_to_aa_map.insert(std::make_pair("PTR",PTR)); //phosphorylated TYR
          str_to_aa_map.insert(std::make_pair("SAC",SAC));
          str_to_aa_map.insert(std::make_pair("SAR",SAR));
          str_to_aa_map.insert(std::make_pair("SCH",SCH));
          str_to_aa_map.insert(std::make_pair("SCS",SCS));
          str_to_aa_map.insert(std::make_pair("SCY",SCY));
          str_to_aa_map.insert(std::make_pair("SEL",SEL));
          str_to_aa_map.insert(std::make_pair("SEP",SEP)); //phosphorylated SER
          str_to_aa_map.insert(std::make_pair("SET",SET));
          str_to_aa_map.insert(std::make_pair("SHC",SHC));
          str_to_aa_map.insert(std::make_pair("SHR",SHR));
          str_to_aa_map.insert(std::make_pair("SMC",SMC));
          str_to_aa_map.insert(std::make_pair("SOC",SOC));
          str_to_aa_map.insert(std::make_pair("STY",STY));
          str_to_aa_map.insert(std::make_pair("SVA",SVA));
          str_to_aa_map.insert(std::make_pair("TIH",TIH));
          str_to_aa_map.insert(std::make_pair("TPL",TPL));
          str_to_aa_map.insert(std::make_pair("TPO",TPO)); //phosphorylated THR
          str_to_aa_map.insert(std::make_pair("TPQ",TPQ));
          str_to_aa_map.insert(std::make_pair("TRG",TRG));
          str_to_aa_map.insert(std::make_pair("TRO",TRO));
          str_to_aa_map.insert(std::make_pair("TYB",TYB));
          str_to_aa_map.insert(std::make_pair("TYQ",TYQ));
          str_to_aa_map.insert(std::make_pair("TYS",TYS));
          str_to_aa_map.insert(std::make_pair("TYY",TYY));
          str_to_aa_map.insert(std::make_pair("UNK",UNK));
     }

     //! Overload [] indexing operator (using key string)
     int operator[](std::string key) const {

          std::map<std::string,int>::const_iterator it = str_to_aa_map.find(boost::to_upper_copy(key));
          if (it != str_to_aa_map.end()) {
               return it->second;
          } else {
               std::cerr << "WARNING: Unknown amino acid type: " << key << "\n";
               return -1;
          }          
     }

     //! Overload () indexing operator
     int operator()(std::string key) const {
          return operator[](key);
     }
} str_to_aa;


//! Traversal enum - used by get_neighbour
enum IterateEnum {BACKBONE=0, CA_ONLY, SC_ONLY, ALL, POSITIONING, ITERATE_ENUM_SIZE};
static const std::string iterateEnumNames[] = {"BACKBONE", "CA_ONLY", "SC_ONLY", "ALL", "POSITIONING"};

//! Input IterateEnum from string
inline std::istream &operator>>(std::istream &input, IterateEnum &ie) {
     std::string raw_string;
     input >> raw_string;

     for (unsigned int i=0; i<ITERATE_ENUM_SIZE; ++i) {
          if (raw_string == iterateEnumNames[i]) {
               ie = IterateEnum(i);
          }
     }     
     return input;
}
//! Output IterateEnum
inline std::ostream &operator<<(std::ostream &o, const IterateEnum &ie) {
     o << iterateEnumNames[static_cast<unsigned int>(ie)];
     return o;
}


//! Terminal enum - indicates whether a residue is at the C-terminal, N-terminal, or internal in the protein
enum TerminalEnum {INTERNAL=0, CTERM, NTERM};




//! Atom selection type
enum AtomSelectionEnum {
    NO_ATOMS=1,  // nonsense: here just for initializations in loops
    BACKBONE_ATOMS=2,
    BACKBONE_O_ATOMS=4,
    BACKBONE_H_ATOMS=8,
    CB_ATOMS=16,
    SIDECHAIN_ATOMS=32,
    NON_BACKBONE_H_ATOMS=64,
    PSEUDO_SIDECHAIN_ATOMS=128,
    // Add non 2-powers last - otherwise the string conversion is off
    ALL_PHYSICAL_ATOMS=127,
    ALL_ATOMS=255,
    ATOM_SELECTION_ENUM_SIZE = 10,
    };

//! Atom selection array (used for output)
const static AtomSelectionEnum atom_selection_array[ATOM_SELECTION_ENUM_SIZE]={
     NO_ATOMS,
     BACKBONE_ATOMS,
     BACKBONE_O_ATOMS,
     BACKBONE_H_ATOMS,
     CB_ATOMS,
     SIDECHAIN_ATOMS,
     NON_BACKBONE_H_ATOMS,
     ALL_PHYSICAL_ATOMS,
     PSEUDO_SIDECHAIN_ATOMS,
     ALL_ATOMS
};     

//! Atom selection names
const static char *atom_selection_name[ATOM_SELECTION_ENUM_SIZE]={
     "NO_ATOMS",
     "BACKBONE_ATOMS",
     "BACKBONE_O_ATOMS",
     "BACKBONE_H_ATOMS",
     "CB_ATOMS",
     "SIDECHAIN_ATOMS",
     "NON_BACKBONE_H_ATOMS",
     "ALL_PHYSICAL_ATOMS",
     "PSEUDO_SIDECHAIN_ATOMS",
     "ALL_ATOMS"};

// Overload output operator for atomSelectionAtom
inline std::ostream &operator<<(std::ostream &o, AtomSelectionEnum a) {
     int power = 1;
     for (unsigned int i=1; i<ATOM_SELECTION_ENUM_SIZE; i++) {
          power*=2;
          if (a & power)
               o<<atom_selection_name[i] << " ";
     }
     return o;
}

// namespace std {
//! Overload input operator for atomSelectionAtom
inline std::istream &operator>>(std::istream &input, AtomSelectionEnum &a) {
     std::string raw_string;
     std::getline(input, raw_string);

     std::vector<std::string> split_raw_string;
     boost::split(split_raw_string, 
                  raw_string, 
                  boost::is_any_of(" "));

     AtomSelectionEnum atom_selection = NO_ATOMS;
     for (unsigned int i=0; i<split_raw_string.size(); ++i) {
          for (unsigned int j=0; j<ATOM_SELECTION_ENUM_SIZE; ++j) {
               if (split_raw_string[i] == atom_selection_name[j]) {
                    atom_selection = AtomSelectionEnum((int)atom_selection | 
                                                       (int)atom_selection_array[j]);
                    break;
               }
          }
     }

     if (atom_selection != NO_ATOMS)
          a = atom_selection;

     return input;
}
// }

//! Overload + operator for AtomSelectionEnum
inline AtomSelectionEnum operator+(AtomSelectionEnum v1, AtomSelectionEnum v2) {
     return AtomSelectionEnum((int)v1 + (int)v2);
}

//! Overload - operator for AtomSelectionEnum
inline AtomSelectionEnum operator-(AtomSelectionEnum v1, AtomSelectionEnum v2) {
     return AtomSelectionEnum((int)v1 - (int)v2);
}

//! Overload | operator for AtomSelectionEnum
inline AtomSelectionEnum operator|(AtomSelectionEnum a, AtomSelectionEnum b) {
     return AtomSelectionEnum(int(a) | int(b));
}

//! Overload |= operator for AtomSelectionEnum
inline AtomSelectionEnum& operator|=(AtomSelectionEnum& a, AtomSelectionEnum b) {
     a = a | b;
     return a;
}

//! Overload & operator for AtomSelectionEnum
inline AtomSelectionEnum operator&(AtomSelectionEnum a, AtomSelectionEnum b) {
     return AtomSelectionEnum(int(a) & int(b));
}

//! Overload &= operator for AtomSelectionEnum
inline AtomSelectionEnum& operator&=(AtomSelectionEnum& a, AtomSelectionEnum b) {
     a = a & b;
     return a;
}


//! Move types
enum MoveTypeEnum {NON_LOCAL=1, LOCAL=2, SIDECHAIN=4};


//! Secondary structure enum
enum SecondaryStructureEnum {SS_H=0, SS_E, SS_C, SECONDARY_STRUCTURE_ENUM_SIZE};

//! Secondary structure enum names
const static char *ss_name[SECONDARY_STRUCTURE_ENUM_SIZE] = {"H","E","C"};

//! Translate index to secondary structure names
inline char index_to_ss(int index) {
	return *(ss_name[index]);
}


//! String to secondary structure conversion
static const class StrToSS {
     std::map<std::string, int> str_to_ss_map;
public:
     StrToSS() {
          str_to_ss_map.insert(std::make_pair("H",SS_H));
          str_to_ss_map.insert(std::make_pair("E",SS_E));
          str_to_ss_map.insert(std::make_pair("C",SS_C));

          // Other DSSP input types
          str_to_ss_map.insert(std::make_pair("G",SS_H));
          str_to_ss_map.insert(std::make_pair("I",SS_H));
          str_to_ss_map.insert(std::make_pair("B",SS_E));
          str_to_ss_map.insert(std::make_pair("T",SS_C));
          str_to_ss_map.insert(std::make_pair("S",SS_C));
          str_to_ss_map.insert(std::make_pair(" ",SS_C));
     }

     int operator[](std::string key) const {
          std::map<std::string,int>::const_iterator it = str_to_ss_map.find(key);
          if (it != str_to_ss_map.end()) {
               return it->second;
          } else {
               return uninitialized<int>();
          }
     }

     int operator()(std::string key) const {
          return operator[](key);
     }
} str_to_ss;


//! Angle types
enum AngleEnum{DIHEDRAL=0, ANGLE, CHI1, CHI2, CHI3, CHI4};

}

}

#endif
