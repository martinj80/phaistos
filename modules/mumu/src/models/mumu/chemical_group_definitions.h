// chemical_group_definitions.h --- Various definitions used by the MUMU energy term
// Copyright (C) 2012-2013 Kristoffer Enøe Johansson
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


#ifndef CHEMICAL_GROUP_DEFINITIONS_H
#define CHEMICAL_GROUP_DEFINITIONS_H

#include "protein/definitions.h"
#include "utils/vector_matrix_3d.h"

using namespace phaistos::definitions;

// chemical groups enum            0      1    2    3    4    5    6    7    8    9   10   11   12   13  
enum ChemicalGroupEnum {NONE=-1, NCO=0, NH3, COO, ALC,  SH,  C2,  C3, CH3, ARO, IMI, GUA, CSC,  CN,  BB,
                        CHEMICAL_GROUP_ENUM_SIZE};

// chemical group names
const char* chemical_group_name[CHEMICAL_GROUP_ENUM_SIZE] =    {"NCO","NH3","COO","ALC", "SH", "C2", "C3",
                                                                "CH3","ARO","IMI","GUA","CSC", "CN", "BB"};
// disk radius of the chemical groups
const double chemical_group_radius[CHEMICAL_GROUP_ENUM_SIZE] = { 2.5,  1.8,  2.0,  1.3,  1.4,  1.7,  2.5, 
                                                                 1.4,  3.0,  2.8,  2.5,  2.5,  1.8,  2.5};

//                         A C D E F G H I K L M N P Q R S T V W Y
int sc_group_number[20] = {1,2,2,2,2,0,2,3,3,3,2,2,2,2,3,2,3,3,3,3};

ChemicalGroupEnum sc_group[20][3] = {{ CH3, NONE, NONE}, //  0 Ala - missing Ca
                                     {  C2,   SH, NONE}, //  1 Cys
                                     {  C2,  COO, NONE}, //  2 Asp
                                     {  C3,  COO, NONE}, //  3 Glu
                                     {  C2,  ARO, NONE}, //  4 Phe
                                     {NONE, NONE, NONE}, //  5 Gly - missing Ca
                                     {  C2,  IMI, NONE}, //  6 His
                                     {  C3,  CH3,  CH3}, //  7 Ile
                                     {  C3,   C2,  NH3}, //  8 Lys
                                     {  C3,  CH3,  CH3}, //  9 Leu
                                     {  C2,  CSC, NONE}, // 10 Met
                                     {  C2,  NCO, NONE}, // 11 Asn
                                     {  C2,   C2, NONE}, // 12 Pro
                                     {  C3,  NCO, NONE}, // 13 Gln
                                     {  C2,   C2,  GUA}, // 14 Arg
                                     {  C2,  ALC, NONE}, // 15 Ser
                                     {  C2,  CH3,  ALC}, // 16 Thr
                                     {  C2,  CH3,  CH3}, // 17 Val
                                     {  C3,   CN,  ARO}, // 18 Trp
                                     {  C2,  ARO,  ALC}};// 19 Tyr

int group_atom_number[20][3] = {{1, 0, 0}, //  0 Ala - CH3
                                {2, 1, 0}, //  1 Cys -  C2, SH
                                {2, 2, 0}, //  2 Asp -  C2, COO
                                {3, 2, 0}, //  3 Glu -  C3, COO
                                {2, 3, 0}, //  4 Phe -  C2, ARO
                                {0, 0, 0}, //  5 Gly - 
                                {2, 5, 0}, //  6 His -  C2, IMI
                                {3, 1, 1}, //  7 Ile -  C3, CH3, CH3
                                {3, 2, 1}, //  8 Lys -  C3,  C2, NH3
                                {3, 1, 1}, //  9 Leu -  C3, CH3, CH3
                                {2, 3, 0}, // 10 Met -  C2, CSC
                                {2, 2, 0}, // 11 Asn -  C2, NCO
                                {2, 2, 0}, // 12 Pro -  C2,  C2
                                {3, 2, 0}, // 13 Gln -  C3, NCO
                                {2, 2, 1}, // 14 Arg -  C2,  C2, GUA
                                {2, 1, 0}, // 15 Ser -  C2, ALC
                                {2, 1, 1}, // 16 Thr -  C2, CH3, ALC
                                {2, 1, 1}, // 17 Val -  C2, CH3, CH3
                                {3, 2, 3}, // 18 Trp -  C3,  CN, ARO
                                {2, 3, 1}};// 19 Tyr -  C2, ARO, ALC

AtomEnum group_atom[20][3][6] = {{{CB}      , {XX}        , {XX}}, //  0 Ala - CH3
                                 {{CA,CB}   , {SG}        , {XX}}, //  1 Cys -  C2,  SH
                                 {{CA,CB}   , {OD1,OD2}   , {XX}}, //  2 Asp -  C2, COO
                                 {{CA,CB,CG}, {OE1,OE2}   , {XX}}, //  3 Glu -  C3, COO
                                 {{CA,CB}   , {CG,CE1,CE2}, {XX}}, //  4 Phe -  C2, ARO
                                 {{XX}      , {XX}        , {XX}}, //  5 Gly - 
                                 {{CA,CB}   , {CG,ND1,CD2,CE1,NE2}, {XX}}, //  6 His -  C2, IMI
                                 {{CA,CB,CG1}, {CG2}      ,{CD1}}, //  7 Ile -  C3, CH3, CH3
                                 {{CA,CB,CG}, {CD,CE}     , {NZ}}, //  8 Lys -  C3,  C2, NH3
                                 {{CA,CB,CG}, {CD1}       ,{CD2}}, //  9 Leu -  C3, CH3, CH3
                                 {{CA,CB}   , {CG,SD,CE}  , {XX}}, // 10 Met -  C2, CSC
                                 {{CA,CB}   , {CG,ND2}    , {XX}}, // 11 Asn -  C2, NCO
                                 {{CA,CB}   , {CG,CD}     , {XX}}, // 12 Pro -  C2,  C2
                                 {{CA,CB,CG}, {CD,NE2}    , {XX}}, // 13 Gln -  C3, NCO
                                 {{CA,CB}   , {CG,CD}     , {CZ}}, // 14 Arg -  C2,  C2, GUA
                                 {{CA,CB}   , {OG}        , {XX}}, // 15 Ser -  C2, ALC
                                 {{CA,CB}   , {CG2}       ,{OG1}}, // 16 Thr -  C2, CH3, ALC
                                 {{CA,CB}   , {CG1}       ,{CG2}}, // 17 Val -  C2, CH3, CH3
                                 {{CA,CB,CG}, {CD1,NE1} , {CD2,CZ2,CZ3}}, // 18 Trp - C3, CN, ARO
                                 {{CA,CB}   , {CG,CE1,CE2},{OH}}};// 19 Tyr -  C2, ARO, ALC

#endif
