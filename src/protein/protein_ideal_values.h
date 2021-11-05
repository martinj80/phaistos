// enghHuberConstants.h --- Bond length and bond angle constants from Engh, Huber, 1991
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

#ifndef ENGHHUBERCONSTANTS_H
#define ENGHHUBERCONSTANTS_H

#include "protein/pdb_input.h"
#include "protein/definitions.h"

namespace phaistos {

//! Bond length constants - primarily from Engh, Huber, 1991
//!
//! \param atom1 First atom
//! \param atom2 Second atom
//! \param res Residue containing atoms
//! \param std_dev Optional output variable for standard deviation
//! \return bond length
inline double bond_length_constants(enum definitions::AtomEnum atom1, enum definitions::AtomEnum atom2,
                                    enum definitions::ResidueEnum res=definitions::AA_UNDEF, double *std_dev=NULL) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Set std_dev pointer to point to local variable if not provided
     if (std_dev==NULL) {
          double std_dev_tmp;
          std_dev = &std_dev_tmp;
     }
     
     double bond_length = UNINITIALIZED;
     if ((atom1==N && atom2==CA) || (atom1==CA && atom2==N)) {
          if (res==GLY) {
               bond_length = 1.451;
               *std_dev = 0.016;
          } else if (res==PRO) {
               bond_length = 1.466;
               *std_dev = 0.015;
          } else {
               bond_length = 1.458;
               *std_dev = 0.019;
          }
     } else if ((atom1==CA && atom2==C) || (atom1==C && atom2==CA)) {
          if (res==GLY) {
               bond_length = 1.516;
               *std_dev = 0.018;
          } else {
               bond_length = 1.525;
               *std_dev = 0.021;
          }
     } else if ((atom1==C && atom2==N) || (atom1==N && atom2==C)) {
          if (res==PRO) {
               bond_length = 1.341;
               *std_dev = 0.016;
          } else {
               bond_length = 1.329;
               *std_dev = 0.014;
          }
     } else if ((atom1==CA && atom2==CB) || (atom1==CB && atom2==CA)) {
          if (res==ALA) {
               bond_length = 1.521;
               *std_dev = 0.033;
          } else if (res==ILE || res==THR || res==VAL) {
               bond_length = 1.540;
               *std_dev = 0.027;
          } else {
               bond_length = 1.530;
               *std_dev = 0.020;
          }
     } else if ((atom1==C && atom2==O) || (atom1==O && atom2==C) || (atom1==C && atom2==OXT) || (atom1==OXT && atom2==C)) {
          bond_length = 1.231;
          *std_dev = 0.020;
          
     
     // Hydrogen atoms. Not from Engh, Huber, 1991
     } else if (atom2>=H && atom2<=H3) {
          bond_length = 1.0;
          *std_dev = 0.020;
          
     // for Calpha chain
     } else if (atom1==CA && atom2==CA) {
          bond_length = 3.8;
          *std_dev = 0.020;
     }

     
     // bond lengths for sidechain bonds 
     if (!is_initialized(bond_length)) {
          switch(res) {
          case ALA:
               break;
          case CYS:                                       // Cysteine
               if ((atom1==CB && atom2==SG) || (atom1==SG && atom2==CB)) {
                    // "ca-cb-sg") {
                    // CH1E-CH2E-S
                    bond_length = 1.808;
                    *std_dev = 0.020;
               }
               break;
          case ASP:                                       // Aspartic Acid
               if ( (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    // ca-cb-cg
                    // C-CH2E-CH1E
                    bond_length = 1.516;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==OD1) || (atom1==OD1 && atom2==CG) ) {
                    // cb-cg-od
                    // CH2E-C-OC
                    bond_length = 1.249;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==OD2) || (atom1==OD2 && atom2==CG) ) {
                    // cb-cg-od
                    // CH2E-C-OC
                    bond_length = 1.249;
                    *std_dev = 0.020;
               }
               break;
          case GLU:                                        // Glutamic Acid
               if  ((atom1==CG && atom2==CB) || (atom1==CB && atom2==CG) ) {
                    // ca-cb-cg
                    // CH1E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if (  (atom1==CG && atom2==CD) || (atom1==CD && atom2==CG) ) {
                    // cb-cg-cd
                    // C-CH2E-CH2E
                    bond_length = 1.516;
                    *std_dev = 0.020;
               } else if ( (atom1==CD && atom2==OE1) || (atom1==OE1 && atom2==CD) ) {
                    // cg-cd-oe
                    // CH2E-C-OC
                    bond_length = 1.249;
                    *std_dev = 0.020;
               } else if ( (atom1==CD && atom2==OE2) || (atom1==OE2 && atom2==CD) ) {
                    // cg-cd-oe
                    // CH2E-C-OC
                    bond_length = 1.249;
                    *std_dev = 0.020;
               }
               break;
          case PHE:                                        // Phenylalanine
               if (  (atom1==CG && atom2==CB) || (atom1==CB && atom2==CG) ) {
                    // ca-cb-cg") {
                    // CF-CH2E-CH1E
                    bond_length = 1.502;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD1) || (atom1==CD1 && atom2==CG) ) {
                    // cb-cg-cd") {
                    // CH2E-CF-CR1E
                    bond_length = 1.384;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD2) || (atom1==CD2 && atom2==CG) ) {
                    // cb-cg-cd") {
                    // CH2E-CF-CR1E
                    bond_length = 1.384;
                    *std_dev = 0.020;
               } else if ( (atom1==CD1 && atom2==CE1) || (atom1==CE1 && atom2==CD1) ) {
                    // bond == "cg-cd-ce") {
                    // CF-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               } else if ( (atom1==CD2 && atom2==CE2) || (atom1==CE2 && atom2==CD2) ) {
                    // bond == "cg-cd-ce") {
                    // CF-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               } else if ( (atom1==CZ && atom2==CE1) || (atom1==CE1 && atom2==CZ) )  {
                    // bond == "cd-ce-cz") {
                    // CR1E-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               } else if ( (atom1==CZ && atom2==CE2) || (atom1==CE2 && atom2==CZ) )  {
                    // bond == "cd-ce-cz") {
                    // CR1E-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               }
               break;
          case GLY:                                        // Glycine
               break;
          case HIS:                                        // Histidine
               if ( (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // C5-CH2E-CH1E
                    bond_length = 1.497;
                    *std_dev = 0.020;
               } else if (  (atom1==CG && atom2==ND1) || (atom1==ND1 && atom2==CG)) {
                    //bond == "cb-cg-nd1") {
                    // CH2E-C5-NR
                    bond_length = 1.371;
                    *std_dev = 0.020;
               } else if (  (atom1==CG && atom2==CD2) || (atom1==CD2 && atom2==CG)) {
                    //bond == "cb-cg-cd2") {
                    // CH2E-C5-CR1E
                    bond_length = 1.356;
                    *std_dev = 0.020;
               } else if ( (atom1==ND1 && atom2==CE1) || (atom1==CE1 && atom2==ND1)) {
                    //bond == "cg-nd1-ce1") {
                    // C5-NR-CR1E
                    bond_length = 1.319;
                    *std_dev = 0.020;
               } else if (  (atom1==CD2 && atom2==NE2) || (atom1==NE2 && atom2==CD2)) {
                    //bond == "cg-cd2-ne2") {
                    // C5-CR1E-NH1
                    bond_length = 1.374;
                    *std_dev = 0.020;
               } else if (  (atom1==CE1 && atom2==NE2) || (atom1==NE2 && atom2==CE1)) {
                    //bond == "nd1-ce1-ne2") {
                    // NR-CR1E-NH1
                    bond_length = 1.374;
                    *std_dev = 0.020;
               }
               break;
          case ILE:                                        // Isoleucine
               if (  (atom1==CB && atom2==CG1) || (atom1==CG1 && atom2==CB)) {
                    //bond == "ca-cb-cg1") {
                    // CH1E-CH1E-CH2E
                    bond_length = 1.530;
                    *std_dev = 0.020;
               } else if ( (atom1==CB && atom2==CG2) || (atom1==CG2 && atom2==CB)) {
                    //bond == "ca-cb-cg2") {
                    // CH1E-CH1E-CH3E
                    bond_length = 1.521;
                    *std_dev = 0.020;
               } else if ((atom1==CG1 && atom2==CD1) || (atom1==CD1 && atom2==CG1)) {
                    //bond == "cb-cg1-cd1") {
                    // CH1E-CH2E-CH3E
                    bond_length = 1.513;
                    *std_dev = 0.020;
               }
               break;
          case LYS:                                        // Lysine
               if ( (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD) || (atom1==CD && atom2==CG)) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ((atom1==CD && atom2==CE) || (atom1==CE && atom2==CD)) {
                    //bond == "cg-cd-ce") {
                    // CH2E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ((atom1==CE && atom2==NZ) || (atom1==NZ && atom2==CE)) {
                    //bond == "cd-ce-nz") {
                    // CH2E-CH2E-NH3
                    bond_length = 1.489;
                    *std_dev = 0.020;
               }
               break;
          case LEU:                                        // Leucine
               if (  (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH1E
                    bond_length = 1.530;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD1) || (atom1==CD1 && atom2==CG)) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH1E-CH3E
                    bond_length = 1.521;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD2) || (atom1==CD2 && atom2==CG)) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH1E-CH3E
                    bond_length = 1.521;
                    *std_dev = 0.020;
               }
               break;
          case MET:                                        // Methionine
               if (  (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==SD) || (atom1==SD && atom2==CG)) {
                    //bond == "cb-cg-sd")
                    // CH2E-CH2E-SM
                    bond_length = 1.803;
                    *std_dev = 0.020;
               } else if ((atom1==SD && atom2==CE) || (atom1==CE && atom2==SD)) {
                    //bond == "cg-sd-ce") {
                    // CH2E-SM-CH3E
                    bond_length = 1.791;
                    *std_dev = 0.020;
               }
               break;
          case ASN:                                        // Asparagine
               if (  (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // C-CH2E-CH1E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==OD1) || (atom1==OD1 && atom2==CG)) {
                    //bond == "cb-cg-od") {
                    // CH2E-C-O
                    bond_length = 1.231;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==ND2) || (atom1==ND2 && atom2==CG)) {
                    //bond == "cb-cg-nd") {
                    // CH2E-C-NH2
                    bond_length = 1.328;
                    *std_dev = 0.020;
               }
               break;
          case PRO:                                        // Proline
               if (  (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2P
                    bond_length = 1.492;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD) || (atom1==CD && atom2==CG)) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH2P-CH2P
                    bond_length = 1.503;
                    *std_dev = 0.020;
               }
               break;
          case GLN:                                        // Glutamine
               if (  (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD) || (atom1==CD && atom2==CG)) {
                    //bond == "cb-cg-cd") {
                    // C-CH2E-CH2E
                    bond_length = 1.516;
                    *std_dev = 0.020;
               } else if ((atom1==CD && atom2==OE1) || (atom1==OE1 && atom2==CD)) {
                    //bond == "cg-cd-oe") {
                    // CH2E-C-O
                    bond_length = 1.231;
                    *std_dev = 0.020;
               } else if ((atom1==CD && atom2==NE2) || (atom1==NE2 && atom2==CD)) {
                    //bond == "cg-cd-ne") {
                    // CH2E-C-NH2
                    bond_length = 1.328;
                    *std_dev = 0.020;
               }
               break;
          case ARG:                                        // Arginine
               if (  (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ((atom1==CG && atom2==CD) || (atom1==CD && atom2==CG)) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH2E-CH2E
                    bond_length = 1.520;
                    *std_dev = 0.020;
               } else if ( (atom1==CD && atom2==NE) || (atom1==CD && atom2==NE)) {
                    //bond == "cg-cd-ne") {
                    // CH2E-CH2E-NH1
                    bond_length = 1.460;
                    *std_dev = 0.020;
               } else if ((atom1==NE && atom2==CZ) || (atom1==CZ && atom2==NE)) {
                    //bond == "cd-ne-cz") {
                    // C-NH1-CH2E
                    bond_length = 1.329;
                    *std_dev = 0.020;
               } else if ((atom1==CZ && atom2==NH1) || (atom1==NH1 && atom2==CZ)) {
                    //bond == "ne-cz-nh") {
                    // NC2-C-NH1
                    bond_length = 1.326;
                    *std_dev = 0.020;
               } else if ((atom1==CZ && atom2==NH2) || (atom1==NH2 && atom2==CZ)) {
                    //bond == "ne-cz-nh") {
                    // NC2-C-NH1
                    bond_length = 1.326;
                    *std_dev = 0.020;
               }
               break;
          case SER:                                        // Serine
               if (  (atom1==CB && atom2==OG) || (atom1==OG && atom2==CB)) {
                    //bond == "ca-cb-og") {
                    // CH1E-CH2E-OH1
                    bond_length = 1.417;
                    *std_dev = 0.020;
               }
               break;
          case THR:                                        // Threonine
               if (  (atom1==CB && atom2==OG1) || (atom1==OG1 && atom2==CB)) {
                    //bond == "ca-cb-og") {
                    // CH1E-CH1E-OH1
                    bond_length = 1.433;
                    *std_dev = 0.020;
               } else if ((atom1==CB && atom2==CG2) || (atom1==CG2 && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH1E-CH3E
                    bond_length = 1.521;
                    *std_dev = 0.020;
               }
               break;
          case VAL:                                        // Valine
               if (  (atom1==CB && atom2==CG1) || (atom1==CG1 && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH1E-CH3E
                    bond_length = 1.521;
                    *std_dev = 0.020;
               } else if (  (atom1==CB && atom2==CG2) || (atom1==CG2 && atom2==CB)) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH1E-CH3E
                    bond_length = 1.521;
                    *std_dev = 0.020;
               }
               break;
          case TRP:                                        // Tryptophan
               if ( (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    // "ca-cb-cg") {
                    // C5W-CH2E-CH1E
                    bond_length = 1.498;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==CD1) || (atom1==CD1 && atom2==CG)) {
                    // "cb-cg-cd1") {
                    // CW-C5W-CH2E
                    bond_length = 1.433;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==CD2) || (atom1==CD2 && atom2==CG)) {
                    // "cb-cg-cd2") {
                    // CH2E-C5W-CR1E
                    bond_length = 1.365;
                    *std_dev = 0.020;
               } else if ( (atom1==CD1 && atom2==NE1) || (atom1==NE1 && atom2==CD1)) {
                    // "cg-cd-ne") {
                    // C5W-CR1E-NH1
                    bond_length = 1.374;
                    *std_dev = 0.020;
               } else if ( (atom1==CE2 && atom2==NE1) || (atom1==NE1 && atom2==CE2)) {
                    // "cd1-ne1-CE2") {
                    // CR1E-NH1-CW
                    bond_length = 1.370;
                    *std_dev = 0.020;
               } else if ( (atom1==CD2 && atom2==CE2) || (atom1==CE2 && atom2==CD2)) {
                    // "cg-cd-ce") {
                    // C5W-CW-CW
                    bond_length = 1.409;
                    *std_dev = 0.020;
               } else if ( (atom1==CD2 && atom2==CE3) || (atom1==CE3 && atom2==CD2)) {
                    // "cg-cd-ce3") {
                    // C5W-CW-CR1E
                    bond_length = 1.398;
                    *std_dev = 0.020;
               } else if ( (atom1==CE2 && atom2==CZ2) || (atom1==CZ2 && atom2==CE2)) {
                    // "cd-ce-cz2") {
                    // CW-CW-CR1W
                    bond_length = 1.394;
                    *std_dev = 0.020;
               } else if ( (atom1==CZ2 && atom2==CH2) || (atom1==CH2 && atom2==CZ2)) {
                    // "ce-cz2-ch2") {
                    // CW-CR1W-CR1W
                    bond_length = 1.368;
                    *std_dev = 0.020;
               } else if ( (atom1==CE3 && atom2==CZ3) || (atom1==CZ3 && atom2==CE3)) {
                    // "cd-ce-cz3") {
                    // CW-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               } else if ( (atom1==CZ3 && atom2==CH2) || (atom1==CH2 && atom2==CZ3)) {
                    // "ce3-cz3-ch") {
                    // CR1E-CR1E-CR1W
                    bond_length = 1.400;
                    *std_dev = 0.020;
               }
               break;
          case TYR:                                        // Tyrosine
               if ( (atom1==CB && atom2==CG) || (atom1==CG && atom2==CB)) {
                    // "ca-cb-cg") {
                    // CY-CH2E-CH1E
                    bond_length = 1.512;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==CD1) || (atom1==CD1 && atom2==CG)) {
                    // "cb-cg-cd") {
                    // CH2E-CY-CR1E
                    bond_length = 1.389;
                    *std_dev = 0.020;
               } else if ( (atom1==CG && atom2==CD2) || (atom1==CD2 && atom2==CG)) {
                    // "cb-cg-cd") {
                    // CH2E-CY-CR1E
                    bond_length = 1.389;
                    *std_dev = 0.020;
               } else if ( (atom1==CD1 && atom2==CE1) || (atom1==CE1 && atom2==CD1) ) {
                    // "cg-cd-ce") {
                    // CY-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               } else if ( (atom1==CD2 && atom2==CE2) || (atom1==CE2 && atom2==CD2) ) {
                    // "cg-cd-ce") {
                    // CY-CR1E-CR1E
                    bond_length = 1.382;
                    *std_dev = 0.020;
               } else if ( (atom1==CE1 && atom2==CZ) || (atom1==CZ && atom2==CE1)) {
                    // "cd-ce-cz") {
                    // CY2-CR1E-CR1E
                    bond_length = 1.378;
                    *std_dev = 0.020;
               } else if ( (atom1==CE2 && atom2==CZ) || (atom1==CZ && atom2==CE2)) {
                    // "cd-ce-cz") {
                    // CY2-CR1E-CR1E
                    bond_length = 1.378;
                    *std_dev = 0.020;
               } else if ( (atom1==CZ && atom2==OH) || (atom1==OH && atom2==CZ)) {
                    // "ce-cz-oh") {
                    // CR1E-CY2-OH1
                    bond_length = 1.376;
                    *std_dev = 0.020;
               }
               break;
          default:
               break;
          }
     }

     // bond lengths for pseudo sidechain bonds
     // These values are not actually fixed - they serve as initial values
     if (!is_initialized(bond_length)) {
          if (atom1==PS || atom2==PS) {
               switch(res) {
               case ALA:                                        // Alanine
                    bond_length = 1.54;
                    *std_dev = 0.0;
                    break;
               case CYS:                                        // Cysteine Acid
                    bond_length = 2.8;
                    *std_dev = 0.0;
                    break;
               case ASP:                                        // Aspartic Acid
                    bond_length = 2.92;
                    *std_dev = 0.0;
                    break;
               case GLU:                                        // Glutamic Acid
                    bond_length = 3.125;
                    *std_dev = 0.0;
                    break;
               case PHE:                                        // Phenylalanine
                    bond_length = 3.79;
                    *std_dev = 0.0;
                    break;
               case GLY:                                        // Glycine
                    bond_length = 1.0;
                    *std_dev = 0.0;
                    break;
               case HIS:                                        // Histidine
                    bond_length = 3.57;
                    *std_dev = 0.0;
                    break;
               case ILE:                                        // Isoleucine
                    bond_length = 2.7;
                    *std_dev = 0.0;
                    break;
               case LYS:                                        // Lysine
                    bond_length = 4.6;
                    *std_dev = 0.0;
                    break;
               case LEU:                                        // Leucine
                    bond_length = 3.05;
                    *std_dev = 0.0;
                    break;
               case MET:                                        // Methionine
                    bond_length = 3.185;
                    *std_dev = 0.0;
                    break;
               case ASN:                                        // Asparagine
                    bond_length = 2.91;
                    *std_dev = 0.0;
                    break;
               case PRO:                                        // Proline
                    bond_length = 2.29;
                    *std_dev = 0.0;
                    break;
               case GLN:                                        // Glutamine
                    bond_length = 3.875;
                    *std_dev = 0.0;
                    break;
               case ARG:                                        // Arginine
                    bond_length = 4.8;
                    *std_dev = 0.0;
                    break;
               case SER:                                        // Serine
                    bond_length = 2.43;
                    *std_dev = 0.0;
                    break;
               case THR:                                        // Threonine
                    bond_length = 2.17;
                    *std_dev = 0.0;
                    break;
               case VAL:                                        // Valine
                    bond_length = 2.19;
                    *std_dev = 0.0;
                    break;
               case TRP:                                        // Tryptophan
                    bond_length = 4.2;
                    *std_dev = 0.0;
                    break;
               case TYR:                                        // Tyrosine
                    bond_length = 4.27;
                    *std_dev = 0.0;
                    break;
               default:
                    break;
               }
          }
     }
     
     assert(is_initialized(bond_length));
     return bond_length;
}


//! Dihedral angle constants
//!
//! Dihedrals are either specified in the normal sense:,
//!     e.g. phi(i) is defined by: C(i-1), N(i), CA(i), C(i+1)
//! To allow the position of for instance H, and O atoms to be specified by
//! a dihedral and an angle, the dihedral is sometimes defined differently:
//!     e.g. the dihedral to specify H(i) is defined by:
//!          N(i), C(i), CA(i), H(i)
//! HOWEVER, WHEN THIS SCHEME IS USED, THE POSITION OF THE TARGET ATOM IS FIRST TRANSLATED
//!                                                       //
//!            o H                  o H                   //
//!            |                    |\                    //
//!            |                    | \                   //
//!            o N                  o  o H'               //
//!           / \                  / \                    //
//!          /   \                /   \                   //
//!       C o     o CA         C o     o CA               //
//!                                                       //
//! This corresponds to having an C-CA directed axis going through N, with the dihedral
//! defined as the rotation around this axis away from the plane.
//!
//! \param atom1 First atom
//! \param atom2 Second atom
//! \param atom3 Third atom
//! \param atom4 Fourth atom
//! \param res Residue containing atoms
//! \param std_dev Optional output variable for standard deviation
//! \param trans Whether residue is in a trans state
//! \return dihedral angle
inline double dihedral_constants(enum definitions::AtomEnum atom1, 
                                 enum definitions::AtomEnum atom2,
                                 enum definitions::AtomEnum atom3,
                                 enum definitions::AtomEnum atom4,
                                 enum definitions::ResidueEnum res=definitions::AA_UNDEF,
                                 double *std_dev=NULL,
                                 bool trans=true) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Set std_dev pointer to point to local variable if not provided
     if (std_dev==NULL) {
          double std_dev_tmp;
          std_dev = &std_dev_tmp;
     }
     
     double dihedral = UNINITIALIZED;

     // Omega angle
     // From J. P. Priestle: Improved dihedral-angle restraints for protein structure refinement
     if (atom1==CA && atom2==C && atom3==N && atom4==CA) {
          dihedral = 179.3;
          if (!trans)
               dihedral = 0.5;
          *std_dev = 6.2;

          
     // Residue unspecific constants
     } else if (atom1==CA && atom2==N && atom3==C && atom4==CB) {    // CB atom
          dihedral = 52.25;
          *std_dev = 0.041;
     } else if (atom1==CA && atom2==CA && atom3==CA && atom4==CB) {  // CB atom - in CA-only representation
          dihedral = 39.31;
          *std_dev = 0.0;
     } else if (atom1==C && atom2==CA && atom3==N && atom4==O){      // O atom
          dihedral = 0.0;
          *std_dev = 0.0;
     } else if (atom1==N && atom2==CA && atom3==C && atom4==O){      // O atom - terminal
          dihedral = 0.0;
          *std_dev = 0.0;
     } else if (atom1==O && atom2==CA && atom3==C && atom4==OXT) {   // OXT atom
          dihedral = 180.0;
          *std_dev = 0.0;
     } else if (atom1==N && atom2==C && atom3==CA && atom4==H) {     // H atom
          dihedral = 0.0;
          *std_dev = 0.0;
     } else if (atom1==C && atom2==CA && atom3==N && atom4==H1) {     // H1 atom
          dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral;
          *std_dev = 0.0;
     } else if (atom1==H1 && atom2==CA && atom3==N && atom4==H2) {   // H2 atom
          dihedral = 120.0;
          *std_dev = 0.0;
     } else if (atom1==H1 && atom2==CA && atom3==N && atom4==H3) {   // H3 atom
          dihedral = 240.0;
          *std_dev = 0.0;
     } else if (atom1==CA && atom2==N && atom3==C && atom4==HA) {    // HA atom
          dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral;
          *std_dev = 0.0;
     }


     // Residue specific constants - Sidechain heavy atoms
     if (!is_initialized(dihedral)) {
          switch (res) {
          case ASP:                                        // Aspartic Acid
               if (atom1==OD1 && atom2==CB && atom3==CG && atom4==OD2) {            // OD2
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          case GLU:                                        // Glutamic Acid
               if (atom1==OE1 && atom2==CG && atom3==CD && atom4==OE2) {            // OE2
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          case PHE:                                        // Phenylalanine
               // if (atom1==CD1 && atom2==CB && atom3==CG && atom4==CD2) {            // CD2
               if (atom1==CB && atom2==CD1 && atom3==CG && atom4==CD2) {            // CD2
                    dihedral = 180;
                    *std_dev = 0.0;
               // } else if (atom1==CB && atom2==CG && atom3==CD1 && atom4==CE1) {     // CE1
               } else if (atom1==CD2 && atom2==CG && atom3==CD1 && atom4==CE1) {     // CE1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CD1 && atom2==CG && atom3==CD2 && atom4==CE2) {     // CE2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==CD1 && atom3==CE1 && atom4==CZ) {     // CZ
                    dihedral = 0.0;
                    *std_dev = 0.0;
               }
               break;
          case GLY:                                        // Glycine
               break;
          case HIS:                                        // Histidine
               // if (atom1==ND1 && atom2==CB && atom3==CG && atom4==CD2) {            // CD2
               if (atom1==CB && atom2==ND1 && atom3==CG && atom4==CD2) {            // CD2
                    dihedral = 180.0;
                    *std_dev = 0.0;
               // } else if (atom1==CB && atom2==CG && atom3==CD2 && atom4==NE2) {     // NE2
               } else if (atom1==ND1 && atom2==CG && atom3==CD2 && atom4==NE2) {     // NE2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               // } else if (atom1==CB && atom2==CG && atom3==ND1 && atom4==CE1) {     // CE1
               } else if (atom1==CD2 && atom2==CG && atom3==ND1 && atom4==CE1) {     // CE1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               }
               break;
          case ILE:                                        // Isoleucine
               if (atom1==CG1 && atom2==CA && atom3==CB && atom4==CG2) {            // CG2
                    dihedral = 4.0*180/3.0;
                    *std_dev = 0.0;
               }
               break;
          case LYS:                                        // Lysine
               break;
          case LEU:                                        // Leucine
               if (atom1==CD1 && atom2==CB && atom3==CG && atom4==CD2) {            // CD2
                    dihedral = 2.0*180/3.0;
                    *std_dev = 0.0;
               }
               break;
          case MET:                                        // Methionine
               break;
          case ASN:                                        // Asparagine
               if (atom1==OD1 && atom2==CB && atom3==CG && atom4==ND2) {            // ND2
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          case PRO:                                        // Proline
               break;
          case GLN:                                        // Glutamine
               if (atom1==OE1 && atom2==CG && atom3==CD && atom4==NE2) {           // NE2
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          case ARG:                                        // Arginine
               if (atom1==CD && atom2==NE && atom3==CZ && atom4==NH1) {            // NH1
                    dihedral = 0;
                    *std_dev = 0.0;
               } else if (atom1==CD && atom2==NE && atom3==CZ && atom4==NH2) {     // NH2
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          case SER:                                        // Serine
               break;
          case THR:                                        // Threonine
               if (atom1==OG1 && atom2==CA && atom3==CB && atom4==CG2) {           // CG2
                    dihedral = 4.0*180/3.0;
                    *std_dev = 0.0;
               }
               break;
          case VAL:                                        // Valine
               if (atom1==CG1 && atom2==CA && atom3==CB && atom4==CG2) {           // CG2
                    dihedral = 2.0*180/3.0;
                    *std_dev = 0.0;
               }
               break;
          case TRP:                                        // Tryptophan
               if (atom1==CB && atom2==CD1 && atom3==CG && atom4==CD2) {          // CD2
                    dihedral = 180;
                    *std_dev = 0.0;
               } else if (atom1==CD2 && atom2==CG && atom3==CD1 && atom4==NE1) {   // NE1
                    dihedral = 0;
                    *std_dev = 0.0;
               } else if (atom1==CD1 && atom2==CG && atom3==CD2 && atom4==CE2) {   // CE2
                    dihedral = 0;
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==CE2 && atom3==CD2 && atom4==CE3) {   // CE3
                    dihedral = 180;
                    *std_dev = 0.0;
               } else if (atom1==CE3 && atom2==CD2 && atom3==CE2 && atom4==CZ2) {  // CZ2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE2 && atom2==CD2 && atom3==CE3 && atom4==CZ3) {  // CZ3
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CD2 && atom2==CE3 && atom3==CZ3 && atom4==CH2) { // CH2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               }
               break;
          case TYR:                                        // Tyrosine
               if (atom1==CB && atom2==CD1 && atom3==CG && atom4==CD2) {          // CD2
                    dihedral = 180;
                    *std_dev = 0.0;
               } else if (atom1==CD2 && atom2==CG && atom3==CD1 && atom4==CE1) {   // CE1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CD1 && atom2==CG && atom3==CD2 && atom4==CE2) {   // CE2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==CD1 && atom3==CE1 && atom4==CZ) {   // CZ
                    dihedral = 0;
                    *std_dev = 0.0;
               } else if (atom1==CD1 && atom2==CE1 && atom3==CZ && atom4==OH) {   // CZ
                    dihedral = 180;
                    *std_dev = 0.0;
               }
          default:
               break;
          }
     }
          
     // Sidechain dihedral values for default rotamers.
     // These are not actual constants but serve as initial values
     // The values are the most probable rotamer states taken from:
     // The Penultimate Rotamer Library. Lovell S, Richardson DC et al.,
     // PROTEINS: Structure, Function and Genetics 40:389-408 (2000)
     if (!is_initialized(dihedral)) {

          switch (res) {
          case CYS:
               if (atom1==N && atom2==CA && atom3==CB && atom4==SG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               }
               break;
          case ASP:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -70.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==OD1) {
                    dihedral = -15.;
                    *std_dev = 0.0;
               }
               break;
          case GLU:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==CD && atom4==OE1) {
                    dihedral = -40.;
                    *std_dev = 0.0;
               }
               break;
          case PHE:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD1) {
                    dihedral = -85.;
                    *std_dev = 0.0;
               }
               break;
          case HIS:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==ND1) {
                    dihedral = -70.;
                    *std_dev = 0.0;
               }
               break;
          case ILE:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG1) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG1 && atom4==CD1) {
                    dihedral = 170.;
                    *std_dev = 0.0;
               }
               break;
          case LYS:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -67.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==CD && atom4==CE) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==CD && atom3==CE && atom4==NZ) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               }
               break;
          case LEU:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -177.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD1) {
                    dihedral = 65.;
                    *std_dev = 0.0;
               }
               break;
          case MET:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==SD) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==SD && atom4==CE) {
                    dihedral = -70.;
                    *std_dev = 0.0;
               }
               break;
          case ASN:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==OD1) {
                    dihedral = -20.;
                    *std_dev = 0.0;
               }
               break;
          case PRO:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = 30.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD) {
                    dihedral = -40.;
                    *std_dev = 0.0;
               }
               break;
          case GLN:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -67.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==CD && atom4==OE1) {
                    dihedral = -25.;
                    *std_dev = 0.0;
               }               
               break;
          case ARG:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -67.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==CD && atom4==NE) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==CD && atom3==NE && atom4==CZ) {
                    dihedral = 180.;
                    *std_dev = 0.0;
               }
               break;
          case SER:
               if (atom1==N && atom2==CA && atom3==CB && atom4==OG) {
                    dihedral = 62.;
                    *std_dev = 0.0;
               }
               break;
          case THR:
               if (atom1==N && atom2==CA && atom3==CB && atom4==OG1) {
                    dihedral = 62.;
                    *std_dev = 0.0;
               }
               break;
          case VAL:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG1) {
                    dihedral = 175.;
                    *std_dev = 0.0;
               }
               break;
          case TRP:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD1) {
                    dihedral = 95.;
                    *std_dev = 0.0;
               }
          case TYR:
               if (atom1==N && atom2==CA && atom3==CB && atom4==CG) {
                    dihedral = -65.;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG && atom4==CD1) {
                    dihedral = -85.;
                    *std_dev = 0.0;
               }               
               break;
          default:
               break;
          }
     }


     // Pseudo-sidechain dihedral values.
     // These are not actual constants but serve as initial values
     // The values reflect the most probable rotamer states
     if (!is_initialized(dihedral)) {
          if (atom4==PS) {
               switch (res) {
               case ALA:           
                    dihedral = 52.25;
                    *std_dev = 0.0;
                    break;
               case CYS:           
                    dihedral = 28.6;
                    *std_dev = 0.0;
                    break;
               case ASP:           
                    dihedral = 28.6;
                    *std_dev = 0.0;
                    break;
               case GLU:           
                    dihedral = 40.1;
                    *std_dev = 0.0;
                    break;
               case PHE:
                    dihedral = 22.9;
                    *std_dev = 0.0;
                    break;
               case GLY:           
                    dihedral = 52.25;
                    *std_dev = 0.0;
                    break;
               case HIS:           
                    dihedral = 22.9;
                    *std_dev = 0.0;
                    break;
               case ILE:
                    dihedral = 37.2;
                    *std_dev = 0.0;
                    break;
               case LYS:
                    dihedral = 34.4;
                    *std_dev = 0.0;
                    break;
               case LEU:
                    dihedral = 28.6;
                    *std_dev = 0.0;
                    break;
               case MET:
                    dihedral = 40.1;
                    *std_dev = 0.0;
                    break;
               case ASN:
                    dihedral = 25.8;
                    *std_dev = 0.0;
                    break;
               case PRO:
                    dihedral = 100.3;
                    *std_dev = 0.0;
                    break;
               case GLN:
                    dihedral = 37.2;
                    *std_dev = 0.0;
                    break;
               case ARG:
                    dihedral = 28.6;
                    *std_dev = 0.0;
                    break;
               case SER:
                    dihedral = 85.9;
                    *std_dev = 0.0;
                    break;
               case THR:
                    dihedral = 34.4;
                    *std_dev = 0.0;
                    break;
               case VAL:
                    dihedral = 34.4;
                    *std_dev = 0.0;
                    break;
               case TRP:
                    dihedral = 5.7;
                    *std_dev = 0.0;
                    break;
               case TYR:
                    dihedral = 22.9;
                    *std_dev = 0.0;
                    break;
               default:
                    break;
               }
          }
     }


     // Semi residue-specific constants
     if (!is_initialized(dihedral)) {
          if (atom1==CB && atom2==CA && atom3==CG && atom4==HB2) {           // HB2 
               dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CB && atom2==CA && atom3==CG && atom4==HB3) {    // HB3 
               dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CG && atom2==CB && atom3==CD && atom4==HG2) {    // HG2 
               dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CG && atom2==CB && atom3==CD && atom4==HG3) {    // HG3 
               dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CD && atom2==CG && atom3==CE && atom4==HD2) {    // HD2 
               dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CD && atom2==CG && atom3==CE && atom4==HD3) {    // HD3 
               dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          }
     }

          
     // Residue specific constants - hydrogens
     if (!is_initialized(dihedral)) {
          switch (res) {
          case ALA:                                        // Alanine
               if (atom1==N && atom2==CA && atom3==CB && atom4==HB1) {            // HB1 
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HB1 && atom2==CA && atom3==CB && atom4==HB2) {   // HB2
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HB1 && atom2==CA && atom3==CB && atom4==HB3) {   // HB3
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case CYS:                                        // Cysteine
               if (atom1==CB && atom2==CA && atom3==SG && atom4==HB2) {           // HB2
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CA && atom3==SG && atom4==HB3) {    // HB3
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==SG && atom4==HG) {     // HG 
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          case PHE:                                       // Phenylalanine
               if (atom1==CD1 && atom2==CG && atom3==CE1 && atom4==HD1) {         // HD1 
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CD2 && atom2==CG && atom3==CE2 && atom4==HD2) {  // HD2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE1 && atom2==CD1 && atom3==CZ && atom4==HE1) {  // HE1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE2 && atom2==CD2 && atom3==CZ && atom4==HE2) {  // HE2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CZ && atom2==CE1 && atom3==CE2 && atom4==HZ) {   // HZ
                    dihedral = 0.0;
                    *std_dev = 0.0;
               }
               break;
          case GLY:                                        // Glycine
               if (atom1==CA && atom2==N && atom3==C && atom4==HA2) {             // HA2
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==N && atom3==C && atom4==HA3) {      // HA3
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180;
                    *std_dev = 0.0;
               }
               break;
          case HIS:                                        // Histidine
               if (atom1==ND1 && atom2==CG && atom3==CE1 && atom4==HD1) {          // HD1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CD2 && atom2==CG && atom3==NE2 && atom4==HD2) {   // HD2 
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE1 && atom2==ND1 && atom3==NE2 && atom4==HE1) {  // HE1 
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==NE2 && atom2==CD2 && atom3==CE1 && atom4==HE2) {  // HE2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               }
               break;
          case ILE:                                        // Isoleucine
               if (atom1==CB && atom2==CA && atom3==CG1 && atom4==HB) {              // HB
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CG1 && atom2==CB && atom3==CD1 && atom4==HG12) {    // HG12
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CG1 && atom2==CB && atom3==CD1 && atom4==HG13) {    // HG13
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG2 && atom4==HG21) {     // HG21 
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HG21 && atom2==CB && atom3==CG2 && atom4==HG22) {   // HG22
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HG21 && atom2==CB && atom3==CG2 && atom4==HG23) {  // HG23
                    dihedral = 240.0;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG1 && atom3==CD1 && atom4==HD11) {    // HD11
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HD11 && atom2==CG1 && atom3==CD1 && atom4==HD12) {  // HD12
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HD11 && atom2==CG1 && atom3==CD1 && atom4==HD13) {  // HD13
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case LYS:                                        // Lysine
               if (atom1==CE && atom2==CD && atom3==NZ && atom4==HE2) {               // HE2
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CE && atom2==CD && atom3==NZ && atom4==HE3) {        // HE3
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CD && atom2==CE && atom3==NZ && atom4==HZ1) {        // HZ1 
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HZ1 && atom2==CE && atom3==NZ && atom4==HZ2) {       // HZ2 
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HZ1 && atom2==CE && atom3==NZ && atom4==HZ3) {       // HZ3 
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case LEU:                                        // Leucine
               if (atom1==CG && atom2==CB && atom3==CD1 && atom4==HG) {              // HG
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==CD1 && atom4==HD11) {     // HD11
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HD11 && atom2==CG && atom3==CD1 && atom4==HD12) {   // HD12
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HD11 && atom2==CG && atom3==CD1 && atom4==HD13) {   // 3HD1
                    dihedral = 240.0;
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CG && atom3==CD2 && atom4==HD21) {     // HD21
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HD21 && atom2==CG && atom3==CD2 && atom4==HD22) {   // HD22
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HD21 && atom2==CG && atom3==CD2 && atom4==HD23) {   // HD23
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case MET:                                        // Methionine
               if (atom1==CG && atom2==CB && atom3==SD && atom4==HG2) {               // HG2
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==CB && atom3==SD && atom4==HG3) {        // HG3
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CG && atom2==SD && atom3==CE && atom4==HE1) {        // HE1
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HE1 && atom2==SD && atom3==CE && atom4==HE2) {       // HE2
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HE1 && atom2==SD && atom3==CE && atom4==HE3) {       // HE3
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case ASN:                                        // Asparagine
               if (atom1==CB && atom2==CG && atom3==ND2 && atom4==HD21) {            // HD21
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HD21 && atom2==CG && atom3==ND2 && atom4==HD22) {   // HD22
                    dihedral = 180.0;
                    *std_dev = 0.0;
               }
               break;
          case PRO:                                        // Proline
               if (atom1==CD && atom2==CG && atom3==N && atom4==HD2) {                // HD2 
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CD && atom2==CG && atom3==N && atom4==HD3) {        // HD3 
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               }
               break;
          case GLN:                                        // Glutamine
               if (atom1==CG && atom2==CD && atom3==NE2 && atom4==HE21) {             // HE21
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HE21 && atom2==CD && atom3==NE2 && atom4==HE22) {    // HE22
                    dihedral = 180.0;
                    *std_dev = 0.0;
               }
               break;
          case ARG:                                        // Arginine
               if (atom1==CD && atom2==CG && atom3==NE && atom4==HD2) {               // HD2 
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CD && atom2==CG && atom3==NE && atom4==HD3) {        // HD3 
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==NE && atom2==CD && atom3==CZ && atom4==HE) {         // HE
                    dihedral = 0.0; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==NE && atom2==CZ && atom3==NH1 && atom4==HH11) {      // HH11
                    dihedral = 0.0; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==HH11 && atom2==CZ && atom3==NH1 && atom4==HH12) {    // HH12
                    dihedral = 180.0; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==NE && atom2==CZ && atom3==NH2 && atom4==HH21) {      // HH21
                    dihedral = 0.0; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==HH21 && atom2==CZ && atom3==NH2 && atom4==HH22) {    // HH22
                    dihedral = 180.0; // tetrahedral
                    *std_dev = 0.0;
               }
               break;
          case SER:                                        // Serine
               if (atom1==CB && atom2==CA && atom3==OG && atom4==HB2) {           // HB2
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CB && atom2==CA && atom3==OG && atom4==HB3) {    // HB3
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==OG && atom4==HG) {     // HG 
                    dihedral = 180; 
                    *std_dev = 0.0;
               }               
               break;
          case THR:                                        // Threonine
               if (atom1==CB && atom2==CA && atom3==OG1 && atom4==HB) {              // HB
                    dihedral = -0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==OG1 && atom4==HG1) {      // HG1 
                    dihedral = 180; 
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG2 && atom4==HG21) {     // HG21 
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HG21 && atom2==CB && atom3==CG2 && atom4==HG22) {   // HG22
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HG21 && atom2==CB && atom3==CG2 && atom4==HG23) {   // HG23
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case VAL:                                        // Valine
               if (atom1==CB && atom2==CA && atom3==CG1 && atom4==HB) {              // HB
                    dihedral = 0.5*acos(-1.0/3.0)/M_PI*180; // tetrahedral
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG1 && atom4==HG11) {     // HG11 
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HG11 && atom2==CB && atom3==CG1 && atom4==HG12) {   // HG12
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HG11 && atom2==CB && atom3==CG1 && atom4==HG13) {   // HG13
                    dihedral = 240.0;
                    *std_dev = 0.0;
               } else if (atom1==CA && atom2==CB && atom3==CG2 && atom4==HG21) {     // HG21 
                    dihedral = 180.0;
                    *std_dev = 0.0;
               } else if (atom1==HG21 && atom2==CB && atom3==CG2 && atom4==HG22) {   // HG22
                    dihedral = 120.0;
                    *std_dev = 0.0;
               } else if (atom1==HG21 && atom2==CB && atom3==CG2 && atom4==HG23) {   // HG23
                    dihedral = 240.0;
                    *std_dev = 0.0;
               }
               break;
          case TRP:                                        // Tryptophane
               if (atom1==CD1 && atom2==CG && atom3==NE1 && atom4==HD1) {            // HD1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==NE1 && atom2==CD1 && atom3==CE2 && atom4==HE1) {    // HE1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE3 && atom2==CD2 && atom3==CZ3 && atom4==HE3) {    // HE3
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CZ2 && atom2==CE2 && atom3==CH2 && atom4==HZ2) {    // HZ2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CZ3 && atom2==CE3 && atom3==CH2 && atom4==HZ3) {    // HZ3
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CH2 && atom2==CZ2 && atom3==CZ3 && atom4==HH2) {    // HH2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               }               
               break;
          case TYR:                                        // Tyrosine
               if (atom1==CD1 && atom2==CG && atom3==CE1 && atom4==HD1) {         // HD1 
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CD2 && atom2==CG && atom3==CE2 && atom4==HD2) {  // HD2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE1 && atom2==CD1 && atom3==CZ && atom4==HE1) {  // HE1
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE2 && atom2==CD2 && atom3==CZ && atom4==HE2) {  // HE2
                    dihedral = 0.0;
                    *std_dev = 0.0;
               } else if (atom1==CE1 && atom2==CZ && atom3==OH && atom4==HH) {    // HH
                    dihedral = 180;
                    *std_dev = 0.0;
               }
               break;
          default:
               break;
          }
     }


     
     assert(is_initialized(dihedral));

     double conversion_factor = M_PI/180;
     dihedral *= conversion_factor;
     *std_dev *= conversion_factor;
     return dihedral;

}

//! Bond angle constants
//!
//! \param atom1 First atom
//! \param atom2 Second atom
//! \param atom3 Third atom
//! \param res Residue containing atoms
//! \param std_dev Optional output variable for standard deviation
//! \param terminal_status Whether residue is located at the N-terminal, C-terminal or internally in the chain
//! \return bond angle
inline double bond_angle_constants(enum definitions::AtomEnum atom1, 
                                   enum definitions::AtomEnum atom2,
                                   enum definitions::AtomEnum atom3,
                                   enum definitions::ResidueEnum res=definitions::AA_UNDEF,
                                   double *std_dev=NULL,
                                   definitions::TerminalEnum terminal_status=definitions::INTERNAL) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Set std_dev pointer to point to local variable if not provided
     if (std_dev==NULL) {
          double std_dev_tmp;
          std_dev = &std_dev_tmp;
     }

     double bond_angle = UNINITIALIZED;

     // Backbone contants from Engh-Huber
     if (atom1==C && atom2==N && atom3==CA) {            // C - N - CA 
          if (res==GLY) {
               bond_angle = 120.6;
               *std_dev = 1.7;
          } else if (res==PRO) {
               bond_angle = 122.6;
               *std_dev = 5.0;
          } else {
               bond_angle = 121.7;
               *std_dev = 1.8;
          }
     } else if (atom1==N && atom2==CA && atom3==C) {     // N - CA - C
          if (res==GLY) {
               bond_angle = 112.5;
               *std_dev = 2.9;
          } else if (res==PRO) {
               bond_angle = 111.8;
               *std_dev = 2.5;
          } else {
               bond_angle = 111.2;
               *std_dev = 2.8;
          }
     } else if (atom1==CA && atom2==C && atom3==N) {     // CA - C - N 
          if (res==GLY) {
               bond_angle = 116.4;
               *std_dev = 2.1;
          } else if (res==PRO) {
               bond_angle = 116.9;
               *std_dev = 1.5;
          } else {
               bond_angle = 116.2;
               *std_dev = 2.0;
          }
     } else if (atom1==CA && atom2==C && atom3==O) {    // CA - C - O 
          if (terminal_status==CTERM) {
               bond_angle = 117.0;
               *std_dev = 2.5;
          } else {
               if (res==GLY) {
                    bond_angle = 120.8;
                    *std_dev = 2.1;
               } else {
                    bond_angle = 120.8;
                    *std_dev = 1.7;
               }               
          }
     } else if (atom1==CA && atom2==C && atom3==OXT) {   // CA - C - OXT 
          bond_angle = 117.0;
          *std_dev = 2.5;
     }

     

     // Side chain contants from Engh-Huber
     if (!is_initialized(bond_angle)) {
          switch (res) {
          case ALA:                                        // Alanine
               break;
          case CYS:                                        // Cysteine
               if (atom1==CA && atom2==CB && atom3==SG) {
                    // "ca-cb-sg") {
                    // CH1E-CH2E-S
                    bond_angle  = 114.4;
                    *std_dev = 0.020;
               }
               break;
          case ASP:                                        // Aspartic Acid
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    // ca-cb-cg
                    // C-CH2E-CH1E
                    bond_angle  = 112.6;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==OD1) {
                    // cb-cg-od
                    // CH2E-C-OC
                    bond_angle  = 118.4;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==OD2) {
                    // cb-cg-od
                    // CH2E-C-OC
                    bond_angle  = 118.4;
                    *std_dev = 0.020;
               }
               break;
          case GLU:                                        // Glutamic Acid
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    // ca-cb-cg
                    // CH1E-CH2E-CH2E
                    bond_angle  = 114.1;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD) {
                    // cb-cg-cd
                    // C-CH2E-CH2E
                    bond_angle  = 112.6;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD && atom3==OE1) {
                    // cg-cd-oe
                    // CH2E-C-OC
                    bond_angle  = 118.4;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD && atom3==OE2) {
                    // cg-cd-oe
                    // CH2E-C-OC
                    bond_angle  = 118.4;
                    *std_dev = 0.020;
               }
               break;
          case PHE:                                        // Phenylalanine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    // ca-cb-cg") {
                    // CF-CH2E-CH1E
                    bond_angle  = 113.8;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD1) {
                    // cb-cg-cd") {
                    // CH2E-CF-CR1E
                    bond_angle  = 120.7;
                    *std_dev = 0.020;
               } else if (atom1==CD1 && atom2==CG && atom3==CD2) {
                    // cb-cg-cd") {
                    // CH2E-CF-CR1E
                    // bondAngle  = 120.7;
                    bond_angle  = 118.6;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD1 && atom3==CE1) {
                    // bond == "cg-cd-ce") {
                    // CF-CR1E-CR1E
                    bond_angle  = 120.7;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD2 && atom3==CE2) {
                    // bond == "cg-cd-ce") {
                    // CF-CR1E-CR1E
                    bond_angle  = 120.7;
                    *std_dev = 0.020;
               } else if (atom1==CD1 && atom2==CE1 && atom3==CZ) {
                    // bond == "cd-ce-cz") {
                    // CR1E-CR1E-CR1E
                    bond_angle  = 120.0;
                    *std_dev = 0.020;
               } else if (atom1==CD2 && atom2==CE2 && atom3==CZ) {
                    // bond == "cd-ce-cz") {
                    // CR1E-CR1E-CR1E
                    bond_angle  = 120.0;
                    *std_dev = 0.020;
               } else if (atom1==CE2 && atom2==CZ && atom3==CE1) {
                    // bond == "ce2-cz-ce1") {
                    // CR1E-CR1E-CR1E
                    bond_angle  = 120.0;
                    *std_dev = 0.020;
               }
               break;
          case GLY:                                        // Glycine
               break;
          case HIS:                                        // Histidine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // C5-CH2E-CH1E
                    bond_angle  = 113.8;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==ND1) {
                    //bond == "cb-cg-nd1") {
                    // CH2E-C5-NR
                    bond_angle  = 121.6;
                    *std_dev = 0.020;
               // } else if (atom1==CB && atom2==CG && atom3==CD2) {
               } else if (atom1==ND1 && atom2==CG && atom3==CD2) {
                    //bond == "cb-cg-cd2") {
                    // CH2E-C5-CR1E
                    // bondAngle  = 129.1;
                    bond_angle  = 109.3;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==ND1 && atom3==CE1) {
                    //bond == "cg-nd1-ce1") {
                    // C5-NR-CR1E
                    bond_angle  = 105.6;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD2 && atom3==NE2) {
                    //bond == "cg-cd2-ne2") {
                    // C5-CR1E-NH1
                    bond_angle  = 106.5;
                    *std_dev = 0.020;
               } else if (atom1==ND1 && atom2==CE1 && atom3==NE2) {
                    //bond == "nd1-ce1-ne2") {
                    // NR-CR1E-NH1
                    bond_angle  = 111.7;
                    *std_dev = 0.020;
               } else if (atom1==CD2 && atom2==NE2 && atom3==CE1 ) {
                    //bond == "ce1-ne2-CD2") {
                    // CR1E-NH1-CR1E
                    bond_angle  = 107.0;
                    *std_dev = 0.020;
               }
               break;
          case ILE:                                        // Isoleucine
               if (atom1 == CA && atom2==CB && atom3==CG1) {
                    //bond == "ca-cb-cg1") {
                    // CH1E-CH1E-CH2E
                    bond_angle  = 110.4;
                    *std_dev = 0.020;
               } else if (atom1==CA && atom2==CB && atom3==CG2) {
                    //bond == "ca-cb-cg2") {
                    // CH1E-CH1E-CH3E
                    bond_angle  = 110.5;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG1 && atom3==CD1) {
                    //bond == "cb-cg1-cd1") {
                    // CH1E-CH2E-CH3E
                    bond_angle  = 113.8;
                    *std_dev = 0.020;
               }
               break;
          case LYS:                                        // Lysine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_angle  = 114.1;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH2E-CH2E
                    bond_angle  = 111.3;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD && atom3==CE) {
                    //bond == "cg-cd-ce") {
                    // CH2E-CH2E-CH2E
                    bond_angle  = 111.3;
                    *std_dev = 0.020;
               } else if (atom1==CD && atom2==CE && atom3==NZ) {
                    //bond == "cd-ce-nz") {
                    // CH2E-CH2E-NH3
                    bond_angle  = 111.9;
                    *std_dev = 0.020;
               }
               break;
          case LEU:                                        // Leucine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH1E
                    bond_angle  = 116.3;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD1) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH1E-CH3E
                    bond_angle  = 110.7;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD2) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH1E-CH3E
                    bond_angle  = 110.7;
                    *std_dev = 0.020;
               }
               break;
          case MET:                                        // Methionine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_angle  = 114.1;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==SD) {
                    //bond == "cb-cg-sd")
                    // CH2E-CH2E-SM
                    bond_angle  = 112.7;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==SD && atom3==CE) {
                    //bond == "cg-sd-ce") {
                    // CH2E-SM-CH3E
                    bond_angle  = 100.9;
                    *std_dev = 0.020;
               }
               break;
          case ASN:                                        // Asparagine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // C-CH2E-CH1E
                    bond_angle  = 112.6;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==OD1) {
                    //bond == "cb-cg-od") {
                    // CH2E-C-O
                    bond_angle  = 120.8;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==ND2) {
                    //bond == "cb-cg-nd") {
                    // CH2E-C-NH2
                    bond_angle  = 116.4;
                    *std_dev = 0.020;
               }
               break;
          case PRO:                                        // Proline
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2P
                    bond_angle  = 104.5;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH2P-CH2P
                    bond_angle  = 106.1;
                    *std_dev = 0.020;
               }
               break;
          case GLN:                                        // Glutamine
               if (atom1==CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_angle  = 114.1;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD) {
                    //bond == "cb-cg-cd") {
                    // C-CH2E-CH2E
                    bond_angle  = 112.6;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD && atom3==OE1) {
                    //bond == "cg-cd-oe") {
                    // CH2E-C-O
                    bond_angle  = 120.8;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD && atom3==NE2) {
                    //bond == "cg-cd-ne") {
                    // CH2E-C-NH2
                    bond_angle  = 116.4;
                    *std_dev = 0.020;
               }
               break;
          case ARG:                                        // Arginine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH2E-CH2E
                    bond_angle  = 114.1;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD) {
                    //bond == "cb-cg-cd") {
                    // CH2E-CH2E-CH2E
                    bond_angle  = 111.3;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD && atom3==NE) {
                    //bond == "cg-cd-ne") {
                    // CH2E-CH2E-NH1
                    bond_angle  = 112.0;
                    *std_dev = 0.020;
               } else if (atom1==CD && atom2==NE && atom3==CZ) {
                    //bond == "cd-ne-cz") {
                    // C-NH1-CH2E
                    bond_angle  = 124.2;
                    *std_dev = 0.020;
               } else if (atom1==NE && atom2==CZ && atom3==NH1) {
                    //bond == "ne-cz-nh") {
                    // NC2-C-NH1
                    bond_angle  = 120.0;
                    *std_dev = 0.020;
               } else if (atom1==NE && atom2==CZ && atom3==NH2) {
                    //bond == "ne-cz-nh") {
                    // NC2-C-NH1
                    bond_angle  = 120.0;
                    *std_dev = 0.020;
               }
               break;
          case SER:                                        // Serine
               if (atom1 == CA && atom2==CB && atom3==OG) {
                    //bond == "ca-cb-og") {
                    // CH1E-CH2E-OH1
                    bond_angle  = 111.1;
                    *std_dev = 0.020;
               }
               break;
          case THR:                                        // Threonine
               if (atom1 == CA && atom2==CB && atom3==OG1) {
                    //bond == "ca-cb-og") {
                    // CH1E-CH1E-OH1
                    bond_angle  = 109.6;
                    *std_dev = 0.020;
               } else if (atom1==CA && atom2==CB && atom3==CG2) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH1E-CH3E
                    bond_angle  = 110.5;
                    *std_dev = 0.020;
               }
               break;
          case VAL:                                        // Valine
               if (atom1 == CA && atom2==CB && atom3==CG1) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH1E-CH3E
                    bond_angle  = 110.5;
                    *std_dev = 0.020;
               } else if (atom1==CA && atom2==CB && atom3==CG2) {
                    //bond == "ca-cb-cg") {
                    // CH1E-CH1E-CH3E
                    bond_angle  = 110.5;
                    *std_dev = 0.020;
               }
               break;
          case TRP:                                        // Tryptophan
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    // "ca-cb-cg") {
                    // C5W-CH2E-CH1E
                    bond_angle  = 113.6;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD1) {
                    // "cb-cg-cd1") {
                    // CW-C5W-CH2E
                    bond_angle  = 126.8;
                    *std_dev = 0.020;
               // } else if (atom1==CB && atom2==CG && atom3==CD2) {
               } else if (atom1==CD1 && atom2==CG && atom3==CD2) {
                    // "cb-cg-cd2") {
                    // CH2E-C5W-CR1E
                    // bondAngle  = 126.9;
                    bond_angle  = 106.3;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD1 && atom3==NE1) {
                    // "cg-cd-ne") {
                    // C5W-CR1E-NH1
                    bond_angle  = 110.2;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD2 && atom3==CE2) {
                    // "cg-cd-ce") {
                    // C5W-CW-CW
                    bond_angle  = 107.2;
                    *std_dev = 0.020;
               // } else if (atom1==CG && atom2==CD2 && atom3==CE3) {
               } else if (atom1==CE2 && atom2==CD2 && atom3==CE3) {
                    // "cg-cd-ce3") {
                    // C5W-CW-CR1E
                    // bondAngle  = 133.9;
                    bond_angle  = 118.9; // 360 - 133.0 - 107.2
                    *std_dev = 0.020;
               } else if (atom1==CD2 && atom2==CE2 && atom3==CZ2) {
                    // "cd-ce-cz2") {
                    // CW-CW-CR1W
                    bond_angle  = 122.4;
                    *std_dev = 0.020;
               }  else if (atom1==CD2 && atom2==CE3 && atom3==CZ3) {
                    // "cd-ce-cz3") {
                    // CW-CR1E-CR1E
                    bond_angle  = 118.6;
                    *std_dev = 0.020;
               } else if (atom1==CE3 && atom2==CZ3 && atom3==CH2) {
                    // "ce3-cz3-ch") {
                    // CR1E-CR1E-CR1W
                    bond_angle  = 121.1;
                    *std_dev = 0.020;
               }
               break;
          case TYR:                                        // Tyrosine
               if (atom1 == CA && atom2==CB && atom3==CG) {
                    // "ca-cb-cg") {
                    // CY-CH2E-CH1E
                    bond_angle  = 113.9;
                    *std_dev = 0.020;
               } else if (atom1==CB && atom2==CG && atom3==CD1) {
                    // "cb-cg-cd") {
                    // CH2E-CY-CR1E
                    bond_angle  = 120.8;
                    *std_dev = 0.020;
               // } else if (atom1==CB && atom2==CG && atom3==CD2) {
               } else if (atom1==CD1 && atom2==CG && atom3==CD2) {
                    // "cb-cg-cd") {
                    // CH2E-CY-CR1E
                    // bondAngle  = 120.8;
                    bond_angle  = 118.4;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD1 && atom3==CE1) {
                    // "cg-cd-ce") {
                    // CY-CR1E-CR1E
                    bond_angle  = 121.2;
                    *std_dev = 0.020;
               } else if (atom1==CG && atom2==CD2 && atom3==CE2) {
                    // "cg-cd-ce") {
                    // CY-CR1E-CR1E
                    bond_angle  = 121.2;
                    *std_dev = 0.020;
               } else if (atom1==CD1 && atom2==CE1 && atom3==CZ) {
                    // "cd-ce-cz") {
                    // CY2-CR1E-CR1E
                    bond_angle  = 119.6;
                    *std_dev = 0.020;
               } else if (atom1==CD2 && atom2==CE2 && atom3==CZ) {
                    // "cd-ce-cz") {
                    // CY2-CR1E-CR1E
                    bond_angle  = 119.6;
                    *std_dev = 0.020;
               } else if (atom1==CE1 && atom2==CZ && atom3==OH) {
                    // "ce-cz-oh") {
                    // CR1E-CY2-OH1
                    bond_angle  = 119.9;
                    *std_dev = 0.020;
               }
               break;
          default:
               break;
          }
     }

     // Pseudo-sidechain dihedral values.
     // These are not actual constants but serve as initial values
     // The values reflect the most probable rotamer states
     if (!is_initialized(bond_angle)) {
          if (atom3==PS) {
               switch (res) {  
               case ALA:           
                    bond_angle = 89.8283;
                    *std_dev = 0.0;
                    break;
               case CYS:           
                    bond_angle = 118.8;
                    *std_dev = 0.0;
                    break;
               case ASP:           
                    bond_angle = 117.5;
                    *std_dev = 0.0;
                    break;
               case GLU:           
                    bond_angle = 108.9;
                    *std_dev = 0.0;
                    break;
               case PHE:
                    bond_angle = 126.1;
                    *std_dev = 0.0;
                    break;
               case GLY:           
                    bond_angle = 89.8283;
                    *std_dev = 0.0;
                    break;
               case HIS:           
                    bond_angle = 126.1;
                    *std_dev = 0.0;
                    break;
               case ILE:
                    bond_angle = 100.3;
                    *std_dev = 0.0;
                    break;
               case LYS:
                    bond_angle = 111.7;
                    *std_dev = 0.0;
                    break;
               case LEU:
                    bond_angle = 108.9;
                    *std_dev = 0.0;
                    break;
               case MET:
                    bond_angle = 103.1;
                    *std_dev = 0.0;
                    break;
               case ASN:
                    bond_angle = 118.8;
                    *std_dev = 0.0;
                    break;
               case PRO:
                    bond_angle = 126.1;
                    *std_dev = 0.0;
                    break;
               case GLN:
                    bond_angle = 106.0;
                    *std_dev = 0.0;
                    break;
               case ARG:
                    bond_angle = 114.6;
                    *std_dev = 0.0;
                    break;
               case SER:
                    bond_angle = 91.7;
                    *std_dev = 0.0;
                    break;
               case THR:
                    bond_angle = 94.5;
                    *std_dev = 0.0;
                    break;
               case VAL:
                    bond_angle = 88.8;
                    *std_dev = 0.0;
                    break;
               case TRP:
                    bond_angle = 108.9;
                    *std_dev = 0.0;
                    break;
               case TYR:
                    bond_angle = 126.1;
                    *std_dev = 0.0;
                    break;
               default:
                    break;
               
               }
          }
     }
     
     // Geometric constants for Phaistos-specific atom placement strategy

     // Residue-unspecific constants
     if (!is_initialized(bond_angle)) {
          if (atom1==CA && atom2==N && atom3==O) {             // O 
               bond_angle = 90.0;
               *std_dev = 0.0;
          } else if (atom1==N && atom2==C && atom3==CB) {      // CB
               bond_angle = 89.8283;
               *std_dev = 1.3006;
          } else if (atom1==CA && atom2==CA && atom3==CB) {    // CB - in CA-only representation
               bond_angle = 96.1;
               *std_dev = 0.0;
          }
     }


     // Hydrogens
     if (!is_initialized(bond_angle)) {

          if (atom1==CA && atom2==N && (atom3==H1 || atom3==H2 || atom3==H3)) {                     // H1,H2,H3
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CA && atom2==CB && (atom3==HB1 || atom3==HB2 || atom3==HB3)) {          // ALA: HB1, HB2, HB3
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CB && atom2==SG && atom3==HG) {                                         // CYS: HG
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CG1 && atom2==CD1 && (atom3==HD11 || atom3==HD12 || atom3==HD13)) {     // ILE: HD11, HD12, HD13
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CB && atom2==CG2 && (atom3==HG21 || atom3==HG22 || atom3==HG23)) {      // ILE,THR,VAL: HG21, HG22, HG23
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CE && atom2==NZ && (atom3==HZ1 || atom3==HZ2 || atom3==HZ3)) {          // LYS: HZ1, HZ2, HZ3
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CG && atom2==CD1 && (atom3==HD11 || atom3==HD12 || atom3==HD13)) {      // LEU: HD11, HD12, HD13
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CG && atom2==CD2 && (atom3==HD21 || atom3==HD22 || atom3==HD23)) {      // LEU: HD21, HD22, HD23
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==SD && atom2==CE && (atom3==HE1 || atom3==HE2 || atom3==HE3)) {          // MET: HE1, HE2, HE3
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CG && atom2==ND2 && (atom3==HD21 || atom3==HD22)) {                     // ASN: HD21, HD22
               bond_angle = 120.0;
               *std_dev = 0.0;
          } else if (atom1==CD && atom2==NE2 && (atom3==HE21 || atom3==HE22)) {                     // GLN: HE21, HE22
               bond_angle = 120.0;
               *std_dev = 0.0;
          } else if (atom1==CZ && atom2==NH1 && (atom3==HH11 || atom3==HH12)) {                     // ARG: HH11, HH12
               bond_angle = 120.0;
               *std_dev = 0.0;
          } else if (atom1==CZ && atom2==NH2 && (atom3==HH21 || atom3==HH22)) {                     // ARG: HH21, HH22
               bond_angle = 120.0;
               *std_dev = 0.0;
          } else if (atom1==CB && atom2==OG && atom3==HG) {                                         // SER: HG
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CB && atom2==OG1 && atom3==HG1) {                                       // THR: HG1
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CB && atom2==CG1 && (atom3==HG11 || atom3==HG12 || atom3==HG13)) {      // VAL: HG11, HG12, HG13
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;
          } else if (atom1==CZ && atom2==OH && atom3==HH) {                                         // TYR: HH
               bond_angle = acos(-1.0/3.0)/M_PI*180; // tetrahedral
               *std_dev = 0.0;

          // Default hydrogen bondAngle is 90 degrees 
          } else if (atom3>=H and atom3<=H3) {
               bond_angle = 90.0;
               *std_dev = 0.0;
          }            
     }
          

     double conversion_factor = M_PI/180;

     bond_angle *= conversion_factor;
     *std_dev *= conversion_factor;

     assert(is_initialized(bond_angle));
     
     return bond_angle;
}

}

#endif
