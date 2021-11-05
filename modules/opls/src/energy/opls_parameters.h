// opls_parameters.h --- Default parameters for OPLS energy class
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson, Wouter Boomsma
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

// The original version of this parameter file was developed in the
// group of Prof. Karl F. Freed, Department of Chemistry, University
// of Chicago by modifying the TINKER OPLS-AA parameters as per the
// Kaminsky, et al. article cited below

#ifndef PARAMETERS_OPLS_H
#define PARAMETERS_OPLS_H

const static std::string parameters_opls = " \n\
\n\
      ##############################\n\
      ##                          ##\n\
      ##  Force Field Definition  ##\n\
      ##                          ##\n\
      ##############################\n\
\n\
\n\
forcefield              OPLS-AA/L\n\
\n\
vdwtype                 LENNARD-JONES\n\
radiusrule              GEOMETRIC\n\
radiustype              SIGMA\n\
radiussize              DIAMETER\n\
epsilonrule             GEOMETRIC\n\
torsionunit             0.5\n\
vdw-14-scale            2.0\n\
chg-14-scale            2.0\n\
dielectric              1.0\n\
\n\
\n\
      #############################\n\
      ##                         ##\n\
      ##  Literature References  ##\n\
      ##                         ##\n\
      #############################\n\
\n\
\n\
G. A. Kaminsky, R. A. Friesner, J. Tirado-Rives and W. L. Jorgensen,\n\
\"Evaluation and Reparametrization of the OPLS-AA Force Field for\n\
Proteins via Comparison with Accurate Quantum Chemical Calculations\n\
on Peptides\", J. Phys. Chem. B, 105, 6474-6487 (2001)\n\
\n\
W. L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, \"Development\n\
and Testing of the OPLS All-Atom Force Field on Conformational\n\
Energetics and Properties of Organic Liquids\", J. Am. Chem. Soc.,\n\
117, 11225-11236 (1996)\n\
\n\
The original version of this parameter file was developed in the\n\
group of Prof. Karl F. Freed, Department of Chemistry, University\n\
of Chicago by modifying the TINKER OPLS-AA parameters as per the\n\
Kaminsky, et al. article cited above\n\
\n\
\n\
      #############################\n\
      ##                         ##\n\
      ##  Atom Type Definitions  ##\n\
      ##                         ##\n\
      #############################\n\
\n\
\n\
atom      1     1    CT      \"RCH3 Alkane\"               6     12.000     4\n\
atom      2     1    CT      \"R2CH2 Alkane\"              6     12.000     4\n\
atom      3     1    CT      \"R3CH Alkane\"               6     12.000     4\n\
atom      4     1    CT      \"CH4 Methane\"               6     12.000     4\n\
atom      5     1    CT      \"R4C Alkane\"                6     12.000     4\n\
atom      6     2    HC      \"HR Alkane\"                 1      1.008     1\n\
atom      7     3    CM      \"R2C=C Alkene\"              6     12.000     3\n\
atom      8     3    CM      \"RHC=C Alkene\"              6     12.000     3\n\
atom      9     3    CM      \"H2C=C Alkene\"              6     12.000     3\n\
atom     10     4    HC      \"H-C=C Alkene\"              1      1.008     1\n\
atom     11     5    CA      \"CH Benzene\"                6     12.000     3\n\
atom     12     6    HA      \"H Benzene\"                 1      1.008     1\n\
atom     13     5    CA      \"C Aromatic Fusion\"         6     12.000     3\n\
atom     14     1    CT      \"CH3 Toluene\"               6     12.000     4\n\
atom     15     1    CT      \"CH2 Ethyl Benzene\"         6     12.000     4\n\
atom     16     7    OH      \"OH Alcohol\"                8     15.999     2\n\
atom     17     8    HO      \"HO Alcohol\"                1      1.008     1\n\
atom     18     2    HC      \"CH3 Methanol\"              1      1.008     1\n\
atom     19     1    CT      \"CH2 & CH3 Alcohol\"         6     12.000     4\n\
atom     20     1    CT      \"CH i-Pr Alcohol\"           6     12.000     4\n\
atom     21     1    CT      \"C t-Bu Alcohol\"            6     12.000     4\n\
atom     22     1    CT      \"CH2 Trifluoroethanol\"      6     12.000     4\n\
atom     23     9    CT      \"CF3 Trifluoroethanol\"      6     12.000     4\n\
atom     24    10    OH      \"OH Trifluoroethanol\"       8     15.999     2\n\
atom     25     8    HO      \"HO Trifluoroethanol\"       1      1.008     1\n\
atom     26    11    F       \"F Trifluoroethanol\"        9     18.998     1\n\
atom     27     2    HC      \"CH2 Trifluoroethanol\"      1      1.008     1\n\
atom     28     5    CA      \"COH Phenol\"                6     12.000     3\n\
atom     29    10    OH      \"OH Phenol\"                 8     15.999     2\n\
atom     30     8    HO      \"HO Phenol\"                 1      1.008     1\n\
atom     31    10    OH      \"OH Diols\"                  8     15.999     2\n\
atom     32     8    HO      \"HO Diols\"                  1      1.008     1\n\
atom     33    12    OS      \"O Ether\"                   8     15.999     2\n\
atom     34     1    CT      \"CH3OR Me Ether\"            6     12.000     4\n\
atom     35     1    CT      \"RCH2OR Et Ether\"           6     12.000     4\n\
atom     36     1    CT      \"R2CHOR i-Pr Ether\"         6     12.000     4\n\
atom     37     1    CT      \"R3COR t-Bu Ether\"          6     12.000     4\n\
atom     38     2    HC      \"H Ether Calpha\"            1      1.008     1\n\
atom     39    13    OS      \"O Acetal\"                  8     15.999     2\n\
atom     40    10    OH      \"OH Hemiacetal\"             8     15.999     2\n\
atom     41     8    HO      \"HO Hemiacetal\"             1      1.008     1\n\
atom     42    14    CO      \"CH2O2 Acetal\"              6     12.000     4\n\
atom     43     2    HC      \"CH2O2 Acetal\"              1      1.008     1\n\
atom     44    14    CO      \"CH2O2 Hemiacetal\"          6     12.000     4\n\
atom     45     2    HC      \"CH2O2 Hemiacetal\"          1      1.008     1\n\
atom     46    14    CO      \"CHRO2 Acetal\"              6     12.000     4\n\
atom     47     2    HC      \"CHRO2 Acetal\"              1      1.008     1\n\
atom     48    14    CO      \"CHRO2 Hemiacetal\"          6     12.000     4\n\
atom     49     2    HC      \"CHRO2 Hemiacetal\"          1      1.008     1\n\
atom     50    14    CO      \"CR2O2 Acetal\"              6     12.000     4\n\
atom     51    14    CO      \"CR2O2 Hemiacetal\"          6     12.000     4\n\
atom     52    15    SH      \"SH Thiol\"                 16     32.066     2\n\
atom     53    16    SH      \"H2S Hydrogen Sulfide\"     16     32.066     2\n\
atom     54    17    S       \"S Sulfide\"                16     32.066     2\n\
atom     55    17    S       \"S Disulfide\"              16     32.066     2\n\
atom     56    18    HS      \"HS Thiol\"                  1      1.008     1\n\
atom     57    18    HS      \"H2S Hydrogen Sulfide\"      1      1.008     1\n\
atom     58     1    CT      \"CH3SH Me Thiol\"            6     12.000     4\n\
atom     59     1    CT      \"RCH2SH Et Thiol\"           6     12.000     4\n\
atom     60     1    CT      \"R2CHSH i-Pr Thiol\"         6     12.000     4\n\
atom     61     1    CT      \"R3CSH t-Bu Thiol\"          6     12.000     4\n\
atom     62     1    CT      \"CH3SR Me Sulfide\"          6     12.000     4\n\
atom     63     1    CT      \"RCH2SR Et Sulfide\"         6     12.000     4\n\
atom     64     1    CT      \"R2CHSR i-Pr Sulfide\"       6     12.000     4\n\
atom     65     1    CT      \"R3CSR t-Bu Sulfide\"        6     12.000     4\n\
atom     66     1    CT      \"CH3SSR Me Disulf\"          6     12.000     4\n\
atom     67     1    CT      \"RCH2SSR Et Disulf\"         6     12.000     4\n\
atom     68     1    CT      \"R2CHSSR i-Pr Disulf\"       6     12.000     4\n\
atom     69     1    CT      \"R3CSSR t-Bu Disulf\"        6     12.000     4\n\
atom     70    19    NT      \"NH2 Primary Amine\"         7     14.007     3\n\
atom     71    20    H2      \"HN Primary Amine\"          1      1.008     1\n\
atom     72     1    CT      \"CH3NH2 Amine\"              6     12.000     4\n\
atom     73     1    CT      \"RCH2NH2 & GLY CA\"          6     12.000     4\n\
atom     74     1    CT      \"R2CHNH2 & ALA CA\"          6     12.000     4\n\
atom     75     1    CT      \"R3CNH2 & AIB CA\"           6     12.000     4\n\
atom     76    19    NT      \"NH2 Aniline\"               7     14.007     3\n\
atom     77    20    H2      \"HN Aniline\"                1      1.008     1\n\
atom     78     5    CA      \"CNH2 Aniline\"              6     12.000     4\n\
atom     79    19    NT      \"NR3 Tertiary Amine\"        7     14.007     3\n\
atom     80     1    CT      \"CH3NR2 Tert Amine\"         6     12.000     4\n\
atom     81     1    CT      \"RCH2NR2 Tert Amine\"        6     12.000     4\n\
atom     82    21    C       \"NC=O Amide\"                6     12.000     3\n\
atom     83    22    O       \"NC=O Amide\"                8     15.999     1\n\
atom     84    23    N       \"H2NC=O Amide\"              7     14.007     3\n\
atom     85    23    N       \"RHNC=O Amide\"              7     14.007     3\n\
atom     86    23    N       \"R2NC=O Amide\"              7     14.007     3\n\
atom     87    24    H       \"H2NC=O Amide\"              1      1.008     1\n\
atom     88    24    H       \"RHNC=O Amide\"              1      1.008     1\n\
atom     89     1    CT      \"O=CNH-Me Amide\"            6     12.000     4\n\
atom     90     1    CT      \"O=CNR-Me Amide\"            6     12.000     4\n\
atom     91     1    CT      \"O=CNH-CH2R Amide\"          6     12.000     4\n\
atom     92     1    CT      \"O=CNR-CH2R Amide\"          6     12.000     4\n\
atom     93     1    CT      \"O=CNR-CHR2 Amide\"          6     12.000     4\n\
atom     94    21    C       \"C Urea\"                    6     12.000     3\n\
atom     95    22    O       \"O Urea\"                    8     15.999     1\n\
atom     96    23    N       \"NH2 Urea\"                  7     14.007     3\n\
atom     97    24    H       \"HN Urea\"                   1      1.008     1\n\
atom     98    21    C       \"O=CNHC=O Imide\"            6     12.000     3\n\
atom     99    22    O       \"O=CNHC=O Imide\"            8     15.999     1\n\
atom    100    23    N       \"O=CNHC=O Imide\"            7     14.007     3\n\
atom    101    24    H       \"O=CNHC=O Imide\"            1      1.008     1\n\
atom    102    25    HC      \"H-C Formimide\"             1      1.008     1\n\
atom    103     1    CT      \"CH3 Imide\"                 6     12.000     4\n\
atom    104     1    CT      \"RCH2 Imide\"                6     12.000     4\n\
atom    105     1    CT      \"R2CH Imide\"                6     12.000     4\n\
atom    106     1    CT      \"R3C Imide\"                 6     12.000     4\n\
atom    107    21    C       \"COOH Carboxylic Acid\"      6     12.000     3\n\
atom    108    26    OH      \"OH Carboxylic Acid\"        8     15.999     2\n\
atom    109    22    O       \"C=O Carboxylic Acid\"       8     15.999     1\n\
atom    110     8    HO      \"COOH Carboxylic Acid\"      1      1.008     1\n\
atom    111    27    C       \"COO- Carboxylate\"          6     12.000     3\n\
atom    112    28    O2      \"COO- Carboxylate\"          8     15.999     1\n\
atom    113     1    CT      \"CH3COO- Carboxylate\"       6     12.000     4\n\
atom    114     1    CT      \"RCH2COO- Carboxylate\"      6     12.000     4\n\
atom    115     1    CT      \"R2CHCOO- Carboxylate\"      6     12.000     4\n\
atom    116     1    CT      \"R3CCOO- Carboxylate\"       6     12.000     4\n\
atom    117    21    C       \"HC=O Aldehyde\"             6     12.000     3\n\
atom    118    22    O       \"HC=O Aldehyde\"             8     15.999     1\n\
atom    119    29    HC      \"HC=O Aldehyde\"             1      1.008     1\n\
atom    120    21    C       \"C=O Ketone\"                6     12.000     3\n\
atom    121    22    O       \"C=O Ketone\"                8     15.999     1\n\
atom    122    29    HC      \"HC-C=O Halpha\"             1      1.008     1\n\
atom    123    30    N3      \"NH4+ Ammonium\"             7     14.007     4\n\
atom    124    30    N3      \"RNH3+ Ammonium\"            7     14.007     4\n\
atom    125    30    N3      \"R2NH2+ Ammonium\"           7     14.007     4\n\
atom    126    30    N3      \"R3NH+ Ammonium\"            7     14.007     4\n\
atom    127    30    N3      \"R4N+ Ammonium\"             7     14.007     4\n\
atom    128    31    H3      \"NH4+ Ammonium\"             1      1.008     1\n\
atom    129    31    H3      \"RNH3+ Ammonium\"            1      1.008     1\n\
atom    130    31    H3      \"R2NH2+ Ammonium\"           1      1.008     1\n\
atom    131    31    H3      \"R3NH+ Ammonium\"            1      1.008     1\n\
atom    132     1    CT      \"CH3NH3+ Ammonium\"          6     12.000     4\n\
atom    133     1    CT      \"RCH2NH3+ Ammonium\"         6     12.000     4\n\
atom    134     1    CT      \"R2CHNH3+ Ammonium\"         6     12.000     4\n\
atom    135     1    CT      \"R3CNH3+ Ammonium\"          6     12.000     4\n\
atom    136    32    N2      \"NH2 Guanidinium\"           7     14.007     3\n\
atom    137    32    N2      \"NHR Guanidinium\"           7     14.007     3\n\
atom    138    31    H3      \"NH2 Guanidinium\"           1      1.008     1\n\
atom    139    31    H3      \"NHR Guanidinium\"           1      1.008     1\n\
atom    140    33    CA      \"C+ Guanidinium\"            6     12.000     3\n\
atom    141     1    CT      \"CH3 Me Guanidinium\"        6     12.000     4\n\
atom    142     1    CT      \"CH3 Et Guanidinium\"        6     12.000     4\n\
atom    143     1    CT      \"CH2 Et Guanidinium\"        6     12.000     4\n\
atom    144     1    CT      \"CH2 ARG CG\"                6     12.000     4\n\
atom    145    34    C*      \"C TRP CG\"                  6     12.000     3\n\
atom    146    35    CB      \"C TRP CD2\"                 6     12.000     3\n\
atom    147    36    CN      \"C TRP CE2\"                 6     12.000     3\n\
atom    148    37    NA      \"NH TRP NE1 & HID/E\"        7     14.007     3\n\
atom    149    24    H       \"HN TRP HE1 & HID/E\"        1      1.008     1\n\
atom    150     1    CT      \"CH2 HIS CB\"                6     12.000     4\n\
atom    151    38    CP      \"CH HID/HIE CE1\"            6     12.000     3\n\
atom    152    39    CV      \"C HID CD2 & HIE CG\"        6     12.000     3\n\
atom    153    40    CW      \"C HID CG & HIE CD2\"        6     12.000     3\n\
atom    154    38    CP      \"CH HIP CE1\"                6     12.000     3\n\
atom    155    40    CW      \"C HIP CG & CD2\"            6     12.000     3\n\
atom    156    41    NB      \"N HID NE2 & HIE ND1\"       7     14.007     2\n\
atom    157    37    NA      \"NH HIP ND1 & NE2\"          7     14.007     3\n\
atom    158    24    H       \"HN HIP ND1 & NE2\"          1      1.008     1\n\
atom    159    42    CT      \"CH2 PRO CD\"                6     12.000     4\n\
atom    160     1    CT      \"CH2 C-term GLY CA\"         6     12.000     4\n\
atom    161     1    CT      \"CH C-term ALA CA\"          6     12.000     4\n\
atom    162     1    CT      \"C C-term AIB CA\"           6     12.000     4\n\
atom    163     1    CT      \"CH C-term PRO CA\"          6     12.000     4\n\
atom    164    43    OW      \"O Water (TIP3P)\"           8     15.999     2\n\
atom    165    44    HW      \"H Water (TIP3P)\"           1      1.008     1\n\
atom    166    45    OW      \"O Water (SPC)\"             8     15.999     2\n\
atom    167    44    HW      \"H Water (SPC)\"             1      1.008     1\n\
atom    168    46    NT      \"N Ammonia\"                 7     14.007     3\n\
atom    169    20    H2      \"H Ammonia\"                 1      1.008     1\n\
atom    170    47    O-      \"O- Hydroxide Ion\"          8     15.999     1\n\
atom    171     8    HO      \"H Hydroxide Ion\"           1      1.008     1\n\
atom    172    48    Li+     \"Li+ Lithium Ion\"           3      6.941     0\n\
atom    173    49    Na+     \"Na+ Sodium Ion\"           11     22.990     0\n\
atom    174    50    K+      \"K+ Potassium Ion\"         19     39.098     0\n\
atom    175    51    Rb+     \"Rb+ Rubidium Ion\"         37     85.468     0\n\
atom    176    52    Cs+     \"Cs+ Cesium Ion\"           55    132.905     0\n\
atom    177    53    Mg+     \"Mg+2 Magnesium Ion\"       12     24.305     0\n\
atom    178    54    Ca+     \"Ca+2 Calcium Ion\"         20     40.078     0\n\
atom    179    55    Sr+     \"Sr+2 Strontium Ion\"       38     87.620     0\n\
atom    180    56    Ba+     \"Ba+2 Barium Ion\"          56    137.327     0\n\
atom    181    57    F-      \"F- Fluoride Ion\"           9     18.998     0\n\
atom    182    58    Cl-     \"Cl- Chloride Ion\"         17     35.453     0\n\
atom    183    59    Br-     \"Br- Bromide Ion\"          35     79.904     0\n\
atom    184    60    He      \"Helium Atom\"               2      4.003     0\n\
atom    185    61    Ne      \"Neon Atom\"                10     20.179     0\n\
atom    186    62    Ar      \"Argon Atom\"               18     39.948     0\n\
atom    187    63    Kr      \"Krypton Atom\"             36     83.800     0\n\
atom    188    64    Xe      \"Xenon Atom\"               54    131.300     0\n\
atom    189    65    CT      \"CB of Serine\"              6     12.000     4\n\
atom    190    66    CT      \"CB of Cysteine\"            6     12.000     4\n\
atom    191    67    CT      \"CB of Asparagine\"          6     12.000     4\n\
atom    192    68    CT      \"CB of Glutamine\"           6     12.000     4\n\
atom    193    69    CT      \"CB of Leucine\"             6     12.000     4\n\
atom    194    69    CT      \"CB of Valine\"              6     12.000     4\n\
atom    195    70    CT      \"CB of Isoleucine\"          6     12.000     4\n\
atom    196    71    CT      \"CB of Methionine\"          6     12.000     4\n\
atom    197    72    CT      \"CB of Tryptophan\"          6     12.000     4\n\
atom    198    73    CT      \"CB of Aspartic Acid\"       6     12.000     4\n\
atom    199    74    CT      \"CB of Glutamic Acid\"       6     12.000     4\n\
atom    200    75    CT      \"CB of Lysine\"              6     12.000     4\n\
atom    201    76    CT      \"CB of Histidine (+)\"       6     12.000     4\n\
atom    202    77    CT      \"CB of Arginine\"            6     12.000     4\n\
\n\
\n\
      ################################\n\
      ##                            ##\n\
      ##  Van der Waals Parameters  ##\n\
      ##                            ##\n\
      ################################\n\
\n\
\n\
vdw          1              3.5000     0.0660\n\
vdw          2              2.5000     0.0300\n\
vdw          3              3.5500     0.0760\n\
vdw          4              2.4200     0.0300\n\
vdw          5              3.5500     0.0700\n\
vdw          6              2.4200     0.0300\n\
vdw          7              3.1200     0.1700\n\
vdw          8              0.0000     0.0000\n\
vdw          9              3.2500     0.0620\n\
vdw         10              3.0700     0.1700\n\
vdw         11              2.9400     0.0610\n\
vdw         12              2.9000     0.1400\n\
vdw         13              3.0000     0.1400\n\
vdw         14              3.5000     0.0660\n\
vdw         15              3.6000     0.4250\n\
vdw         16              3.7000     0.2500\n\
vdw         17              3.6000     0.3550\n\
vdw         18              0.0000     0.0000\n\
vdw         19              3.2500     0.1700\n\
vdw         20              0.0000     0.0000\n\
vdw         21              3.7500     0.1050\n\
vdw         22              2.9600     0.2100\n\
vdw         23              3.2500     0.1700\n\
vdw         24              0.0000     0.0000\n\
vdw         25              2.5000     0.0200\n\
vdw         26              3.0000     0.1700\n\
vdw         27              3.7500     0.1050\n\
vdw         28              2.9600     0.2100\n\
vdw         29              2.4200     0.0150\n\
vdw         30              3.2500     0.1700\n\
vdw         31              0.0000     0.0000\n\
vdw         32              3.2500     0.1700\n\
vdw         33              2.2500     0.0500\n\
vdw         34              3.5500     0.0700\n\
vdw         35              3.5500     0.0700\n\
vdw         36              3.5500     0.0700\n\
vdw         37              3.2500     0.1700\n\
vdw         38              3.5500     0.0700\n\
vdw         39              3.5500     0.0700\n\
vdw         40              3.5500     0.0700\n\
vdw         41              3.2500     0.1700\n\
vdw         42              3.5000     0.0660\n\
vdw         43         3.150656111     0.152072595\n\
vdw         44              0.0000     0.0000\n\
vdw         45         3.165555296     0.155406042\n\
vdw         46              3.3600     0.2100\n\
vdw         47              3.2000     0.2500\n\
vdw         48            2.126452     0.018279\n\
vdw         49            3.330445     0.002772\n\
vdw         50            4.934628     0.000328\n\
vdw         51            5.621773     0.000171\n\
vdw         52            6.715999     0.000081\n\
vdw         53            1.644471     0.875044\n\
vdw         54            2.412031     0.449657\n\
vdw         55            3.102688     0.118226\n\
vdw         56            3.816610     0.047096\n\
vdw         57             2.73295     0.7200\n\
vdw         58             4.41724     0.1180\n\
vdw         59             4.62376     0.0900\n\
vdw         60              2.5560     0.0200\n\
vdw         61              2.7800     0.0690\n\
vdw         62              3.4010     0.2339\n\
vdw         63              3.6240     0.3170\n\
vdw         64              3.9350     0.4330\n\
vdw         65              3.5000     0.0660\n\
vdw         66              3.5000     0.0660\n\
vdw         67              3.5000     0.0660\n\
vdw         68              3.5000     0.0660\n\
vdw         69              3.5000     0.0660\n\
vdw         70              3.5000     0.0660\n\
vdw         71              3.5000     0.0660\n\
vdw         72              3.5000     0.0660\n\
vdw         73              3.5000     0.0660\n\
vdw         74              3.5000     0.0660\n\
vdw         75              3.5000     0.0660\n\
vdw         76              3.5000     0.0660\n\
vdw         77              3.5000     0.0660\n\
\n\
\n\
      ##################################\n\
      ##                              ##\n\
      ##  Bond Stretching Parameters  ##\n\
      ##                              ##\n\
      ##################################\n\
\n\
\n\
bond         1    1          268.0     1.5290\n\
bond         1    2          340.0     1.0900\n\
bond         1    5          317.0     1.5100\n\
bond         1    7          320.0     1.4100\n\
bond         1   15          222.0     1.8100\n\
bond         1   17          222.0     1.8100\n\
bond         1   21          317.0     1.5220\n\
bond         1   23          337.0     1.4490\n\
bond         1   27          317.0     1.5220\n\
bond         1   30          367.0     1.4710\n\
bond         1   32          337.0     1.4630\n\
bond         1   34          317.0     1.4950\n\
bond         1   39          317.0     1.5040\n\
bond         1   40          317.0     1.5040\n\
bond         1   42          268.0     1.5290\n\
bond         1   65          268.0     1.5290\n\
bond         1   66          268.0     1.5290\n\
bond         1   67          268.0     1.5290\n\
bond         1   68          268.0     1.5290\n\
bond         1   69          268.0     1.5290\n\
bond         1   70          268.0     1.5290\n\
bond         1   71          268.0     1.5290\n\
bond         1   72          268.0     1.5290\n\
bond         1   73          268.0     1.5290\n\
bond         1   74          268.0     1.5290\n\
bond         1   75          268.0     1.5290\n\
bond         1   76          268.0     1.5290\n\
bond         1   77          268.0     1.5290\n\
bond         2   42          340.0     1.0900\n\
bond         2   65          340.0     1.0900\n\
bond         2   66          340.0     1.0900\n\
bond         2   67          340.0     1.0900\n\
bond	     2   68          340.0     1.0900\n\
bond         2   69          340.0     1.0900\n\
bond         2   70          340.0     1.0900\n\
bond         2   71          340.0     1.0900\n\
bond         2   72          340.0     1.0900\n\
bond         2   73          340.0     1.0900\n\
bond         2   74          340.0     1.0900\n\
bond         2   75          340.0     1.0900\n\
bond         2   76          340.0     1.0900\n\
bond         2   77          340.0     1.0900\n\
bond         5    5          469.0     1.4000\n\
bond         5    6          367.0     1.0800\n\
bond         5   10          450.0     1.3640\n\
bond         5   34          546.0     1.3520\n\
bond         5   35          469.0     1.4000\n\
bond         5   36          469.0     1.4000\n\
bond         5   37          427.0     1.3810\n\
bond         6   38          367.0     1.0800\n\
bond         6   39          367.0     1.0800\n\
bond         6   40          367.0     1.0800\n\
bond         7    8          553.0     0.9450\n\
bond         7   65          320.0     1.4100\n\
bond         8   10          553.0     0.9450\n\
bond        15   18          274.0     1.3360\n\
bond        15   66          222.0     1.8100\n\
bond        17   17          166.0     2.0380\n\
bond        21   22          570.0     1.2290\n\
bond        21   23          490.0     1.3350\n\
bond        21   29          340.0     1.0900\n\
bond        21   67          317.0     1.5220\n\
bond        23   24          434.0     1.0100\n\
bond        23   42          337.0     1.4490\n\
bond        24   37          434.0     1.0100\n\
bond        27   28          656.0     1.2500\n\
bond        27   73          317.0     1.5220\n\
bond        30   31          434.0     1.0100\n\
bond        31   32          434.0     1.0100\n\
bond        32   33          481.0     1.3400\n\
bond        34   35          388.0     1.4590\n\
bond        34   72          317.0     1.4950\n\
bond        35   36          447.0     1.4190\n\
bond        36   37          428.0     1.3800\n\
bond        37   38          477.0     1.3430\n\
bond        37   40          427.0     1.3810\n\
bond        38   41          488.0     1.3350\n\
bond        39   40          518.0     1.3710\n\
bond        39   41          410.0     1.3940\n\
bond        40   40          518.0     1.3710\n\
bond        40   76          317.0     1.5040\n\
bond        43   44          529.6     0.9572\n\
bond        44   45          527.2     1.0000\n\
\n\
\n\
      ################################\n\
      ##                            ##\n\
      ##  Angle Bending Parameters  ##\n\
      ##                            ##\n\
      ################################\n\
\n\
\n\
angle        1    1    1     58.35     112.70\n\
angle        1    1    2     37.50     110.70\n\
angle        1    1    5     63.00     114.00\n\
angle        1    1    7     50.00     109.50\n\
angle        1    1   15     50.00     108.60\n\
angle        1    1   17     50.00     114.70\n\
angle        1    1   21     63.00     111.10\n\
angle        1    1   23     80.00     109.70\n\
angle        1    1   27     63.00     111.10\n\
angle        1    1   30     80.00     111.20\n\
angle        1    1   32     80.00     111.20\n\
angle        1    1   34     63.00     115.60\n\
angle        1    1   39     63.00     114.00\n\
angle        1    1   40     63.00     114.00\n\
angle        1    1   42     58.35     112.70\n\
angle        1    1   69     58.35     112.70\n\
angle        1    1   70     58.35     112.70\n\
angle        1    1   75     58.35     112.70\n\
angle        1    1   77     58.35     112.70\n\
angle        2    1    2     33.00     107.80\n\
angle        2    1    5     35.00     109.50\n\
angle        2    1    7     35.00     109.50\n\
angle        2    1   15     35.00     109.50\n\
angle        2    1   17     35.00     109.50\n\
angle        2    1   21     35.00     109.50\n\
angle        2    1   23     35.00     109.50\n\
angle        2    1   27     35.00     109.50\n\
angle        2    1   30     35.00     109.50\n\
angle        2    1   32     35.00     109.50\n\
angle        2    1   34     35.00     109.50\n\
angle        2    1   39     35.00     109.50\n\
angle        2    1   40     35.00     109.50\n\
angle        2    1   42     37.50     110.70\n\
angle        2    1   65     37.50     110.70\n\
angle        2    1   66     37.50     110.70\n\
angle        2    1   67     37.50     110.70\n\
angle        2    1   68     37.50     110.70\n\
angle        2    1   69     37.50     110.70\n\
angle        2    1   70     37.50     110.70\n\
angle        2    1   71     37.50     110.70\n\
angle        2    1   72     37.50     110.70\n\
angle        2    1   73     37.50     110.70\n\
angle        2    1   74     37.50     110.70\n\
angle        2    1   75     37.50     110.70\n\
angle        2    1   76     37.50     110.70\n\
angle        2    1   77     37.50     110.70\n\
angle       17    1   71     50.00     114.70\n\
angle       21    1   23     63.00     110.10\n\
angle       21    1   30     80.00     111.20\n\
angle       21    1   65     63.00     111.10\n\
angle       21    1   66     63.00     111.10\n\
angle       21    1   67     63.00     111.10\n\
angle       21    1   68     63.00     111.10\n\
angle       21    1   69     63.00     111.10\n\
angle       21    1   70     63.00     111.10\n\
angle       21    1   71     63.00     111.10\n\
angle       21    1   72     63.00     111.10\n\
angle       21    1   73     63.00     111.10\n\
angle       21    1   74     63.00     111.10\n\
angle       21    1   75     63.00     111.10\n\
angle       21    1   76     63.00     111.10\n\
angle       21    1   77     63.00     111.10\n\
angle       23    1   27     63.00     110.10\n\
angle       23    1   65     80.00     109.70\n\
angle       23    1   66     80.00     109.70\n\
angle       23    1   67     80.00     109.70\n\
angle       23    1   68     80.00     109.70\n\
angle       23    1   69     80.00     109.70\n\
angle       23    1   70     80.00     109.70\n\
angle       23    1   71     80.00     109.70\n\
angle       23    1   72     80.00     109.70\n\
angle       23    1   73     80.00     109.70\n\
angle       23    1   74     80.00     109.70\n\
angle       23    1   75     80.00     109.70\n\
angle       23    1   76     80.00     109.70\n\
angle       23    1   77     80.00     109.70\n\
angle       27    1   30     80.00     111.20\n\
angle       27    1   65     63.00     111.10\n\
angle       27    1   66     63.00     111.10\n\
angle       27    1   67     63.00     111.10\n\
angle       27    1   68     63.00     111.10\n\
angle       27    1   69     63.00     111.10\n\
angle       27    1   70     63.00     111.10\n\
angle       27    1   71     63.00     111.10\n\
angle       27    1   72     63.00     111.10\n\
angle       27    1   73     63.00     111.10\n\
angle       27    1   74     63.00     111.10\n\
angle       27    1   75     63.00     111.10\n\
angle       27    1   76     63.00     111.10\n\
angle       27    1   77     63.00     111.10\n\
angle       30    1   65     80.00     111.20\n\
angle       30    1   66     80.00     111.20\n\
angle       30    1   67     80.00     111.20\n\
angle       30    1   68     80.00     111.20\n\
angle       30    1   69     80.00     111.20\n\
angle       30    1   70     80.00     111.20\n\
angle       30    1   71     80.00     111.20\n\
angle       30    1   72     80.00     111.20\n\
angle       30    1   73     80.00     111.20\n\
angle       30    1   74     80.00     111.20\n\
angle       30    1   75     80.00     111.20\n\
angle       30    1   76     80.00     111.20\n\
angle       30    1   77     80.00     111.20\n\
angle        1    5    5     70.00     120.00\n\
angle        5    5    5     63.00     120.00\n\
angle        5    5    6     35.00     120.00\n\
angle        5    5   10     70.00     120.00\n\
angle        5    5   35     63.00     120.00\n\
angle        5    5   36     85.00     120.00\n\
angle        6    5   34     35.00     120.00\n\
angle        6    5   35     35.00     120.00\n\
angle        6    5   36     35.00     120.00\n\
angle        6    5   37     35.00     121.60\n\
angle       34    5   37     70.00     108.70\n\
angle        1    7    8     55.00     108.50\n\
angle        8    7   65     55.00     108.50\n\
angle        5   10    8     35.00     113.00\n\
angle        1   15   18     44.00      96.00\n\
angle       18   15   66     44.00      96.00\n\
angle        1   17    1     62.00      98.90\n\
angle        1   17   17     68.00     103.70\n\
angle        1   21   22     80.00     120.40\n\
angle        1   21   23     70.00     116.60\n\
angle       22   21   23     80.00     122.90\n\
angle       22   21   29     35.00     120.00\n\
angle       22   21   67     80.00     120.40\n\
angle       23   21   29     40.00     115.00\n\
angle       23   21   67     70.00     116.60\n\
angle        1   23    1     50.00     118.00\n\
angle        1   23   21     50.00     121.90\n\
angle        1   23   24     38.00     118.40\n\
angle        1   23   42     50.00     118.00\n\
angle       21   23   24     35.00     119.80\n\
angle       21   23   42     50.00     121.90\n\
angle       24   23   24     35.00     120.00\n\
angle        1   27   28     70.00     117.00\n\
angle       28   27   28     80.00     126.00\n\
angle       28   27   73     70.00     117.00\n\
angle        1   30    1     50.00     113.00\n\
angle        1   30   31     35.00     109.50\n\
angle       31   30   31     35.00     109.50\n\
angle        1   32   31     35.00     118.40\n\
angle        1   32   33     50.00     123.20\n\
angle       31   32   31     35.00     120.00\n\
angle       31   32   33     35.00     120.00\n\
angle       32   33   32     70.00     120.00\n\
angle        1   34    5     70.00     125.00\n\
angle        1   34   35     70.00     128.60\n\
angle        5   34   35     85.00     106.40\n\
angle        5   34   72     70.00     125.00\n\
angle       35   34   72     70.00     128.60\n\
angle        5   35   34     85.00     134.90\n\
angle        5   35   36     85.00     116.20\n\
angle       34   35   36     85.00     108.80\n\
angle        5   36   35     85.00     122.70\n\
angle        5   36   37     70.00     132.80\n\
angle       35   36   37     70.00     104.40\n\
angle        5   37   24     35.00     120.00\n\
angle        5   37   36     70.00     111.60\n\
angle       24   37   36     35.00     123.10\n\
angle       24   37   38     35.00     126.35\n\
angle       24   37   40     35.00     126.35\n\
angle       38   37   40     70.00     107.30\n\
angle        6   38   37     35.00     120.00\n\
angle        6   38   41     35.00     120.00\n\
angle       37   38   37     70.00     110.75\n\
angle       37   38   41     70.00     111.60\n\
angle        1   39   40     70.00     130.70\n\
angle        1   39   41     70.00     124.50\n\
angle        6   39   40     35.00     128.20\n\
angle        6   39   41     35.00     120.00\n\
angle       40   39   41     70.00     108.30\n\
angle        1   40   37     70.00     121.60\n\
angle        1   40   39     70.00     130.70\n\
angle        1   40   40     70.00     130.70\n\
angle        6   40   37     35.00     121.60\n\
angle        6   40   39     35.00     130.70\n\
angle        6   40   40     35.00     130.70\n\
angle       37   40   39     70.00     108.70\n\
angle       37   40   40     70.00     106.30\n\
angle       37   40   76     70.00     121.60\n\
angle       40   40   76     70.00     130.70\n\
angle       38   41   39     70.00     105.30\n\
angle        1   42    2     37.50     110.70\n\
angle        1   42   23     80.00     109.70\n\
angle        2   42    2     33.00     107.80\n\
angle        2   42   23     35.00     109.50\n\
angle       44   43   44     34.05     104.52\n\
angle       44   45   44     37.95     109.47\n\
angle        1   65    2     37.50     110.70\n\
angle        1   65    7     50.00     109.50\n\
angle        2   65    2     33.00     107.80\n\
angle        2   65    7     35.00     109.50\n\
angle        1   66    2     37.50     110.70\n\
angle        1   66   15     50.00     108.60\n\
angle        2   66    2     33.00     107.80\n\
angle        2   66   15     35.00     109.50\n\
angle        1   67    2     37.50     110.70\n\
angle        1   67   21     63.00     111.10\n\
angle        2   67    2     33.00     107.80\n\
angle        2   67   21     35.00     109.50\n\
angle        1   68    1     58.35     112.70\n\
angle        1   68    2     37.50     110.70\n\
angle        2   68    2     33.00     107.80\n\
angle        1   69    1     58.35     112.70\n\
angle        1   69    2     37.50     110.70\n\
angle        2   69    2     33.00     107.80\n\
angle        1   70    1     58.35     112.70\n\
angle        1   70    2     37.50     110.70\n\
angle        1   71    1     58.35     112.70\n\
angle        1   71    2     37.50     110.70\n\
angle        2   71    2     33.00     107.80\n\
angle        1   72    2     37.50     110.70\n\
angle        1   72   34     63.00     115.60\n\
angle        2   72    2     33.00     107.80\n\
angle        2   72   34     35.00     109.50\n\
angle        1   73    2     37.50     110.70\n\
angle        1   73   27     63.00     111.10\n\
angle        2   73    2     33.00     107.80\n\
angle        2   73   27     35.00     109.50\n\
angle        1   74    1     58.35     112.70\n\
angle        1   74    2     37.50     110.70\n\
angle        2   74    2     33.00     107.80\n\
angle        1   75    1     58.35     112.70\n\
angle        1   75    2     37.50     110.70\n\
angle        2   75    2     33.00     107.80\n\
angle        1   76    2     37.50     110.70\n\
angle        1   76   40     63.00     114.00\n\
angle        2   76    2     33.00     107.80\n\
angle        2   76   40     35.00     109.50\n\
angle        1   77    1     58.35     112.70\n\
angle        1   77    2     37.50     110.70\n\
angle        2   77    2     33.00     107.80\n\
\n\
\n\
      ###############################\n\
      ##                           ##\n\
      ##  Urey-Bradley Parameters  ##\n\
      ##                           ##\n\
      ###############################\n\
\n\
\n\
ureybrad    44   43   44     38.25     1.5139\n\
ureybrad    44   45   44      39.9     1.6330\n\
\n\
\n\
      ########################################\n\
      ##                                    ##\n\
      ##  Improper Torsional Parameters     ##\n\
      ##   I belive this version has 3      ##\n\
      ##   errors:                          ##\n\
      ##   37 40 40  1 -> 37 40 40 76       ##\n\
      ##   23 29 21 22 -> 23 67 21 22       ##\n\
      ##   35  5 34  1 -> 35  5 34 72       ##\n\
      ##   and one missing parameter:       ##\n\
      ##   28 73 27 28                      ##\n\
      ##   enoee@binf.ku.dk                 ##\n\
      ########################################\n\
\n\
\n\
imptors      5    5    5    1                     2.200 180.0 2\n\
imptors      5    5    5    6                     2.200 180.0 2\n\
imptors      5    5    5   10                     2.200 180.0 2\n\
imptors      5   36    5    6                     2.200 180.0 2\n\
imptors     35    5    5    6                     2.200 180.0 2\n\
imptors     37   34    5    6                     2.200 180.0 2\n\
imptors     23    1   21   22                    21.000 180.0 2\n\
imptors     23   29   21   22                    21.000 180.0 2\n\
imptors      1   21   23   24                     2.000 180.0 2\n\
imptors      1   21   23   42                     2.000 180.0 2\n\
imptors     24   21   23   24                     2.000 180.0 2\n\
imptors     28    1   27   28                    21.000 180.0 2\n\
imptors      1   33   32   31                     2.000 180.0 2\n\
imptors     31   33   32   31                     2.000 180.0 2\n\
imptors     32   32   33   32                    21.000 180.0 2\n\
imptors     35    5   34    1                     2.200 180.0 2\n\
imptors     34    5   35   36                     2.200 180.0 2\n\
imptors     37   35   36    5                     2.200 180.0 2\n\
imptors     36    5   37   24                     2.000 180.0 2\n\
imptors     40   38   37   24                     2.000 180.0 2\n\
imptors     37   37   38    6                     2.200 180.0 2\n\
imptors     37   41   38    6                     2.200 180.0 2\n\
imptors     40   41   39    1                     2.200 180.0 2\n\
imptors     40   41   39    6                     2.200 180.0 2\n\
imptors     37   39   40    1                     2.200 180.0 2\n\
imptors     37   39   40    6                     2.200 180.0 2\n\
imptors     37   40   40    1                     2.200 180.0 2\n\
imptors     37   40   40    6                     2.200 180.0 2\n\
\n\
\n\
      ############################\n\
      ##                        ##\n\
      ##  Torsional Parameters  ##\n\
      ##                        ##\n\
      ############################\n\
\n\
\n\
torsion      1    1    1    1      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1    1    2                                      0.366 0.0 3\n\
torsion      1    1    1   17      2.619 0.0 1   -0.620 180.0 2    0.258 0.0 3\n\
torsion      1    1    1   21     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion      1    1    1   23      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion      1    1    1   27     -3.185 0.0 1   -0.825 180.0 2    0.493 0.0 3\n\
torsion      1    1    1   30      2.732 0.0 1   -0.229 180.0 2    0.485 0.0 3\n\
torsion      1    1    1   32      1.964 0.0 1                     0.659 0.0 3\n\
torsion      1    1    1   42      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1    1   75      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      2    1    1    2                                      0.318 0.0 3\n\
torsion      2    1    1    5                                      0.462 0.0 3\n\
torsion      2    1    1    7                                      0.468 0.0 3\n\
torsion      2    1    1   15                                      0.452 0.0 3\n\
torsion      2    1    1   17                                      0.452 0.0 3\n\
torsion      2    1    1   21                                     -0.076 0.0 3\n\
torsion      2    1    1   23                                      0.464 0.0 3\n\
torsion      2    1    1   27                                     -0.225 0.0 3\n\
torsion      2    1    1   30                                      0.384 0.0 3\n\
torsion      2    1    1   32                                      0.464 0.0 3\n\
torsion      2    1    1   34                                      0.462 0.0 3\n\
torsion      2    1    1   39                                      0.462 0.0 3\n\
torsion      2    1    1   40                                      0.462 0.0 3\n\
torsion      2    1    1   42                                      0.366 0.0 3\n\
torsion      2    1    1   69                                      0.366 0.0 3\n\
torsion      2    1    1   70                                      0.366 0.0 3\n\
torsion      2    1    1   75                                      0.366 0.0 3\n\
torsion      2    1    1   77                                      0.366 0.0 3\n\
torsion      5    1    1   21     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion      5    1    1   23      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion      5    1    1   27     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion      5    1    1   30      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion      7    1    1   21     -6.180 0.0 1\n\
torsion      7    1    1   23      6.280 0.0 1   -1.467 180.0 2    2.030 0.0 3\n\
torsion      7    1    1   27     -6.180 0.0 1\n\
torsion      7    1    1   30      6.280 0.0 1   -1.467 180.0 2    2.030 0.0 3\n\
torsion     15    1    1   21     -4.214 0.0 1   -2.114 180.0 2    0.969 0.0 3\n\
torsion     15    1    1   23      0.583 0.0 1   -1.163 180.0 2    0.141 0.0 3\n\
torsion     15    1    1   27     -4.214 0.0 1   -2.114 180.0 2    0.969 0.0 3\n\
torsion     15    1    1   30      0.583 0.0 1   -1.163 180.0 2    0.141 0.0 3\n\
torsion     17    1    1   21     -4.214 0.0 1   -2.114 180.0 2    0.969 0.0 3\n\
torsion     17    1    1   23      0.583 0.0 1   -1.163 180.0 2    0.141 0.0 3\n\
torsion     17    1    1   27     -4.214 0.0 1   -2.114 180.0 2    0.969 0.0 3\n\
torsion     17    1    1   30      0.583 0.0 1   -1.163 180.0 2    0.141 0.0 3\n\
torsion     21    1    1   21     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion     21    1    1   23      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     21    1    1   27     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion     21    1    1   30      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     21    1    1   34     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion     21    1    1   39     -1.607 0.0 1    0.046 180.0 2\n\
torsion     21    1    1   40     -1.607 0.0 1    0.046 180.0 2\n\
torsion     23    1    1   27      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     23    1    1   34      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     23    1    1   39     -0.713 0.0 1    0.502 180.0 2    0.289 0.0 3\n\
torsion     23    1    1   40     -0.713 0.0 1    0.502 180.0 2    0.289 0.0 3\n\
torsion     27    1    1   27      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     27    1    1   30      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     27    1    1   34     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion     27    1    1   39     -1.607 0.0 1    0.046 180.0 2\n\
torsion     30    1    1   34      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion     32    1    1   77      1.964 0.0 1                     0.659 0.0 3\n\
torsion      1    1    5    5\n\
torsion      2    1    5    5\n\
torsion      1    1    7    8     -0.356 0.0 1   -0.174 180.0 2    0.492 0.0 3\n\
torsion      2    1    7    8                                      0.450 0.0 3\n\
torsion      1    1   15   18     -0.759 0.0 1   -0.282 180.0 2    0.603 0.0 3\n\
torsion      2    1   15   18                                      0.451 0.0 3\n\
torsion      1    1   17    1      0.925 0.0 1   -0.576 180.0 2    0.677 0.0 3\n\
torsion      1    1   17   17      1.941 0.0 1   -0.836 180.0 2    0.935 0.0 3\n\
torsion      2    1   17    1                                      0.647 0.0 3\n\
torsion      2    1   17   17                                      0.558 0.0 3\n\
torsion     71    1   17    1      0.925 0.0 1   -0.576 180.0 2    0.677 0.0 3\n\
torsion      1    1   21   22\n\
torsion      1    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion      2    1   21   22\n\
torsion      2    1   21   23\n\
torsion     23    1   21   22\n\
torsion     23    1   21   23      0.743 0.0 1    2.508 180.0 2   -0.805 0.0 3\n\
torsion     30    1   21   22\n\
torsion     30    1   21   23      1.816 0.0 1    1.222 180.0 2    1.581 0.0 3\n\
torsion     65    1   21   22\n\
torsion     65    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     66    1   21   22\n\
torsion     66    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     67    1   21   22\n\
torsion     67    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     68    1   21   22\n\
torsion     68    1   21   23     -0.290 0.0 1    2.621 180.0 2   -1.778 0.0 3\n\
torsion     69    1   21   22\n\
torsion     69    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     70    1   21   22\n\
torsion     70    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     71    1   21   22\n\
torsion     71    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     72    1   21   22\n\
torsion     72    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     73    1   21   22\n\
torsion     73    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     74    1   21   22\n\
torsion     74    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     75    1   21   22\n\
torsion     75    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     76    1   21   22\n\
torsion     76    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion     77    1   21   22\n\
torsion     77    1   21   23      1.865 0.0 1    0.089 180.0 2    0.351 0.0 3\n\
torsion      1    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion      1    1   23   24\n\
torsion      1    1   23   42      4.753 0.0 1   -0.734 180.0 2\n\
torsion      2    1   23    1\n\
torsion      2    1   23   21\n\
torsion      2    1   23   24\n\
torsion      2    1   23   42\n\
torsion     21    1   23   21     -0.596 0.0 1    0.279 180.0 2   -4.913 0.0 3\n\
torsion     21    1   23   24\n\
torsion     21    1   23   42     -1.737 0.0 1    1.251 180.0 2   -3.501 0.0 3\n\
torsion     27    1   23   21     -2.365 0.0 1    0.912 180.0 2   -0.850 0.0 3\n\
torsion     27    1   23   24\n\
torsion     27    1   23   42     -1.737 0.0 1    1.251 180.0 2   -3.501 0.0 3\n\
torsion     65    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     65    1   23   24\n\
torsion     66    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     66    1   23   24\n\
torsion     67    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     67    1   23   24\n\
torsion     68    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     68    1   23   24\n\
torsion     69    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     69    1   23   24\n\
torsion     70    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     70    1   23   24\n\
torsion     71    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     71    1   23   24\n\
torsion     72    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     72    1   23   24\n\
torsion     73    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     73    1   23   24\n\
torsion     74    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     74    1   23   24\n\
torsion     75    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     75    1   23   24\n\
torsion     76    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     76    1   23   24\n\
torsion     77    1   23   21      0.519 0.0 1    0.877 180.0 2    5.233 0.0 3\n\
torsion     77    1   23   24\n\
torsion      1    1   27   28                     0.820 180.0 2\n\
torsion      2    1   27   28\n\
torsion     23    1   27   28\n\
torsion     30    1   27   28\n\
torsion     65    1   27   28                     0.820 180.0 2\n\
torsion     66    1   27   28                     0.820 180.0 2\n\
torsion     67    1   27   28                     0.820 180.0 2\n\
torsion     68    1   27   28                     0.820 180.0 2\n\
torsion     69    1   27   28                     0.820 180.0 2\n\
torsion     70    1   27   28                     0.820 180.0 2\n\
torsion     71    1   27   28                     0.820 180.0 2\n\
torsion     72    1   27   28                     0.820 180.0 2\n\
torsion     73    1   27   28                     0.820 180.0 2\n\
torsion     74    1   27   28                     0.820 180.0 2\n\
torsion     75    1   27   28                     0.820 180.0 2\n\
torsion     76    1   27   28                     0.820 180.0 2\n\
torsion     77    1   27   28                     0.820 180.0 2\n\
torsion      1    1   30    1      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1   30   31                                      0.347 0.0 3\n\
torsion      2    1   30    1                                      0.366 0.0 3\n\
torsion      2    1   30   31                                      0.261 0.0 3\n\
torsion     21    1   30    1     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion     21    1   30   31                                      0.347 0.0 3\n\
torsion     27    1   30   31                                      0.347 0.0 3\n\
torsion     65    1   30   31                                      0.347 0.0 3\n\
torsion     66    1   30   31                                      0.347 0.0 3\n\
torsion     67    1   30   31                                      0.347 0.0 3\n\
torsion     68    1   30   31                                      0.347 0.0 3\n\
torsion     69    1   30   31                                      0.347 0.0 3\n\
torsion     70    1   30   31                                      0.347 0.0 3\n\
torsion     71    1   30   31                                      0.347 0.0 3\n\
torsion     72    1   30   31                                      0.347 0.0 3\n\
torsion     73    1   30   31                                      0.347 0.0 3\n\
torsion     74    1   30   31                                      0.347 0.0 3\n\
torsion     75    1   30   31                                      0.347 0.0 3\n\
torsion     76    1   30   31                                      0.347 0.0 3\n\
torsion     77    1   30   31                                      0.347 0.0 3\n\
torsion      1    1   32   31\n\
torsion      1    1   32   33      1.829 0.0 1    0.243 180.0 2   -0.498 0.0 3\n\
torsion      2    1   32   31\n\
torsion      2    1   32   33                                      0.177 0.0 3\n\
torsion      1    1   34    5     -0.714 0.0 1\n\
torsion      1    1   34   35\n\
torsion      2    1   34    5                                     -0.480 0.0 3\n\
torsion      2    1   34   35\n\
torsion      1    1   39   40\n\
torsion      1    1   39   41     -0.543 0.0 1    0.014 180.0 2    0.700 0.0 3\n\
torsion      2    1   39   40\n\
torsion      2    1   39   41                                      0.419 0.0 3\n\
torsion      1    1   40   37     -0.543 0.0 1    0.014 180.0 2    0.700 0.0 3\n\
torsion      1    1   40   39\n\
torsion      1    1   40   40\n\
torsion      2    1   40   37                                      0.419 0.0 3\n\
torsion      2    1   40   39\n\
torsion      2    1   40   40\n\
torsion      1    1   42    2                                      0.366 0.0 3\n\
torsion      1    1   42   23      0.845 0.0 1   -0.962 180.0 2    0.713 0.0 3\n\
torsion      2    1   42    2                                      0.318 0.0 3\n\
torsion      2    1   42   23                                      0.464 0.0 3\n\
torsion      2    1   65    2                                      0.318 0.0 3\n\
torsion      2    1   65    7                                      0.468 0.0 3\n\
torsion     21    1   65    2                                     -0.076 0.0 3\n\
torsion     21    1   65    7     -5.654 0.0 1   -0.872 180.0 2\n\
torsion     23    1   65    2                                      0.464 0.0 3\n\
torsion     23    1   65    7      5.429 0.0 1   -0.879 180.0 2    1.058 0.0 3\n\
torsion     27    1   65    2                                     -0.225 0.0 3\n\
torsion     27    1   65    7     -5.654 0.0 1   -0.872 180.0 2\n\
torsion     30    1   65    2                                      0.384 0.0 3\n\
torsion     30    1   65    7      5.429 0.0 1   -0.879 180.0 2    1.058 0.0 3\n\
torsion      2    1   66    2                                      0.318 0.0 3\n\
torsion      2    1   66   15                                      0.452 0.0 3\n\
torsion     21    1   66    2                                     -0.076 0.0 3\n\
torsion     21    1   66   15     -4.344 0.0 1   -1.714 180.0 2\n\
torsion     23    1   66    2                                      0.464 0.0 3\n\
torsion     23    1   66   15      1.428 0.0 1    0.086 180.0 2    0.029 0.0 3\n\
torsion     27    1   66    2                                     -0.225 0.0 3\n\
torsion     27    1   66   15     -4.344 0.0 1   -1.714 180.0 2\n\
torsion     30    1   66    2                                      0.384 0.0 3\n\
torsion     30    1   66   15      1.428 0.0 1    0.086 180.0 2    0.029 0.0 3\n\
torsion      2    1   67    2                                      0.318 0.0 3\n\
torsion      2    1   67   21                                     -0.076 0.0 3\n\
torsion     21    1   67    2                                     -0.076 0.0 3\n\
torsion     21    1   67   21      1.045 0.0 1   -1.603 180.0 2\n\
torsion     21    1   67   27      1.045 0.0 1   -1.603 180.0 2\n\
torsion     23    1   67    2                                      0.464 0.0 3\n\
torsion     23    1   67   21     -3.467 0.0 1   -0.677 180.0 2    1.465 0.0 3\n\
torsion     27    1   67    2                                     -0.225 0.0 3\n\
torsion     27    1   67   21     -1.697 0.0 1   -0.456 180.0 2    0.585 0.0 3\n\
torsion     30    1   67    2                                      0.384 0.0 3\n\
torsion     30    1   67   21     -3.467 0.0 1   -0.677 180.0 2    1.465 0.0 3\n\
torsion      2    1   68    1                                      0.366 0.0 3\n\
torsion      2    1   68    2                                      0.318 0.0 3\n\
torsion     21    1   68    1     -2.692 0.0 1    0.836 180.0 2\n\
torsion     21    1   68    2                                     -0.076 0.0 3\n\
torsion     23    1   68    1      1.167 0.0 1   -0.526 180.0 2    0.805 0.0 3\n\
torsion     23    1   68    2                                      0.464 0.0 3\n\
torsion     27    1   68    1     -2.692 0.0 1    0.836 180.0 2\n\
torsion     27    1   68    2                                     -0.225 0.0 3\n\
torsion     30    1   68    1      1.167 0.0 1   -0.526 180.0 2    0.805 0.0 3\n\
torsion     30    1   68    2                                      0.384 0.0 3\n\
torsion      1    1   69    1      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1   69    2                                      0.366 0.0 3\n\
torsion      2    1   69    1                                      0.366 0.0 3\n\
torsion      2    1   69    2                                      0.318 0.0 3\n\
torsion     21    1   69    1     -0.916 0.0 1    0.411 180.0 2\n\
torsion     21    1   69    2                                     -0.076 0.0 3\n\
torsion     23    1   69    1      0.528 0.0 1    0.323 180.0 2    0.517 0.0 3\n\
torsion     23    1   69    2                                      0.464 0.0 3\n\
torsion     27    1   69    1     -0.916 0.0 1    0.411 180.0 2\n\
torsion     27    1   69    2                                     -0.225 0.0 3\n\
torsion     30    1   69    1      0.528 0.0 1    0.323 180.0 2    0.517 0.0 3\n\
torsion     30    1   69    2                                      0.384 0.0 3\n\
torsion      1    1   70    1      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1   70    2                                      0.366 0.0 3\n\
torsion      2    1   70    1                                      0.366 0.0 3\n\
torsion      2    1   70    2                                      0.318 0.0 3\n\
torsion     21    1   70    1      0.678 0.0 1    1.276 180.0 2\n\
torsion     21    1   70    2                                     -0.076 0.0 3\n\
torsion     23    1   70    1      2.647 0.0 1    0.954 180.0 2    0.522 0.0 3\n\
torsion     23    1   70    2                                      0.464 0.0 3\n\
torsion     27    1   70    1      0.678 0.0 1    1.276 180.0 2\n\
torsion     27    1   70    2                                     -0.225 0.0 3\n\
torsion     30    1   70    1      2.647 0.0 1    0.954 180.0 2    0.522 0.0 3\n\
torsion     30    1   70    2                                      0.384 0.0 3\n\
torsion      2    1   71    1                                      0.366 0.0 3\n\
torsion      2    1   71    2                                      0.318 0.0 3\n\
torsion     17    1   71    1      2.619 0.0 1   -0.620 180.0 2    0.258 0.0 3\n\
torsion     17    1   71    2                                      0.452 0.0 3\n\
torsion     21    1   71    1      0.384 0.0 1   -0.762 180.0 2\n\
torsion     21    1   71    2                                     -0.076 0.0 3\n\
torsion     23    1   71    1      1.493 0.0 1    0.199 180.0 2   -0.098 0.0 3\n\
torsion     23    1   71    2                                      0.464 0.0 3\n\
torsion     27    1   71    1      0.384 0.0 1   -0.762 180.0 2\n\
torsion     27    1   71    2                                     -0.225 0.0 3\n\
torsion     30    1   71    1      1.493 0.0 1    0.199 180.0 2   -0.098 0.0 3\n\
torsion     30    1   71    2                                      0.384 0.0 3\n\
torsion      2    1   72    2                                      0.318 0.0 3\n\
torsion      2    1   72   34                                      0.462 0.0 3\n\
torsion     21    1   72    2                                     -0.076 0.0 3\n\
torsion     21    1   72   34     -1.058 0.0 1   -0.625 180.0 2\n\
torsion     23    1   72    2                                      0.464 0.0 3\n\
torsion     23    1   72   34     -1.294 0.0 1    0.562 180.0 2    0.094 0.0 3\n\
torsion     27    1   72    2                                     -0.225 0.0 3\n\
torsion     27    1   72   34     -1.058 0.0 1   -0.625 180.0 2\n\
torsion     30    1   72    2                                      0.384 0.0 3\n\
torsion     30    1   72   34     -1.294 0.0 1    0.562 180.0 2    0.094 0.0 3\n\
torsion      2    1   73    2                                      0.318 0.0 3\n\
torsion      2    1   73   27                                     -0.225 0.0 3\n\
torsion     21    1   73    2                                     -0.076 0.0 3\n\
torsion     21    1   73   27     -2.059 0.0 1    3.103 180.0 2\n\
torsion     23    1   73    2                                      0.464 0.0 3\n\
torsion     23    1   73   27     -4.500 0.0 1    0.958 180.0 2   -0.293 0.0 3\n\
torsion     27    1   73    2                                     -0.225 0.0 3\n\
torsion     27    1   73   27     -2.059 0.0 1    3.103 180.0 2\n\
torsion     30    1   73    2                                      0.384 0.0 3\n\
torsion     30    1   73   27     -4.500 0.0 1    0.958 180.0 2   -0.293 0.0 3\n\
torsion      2    1   74    1                                      0.366 0.0 3\n\
torsion      2    1   74    2                                      0.318 0.0 3\n\
torsion     21    1   74    1     -1.618 0.0 1   -0.571 180.0 2\n\
torsion     21    1   74    2                                     -0.076 0.0 3\n\
torsion     23    1   74    1      4.952 0.0 1   -0.257 180.0 2   -0.235 0.0 3\n\
torsion     23    1   74    2                                      0.464 0.0 3\n\
torsion     27    1   74    1     -1.618 0.0 1   -0.571 180.0 2\n\
torsion     27    1   74    2                                     -0.225 0.0 3\n\
torsion     30    1   74    1      4.952 0.0 1   -0.257 180.0 2   -0.235 0.0 3\n\
torsion     30    1   74    2                                      0.384 0.0 3\n\
torsion      1    1   75    1      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1   75    2                                      0.366 0.0 3\n\
torsion      2    1   75    1                                      0.366 0.0 3\n\
torsion      2    1   75    2                                      0.318 0.0 3\n\
torsion     21    1   75    1     -2.235 0.0 1   -0.536 180.0 2\n\
torsion     21    1   75    2                                     -0.076 0.0 3\n\
torsion     23    1   75    1      0.639 0.0 1   -0.214 180.0 2    0.399 0.0 3\n\
torsion     23    1   75    2                                      0.464 0.0 3\n\
torsion     27    1   75    1     -2.235 0.0 1   -0.536 180.0 2\n\
torsion     27    1   75    2                                     -0.225 0.0 3\n\
torsion     30    1   75    1      0.639 0.0 1   -0.214 180.0 2    0.399 0.0 3\n\
torsion     30    1   75    2                                      0.384 0.0 3\n\
torsion      2    1   76    2                                      0.318 0.0 3\n\
torsion      2    1   76   40                                      0.462 0.0 3\n\
torsion     21    1   76    2                                     -0.076 0.0 3\n\
torsion     21    1   76   40      1.679 0.0 1    0.082 180.0 2\n\
torsion     23    1   76    2                                      0.464 0.0 3\n\
torsion     23    1   76   40      0.347 0.0 1   -0.350 180.0 2    1.468 0.0 3\n\
torsion     27    1   76    2                                     -0.225 0.0 3\n\
torsion     27    1   76   40      1.679 0.0 1    0.082 180.0 2\n\
torsion     30    1   76    2                                      0.384 0.0 3\n\
torsion     30    1   76   40      0.347 0.0 1   -0.350 180.0 2    1.468 0.0 3\n\
torsion      1    1   77    1      1.740 0.0 1   -0.157 180.0 2    0.279 0.0 3\n\
torsion      1    1   77    2                                      0.366 0.0 3\n\
torsion      2    1   77    1                                      0.366 0.0 3\n\
torsion      2    1   77    2                                      0.318 0.0 3\n\
torsion     21    1   77    1     -0.678 0.0 1    1.653 180.0 2\n\
torsion     21    1   77    2                                     -0.076 0.0 3\n\
torsion     23    1   77    1      1.821 0.0 1    0.951 180.0 2    1.168 0.0 3\n\
torsion     23    1   77    2                                      0.464 0.0 3\n\
torsion     27    1   77    1     -0.678 0.0 1    1.653 180.0 2\n\
torsion     27    1   77    2                                     -0.225 0.0 3\n\
torsion     30    1   77    1      1.821 0.0 1    0.951 180.0 2    1.168 0.0 3\n\
torsion     30    1   77    2                                      0.384 0.0 3\n\
torsion      1    5    5    5                     7.250 180.0 2\n\
torsion      1    5    5    6                     7.250 180.0 2\n\
torsion      5    5    5    5                     7.250 180.0 2\n\
torsion      5    5    5    6                     7.250 180.0 2\n\
torsion      5    5    5   10                     7.250 180.0 2\n\
torsion      5    5    5   35                     7.250 180.0 2\n\
torsion      5    5    5   36                     7.250 180.0 2\n\
torsion      6    5    5    6                     7.250 180.0 2\n\
torsion      6    5    5   10                     7.250 180.0 2\n\
torsion      6    5    5   35                     7.250 180.0 2\n\
torsion      6    5    5   36                     7.250 180.0 2\n\
torsion      5    5   10    8                     1.682 180.0 2\n\
torsion      6    5   34    1                    13.050 180.0 2\n\
torsion      6    5   34   35                    13.050 180.0 2\n\
torsion      6    5   34   72                    13.050 180.0 2\n\
torsion     37    5   34    1                    13.050 180.0 2\n\
torsion     37    5   34   35                    13.050 180.0 2\n\
torsion     37    5   34   72                    13.050 180.0 2\n\
torsion      5    5   35   34                     7.000 180.0 2\n\
torsion      5    5   35   36                     7.000 180.0 2\n\
torsion      6    5   35   34                     7.000 180.0 2\n\
torsion      6    5   35   36                     7.000 180.0 2\n\
torsion      5    5   36   35                     7.250 180.0 2\n\
torsion      5    5   36   37                     7.250 180.0 2\n\
torsion      6    5   36   35                     7.250 180.0 2\n\
torsion      6    5   36   37                     7.250 180.0 2\n\
torsion      6    5   37   24                     3.000 180.0 2\n\
torsion      6    5   37   36                     3.000 180.0 2\n\
torsion     34    5   37   24                     3.000 180.0 2\n\
torsion     34    5   37   36                     3.000 180.0 2\n\
torsion      8    7   65    1     -0.991 0.0 1   -0.869 180.0 2    0.739 0.0 3\n\
torsion      8    7   65    2                                      0.450 0.0 3\n\
torsion     18   15   66    1     -0.759 0.0 1   -0.282 180.0 2    0.603 0.0 3\n\
torsion     18   15   66    2                                      0.451 0.0 3\n\
torsion      1   17   17    1                    -7.414 180.0 2    1.705 0.0 3\n\
torsion      1   21   23    1      2.300 0.0 1    6.089 180.0 2\n\
torsion      1   21   23   24                     4.900 180.0 2\n\
torsion      1   21   23   42      2.300 0.0 1    6.089 180.0 2\n\
torsion     22   21   23    1                     6.089 180.0 2\n\
torsion     22   21   23   24                     4.900 180.0 2\n\
torsion     22   21   23   42                     6.089 180.0 2\n\
torsion     29   21   23    1                     6.089 180.0 2\n\
torsion     29   21   23   24                     4.900 180.0 2\n\
torsion     29   21   23   42                     6.089 180.0 2\n\
torsion     67   21   23   24                     4.900 180.0 2\n\
torsion     22   21   67    1\n\
torsion     22   21   67    2\n\
torsion     23   21   67    1      0.546 0.0 1   -2.127 180.0 2   -0.832 0.0 3\n\
torsion     23   21   67    2\n\
torsion      1   23   42    1      2.859 0.0 1    2.058 180.0 2  -11.266 0.0 3\n\
torsion      1   23   42    2\n\
torsion     21   23   42    1                     0.462 180.0 2\n\
torsion     21   23   42    2\n\
torsion     28   27   73    1                     0.820 180.0 2\n\
torsion     28   27   73    2\n\
torsion      1   32   33   32                     7.936 180.0 2\n\
torsion     31   32   33   32                     3.900 180.0 2\n\
torsion      1   34   35    5                     3.350 180.0 2\n\
torsion      1   34   35   36                     3.350 180.0 2\n\
torsion      5   34   35    5                     3.350 180.0 2\n\
torsion      5   34   35   36                     3.350 180.0 2\n\
torsion     72   34   35    5                     3.350 180.0 2\n\
torsion     72   34   35   36                     3.350 180.0 2\n\
torsion      5   34   72    1     -0.714 0.0 1\n\
torsion      5   34   72    2                                     -0.480 0.0 3\n\
torsion     35   34   72    1\n\
torsion     35   34   72    2\n\
torsion      5   35   36    5                     6.000 180.0 2\n\
torsion      5   35   36   37                     6.000 180.0 2\n\
torsion     34   35   36    5                     6.000 180.0 2\n\
torsion     34   35   36   37                     6.000 180.0 2\n\
torsion      5   36   37    5                     3.050 180.0 2\n\
torsion      5   36   37   24                     3.050 180.0 2\n\
torsion     35   36   37    5                     3.050 180.0 2\n\
torsion     35   36   37   24                     3.050 180.0 2\n\
torsion     24   37   38    6                     4.650 180.0 2\n\
torsion     24   37   38   37                     4.650 180.0 2\n\
torsion     24   37   38   41                     4.650 180.0 2\n\
torsion     40   37   38    6                     4.650 180.0 2\n\
torsion     40   37   38   37                     4.650 180.0 2\n\
torsion     40   37   38   41                     4.650 180.0 2\n\
torsion     24   37   40    1                     2.800 180.0 2\n\
torsion     24   37   40    6                     3.200 180.0 2\n\
torsion     24   37   40   39                     2.800 180.0 2\n\
torsion     24   37   40   40                     2.800 180.0 2\n\
torsion     24   37   40   76                     2.800 180.0 2\n\
torsion     38   37   40    1                     2.800 180.0 2\n\
torsion     38   37   40    6                     3.200 180.0 2\n\
torsion     38   37   40   39                     2.800 180.0 2\n\
torsion     38   37   40   40                     2.800 180.0 2\n\
torsion     38   37   40   76                     2.800 180.0 2\n\
torsion      6   38   41   39                    10.000 180.0 2\n\
torsion     37   38   41   39                    10.000 180.0 2\n\
torsion      1   39   40    6                    10.750 180.0 2\n\
torsion      1   39   40   37                    10.750 180.0 2\n\
torsion      6   39   40    1                    10.750 180.0 2\n\
torsion      6   39   40   37                    10.750 180.0 2\n\
torsion     41   39   40    1                    10.750 180.0 2\n\
torsion     41   39   40    6                    10.750 180.0 2\n\
torsion     41   39   40   37                    10.750 180.0 2\n\
torsion      1   39   41   38                     4.800 180.0 2\n\
torsion      6   39   41   38                     4.800 180.0 2\n\
torsion     40   39   41   38                     4.800 180.0 2\n\
torsion      1   40   40    6                    10.750 180.0 2\n\
torsion      1   40   40   37                    10.750 180.0 2\n\
torsion      6   40   40   37                    10.750 180.0 2\n\
torsion      6   40   40   76                    10.750 180.0 2\n\
torsion     37   40   40   37                    10.750 180.0 2\n\
torsion     37   40   40   76                    10.750 180.0 2\n\
torsion     37   40   76    1     -0.543 0.0 1    0.014 180.0 2    0.700 0.0 3\n\
torsion     37   40   76    2                                      0.419 0.0 3\n\
torsion     40   40   76    1\n\
torsion     40   40   76    2\n\
\n\
\n\
      ########################################\n\
      ##                                    ##\n\
      ##  Atomic Partial Charge Parameters  ##\n\
      ##                                    ##\n\
      ########################################\n\
\n\
\n\
charge       1       -0.1800\n\
charge       2       -0.1200\n\
charge       3       -0.0600\n\
charge       4       -0.2400\n\
charge       5        0.0000\n\
charge       6        0.0600\n\
charge       7        0.0000\n\
charge       8       -0.1150\n\
charge       9       -0.2300\n\
charge      10        0.1150\n\
charge      11       -0.1150\n\
charge      12        0.1150\n\
charge      13        0.0000\n\
charge      14       -0.0650\n\
charge      15       -0.0050\n\
charge      16       -0.6830\n\
charge      17        0.4180\n\
charge      18        0.0400\n\
charge      19        0.1450\n\
charge      20        0.2050\n\
charge      21        0.2650\n\
charge      22        0.1263\n\
charge      23        0.5323\n\
charge      24       -0.6351\n\
charge      25        0.4286\n\
charge      26       -0.2057\n\
charge      27        0.0825\n\
charge      28        0.1500\n\
charge      29       -0.5850\n\
charge      30        0.4350\n\
charge      31        0.7000\n\
charge      32        0.4350\n\
charge      33       -0.4000\n\
charge      34        0.1100\n\
charge      35        0.1400\n\
charge      36        0.1700\n\
charge      37        0.2000\n\
charge      38        0.0300\n\
charge      39       -0.4000\n\
charge      40       -0.7000\n\
charge      41        0.4350\n\
charge      42        0.2000\n\
charge      43        0.1000\n\
charge      44        0.2650\n\
charge      45        0.1000\n\
charge      46        0.3000\n\
charge      47        0.1000\n\
charge      48        0.3650\n\
charge      49        0.1000\n\
charge      50        0.4000\n\
charge      51        0.4650\n\
charge      52       -0.3350\n\
charge      53       -0.4700\n\
charge      54       -0.3350\n\
charge      55       -0.2175\n\
charge      56        0.1550\n\
charge      57        0.2350\n\
charge      58        0.0000\n\
charge      59        0.0600\n\
charge      60        0.1200\n\
charge      61        0.1800\n\
charge      62       -0.0125\n\
charge      63        0.0475\n\
charge      64        0.1075\n\
charge      65        0.1675\n\
charge      66        0.0375\n\
charge      67        0.0975\n\
charge      68        0.1575\n\
charge      69        0.2175\n\
charge      70       -0.9000\n\
charge      71        0.3500\n\
charge      72        0.0200\n\
charge      73        0.0800\n\
charge      74        0.1400\n\
charge      75        0.2000\n\
charge      76       -0.9000\n\
charge      77        0.4000\n\
charge      78        0.1000\n\
charge      79       -0.6600\n\
charge      80        0.0400\n\
charge      81        0.1000\n\
charge      82        0.5000\n\
charge      83       -0.5000\n\
charge      84       -0.7600\n\
charge      85       -0.5000\n\
charge      86       -0.1400\n\
charge      87        0.3800\n\
charge      88        0.3000\n\
charge      89        0.0200\n\
charge      90       -0.1100\n\
charge      91        0.0800\n\
charge      92       -0.0500\n\
charge      93        0.0100\n\
charge      94        0.1420\n\
charge      95       -0.3900\n\
charge      96       -0.5420\n\
charge      97        0.3330\n\
charge      98        0.4200\n\
charge      99       -0.4200\n\
charge     100       -0.4900\n\
charge     101        0.3700\n\
charge     102        0.0600\n\
charge     103       -0.1200\n\
charge     104       -0.0600\n\
charge     105        0.0000\n\
charge     106        0.0600\n\
charge     107        0.5200\n\
charge     108       -0.5300\n\
charge     109       -0.4400\n\
charge     110        0.4500\n\
charge     111        0.7000\n\
charge     112       -0.8000\n\
charge     113       -0.2800\n\
charge     114       -0.2200\n\
charge     115       -0.1600\n\
charge     116       -0.1000\n\
charge     117        0.4500\n\
charge     118       -0.4500\n\
charge     119        0.0000\n\
charge     120        0.4700\n\
charge     121       -0.4700\n\
charge     122        0.0600\n\
charge     123       -0.4000\n\
charge     124       -0.3000\n\
charge     125       -0.2400\n\
charge     126       -0.2200\n\
charge     127        0.0000\n\
charge     128        0.3500\n\
charge     129        0.3300\n\
charge     130        0.3100\n\
charge     131        0.2900\n\
charge     132        0.1300\n\
charge     133        0.1900\n\
charge     134        0.2500\n\
charge     135        0.3100\n\
charge     136       -0.8000\n\
charge     137       -0.7000\n\
charge     138        0.4600\n\
charge     139        0.4400\n\
charge     140        0.6400\n\
charge     141        0.2000\n\
charge     142       -0.1100\n\
charge     143        0.1900\n\
charge     144       -0.0500\n\
charge     145        0.0750\n\
charge     146       -0.0550\n\
charge     147        0.1300\n\
charge     148       -0.5700\n\
charge     149        0.4200\n\
charge     150       -0.0050\n\
charge     151        0.2950\n\
charge     152       -0.0150\n\
charge     153        0.0150\n\
charge     154        0.3850\n\
charge     155        0.2150\n\
charge     156       -0.4900\n\
charge     157       -0.5400\n\
charge     158        0.4600\n\
charge     159       -0.0500\n\
charge     160       -0.0200\n\
charge     161        0.0400\n\
charge     162        0.1000\n\
charge     163       -0.0900\n\
charge     164       -0.8340\n\
charge     165        0.4170\n\
charge     166       -0.8200\n\
charge     167        0.4100\n\
charge     168       -1.0260\n\
charge     169        0.3420\n\
charge     170       -1.3000\n\
charge     171        0.3000\n\
charge     172        1.0000\n\
charge     173        1.0000\n\
charge     174        1.0000\n\
charge     175        1.0000\n\
charge     176        1.0000\n\
charge     177        2.0000\n\
charge     178        2.0000\n\
charge     179        2.0000\n\
charge     180        2.0000\n\
charge     181       -1.0000\n\
charge     182       -1.0000\n\
charge     183       -1.0000\n\
charge     184        0.0000\n\
charge     185        0.0000\n\
charge     186        0.0000\n\
charge     187        0.0000\n\
charge     188        0.0000\n\
charge     189        0.1450\n\
charge     190        0.0600\n\
charge     191       -0.1200\n\
charge     192       -0.1200\n\
charge     193       -0.1200\n\
charge     194       -0.0600\n\
charge     195       -0.0600\n\
charge     196       -0.1200\n\
charge     197       -0.1200\n\
charge     198       -0.2200\n\
charge     199       -0.1200\n\
charge     200       -0.1200\n\
charge     201       -0.0050\n\
charge     202       -0.1200\n\
\n\
\n\
      ########################################\n\
      ##                                    ##\n\
      ##  Biopolymer Atom Type Conversions  ##\n\
      ##                                    ##\n\
      ########################################\n\
\n\
\n\
biotype      1    N       \"Glycine\"                  85\n\
biotype      2    CA      \"Glycine\"                  73\n\
biotype      3    C       \"Glycine\"                  82\n\
biotype      4    HN      \"Glycine\"                  88\n\
biotype      5    O       \"Glycine\"                  83\n\
biotype      6    HA      \"Glycine\"                   6\n\
biotype      7    N       \"Alanine\"                  85\n\
biotype      8    CA      \"Alanine\"                  74\n\
biotype      9    C       \"Alanine\"                  82\n\
biotype     10    HN      \"Alanine\"                  88\n\
biotype     11    O       \"Alanine\"                  83\n\
biotype     12    HA      \"Alanine\"                   6\n\
biotype     13    CB      \"Alanine\"                   1\n\
biotype     14    HB      \"Alanine\"                   6\n\
biotype     15    N       \"Valine\"                   85\n\
biotype     16    CA      \"Valine\"                   74\n\
biotype     17    C       \"Valine\"                   82\n\
biotype     18    HN      \"Valine\"                   88\n\
biotype     19    O       \"Valine\"                   83\n\
biotype     20    HA      \"Valine\"                    6\n\
biotype     21    CB      \"Valine\"                  194\n\
biotype     22    HB      \"Valine\"                    6\n\
biotype     23    CG1     \"Valine\"                    1\n\
biotype     24    HG1     \"Valine\"                    6\n\
biotype     25    CG2     \"Valine\"                    1\n\
biotype     26    HG2     \"Valine\"                    6\n\
biotype     27    N       \"Leucine\"                  85\n\
biotype     28    CA      \"Leucine\"                  74\n\
biotype     29    C       \"Leucine\"                  82\n\
biotype     30    HN      \"Leucine\"                  88\n\
biotype     31    O       \"Leucine\"                  83\n\
biotype     32    HA      \"Leucine\"                   6\n\
biotype     33    CB      \"Leucine\"                 193\n\
biotype     34    HB      \"Leucine\"                   6\n\
biotype     35    CG      \"Leucine\"                   3\n\
biotype     36    HG      \"Leucine\"                   6\n\
biotype     37    CD1     \"Leucine\"                   1\n\
biotype     38    HD1     \"Leucine\"                   6\n\
biotype     39    CD2     \"Leucine\"                   1\n\
biotype     40    HD2     \"Leucine\"                   6\n\
biotype     41    N       \"Isoleucine\"               85\n\
biotype     42    CA      \"Isoleucine\"               74\n\
biotype     43    C       \"Isoleucine\"               82\n\
biotype     44    HN      \"Isoleucine\"               88\n\
biotype     45    O       \"Isoleucine\"               83\n\
biotype     46    HA      \"Isoleucine\"                6\n\
biotype     47    CB      \"Isoleucine\"              195\n\
biotype     48    HB      \"Isoleucine\"                6\n\
biotype     49    CG1     \"Isoleucine\"                2\n\
biotype     50    HG1     \"Isoleucine\"                6\n\
biotype     51    CG2     \"Isoleucine\"                1\n\
biotype     52    HG2     \"Isoleucine\"                6\n\
biotype     53    CD      \"Isoleucine\"                1\n\
biotype     54    HD      \"Isoleucine\"                6\n\
biotype     55    N       \"Serine\"                   85\n\
biotype     56    CA      \"Serine\"                   74\n\
biotype     57    C       \"Serine\"                   82\n\
biotype     58    HN      \"Serine\"                   88\n\
biotype     59    O       \"Serine\"                   83\n\
biotype     60    HA      \"Serine\"                    6\n\
biotype     61    CB      \"Serine\"                  189\n\
biotype     62    HB      \"Serine\"                    6\n\
biotype     63    OG      \"Serine\"                   16\n\
biotype     64    HG      \"Serine\"                   17\n\
biotype     65    N       \"Threonine\"                85\n\
biotype     66    CA      \"Threonine\"                74\n\
biotype     67    C       \"Threonine\"                82\n\
biotype     68    HN      \"Threonine\"                88\n\
biotype     69    O       \"Threonine\"                83\n\
biotype     70    HA      \"Threonine\"                 6\n\
biotype     71    CB      \"Threonine\"                20\n\
biotype     72    HB      \"Threonine\"                 6\n\
biotype     73    OG1     \"Threonine\"                16\n\
biotype     74    HG1     \"Threonine\"                17\n\
biotype     75    CG2     \"Threonine\"                 1\n\
biotype     76    HG2     \"Threonine\"                 6\n\
biotype     77    N       \"Cysteine (-SH)\"           85\n\
biotype     78    CA      \"Cysteine (-SH)\"           74\n\
biotype     79    C       \"Cysteine (-SH)\"           82\n\
biotype     80    HN      \"Cysteine (-SH)\"           88\n\
biotype     81    O       \"Cysteine (-SH)\"           83\n\
biotype     82    HA      \"Cysteine (-SH)\"            6\n\
biotype     83    CB      \"Cysteine (-SH)\"          190\n\
biotype     84    HB      \"Cysteine (-SH)\"            6\n\
biotype     85    SG      \"Cysteine (-SH)\"           52\n\
biotype     86    HG      \"Cysteine (-SH)\"           56\n\
biotype     87    N       \"Cystine (-SS-)\"           85\n\
biotype     88    CA      \"Cystine (-SS-)\"           74\n\
biotype     89    C       \"Cystine (-SS-)\"           82\n\
biotype     90    HN      \"Cystine (-SS-)\"           88\n\
biotype     91    O       \"Cystine (-SS-)\"           83\n\
biotype     92    HA      \"Cystine (-SS-)\"            6\n\
biotype     93    CB      \"Cystine (-SS-)\"           67\n\
biotype     94    HB      \"Cystine (-SS-)\"            6\n\
biotype     95    SG      \"Cystine (-SS-)\"           55\n\
biotype     96    N       \"Proline\"                  86\n\
biotype     97    CA      \"Proline\"                  93\n\
biotype     98    C       \"Proline\"                  82\n\
biotype     99    O       \"Proline\"                  83\n\
biotype    100    HA      \"Proline\"                   6\n\
biotype    101    CB      \"Proline\"                   2\n\
biotype    102    HB      \"Proline\"                   6\n\
biotype    103    CG      \"Proline\"                   2\n\
biotype    104    HG      \"Proline\"                   6\n\
biotype    105    CD      \"Proline\"                 159\n\
biotype    106    HD      \"Proline\"                   6\n\
biotype    107    N       \"Phenylalanine\"            85\n\
biotype    108    CA      \"Phenylalanine\"            74\n\
biotype    109    C       \"Phenylalanine\"            82\n\
biotype    110    HN      \"Phenylalanine\"            88\n\
biotype    111    O       \"Phenylalanine\"            83\n\
biotype    112    HA      \"Phenylalanine\"             6\n\
biotype    113    CB      \"Phenylalanine\"            15\n\
biotype    114    HB      \"Phenylalanine\"             6\n\
biotype    115    CG      \"Phenylalanine\"            11\n\
biotype    116    CD      \"Phenylalanine\"            11\n\
biotype    117    HD      \"Phenylalanine\"            12\n\
biotype    118    CE      \"Phenylalanine\"            11\n\
biotype    119    HE      \"Phenylalanine\"            12\n\
biotype    120    CZ      \"Phenylalanine\"            11\n\
biotype    121    HZ      \"Phenylalanine\"            12\n\
biotype    122    N       \"Tyrosine\"                 85\n\
biotype    123    CA      \"Tyrosine\"                 74\n\
biotype    124    C       \"Tyrosine\"                 82\n\
biotype    125    HN      \"Tyrosine\"                 88\n\
biotype    126    O       \"Tyrosine\"                 83\n\
biotype    127    HA      \"Tyrosine\"                  6\n\
biotype    128    CB      \"Tyrosine\"                 15\n\
biotype    129    HB      \"Tyrosine\"                  6\n\
biotype    130    CG      \"Tyrosine\"                 11\n\
biotype    131    CD      \"Tyrosine\"                 11\n\
biotype    132    HD      \"Tyrosine\"                 12\n\
biotype    133    CE      \"Tyrosine\"                 11\n\
biotype    134    HE      \"Tyrosine\"                 12\n\
biotype    135    CZ      \"Tyrosine\"                 28\n\
biotype    136    OH      \"Tyrosine\"                 29\n\
biotype    137    HH      \"Tyrosine\"                 30\n\
biotype    138    N       \"Tryptophan\"               85\n\
biotype    139    CA      \"Tryptophan\"               74\n\
biotype    140    C       \"Tryptophan\"               82\n\
biotype    141    HN      \"Tryptophan\"               88\n\
biotype    142    O       \"Tryptophan\"               83\n\
biotype    143    HA      \"Tryptophan\"                6\n\
biotype    144    CB      \"Tryptophan\"              197\n\
biotype    145    HB      \"Tryptophan\"                6\n\
biotype    146    CG      \"Tryptophan\"              145\n\
biotype    147    CD1     \"Tryptophan\"               11\n\
biotype    148    HD1     \"Tryptophan\"               12\n\
biotype    149    CD2     \"Tryptophan\"              146\n\
biotype    150    NE1     \"Tryptophan\"              148\n\
biotype    151    HE1     \"Tryptophan\"              149\n\
biotype    152    CE2     \"Tryptophan\"              147\n\
biotype    153    CE3     \"Tryptophan\"               11\n\
biotype    154    HE3     \"Tryptophan\"               12\n\
biotype    155    CZ2     \"Tryptophan\"               11\n\
biotype    156    HZ2     \"Tryptophan\"               12\n\
biotype    157    CZ3     \"Tryptophan\"               11\n\
biotype    158    HZ3     \"Tryptophan\"               12\n\
biotype    159    CH2     \"Tryptophan\"               11\n\
biotype    160    HH2     \"Tryptophan\"               12\n\
biotype    161    N       \"Histidine (+)\"            85\n\
biotype    162    CA      \"Histidine (+)\"            74\n\
biotype    163    C       \"Histidine (+)\"            82\n\
biotype    164    HN      \"Histidine (+)\"            88\n\
biotype    165    O       \"Histidine (+)\"            83\n\
biotype    166    HA      \"Histidine (+)\"             6\n\
biotype    167    CB      \"Histidine (+)\"           201\n\
biotype    168    HB      \"Histidine (+)\"             6\n\
biotype    169    CG      \"Histidine (+)\"           155\n\
biotype    170    ND1     \"Histidine (+)\"           157\n\
biotype    171    HD1     \"Histidine (+)\"           158\n\
biotype    172    CD2     \"Histidine (+)\"           155\n\
biotype    173    HD2     \"Histidine (+)\"            12\n\
biotype    174    CE1     \"Histidine (+)\"           154\n\
biotype    175    HE1     \"Histidine (+)\"            12\n\
biotype    176    NE2     \"Histidine (+)\"           157\n\
biotype    177    HE2     \"Histidine (+)\"           158\n\
biotype    178    N       \"Histidine (HD)\"           85\n\
biotype    179    CA      \"Histidine (HD)\"           74\n\
biotype    180    C       \"Histidine (HD)\"           82\n\
biotype    181    HN      \"Histidine (HD)\"           88\n\
biotype    182    O       \"Histidine (HD)\"           83\n\
biotype    183    HA      \"Histidine (HD)\"            6\n\
biotype    184    CB      \"Histidine (HD)\"          150\n\
biotype    185    HB      \"Histidine (HD)\"            6\n\
biotype    186    CG      \"Histidine (HD)\"          153\n\
biotype    187    ND1     \"Histidine (HD)\"          148\n\
biotype    188    HD1     \"Histidine (HD)\"          149\n\
biotype    189    CD2     \"Histidine (HD)\"          152\n\
biotype    190    HD2     \"Histidine (HD)\"           12\n\
biotype    191    CE1     \"Histidine (HD)\"          151\n\
biotype    192    HE1     \"Histidine (HD)\"           12\n\
biotype    193    NE2     \"Histidine (HD)\"          156\n\
biotype    194    N       \"Histidine (HE)\"           85\n\
biotype    195    CA      \"Histidine (HE)\"           74\n\
biotype    196    C       \"Histidine (HE)\"           82\n\
biotype    197    HN      \"Histidine (HE)\"           88\n\
biotype    198    O       \"Histidine (HE)\"           83\n\
biotype    199    HA      \"Histidine (HE)\"            6\n\
biotype    200    CB      \"Histidine (HE)\"          150\n\
biotype    201    HB      \"Histidine (HE)\"            6\n\
biotype    202    CG      \"Histidine (HE)\"          152\n\
biotype    203    ND1     \"Histidine (HE)\"          156\n\
biotype    204    CD2     \"Histidine (HE)\"          153\n\
biotype    205    HD2     \"Histidine (HE)\"           12\n\
biotype    206    CE1     \"Histidine (HE)\"          151\n\
biotype    207    HE1     \"Histidine (HE)\"           12\n\
biotype    208    NE2     \"Histidine (HE)\"          148\n\
biotype    209    HE2     \"Histidine (HE)\"          149\n\
biotype    210    N       \"Aspartic Acid\"            85\n\
biotype    211    CA      \"Aspartic Acid\"            74\n\
biotype    212    C       \"Aspartic Acid\"            82\n\
biotype    213    HN      \"Aspartic Acid\"            88\n\
biotype    214    O       \"Aspartic Acid\"            83\n\
biotype    215    HA      \"Aspartic Acid\"             6\n\
biotype    216    CB      \"Aspartic Acid\"           198\n\
biotype    217    HB      \"Aspartic Acid\"             6\n\
biotype    218    CG      \"Aspartic Acid\"           111\n\
biotype    219    OD      \"Aspartic Acid\"           112\n\
biotype    220    N       \"Asparagine\"               85\n\
biotype    221    CA      \"Asparagine\"               74\n\
biotype    222    C       \"Asparagine\"               82\n\
biotype    223    HN      \"Asparagine\"               88\n\
biotype    224    O       \"Asparagine\"               83\n\
biotype    225    HA      \"Asparagine\"                6\n\
biotype    226    CB      \"Asparagine\"              191\n\
biotype    227    HB      \"Asparagine\"                6\n\
biotype    228    CG      \"Asparagine\"               82\n\
biotype    229    OD1     \"Asparagine\"               83\n\
biotype    230    ND2     \"Asparagine\"               84\n\
biotype    231    HD2     \"Asparagine\"               87\n\
biotype    232    N       \"Glutamic Acid\"            85\n\
biotype    233    CA      \"Glutamic Acid\"            74\n\
biotype    234    C       \"Glutamic Acid\"            82\n\
biotype    235    HN      \"Glutamic Acid\"            88\n\
biotype    236    O       \"Glutamic Acid\"            83\n\
biotype    237    HA      \"Glutamic Acid\"             6\n\
biotype    238    CB      \"Glutamic Acid\"           199\n\
biotype    239    HB      \"Glutamic Acid\"             6\n\
biotype    240    CG      \"Glutamic Acid\"           114\n\
biotype    241    HG      \"Glutamic Acid\"             6\n\
biotype    242    CD      \"Glutamic Acid\"           111\n\
biotype    243    OE      \"Glutamic Acid\"           112\n\
biotype    244    N       \"Glutamine\"                85\n\
biotype    245    CA      \"Glutamine\"                74\n\
biotype    246    C       \"Glutamine\"                82\n\
biotype    247    HN      \"Glutamine\"                88\n\
biotype    248    O       \"Glutamine\"                83\n\
biotype    249    HA      \"Glutamine\"                 6\n\
biotype    250    CB      \"Glutamine\"               192\n\
biotype    251    HB      \"Glutamine\"                 6\n\
biotype    252    CG      \"Glutamine\"                 2\n\
biotype    253    HG      \"Glutamine\"                 6\n\
biotype    254    CD      \"Glutamine\"                82\n\
biotype    255    OE1     \"Glutamine\"                83\n\
biotype    256    NE2     \"Glutamine\"                84\n\
biotype    257    HE2     \"Glutamine\"                87\n\
biotype    258    N       \"Methionine\"               85\n\
biotype    259    CA      \"Methionine\"               74\n\
biotype    260    C       \"Methionine\"               82\n\
biotype    261    HN      \"Methionine\"               88\n\
biotype    262    O       \"Methionine\"               83\n\
biotype    263    HA      \"Methionine\"                6\n\
biotype    264    CB      \"Methionine\"              196\n\
biotype    265    HB      \"Methionine\"                6\n\
biotype    266    CG      \"Methionine\"               63\n\
biotype    267    HG      \"Methionine\"                6\n\
biotype    268    SD      \"Methionine\"               54\n\
biotype    269    CE      \"Methionine\"               62\n\
biotype    270    HE      \"Methionine\"                6\n\
biotype    271    N       \"Lysine\"                   85\n\
biotype    272    CA      \"Lysine\"                   74\n\
biotype    273    C       \"Lysine\"                   82\n\
biotype    274    HN      \"Lysine\"                   88\n\
biotype    275    O       \"Lysine\"                   83\n\
biotype    276    HA      \"Lysine\"                    6\n\
biotype    277    CB      \"Lysine\"                  200\n\
biotype    278    HB      \"Lysine\"                    6\n\
biotype    279    CG      \"Lysine\"                    2\n\
biotype    280    HG      \"Lysine\"                    6\n\
biotype    281    CD      \"Lysine\"                    2\n\
biotype    282    HD      \"Lysine\"                    6\n\
biotype    283    CE      \"Lysine\"                  133\n\
biotype    284    HE      \"Lysine\"                    6\n\
biotype    285    NZ      \"Lysine\"                  124\n\
biotype    286    HZ      \"Lysine\"                  129\n\
biotype    287    N       \"Arginine\"                 85\n\
biotype    288    CA      \"Arginine\"                 74\n\
biotype    289    C       \"Arginine\"                 82\n\
biotype    290    HN      \"Arginine\"                 88\n\
biotype    291    O       \"Arginine\"                 83\n\
biotype    292    HA      \"Arginine\"                  6\n\
biotype    293    CB      \"Arginine\"                202\n\
biotype    294    HB      \"Arginine\"                  6\n\
biotype    295    CG      \"Arginine\"                144\n\
biotype    296    HG      \"Arginine\"                  6\n\
biotype    297    CD      \"Arginine\"                143\n\
biotype    298    HD      \"Arginine\"                  6\n\
biotype    299    NE      \"Arginine\"                137\n\
biotype    300    HE      \"Arginine\"                139\n\
biotype    301    CZ      \"Arginine\"                140\n\
biotype    302    NH      \"Arginine\"                136\n\
biotype    303    HH      \"Arginine\"                138\n\
biotype    304    N       \"Ornithine\"                85\n\
biotype    305    CA      \"Ornithine\"                74\n\
biotype    306    C       \"Ornithine\"                82\n\
biotype    307    HN      \"Ornithine\"                88\n\
biotype    308    O       \"Ornithine\"                83\n\
biotype    309    HA      \"Ornithine\"                 6\n\
biotype    310    CB      \"Ornithine\"                 2\n\
biotype    311    HB      \"Ornithine\"                 6\n\
biotype    312    CG      \"Ornithine\"                 2\n\
biotype    313    HG      \"Ornithine\"                 6\n\
biotype    314    CD      \"Ornithine\"               133\n\
biotype    315    HD      \"Ornithine\"                 6\n\
biotype    316    NE      \"Ornithine\"               124\n\
biotype    317    HE      \"Ornithine\"               129\n\
biotype    318    N       \"MethylAlanine (AIB)\"      85\n\
biotype    319    CA      \"MethylAlanine (AIB)\"      75\n\
biotype    320    C       \"MethylAlanine (AIB)\"      82\n\
biotype    321    HN      \"MethylAlanine (AIB)\"      88\n\
biotype    322    O       \"MethylAlanine (AIB)\"      83\n\
biotype    323    CB      \"MethylAlanine (AIB)\"       1\n\
biotype    324    HB      \"MethylAlanine (AIB)\"       6\n\
biotype    325    N       \"Pyroglutamic Acid\"        85\n\
biotype    326    CA      \"Pyroglutamic Acid\"        74\n\
biotype    327    C       \"Pyroglutamic Acid\"        82\n\
biotype    328    HN      \"Pyroglutamic Acid\"        88\n\
biotype    329    O       \"Pyroglutamic Acid\"        83\n\
biotype    330    HA      \"Pyroglutamic Acid\"         6\n\
biotype    331    CB      \"Pyroglutamic Acid\"         2\n\
biotype    332    HB      \"Pyroglutamic Acid\"         6\n\
biotype    333    CG      \"Pyroglutamic Acid\"         2\n\
biotype    334    HG      \"Pyroglutamic Acid\"         6\n\
biotype    335    CD      \"Pyroglutamic Acid\"        82\n\
biotype    336    OE      \"Pyroglutamic Acid\"        83\n\
biotype    337    C       \"Formyl N-Terminus\"        82\n\
biotype    338    H       \"Formyl N-Terminus\"       119\n\
biotype    339    O       \"Formyl N-Terminus\"        83\n\
biotype    340    CH3     \"Acetyl N-Terminus\"         1\n\
biotype    341    H       \"Acetyl N-Terminus\"         6\n\
biotype    342    C       \"Acetyl N-Terminus\"        82\n\
biotype    343    O       \"Acetyl N-Terminus\"        83\n\
biotype    344    N       \"Amide C-Terminus\"         84\n\
biotype    345    HN      \"Amide C-Terminus\"         87\n\
biotype    346    N       \"N-MeAmide C-Terminus\"     85\n\
biotype    347    HN      \"N-MeAmide C-Terminus\"     88\n\
biotype    348    CH3     \"N-MeAmide C-Terminus\"     89\n\
biotype    349    H       \"N-MeAmide C-Terminus\"      6\n\
biotype    350    N       \"N-Terminal GLY\"          124\n\
biotype    351    CA      \"N-Terminal GLY\"          133\n\
biotype    352    C       \"N-Terminal GLY\"           82\n\
biotype    353    HN      \"N-Terminal GLY\"          129\n\
biotype    354    O       \"N-Terminal GLY\"           83\n\
biotype    355    HA      \"N-Terminal GLY\"            6\n\
biotype    356    N       \"N-Terminal ALA\"          124\n\
biotype    357    CA      \"N-Terminal ALA\"          134\n\
biotype    358    C       \"N-Terminal ALA\"           82\n\
biotype    359    HN      \"N-Terminal ALA\"          129\n\
biotype    360    O       \"N-Terminal ALA\"           83\n\
biotype    361    HA      \"N-Terminal ALA\"            6\n\
biotype    362    N       \"N-Terminal VAL\"          124\n\
biotype    363    CA      \"N-Terminal VAL\"          134\n\
biotype    364    C       \"N-Terminal VAL\"           82\n\
biotype    365    HN      \"N-Terminal VAL\"          129\n\
biotype    366    O       \"N-Terminal VAL\"           83\n\
biotype    367    HA      \"N-Terminal VAL\"            6\n\
biotype    368    N       \"N-Terminal LEU\"          124\n\
biotype    369    CA      \"N-Terminal LEU\"          134\n\
biotype    370    C       \"N-Terminal LEU\"           82\n\
biotype    371    HN      \"N-Terminal LEU\"          129\n\
biotype    372    O       \"N-Terminal LEU\"           83\n\
biotype    373    HA      \"N-Terminal LEU\"            6\n\
biotype    374    N       \"N-Terminal ILE\"          124\n\
biotype    375    CA      \"N-Terminal ILE\"          134\n\
biotype    376    C       \"N-Terminal ILE\"           82\n\
biotype    377    HN      \"N-Terminal ILE\"          129\n\
biotype    378    O       \"N-Terminal ILE\"           83\n\
biotype    379    HA      \"N-Terminal ILE\"            6\n\
biotype    380    N       \"N-Terminal SER\"          124\n\
biotype    381    CA      \"N-Terminal SER\"          134\n\
biotype    382    C       \"N-Terminal SER\"           82\n\
biotype    383    HN      \"N-Terminal SER\"          129\n\
biotype    384    O       \"N-Terminal SER\"           83\n\
biotype    385    HA      \"N-Terminal SER\"            6\n\
biotype    386    N       \"N-Terminal THR\"          124\n\
biotype    387    CA      \"N-Terminal THR\"          134\n\
biotype    388    C       \"N-Terminal THR\"           82\n\
biotype    389    HN      \"N-Terminal THR\"          129\n\
biotype    390    O       \"N-Terminal THR\"           83\n\
biotype    391    HA      \"N-Terminal THR\"            6\n\
biotype    392    N       \"N-Terminal CYS (-SH)\"    124\n\
biotype    393    CA      \"N-Terminal CYS (-SH)\"    134\n\
biotype    394    C       \"N-Terminal CYS (-SH)\"     82\n\
biotype    395    HN      \"N-Terminal CYS (-SH)\"    129\n\
biotype    396    O       \"N-Terminal CYS (-SH)\"     83\n\
biotype    397    HA      \"N-Terminal CYS (-SH)\"      6\n\
biotype    398    N       \"N-Terminal CYS (-SS)\"    124\n\
biotype    399    CA      \"N-Terminal CYS (-SS)\"    134\n\
biotype    400    C       \"N-Terminal CYS (-SS)\"     82\n\
biotype    401    HN      \"N-Terminal CYS (-SS)\"    129\n\
biotype    402    O       \"N-Terminal CYS (-SS)\"     83\n\
biotype    403    HA      \"N-Terminal CYS (-SS)\"      6\n\
biotype    404    N       \"N-Terminal PRO\"          125\n\
biotype    405    CA      \"N-Terminal PRO\"          134\n\
biotype    406    C       \"N-Terminal PRO\"           82\n\
biotype    407    HN      \"N-Terminal PRO\"          130\n\
biotype    408    O       \"N-Terminal PRO\"           83\n\
biotype    409    HA      \"N-Terminal PRO\"            6\n\
biotype    410    CD      \"N-Terminal PRO\"          133\n\
biotype    411    HD      \"N-Terminal PRO\"            6\n\
biotype    412    N       \"N-Terminal PHE\"          124\n\
biotype    413    CA      \"N-Terminal PHE\"          134\n\
biotype    414    C       \"N-Terminal PHE\"           82\n\
biotype    415    HN      \"N-Terminal PHE\"          129\n\
biotype    416    O       \"N-Terminal PHE\"           83\n\
biotype    417    HA      \"N-Terminal PHE\"            6\n\
biotype    418    N       \"N-Terminal TYR\"          124\n\
biotype    419    CA      \"N-Terminal TYR\"          134\n\
biotype    420    C       \"N-Terminal TYR\"           82\n\
biotype    421    HN      \"N-Terminal TYR\"          129\n\
biotype    422    O       \"N-Terminal TYR\"           83\n\
biotype    423    HA      \"N-Terminal TYR\"            6\n\
biotype    424    N       \"N-Terminal TRP\"          124\n\
biotype    425    CA      \"N-Terminal TRP\"          134\n\
biotype    426    C       \"N-Terminal TRP\"           82\n\
biotype    427    HN      \"N-Terminal TRP\"          129\n\
biotype    428    O       \"N-Terminal TRP\"           83\n\
biotype    429    HA      \"N-Terminal TRP\"            6\n\
biotype    430    N       \"N-Terminal HIS (+)\"      124\n\
biotype    431    CA      \"N-Terminal HIS (+)\"      134\n\
biotype    432    C       \"N-Terminal HIS (+)\"       82\n\
biotype    433    HN      \"N-Terminal HIS (+)\"      129\n\
biotype    434    O       \"N-Terminal HIS (+)\"       83\n\
biotype    435    HA      \"N-Terminal HIS (+)\"        6\n\
biotype    436    N       \"N-Terminal HIS (HD)\"     124\n\
biotype    437    CA      \"N-Terminal HIS (HD)\"     134\n\
biotype    438    C       \"N-Terminal HIS (HD)\"      82\n\
biotype    439    HN      \"N-Terminal HIS (HD)\"     129\n\
biotype    440    O       \"N-Terminal HIS (HD)\"      83\n\
biotype    441    HA      \"N-Terminal HIS (HD)\"       6\n\
biotype    442    N       \"N-Terminal HIS (HE)\"     124\n\
biotype    443    CA      \"N-Terminal HIS (HE)\"     134\n\
biotype    444    C       \"N-Terminal HIS (HE)\"      82\n\
biotype    445    HN      \"N-Terminal HIS (HE)\"     129\n\
biotype    446    O       \"N-Terminal HIS (HE)\"      83\n\
biotype    447    HA      \"N-Terminal HIS (HE)\"       6\n\
biotype    448    N       \"N-Terminal ASP\"          124\n\
biotype    449    CA      \"N-Terminal ASP\"          134\n\
biotype    450    C       \"N-Terminal ASP\"           82\n\
biotype    451    HN      \"N-Terminal ASP\"          129\n\
biotype    452    O       \"N-Terminal ASP\"           83\n\
biotype    453    HA      \"N-Terminal ASP\"            6\n\
biotype    454    N       \"N-Terminal ASN\"          124\n\
biotype    455    CA      \"N-Terminal ASN\"          134\n\
biotype    456    C       \"N-Terminal ASN\"           82\n\
biotype    457    HN      \"N-Terminal ASN\"          129\n\
biotype    458    O       \"N-Terminal ASN\"           83\n\
biotype    459    HA      \"N-Terminal ASN\"            6\n\
biotype    460    N       \"N-Terminal GLU\"          124\n\
biotype    461    CA      \"N-Terminal GLU\"          134\n\
biotype    462    C       \"N-Terminal GLU\"           82\n\
biotype    463    HN      \"N-Terminal GLU\"          129\n\
biotype    464    O       \"N-Terminal GLU\"           83\n\
biotype    465    HA      \"N-Terminal GLU\"            6\n\
biotype    466    N       \"N-Terminal GLN\"          124\n\
biotype    467    CA      \"N-Terminal GLN\"          134\n\
biotype    468    C       \"N-Terminal GLN\"           82\n\
biotype    469    HN      \"N-Terminal GLN\"          129\n\
biotype    470    O       \"N-Terminal GLN\"           83\n\
biotype    471    HA      \"N-Terminal GLN\"            6\n\
biotype    472    N       \"N-Terminal MET\"          124\n\
biotype    473    CA      \"N-Terminal MET\"          134\n\
biotype    474    C       \"N-Terminal MET\"           82\n\
biotype    475    HN      \"N-Terminal MET\"          129\n\
biotype    476    O       \"N-Terminal MET\"           83\n\
biotype    477    HA      \"N-Terminal MET\"            6\n\
biotype    478    N       \"N-Terminal LYS\"          124\n\
biotype    479    CA      \"N-Terminal LYS\"          134\n\
biotype    480    C       \"N-Terminal LYS\"           82\n\
biotype    481    HN      \"N-Terminal LYS\"          129\n\
biotype    482    O       \"N-Terminal LYS\"           83\n\
biotype    483    HA      \"N-Terminal LYS\"            6\n\
biotype    484    N       \"N-Terminal ARG\"          124\n\
biotype    485    CA      \"N-Terminal ARG\"          134\n\
biotype    486    C       \"N-Terminal ARG\"           82\n\
biotype    487    HN      \"N-Terminal ARG\"          129\n\
biotype    488    O       \"N-Terminal ARG\"           83\n\
biotype    489    HA      \"N-Terminal ARG\"            6\n\
biotype    490    N       \"N-Terminal ORN\"          124\n\
biotype    491    CA      \"N-Terminal ORN\"          134\n\
biotype    492    C       \"N-Terminal ORN\"           82\n\
biotype    493    HN      \"N-Terminal ORN\"          129\n\
biotype    494    O       \"N-Terminal ORN\"           83\n\
biotype    495    HA      \"N-Terminal ORN\"            6\n\
biotype    496    N       \"N-Terminal AIB\"          124\n\
biotype    497    CA      \"N-Terminal AIB\"          135\n\
biotype    498    C       \"N-Terminal AIB\"           82\n\
biotype    499    HN      \"N-Terminal AIB\"          129\n\
biotype    500    O       \"N-Terminal AIB\"           83\n\
biotype    501    N       \"C-Terminal GLY\"           85\n\
biotype    502    CA      \"C-Terminal GLY\"          160\n\
biotype    503    C       \"C-Terminal GLY\"          111\n\
biotype    504    HN      \"C-Terminal GLY\"           88\n\
biotype    505    OXT     \"C-Terminal GLY\"          112\n\
biotype    506    HA      \"C-Terminal GLY\"            6\n\
biotype    507    N       \"C-Terminal ALA\"           85\n\
biotype    508    CA      \"C-Terminal ALA\"          161\n\
biotype    509    C       \"C-Terminal ALA\"          111\n\
biotype    510    HN      \"C-Terminal ALA\"           88\n\
biotype    511    OXT     \"C-Terminal ALA\"          112\n\
biotype    512    HA      \"C-Terminal ALA\"            6\n\
biotype    513    N       \"C-Terminal VAL\"           85\n\
biotype    514    CA      \"C-Terminal VAL\"          161\n\
biotype    515    C       \"C-Terminal VAL\"          111\n\
biotype    516    HN      \"C-Terminal VAL\"           88\n\
biotype    517    OXT     \"C-Terminal VAL\"          112\n\
biotype    518    HA      \"C-Terminal VAL\"            6\n\
biotype    519    N       \"C-Terminal LEU\"           85\n\
biotype    520    CA      \"C-Terminal LEU\"          161\n\
biotype    521    C       \"C-Terminal LEU\"          111\n\
biotype    522    HN      \"C-Terminal LEU\"           88\n\
biotype    523    OXT     \"C-Terminal LEU\"          112\n\
biotype    524    HA      \"C-Terminal LEU\"            6\n\
biotype    525    N       \"C-Terminal ILE\"           85\n\
biotype    526    CA      \"C-Terminal ILE\"          161\n\
biotype    527    C       \"C-Terminal ILE\"          111\n\
biotype    528    HN      \"C-Terminal ILE\"           88\n\
biotype    529    OXT     \"C-Terminal ILE\"          112\n\
biotype    530    HA      \"C-Terminal ILE\"            6\n\
biotype    531    N       \"C-Terminal SER\"           85\n\
biotype    532    CA      \"C-Terminal SER\"          161\n\
biotype    533    C       \"C-Terminal SER\"          111\n\
biotype    534    HN      \"C-Terminal SER\"           88\n\
biotype    535    OXT     \"C-Terminal SER\"          112\n\
biotype    536    HA      \"C-Terminal SER\"            6\n\
biotype    537    N       \"C-Terminal THR\"           85\n\
biotype    538    CA      \"C-Terminal THR\"          161\n\
biotype    539    C       \"C-Terminal THR\"          111\n\
biotype    540    HN      \"C-Terminal THR\"           88\n\
biotype    541    OXT     \"C-Terminal THR\"          112\n\
biotype    542    HA      \"C-Terminal THR\"            6\n\
biotype    543    N       \"C-Terminal CYS (-SH)\"     85\n\
biotype    544    CA      \"C-Terminal CYS (-SH)\"    161\n\
biotype    545    C       \"C-Terminal CYS (-SH)\"    111\n\
biotype    546    HN      \"C-Terminal CYS (-SH)\"     88\n\
biotype    547    OXT     \"C-Terminal CYS (-SH)\"    112\n\
biotype    548    HA      \"C-Terminal CYS (-SH)\"      6\n\
biotype    549    N       \"C-Terminal CYS (-SS)\"     85\n\
biotype    550    CA      \"C-Terminal CYS (-SS)\"    161\n\
biotype    551    C       \"C-Terminal CYS (-SS)\"    111\n\
biotype    552    HN      \"C-Terminal CYS (-SS)\"     88\n\
biotype    553    OXT     \"C-Terminal CYS (-SS)\"    112\n\
biotype    554    HA      \"C-Terminal CYS (-SS)\"      6\n\
biotype    555    N       \"C-Terminal PRO\"           86\n\
biotype    556    CA      \"C-Terminal PRO\"          163\n\
biotype    557    C       \"C-Terminal PRO\"          111\n\
biotype    558    OXT     \"C-Terminal PRO\"          112\n\
biotype    559    HA      \"C-Terminal PRO\"            6\n\
biotype    560    N       \"C-Terminal PHE\"           85\n\
biotype    561    CA      \"C-Terminal PHE\"          161\n\
biotype    562    C       \"C-Terminal PHE\"          111\n\
biotype    563    HN      \"C-Terminal PHE\"           88\n\
biotype    564    OXT     \"C-Terminal PHE\"          112\n\
biotype    565    HA      \"C-Terminal PHE\"            6\n\
biotype    566    N       \"C-Terminal TYR\"           85\n\
biotype    567    CA      \"C-Terminal TYR\"          161\n\
biotype    568    C       \"C-Terminal TYR\"          111\n\
biotype    569    HN      \"C-Terminal TYR\"           88\n\
biotype    570    OXT     \"C-Terminal TYR\"          112\n\
biotype    571    HA      \"C-Terminal TYR\"            6\n\
biotype    572    N       \"C-Terminal TRP\"           85\n\
biotype    573    CA      \"C-Terminal TRP\"          161\n\
biotype    574    C       \"C-Terminal TRP\"          111\n\
biotype    575    HN      \"C-Terminal TRP\"           88\n\
biotype    576    OXT     \"C-Terminal TRP\"          112\n\
biotype    577    HA      \"C-Terminal TRP\"            6\n\
biotype    578    N       \"C-Terminal HIS (+)\"       85\n\
biotype    579    CA      \"C-Terminal HIS (+)\"      161\n\
biotype    580    C       \"C-Terminal HIS (+)\"      111\n\
biotype    581    HN      \"C-Terminal HIS (+)\"       88\n\
biotype    582    OXT     \"C-Terminal HIS (+)\"      112\n\
biotype    583    HA      \"C-Terminal HIS (+)\"        6\n\
biotype    584    N       \"C-Terminal HIS (HD)\"      85\n\
biotype    585    CA      \"C-Terminal HIS (HD)\"     161\n\
biotype    586    C       \"C-Terminal HIS (HD)\"     111\n\
biotype    587    HN      \"C-Terminal HIS (HD)\"      88\n\
biotype    588    OXT     \"C-Terminal HIS (HD)\"     112\n\
biotype    589    HA      \"C-Terminal HIS (HD)\"       6\n\
biotype    590    N       \"C-Terminal HIS (HE)\"      85\n\
biotype    591    CA      \"C-Terminal HIS (HE)\"     161\n\
biotype    592    C       \"C-Terminal HIS (HE)\"     111\n\
biotype    593    HN      \"C-Terminal HIS (HE)\"      88\n\
biotype    594    OXT     \"C-Terminal HIS (HE)\"     112\n\
biotype    595    HA      \"C-Terminal HIS (HE)\"       6\n\
biotype    596    N       \"C-Terminal ASP\"           85\n\
biotype    597    CA      \"C-Terminal ASP\"          161\n\
biotype    598    C       \"C-Terminal ASP\"          111\n\
biotype    599    HN      \"C-Terminal ASP\"           88\n\
biotype    600    OXT     \"C-Terminal ASP\"          112\n\
biotype    601    HA      \"C-Terminal ASP\"            6\n\
biotype    602    N       \"C-Terminal ASN\"           85\n\
biotype    603    CA      \"C-Terminal ASN\"          161\n\
biotype    604    C       \"C-Terminal ASN\"          111\n\
biotype    605    HN      \"C-Terminal ASN\"           88\n\
biotype    606    OXT     \"C-Terminal ASN\"          112\n\
biotype    607    HA      \"C-Terminal ASN\"            6\n\
biotype    608    N       \"C-Terminal GLU\"           85\n\
biotype    609    CA      \"C-Terminal GLU\"          161\n\
biotype    610    C       \"C-Terminal GLU\"          111\n\
biotype    611    HN      \"C-Terminal GLU\"           88\n\
biotype    612    OXT     \"C-Terminal GLU\"          112\n\
biotype    613    HA      \"C-Terminal GLU\"            6\n\
biotype    614    N       \"C-Terminal GLN\"           85\n\
biotype    615    CA      \"C-Terminal GLN\"          161\n\
biotype    616    C       \"C-Terminal GLN\"          111\n\
biotype    617    HN      \"C-Terminal GLN\"           88\n\
biotype    618    OXT     \"C-Terminal GLN\"          112\n\
biotype    619    HA      \"C-Terminal GLN\"            6\n\
biotype    620    N       \"C-Terminal MET\"           85\n\
biotype    621    CA      \"C-Terminal MET\"          161\n\
biotype    622    C       \"C-Terminal MET\"          111\n\
biotype    623    HN      \"C-Terminal MET\"           88\n\
biotype    624    OXT     \"C-Terminal MET\"          112\n\
biotype    625    HA      \"C-Terminal MET\"            6\n\
biotype    626    N       \"C-Terminal LYS\"           85\n\
biotype    627    CA      \"C-Terminal LYS\"          161\n\
biotype    628    C       \"C-Terminal LYS\"          111\n\
biotype    629    HN      \"C-Terminal LYS\"           88\n\
biotype    630    OXT     \"C-Terminal LYS\"          112\n\
biotype    631    HA      \"C-Terminal LYS\"            6\n\
biotype    632    N       \"C-Terminal ARG\"           85\n\
biotype    633    CA      \"C-Terminal ARG\"          161\n\
biotype    634    C       \"C-Terminal ARG\"          111\n\
biotype    635    HN      \"C-Terminal ARG\"           88\n\
biotype    636    OXT     \"C-Terminal ARG\"          112\n\
biotype    637    HA      \"C-Terminal ARG\"            6\n\
biotype    638    N       \"C-Terminal ORN\"           85\n\
biotype    639    CA      \"C-Terminal ORN\"          161\n\
biotype    640    C       \"C-Terminal ORN\"          111\n\
biotype    641    HN      \"C-Terminal ORN\"           88\n\
biotype    642    OXT     \"C-Terminal ORN\"          112\n\
biotype    643    HA      \"C-Terminal ORN\"            6\n\
biotype    644    N       \"C-Terminal AIB\"           85\n\
biotype    645    CA      \"C-Terminal AIB\"          162\n\
biotype    646    C       \"C-Terminal AIB\"          111\n\
biotype    647    HN      \"C-Terminal AIB\"           88\n\
biotype    648    OXT     \"C-Terminal AIB\"          112\n\
";

#endif
