// rmsd.h --- rmsd energy term
// Copyright (C) 2008-2011 Wouter Boomsma
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


#ifndef RMSD_H
#define RMSD_H

#include "utils/vector_matrix_3d.h"
#include "energy/energy_term.h"
#include "protein/iterators/atom_iterator.h"

namespace phaistos {

//! Calculate RMSD within range specified by iterators.
//! \param begin1 Iterator start-point chain1
//! \param end1 Iterator end-point chain1
//! \param begin2 Iterator start-point chain2
//! \param end2 Iterator end-point chain2
template <class ITERATORTYPE>
double calc_rmsd_from_iterators(ITERATORTYPE begin1, ITERATORTYPE end1,
                                ITERATORTYPE begin2, ITERATORTYPE end2) {
     Matrix_3D rotationMatrix;
     Vector_3D singularValues;
     Vector_3D CM1, CM2;
     calc_superimpose_rotation_matrix(begin1, end1, begin2, end2, CM1, CM2, rotationMatrix, singularValues);

     // Calculate rmsd directly from singular values
     double sum = 0;
     int length = 0;
     ITERATORTYPE it1=begin1, it2=begin2;
     for (; it1 != end1 && it2 != end2; ++it1, ++it2) {
	  
          length++;
          for (int j=0; j<3; j++) {
               double value1 = (*it1)[j] - CM1[j];
               double value2 = (*it2)[j] - CM2[j];
               sum += value1*value1 + value2*value2;
          }
     }

     // Check that compared sequences are of equal length
     if (it1!=end1 || it2!=end2) {
	  std::cerr << "WARNING (RMSD): compared sequences are not of the same length.\n";
     }
     
     double singularValueSum = 0;
     for (int i=0; i<3; i++) {
          singularValueSum += singularValues[i];
     }

     double rmsd = sqrt(std::fabs(sum - 2*singularValueSum) / length);

     return rmsd;
}


//! Calculate RMSD within range specified by iterator
//! \param it1 Iterator chain1 (iterator contains both start and end information)
//! \param it2 Iterator chain2 (iterator contains both start and end information)
template<typename ITERATORTYPE>
double calc_rmsd_from_iterators(ITERATORTYPE it1, ITERATORTYPE it2) {

     Matrix_3D rotationMatrix;
     Vector_3D singularValues;
     Vector_3D CM1, CM2;
     calc_superimpose_rotation_matrix(it1, it2, CM1, CM2, rotationMatrix, singularValues);

     // Calculate rmsd directly from singular values
     double sum = 0;
     int length = 0;
     
     for (; !it1.end() && !it2.end(); ++it1, ++it2) {

          length++;
          for (int j=0; j<3; j++) {
               double value1 = (*it1)[j] - CM1[j];
               double value2 = (*it2)[j] - CM2[j];
               sum += value1*value1 + value2*value2;
          }
     }

     // Check that compared sequences are of equal length
     if (!it1.end() || !it2.end()) {
	  std::cerr << "WARNING (RMSD): compared sequences are not of the same length.\n";
     }
     
     double singularValueSum = 0;
     for (int i=0; i<3; i++) {
          singularValueSum += singularValues[i];
     }

     double rmsd = sqrt(std::fabs(sum - 2*singularValueSum) / length);

     return rmsd;
}

//! Calculate RMSD between chains
template <definitions::IterateEnum ITERATION_MODE>
double calc_rmsd(const ChainFB &chain1, const ChainFB &chain2, int start_index=0, int end_index=-1) {
     return calc_rmsd_from_iterators(AtomIterator<ChainFB,ITERATION_MODE,Vector_3D>(chain1, start_index, end_index),
                                     AtomIterator<ChainFB,ITERATION_MODE,Vector_3D>(chain2, start_index, end_index));
}

//! Calculate RMSD between chains
template <definitions::IterateEnum ITERATION_MODE>
double calc_rmsd(const ChainCA &chain1, const ChainCA &chain2, int start_index=0, int end_index=-1) {
     return calc_rmsd_from_iterators(AtomIterator<ChainCA,ITERATION_MODE,Vector_3D>(chain1, start_index, end_index),
                                     AtomIterator<ChainCA,ITERATION_MODE,Vector_3D>(chain2, start_index, end_index));
}

//! Calculate RMSD within range specified by iterators (without superimposing).
//! \param begin1 Iterator start-point chain1
//! \param end1 Iterator end-point chain1
//! \param begin2 Iterator start-point chain2
//! \param end2 Iterator end-point chain2
template <class ITERATORTYPE>
double calc_rmsd_no_superimpose(ITERATORTYPE begin1, ITERATORTYPE end1,
                                ITERATORTYPE begin2, ITERATORTYPE end2){

     double sum=0.0;
     int length = 0;
     for (ITERATORTYPE it1(begin1),it2(begin2);
          it1 != end1 && it2 != end2; ++it1, ++it2) {
          length++;
	  sum += (*it1-*it2)*(*it1-*it2);
     }
     return sqrt(sum/length);     
}

// Calculate RMSD within range specified by iterators (without superimposing).
//! \param it1 Iterator chain1 (iterator contains both start and end information)
//! \param it2 Iterator chain2 (iterator contains both start and end information)
template <class ITERATORTYPE>
double calc_rmsd_no_superimpose(ITERATORTYPE it1, ITERATORTYPE it2){
     double sum=0.0;
     int length = 0;
     for (; !it1.end() && !it2.end(); ++it1, ++it2) {
          length++;
	  sum += (*it1-*it2)*(*it1-*it2);
     }
     return sqrt(sum/length);     
}

//! Calculate RMSD between two chains (without superimposing).
//! \param chain1 Molecule chain1
//! \param chain2 Molecule chain2
template <definitions::IterateEnum ITERATION_MODE>
double calc_rmsd_no_superimpose(const ChainFB &chain1, const ChainFB &chain2) {
     return calc_rmsd_no_superimpose(AtomIterator<ChainFB,ITERATION_MODE,Vector_3D>(chain1),
                                     AtomIterator<ChainFB,ITERATION_MODE,Vector_3D>(chain2));
}

//! Calculate RMSD between two chains (without superimposing).
//! \param chain1 Molecule chain1
//! \param chain2 Molecule chain2
template <definitions::IterateEnum ITERATION_MODE>
double calc_rmsd_no_superimpose(const ChainCA &chain1, const ChainCA &chain2) {
     return calc_rmsd_no_superimpose(AtomIterator<ChainCA,ITERATION_MODE,Vector_3D>(chain1),
                                     AtomIterator<ChainCA,ITERATION_MODE,Vector_3D>(chain2));
}



//! RMSD energy term
template <typename CHAIN_TYPE>
class TermRmsd: public EnergyTermCommon<TermRmsd<CHAIN_TYPE>, CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermRmsd<CHAIN_TYPE>, CHAIN_TYPE> EnergyTermCommon;

     //! Reference chain to which RMSD is measured 
     CHAIN_TYPE *reference_chain;

     //! Function pointer used to select iteration mode     
     double (*energy_function)(const CHAIN_TYPE &chain1, const CHAIN_TYPE &chain2, int START, int END);
     
public:

     //! Local settings class.     
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          
          //! PDB of Reference chain to which RMSD is measured 
          std::string reference_pdb_file;

          //! Whether to calculate CA-only RMSD          
          bool ca_only;
          
          //! Start residue:
          int residue_start;

          //! End residue:
          int residue_end;

          //! Constructor. Defines default values for settings object.          
          Settings(std::string reference_pdb_file="",
                   bool ca_only=true,
                   int residue_start=0,
                   int residue_end=-1)
               : reference_pdb_file(reference_pdb_file),
                 ca_only(ca_only),
                 residue_start(residue_start),
                 residue_end(residue_end) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "reference-pdb-file:" << settings.reference_pdb_file << "\n";
               o << "ca-only:" << settings.ca_only << "\n";
               o << "residue-start:" << settings.residue_start << "\n";
               o << "residue-end:" << settings.residue_end << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }                    
     } settings;   //!< Local settings object


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object     
     //! \param random_number_engine Object from which random number generators can be created.
     TermRmsd(CHAIN_TYPE *chain, const Settings &settings=Settings(),
              RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "rmsd", settings, random_number_engine),
            settings(settings) {

          reference_chain = new CHAIN_TYPE(settings.reference_pdb_file);
      
	  if (settings.ca_only) {
	       energy_function = &calc_rmsd<definitions::CA_ONLY>;
	  } else {
	       energy_function = &calc_rmsd<definitions::ALL>;
	  }
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermRmsd(const TermRmsd &other, 
              RandomNumberEngine *random_number_engine,
              int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            energy_function(other.energy_function),
            settings(other.settings) {

          reference_chain = new CHAIN_TYPE(settings.reference_pdb_file);

     }          

     //! Destructor
     ~TermRmsd() {
          delete reference_chain;
     }

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {
	  return (*energy_function)(*this->chain, *reference_chain,
                                    settings.residue_start,
                                    settings.residue_end);
     }
};

}

#endif
