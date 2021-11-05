// move_info.h --- Information about last move
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

#include "protein/definitions.h"

#ifndef MOVE_INFO_H
#define MOVE_INFO_H

namespace phaistos {

//! Information about the last executed move
class MoveInfo {
public:

     //! Flag indicating whether move was successful
     bool success;
     
     //! Which angles are modified by the current move - vector of (start,end) pairs
     std::vector<std::pair<int,int> > modified_angles;

     //! Total range for modified angles - start index
     int modified_angles_start;

     //! Total range for modified angles - end index
     int modified_angles_end;

     //! Which positions are modified by the current move - vector of (start,end) pairs
     std::vector<std::pair<int,int> > modified_positions;

     //! Total range for modified positions - start index
     int modified_positions_start;

     //! Total range for modified positions - end index
     int modified_positions_end;

     //! Type of move
     definitions::MoveTypeEnum move_type;

     //! Constructor
     //! \param move_type Type of move
     MoveInfo(definitions::MoveTypeEnum move_type=definitions::NON_LOCAL)
          : success(true),
            modified_angles_start(uninitialized<int>()),
            modified_angles_end(uninitialized<int>()),
            move_type(move_type){}


     //! Add new information
     //! \param modified_angle_range [start_index,end_index] pair specifying for which residues backbone angles have changed
     //! \param modified_position_range [start_index,end_index] pair specifying for which residues positions have changed
     //! \returns The object itself (this)
     MoveInfo *add_info(const std::pair<int, int> &modified_angle_range,
                        const std::pair<int, int> &modified_position_range) {

          // Set overall start end values (combining values from vector)
          if (modified_angles.size() == 0) {
               modified_angles_start = modified_angle_range.first;
               modified_angles_end   = modified_angle_range.second;
               modified_positions_start = modified_position_range.first;
               modified_positions_end   = modified_position_range.second;
          } else {
               modified_angles_start = std::min(modified_angles_start, modified_angle_range.first);
               modified_angles_end   = std::max(modified_angles_end, modified_angle_range.second);
               modified_positions_start = std::min(modified_positions_start, modified_position_range.first);
               modified_positions_end   = std::max(modified_positions_end, modified_position_range.second);
          }

          modified_angles.push_back(modified_angle_range);
          modified_positions.push_back(modified_position_range);

          return this;
     }


     // Add new information using a vector of modified pairs
     //! \param modified_angles vector of [start_index,end_index] pairs specifying for which residues backbone angles have changed
     //! \param modified_positions vector of [start_index,end_index] pairs specifying for which residues positions have changed
     //! \returns The object itself (this)
     MoveInfo *add_info(const std::vector<std::pair<int, int> > &modified_angles,
                        const std::vector<std::pair<int, int> > &modified_positions) {

          assert(modified_angles.size() == modified_positions.size());

          for (unsigned int i=0; i<modified_angles.size(); ++i) {
               this->add_info(modified_angles[i],
                              modified_positions[i]);
          }
          return this;
     }

     //! Destructor
     virtual ~MoveInfo(){}

     //! Output operator
     friend std::ostream &operator<<(std::ostream &o, const MoveInfo &mi) {
          o << "Modified angles: " << mi.modified_angles << "\n";
          o << "Modified angles - total range: " << mi.modified_angles_start << "-" << mi.modified_angles_end << "\n";
          o << "Modified positions: " << mi.modified_positions << "\n";
          o << "Move type: " << mi.move_type << "\n";
          o << "Success: " << mi.success << "\n";
          return o;
     }                                             
};


//! Information about the last executed move - derived class used by moves
//! which contain a probabilistic model. This makes it possible to distinguish
//! between move types based on the associated move_info object, a feature
//! which is used by Move::notify()
template <typename DBN_TYPE>
class MoveInfoDbn: public MoveInfo {
public:

     //! Dynamic Bayesian network model
     DBN_TYPE *dbn;

     //! Constructor
     //! \param dbn Dynamic Bayesian network model
     //! \param move_type Type of move
     MoveInfoDbn(DBN_TYPE *dbn,
                 definitions::MoveTypeEnum move_type=definitions::NON_LOCAL)
          : MoveInfo(move_type), dbn(dbn) {}
};

}

#endif
