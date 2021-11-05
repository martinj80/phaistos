// move_none.h --- None move: No structural resampling
// Copyright (C) 2011 Simon Olsson, Wouter Boomsma
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


#ifndef MOVE_NONE_H
#define MOVE_NONE_H

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"
#include "move.h"

namespace phaistos {

//! Move that does nothing. Used for testing purposes and in cases
//! where extended state variables are resampled conditioned on the current
//! structural state
template <typename CHAIN_TYPE>
class MoveNone: public MoveCommon<MoveNone<CHAIN_TYPE>,
                                          CHAIN_TYPE> {

     //! For convenience, define local MoveCommon
     typedef ::phaistos::MoveCommon<MoveNone<CHAIN_TYPE>,CHAIN_TYPE> MoveCommon;

public:

     //! Use base class Settings class
     typedef typename Move<CHAIN_TYPE>::Settings Settings;

     //! Local settings object
     const Settings settings;
     
     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveNone(CHAIN_TYPE *chain, const Settings &settings = Settings(),
              RandomNumberEngine *random_number_engine = &random_global)
          : MoveCommon(chain, "none", settings, random_number_engine), 
            settings(settings) {
     }
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveNone(const MoveNone &other)
          : MoveCommon(other),
            settings(other.settings) {
     }

     //! Copy constructor. Random number engine, thread index and chain specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Which thread the copy will run in
     //! \param chain Molecule chain object
     MoveNone(const MoveNone &other, 
              RandomNumberEngine *random_number_engine, int thread_index, CHAIN_TYPE *chain)
          : MoveCommon(other, random_number_engine, chain),
            settings(other.settings) {
     }

     //! Apply move
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether move was succesful
     bool apply(int start_index=-1, int end_index=-1) {
          // Call base class apply method
          Move<CHAIN_TYPE>::apply(start_index, end_index);
          this->move_info = new MoveInfo(definitions::NON_LOCAL);

          return true;
     }

     //! Calculate the log-bias that should be included when this move is accepted/rejected
     //! \return log-probability
     double get_log_bias() {
          return 0.0;
     }


     //! Accept last move
     void accept() {

          // Call base class accept method
          Move<CHAIN_TYPE>::accept();

     }


     //! Reject last move
     void reject() {

          // Call base class accept method
          Move<CHAIN_TYPE>::reject();

     }

};

}

#endif
