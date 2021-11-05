// move_collection.h --- Move collection class
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


#ifndef MOVE_COLLECTION_H
#define MOVE_COLLECTION_H

#include <boost/lexical_cast.hpp>

#include "utils/indented_stream.h"

#include "move.h"
#include "move_statistics.h"
#include "move_info.h"

namespace phaistos {

//! Collection of move objects
template <typename CHAIN_TYPE>
class MoveCollection {
private:

     //! Random number engine from which random number generators can be created.
     RandomNumberEngine *random_number_engine;

     //! Iteration number
     PHAISTOS_LONG_LONG iteration;

     //! Vector of move objects
     std::vector<Move<CHAIN_TYPE> *> moves;

     //! Vector of weights associated with the move objects
     std::vector<double> weights;

     //! Index of currently used move
     int current_move_index;

     //! Special move used for reinitialization
     Move<CHAIN_TYPE> *initializing_move;

     //! Move statistics - number of moves attempted
     PHAISTOS_LONG_LONG moves_attempted;

     //! Move statistics - number of moves accepted
     PHAISTOS_LONG_LONG moves_accepted;

public:

     //! Molecule chain object     
     CHAIN_TYPE *chain;


     //! Constructor
     //! \param chain Molecule chain
     //! \param random_number_engine Object from which random number generators can be created.
     MoveCollection(CHAIN_TYPE *chain,
                    RandomNumberEngine *random_number_engine = &random_global):
          random_number_engine(random_number_engine),
          iteration(0), initializing_move(NULL), 
          moves_attempted(0), moves_accepted(0), chain(chain) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     MoveCollection(const MoveCollection &other):
          random_number_engine(other.random_number_engine),
          iteration(other.iteration),
          moves_attempted(other.moves_attempted),
          moves_accepted(other.moves_accepted),
          chain(other.chain) {
          
          initializing_move = other.initializing_move->clone();
          
          for (unsigned int i=0; i<other.moves.size(); i++) {
               moves.push_back(other.moves[i]->clone());
               weights.push_back(other.weights[i]);
          }
     }

     //! Copy constructor - with different chain and different random number generator specified.
     //! \param other Source object from which copy is made
     //! \param chain Molecule chain object
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Optionally provide a thread number for which the copied move will be used
     MoveCollection(const MoveCollection &other, CHAIN_TYPE *chain,
                    RandomNumberEngine *random_number_engine = &random_global, 
                    int thread_index=0):
          random_number_engine(random_number_engine),
          iteration(other.iteration), 
          moves_attempted(other.moves_attempted), moves_accepted(other.moves_accepted), 
          chain(chain) {

          initializing_move = other.initializing_move->clone(random_number_engine, thread_index, chain);

          for (unsigned int i=0; i<other.moves.size(); i++) {
               moves.push_back(other.moves[i]->clone(random_number_engine, thread_index, chain));
               weights.push_back(other.weights[i]);
          }          
     }

     //! Destructor
     ~MoveCollection() {

          // Delete individual moves
          for (unsigned int i=0; i<moves.size(); i++) {
               delete moves[i];
          }

          // Delete initializing_move
          if (initializing_move) {
               delete initializing_move;
          }
     }
     
     //! Set special move that initializes chain
     //! \param move Move object
     void set_initializer(Move<CHAIN_TYPE> *move) {
          if (initializing_move)
               delete initializing_move;
	  initializing_move = move;
     }

     //! Add move object
     //! \param move Move object
     void add_move(Move<CHAIN_TYPE> *move) {

          // Duplicate names are made unique using number postfix
          bool id_match_found;
          std::string move_id = move->id;
          for (unsigned int i=0; i==0 || id_match_found; ++i) {
               id_match_found = false;
               for (unsigned int j=0; j<moves.size(); ++j) {
                    if (moves[j]->id == move_id) {
                         id_match_found = true;
                         move_id = move->id + "-" + boost::lexical_cast<std::string>(i+1);
                         break;
                    }
               }
          }
          move->id = move_id;

	  moves.push_back(move);
	  weights.push_back(move->settings.weight);
     }

     //! Apply a random move from the collection (using the associated weights)
     //! \param start_index Sequence index at which move begins
     //! \param end_index Sequence index at which move ends (not included)
     //! \param forced_move_index Optionally, a specific move index can be chosen (rather than a random one)
     bool apply(int start_index=0, int end_index=0, int forced_move_index=-1) {

          if (moves.size()==0) {
               current_move_index = -1;
               iteration++;
               moves_attempted++;
               return true;
          }

	  // Pick move
          if (forced_move_index == -1) {
               current_move_index = DiscreteSampler::sample(weights, NULL, random_number_engine);
          } else {
               current_move_index = forced_move_index;
          }
	  Move<CHAIN_TYPE> *move = moves[current_move_index];

	  // Apply move
	  bool success = move->apply(start_index, end_index);

	  iteration++;
          moves_attempted++;

          if (success) {

               // Update chaintree if present
               if (chain->chain_tree && get_last_move()->move_info->modified_angles.size() > 0) {
                    chain->chain_tree->update(*move->move_info);
               }
          }
          
	  return success;
     }

     //! Return last applied move
     //! \return Point to move object
     Move<CHAIN_TYPE> *get_last_move() {
          if (current_move_index < 0)
               return NULL;
          return moves[current_move_index];
     }

     //! Calculate log bias of last move (for MCMC acceptance probability)
     //! \return log-probability
     double get_log_bias() {
          if (current_move_index < 0)
               return 0;
          return moves[current_move_index]->get_log_bias();
     }
     
     //! Accept last move
     void accept() {
          if (current_move_index < 0)
               return;

	  moves[current_move_index]->accept();
          moves_accepted++;

          // Notify entire move collection of last move
          for (unsigned int i=0; i<moves.size(); i++) {
               if ((int)i != current_move_index) {
                    moves[i]->notify(this->get_last_move()->move_info);
               }
          }
     }

     //! Reject last move
     void reject() {
          if (!get_last_move()->move_info) {
               return;
          }

	  moves[current_move_index]->reject();

          // Update chaintree if present
          if (chain->chain_tree && get_last_move()->move_info->success && !get_last_move()->move_info->modified_angles.empty()) {
               chain->chain_tree->undo();
          }          
     }

     //! Reinitialize chain either using initialization move.
     void reinitialize() {

	  if (initializing_move) {
	       chain->check_consistency();

	       while (!initializing_move->apply()) {
		    initializing_move->reject();
	       }
	       initializing_move->accept();

               // Notify entire move collection of last move
               for (unsigned int i=0; i<moves.size(); i++) {
                    moves[i]->notify(initializing_move->move_info);
               }

	       iteration++;

               // Update chaintree if present
               if (chain->chain_tree) {
                    // chain->chain_tree->update(*move_info);
                    chain->chain_tree->reset();
               }
	       chain->check_consistency();  

               // After reinitialize, chain time_stamp=0
               chain->time_stamp = 0;
	  }
     }

     //! Retrieve information on previous move
     //! \return MoveInfo object 
     MoveInfo *get_move_info() {
          if (current_move_index < 0)
               return NULL;
          return this->get_last_move()->move_info;
     }

     //! Retrieve move statistics
     //! \return MoveStatisticsSummary object
     typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary get_statistics() {
          typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary summary;
          for (unsigned int i=0; i<moves.size(); i++) {
               typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary summary_tmp;
               summary_tmp.add_pair("Move usage", ((double)(moves[i]->statistics->moves_attempted)/
                                                     moves_attempted));
               summary_tmp += moves[i]->statistics->get_summary();

               summary.add_pair(moves[i]->id, summary_tmp);
          }
          return summary;
     }

     //! Display the settings used in the moves
     void display_settings() {

          // Create output stream where everything is indented
          boost::iostreams::filtering_ostream out;
          out.push(IndentedStreamFilter(8));
          out.push(std::cout);           

          std::cout << "Move settings:\n";
          for (unsigned int i=0; i<moves.size(); i++) {
               std::cout << "    " << moves[i]->id << "\n";
               std::cout.flush();
               out << moves[i]->display_settings();
               out.flush();
          }
     }

};

}

#endif
