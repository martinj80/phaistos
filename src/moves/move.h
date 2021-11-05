// move.h --- Move base class
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


#ifndef MOVE_H
#define MOVE_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

#include "utils/random.h"

#include "nodes/discrete.h"

namespace phaistos {

// Forward declarations
template <typename CHAIN_TYPE> class MoveCollection;
template <typename CHAIN_TYPE> class MoveStatistics;
template <typename CHAIN_TYPE, typename MOVE_TYPE> class MoveStatistics_localmove;

//! Move base class
template <typename CHAIN_TYPE>
class Move {

     //! Make MoveStatistics a friend
     friend class MoveStatistics<CHAIN_TYPE>;

public:

     //! Type of chain molecule associated with move
     typedef CHAIN_TYPE ChainType;
     
protected:

     //! Chain object
     CHAIN_TYPE *chain;

     //! Temporary backup chain used by move for acceptance and rejection
     CHAIN_TYPE *chain_backup;

     //! Random number engine from which random number generators can be created.
     RandomNumberEngine *random_number_engine;

     //! Random number generator - uniform
     boost::variate_generator<RandomNumberEngine&, 
                              boost::uniform_real<> > random_generator_uniform_01;


     //! Sample random region from region vector of (start,end) indices
     //! \return [start_index, end_index] pair
     const std::pair<int,int> &random_move_region() {

          if (settings.regions.size()==1) {
               return settings.regions[0];
          }
          // Ranges must be sampled weighed by their length 
          std::vector<double> probability_vector;
          for (unsigned int i=0; i<settings.regions.size(); ++i) {
               if ((settings.regions[i].second - settings.regions[i].first) < static_cast<int>(settings.move_length_min)) {
                    probability_vector.push_back(0);
               } else {
                    probability_vector.push_back(settings.regions[i].second - settings.regions[i].first);
               }
          }
          return settings.regions[DiscreteSampler::sample(probability_vector,NULL,random_number_engine)];
     } 


     //! Get random move length
     //! \return length of move
     int random_move_length() {

          boost::variate_generator<RandomNumberEngine&,
                                   boost::uniform_int<> >
               random_generator_uniform_int(*random_number_engine,
                                            boost::uniform_int<>(settings.move_length_min,
                                                                 settings.move_length_max));

          unsigned int length = random_generator_uniform_int();

          if (length > (unsigned int)chain->size()) {
               std::cerr << "Error: Attempting move with length " << length << " on chain of length " << chain->size() << "\n";
               exit(1);
          }

          return length;
     }


     //! Get random move length
     //! \param region Region in which move is allowed
     //! \return length of move
     int random_move_length(const std::pair<int,int> &region) {

          unsigned int region_length = region.second - region.first;
          if ((region_length) < settings.move_length_min) {
               std::cerr << "Error: Move - Attempting move with minimum length " 
                         << settings.move_length_min << " in region of length " << region_length << "\n";
          }
          boost::variate_generator<RandomNumberEngine&, 
                                   boost::uniform_int<> > 
               random_generator_uniform_int(*random_number_engine, 
                                            boost::uniform_int<>(settings.move_length_min, 
                                                                 settings.move_length_max));          
	  unsigned int length = random_generator_uniform_int();

          if (length > (unsigned int)chain->size()) {
               std::cerr << "Error: Attempting move with length " << length << " on chain of length " << chain->size() << "\n";
               exit(1);
          }

          return std::min(length, region_length); 
     }


     //! Get random move range given move length and a range minimum and maximum
     //! \param length Length of move
     //! \param move_range_min Minimum index of move
     //! \param move_range_max Maximum index of move
     //! \param start_index Pointer to variable where start index is written
     //! \param end_index Pointer to variable where end index is written
     void random_move_range(int length, int move_range_min, int move_range_max, int *start_index, int *end_index) {

          boost::variate_generator<RandomNumberEngine&, boost::uniform_int<> >
               random_generator_uniform_int(*random_number_engine, boost::uniform_int<>(move_range_min, move_range_max - length));
          *start_index = random_generator_uniform_int();
          *end_index = *start_index + length;
     }


     //! Get random move range given a region a length and a range of
     //! indices in which the move is allowed to operatate.  This
     //! method will return a range compatible with the given region.
     //! \param region A region of the chain typically specified by the user
     //! \param length Length of the move
     //! \param move_range_min The minimum of the range of the chain in which the move (by design) can operate.
     //! \param move_range_max The maximum of the range of the chain in which the move (by design) can operate.
     //! \param start_index Pointer to variable where start index is written
     //! \param end_index Pointer to variable where end index is written
     void random_move_range(const std::pair<int,int> &region, int length, 
                          int move_range_min, int move_range_max, int *start_index, int *end_index) {

          // Make sure that the move range is within region
          move_range_min = std::max(region.first, move_range_min);
          move_range_max = std::min(region.second, move_range_max);

          // Choose interval in range [move_range_min, move_range_max] of length 'length'.
          // This corresponds to the interval [-(length-1), chainsize+(length-1)]
          // for move_range_min=-(length-1) and move_range_max=chainsize+length-1
          boost::variate_generator<RandomNumberEngine&, 
                                   boost::uniform_int<> > 
               random_generator_uniform_int(*random_number_engine, 
                                             boost::uniform_int<>(move_range_min, move_range_max-length));          
          *start_index = random_generator_uniform_int();
          *end_index = *start_index+length;	  
     }


     //! Get random move range - automatically determines region.
     //! \param move_range_min The minimum of the range of the chain in which the move (by design) can operate.
     //! \param move_range_max The maximum of the range of the chain in which the move (by design) can operate.
     //! \param start_index Pointer to variable where start index is written
     //! \param end_index Pointer to variable where end index is written
     void random_move_range(int move_range_min, int move_range_max, int *start_index, int *end_index) {
          const std::pair<int,int> &region = random_move_region();
          int length = random_move_length(region);
          return random_move_range(region, length, move_range_min, move_range_max, start_index, end_index);
     }

     //! Get random range - automatically determines region and chooses among all regions positions
     //! \param start_index Pointer to variable where start index is written
     //! \param end_index Pointer to variable where end index is written
     void random_move_range(int *start_index, int *end_index) {
          const std::pair<int,int> &region = random_move_region();
          int length = random_move_length(region);
          return random_move_range(region, length, 0, chain->size(), start_index, end_index);
     }

public:

     //! Local Settings class - common to all moves.
     const class Settings: public ::Settings {
     public:

          //! Minimum move length
          unsigned int move_length_min;

          //! Maximum move length
          unsigned int move_length_max;

          //! Region of chain in which move will be applied (start included, end excluded)
          std::vector<std::pair<int,int> > regions;

          //! Weight of move
          double weight;

          //! Level of debug information
          int debug;

          //! Constructor
          Settings(int move_length_min=1, int move_length_max=1,
                   const std::vector<std::pair<int,int> > &regions = 
                   (std::vector<std::pair<int,int> >(1,std::make_pair(std::numeric_limits<int>::min(),
                                                                      std::numeric_limits<int>::max()))),
                   double weight=1.0,
                   int debug=0)
               : move_length_min(move_length_min), 
                 move_length_max(move_length_max),
                 regions(regions),
                 weight(weight), 
                 debug(debug) {}


          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "move-length-min:" << settings.move_length_min << "\n";
               o << "move-length-max:" << settings.move_length_max << "\n";
               o << "regions:" << replace(replace(stringify(settings.regions), "(","["), ")","]") << "\n";
               o << "weight:" << settings.weight << "\n";
               o << "debug:" << settings.debug << "\n";
               return o;
          }
     } settings;    //!< Local settings object


     //! Display the setting object as a string
     virtual std::string display_settings() {
          return stringify(settings);
     }

     //! Information about last move
     MoveInfo *move_info;

     //! Move statistics
     MoveStatistics<CHAIN_TYPE> *statistics;

     //! Move identifier
     std::string id;            // Name of move


     //! Constructor
     //! \param chain Molecule chain
     //! \param name Name of move
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     Move(CHAIN_TYPE *chain, 
          std::string name,
          const Settings &settings,
          RandomNumberEngine *random_number_engine = &random_global);

     //! Copy constructor
     //! \param other Source object from which copy is made
     Move(const Move &other);

     //! Copy constructor. Different random_number_generator specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param chain Optionally specify different chain object
     Move(const Move &other, RandomNumberEngine *random_number_engine, CHAIN_TYPE *chain=NULL);

     //! Clone: Corresponds to a virtual copy constructor
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Optionally provide a thread number for which the copied move will be used
     //! \param chain Optionally specify different chain object
     virtual Move *clone(RandomNumberEngine *random_number_engine=NULL, int thread_index=0, CHAIN_TYPE *chain=NULL)=0;

     //! Destructor
     virtual ~Move();

     //! Apply move
     //! The active region is from start_index (included) to end_index (not included).
     //! \param start_index Start index in sequence
     //! \param end_index End index in sequence
     //! \return Boolean indicating whether move was succesful
     virtual bool apply(int start_index=-1, int end_index=-1);
     
     //! Calculate log-bias of last move (for MCMC acceptance probability)
     //! \return log probability
     virtual double get_log_bias()=0;

     //! Accept last move
     virtual void accept();

     //! Reject last move
     virtual void reject();

     //! Defines how move reacts to the execution of other moves
     virtual void notify(MoveInfo *move_info);
};

}

#include "move_statistics.h"

namespace phaistos {

// Constructor
template <typename CHAIN_TYPE>
Move<CHAIN_TYPE>::Move(CHAIN_TYPE *chain, 
                       std::string name,
                       const Settings &settings,
                       RandomNumberEngine *random_number_engine)
     : chain(chain),
       random_number_engine(random_number_engine),
       random_generator_uniform_01(*random_number_engine, boost::uniform_real<>(0,1)),
       settings(settings),
       move_info(NULL) {

     id = name;
          
     this->chain_backup = NULL;

     this->statistics = new MoveStatistics<CHAIN_TYPE>(this);
}



// Copy constructor
template <typename CHAIN_TYPE>
Move<CHAIN_TYPE>::Move(const Move &other)
     : chain(other.chain),
       random_number_engine(other.random_number_engine),
       random_generator_uniform_01(*(other.random_number_engine), boost::uniform_real<>(0,1)),
       settings(other.settings),
       move_info(NULL),
       id(other.id) {

     this->chain_backup = NULL;

     // Create new moveStatistics
     this->statistics = new MoveStatistics<CHAIN_TYPE>(this);
}

// Copy constructor - different random_number_engine
template <typename CHAIN_TYPE>
Move<CHAIN_TYPE>::Move(const Move &other, RandomNumberEngine *random_number_engine, CHAIN_TYPE *chain)
     : random_number_engine(random_number_engine),
       random_generator_uniform_01(*random_number_engine, boost::uniform_real<>(0,1)),
       settings(other.settings),
       move_info(NULL),
       id(other.id) {

     if (chain)
          this->chain = chain;
     else
          this->chain = other.chain;      

     this->chain_backup = NULL;

     if (random_number_engine==NULL) {
          this->random_number_engine = other.random_number_engine;
     } else {
          this->random_number_engine = random_number_engine;
     }

     // Create new moveStatistics
     this->statistics = new MoveStatistics<CHAIN_TYPE>(this);
}

// Destructor
template <typename CHAIN_TYPE>
Move<CHAIN_TYPE>::~Move(){
     if (chain_backup)
          delete chain_backup;
     if (move_info) {
          delete move_info;
          move_info = NULL;
     }
     delete statistics;
}     

// Apply move
template <typename CHAIN_TYPE>
bool Move<CHAIN_TYPE>::apply(int start_index, int end_index) {
     statistics->update_on_apply();

     // Clear modified angles and positions
     if (move_info) {
          delete move_info;
          move_info = NULL;
     }

     // Increase time stamp
     chain->time_stamp+=1;

     return true;
}

// Accept last move
template <typename CHAIN_TYPE>
void Move<CHAIN_TYPE>::accept() {
     statistics->update_on_accept();
}

// Reject last move
template <typename CHAIN_TYPE>
void Move<CHAIN_TYPE>::reject() {
     statistics->update_on_reject();

     // Decrease time stamp
     chain->time_stamp-=1;
}

// Defines how move reacts to the execution of other moves
// Default: do nothing
template <typename CHAIN_TYPE>
void Move<CHAIN_TYPE>::notify(MoveInfo *move_info) {
}



//! Layer between user-defined moves and the Move base class.
//! This class takes the derived class as a template argument, and can 
//! therefore define common methods for all energy terms that would normally
//! have to be defined explicitly in each term.
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class MoveCommon: public Move<CHAIN_TYPE> {
protected: 
     
     //! Local settings definition.
     typedef typename Move<CHAIN_TYPE>::Settings Settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param name Name of move
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     MoveCommon(CHAIN_TYPE *chain, std::string name, const Settings &settings,
                RandomNumberEngine *random_number_engine)
          : Move<CHAIN_TYPE>(chain, name, settings, random_number_engine) {}

     //! Copy constructor. Different random_number_generator specified.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param chain Optionally specify different chain object
     MoveCommon(const MoveCommon &other, 
                RandomNumberEngine *random_number_engine, 
                CHAIN_TYPE *chain)
          : Move<CHAIN_TYPE>(other, random_number_engine, chain){}
     
     //! Clone: Corresponds to a virtual copy constructor
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Optionally provide a thread number for which the copied move will be used
     //! \param chain Optionally specify different chain object
     Move<CHAIN_TYPE> *clone(RandomNumberEngine *random_number_engine=NULL, int thread_index=0, CHAIN_TYPE *chain=NULL) {
          return new DERIVED_CLASS(dynamic_cast<DERIVED_CLASS const&>(*this), random_number_engine, thread_index, chain);
     }

     //! Display the setting object as a string
     virtual std::string display_settings() {
          return stringify(((DERIVED_CLASS*)(this))->settings);
     }
};

}

#endif
