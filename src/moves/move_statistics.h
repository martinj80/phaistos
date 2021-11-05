// move_statistics.h --- Move statistics class
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

#ifndef MOVE_STATISTICS_H
#define MOVE_STATISTICS_H

#include "protein/iterators/dof_iterator.h"

namespace phaistos {

//! Class for standard move statistics
template <typename CHAIN_TYPE>
class MoveStatistics {
protected:

     //! Pointer to move class
     Move<CHAIN_TYPE> *move;

public:

     //! Statistics summary class
     class MoveStatisticsSummary {
     protected:
          //! keys (names)
          std::vector<std::string> keys;

          //! map between strings and summary objects
          std::map<std::string, MoveStatisticsSummary> map;

          //! statistic values
          std::vector<double> values;

     public:

          //! Constructor (default)
          MoveStatisticsSummary() {}

          //! Constructor (with double value)
          MoveStatisticsSummary(double value)
               : values(std::vector<double>(1,value)) {}

          //! Constructor (with vector of doubles)
          MoveStatisticsSummary(std::vector<double> values)
               : values(values) {}

          //! Constructor making it possible to gather a vector of summaries into one
          MoveStatisticsSummary(std::vector<MoveStatisticsSummary> &summaries) {

               // Assume key vector of all summaries are the same
               for (unsigned int i=0; i<summaries[0].keys.size(); i++) {

                    std::string key = summaries[0].keys[i];

                    std::vector<MoveStatisticsSummary> sub_summaries;
                    for (unsigned int j=0; j<summaries.size(); j++) {
                         sub_summaries.push_back(summaries[j][key]);
                    }
                    add_pair(key, sub_summaries);
               }
               for (unsigned int j=0; j<summaries.size(); j++) {
                    for (unsigned int k=0; k<summaries[j].values.size(); k++) {
                         values.push_back(summaries[j].values[k]);
                    }
               }
          }

          //! Add a (key, summary) pair
          //! \param key Name
          //! \param move_statistics_summary Statistics summary object
          void add_pair(const std::string key, const MoveStatisticsSummary &move_statistics_summary) {
               map.insert(std::make_pair(key, move_statistics_summary));
               keys.push_back(key);
          }

          //! Merge two summary objects
          //! \param other Source object from which copy is made
          //! \return Merged summary object
          MoveStatisticsSummary &operator+=(const MoveStatisticsSummary &other) {
               for (unsigned int i=0; i<other.keys.size(); i++) {
                    const std::string key = other.keys[i];
                    const MoveStatisticsSummary value = other.map.find(key)->second;
                    add_pair(key, value);
               }
               for (unsigned int i=0; i<other.values.size(); i++) {
                    this->values.push_back(other.values[i]);
               }
               return *this;
          }

          //! Overload [] operator - lookup up summary by key
          MoveStatisticsSummary &operator[](std::string key) {
               return map[key];
          }

          //! Output summary data (recursive)
          friend std::string output(const MoveStatisticsSummary &mss, int indentation=0) {

               std::string str = "";
               for (unsigned int i=0; i<mss.keys.size(); i++) {
                    typename std::map<std::string, MoveStatisticsSummary>::const_iterator iter = mss.map.find(mss.keys[i]);
                    if( iter != mss.map.end() ) {
                         if (!(indentation==0 && i==0))
                              str += "\n";
                         for (int j=0; j<indentation; j++) {
                              str +="    ";
                         }
                         str += iter->first + ": ";
                         str += output(iter->second, indentation+1);
                    }
               }
               if (mss.values.size() ==1 ) {
                    str += stringify(mss.values[0]);
               } else if (mss.values.size() > 1 ) {
                    str += stringify(mss.values);
               }
               return str;
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &o, const MoveStatisticsSummary &mss) {
               o << output(mss);
               return o;
          }          
     };

     //! Number of moves attempted
     PHAISTOS_LONG_LONG moves_attempted;
     
     //! Number of moves failed
     PHAISTOS_LONG_LONG moves_failed;

     //! Number of moves accepted
     PHAISTOS_LONG_LONG moves_accepted;

     //! Number of moves rejected
     PHAISTOS_LONG_LONG moves_rejected;


     //! Constructor
     //! \param move Move object
     MoveStatistics(Move<CHAIN_TYPE> *move):
          move(move),
          moves_attempted(0),
          moves_failed(0),
          moves_accepted(0) {}

     //! Destructor
     virtual ~MoveStatistics() {}

     //! Action when accepting a move
     virtual void update_on_accept() {
          moves_attempted++;
          moves_accepted++;
     }

     //! Action when rejecting a move
     virtual void update_on_reject() {
          moves_attempted++;
          if (!this->move->move_info->success)
               moves_failed++;
     }

     //! Action when applying a move     
     virtual void update_on_apply() {}


     //! Retrieve statistics summary
     //! \return Statistics Summary object
     virtual MoveStatisticsSummary get_summary() const {
          MoveStatisticsSummary summary = MoveStatisticsSummary();
          summary.add_pair("Acceptance rate", moves_accepted/double(moves_attempted-moves_failed));
          summary.add_pair("Failure rate", moves_failed/double(moves_attempted));
          return summary;
     }

     //! Overload output operator
     friend std::ostream &operator<<(std::ostream &o, const MoveStatistics &ms) {
          o << ms.get_summary();
          return o;
     }          
};


//! Statistics class for local move statistics
//! This class can be used in certain local move classes to get additional statistics
template <typename CHAIN_TYPE, typename MOVE_TYPE>
class MoveStatisticsLocalMove: public MoveStatistics<CHAIN_TYPE> {
protected:

     // Move object
     MOVE_TYPE *move;     
     
public:

     //! Number of attempted end moves
     PHAISTOS_LONG_LONG end_moves_attempted;

     //! Number of failed end moves
     PHAISTOS_LONG_LONG end_moves_failed;

     //! Number of accepted end moves
     PHAISTOS_LONG_LONG end_moves_accepted;

     //! Stepsize of endmoves
     double end_moves_stepsize;

     //! Stepsize of internal moves
     double internal_moves_stepsize;

     //! Constructor
     //! \param move Move object
     MoveStatisticsLocalMove(MOVE_TYPE *move):
          MoveStatistics<CHAIN_TYPE>(move),
          end_moves_attempted(0),
          end_moves_failed(0),
          end_moves_accepted(0),
          end_moves_stepsize(0.0),
          internal_moves_stepsize(0.0) {

          // Overwrite base class move with derived move
          this->move = move;
     }

     //! Action when accepting a move
     virtual void update_on_accept() {
          MoveStatistics<CHAIN_TYPE>::update_on_accept(); // base class call
    
          double mean = 0;
          int dihedral_counter = 0;
          for (DofIterator<CHAIN_TYPE> it1=*move->begin, it2=*move->backup_begin;
               (it1!= *move->end && it2 != *move->backup_end); ++it1,++it2) {
               if(it1.get_dof_type() == definitions::DIHEDRAL){
                    if(it1.get_atom()->atom_type != definitions::N) {
                         dihedral_counter++;
                         double diff_1  = std::fabs(*it1-*it2);
                         double diff_2 = 2*M_PI-std::fabs(*it1-*it2);
                         assert(diff_2>=0);
                         if (diff_1<diff_2)
                              mean += diff_1;
                         else
                              mean += diff_2;
                    }
               }
          }
          if (dihedral_counter > 0)
               mean /= (double)dihedral_counter;
          else
               mean = 0.0;

          if (this->move->move_region!=definitions::INTERNAL) {
               end_moves_attempted++;
               end_moves_accepted++;
               end_moves_stepsize += mean;
          } else {
               internal_moves_stepsize += mean;
          }
     }

     //! Action when rejecting a move
     virtual void update_on_reject() {
          MoveStatistics<CHAIN_TYPE>::update_on_reject(); // base class call
          
          if (this->move->move_region!=definitions::INTERNAL) {
               end_moves_attempted++;
               end_moves_stepsize += 0.0;
               if (!this->move->move_info->success)
                    end_moves_failed++;
          } else {
               internal_moves_stepsize += 0.0;
          }
     }


     //! Action when applying a move
     virtual void update_on_apply() {
          MoveStatistics<CHAIN_TYPE>::update_on_apply(); // base class call
     }

     //! Retrieve statistics summary
     //! \return Statistics Summary object
     typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary get_summary() {
          typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary summary = MoveStatistics<CHAIN_TYPE>::get_summary();
          
          typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary summary_internalmove;
          summary_internalmove.add_pair("Move usage", ((this->moves_attempted-end_moves_attempted)/
                                                         double(this->moves_attempted)));
          summary_internalmove.add_pair("Acceptance rate", ((this->moves_accepted-end_moves_accepted)/
                                                   double(this->moves_attempted-end_moves_attempted)));
          summary_internalmove.add_pair("Failure rate", ((this->moves_failed-end_moves_failed)/
                                                double(this->moves_attempted-end_moves_attempted)));
          summary_internalmove.add_pair("Stepsize (degr.)", ((internal_moves_stepsize/
                                                     double(this->moves_attempted-end_moves_attempted))/M_PI*180));

          typename MoveStatistics<CHAIN_TYPE>::MoveStatisticsSummary summary_endmove;
          summary_endmove.add_pair("Move usage", ((end_moves_attempted)/(double(this->moves_attempted))));
          summary_endmove.add_pair("Acceptance rate", (end_moves_accepted/
                                                             double(end_moves_attempted)));
          summary_endmove.add_pair("Failure rate", (end_moves_failed/
                                                          double(end_moves_attempted)));
          summary_endmove.add_pair("Stepsize (degr.)", ((end_moves_stepsize/
                                                               double(end_moves_attempted))/M_PI*180));

          summary.add_pair("Internal move", summary_internalmove);
          summary.add_pair("End move", summary_endmove);
          return summary;
     }

};

}

#endif
