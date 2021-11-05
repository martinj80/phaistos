// energy.h --- Main energy class
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


#ifndef ENERGY_H
#define ENERGY_H

#include "protein/residue.h"
#include "protein/atom.h"
#include "utils/utils.h"
#include "utils/random.h"
#include "energy_term.h"
#include "limits.h"
#include "utils/indented_stream.h"

namespace phaistos {

//! Main energy class - a collection of energy terms of which the sum
//! constitutes the total energy value.
//! Energies in Phaistos are actually -log probabilities and can
//! therefore also represent probabilistic models
template <typename CHAIN_TYPE>
class Energy {
protected:

     //! Chain molecule
     CHAIN_TYPE *chain;

     //! Vector of booleans specifying which terms are "owned" by the
     //! energy, in the sense that they should be deleted upon
     //! deconstruction
     std::vector<bool> terms_ownership;

     //! Internal table for fast lookup of terms by their name
     std::map<std::string, int> name_index_map;

     //! For timing purposes: time spent by each term
     std::vector<double> term_times;

     //! For timing purposes: number of evaluations by each term
     std::vector<PHAISTOS_LONG_LONG> term_evaluations;

     //! Total number of energy evaluations
     PHAISTOS_LONG_LONG evaluations;

public:

     //! Vector of most recently calculated values for energy terms
     std::vector<double> current_values;

     //! Backup vector of most recently calculated values for energy
     //! terms (used in case of reject)
     std::vector<double> current_values_backup;

     //! Sum of most recently calculated values for energy terms
     double sum;

     //! Backup value of sum of most recently calculated values for energy terms
     double sum_backup;

     //! Vector of energy terms
     std::vector<EnergyTerm<CHAIN_TYPE> *> terms;

     //! Representation of energy evaluation data. 
     //! Used to generate a report of the time spent on each term.
     class EnergyData {

          //! Optional output header
          std::string title;

          //! vector of term names
          std::vector<std::string> names;

          //! Internal table for fast lookup of terms by their name
          std::map<std::string, int> name_index_map;

          //! Vector of most recently calculated values for energy terms
          std::vector<std::vector<double> > values;

          //! Summed total of energy terms
          std::vector<double> value_total;

          //! Average time spent on each terms
          std::vector<std::vector<double> > term_time_avg;

          // Average time spent on all terms
          std::vector<double> time_avg;
     public:

          //! Default constructor
          EnergyData(){}

          //! Constructor.
          //! \param title optional output header
          //! \param names names of energy terms
          //! \param values most recent values of energy terms
          //! \param value_total most recent energy total
          //! \param term_times total time spent on each term
          //! \param term_evaluations total number of energy evaluations of each term
          //! \param evaluations total number of energy evaluations of any term
          EnergyData(std::string title,
                     std::vector<std::string> &names,
                     std::vector<double> &values,
                     double value_total,
                     std::vector<double> term_times,
                     std::vector<PHAISTOS_LONG_LONG> term_evaluations,
                     PHAISTOS_LONG_LONG evaluations)
               : title(title) {

               int evaluations_total = 0;
               double time_total = 0;

               for (unsigned int i=0; i<names.size(); ++i) {
                    this->names.push_back(names[i]);
                    this->values.push_back(std::vector<double>(1, values[i]));
                    this->value_total = std::vector<double>(1, value_total);
                    this->name_index_map.insert(make_pair(names[i], i));

                    // timing
                    this->term_time_avg.push_back(std::vector<double>(1, term_times[i]/((double)term_evaluations[i])));
                    time_total += term_times[i];
                    evaluations_total += term_evaluations[i];
               }

               this->time_avg.push_back(time_total/((double)evaluations));
          }

          //! Constructor. Makes it possible to gather a vector of data objects into one
          //! \param data vector of data objects
          EnergyData(const std::vector<EnergyData> &data) {

               // Assume name and weight vector of all data objects are the same
               for (unsigned int i=0; i<data[0].names.size(); ++i) {
                    this->names.push_back(data[0].names[i]);
                    // this->weights.push_back(data[0].weights[i]);
                    this->name_index_map.insert(make_pair(names[i], i));
                    this->values.push_back(std::vector<double>());
                    this->term_time_avg.push_back(std::vector<double>());
               }

               for (unsigned int i=0; i<data.size(); i++) {
                    for (unsigned int j=0; j< data[i].value_total.size(); j++) {
                         value_total.push_back(data[i].value_total[j]);
                         time_avg.push_back(data[i].time_avg[j]);
                    }

                    for (unsigned int j=0; j<data[i].values.size(); j++) {
                         for (unsigned int k=0; k<data[i].values[j].size(); k++) {
                              values[j].push_back(data[i].values[j][k]);
                              term_time_avg[j].push_back(data[i].term_time_avg[j][k]);
                         }
                    }
               }
          }

          //! Copy constructor.
          //! \param other Source object from which copy is made.
          //! \param indices Vector of indices to include
          EnergyData(const EnergyData &other, std::vector<int> indices)
               : title(other.title),
                 value_total(other.value_total),
                 time_avg(other.time_avg) {

               for (unsigned int i=0; i<indices.size(); ++i) {
                    int index = indices[i];
                    names.push_back(other.names[index]);
                    values.push_back(other.values[index]);
                    term_time_avg.push_back(other.term_time_avg[index]);
               }
               for (unsigned int i=0; i<other.names.size(); ++i) {
                    if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
                         for (unsigned int j=0; j<value_total.size(); ++j) {
                              value_total[j] -= other.values[i][j];
                              
                              // Approximate expression for time_avg
                              time_avg[j] -= other.term_time_avg[i][j];
                         }
                    }
               }
          }

          //! Return total energy
          std::vector<double> get_total() {
               return value_total;
          }

          //! Overload [] operator
          //! \param term_name name of energy term
          const std::vector<double> &operator[](const std::string &term_name) {
               return this->values[name_index_map[term_name]];
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &o, const EnergyData &ed) {
               
               std::string tot_str;
               if (ed.title != "")
                    tot_str += ed.title + ":\n";
               tot_str += "    Total: ";
               if (ed.value_total.size()==0)
                    tot_str += "0.0";
               else if (ed.value_total.size()==1)
                    tot_str += stringify(ed.value_total[0]);
               else
                    tot_str += stringify(ed.value_total);

               int max_width = tot_str.size();

               std::vector<std::string> str_vec;
               for (unsigned int i=0; i<ed.values.size(); ++i) {
                    // str_vec.push_back("\n    " + ed.names[i] + "(" + stringify(ed.weights[i]) + "): ");
                    str_vec.push_back("\n        " + ed.names[i] + ": ");
                    if (ed.values[i].size() == 1) { 
                         str_vec.back() += stringify(ed.values[i][0]);
                    } else {
                         str_vec.back() += stringify(ed.values[i]);
                    }
                    max_width = (int) fmax(max_width, str_vec.back().size());
               }

               o << tot_str;
               for (unsigned int j=0; j<((max_width-tot_str.size())); ++j) {
                    o << " ";
               }               
               o << "    \t(avg time/eval: " << ed.time_avg << " sec)\n";

               o << "    Terms:";
               for(unsigned int i=0; i<str_vec.size(); i++) {
                    o << str_vec[i];
                    for (unsigned int j=0; j<((max_width-str_vec[i].size())); ++j) {
                         o << " ";
                    }
                    o << "    \t(avg time/eval: " << ed.term_time_avg[i] << " sec)";
               }
               return o;
          }          

     };

     //! A collection of several energy data objects
     class EnergyDataCollection {
          std::vector<EnergyData> data_vector;
     public:
          //! Default constructor
          EnergyDataCollection(){}

          //! Constructor
          //! \param data single EnergyData object
          EnergyDataCollection(const EnergyData &data)
               : data_vector(1, data) {}

          //! Add a data object
          //! \param data single EnergyData object
          void add(const EnergyData &data) {
               data_vector.push_back(data);
          }

          //! Resize collection
          //! \param size New size of collection
          void resize(unsigned int size) {
               data_vector.resize(size);
          }

          //! Return size of collection
          //! \return size of collection
          unsigned int size() const {
               return data_vector.size();
          }

          //! Overload [] operator
          EnergyData &operator[](const unsigned int index) {
               return data_vector[index];
          }

          //! Overload [] operator - const
          EnergyData operator[](const unsigned int index) const {
               return data_vector[index];
          }

          //! Overload output operator
          friend std::ostream &operator<<(std::ostream &o, const EnergyDataCollection &edc) {
               for (unsigned int i=0; i<edc.size(); ++i) {
                    o << edc[i];
                    if (i<edc.size()-1) 
                         o << "\n";
               }
               return o;
          }
     };

     
     //! Constructor.
     //! \param chain Molecule chain.
     //! \param random_number_engine Object from which random number generators can be created.
     Energy(CHAIN_TYPE *chain,
            RandomNumberEngine *random_number_engine = &random_global)
          : chain(chain),
            evaluations(0) {
     }

     //! Copy constructor helper function.
     void copy(const Energy &other, CHAIN_TYPE *chain=NULL, 
               RandomNumberEngine *random_number_engine = &random_global,                
               int thread_index=0) {
          
          for (unsigned int i=0; i<other.terms.size(); i++) {

               // These terms are newly allocated - so should be deleted upon deconstruction
               bool ownership = true;

               (*this).add_term(other.terms[i]->clone(random_number_engine, thread_index, chain), ownership);
          }
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made.
     Energy(const Energy &other)
          : chain(other.chain),
            evaluations(0) {
          copy(other, other.chain);
     }

     //! Copy constructor - using different chain.
     //! \param other Source object from which copy is made.
     //! \param chain Molecule chain.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists.
     Energy(const Energy &other, CHAIN_TYPE *chain, 
            RandomNumberEngine *random_number_engine = &random_global, 
            int thread_index=0)
          : chain(chain),
            evaluations(0) {
          copy(other, chain, random_number_engine, thread_index);
     }
     
     //! Destructor.
     virtual ~Energy() {
          for (unsigned int i=0; i<terms.size(); i++) {
               if (terms_ownership[i]) {
                    delete terms[i];
               }
          }

     }

     //! Overload indexing operator (using index).
     //! \param index index of energy term in terms vector
     EnergyTerm<CHAIN_TYPE> *operator[](int index) {
          return terms.at(index);
     }

     //! Overload indexing operator.
     //! \param name name of energy term.
     EnergyTerm<CHAIN_TYPE> *operator[](std::string name) {
          return terms.at(map_lookup(name_index_map,name,-1));
     }
     
     //! Add energy term to collection. 
     //!
     //! When ownership==true, the energy term pointer will be freed
     //! upon destruction of the Energy object. This makes it possible
     //! register both energy terms that are allocated on the stack
     //! and on the heap (using new) and on the stack.
     //!
     //! \param term Pointer to energy term
     //! \param ownership Whether the pointer belongs to the collection.
     void add_term(EnergyTerm<CHAIN_TYPE> *term, bool ownership=true) {

          // Duplicate names are made unique using number postfix
          bool name_match_found;
          std::string term_name = term->name;
          for (unsigned int i=0; i==0 || name_match_found; ++i) {
               name_match_found = false;
               for (unsigned int j=0; j<terms.size(); ++j) {
                    if (terms[j]->name == term_name) {
                         name_match_found = true;
                         term_name = term->name + "-" + boost::lexical_cast<std::string>(i+1);
                         break;
                    }
               }
          }
          term->name = term_name;

          terms.push_back(term);
          terms_ownership.push_back(ownership);
          name_index_map.insert(std::pair<std::string,int>(term->name,terms.size()-1));

          current_values.push_back(-1.0);
          current_values_backup.push_back(-1.0);

          // Set up timing
          term_times.push_back(0.0);
          term_evaluations.push_back(0);
          evaluations = 0;
     }


     //! Accept energy due to last move. 
     //! Since energy terms are allowed to maintain a state (for
     //! instance a cache), it is important that the energy terms are
     //! notified about whether the last state was accepted or
     //! rejected.
     void accept() {
          for (unsigned int i=0; i<terms.size(); ++i) {
               terms[i]->accept();
          }
          for (unsigned int i=0; i<current_values_backup.size(); ++i) {
               current_values_backup[i] = current_values[i];
               sum_backup = sum;
          }
     }

     //! Reject energy due to last move.
     //! Since energy terms are allowed to maintain a state (for
     //! instance a cache), it is important that the energy terms are
     //! notified about whether the last state was accepted or
     //! rejected.
     void reject() {
          for (unsigned int i=0; i<terms.size(); ++i) {
               terms[i]->reject();
          }
          for (unsigned int i=0; i<current_values.size(); ++i) {
               current_values[i] = current_values_backup[i];
               sum = sum_backup;
          }
     }

     //! Evaluate energy using all terms.
     //! A MoveInfo object can optionally be passed along, making it
     //! possible for the individual energy terms to speed up
     //! calculation by only recalculating energy contributions from
     //! modified regions of the chain.
     //! \param move_info object containing information about last move
     double evaluate(MoveInfo *move_info=NULL) {

          // Evaluate energy
          sum = 0.0;
          for (unsigned int i=0; i<terms.size(); ++i) {
               
               double time1 = get_time();

               current_values_backup[i] = current_values[i];
               current_values[i] = terms[i]->evaluate_weighted(move_info);
               sum += current_values[i];

               double time2 = get_time();
               term_times[i] += time2-time1;
               term_evaluations[i]++;

               // Exit immediately if infinite energy has been returned
               // By putting a boolean clash detection as the first term, this can speed
               // up calculations 
               if (!std::isfinite(sum)) {
                    // In case it is a Nan, we set it to inf (guaranteed rejection)
                    current_values[i] = std::numeric_limits<double>::infinity();
                    sum = std::numeric_limits<double>::infinity();
                    break;
               }
          }
          evaluations++;
          return sum;
     }

     //! Evaluate all terms - return vector of values.
     //! A MoveInfo object can optionally be passed along, making it
     //! possible for the individual energy terms to speed up
     //! calculation by only recalculating energy contributions from
     //! modified regions of the chain.
     //! \param move_info object containing information about last move
     std::vector<double> &evaluate_terms(MoveInfo *move_info=NULL) {
          evaluate(move_info);
          return current_values;
     }


     //! Get energy data. Returns a summary of energy evaluation times
     //! for each term.
     EnergyData get_data(std::string title="") {

          std::vector<std::string> names;
          double energy_total = 0.0;
          for (unsigned int i=0; i<terms.size(); i++) {
               names.push_back(terms[i]->name);
               energy_total += current_values[i];
          }
          return EnergyData(title,
                            names, current_values, energy_total, 
                            term_times, term_evaluations, evaluations);  
     }

     //! Returns a bias associated with the energy terms.
     //! In (rare) cases where an energy term adds extra parameters to the
     //! state of the simulation, the resampling of these values might
     //! be associated with a bias. 
     double get_log_bias(MoveInfo *moveInfo = NULL) {

          double log_bias = 0.0;
          for (unsigned int i=0; i<terms.size(); ++i) {
               log_bias += terms[i]->get_log_bias(moveInfo);
          }         
          return log_bias;
     }

     //! Display the settings used in the moves.
     void display_settings() {
          // Create output stream where everything is indented
          boost::iostreams::filtering_ostream out;
          out.push(IndentedStreamFilter(8));
          out.push(std::cout);           

          std::cout << "Energy settings:\n";
          for (unsigned int i=0; i<terms.size(); i++) {
               std::cout << "    " << terms[i]->name << "\n";
               std::cout.flush();
               out << terms[i]->display_settings();
               out.flush();               
          }
     }

     //! Output last evaluated energy.
     friend std::ostream &operator<<(std::ostream &o, Energy &e) {

          std::string out = "";
          out += "    Energy terms:";
          double tot_energy = 0.0;
          for (unsigned int i=0; i<e.current_values.size(); i++) {
               tot_energy += e.current_values[i];
               out += std::string("\n        ") + e.terms[i]->name + ": \t" + stringify(e.current_values[i]);
          }
          out = std::string("Total Energy: ") + stringify(tot_energy) + "\n" + out;
          o << out;
          return o;
     }     
};

}

#endif    
