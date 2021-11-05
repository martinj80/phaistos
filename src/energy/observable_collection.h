// observable_collection.h --- Collection of observables. Similar to an energy, but with observables instead of energy terms.
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

#ifndef OBSERVABLE_COLLECTION_H
#define OBSERVABLE_COLLECTION_H

namespace phaistos {

//! Output file stream that switches
//! target filename whenever update is called
class MultiFileStream: public std::ofstream {

     //! Filename template: %i will be replaced with the id sent along to update  
     std::string filename_template;

     //! Which mode to open file in
     std::ios_base::openmode mode;

public:

     //! Constructor
     //! \param filename_template String optionally containing a %i format character
     //! \param mode In which output mode to open files
     MultiFileStream(std::string filename_template, std::ios_base::openmode mode)
          : filename_template(filename_template),
            mode(mode){}

     //! Create new active filename based on filename_template and current id
     //! and open corresponding stream
     void update(std::string id) {
          this->flush();
          this->close();
          std::string log_filename = replace(filename_template, "%i", boost::lexical_cast<std::string>(id));
          this->open(log_filename.c_str(), mode);
     }     
};

//! Base class for all output styles
class OutputStyleBase {
public:

     //! Internal stream object 
     std::ostream* output_stream;

     //! Which term indices use this output style
     std::vector<int> selected_term_indices;
     
     //! Constructor
     //! \param output_stream Stream object
     OutputStyleBase(std::ostream* output_stream)
          : output_stream(output_stream) {}

     //! Destructor
     virtual ~OutputStyleBase() {
          std::ofstream *output_file_stream = dynamic_cast<std::ofstream*>(output_stream);
          if (output_file_stream) {
               output_file_stream->close();
          }
          std::stringstream *output_string_stream = dynamic_cast<std::stringstream*>(output_stream);
          if (output_file_stream || output_string_stream) {
               delete output_stream;                    
          }               
     }

     //! Register term as using this style
     //! \param index Term index
     void add_term_index(int index) {
          selected_term_indices.push_back(index);
     }

     //! Return stream object
     std::ostream* get_stream() {
          return output_stream;
     }

     //! Update stream.
     //! This function makes it possible to update the stream with a given id.
     void update(std::string id) {

          // In case of MultiFileStream, replace %i in filename with current id
          MultiFileStream* output_file_stream = dynamic_cast<MultiFileStream*>(output_stream);
          if (output_file_stream) {
               output_file_stream->update(id);
          }
     }


};

class ObservableCollectionBase {
public:

     // Constructor
     ObservableCollectionBase()
          : label_to_style_indices(OBSERVABLE_STREAM_LABEL_SIZE) {}

     //! Labels that can be used to fetch stream output
     enum ObservableStreamLabel {NONE=0, PDB_HEADER, PDB_B_FACTOR, OBSERVABLE_STREAM_LABEL_SIZE};

     //! Vector of output styles
     std::vector<OutputStyleBase*> output_styles;

     //! Map from label(enum) to style indices
     std::vector<std::vector<int> > label_to_style_indices;

     //! Which styles were active in this iteration
     std::vector<bool> style_activity;

     //! Gather all string streams with a given label
     //! \param stream_label Query label (enum)
     std::string gather_string_streams(ObservableStreamLabel stream_label) {
          std::string output = "";
          std::vector<int> &stream_indices = label_to_style_indices[(int)stream_label];
          for (unsigned int i=0; i<stream_indices.size(); ++i) {
               int stream_index = stream_indices[i];
               if (style_activity[stream_index]) {
                    std::ostream *stream = output_styles[stream_index]->get_stream();
                    std::stringstream *string_stream = dynamic_cast<std::stringstream*>(stream);
                    if (string_stream) {
                         output += string_stream->str();
                         string_stream->str("");
                    }
               }
          }
          return output;
     }
};

//! Collection of Observables. Derived from Energy.
template <typename CHAIN_TYPE>
class ObservableCollection: public Energy<CHAIN_TYPE>, public ObservableCollectionBase {

     //! Style functor defining how observables are output: Base class
     class OutputStyle: public OutputStyleBase {
     public:

          //! Constructor
          //! \param output_stream Stream object
          OutputStyle(std::ostream* output_stream)
               : OutputStyleBase(output_stream) {}

          //! Output observable (defined by derived classes)
          //! \param observable_collection The Observables to be output
          //! \param id Identification of this observation
          //! \param initialized Whether stream has been initialized.
          virtual void output(ObservableCollection<CHAIN_TYPE> &observable_collection,
                              std::string id, bool initialized=false)=0;
     };

     //! Style functor defining how observables are output: Verbose annotated with time spent.
     //! Uses Energy::get_data.
     class OutputStyleVerboseWithTime: public OutputStyle {
     public:
          //! Constructor
          //! \param output_stream Stream object
          OutputStyleVerboseWithTime(std::ostream* stream)
               : OutputStyle(stream){}

          //! Output observable (defined by derived classes)
          //! \param observable_collection The Observables to be output
          //! \param id Identification of this observation
          //! \param initialized Whether stream has been initialized.
          void output(ObservableCollection<CHAIN_TYPE> &observable_collection,
                      std::string id, bool initialized=false) {
               *this->output_stream << id << "\nObservables:\n";
               *this->output_stream << typename Energy<CHAIN_TYPE>::EnergyData(
                    observable_collection.get_data(), this->selected_term_indices) << "\n";
               this->output_stream->flush();
          }
     };

     //! Style functor defining how observables are output: Verbose.
     class OutputStyleVerbose: public OutputStyle {
     public:
          //! Constructor
          //! \param output_stream Stream object
          OutputStyleVerbose(std::ostream* stream)
               : OutputStyle(stream){}

          //! Output observable (defined by derived classes)
          //! \param observable_collection The Observables to be output
          //! \param id Identification of this observation
          //! \param initialized Whether stream has been initialized.
          void output(ObservableCollection<CHAIN_TYPE> &observable_collection,
                      std::string id, bool initialized=false) {
               *this->output_stream << "Observables:\n";
               for (unsigned int i=0; i<this->selected_term_indices.size(); ++i) {
                    *this->output_stream << "    " << 
                         observable_collection.terms[this->selected_term_indices[i]]->name << ": ";
                    *this->output_stream << boost::lexical_cast<std::string>(
                         observable_collection.current_values[this->selected_term_indices[i]]) << "\n";;
               }
               this->output_stream->flush();
          }
     };

     //! Style functor defining how observables are output: Verbose with ID.
     class OutputStyleVerboseWithId: public OutputStyle {
     public:
          //! Constructor
          //! \param output_stream Stream object
          OutputStyleVerboseWithId(std::ostream* stream)
               : OutputStyle(stream){}

          //! Output observable (defined by derived classes)
          //! \param observable_collection The Observables to be output
          //! \param id Identification of this observation
          //! \param initialized Whether stream has been initialized.
          void output(ObservableCollection<CHAIN_TYPE> &observable_collection,
                      std::string id, bool initialized=false) {
               *this->output_stream << id << "\nObservables:\n";
               for (unsigned int i=0; i<this->selected_term_indices.size(); ++i) {
                    *this->output_stream << "    " << 
                         observable_collection.terms[this->selected_term_indices[i]]->name << ": ";
                    *this->output_stream << boost::lexical_cast<std::string>(
                         observable_collection.current_values[this->selected_term_indices[i]]) << "\n";;
               }
               this->output_stream->flush();
          }
     };

     //! Style functor defining how observables are output: Compact
     class OutputStyleCompact: public OutputStyle {
     public:
          //! Whether stream has been initialized
          bool initialized;

          //! Constructor
          //! \param output_stream Stream object
          OutputStyleCompact(std::ostream* stream)
               : OutputStyle(stream),
                 initialized(false) {}

          //! Output observable (defined by derived classes)
          //! \param observable_collection The Observables to be output
          //! \param id Identification of this observation
          //! \param initialized Whether stream has been initialized.
          void output(ObservableCollection<CHAIN_TYPE> &observable_collection,
                      std::string id, bool initialized=false) {

               //! Output ID if not initialized
               if (!initialized && !this->initialized) {
                    *this->output_stream << "# ID\t";
                    for (unsigned int i=0; i<this->selected_term_indices.size(); ++i) {
                         *this->output_stream << observable_collection.terms[this->selected_term_indices[i]]->name << "\t";
                    }
                    *this->output_stream << "\n";
                    this->initialized = true;
               }
               *this->output_stream << id << "\t";
               for (unsigned int i=0; i<this->selected_term_indices.size(); ++i) {
                    std::string output_string = observable_collection.current_values[this->selected_term_indices[i]];
                    *this->output_stream << output_string << "\t";
               }
               *this->output_stream << "\n";
               this->output_stream->flush();
          }
     };

     //! Style functor defining how observables are output: Minimal
     class OutputStyleMinimal: public OutputStyle {
     public:
          //! Constructor
          //! \param output_stream Stream object
          OutputStyleMinimal(std::ostream* stream)
               : OutputStyle(stream) {}

          //! Output observable (defined by derived classes)
          //! \param observable_collection The Observables to be output
          //! \param id Identification of this observation
          //! \param initialized Whether stream has been initialized.
          void output(ObservableCollection<CHAIN_TYPE> &observable_collection,
                      std::string id, bool initialized=false) {

               for (unsigned int i=0; i<this->selected_term_indices.size(); ++i) {
                    std::string output_string = boost::lexical_cast<std::string>(observable_collection.current_values[this->selected_term_indices[i]]);
                    *this->output_stream << output_string << "\t";
               }
               *this->output_stream << "\n";
               this->output_stream->flush();
          }
     };


     //! Create output style from string specification
     //! \param output_style_string String specification of output style
     //! \param output_stream Stream used by output style.
     //! \return New OutputStyle object
     OutputStyle *get_output_style(std::string output_style_string,
                                   std::ostream* output_stream) {

          OutputStyle *output_style = NULL;
          if (output_style_string == "verbose") {
               output_style = new OutputStyleVerbose(output_stream);
          } else if (output_style_string == "verbose-with-id") {
               output_style = new OutputStyleVerboseWithId(output_stream);
          } else if (output_style_string == "verbose-with-time") {
               output_style = new OutputStyleVerboseWithTime(output_stream);
          } else if (output_style_string == "compact") {
               output_style = new OutputStyleCompact(output_stream);
          } else if (output_style_string == "minimal") {
               output_style = new OutputStyleMinimal(output_stream);
          } else {
               std::cerr << "Error (ObservableCollection): Unknown output style " << output_style_string << "\n";
               assert(false);
          }
          return output_style;
     }

     // //! Vector of output styles
     // std::vector<OutputStyleBase*> output_styles;

     //! Map from term index to style indices
     std::vector<std::vector<int> > term_index_to_style_indices;

     //! Map from style name to style index
     std::map<std::string,std::map<std::string,int> > style_name_to_style_index;

     //! Settings objects used by the observable terms
     std::vector<const ObservableBase::Settings*> settings_objects;

     // //! Map from label(enum) to style indices
     // std::vector<std::vector<int> > label_to_style_indices;

     //! Which styles require update() to be called
     std::vector<int> update_requiring_style_indices;

     //! Current observed output strings 
     std::vector<std::string> current_values;

     //! Current iteration number
     PHAISTOS_LONG_LONG current_iteration;

     // //! Which styles were active in this iteration
     // std::vector<bool> style_activity;

     //! Which terms were active in this iteration
     std::vector<bool> term_activity;

     //! Whether existing output files should be overwritten
     bool overwrite_existing_files;

public: 

     //! Constructor
     //! \param chain Molecule chain.
     //! \param overwrite_existing_files Whether existing files should be overwritten
     ObservableCollection(CHAIN_TYPE *chain, bool overwrite_existing_files=true)
          : Energy<CHAIN_TYPE>(chain),
            current_iteration(0),
            overwrite_existing_files(overwrite_existing_files) {
     }

     //! Copy constructor - using different chain.
     //! \param other Source object from which copy is made.
     //! \param chain Molecule chain.
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists.
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     ObservableCollection(const ObservableCollection &other, 
                          CHAIN_TYPE *chain, 
                          RandomNumberEngine *random_number_engine = &random_global, 
                          int thread_index=0,
                          Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : Energy<CHAIN_TYPE>(other,chain, random_number_engine, thread_index),
            current_values(other.current_values),
            current_iteration(0),
            overwrite_existing_files(other.overwrite_existing_files) {

          for (unsigned int i=0; i<this->terms.size(); ++i) {
               convert_term_to_observable(i, *other.settings_objects[i], thread_index, reference_energy_function);
          }
     }

     //! Desctructor
     ~ObservableCollection() {
          for (unsigned int i=0; i<output_styles.size(); ++i) {
               delete output_styles[i];
          }
     }

     void initialize(unsigned int index, std::string output_target, int thread_index=0) {

          // Allocate space in current values vector
          if (current_values.size() < this->terms.size())
               current_values.resize(this->terms.size(), "");

          // observables with base class IndependentObservable control their own I/O, and
          // are not registered here
          IndependentObservable *independent_observable = dynamic_cast<IndependentObservable*>(this->terms[index]);
          if (independent_observable) {
               return;
          }

          // Parse output_target string (a term can have several output styles)
          std::vector<std::string> output_target_vector;
          boost::split(output_target_vector, output_target, boost::is_any_of("+"));
          for (unsigned int i=0; i<output_target_vector.size(); ++i) {

               std::string output_target = boost::trim_copy(output_target_vector[i]);

               std::string output_style_string = "";
               if (output_target.find('#') != std::string::npos) {
                    std::vector<std::string> output_target_style;
                    boost::split(output_target_style, output_target, boost::is_any_of("#"));
                    output_target = output_target_style[0];
                    output_style_string = output_target_style[1];
               }

               // Choose style 
               OutputStyle *output_style;
               bool single_output_only = false;
               bool stream_requires_update = false;
               ObservableStreamLabel label = NONE;
               if (output_target == "stdout" || output_target == "cout") {
                    if (output_style_string == "")
                         output_style_string = "verbose";
                    output_style = get_output_style(output_style_string, &std::cout);
               } else if (output_target == "stderr" || output_target == "cerr") {
                    if (output_style_string == "")
                         output_style_string = "verbose";
                    output_style = get_output_style(output_style_string, &std::cerr);
               } else if (output_target == "pdb-header") {
                    if (output_style_string == "")
                         output_style_string = "verbose";
                    output_style = get_output_style(output_style_string, new std::stringstream());
                    label = PDB_HEADER;
               } else if (output_target == "pdb-b-factor" || output_target == "pdb-b-factors") {
                    if (output_style_string == "")
                         output_style_string = "minimal";
                    output_style = get_output_style(output_style_string, new std::stringstream());
                    single_output_only = true;
                    label = PDB_B_FACTOR;
               } else {
                    if (output_target.find("%i") != std::string::npos) {
                         if (output_style_string == "")
                              output_style_string = "verbose";
                         std::string log_filename = generate_log_filename(output_target, thread_index);
                         std::ios_base::openmode mode = std::ios_base::out;
                         if (!overwrite_existing_files)
                              mode = std::ios_base::app;
                         output_style = get_output_style(output_style_string, new MultiFileStream(log_filename.c_str(), mode));
                         stream_requires_update = true;
                    } else {
                         if (output_style_string == "")
                              output_style_string = "compact";
                         std::ios_base::openmode mode = std::ios_base::out;
                         if (!overwrite_existing_files)
                              mode = std::ios_base::app;
                         std::string log_filename = generate_log_filename(output_target, thread_index);
                         output_style = get_output_style(output_style_string, new std::ofstream(log_filename.c_str(), mode));
                    }
               }

               // Check if stream is already present
               int stream_index = -1;
               if ((stream_index = map_lookup(map_lookup(style_name_to_style_index, 
                                                         output_target), 
                                              output_style_string, -1)) != -1) {

                    if (single_output_only) {
                         std::cerr << "Error: Only single output stream allowed for " << output_style_string << ". Aborting\n";
                         assert(false);
                    }

                    output_styles[stream_index]->add_term_index(index);

                    if (index >= term_index_to_style_indices.size())
                         term_index_to_style_indices.resize(index+1);
                    term_index_to_style_indices[index].push_back(stream_index);

               } else {
                    output_styles.push_back(output_style);

                    output_styles.back()->add_term_index(index);
                    stream_index = output_styles.size()-1;
                    if (style_name_to_style_index.count(output_target) == 0)
                         style_name_to_style_index[output_target] = std::map<std::string,int>();
                    style_name_to_style_index[output_target][output_style_string] = stream_index;
                    if (index >= term_index_to_style_indices.size())
                         term_index_to_style_indices.resize(index+1);
                    term_index_to_style_indices[index].push_back(stream_index);

                    label_to_style_indices[(int)label].push_back(stream_index);

                    style_activity.push_back(false);

                    if (stream_requires_update)
                         update_requiring_style_indices.push_back(stream_index);
               }
          }
     }

     //! Check whether stream with label is present in collection
     //! \param stream_label Query label (enum)
     bool has_stream(ObservableStreamLabel stream_label) const {
          const std::vector<int> &stream_indices = label_to_style_indices[(int)stream_label];
          if (!stream_indices.empty()) {
               return true;
          } else {
               return false;
          }
     }

     //! Convert a term to an observable
     //! \param index Index where energy term is found
     //! \param settings ObservableBase Settings object
     //! \param thread_index Optional thread index
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     void convert_term_to_observable(unsigned int index, 
                                     const ObservableBase::Settings &settings, 
                                     int thread_index=0,
                                     Energy<CHAIN_TYPE> *reference_energy_function=NULL) {

          // Transform energy term to observable
          this->terms[index] = this->terms[index]->clone_to_observable(settings, reference_energy_function);

          // Save settings object
          if (index >= settings_objects.size())
               settings_objects.resize(index+1);
          settings_objects[index] = &settings;

          initialize(index, settings.output_target, thread_index);
     }


     //! Convert a term to an observable - simple version for when no settings object is available
     //! \param index Index where energy term is found
     //! \param thread_index Optional thread index
     void convert_term_to_observable(unsigned int index, int thread_index=0) {

          std::cout << "WARNING (observable_collection): Using observable without associated observable Settings object.\n";

          // Allocate space in current values vector
          if (current_values.size() < this->terms.size())
               current_values.resize(this->terms.size(), "");

          initialize(index, "stdout", thread_index);
     }



     //! Evaluate energy using all terms.
     //! A MoveInfo object can optionally be passed along, making it
     //! possible for the individual energy terms to speed up
     //! calculation by only recalculating energy contributions from
     //! modified regions of the chain.
     //! \param move_info object containing information about last move
     double evaluate(MoveInfo *move_info=NULL) {

          // Reset style activity
          this->style_activity = std::vector<bool>(style_activity.size(), false);

          // Reset term activity
          this->term_activity = std::vector<bool>(this->terms.size(), false);

          // Evaluate energy
          for (unsigned int i=0; i<this->terms.size(); ++i) {

               double time1 = get_time();

               ObservableBase *observable_base = dynamic_cast<ObservableBase*>(this->terms[i]);
               if (observable_base) {
                    const ObservableBase::Settings *observable_base_settings = settings_objects[i];
                    if (((this->current_iteration % observable_base_settings->register_interval)==0) && 
                        (this->current_iteration >= observable_base_settings->register_burnin)) {

                         term_activity[i] = true;

                         if ((this->current_iteration % observable_base_settings->output_interval)==0) {
                              this->current_values[i] = observable_base->observe(move_info, this->current_iteration);

                              IndependentObservable *independent_observable = 
                                   dynamic_cast<IndependentObservable*>(this->terms[i]);

                              // Observables with base class IndependentObservable control their own I/O
                              // and are not registered here
                              if (!independent_observable) {
                                   std::vector<int> &output_indices = term_index_to_style_indices[i];
                                   for (unsigned int j=0; j<output_indices.size(); ++j) {
                                        this->style_activity[output_indices[j]] = true;
                                   }
                              }
                         } else {
                              bool register_only=true;
                              this->current_values[i] = observable_base->observe(move_info, this->current_iteration, register_only);
                         }
                    } else {
                         this->current_values[i] = "0";
                    }
               } else {
                    this->current_values[i] = boost::lexical_cast<std::string>(this->terms[i]->evaluate_weighted(move_info));
                    std::vector<int> &output_indices = term_index_to_style_indices[i];
                    for (unsigned int j=0; j<output_indices.size(); ++j) {
                         this->style_activity[output_indices[j]] = true;
                    }
               }

               double time2 = get_time();
               this->term_times[i] += time2-time1;
               this->term_evaluations[i]++;
          }

          this->evaluations++;
          return UNINITIALIZED;
     }
     

     //! Make an observation
     //! \param id Identification string for observation
     //! \param iteration Iteration number (used to decide when to record and output)
     //! \param initialized Whether stream has been initialized
     void observe(std::string id, PHAISTOS_LONG_LONG iteration=0, bool initialized=false) {
          this->current_iteration = iteration;
          this->evaluate();

          for (unsigned int i=0; i<update_requiring_style_indices.size(); ++i) {
               if (style_activity[i]) {
                    output_styles[update_requiring_style_indices[i]]->update(id);
               }
          }

          for (unsigned int i=0; i<output_styles.size(); ++i) {
               if (style_activity[i]) {
                    OutputStyle *output_style = dynamic_cast<OutputStyle*>(output_styles[i]);
                    output_style->output(*this, id, initialized);
               }
          }

          for (unsigned int i=0; i<this->terms.size(); ++i) {
               if (term_activity[i]) {
                    ObservableBase *observable_base = dynamic_cast<ObservableBase*>(this->terms[i]);
                    observable_base->observe_finalize(this);
               }
          }
     }

     //! Get energy data. Returns a summary of energy evaluation times
     //! for each term.
     //! \param title Title to be inserted in output
     typename Energy<CHAIN_TYPE>::EnergyData get_data(std::string title="") {

          std::vector<std::string> names;
          std::vector<double> current_values_numeric;
          double energy_total = 0.0;
          for (unsigned int i=0; i<this->terms.size(); i++) {
               names.push_back(this->terms[i]->name);
               current_values_numeric.push_back(boost::lexical_cast<double>(current_values[i]));
               energy_total += current_values_numeric[i];
          }
          return typename Energy<CHAIN_TYPE>::EnergyData(title,
                                                         names, current_values_numeric, energy_total, 
                                                         this->term_times, this->term_evaluations, this->evaluations);  
     }

};

}

#endif
