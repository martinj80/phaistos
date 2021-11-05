// contact_map.h
// Copyright (C) 2010 Mikael Borg
// 
// This file is part of Phaistos

// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#ifndef CONTACTMAP
#define CONTACTMAP

#include <iostream>
#include <vector>
#include <cmath>
#include <string.h>

#include "protein/iterators/atom_iterator.h"
#include "energy/energy_term.h"
#include "protein/definitions.h"

namespace phaistos {

//! Base class for energies using contact maps
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class TermContactMapBase: public EnergyTermCommon<DERIVED_CLASS, CHAIN_TYPE> {

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS,CHAIN_TYPE> EnergyTermCommon;

public:

     //! Contact.
     //! A contact between two atoms is specified by the residue numbers and
     //! AtomEnums, a contact distance, a width of the potential well and a
     //! weight.
     //!
     //! The energy of a contact is calculated as either
     //! \f[ -e^{-\left(\frac{r-r_0}{width}\right)^2}*weight \f]
     //! or
     //! \f[ \left(\frac{r-r_0}{width}\right)^2*weight \f]
     //! depending on the argument dist_squared. The former is the default.
     //!
     class Contact {
     public:

          //! Index in chain of residue1
          int residue_index1;
          
          //! Atom type of atom1
          definitions::AtomEnum atom_type1;

          //! Index in chain of residue2
          int residue_index2;

          //! Atom type of atom1
          definitions::AtomEnum atom_type2;

          //! Distance between atom1 and atom2
          double distance;

          //! Weight of this contact
          double weight;

          //! Default constructor
          Contact(){}
     
          //! Constructor
          //! \param residue_index1 index of residue 1
          //! \param atom_type1 definitions::AtomEnum of atom in residue 1
          //! \param residue_index2 index of residue 2
          //! \param atom_type2 definitions::AtomEnum of atom in residue 2
          //! \param distance ideal distance between a1 and a2
          //! \param weight weight for this contact
          Contact(int residue_index1, definitions::AtomEnum atom_type1, 
                  int residue_index2, definitions::AtomEnum atom_type2, 
                  double distance, double weight=1.0)
               : residue_index1(residue_index1),
                 atom_type1(atom_type1),
                 residue_index2(residue_index2),
                 atom_type2(atom_type2),
                 distance(distance),
                 weight(weight) {
          }

          //! Output as string - verbose version
          std::string to_string_verbose() const {
               std::string output = "";
               output += (boost::lexical_cast<std::string>(residue_index1) + " " +  
                          boost::lexical_cast<std::string>(atom_type1) + " " + 
                          boost::lexical_cast<std::string>(residue_index2) + " " + 
                          boost::lexical_cast<std::string>(atom_type2) + " " + 
                          boost::lexical_cast<std::string>(distance) + " " + 
                          boost::lexical_cast<std::string>(weight));
               return output;
          }

          //! Overload << operator for Contact (compact output)
          friend std::ostream &operator<<(std::ostream &o, const Contact &c) {
               o << "("
                 << c.residue_index1 << ","
                 << c.atom_type1 << ","
                 << c.residue_index2 << ","
                 << c.atom_type2 << ","
                 << c.distance << ","
                 << c.weight << ")";
               return o;
          }

          //! Overload >> operator for Contact pointer (compact input)
          friend std::istream &operator>>(std::istream &input, Contact &c) {
               input.ignore(1);
               input >> c.residue_index1;
               input.ignore(1);
               input >> c.atom_type1;
               input.ignore(1);
               input >> c.residue_index2;
               input.ignore(1);
               input >> c.atom_type2;
               input.ignore(1);
               input >> c.distance;
               if (input.peek() != ')') {
                    input.ignore(1);
                    input >> c.weight;
               }
               input.ignore(1);
               return input;
          }
     };

protected:

     //! vector of Contact objects
     std::vector<Contact> contact_map;

     //! normalizing factor (sum of weights)
     double normalization_factor;
     
public:

     //! Local settings class.
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! File containing contacts, two supported file formats: 
          //! 1: One line: [(residue_index1,atom_type1,residue_index2,atom_type2,distance,weight),(...),...]
          //!    Example: [(0,CA,14,CA,8.35098,1),(0,CA,15,CA,5.72554,1),...]
          //! 2: One contact per line:
          //!     residueno1 atomname1 residueno2 atomname2 idealdistance [width [weight]].          
          //!    Example:
          //!    5 CA 1 CA 13.029
          //!    6 CA 2 CA 13.168 2.0 1.5
          std::string contact_map_filename;

          //! String containing contacts.
          std::string contact_map_string;

          //! Name of PDB file - to initialize contact map from structure
          std::string pdb_filename;
          
          //! maximum contact distance if constructing contactmap from a chain
          double cutoff;

          //! Minimum distance between residues in contacts
          int minimum_residue_distance;

          //! Which atoms to include when constructing contactmap from chain
          definitions::IterateEnum iteration_type;

          //! Constructor
          Settings(std::string contact_map_filename="",
                   std::string contact_map_string="",
                   std::string pdb_filename="",
                   double cutoff = std::numeric_limits<double>::infinity(), 
                   int minimum_residue_distance = 7,
                   definitions::IterateEnum iteration_type = definitions::ALL)
               : contact_map_filename(contact_map_filename),
                 contact_map_string(contact_map_string),
                 pdb_filename(pdb_filename),
                 cutoff(cutoff),
                 minimum_residue_distance(minimum_residue_distance),
                 iteration_type(iteration_type) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "contact-map-file: " << settings.contact_map_filename << "\n";
               o << "contact-map-string: " << settings.contact_map_string << "\n";
               o << "pdb-file: " << settings.pdb_filename << "\n";
               o << "cutoff: " << settings.cutoff << "\n";
               o << "minimum-residue-distance: " << settings.minimum_residue_distance << "\n";
               o << "iteration-type: " << settings.iteration_type << "\n";
               o << (typename EnergyTerm<CHAIN_TYPE>::Settings &)settings;
               return o;
          }                    
     } settings;  //!< Local settings object 


     //! contactmap constructor
     //! \param chain chain for contact map 
     //! \param name Name of energy term
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermContactMapBase(CHAIN_TYPE *chain, std::string name,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine),
            settings(settings) {
          
          if(settings.contact_map_filename.size()>0) {
               std::ifstream input_stream(settings.contact_map_filename.c_str());
               if (!input_stream.is_open()) {          
                    std::cerr << "# Error: Cannot open contact map file " << settings.contact_map_filename << " .\n";
                    exit(EXIT_FAILURE);
               }
               read_contact_map(chain, input_stream, settings);
          } else if(settings.contact_map_string.size()>0) {
               std::stringstream input_stream(settings.contact_map_string.c_str());
               read_contact_map(chain, input_stream, settings);
          } else if (settings.pdb_filename.size()>0) {
               if(settings.debug)
                    std::cout << "# Constructing contact map from " << settings.pdb_filename << "\n";
               CHAIN_TYPE contact_map_chain(settings.pdb_filename);

               construct_contact_map(chain, &contact_map_chain);
          }

          if (settings.debug) {
               std::cout << "Initial ContactMap: \n" << contact_map << "\n";
          }
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermContactMapBase(const TermContactMapBase &other, 
                        RandomNumberEngine *random_number_engine,
                        int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            contact_map(other.contact_map),
            settings(other.settings) {}


     //! Read contact map from file
     //! \param chain Chain object for contact map 
     //! \param input_stream Stream from which data is read.
     //! \param settings settings object
     void read_contact_map(CHAIN_TYPE *chain, std::istream &input_stream, const Settings &settings=Settings()) {

          while (input_stream.good()) {
          
               std::string line;
               std::getline(input_stream, line);

               boost::trim(line);

               if (line.size()==0 || line[0] == '#') {
                    continue;
               }
                
               std::vector<std::string> split_line;
               boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);
               
               // Parse single-line entry (compact format)
               if (split_line.size() == 1) {
                    std::istringstream buffer(line);
                    buffer >> contact_map;
                    std::cout << "1: " << contact_map << "\n";
                    return;
                    // contact_map = boost::lexical_cast<Contact>(line);
               } else {
                    int residue_index1 = boost::lexical_cast<int>(split_line[0]);
                    definitions::AtomEnum atom_type1 = boost::lexical_cast<definitions::AtomEnum>(split_line[1]);
                    int residue_index2 = boost::lexical_cast<int>(split_line[2]);
                    definitions::AtomEnum atom_type2 = boost::lexical_cast<definitions::AtomEnum>(split_line[3]);
                    double distance = boost::lexical_cast<double>(split_line[4]);;
                    double weight = 1.0;
                    if (split_line.size() == 6) {
                         weight = boost::lexical_cast<double>(split_line[5]);;
                    } 
                    Contact contact(residue_index1,
                                    atom_type1,
                                    residue_index2,
                                    atom_type2,
                                    distance,
                                    weight);
                    contact_map.push_back(contact);
               }
          }          
     }

     //! Construct a contact map from chain.
     //! Wrapper for templated version below.
     void construct_contact_map(CHAIN_TYPE *chain, CHAIN_TYPE *cm_chain) {
          switch (settings.iteration_type) {
          case definitions::BACKBONE:
               construct_contact_map<definitions::BACKBONE>(chain, cm_chain);
               break;
          case definitions::CA_ONLY:
               construct_contact_map<definitions::CA_ONLY>(chain, cm_chain);
               break;
          case definitions::SC_ONLY:
               construct_contact_map<definitions::SC_ONLY>(chain, cm_chain);
               break;
          case definitions::ALL:
               construct_contact_map<definitions::ALL>(chain, cm_chain);
               break;
          default:
               std::cerr << "ContactMap: Unsupported iteration type: " << settings.iteration_type << "\n";
               assert(false);
          }
     }


     //! Construct contact map based on chain structure
     //! \param chain chain for contact map 
     //! \param cm_chain the structure that defines the contact map
     template <definitions::IterateEnum iteration_type>
     void construct_contact_map(CHAIN_TYPE *chain,  CHAIN_TYPE *cm_chain) {

          contact_map.clear();
          for (ResidueIterator<CHAIN_TYPE> res_it1(*cm_chain); !res_it1.end(); ++res_it1) {
               for (ResidueIterator<CHAIN_TYPE> res_it2(res_it1+settings.minimum_residue_distance); !res_it2.end(); ++res_it2) {
                    for (AtomIterator<CHAIN_TYPE, iteration_type> it1(*res_it1); !it1.end(); ++it1) {
                         for (AtomIterator<CHAIN_TYPE, iteration_type> it2(*res_it2); !it2.end(); ++it2) {
                              double distance = (it1->position - it2->position).norm();
                              if (distance < settings.cutoff) {
                                   contact_map.push_back(Contact(res_it1->index, it1->atom_type,
                                                                 res_it2->index, it2->atom_type,
                                                                 distance));
                              }                              
                         }
                    }
               }
          }

          normalization_factor = contact_map.size();

          if(settings.debug)
               std::cout << "Created contact map w/ " << contact_map.size() << " contacts.\n";
     }


     //! Return string representation of contact map
     static std::string contact_map_to_string_verbose(const std::vector<Contact> &contact_map) {
          std::string output = "";
          for(uint i = 0; i < contact_map.size(); i++)
               output += contact_map[i].to_string_verbose() + "\n";
          return output;
     }

     //! Return string representation of contact map - compact version
     static std::string contact_map_to_string(const std::vector<Contact> &contact_map) {
          std::string output = "";
          output += boost::lexical_cast<std::string>(contact_map);
          // for(uint i = 0; i < contact_map.size(); i++) {
          //      if (i>0)
          //           output += ",";
          //      output += contact_map[i].to_string_compact();
          // }
          return output;
     }

};

//! Contact map energy term.
//! A contact map consists of a vector of Contact objects
//! the energy is the sum of the contact energies.
//! 
//! The functional form of the contact energies is either exp(-((d-d0)/width)^2)*w (default) 
//! or ((d-d0)/width)^2*w and it is normalized so that
//! when all contacts in the contactmap are at ideal distances, the total
//! energy is -1 (or 0 for the alternative functional form).
//!
template <typename CHAIN_TYPE>
class TermContactMap: public TermContactMapBase<TermContactMap<CHAIN_TYPE>, CHAIN_TYPE> {
private:

     //! For convenience, define local TermContactMap
     typedef phaistos::TermContactMapBase<TermContactMap<CHAIN_TYPE>,CHAIN_TYPE> TermContactMapBase;

protected:

     //! For convenience, define local Contact
     typedef typename TermContactMapBase::Contact Contact;

public:

     //! Local settings class.
     const class Settings: public TermContactMapBase::Settings {
     public:

          //! Whether to use the distance^2 version of the energy
          bool dist_squared;

          //! Width of energy potential
          double potential_width;

          //! Constructor
          Settings(bool dist_squared = false,
                   double potential_width = 1.0)
               :dist_squared(dist_squared),
                potential_width(potential_width) {}
          
          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "dist-squared: " << settings.dist_squared << "\n";
               o << "potential-width: " << settings.potential_width << "\n";
               o << (typename EnergyTerm<CHAIN_TYPE>::Settings &)settings;
               return o;
          }                    
     } settings;
     
     //! Constructor
     //! \param chain chain for contact map 
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermContactMap(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : TermContactMapBase(chain, "contact-map", settings, random_number_engine),
            settings(settings) {
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermContactMap(const TermContactMap &other, 
                    RandomNumberEngine *random_number_engine,
                    int thread_index, CHAIN_TYPE *chain)
          : TermContactMapBase(other, random_number_engine, thread_index, chain),
            settings(other.settings) {}


     //! Calculate contact energy given two residues.
     //! \param contact contact data
     //! \param r1 Residue 1
     //! \param r2 Residue 2
     //! \param dist_squared flag for energy ~ distance^2
     //! \return contact energy
     double contact_map_energy(const Contact &contact, 
                               Residue *r1, Residue *r2, bool dist_squared = false){
          double d=((*r1)[contact.a1]->position-(*r2)[contact.a2]->position).norm();
          return contact_map_energy(d, dist_squared);
     }

     //! Energy contribution for one contact
     //! \param contact contact data
     //! \param r distance between atom 1 and atom 2
     //! \param dist_squared flag for energy ~ distance^2
     //! \return energy
     double contact_map_energy(const Contact &contact, 
                               double r, bool dist_squared = false){
          //std::cout << r << " " << distance << " " << width << " " << dist_squared << "\n";
          return dist_squared ? pow((r-contact.distance)/settings.potential_width, 2)*contact.weight : 
               -std::exp(-pow((r-contact.distance)/settings.potential_width, 2)*contact.weight);
     }

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move
     double evaluate(MoveInfo *move_info=NULL) {

          double res=0.0;
          int no=0;
          for(unsigned int i=0; i<this->contact_map.size(); i++){
               if(this->contact_map[i].residue_index1 >= (this->chain)->size() || 
                  this->contact_map[i].residue_index2 >= (this->chain)->size())
                    continue;
               if(!(*(this->chain))[this->contact_map[i].residue_index1].has_atom(this->contact_map[i].atom_type1))
                    continue;
               if(!(*(this->chain))[this->contact_map[i].residue_index2].has_atom(this->contact_map[i].atom_type2))
                    continue;
               Vector_3D r1 = (*(this->chain))(this->contact_map[i].residue_index1, this->contact_map[i].atom_type1)->position;
               Vector_3D r2 = (*(this->chain))(this->contact_map[i].residue_index2, this->contact_map[i].atom_type2)->position;
               double d = (r1 - r2).norm();
               res += contact_map_energy(this->contact_map[i], d, settings.dist_squared);
               no++;
          }
          return res/this->normalization_factor;
     }
};



//! Observable specialization for TermContactMap
template <typename CHAIN_TYPE>
class Observable<TermContactMap<CHAIN_TYPE> >: public TermContactMap<CHAIN_TYPE>, public ObservableBase {

     std::vector<std::vector<std::map<std::pair<definitions::AtomEnum, definitions::AtomEnum>,std::pair<double,double> > > > contact_average;

     //! Number of evaluations/total weight (maximum number of occurrences of a contact in average mode)
     double weight_sum_total;

public:

     //! Local settings class.
     const class Settings: public TermContactMap<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Whether to output contact map in verbose mode
          bool verbose;

          //! Whether to calculate persistency of contacts - rather than a new contact map in each iteration
          bool average_mode;

          //! In average_mode: only include contacts which are present more frequently than this cutoff value
          double persistency_cutoff;

          //! In average_mode: whether to use the occurence-counts as weights
          bool counts_as_weights;

          //! Constructor. Defines default values for settings object.
          Settings(bool verbose=false,
                   bool average_mode=false,
                   double persistency_cutoff=0.0,
                   bool counts_as_weights=false)
               : verbose(verbose),
                 average_mode(average_mode),
                 persistency_cutoff(persistency_cutoff),
                 counts_as_weights(counts_as_weights){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "verbose:" << settings.verbose << "\n";
               o << "average-mode:" << settings.average_mode << "\n";
               o << "persistency-cutoff:" << settings.persistency_cutoff << "\n";
               o << "counts-as-weights:" << settings.counts_as_weights << "\n";
               o << static_cast<typename TermContactMap<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     
     //! Constructor.
     //! \param energy_term ContactMap energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermContactMap<CHAIN_TYPE> &energy_term, 
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermContactMap<CHAIN_TYPE>(energy_term),
            weight_sum_total(0),
            settings(dynamic_cast<const Settings&>(settings)) {

          if (this->settings.average_mode) {
               contact_average.resize(
                    this->chain->size(),
                    std::vector<std::map<std::pair<definitions::AtomEnum,definitions::AtomEnum>,std::pair<double,double> > >(this->chain->size()));
          }          

          if (this->settings.pdb_filename != "") {
               std::cerr << "WARNING (ContactMap observable): fixed pdb structures specified using pdb-file option - observable will not reflect current chain\n";
          }
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index, typename TermContactMap<CHAIN_TYPE>::ChainType *chain)
          : TermContactMap<CHAIN_TYPE>(other, thread_index, chain),
            weight_sum_total(0),
            settings(other.settings) {

          if (settings.average_mode) {
               contact_average.resize(
                    chain->size(),
                    std::vector<std::map<std::pair<definitions::AtomEnum,definitions::AtomEnum>,std::pair<double,double> > >(chain->size()));
          }          
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermContactMap<CHAIN_TYPE> *clone(int thread_index=0, typename TermContactMap<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermContactMap<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL, PHAISTOS_LONG_LONG current_iteration=0, bool register_only=false) {

          std::vector<typename TermContactMap<CHAIN_TYPE>::Contact> *contact_map_active = &this->contact_map;

          double weight = this->get_weight();

          // If a specific pdb file is selected from the command line, we repeatedly report the contact map for that
          // Otherwise, it is recalculated
          if (settings.pdb_filename == "") {

               this->construct_contact_map(this->chain, this->chain);

               weight_sum_total += weight;

               if (settings.average_mode) {
                    for (unsigned int i=0; i<this->contact_map.size(); ++i) {
                         typename TermContactMap<CHAIN_TYPE>::Contact &contact = this->contact_map[i];
                         std::map<std::pair<definitions::AtomEnum,definitions::AtomEnum>,std::pair<double,double> > &atom_type_map = 
                              contact_average[contact.residue_index1][contact.residue_index2];
                         std::pair<definitions::AtomEnum,definitions::AtomEnum> key = 
                              std::make_pair(contact.atom_type1,contact.atom_type2);
                         if (atom_type_map.count(key) == 0) {
                              atom_type_map[key] = std::make_pair(0,0);
                         }
                         std::pair<double,double> &entry = atom_type_map[key];
                         double &distance_sum = entry.first;
                         double &weight_sum = entry.second;
                         weight_sum += weight;
                         distance_sum += weight*contact.distance;
                    }

                    if (!register_only) {
                         contact_map_active = new std::vector<typename TermContactMap<CHAIN_TYPE>::Contact>;
                         for (unsigned int i=0; i<contact_average.size(); ++i) {
                              int residue_index1 = i;
                              for (unsigned int j=0; j<contact_average[i].size(); ++j) {
                                   int residue_index2 = j;
                    
                                   std::map<std::pair<definitions::AtomEnum,definitions::AtomEnum>,std::pair<double,double> >::iterator it
                                        = contact_average[i][j].begin();
                                   for (; it != contact_average[i][j].end(); ++it) {
                                        definitions::AtomEnum atom_type1= it->first.first;
                                        definitions::AtomEnum atom_type2= it->first.second;
                                        double distance_sum = it->second.first;
                                        double weight_sum = it->second.second;

                                        double persistency = weight_sum/weight_sum_total;
                                        if (!std::isfinite(persistency))
                                             persistency = 0;

                                        // Omit contacts with low persistency
                                        if (persistency < settings.persistency_cutoff)
                                             continue;

                                        // When using a cutoff, the weights of the reported contacts are set to 1.0
                                        double contact_weight = persistency;
                                        if (settings.persistency_cutoff > 0)
                                             contact_weight = 1.0;
                                        if (settings.counts_as_weights) 
                                             contact_weight = weight_sum;

                                        double distance_average = distance_sum/weight_sum;
                                        if (!std::isfinite(distance_average))
                                             distance_average = 0;

                                        typename TermContactMap<CHAIN_TYPE>::Contact contact(
                                             residue_index1, atom_type1, 
                                             residue_index2, atom_type2, 
                                             distance_average, contact_weight);
                              
                                        contact_map_active->push_back(contact);
                                   }
                              }
                         }
                    }
               }
          }

          std::string output;
          if (!register_only) {
               if (settings.verbose) {
                    output = this->contact_map_to_string_verbose(*contact_map_active);
               } else {
                    output = this->contact_map_to_string(*contact_map_active);
               }
          }
               
          return output;
     }

};


}
#endif
