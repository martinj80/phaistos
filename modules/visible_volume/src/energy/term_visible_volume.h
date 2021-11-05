// term_visible_volume.h --- Visible Volume energy term
// Copyright (C) 2012-2013 Kristoffer Enøe Johansson
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

#ifndef VISIBLE_VOLUME_ENERGY_H
#define VISIBLE_VOLUME_ENERGY_H

#include "models/visible_volume/sphere.h"

#include "protein/iterators/residue_iterator.h"

#include "energy/energy_term.h"

namespace phaistos {
     

//! Parameters for hydrophobic effect of amino acid sidechains solvation from FoldX (kcal/mol):
//! Raphael Guerois, Jens Erik Nielsen and Luis Serrano:
//! "Predicting Changes in the Stability of Proteins and Protein Complexes: A Study of More Than 1000 Mutations"
//! J. Mol. Biol. (2002) 320, 369-387
const double hydrophobic_solvation_param[20] =
// Ala   Cys  Asp  Glu   Phe  Gly  His   Ile  Lys   Leu   Met  Asn   Pro  Gln  Arg  Ser  Thr   Val   Trp   Tyr
{-0.40,-0.59,2.28,1.88,-2.22,0.11,0.08,-1.60,2.10,-1.60,-1.35,1.00,-1.10,0.60,2.50,0.67,0.27,-1.20,-2.32,-1.42};


// //! Definition of neighbor representation
// enum NeighborModeEnum {ALL_ATOM_NEIGHBORS=0, CHEMICAL_GROUPS_NEIGHBORS=1, SIDE_CHAIN_NEIGHBORS=2,
//                         NEIGHBOR_MODE_ENUM_SIZE};

// static const std::string neighbor_mode_enum_name[] = {"all-atom", "chemical-groups", "side-chain"};

// //! Input NeighborModeEnum from string
// inline std::istream &operator>>(std::istream &input, NeighborModeEnum &g) {
//      std::string raw_string;
//      input >> raw_string;

//      for (unsigned int i=0; i<NEIGHBOR_MODE_ENUM_SIZE; ++i) {
//           if (raw_string == neighbor_mode_enum_name[i]) {
//                g = NeighborModeEnum(i);
//           }
//      }     
//      return input;
// }

// //! Output NeighborModeEnum
// inline std::ostream &operator<<(std::ostream &o, const NeighborModeEnum &g) {
//      o << neighbor_mode_enum_name[static_cast<unsigned int>(g)];
//      return o;
// }


//! Visible volume energy term
template <typename CHAIN_TYPE>
class TermVisibleVolume:public EnergyTermCommon<TermVisibleVolume<CHAIN_TYPE>, CHAIN_TYPE> {

public:
     
     //! Calculate the ideal CB position from backbone N, CA and C
     Vector_3D calc_cb_pos(Residue *res) {

          // get Backbone N atom
          assert(res->has_atom(definitions::N));
          Vector_3D &n_pos = (*res)[definitions::N]->position;

          // get Backbone Ca atom
          assert(res->has_atom(definitions::CA));
          Vector_3D &ca_pos = (*res)[definitions::CA]->position;

          // get Backbone C atom
          assert(res->has_atom(definitions::C));
          Vector_3D &c_pos = (*res)[definitions::C]->position;

          // calculate Cb position assuming ideal CA geometry
          double bond_length = 1.53;           //enghHuberConstants.h L71
          double angle = 89.83*M_PI/180.0;    //enghHuberConstants.h L660
          double dihedral = 52.25*M_PI/180.0; //enghHuberConstants.h L1262
          phaistos::Vector_3D D(bond_length * cos(M_PI-angle),
                                bond_length * cos(M_PI-dihedral) * sin(M_PI-angle),
                                bond_length * sin(M_PI-dihedral) * sin(M_PI-angle));
          phaistos::Vector_3D bc = (c_pos-n_pos).normalize();
          phaistos::Vector_3D no = cross_product((c_pos-ca_pos), bc).normalize();
          phaistos::Vector_3D nbc = cross_product(bc,no);
          phaistos::Matrix_3D basis_change(bc, nbc, no, true);
          D = basis_change*D + c_pos;
          D -= (c_pos - ca_pos);
          return D;
     }

     //! Extract CB position from residue
     Vector_3D get_cb_pos(Residue *res) {
          if ( res->has_atom(definitions::CB) ) {
               return (*res)[definitions::CB]->position;
          } else {
               return calc_cb_pos(res);
          }
     }

     //! Find the geometrical center of a residue
     Vector_3D calc_geo_center(Residue *res) {

          // if pseudo sidechain is present use this position
          if (res->has_atom(definitions::PS))
               return (*res)[definitions::PS]->position;

          // use Cb position for alanine
          if (res->residue_type == definitions::ALA)
               return get_cb_pos(res);

          // use Ca position for glysine
          if (res->residue_type == definitions::GLY) {
               assert(res->has_atom(definitions::CA));
               return (*res)[definitions::CA]->position;
          }

          // geometric center to be returned
          Vector_3D gc(0.0, 0.0, 0.0);
          int n_atoms = res->atoms.size();
          int sc_heavy_atoms = 0;

          for (int i=0; i<n_atoms; i++) {
               Atom *atom = res->atoms[i];
               // sum side chain heavy atoms excluding CB
               if (atom->atom_type > definitions::CB && atom->atom_type < definitions::PS) {
                    gc += atom->position;
                    sc_heavy_atoms += 1;
               }
          }
          
          return gc/sc_heavy_atoms;
     }     

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermVisibleVolume<CHAIN_TYPE>, CHAIN_TYPE> EnergyTermCommon;

     void init() {

          // construct sphere object
          sphere = new Sphere(this->settings.sphere_radius, this->settings.min_sphere_points);

          // number of residues
          chain_size = this->chain->size();

          // init cache
          volume = std::vector<double>(chain_size, 0.0);
          angexp = std::vector<double>(chain_size, 0.0);
          fscn = std::vector<int>(chain_size, 0);
     }
          
protected:

     //! Sphere object to calculate visible volume derived measures
     Sphere *sphere;

     //! Number of residues
     int chain_size;

public:

     //! Visible Volume
     std::vector<double> volume;

     //! Solvent accessible surface area
     std::vector<double> angexp;

     //! Contact number
     std::vector<int> fscn;

     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! Sphere radius
          double sphere_radius;

          //! Minimum number of angle points
          int min_sphere_points;

          // //! Include side chain atoms in residue contact center positioning
          // bool design_mode;

          // //! How to represent neighbors
          // NeighborModeEnum neighbor_mode;

          // //! Number of iterations between neighbor list update
          // int neighbor_list_every;
          
          //! Disk radius of atom neighbors
          double atom_radius;

          //! Calculate visible volume for each residue
          bool volume;
          
          //! Calculate surface accessible surface area for each residue
          bool angexp;
          
          //! Calculate contact number, i.e. total number of neighbors in contact
          bool fscn;
          
          //! Minimum number of angle points to define contact
          int fscn_threshold;
          
          //! Constructor
          Settings(double sphere_radius=11.0, int min_sphere_points=1000,
                   // bool design_mode=false, NeighborModeEnum neighbor_mode=ALL_ATOM_NEIGHBORS,
                   // int neighbor_list_every=1,
                   double atom_radius=1.8,
                   bool volume=true, bool angexp=false, bool fscn=false, int fscn_threshold=35)
               : sphere_radius(sphere_radius),
                 min_sphere_points(min_sphere_points),
                 // design_mode(design_mode),
                 // neighbor_mode(neighbor_mode),
                 // neighbor_list_every(neighbor_list_every),
                 atom_radius(atom_radius),
                 volume(volume),
                 angexp(angexp),
                 fscn(fscn),
                 fscn_threshold(fscn_threshold) {
          }
          
          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "sphere-radius:" << settings.sphere_radius << "\n";
               o << "min-sphere-points:" << settings.min_sphere_points << "\n";
               // o << "design-mode:" << settings.design_mode << "\n";
               // o << "neighbor-mode:" << settings.neighbor_mode << "\n";
               // o << "neighbor-list-every:" << settings.neighbor_list_every << "\n";
               // o << "atom-radius:" << settings.atom_radius << "\n";
               o << "volume:" << settings.volume << "\n";
               o << "angexp:" << settings.angexp << "\n";
               o << "fscn:" << settings.fscn << "\n";
               o << "cn-min-counts:" << settings.fscn_threshold << "\n";
               return o;
          }                    
     } settings;
     

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermVisibleVolume(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                       RandomNumberEngine *random_number_engine = &random_global) 
          : EnergyTermCommon(chain, "visible-volume", settings, random_number_engine),
            settings(settings) {

          init();
     };
     
     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain     
     TermVisibleVolume(const TermVisibleVolume &other,
                       RandomNumberEngine *random_number_engine,
                       int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            settings(other.settings) {

          init();          
     };
     
     //! Destructor
     ~TermVisibleVolume() {
          delete sphere;
     };

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move     
     double evaluate(MoveInfo *move_info=NULL) {

          // visible volume energy
          double volume_energy = 0.0;
          
          // SASA energy
          double angexp_energy = 0.0;

          // Contact number
          int fscn_energy = 0;

          // calc visible volume for all residues
          for (int i=0; i<chain_size; i++) {
               ResidueFB *resi = &(*this->chain)[i];
               int resi_type = resi->residue_type;

               // reset sphere
               sphere->reset();

               // find contact center
               Vector_3D cc;
               if (resi_type == definitions::ALA)
                    cc = get_cb_pos(resi);
               else if (resi_type == definitions::GLY)
                    cc = calc_cb_pos(resi);
               else
                    cc = calc_geo_center(resi);
               
               // add neighbors
               for (int j=0; j<chain_size; j++) {

                    // no shadow from self
                    if (i==j)
                         continue;

                    ResidueFB *resj = &(*this->chain)[j];
                    int resj_type = (*this->chain)[j].residue_type;

                    // small optimization
                    Vector_3D cbj = get_cb_pos(resj);
                    double dx = cbj[0]-cc[0];
                    double dy = cbj[1]-cc[1];
                    double dz = cbj[2]-cc[2];
                    double cut = settings.sphere_radius + 5.0;
                    if (dx*dx + dy*dy + dz*dz > cut*cut)
                         continue;
                    
                    int atom_size = resj->atoms.size();
                    for (int a=0; a<atom_size; a++) {
                         Atom *atom = resj->atoms[a];
                         int atom_type = atom->atom_type;
                         
                         // atom position relative to contact center
                         Vector_3D pos = atom->position - cc;

                         // put in sphere object
                         sphere->shade(pos.get_array(), settings.atom_radius, j, a, resj_type, atom_type);
                    }
               }

               if (settings.volume) {
                    // Count exposed angle points
                    std::vector<double> vol = sphere->calc_visible_volume();

                    // Store residue i exposure
                    volume[i] = vol[0];

                    volume_energy += volume[i]; // times some value...
               }

               // if (settings.angexp) { // Currently, only angexp is implemented as an energy so always return angexp energy
                    // Count exposed angle points
                    std::vector<std::vector<int> > exposed_angle_points = sphere->atom_type_histogram(0);

                    // Store residue i exposure
                    angexp[i] = exposed_angle_points[0][0]*1.0/sphere->n_points;

                    // Add angular exposure energy contribution for this residue
                    angexp_energy += (1.0 - angexp[i]) * hydrophobic_solvation_param[resi_type];
               // }

               if (settings.fscn) {
                    // Calculate contacts multinomial
                    std::vector<std::vector<std::vector<int> > > contacts = sphere->get_contacts();

                    int contact_number = 0;
                    if (contacts.size() > 0) {
                         const int c_size = contacts[0].size();
                         for (int c=0; c<c_size; c++) {
                              if (contacts[0][c][2] >= settings.fscn_threshold) {
                                   contact_number += 1;
                              }
                         }
                    }
                    
                    // Store residue i contact number
                    fscn[i] = contact_number;

                    // Add contact number energy contribution for this residue
                    fscn_energy += contact_number; // times some number
               }
          }

          // return volume_energy + angexp_energy + fscn_energy;

          return angexp_energy; // Currently, only angexp is implemented as an energy so always return angexp energy
     }
};


//! Observable specialization for TermVisibleVolume
template <typename CHAIN_TYPE>
class Observable<TermVisibleVolume<CHAIN_TYPE> >: public TermVisibleVolume<CHAIN_TYPE>, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermVisibleVolume<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename TermVisibleVolume<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<const ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     
     //! Constructor.
     //! \param energy_term VisibleVolume energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermVisibleVolume<CHAIN_TYPE> &energy_term,
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermVisibleVolume<CHAIN_TYPE>(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index,
                typename TermVisibleVolume<CHAIN_TYPE>::ChainType *chain)
          : TermVisibleVolume<CHAIN_TYPE>(other, thread_index, chain),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermVisibleVolume<CHAIN_TYPE> *clone(int thread_index=0,
                                          typename TermVisibleVolume<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermVisibleVolume<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     //! Make observation.
     //! \param move_info Object containing information about the last executed move     
     //! \param current_iteration Index of current iteration
     //! \param register_only Whether this call should only register data (not output it)
     virtual std::string observe(MoveInfo *move_info=NULL,
                                 PHAISTOS_LONG_LONG current_iteration=0,
                                 bool register_only=false) {

          // this->evaluate_weighted(move_info);
          this->evaluate(move_info);

          if (settings.angexp+settings.volume+settings.fscn == 0)
               return "No measures selected: Set e.g. --observable-visible-volume-volume=1";
          
          std::string output = "";

          // strip output-target for output-mode
          std::string ot = settings.output_target;
          ot = ot.substr(0, ot.find('#'));

          if (ot == "stdout" || ot == "cout" || ot == "stderr" || ot == "cerr") {

               // screen dump format
               output += "\nResidue  ";
               if (settings.volume)
                    output += "    Volume";
               if (settings.angexp)
                    output += "    AngExp";
               if (settings.fscn  )
                    output += " FSCN";
               output += "\n";
          
               char buffer[128];
               for (int i=0; i<this->chain_size; i++) {
                    int res_type = (*this->chain)[i].residue_type;
                    int seq_index = (*this->chain)[i].index_res_seq;
                    sprintf(buffer,"%3s-%04d:", definitions::residue_name[res_type], seq_index);
                    output += buffer;
                    if (settings.volume) {
                         sprintf(buffer,"  %8.3f", this->volume[i]);
                         output += buffer;
                    }
                    if (settings.angexp) {
                         sprintf(buffer,"  %8.5f", this->angexp[i]);
                         output += buffer;
                    }
                    if (settings.fscn) {
                         sprintf(buffer,"  %3d", this->fscn[i]);
                         output += buffer;
                    }
                    output += "\n";
               }
          } else if (ot == "pdb-b-factor") {
               
               for (ResidueIterator<CHAIN_TYPE> it(*this->chain); !it.end(); ++it) {
                    std::string output_entry = ObservableBase::vector_output_tag(*it);
                    
                    // Only one output measure is meningful as PDB b-factor output
                    assert((settings.volume+settings.angexp+settings.fscn) == 1);

                    if (settings.volume) {
                         // output_entry += boost::lexical_cast<std::string>(this->volume[it->index]);
                         char buffer[12];
                         std::sprintf(buffer, "%.1f", this->volume[it->index]);
                         output_entry += std::string(buffer);
                    } else if (settings.angexp)
                         output_entry += boost::lexical_cast<std::string>(this->angexp[it->index]);
                    else //settings.fscn
                         output_entry += boost::lexical_cast<std::string>(this->fscn[it->index]);
                    
                    if (output != "")
                         output += ",";
                    output += output_entry;
               }
          } else {
               char buffer[128];
               // logfile and pdbfile dump format
               for (ResidueIterator<CHAIN_TYPE> it(*this->chain); !it.end(); ++it) {
                    std::string output_entry = ObservableBase::vector_output_tag(*it);
                    if (settings.volume) {
                         sprintf(buffer,"%.3f",this->volume[it->index]);
                         output_entry += buffer;
                         if (settings.angexp+settings.fscn > 0)
                              output_entry += ":";
                    }
                    if (settings.angexp) {
                         sprintf(buffer,"%.5f",this->angexp[it->index]);
                         output_entry += buffer;
                         if (settings.fscn)
                              output_entry += ":";
                    }
                    if (settings.fscn) {
                         sprintf(buffer,"%d",this->fscn[it->index]);
                         output_entry += buffer;
                    }
                    if (output != "")
                         output += ",";
                    output += output_entry;
               }
          }
          
          return output;
     }
     
};
     
}

#endif
