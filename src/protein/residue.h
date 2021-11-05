// residue.h --- residue class
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


#include "atom.h"

#ifndef RESIDUE_H
#define RESIDUE_H

#include <vector>
#include <queue>
#include "utils/vector_matrix_3d.h"

#include "definitions.h"

namespace phaistos {

// Forward declarations
template <typename RESIDUETYPE> class Chain;
class ChainCA;
class ChainFB;
namespace deprecated {
class AtomIteratorAll;
class AtomIteratorCA;
class AtomIteratorBackbone;
class AtomIteratorSidechain;
}
class Residue;


//! Residue (amino acid) class
class Residue {
public:

     //! Local typedef indicating the chaintype that a given residue type corresponds to
     //! The base class does not correspond to any well-defined chain, and is therefore set
     //! to void
     typedef void ChainType;

     //! Type of residue
     definitions::ResidueEnum residue_type;

     //! Where residue is located in the chain (INTERNAL | CTERM | NTERM)
     definitions::TerminalEnum terminal_status;

     //! Vector of atoms (active)
     std::vector<Atom *> atoms;

     //! Vector of atoms (deactivated)
     std::vector<Atom *> atoms_inactive;

     //! Index of residue in enclosing chain
     int index;

     //! Index of residue in enclosing chain - according to pdb file
     int index_res_seq;

     //! Specifies which atoms are at which index
     int atom_index[definitions::ATOM_ENUM_SIZE];        

     //! Internal variable for keeping track of atoms_length
     //!  - detects if the atoms length suddenly changes     
     int atoms_length;

     //! Internal enum: BEGIN,END
     enum RangeEnum{BEGIN=0, END=1};

     //! Keep track of where each iteration mode begins and ends in current residue
     int iteration_range_indices[definitions::ITERATE_ENUM_SIZE][2];

     //! Pointer to atoms corresponding to chi-angles
     std::vector<std::pair<Atom*,definitions::AngleEnum> > chi_atoms;

     //! Pointer to atoms containing minor degrees of freedom which
     //! can be modified without effecting the backbone     
     std::vector<std::pair<Atom*,definitions::AngleEnum> > minor_dof_atoms;

     //! Status of sidechain (true: sidechain present, false: sidechain absent (or deactivated))
     bool sidechain_status;

     //! Covalent distances between atoms in residue
     std::vector<int> atom_chain_distances;
     


     //! Constructor
     //!
     //! \param residue_type Type of residue
     //! \param index Index of residue in enclosing chain
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file
     Residue(definitions::ResidueEnum residue_type, int index, int index_res_seq);

     //! Copy Constructor
     //!
     //! \param other Source object from which copy is made.
     //! \param index Index of residue in enclosing chain
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file
     Residue(Residue &other, int index, int index_res_seq);

     //! Destructor
     virtual ~Residue();

     //! Initializer
     void init();

     //! Initialize bond lengths
     void init_bond_lengths();

     //! Initialize bond angles and dihedrals
     void init_angles();

     //! Number of atoms in atoms vector
     //!
     //! \return Size of active atoms vector
     int size () const;

     //! Overload indexing operator - by atom_type 
     //! 
     //! \param atom_type Type of atom (e.g. res[N])
     //! \return atom pointer
     Atom *operator[](const definitions::AtomEnum atom_type);

     //! Overload indexing operator - by index
     //!
     //! \param index atom index in residue
     //! \return atom pointer
     Atom *operator[](unsigned int index);

     //! Test for presence of atom_type in residue
     //!
     //! \param atom_type Type of atom
     //! \return True if atom_type is present in residue
     bool has_atom(const definitions::AtomEnum atom_type);

     //! Test for presence of side-chain
     //!
     //! \return True if residue has a sidechain
     bool has_sidechain();

     //! Ensure that atom names adhere to our internal PDB standard
     void standardize_atoms();

     //! Rename atom
     //!
     //! \param atom_type Current type of atom
     //! \param possible_replacements Vector of possible alternative atom types
     //! \param preferred_type Preferred alternative (in case there are several possibilities.
     void rename_atom(definitions::AtomEnum atom_type, 
                      std::vector<definitions::AtomEnum> possible_replacements, 
                      definitions::AtomEnum preferred_type=definitions::XX);
     
     //! Sort atom vector so that they are in correct PDB order and remove NULL entries
     void sort_atoms();

     //! Sort inactive atom vector so that they are in correct PDB order and remove NULL entries
     void sort_atoms_inactive();
     
     //! Update backbone dihedral and bondangle values based on current positions
     //! \param fix_end_points If true, the angle value at the non-well-defined 
     //!        endpoint positions are set to zero.
     void update_angles(bool fix_end_points=true);

     //! Update positions based on current angle values - backbone positions only
     //!
     //! \param direction Iteration direction along the chain
     void update_positions_backbone(int direction);

     //! Update positions based on current angle values - non-backbone positions only
     void update_positions_non_backbone();

     //! Set ideal bond lengths and bond angles, and update positions
     //!
     //! \param atom_selection Which atoms to idealize     
     void idealize(definitions::AtomSelectionEnum atom_selection=definitions::ALL_ATOMS);

     //! Return distance along chain from atom1 to atom2 (in this residue)
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \return distance along chain (number of covalent bonds)
     template <typename CHAIN_TYPE>
     int chain_distance(Atom *atom1, Atom *atom2);

     //! Get current angle values
     //!
     //! \param include_omega Whether omega angles should be included in the output vector
     //! \param include_bond_angles Whether bond angles should be included in the output vector.
     //! \return vector of angles.
     virtual std::vector<double> get_angles(bool include_omega=false, 
                                            bool include_bond_angles=false)=0;

     //! Set angles
     //!
     //! \param angles Vector of angles
     virtual void set_angles(std::vector<double> angles)=0;
     
     //! Retrieve positions for all atoms
     //!
     //! \return vector of 3D-coordinates
     std::vector<Vector_3D> get_positions();

     //! Get neighbouring residue
     //!
     //! \param offset Distance along chain to neighbour (nth neighbour)
     //! \return Residue pointer
     virtual Residue *get_neighbour(int offset)=0;
     
     //! Apply translation operation to all atom positions
     //!
     //! \param translation translation vector
     void translate(Vector_3D translation);

     //! Apply rotation operation to all atom positions
     //!
     //! \param rotation rotation matrix
     void rotate(Matrix_3D rotation);

     //! Check consistency between the two representations of the atom. 
     //! The chain is overspecified since it keeps track of both dihedral angles and positions.
     //! This method checks that the two representations are consistent.
     void check_consistency();

     //! Output residue as string
     //!
     //! \param o Output stream
     void output(std::ostream &o);

     //! Overload output operator for Residue
     friend std::ostream & operator<<(std::ostream &o, Residue &r);

     //! Output residue in PDB format
     //!
     //! \param counter_offset Offset to add to atom index
     //! \param chain_number The chain id
     //! \param b_factor_string Optional vector string to output as b factors
     //! \return Output in pdb format     
     std::string output_as_pdb(int counter_offset, const int chain_number=1, std::string *b_factor_string=NULL);

     //! Clears the cache containing information about neighbouring residues
     void clear_neighbour_cache();

     
     /////////////////////////////
     // Sidechain functionality //
     /////////////////////////////
     
     //! Set chi angles or pseudo-sidechain degrees of freedom
     //!
     //! \param dof_value_vector Vector of input values
     //! \param mode SIDECHAIN_ATOMS -> dofValueVector contains chi angles
     //!             PSEUDO_SIDECHAIN_ATOMS -> dofValueVector contains (bondlength, angle, dihedral)     
     void set_sidechain_dof_values(std::vector<double> &dof_value_vector,
                                   definitions::AtomSelectionEnum mode=definitions::SIDECHAIN_ATOMS);

     //! Get chi angles or pseudo-sidechain degrees of freedom
     //!
     //! \param mode SIDECHAIN_ATOMS -> dofValueVector contains chi angles
     //!              PSEUDO_SIDECHAIN_ATOMS -> dofValueVector contains (bondlength, angle, dihedral)
     //! \return vector of degree-of-freedom values
     std::vector<double> get_sidechain_dof_values(definitions::AtomSelectionEnum mode=definitions::ALL_ATOMS);

     //! Set values of minor degrees of freedom in residue which can be modified freely 
     //! without modifying the backbone
     //!
     //! \param dof_value_vector Vector of input values
     void set_minor_dof_values(const std::vector<double> &dof_value_vector);

     //! Get values of bondangles in residue which can be modified freely 
     //! without modifying the backbone
     //!
     //! \return vector of degree-of-freedom values
     std::vector<double> get_minor_dof_values();

     //! Sets up one of the pointers in the chi-pointer array
     //!
     //! \param atom Atom pointer
     //! \param index Index of chi angle
     void set_chi_atom(Atom *atom, int index=-1);

     //! Sets up one of the pointers in the minor-dof-pointer array
     //!
     //! \param atom Atom pointer
     //! \param angle_type Type of angle
     //! \param index of minor dof
     void set_minor_dof_atom(Atom *atom, definitions::AngleEnum angle_type, int index=-1);

     //! Remove an atom from the chain
     //!
     //! \param atom_type Type of atom
     //! \param execute_sort_atoms Whether to sort atoms after removal
     void remove_atom(definitions::AtomEnum atom_type, bool execute_sort_atoms);

     //! Deactivate an atom (move it to atoms_inactive)
     //!
     //! \param atom_type Type of atom
     //! \param execute_sort_atoms Whether to sort atoms after removal
     void deactivate_atom(definitions::AtomEnum atom_type, bool execute_sort_atoms);

     //! Activate an atom (move it from atoms_inactive to atoms)
     //!
     //! \param atom_type Type of atom
     //! \param execute_sort_atoms Whether to sort atoms after removal
     bool activate_atom(definitions::AtomEnum atom_type, bool execute_sort_atoms);

     //! Remove a selection of atoms from the chain
     //!
     //! \param atom_selection Specifies which atom to remove
     //! \param execute_sort_atoms Whether to sort atoms after removal
     void remove_atoms(definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS+
                                                                      definitions::BACKBONE_H_ATOMS+
                                                                      definitions::CB_ATOMS+
                                                                      definitions::SIDECHAIN_ATOMS+
                                                                      definitions::NON_BACKBONE_H_ATOMS), bool execute_sort_atoms=true);

     //! Activates a selection of atoms in the chain
     //!
     //! \param atom_selection Specifies which atom to activate
     //! \param execute_sort_atoms Whether to sort atoms after removal
     void activate_atoms(definitions::AtomSelectionEnum atom_selection=(definitions::ALL_ATOMS), bool execute_sort_atoms=true);

     //! Deactivates a selection of atoms in the chain
     //!
     //! \param atom_selection Specifies which atom to deactivate
     //! \param execute_sort_atoms Whether to sort atoms after removal
     void deactivate_atoms(definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS+
                                                                          definitions::BACKBONE_H_ATOMS+
                                                                          definitions::CB_ATOMS+
                                                                          definitions::SIDECHAIN_ATOMS+
                                                                          definitions::NON_BACKBONE_H_ATOMS), bool execute_sort_atoms=true);

     //! Add a C_beta atom
     //!
     //! \param atom_selection Specifies which atom to add
     //! \param execute_sort_atoms Whether to sort atoms afterwards
     void add_CB(definitions::AtomSelectionEnum atom_selection=definitions::CB_ATOMS, bool execute_sort_atoms=true);

     //! Add full-atom sidechain
     //!
     //! \param atom_selection Specifies which atom to add
     //! \param chi_angles Values for the new degrees of freedom
     //! \param execute_sort_atoms Whether to sort atoms afterwards
     void add_sidechain(definitions::AtomSelectionEnum atom_selection=definitions::SIDECHAIN_ATOMS,
                        const std::vector<double> &chi_angles=std::vector<double>(), bool execute_sort_atoms=true);

     //! Add a pseudo-sidechain atom
     //!
     //! \param atom_selection Specifies which atom to add
     //! \param bond_length_and_angles Values for the new degrees of freedom
     //! \param execute_sort_atoms Whether to sort atoms afterwards
     void add_pseudo_sidechain(definitions::AtomSelectionEnum atom_selection=definitions::PSEUDO_SIDECHAIN_ATOMS,
                               const std::vector<double> &bond_length_and_angles=std::vector<double>(), bool execute_sort_atoms=true);

     //! Calculates the position of the side chain center of mass
     //!
     //! \return 3D-coordinate of sidechain center of mass
     Vector_3D calc_sidechain_center_of_mass();

     //! Calculates the position of the backbone center of mass
     //!
     //! \return 3D-coordinate of backbone center of mass
     Vector_3D calc_backbone_center_of_mass();

     //! calculates the position of the residue center of mass
     //!
     //! \return 3D-coordinate of residue center of mass
     Vector_3D calc_center_of_mass();

     //! Return pointer to the first atom, given an iteration mode
     //!
     //! \param iteration_mode Type of iteration
     //! \return Atom pointer
     Atom *get_first_atom(definitions::IterateEnum iteration_mode);

     //! Return pointer to the last atom, given an iteration mode
     //!
     //! \param iteration_mode Type of iteration
     //! \return Atom pointer
     Atom *get_last_atom(definitions::IterateEnum iteration_mode);
     
     //! Switch sidechain representation
     void toggle_sidechain();

     //! Switch to pseudo sidechains
     void toggle_sidechain_to_PS();

     //! Switch to full-atom sidechains
     void toggle_sidechain_to_full_atom();


protected:
     ///////////////
     // Iterators //
     ///////////////
     
     //! Internal vector iterator
     typedef std::vector<Atom *>::iterator Iterator;

     //! Internal vector iterator (reverse)
     typedef std::vector<Atom *>::reverse_iterator ReverseIterator;

};



//! C_alpha residue type
class ResidueCA: public Residue {
public:

     //! Local typedef indicating the chaintype that a given residue type corresponds to
     typedef ChainCA ChainType;

     //! Pointer to the chain in which the residue is located
     ChainCA *chain;

     //! Constructor - incomplete initialization
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file
     ResidueCA(definitions::ResidueEnum residue_type, Chain<ResidueCA> *chain, int index, int index_res_seq=uninitialized<int>());

     //! Constructor - incomplete initialization - angles must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param position_CA C_alpha position
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file     
     ResidueCA(definitions::ResidueEnum residue_type, ChainCA *chain, int index, Vector_3D position_CA, int index_res_seq=uninitialized<int>());


     //! Constructor - incomplete initialization - angles must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param positions Vector of (atom type, 3D-coordinate) pairs
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file     
     ResidueCA(definitions::ResidueEnum residue_type, ChainCA *chain, int index,
               std::vector<std::pair<definitions::AtomEnum, Vector_3D> > &positions, int index_res_seq=uninitialized<int>());
     
     
     //! Constructor - incomplete initialization - atom positions must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param theta Theta bond angle
     //! \param tau Tau dihedral angle
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file     
     ResidueCA(definitions::ResidueEnum residue_type, ChainCA *chain, 
               int index, double theta, double tau, int index_res_seq=uninitialized<int>());


     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file
     ResidueCA(ResidueCA &other, ChainCA *chain, int index, int index_res_seq=uninitialized<int>());

     //! Initializer
     void init();

     //! Destructor
     ~ResidueCA();

     //! Set tau dihedral angle
     //!
     //! \param tau Tau dihedral angle
     void set_tau(double tau);

     //! Get tau dihedral angle
     //!
     //! \return Tau dihedral angle
     double get_tau();

     //! Set theta bond angle
     //!
     //! \param theta Theta bond angle     
     void set_theta(double theta);

     //! Get theta bond angle
     //!
     //! \return Theta bond angle     
     double get_theta();

     //! Return angle vector
     //!
     //! \param include_omega Whether to include omega in vector. 
     //! This only makes sense for ResidueFB - but is included here to align the interface of the two)
     //! \param include_bond_angles Whether to include bond angles in vector
     //! \return Vector of angles
     std::vector<double> get_angles(bool include_omega=false, bool include_bond_angles=false);

     //! Set angle vector
     //!
     //! \param angles Vector of angles
     void set_angles(std::vector<double> angles);

     //! Get neighbouring residue
     //!
     //! \param offset Distance along chain to neighbour (nth neighbour)
     //! \return Residue pointer
     ResidueCA *get_neighbour(int offset);

     //! Add atoms to chain
     //!
     //! \param atom_selection Which atoms to add
     //! \param sidechain_dof_values Optionally set values for sidechain degrees of freedom
     //! \param execute_sort_atoms Whether to sort atoms after removal
     void add_atoms(definitions::AtomSelectionEnum atom_selection=definitions::CB_ATOMS,
                    const std::vector<double> &sidechain_dof_values=std::vector<double>(), bool execute_sort_atoms=true);

     //! Output residue as string
     //!
     //! \param o Output stream
     void output(std::ostream &o);

     //! Overload output operator for ResidueCa     
     friend std::ostream & operator<<(std::ostream &o, ResidueCA &r);
};


//! Full-atom backbone residue type
class ResidueFB: public Residue {
public:

     //! Local typedef indicating the chaintype that a given residue type corresponds to
     typedef ChainFB ChainType;

     //! Pointer to the chain in which the residue is located
     ChainFB *chain;

     //! Constructor - incomplete initialization
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file
     ResidueFB(definitions::ResidueEnum residue_type, Chain<ResidueFB> *chain, int index, int index_res_seq=uninitialized<int>());

     //! Constructor - incomplete initialization - angles must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param position_N N position
     //! \param position_CA C_alpha position
     //! \param position_C C' position
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file          
     ResidueFB(definitions::ResidueEnum residue_type, ChainFB *chain, int index, Vector_3D position_N,
               Vector_3D position_CA, Vector_3D position_C, int index_res_seq=uninitialized<int>());

     //! Constructor - incomplete initialization - angles must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param positions Vector of (atom type, 3D-coordinate) pairs
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file          
     ResidueFB(definitions::ResidueEnum residue_type, ChainFB *chain, int index, std::vector<std::pair<definitions::AtomEnum, Vector_3D> > &positions,
               int index_res_seq=uninitialized<int>());

     //! Constructor - incomplete initialization - atom positions must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param phi Phi dihedral angle
     //! \param psi Psi dihedral angle
     //! \param omega Omega dihedral angle
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file          
     ResidueFB(definitions::ResidueEnum residue_type, ChainFB *chain, int index, double phi,
               double psi, double omega=M_PI,
               int index_res_seq=uninitialized<int>());

     //! Constructor - incomplete initialization - atom positions must be initalized later
     //!
     //! \param residue_type Type of residue
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param angles Angle vector
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file               
     ResidueFB(definitions::ResidueEnum residue_type, ChainFB *chain, int index, const std::vector<double> &angles,
               int index_res_seq=uninitialized<int>());

     //! Copy constructor
     //!
     //! \param other Source object from which copy is made.
     //! \param chain Enclosing chain
     //! \param index Index of residue in enclosing chain
     //! \param index_res_seq Index of residue in enclosing chain - according to pdb file
     ResidueFB(ResidueFB &other, ChainFB *chain, int index, int index_res_seq=uninitialized<int>());

     //! Initializer
     void init();

     //! Destructor
     ~ResidueFB();

     //! Set phi dihedral angle
     //!
     //! \param phi Phi dihedral angle
     void set_phi(double phi);

     //! Get phi dihedral angle
     //!
     //! \return Phi dihedral angle
     double &get_phi();

     //! Set psi dihedral angle
     //!
     //! \param psi Psi dihedral angle     
     void set_psi(double psi);

     //! Get psi dihedral angle
     //!
     //! \return Psi dihedral angle     
     double &get_psi();

     //! Set omega dihedral angle
     //!
     //! \param omega Omega dihedral angle     
     void set_omega(double omega);

     //! Get omega dihedral angle
     //!
     //! \return Omega dihedral angle          
     double &get_omega();

     //! Get angles as a vector
     //!
     //! \param include_omega Whether to include omega in vector.
     //! \param include_bond_angles Whether to include bond angles in vector
     //! \return Vector of angles
     std::vector<double> get_angles(bool include_omega=false, bool include_bond_angles=false);

     //! Set angle vector
     //!
     //! \param angles Vector of angles
     void set_angles(std::vector<double> angles);

     //! Get neighbouring residue
     //!
     //! \param offset Distance along chain to neighbour (nth neighbour)
     //! \return Residue pointer
     ResidueFB *get_neighbour(int offset);

     //! Add atoms to chain
     //!
     //! \param atom_selection Which atoms to add
     //! \param sidechain_dof_values Optionally set values for sidechain degrees of freedom
     //! \param execute_sort_atoms Whether to sort atoms afterwards
     void add_atoms(definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_O_ATOMS+
                                                                   definitions::BACKBONE_H_ATOMS+
                                                                   definitions::CB_ATOMS),
                    const std::vector<double> &sidechain_dof_values=std::vector<double>(), bool execute_sort_atoms=true);

     //! Add oxygen atoms to chain
     //!
     //! \param atom_selection Specifies which atom to add
     //! \param execute_sort_atoms Whether to sort atoms afterwards
     void add_oxygens(definitions::AtomSelectionEnum atom_selection=definitions::BACKBONE_O_ATOMS,
                      bool execute_sort_atoms=true);

     //! Add hydrogen atoms to chain
     //!
     //! \param atom_selection Specifies which atom to add
     //! \param execute_sort_atoms Whether to sort atoms afterwards
     void add_hydrogens(definitions::AtomSelectionEnum atom_selection=(definitions::BACKBONE_H_ATOMS+
                                                                       definitions::NON_BACKBONE_H_ATOMS),
                        bool execute_sort_atoms=true);

     //! Returns residue name - potentially dependent on residue 
     //! protonation state
     std::string get_protonation_dependent_name();

     //! Output residue as string
     //!
     //! \param o Output stream
     void output(std::ostream &o);

     //! Overload output operator for ResidueFB     
     friend std::ostream & operator<<(std::ostream &o, ResidueFB &r);
};

}

// Implementations
#include "iterators/covalent_bond_iterator.h"

namespace phaistos {

//! Return distance along chain from atom1 to atom2 (in this residue)
//! \param atom1 First atom
//! \param atom1 Second atom
//! \return distance along chain (number of covalent bonds)
template <typename CHAIN_TYPE>
int Residue::chain_distance(Atom *atom1, Atom *atom2) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

     // Check whether atom1 and atom2 are the same
     if (atom1==atom2)
          return 0;
     
     // Make sure that atoms are within the same residue (use chain_distance in chain otherwise)
     assert(atom1->residue == atom2->residue);

     if (atom_chain_distances.size()==0) {
          atom_chain_distances.resize((atoms.size()*(atoms.size()-1)) >> 1, -1);

          // Initialize ALL-N and ALL-C distances
          // these are frequently used by chain::chain_distance
          // and are faster to compute if done in one go
          CovalentBondIterator<CHAIN_TYPE> it1((*this)[N], CovalentBondIterator<CHAIN_TYPE>::INTRA_RESIDUE_ONLY);
          for (; !it1.end(); ++it1) {
               int atom_index1 = (*this)[N]->index;
               int atom_index2 = it1->index;
               if (atom_index1 > atom_index2) {
                    atom_index1 = it1->index;
                    atom_index2 = (*this)[N]->index;
               }
               int cache_index = ((atom_index2*(atom_index2-1)) >> 1) + atom_index1;
               atom_chain_distances[cache_index] = it1.depth;
          }
          CovalentBondIterator<CHAIN_TYPE> it2((*this)[C], CovalentBondIterator<CHAIN_TYPE>::INTRA_RESIDUE_ONLY);
          for (; !it2.end(); ++it2) {
               int atom_index1 = (*this)[C]->index;
               int atom_index2 = it2->index;
               if (atom_index1 > atom_index2) {
                    atom_index1 = it2->index;
                    atom_index2 = (*this)[C]->index;
               }
               int cache_index = ((atom_index2*(atom_index2-1)) >> 1) + atom_index1;
               atom_chain_distances[cache_index] = it2.depth;
          }
     }

     int atom_index1 = atom1->index;
     int atom_index2 = atom2->index;
     if (atom1->index > atom2->index) {
          atom_index1 = atom2->index;
          atom_index2 = atom1->index;
     }
     int cache_index = ((atom_index2*(atom_index2-1)) >> 1) + atom_index1;

     if (atom_chain_distances[cache_index]==-1) {
          // deprecated::CovalentBondIterator it(atom1, deprecated::CovalentBondIterator::INTRA_RESIDUE_ONLY);
          CovalentBondIterator<CHAIN_TYPE> it(atom1, CovalentBondIterator<CHAIN_TYPE>::INTRA_RESIDUE_ONLY);
          for (; !it.end() && &*it != atom2; ++it) {}
          if (&*it == atom2) {
               atom_chain_distances[cache_index] = it.depth;
          } else {
               std::cerr << "Error - chainDistance: (" << atom1->atom_type << "," << atom2->atom_type <<") not found in residue " << atom1->residue->index << ".\n";
               assert(false);
          }
     }
     return atom_chain_distances[cache_index];
}

}

#endif
