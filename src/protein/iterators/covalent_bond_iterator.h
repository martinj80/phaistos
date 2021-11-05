// covalent_bond_iterator.h --- Iterators over covalently bonded atoms
// Copyright (C) 2006-2010 Wouter Boomsma
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

#ifndef COVALENTBOND_ITERATOR_H
#define COVALENTBOND_ITERATOR_H

#include "iterator_base.h"
#include "protein/atom.h"

namespace phaistos {

//! Iterate recursively over all covalent bond neighbours of an atom 
template <typename CHAIN_TYPE>
class CovalentBondIterator: public IteratorBase<Atom, CovalentBondIterator<CHAIN_TYPE> > {
protected:

     //! Queue of atom,depth pairs used for breadth-first search
     std::queue<std::pair<Atom *, int> > atom_queue;

     //! Vector containing already visited atoms
     std::vector<std::vector<bool> > *visited_atoms;

     //! Indicator determining whether the iterator should free
     //! the vector of visited atoms itself
     bool visited_atoms_owner;

     //! Indicates whether iteration should be limited to current residue
     bool intra_residue_only;
     
public:
     //! Mode of iteration
     enum ModeEnum{ALL_BONDS=0, INTRA_RESIDUE_ONLY, DEPTH_1_ONLY} iteration_mode; //!< Mode of iteration

     //! Current atom     
     Atom *atom_current;

     //! Current depth of search     
     int depth;
     
     //! Constructor
     //!
     //! \param atom_begin Atom at which to start
     //! \param iteration_mode Type of iteration (ALL_BONDS, INTRA_RESIDUE_ONLY, DEPTH_1_ONLY)
     explicit CovalentBondIterator(Atom *atom_begin, 
                                   ModeEnum iteration_mode=ALL_BONDS)
          : visited_atoms(NULL), visited_atoms_owner(false),
            iteration_mode(iteration_mode), 
            atom_current(atom_begin), depth(0) {

          // Add atom to queue 
          atom_queue.push(std::pair<Atom *, int>(atom_begin,0));

          if (iteration_mode != DEPTH_1_ONLY) {

               init();

               // Add atom to visited_atoms vector
               int residue_index = atom_begin->residue->index;
               if (iteration_mode==INTRA_RESIDUE_ONLY)
                    residue_index = 0;
               (*visited_atoms)[residue_index].resize(atom_begin->index+1, false);
               (*visited_atoms)[residue_index][atom_begin->index] = true;
          }

	  // Start by making a single iteration-step 
	  ++(*this);
     }

     //! Constructor.
     //!
     //! \param atom_begin Atom at which to start
     //! \param forbidden_neighbours Specifies a vector of forbidden 
     //!        neighbours which will not be included in the iteration
     //!        This is useful for providing a specific direction.
     //!        For instance a=CB, forbiddenNeighbour=CA, will only iterate
     //!        in the side chain direction.
     //! \param iteration_mode Mode of iteration (ALL_BONDS, INTRA_RESIDUE_ONLY, DEPTH_1_ONLY)
     CovalentBondIterator(Atom *atom_begin, 
                          std::vector<Atom *> forbidden_neighbours,
                          ModeEnum iteration_mode=ALL_BONDS):
	  visited_atoms(NULL), visited_atoms_owner(false),
	  iteration_mode(iteration_mode), 
          atom_current(atom_begin), depth(0) {

          // Add atom to queue 
          atom_queue.push(std::pair<Atom *, int>(atom_begin,0));

          if (iteration_mode != DEPTH_1_ONLY) {

               init();

               // Add atom to visitedAtoms vector
               int residue_index = atom_begin->residue->index;
               if (iteration_mode==INTRA_RESIDUE_ONLY)
                    residue_index = 0;
               (*visited_atoms)[residue_index].resize(atom_begin->index+1, false);
               (*visited_atoms)[residue_index][atom_begin->index] = true;

               // Put forbidden neighbours in visited list
               for (unsigned int i=0; i<forbidden_neighbours.size(); i++) {
                    (*visited_atoms)[residue_index][forbidden_neighbours[i]->index] = true;
               }
          }

          
	  // Start by making a single iteration-step 
	  // ++(*this);
     }     
     
     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     CovalentBondIterator(const CovalentBondIterator &other) {
	  this->atom_queue = other.atom_queue;
	  this->atom_current = other.atom_current;
	  this->depth = other.depth;
	  this->iteration_mode = other.iteration_mode;

	  // Copy visited_atoms pointer - original pointer remains owner
	  this->visited_atoms = other.visited_atoms;
	  this->visited_atoms_owner = false;
     }

     //! Assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     CovalentBondIterator& operator=(const CovalentBondIterator& other) {
	  this->atom_queue = other.atom_queue;
	  this->atom_current = other.atom_current;
	  this->depth = other.depth;
	  this->iteration_mode = other.iteration_mode;

	  // Copy visited_atoms pointer - original pointer remains owner
	  this->visited_atoms = other.visited_atoms;
	  this->visited_atoms_owner = false;
          return (*this);
     }     
     
     //! Destructor
     ~CovalentBondIterator() {
	  if (visited_atoms_owner)
	       delete visited_atoms;
     }

     //! Initializer
     void init() {
	  if (visited_atoms == NULL) {
	       int size = atom_current->residue->index+1;
	       if (iteration_mode==INTRA_RESIDUE_ONLY)
	            size = 1;
	       visited_atoms = new std::vector<std::vector<bool> >(size);
	       visited_atoms_owner = true;
	  }
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     CovalentBondIterator& operator++() {
	  if (!atom_queue.empty()) {

	       // Pop front element of queue
	       std::pair<Atom *, int> front = atom_queue.front(); atom_queue.pop();
	       Atom *front_atom = front.first;
	       int front_depth = front.second;

	       // Iterate over all covalent neighbours of current atom
	       for (unsigned int i=0; i< front_atom->covalent_neighbours.size(); i++) {
		    definitions::AtomEnum neighbourAtomType = front_atom->covalent_neighbours[i].first;
		    int residue_offset = front_atom->covalent_neighbours[i].second;

		    // Disregard neighbours involving other residues if intraResidueOnly flag is set
		    if (residue_offset != 0 && (iteration_mode==INTRA_RESIDUE_ONLY))
			 continue;
		    
		    Residue *residue = front_atom->residue->get_neighbour(residue_offset);
		    if (residue->has_atom(neighbourAtomType)) {
			 Atom *neighbour = (*residue)[neighbourAtomType];
			 int neighbour_atom_index = neighbour->index;
			 int neighbour_residue_index = residue->index;

                         if (iteration_mode==DEPTH_1_ONLY) {
                              if (front_depth == 0)
                                   atom_queue.push(std::pair<Atom *, int>(neighbour, front_depth+1));
                         } else {
                         
                              if (iteration_mode==INTRA_RESIDUE_ONLY)
                                   neighbour_residue_index = 0;

                              // Resize vector if necessary (residue dimension)
                              if (neighbour_residue_index >= (int)visited_atoms->size())
                                   visited_atoms->resize(neighbour_residue_index+1);
                         
                              // Resize vector if necessary (atom dimension)
                              if (neighbour_atom_index >= (int)((*visited_atoms)[neighbour_residue_index].size()))
                                   (*visited_atoms)[neighbour_residue_index].resize(neighbour_atom_index+1, false);
                         
                              if (!(*visited_atoms)[neighbour_residue_index][neighbour_atom_index]) {
                                   atom_queue.push(std::pair<Atom *, int>(neighbour, front_depth+1));
                                   (*visited_atoms)[neighbour_residue_index][neighbour_atom_index] = true;
                              }
                         }
		    }	       
	       }
	  }

	  if (!atom_queue.empty()) {
	       atom_current = atom_queue.front().first;
	       depth = atom_queue.front().second;
	  } else {
	       atom_current = NULL;
	  }
	  return (*this);
     }

     //! Increment with value
     //!
     //! \param v Value to increment with
     //! \return Current iterator (this)
     CovalentBondIterator &operator+=(const int v) {
          for (int i=0; i<v; i++) {
               ++(*this);
          }
          return (*this);
     }

     //! Increment with value - return new iterator
     //!
     //! \param v Value to increment with
     //! \return New iterator
     CovalentBondIterator operator+(const int v) const {
          CovalentBondIterator it(*this);
          it += v;
          return it;
     }

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(const CovalentBondIterator& other) const {
          return (atom_current == other.atom_current);
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     bool operator>(const CovalentBondIterator& other) const {
          return depth > other.depth;
     }

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     bool operator<(const CovalentBondIterator& other) const {
          return depth < other.depth;
     }

     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     // Implemented in different specializations above
     Atom &operator*() const {
	  return *atom_current;
     }
     
     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     bool end() const {
          return atom_current == NULL;
     }

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, const CovalentBondIterator &it) {
	  o << it.atom_current << " " << it.depth;
	  return o;
     }

};

}

#endif
