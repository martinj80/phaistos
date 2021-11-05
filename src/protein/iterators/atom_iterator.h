// atom_iterator.h --- Iterators over atoms
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

#ifndef ATOM_ITERATOR_H
#define ATOM_ITERATOR_H

#include "iterator_base.h"
#include "protein/atom.h"

namespace phaistos {

// Forward declation
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE, typename ENTITY_TYPE> class AtomIterator;
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE> class AtomPairIterator;


//! Helper class used to implement partial specialization for dereference member function
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE, typename ENTITY_TYPE=Atom>
struct AtomIteratorDereferHelper {
     //! Dereference AtomIterator
     static ENTITY_TYPE &dereference(const AtomIterator<CHAIN_TYPE, ITERATION_MODE, ENTITY_TYPE> &atom_iterator) {
          return *atom_iterator.atom_current;
     }
};

//! Helper class used to implement partial specialization for dereference member function
//! Partial specialization for Vector_3D ENTITY_TYPE.
//! Makes it possible to iterate directly over positions
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE>
struct AtomIteratorDereferHelper<CHAIN_TYPE,ITERATION_MODE,Vector_3D> {
     //! Dereference AtomIterator
     static Vector_3D &dereference(const AtomIterator<CHAIN_TYPE, ITERATION_MODE, Vector_3D> &atom_iterator) {
          return atom_iterator.atom_current->position;
     }
};



//! Iterator over atoms in chain
//! the optional ENTITY_TYPE template argument can be set to Vector_3D to create a 
//! iterator that directly returns atom positions
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE, typename ENTITY_TYPE=Atom>
class AtomIterator: public IteratorBase<ENTITY_TYPE,
                                        AtomIterator<CHAIN_TYPE, ITERATION_MODE, ENTITY_TYPE> > {

     //! Current atom pointer
     Atom *atom_current;

     //! End of iteration pointer
     Atom *atom_end;

     // Friends
     friend class AtomIteratorDereferHelper<CHAIN_TYPE, ITERATION_MODE, ENTITY_TYPE>;
     friend class AtomPairIterator<CHAIN_TYPE, ITERATION_MODE>;
public:

     //! Constructor: from atom reference
     //! note: Runs to end of chain
     //!
     //! \param atom_begin Atom at which to start
     AtomIterator(Atom &atom_begin):
          atom_current(&atom_begin), atom_end(NULL) {}

     //! Constructor: from atom to atom
     //!
     //! \param atom_begin Atom at which to start
     //! \param atom_end Atom at which to end
     AtomIterator(Atom &atom_begin, Atom &atom_end):
          atom_current(&atom_begin), atom_end(&atom_end) {}

     //! Constructor: from atom to atom (pointers)
     //!
     //! \param atom_begin Atom at which to start
     //! \param atom_end Atom at which to end
     AtomIterator(Atom *atom_begin, Atom *atom_end=NULL):
          atom_current(atom_begin), atom_end(atom_end) {}

     //! Constructor: from residue
     //! note: Runs from beginning to end of residue
     //!
     //! \param res Residue to iterate over
     AtomIterator(Residue &res) {

          atom_current = res.get_first_atom(ITERATION_MODE);
          atom_end     = res.get_last_atom(ITERATION_MODE)->template get_neighbour<+1>(ITERATION_MODE);
     }

     //! Constructor: from residue (pointer)
     //! note: Runs from beginning to end of residue
     //!
     //! \param res Residue to iterate over
     AtomIterator(Residue *res) {
          if (res) {
               atom_current = res->get_first_atom(ITERATION_MODE);
               atom_end     = res->get_last_atom(ITERATION_MODE)->template get_neighbour<+1>(ITERATION_MODE);
          } else {
               atom_current = NULL;
               atom_end = NULL;
          }
     }

     //! Constructor: from residue
     //! note: Runs from beginning of residue to end of chain
     //!
     //! \param chain Molecule chain
     //! \param res Residue in which to start
     AtomIterator(const CHAIN_TYPE &chain, Residue &res) {
          atom_current = res.get_first_atom(ITERATION_MODE);
          atom_end     = NULL;
     }

     //! Constructor: from residue (pointer)
     //! note: Runs from beginning to end of chain
     //!
     //! \param chain Molecule chain
     //! \param res Residue in which to start
     AtomIterator(const CHAIN_TYPE &chain, Residue *res) {
          if (res) {
               atom_current = res->get_first_atom(ITERATION_MODE);
          } else {
               atom_current = NULL;
          }
          atom_end = NULL;
     }

     //! Constructor: from residue and atom type
     //! note: Runs from specified atom to end of residue
     //!
     //! \param res_start Residue in which to start
     //! \param atom_type_start Type of atom to start at
     AtomIterator(Residue &res_start, definitions::AtomEnum atom_type_start) {

          if (res_start.has_atom(atom_type_start))
               atom_current = res_start[atom_type_start];
          else
               atom_current = NULL;

          atom_end = res_start.get_last_atom(ITERATION_MODE)->template get_neighbour<+1>(ITERATION_MODE);
     }

     //! Constructor: from residue (pointer) and atom type
     //! note: Runs from specified atom to end of residue
     //!
     //! \param res_start Residue in which to start
     //! \param atom_type_start Type of atom to start at
     AtomIterator(Residue *res_start, definitions::AtomEnum atom_type_start) {

          if (res_start->has_atom(atom_type_start))
               atom_current = (*res_start)[atom_type_start];
          else
               atom_current = NULL;

          atom_end = res_start->get_last_atom(ITERATION_MODE)->template get_neighbour<+1>(ITERATION_MODE);
     }

     //! Constructor: from residue and atom type to residue and atom type
     //!
     //! \param res_start Residue in which to start
     //! \param atom_type_start Type of atom to start at
     //! \param res_end Residue in which to end
     //! \param atom_type_end Type of atom to end at
     AtomIterator(Residue &res_start, definitions::AtomEnum atom_type_start,
                  Residue &res_end,   definitions::AtomEnum atom_type_end) {
          if (res_start.has_atom(atom_type_start))
               atom_current = res_start[atom_type_start];
          if (res_end.has_atom(atom_type_end))
               atom_end = res_end[atom_type_end];
     }

     //! Constructor: from residue (pointer) and atom type to residue (pointer) and atom type
     //!
     //! \param res_start Residue in which to start
     //! \param atom_type_start Type of atom to start at
     //! \param res_end Residue in which to end
     //! \param atom_type_end Type of atom to end at
     AtomIterator(Residue *res_start, definitions::AtomEnum atom_type_start,
                  Residue *res_end,   definitions::AtomEnum atom_type_end) {
          if (res_start->has_atom(atom_type_start))
               atom_current = (*res_start)[atom_type_start];
          if (res_end->has_atom(atom_type_end))
               atom_end = (*res_end)[atom_type_end];
     }

     //! Constructor: from chain
     //!
     //! \param chain Molecule chain
     AtomIterator(const CHAIN_TYPE &chain) {
          atom_current = NULL;
          atom_end = NULL;               
          if (chain.size() > 0)
               atom_current = chain[0].get_first_atom(ITERATION_MODE);
     }

     //! Constructor: from chain, residue index and atom types
     //! note: Runs from specified atom to end of chain
     //!
     //! \param chain Molecule chain
     //! \param res_index_start Residue index at which to start
     //! \param atom_type_start Type of atom to start at
     AtomIterator(const CHAIN_TYPE &chain, 
                  int res_index_start, definitions::AtomEnum atom_type_start) {
          atom_current = NULL;
          atom_end = NULL;               
          if (chain.size() > 0) {
               if (res_index_start >= 0) {
                    atom_current = chain(res_index_start, atom_type_start);
               } else {
                    atom_current = chain(0, atom_type_start);
               }
          }
     }

     //! Constructor: from chain and residue indices
     //!
     //! \param chain Molecule chain
     //! \param res_index_start Residue index at which to start
     //! \param res_index_end Residue index at which to end
     AtomIterator(const CHAIN_TYPE &chain, 
                  int res_index_start, int res_index_end) {

          atom_current = NULL;
          atom_end = NULL;               
          if (chain.size() > 0) {
               if (res_index_start >= 0) {
                    atom_current = chain[res_index_start].get_first_atom(ITERATION_MODE);
               } else {
                    atom_current = chain[0].get_first_atom(ITERATION_MODE);
               }

               if (res_index_end >= 0 && res_index_end < chain.size()) {
                    atom_end = chain[res_index_end].get_first_atom(ITERATION_MODE);
               }
          }
     }

     //! Constructor: from chain, residue indices and atom types
     //!
     //! \param chain Molecule chain
     //! \param res_index_start Residue index at which to start
     //! \param atom_type_start Type of atom to start at
     //! \param res_index_end Residue index at which to end
     //! \param atom_type_end Type of atom to end at
     AtomIterator(const CHAIN_TYPE &chain, 
                  int res_index_start, definitions::AtomEnum atom_type_start,
                  int res_index_end,   definitions::AtomEnum atom_type_end) {

          atom_current = NULL;
          atom_end = NULL;               
          if (chain.size() > 0) {
               if (res_index_start >= 0) {
                    atom_current = chain(res_index_start, atom_type_start);
               } else {
                    atom_current = chain(0, atom_type_start);
               }

               if (res_index_end >= 0 && res_index_end < chain.size()) {
                    atom_end =  chain(res_index_end, atom_type_end);
               }
          }
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     AtomIterator& operator=(const AtomIterator& other) {
          atom_current = other.atom_current;
          atom_end = other.atom_end;
               
          return(*this);
     }

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(const AtomIterator& other) const {
          return(atom_current == other.atom_current);
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     bool operator>(const AtomIterator& other) const {
          if (!atom_current) {
               if (!other.atom_current) {
                    return false;
               } else {
                    return true;
               }
          } else {
               if (!other.atom_current) {
                    return false;
               } else if (atom_current->residue->index != other.atom_current->residue->index) {
                    return (atom_current->residue->index > other.atom_current->residue->index);
               } else {
                    return (atom_current->index > other.atom_current->index);                    
               }
          }
     }

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     bool operator<(const AtomIterator& other) const {
          if (!atom_current) {
               if (!other.atom_current) {
                    return false;
               } else {
                    return false;
               }
          } else {
               if (!other.atom_current) {
                    return true;
               } else if (atom_current->residue->index != other.atom_current->residue->index) {
                    return (atom_current->residue->index < other.atom_current->residue->index);
               } else {
                    return (atom_current->index < other.atom_current->index);
               }
          }
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     AtomIterator &operator++() {
          atom_current = atom_current->get_neighbour<+1>(ITERATION_MODE);
          return (*this);
     }


     //! Increment with value
     //!
     //! \param v Value to increment with
     //! \return Current iterator (this)
     AtomIterator &operator+=(const int v) {
          atom_current = atom_current->get_neighbour(+v, ITERATION_MODE);
          return (*this);
     }

     //! Increment with value - return new iterator
     //!
     //! \param v Value to increment with
     //! \return New iterator
     AtomIterator operator+(const int v) const {
          return AtomIterator<CHAIN_TYPE,ITERATION_MODE,ENTITY_TYPE>(*atom_current->get_neighbour(+v, ITERATION_MODE));
     }

     // THE ITERATORS ARE UNIDIRECTIONAL, BECAUSE IT WOULD OTHERWISE NOT BE POSSIBLE
     // TO DEFINE AUTOMATIC BOUNDS (WITHOUT AN ADDITIONAL CHECK). IF WE NEED A REVERSE
     // ITERATOR, IT SHOULD BE DEFINED SEPARATELY
     // // Decrement
     // AtomIterator &operator--() {
     //      atom_current = atom_current->get_neighbour<-1>(ITERATION_MODE);
     //      return (*this);
     // }

     // // Decrement with value
     // AtomIterator &operator-=(const int v) {
     //      atom_current = atom_current->get_neighbour(-v, ITERATION_MODE);
     //      return (*this);
     // }

     // // Subtraction
     // AtomIterator operator-(const int v) const {
     //      return AtomIterator<CHAIN_TYPE,ITERATION_MODE,ENTITY_TYPE>(atom_current->get_neighbour(-v, ITERATION_MODE));
     // }

     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     // Implemented in different specializations above
     ENTITY_TYPE &operator*() const {
          return AtomIteratorDereferHelper<CHAIN_TYPE,ITERATION_MODE, ENTITY_TYPE>::dereference(*this);
     }

     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     bool end() const {
          return atom_current==atom_end;
     }
};

}

#endif
