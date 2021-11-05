// atom_pair_iterator.h --- Iterators over atom pairs
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

#ifndef ATOMPAIR_ITERATOR_H
#define ATOMPAIR_ITERATOR_H

#include "iterator_base.h"

namespace phaistos {

//! Iterator over atom-pairs in chain. ITERATION_MODE can be 
//! ALL, BACKBONE, SC_ONLY or CA_ONLY
template <typename CHAIN_TYPE, definitions::IterateEnum ITERATION_MODE>
class AtomPairIterator: public IteratorBase<const std::pair<Atom*,Atom*>,
                                            AtomPairIterator<CHAIN_TYPE,ITERATION_MODE> > {

     //@{
     //! Inner atom iterators
     AtomIterator<CHAIN_TYPE,ITERATION_MODE> it1;
     AtomIterator<CHAIN_TYPE,ITERATION_MODE> it2;
     //@}

public:

     //! This pair constains references to 
     //! the current atoms in it1 and it2
     std::pair<Atom*, Atom*> atom_pair;

     //! Constructor: from atom pointers 
     //!
     //! \param atom_start1 Atom at which to start first element of pair
     //! \param atom_start2 Atom at which to start second element of pair
     //! \param atom_end1 Atom at which to end first element of pair
     //! \param atom_end2 Atom at which to end second element of pair
     AtomPairIterator(Atom *atom_start1, Atom *atom_start2,
                      Atom *atom_end1=NULL, Atom *atom_end2=NULL)
          : it1(atom_start1, atom_end1), 
            it2(atom_start2, atom_end2),
            atom_pair(it1.atom_current,
                      it2.atom_current) {

     }

     //! Constructor: from chain
     //!
     //! \param chain Molecule chain     
     AtomPairIterator(const CHAIN_TYPE &chain)
          : it1(chain), 
            it2(chain),
            atom_pair(it1.atom_current,
                      it2.atom_current) {

     }

     //! Copy constructor.
     //!
     //! \param other Source object from which copy is made.
     AtomPairIterator(const AtomPairIterator& other)
          : it1(other.it1),
            it2(other.it2),
            atom_pair(it1.atom_current,
                      it2.atom_current){}


     //! Assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     AtomPairIterator& operator=(const AtomPairIterator& other) {
          it1 = other.it1;
          it2 = other.it2;
          atom_pair = other.atom_pair;
          return (*this);
     }     
     
     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(const AtomPairIterator& other) const {
          return (it1 == other.it1 && 
                  it2 == other.it2);
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     bool operator>(const AtomPairIterator& other) const {
          if (it1 == other.it1)
               return (it2 > other.it2);
          else
               return it1 > other.it1;
     }

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     bool operator<(const AtomPairIterator& other) const {
          if (it1 == other.it1)
               return (it2 < other.it2);
          else
               return it1 < other.it1;
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     AtomPairIterator &operator++() {
          ++it2;
          atom_pair.second = it2.atom_current;
          if (it2.end()) {
               ++it1;
               atom_pair.first = it1.atom_current;
               it2 = it1;
               atom_pair.second = atom_pair.first;
          }
          return (*this);
     }


     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     const std::pair<Atom*,Atom*> &operator*() const {
          return atom_pair;
     }


     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     bool end() const {
          return it2.end();
     }

     //! Foreach member function that will apply a function to all entities
     //! This is often slightly faster than iterating manually
     //!
     //! \param func Functor object
     template <typename CALLABLE_TYPE>
     CALLABLE_TYPE &for_each(CALLABLE_TYPE &func) {
          // Make local copy of it2 for extra speed 
          // (lookup of local variables is slightly faster than lookup
          //  of member variables)
          AtomIterator<CHAIN_TYPE,ALL> it2(this->it2);
          std::pair<Atom*, Atom*> atom_pair(it1.atom_current,
                                            it2.atom_current);
          for (; !it1.end(); ++it1) {
               atom_pair.first = it1.atom_current;
               it2.atom_current = it1.atom_current;               
               for (++it2; !it2.end(); ++it2) {
                    atom_pair.second = it2.atom_current;
                    func(atom_pair);
               }
          }
          return func;
     }

};

}

#endif
