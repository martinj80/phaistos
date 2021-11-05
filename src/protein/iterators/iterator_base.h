// iterator_base.h --- Base class for chain iterators
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

#ifndef ITERATOR_BASE_H
#define ITERATOR_BASE_H

namespace phaistos {

//! Base-class/Interface for all chain iterators
//! Note that these iterators DO NOT follow the STL-standard. They are similar
//! to iterators in java - where iterators know about the beginning and end
//! of the container.
template <typename ENTITY_TYPE, typename DERIVED_CLASS>
class IteratorBase {
public:

     //! Local version of Entity template parameter. Used for outside referencing to the type we are iterating over
     typedef ENTITY_TYPE Entity;

     //! Destructor
     virtual ~IteratorBase(){}

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     virtual DERIVED_CLASS& operator=(const DERIVED_CLASS& other) =0;     
     
     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     virtual bool operator==(const DERIVED_CLASS& other) const =0;

     //! Inequality test
     //!
     //! \param other Object to compare with
     //! \return False if objects are identical
     bool operator!=(const DERIVED_CLASS& other) const {
          return (!(*(DERIVED_CLASS*)this == other));
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     virtual bool operator>(const DERIVED_CLASS& other) const =0;

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     virtual bool operator<(const DERIVED_CLASS& other) const =0;

     //! Increment operator
     //!
     //! \return Current iterator (this)
     virtual DERIVED_CLASS &operator++()=0;

     //! Increment with value
     //!
     //! \param v Value to increment with
     //! \return Current iterator (this)
     DERIVED_CLASS &operator+=(const int v) {
          for (int i=0; i<v; i++) {
               ++(*((DERIVED_CLASS*)this));
          }
          return (*((DERIVED_CLASS*)this));          
     }

     //! Increment with value - return new iterator
     //!
     //! \param v Value to increment with
     //! \return New iterator
     DERIVED_CLASS operator+(const int v) {
          DERIVED_CLASS it(*((DERIVED_CLASS*)this));
          it += v;
          return it;
     }

     // THE ITERATORS ARE UNIDIRECTIONAL - THESE OPERATIONS
     // HAVE THEREFORE BEEN REMOVED

     // // Decrement
     // virtual DERIVED_CLASS &operator--()=0;

     // // Decrement with value
     // virtual DERIVED_CLASS &operator-=(const int v)=0;

     // // Subtraction
     // virtual DERIVED_CLASS operator-(const int v) const=0;

     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     virtual ENTITY_TYPE &operator*() const =0;
     
     //! Dereference operator (*)
     //!
     //! \return pointer to underlying entity type
     ENTITY_TYPE *operator->() const {
          return &(*(*((DERIVED_CLASS*)this)));
     }

     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     virtual bool end() const =0;

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, const IteratorBase<ENTITY_TYPE, DERIVED_CLASS> &it) {
	  o << *it;
	  return o;
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
          DERIVED_CLASS it(*(DERIVED_CLASS*)this);
          for (; !it.end(); ++it) {
               func(*it);
          }
          return func;
     }

};


//! Cached iterator base class
//!
//! \tparam CONTRIBUTION_TYPE Type of value added in each step
//! \tparam RETURN_VALUE_TYPE Type of value returned by cached iterator (total value)
template <typename CONTRIBUTION_TYPE, typename RETURN_VALUE_TYPE>
class CachedIteratorBase {
protected: 

     //! Time stamp for testing whether cache is in sync
     long int time_stamp;

     //! Constructor (default)
     CachedIteratorBase():time_stamp(0){}

     //! Destructor
     virtual ~CachedIteratorBase(){}

     //! Register a contribution to the cached iterator
     //!
     //! \param contribution Value to register in cache
     //! \return reference to cache entry in which value is stored
     virtual RETURN_VALUE_TYPE &register_contribution(const CONTRIBUTION_TYPE &contribution=0)=0;

     //! Compute total value for current evaluation
     //!
     //! \return Total value
     virtual RETURN_VALUE_TYPE &compute_total()=0;

     //! Accept last evaluation
     virtual void accept()=0;

     //! Reject last evaluation
     virtual void reject()=0;

     //! Reset chaintree if cache is not in sync with chain
     //!
     //! \param chain Molecule chain
     //! \return True if cache was in sync
     template <typename CHAIN_TYPE>
     bool enforce_cache_sync(CHAIN_TYPE &chain) {
          if (this->time_stamp != chain.time_stamp) {
               std::cout << "Cache out of sync with chain (" << this->time_stamp << "," << chain.time_stamp << "). Resetting to " << chain.time_stamp << ".\n";
               // chain.time_stamp = this->time_stamp;
               this->time_stamp = chain.time_stamp;
               return true;
          } else {
               return false;
          }
     }
};

//! Cached iterator.
//! For each iterator, a specialization of this class should be implemented
template <typename ITERATOR_TYPE, typename CONTRIBUTION_TYPE=double, typename RETURN_VALUE_TYPE=double>
class CachedIterator: public CachedIteratorBase<CONTRIBUTION_TYPE, RETURN_VALUE_TYPE>, 
                      public ITERATOR_TYPE { 
};

}
#endif
