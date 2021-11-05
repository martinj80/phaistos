// residue_iterator.h --- Iterators over residues
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

#ifndef RESIDUE_ITERATOR_H
#define RESIDUE_ITERATOR_H

#include "iterator_base.h"

namespace phaistos {

//! Iterator over residues
template <typename CHAIN_TYPE>
class ResidueIterator: public IteratorBase<typename CHAIN_TYPE::Residue,
                                           ResidueIterator<CHAIN_TYPE> > {
protected:

     //! Local Residue container
     const std::vector<typename CHAIN_TYPE::Residue *> &residues;

     //! Current residue index
     int residue_index_current;

     //! End of iteration pointer     
     int residue_index_end;

     //! Constructor: For internal purposes. Used when copying.
     //!
     //! \param residues Local Residue Container
     //! \param residue_index_start Residue index at which to start
     //! \param residue_index_end Residue index at which to end
     ResidueIterator(const std::vector<typename CHAIN_TYPE::Residue *> &residues,
                     int residue_index_start=0, int residue_index_end=-1)
          : residues(residues),
            residue_index_current(residue_index_start),
            residue_index_end(residue_index_end) {

          if (residue_index_end<0)
               this->residue_index_end = residues.size();

     }

public:

     //! Constructor: from chain, residue index
     //!
     //! \param chain Molecule chain
     //! \param residue_index_start Residue index at which to start
     //! \param residue_index_end Residue index at which to end
     // Note: this function has external requirement not available at
     // this point in the code. The definition is therefore included at the
     // end of the file.
     ResidueIterator(const CHAIN_TYPE &chain,
                     int residue_index_start=0, int residue_index_end=-1);

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     ResidueIterator& operator=(const ResidueIterator& other) {
          residue_index_current = other.residue_index_current;
          residue_index_end = other.residue_index_end;
          return(*this);
     }

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(const ResidueIterator& other) const {
          return(residue_index_current == other.residue_index_current);
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     bool operator>(const ResidueIterator& other) const {
          return (residue_index_current > other.residue_index_current);
     }

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     bool operator<(const ResidueIterator& other) const {
          return (residue_index_current < other.residue_index_current);
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     ResidueIterator &operator++() {
          residue_index_current++;
          return (*this);
     }


     //! Increment with value
     //!
     //! \param v Value to increment with
     //! \return Current iterator (this)
     ResidueIterator &operator+=(const int v) {
          residue_index_current+=v;
          return (*this);
     }

     //! Increment with value - return new iterator
     //!
     //! \param v Value to increment with
     //! \return New iterator
     ResidueIterator operator+(const int v) const {
          return ResidueIterator(residues,
                                 residue_index_current+v,
                                 residue_index_end);
     }

     // THE ITERATORS ARE UNIDIRECTIONAL - THESE OPERATIONS
     // HAVE THEREFORE BEEN REMOVED

     // // Decrement
     // ResidueIterator &operator--() {
     //      residue_index_current--;
     //      return (*this);
     // }


     // // Decrement with value
     // ResidueIterator &operator-=(const int v) {
     //      residue_index_current-=v;
     //      return (*this);
     // }

     // // Subtraction
     // ResidueIterator operator-(const int v) const {
     //      return ResidueIterator(chain,
     //                             residue_index_current-v,
     //                             residue_index_end);
     // }

     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     typename CHAIN_TYPE::Residue &operator*() const {
          return *residues[residue_index_current];
     }

     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     bool end() const {
          return (residue_index_current>=residue_index_end);
     }
};

}

// The code below requires complete types of chains

#include "protein/chain_fb.h"
#include "protein/chain_ca.h"

namespace phaistos {

//! Constructor: from chain, residue index
//!
//! \param chain Molecule chain
//! \param residue_index_start Residue index at which to start
//! \param residue_index_end Residue index at which to end
template <typename CHAIN_TYPE>
ResidueIterator<CHAIN_TYPE>::ResidueIterator(const CHAIN_TYPE &chain,
                                               int residue_index_start,
                                               int residue_index_end)
     : residues(chain.residues),
       residue_index_current(residue_index_start),
       residue_index_end(residue_index_end) {

     if (residue_index_end<0) {
          this->residue_index_end = chain.size();
     }
}





//! CachedIterator specialization for ResidueIterator
template <typename CHAIN_TYPE, typename RETURN_VALUE_TYPE>
class CachedIterator<ResidueIterator<CHAIN_TYPE>, RETURN_VALUE_TYPE, RETURN_VALUE_TYPE>
     : public ResidueIterator<CHAIN_TYPE>,
       public CachedIteratorBase<RETURN_VALUE_TYPE,RETURN_VALUE_TYPE> {

     //! Cache for each residue
     std::vector<RETURN_VALUE_TYPE> cache;

     //! Vector of (index,value) pairs containing backup values
     std::vector<std::pair<int,RETURN_VALUE_TYPE> > cache_backup;

     //! Total sum value
     RETURN_VALUE_TYPE sum;

     //! Flag indicating whether cache has been initialized
     bool initialized;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     CachedIterator(CHAIN_TYPE &chain)
          : ResidueIterator<CHAIN_TYPE>(chain),
            cache(chain.size(), 0),
            sum(UNINITIALIZED),
            initialized(false){
     }

     //! Overload () operator to act as a constructor (has the same interface as the ResidueIterator constructor)
     //!
     //! \param chain Molecule chain
     //! \param residue_index_start Residue index at which to start
     //! \param residue_index_end Residue index at which to end
     //! \param residue_index_full_range_start Optionally specify range to use when reinitializing cache
     //! \param residue_index_full_range_end Optionally specify range to use when reinitializing cache
     void operator()(CHAIN_TYPE &chain,
                     int residue_index_start,
                     int residue_index_end,
                     int residue_index_full_range_start=0,
                     int residue_index_full_range_end=-1) {

          if (residue_index_full_range_end < 0)
               residue_index_full_range_end = chain.size();

          // Ensure that cache is in sync with chain;
          if (this->enforce_cache_sync(chain))
               initialized=false;

          if (initialized) {
               ResidueIterator<CHAIN_TYPE>::operator=(ResidueIterator<CHAIN_TYPE>(chain,
                                                                                  residue_index_start,
                                                                                  residue_index_end));

          } else {

               sum = 0;
               cache = std::vector<RETURN_VALUE_TYPE>(chain.size(), 0.0);

               // If not initialized, do a full iteration
               ResidueIterator<CHAIN_TYPE>::operator=(ResidueIterator<CHAIN_TYPE>(chain,
                                                                                  residue_index_full_range_start,
                                                                                  residue_index_full_range_end));
          }
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     CachedIterator &operator=(const ResidueIterator<CHAIN_TYPE> &other) {
          ResidueIterator<CHAIN_TYPE>::operator=(other);
          return *this;
     }

     //! Register a contribution to the cached iterator
     //!
     //! \param contribution Value to register
     //! \return reference to cache entry in which value is stored
     RETURN_VALUE_TYPE &register_contribution(const RETURN_VALUE_TYPE &contribution=0) {

          // Save old contribution to backup
          cache_backup.push_back(std::make_pair(this->residue_index_current, cache[this->residue_index_current]));

          // Subtract old contribution from sum
          sum -= cache[this->residue_index_current];

          // Add new contribution to sum
          sum += contribution;

          // Write new contribution to cache
          cache[this->residue_index_current] = contribution;

          return cache[this->residue_index_current];
     }

     //! Compute total value
     //!
     //! \return Total sum
     RETURN_VALUE_TYPE &compute_total() {
          return sum;
     }

     //! Accept last evaluation
     void accept() {
          if (is_initialized(sum)) {
               initialized = true;
               cache_backup.clear();
               this->time_stamp++;
          }
     }

     //! Reject last evaluation
     void reject() {
          // Iterate over backup to restore cache
          for (unsigned int i=0; i<cache_backup.size(); ++i) {

               // Subtract current value from sum
               sum -= cache[cache_backup[i].first];

               // Add value from backup to sum
               sum += cache_backup[i].second;

               // Replace value stored in cache
               cache[cache_backup[i].first] = cache_backup[i].second;
          }

          // Clear backup
          cache_backup.clear();
     }
};

}

#endif
