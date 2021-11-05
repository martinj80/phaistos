// parameter_data_structures.h --- Datastructures for storing energy parameters
// Copyright (C) 2008-2011 Kristoffer Enøe Johansson, Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef PARAMETERDATASTRUCTURES_H
#define PARAMETERDATASTRUCTURES_H

#include <map>
#include <iostream>
#include <algorithm>

namespace phaistos {

//! 1D parameter lookup table
template <typename VAL_TYPE>
class ParameterData1D {

     //! int -> VAL_TYPE map data type
     std::map<int, VAL_TYPE> data;

public:

     //! Lookup operator
     VAL_TYPE &operator()(int atom_type) {
	  return data[atom_type];
     }     
};


//! 2D parameter lookup table
template <typename VAL_TYPE>
class ParameterData2D {

private:

     //! Internal class for manipulation two integers
     class Int2 {

	  //! Internal data representation
	  int data[2];
	  
     public:

	  //! Constructor
	  Int2(int atom1, int atom2) {
	       data[0] = atom1;
	       data[1] = atom2;
	  }

	  //! Equality operator
	  friend bool operator==(const Int2 &first, const Int2 &second) {
	       return ((first.data[0] == second.data[0]) &&
		       (first.data[1] == second.data[1]));	       
	  }

	  //! Smaller than operator
	  friend bool operator<(const Int2 &first, const Int2 &second) {
	       return std::lexicographical_compare(first.data, first.data+2,
						   second.data, second.data+2);
	  }
     };
     
     
     //! (int,int) -> VAL_TYPE map data type
     std::map<Int2, VAL_TYPE> data;

public:

     //! Lookup operator
     VAL_TYPE &operator()(int atom_type1, int atom_type2) {
	  return data[Int2(atom_type1, atom_type2)];
     }     
};


//! 3D parameter lookup table
template <typename VAL_TYPE>
class ParameterData3D {

private:

     //! Internal class for manipulation three integers
     class Int3 {

	  //! Internal data representation
	  int data[3];
	  
     public:

	  //! Constructor
	  Int3(int atom1, int atom2, int atom3) {
	       data[0] = atom1;
	       data[1] = atom2;
	       data[2] = atom3;
	  }

	  //! Equality operator
	  friend bool operator==(const Int3 &first, const Int3 &second) {
	       return ((first.data[0] == second.data[0]) &&
		       (first.data[1] == second.data[1]) &&
		       (first.data[2] == second.data[2]));	       
	  }

	  //! Smaller than operator
	  friend bool operator<(const Int3 &first, const Int3 &second) {
	       return std::lexicographical_compare(first.data, first.data+3,
						   second.data, second.data+3);
	  }
     };
     
     
     //! (int,int) -> VAL_TYPE map data type
     std::map<Int3, VAL_TYPE> data;

public:

     //! Lookup operator
     VAL_TYPE &operator()(int atom_type1, int atom_type2, int atom_type3) {
	  return data[Int3(atom_type1, atom_type2, atom_type3)];
     }     
};

//! 4D parameter lookup table
template <typename VAL_TYPE>
class ParameterData4D {
private:

     //! Internal class for manipulation four integers
     class Int4 {

	  //! Internal data representation
	  int data[4];
	  
     public:

	  //! Constructor
	  Int4(int atom1, int atom2, int atom3, int atom4) {
	       data[0] = atom1;
	       data[1] = atom2;
	       data[2] = atom3;
	       data[3] = atom4;
	  }

	  //! Equality operator
	  friend bool operator==(const Int4 &first, const Int4 &second) {
	       return ((first.data[0] == second.data[0]) &&
		       (first.data[1] == second.data[1]) &&
		       (first.data[2] == second.data[2]) &&
		       (first.data[3] == second.data[3]));	       
	  }

	  //! Smaller than operator
	  friend bool operator<(const Int4 &first, const Int4 &second) {
	       return std::lexicographical_compare(first.data, first.data+4,
						   second.data, second.data+4);
	  }
     };
     
     
     //! (int,int) -> VAL_TYPE map data type
     std::map<Int4, VAL_TYPE> data;

public:

     //! Lookup operator
     VAL_TYPE &operator()(int atom_type1, int atom_type2, int atom_type3, int atom_type4) {
	  return data[Int4(atom_type1, atom_type2, atom_type3, atom_type4)];
     }     
};

}

#endif     
