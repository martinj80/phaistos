// git_element.cpp --- Encapsulated GIT vector class for the clustering
// Copyright (C) 2011  Tim Harder, Thomas Hamelryck, Wouter Boomsma
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
#include <vector>

#include "git.h"
#include "git_element.h"

#include "utils/vector_matrix_3d.h"
#include "protein/chain_ca.h"
#include "protein/chain_fb.h"

#include <stdio.h>
#include <string>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

namespace phaistos {


//! Constructor
GitElement::GitElement(std::string &name, int &chain, int &index, int &domain, int &length, std::vector<double> &git) {
     this->name = name;
     this->chain = chain;
     this->index = index;
     this->domain = domain;
     this->length = length;
     this->git = git;
     this->distance_to_mean = -1;
     this->weight = 1;
}

///////////////////////////////////////////////////
//
// Getters
//
double GitElement::get_weight() {
     return this->weight;
}
std::string GitElement::get_name() {
     return this->name;
}
int GitElement::get_chain() {
     return this->chain;
}
int GitElement::get_index() {
     return this->index;
}
int GitElement::get_domain() {
     return this->domain;
}
int GitElement::get_length() {
     return this->length;
}
std::vector<double> GitElement::get_git() {
     return this->git;
}
double GitElement::get_git(int index) {
     return this->git[index];
}
double GitElement::get_distance_to_mean() const {
     return this->distance_to_mean;
}

///////////////////////////////////////////////////
//
// Setters
//
void GitElement::set_weight(double weight) {
     this->weight = weight;
}
void GitElement::set_name(std::string name) {
     this->name = name;
}
void GitElement::set_name(const char *name) {
     this->name = name;
}
void GitElement::set_domain(int domain) {
     this->domain = domain;
}
void GitElement::set_chain(int chain) {
     this->chain = chain;
}
void GitElement::set_index(int index) {
     this->index = index;
}
void GitElement::set_length(int length) {
     this->length = length;
}
void GitElement::set_git(double(&git_arr)[31]) {
     this->git.clear();
     for (int i = 0; i < 31; i++) {
          this->git.push_back(git_arr[i]);
     }
}
void GitElement::set_git(std::vector<double> git) {
     this->git.clear();
     this->git = git;
}

void GitElement::set_distance_to_mean(double distance) {
     this->distance_to_mean = distance;
}
void GitElement::set_distance_to_mean(GitElement &mean) {
     this->distance_to_mean = get_git_distance(this->git, mean.get_git());
}


// Output vector to the standard output stream
void GitElement::print() {
     for (unsigned int i = 0; i < this->git.size(); i++) {
          printf("%.4f ", this->git[i]);
     }
     printf("\n");
}


}
