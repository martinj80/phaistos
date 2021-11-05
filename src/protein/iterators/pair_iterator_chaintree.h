// pair_iterator_chaintree.h --- Iterators over pairs using chain tree
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

#ifndef PAIR_ITERATOR_CHAINTREE_H
#define PAIR_ITERATOR_CHAINTREE_H

#include "protein/chaintree.h"
#include "iterator_base.h"

namespace phaistos {

// NOTE: This iterator is defined in chaintree namespace
namespace chaintree {


// Forward declarations
template <typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2> class PairIteratorOneVsAll;
template <typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2> class PairIterator;
template <typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2> class PairIteratorAtomBase;
template <typename DERIVED_CLASS, typename CHAIN_TYPE, 
          typename ENTITY_TYPE1, typename ENTITY_TYPE2> class PairIteratorOneVsAllBase;

//! Return type used by chaintree pair iterators: a entity pair and a distance
template <typename TYPE1, typename TYPE2>
class PairWithDistance: public std::pair<TYPE1,TYPE2> {
public:
     //! Distance value
     double distance;

     //! Constructor
     //!
     //! \param val1 Value for element 1
     //! \param val2 Value for element 2
     //! \param distance between element 1 and 2
     PairWithDistance(const TYPE1 &val1,
                      const TYPE2 &val2,
                      const double distance):
          std::pair<TYPE1,TYPE2>(val1,val2),
          distance(distance){}

     //! Overload output operator
     friend std::ostream &operator<<(std::ostream &o, PairWithDistance<TYPE1,TYPE2> pair) {
          o << "(" << pair.first << ", " << pair.second << ", " << pair.distance << ")";
          return o;
     }
};


//! Base class for all chaintree pair iterators
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename ENTITY_TYPE1=Atom, typename ENTITY_TYPE2=Atom>
class PairIteratorBase: public IteratorBase<const PairWithDistance<ENTITY_TYPE1*,ENTITY_TYPE2*>,
                                             DERIVED_CLASS> {

protected:
     
     // Friends
     friend class PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>;
     friend class PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,ENTITY_TYPE1>;

     //! Local settings class.
     //! It is useful to declare your settings once and for all, instead
     //! of on every construction of the iterator, because some computations
     //! are involved.
     class Settings {
     public:
          //! Maximum interaction distance between two bounding volumes
          //! this is automatically set to max(maxInteractionDistances)
          double max_interaction_distance_bv;

          //! The maximum interaction distance for atoms not specified explictly in
          //! selected_atoms and maxInteractionDistances
          //! Note: setting this value forces an iteration over all atoms in a node
          double default_max_interaction_distance;
          
          //! Whether only to report pairs that were modified in the last move
          bool only_modified_pairs;

          //! How far the residues should be separated in the chain
          //! before pairs are considered. The Default is 0, meaning all
          //! pairs will be considered.
          //! NOTE that a "residue" here covers the N,H,CA,HA* and sidechain atoms
          //! of the current residue PLUS the C,O of the previous atom (peptide
          //! plane), and EXCLUDES the C,O of the current residue.
          int minimum_residue_distance;

          //! Whether to omit the most fine grained level of the iteration
          bool omit_lowest_level;

          //! Constructor
          //!
          //! \param max_interaction_distance Maximum interaction distance between two bounding volumes
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration
          Settings(double max_interaction_distance=std::numeric_limits<double>::infinity(),
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false):
               max_interaction_distance_bv(max_interaction_distance),
               default_max_interaction_distance(max_interaction_distance),
               only_modified_pairs(only_modified_pairs),
               minimum_residue_distance(minimum_residue_distance),
               omit_lowest_level(omit_lowest_level) {

               // Find max_interaction_distance_bv: max(max_interaction_distance_matrix)
               this->max_interaction_distance_bv = default_max_interaction_distance;                    
          }

     } settings;    //!< Local settings object

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Node corresponding to current first element
     Node<Residue, BvType> *node1;

     //! Node corresponding to current second element
     Node<Residue, BvType> *node2;

     //! A pair with references to the entity pointers (used by operator*)
     PairWithDistance<ENTITY_TYPE1*, ENTITY_TYPE2*> entity_pair;

     //! Currently active entity 1 (Atom|BvType|NodeType)
     ENTITY_TYPE1 *&entity1;

     //! Currently active entity 2 (Atom|BvType|NodeType)
     ENTITY_TYPE2 *&entity2;

     //! Distance between entity 1 and 2
     double &distance;

     //! Stack of nodes to be processed
     //! a NULL pointer on this stack denotes that a nodePair should be popped
     std::stack<Node<typename CHAIN_TYPE::Residue, BvType> *> node_stack;

     //! Stack of node-pairs to be processed
     std::stack<Node<typename CHAIN_TYPE::Residue, BvType> *> node_pair_stack;

     //! Stack of rotation matrices relating the frames of reference
     //! of the two nodes popped from node_pair_stack
     std::stack<Matrix_3D> rotation_stack;

     //! Stack of translation vectors relating the frames of reference
     //! of the two nodes popped from node_pair_stack
     std::stack<Vector_3D> translation_stack;

     //! Stack of booleans specifiying for each popped pair whether
     //! the nodes are "separated". This is used to allow iteration
     //! over only those pairs that were modified in last move 
     std::stack<bool> separation_stack;

     //! Current rotation
     Matrix_3D rotation;

     //! Current translation
     Vector_3D translation;

     //! Current separation
     bool separation;

     //! The associated chain tree 
     typename CHAIN_TYPE::ChainTree *ct;
          
public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIteratorBase(CHAIN_TYPE &chain,
                      const Settings &settings)
          : settings(settings),
            node1(NULL), node2(NULL),
            entity_pair(NULL, NULL, UNINITIALIZED), 
            entity1(entity_pair.first), entity2(entity_pair.second), distance(entity_pair.distance),
            ct(chain.get_chain_tree()) {
     }

     //! Initializer
     void init() {

          // Initialize node stack at root node
          node_stack.push(ct->nodes[ct->nodes.size()-1]);

          // Start at first pair
          ++(*((DERIVED_CLASS*)this));
     }

     //! Determine whether pair of nodes are separated by more
     //! than the maximum interaction distance
     //!
     //! \param node1 First node
     //! \param node2 Second node
     //! \param translation Translation between reference frames of node1 and node2
     //! \param rotation Rotation between reference frames of node1 and node2
     //! \return if distance exceeds cutoff distance
     bool distance_within_limit(Node<Residue, BvType> *node1, 
                                Node<Residue, BvType> *node2, 
                                Vector_3D &translation, 
                                Matrix_3D &rotation) {

          bool test1 = node1->bv->distance_within_limit(*node2->bv, translation, rotation,
                                                        settings.max_interaction_distance_bv);

          // distance = node1->bv->compute_distance(*node2->bv, translation, rotation);
          // std::cout << "distance: " << distance << " " << distance*distance << "\n";
          // bool test2 = !(distance > settings.max_interaction_distance_bv);
          // assert(test1 == test2);
          return test1;
          // // DEBUG TEST: Check that distance is smaller than the distance between all pairs of vertices
          // for (uint j=0; j<node2->atoms.size(); j++) {
          //      for (uint i=0; i<node1->atoms.size(); i++) {
          //           if ((node2->atoms[j]->position - node1->atoms[i]->position).norm() < (distance-0.0001)) {
          //                std::cout << "WARNING: CHAINTREE INCONSISTENT. Distance between a pair of boxes is larger than the distance between a pair of the contained vertices: " << (node2->atoms[j]->position - node1->atoms[i]->position).norm() << " " << distance << " " << node1->bv->radius << " " << node2->bv->radius << " " << node1->atoms.size() << " " << node2->atoms.size() << " " << i << " " << j << "\n";
          //           }
          //      }
          // }

          // return !(distance > settings.max_interaction_distance_bv);
     }

     //! Determine whether pair consisting of a node and a position are separated by more
     //! than the maximum interaction distance
     //!
     //! \param node1 Node value
     //! \param position 3D-coordinate
     //! \param translation Translation between reference frames of node1 and position
     //! \param rotation Rotation between reference frames of node1 and position
     //! \return if distance exceeds cutoff distance
     bool distance_within_limit(Node<Residue, BvType> *node1, 
                                Vector_3D &position, 
                                Vector_3D &translation, 
                                Matrix_3D &rotation) {
          distance = node1->bv->computeDistance(position, translation, rotation);

          return !(distance > settings.max_interaction_distance_bv);
     }
          
     //! Determine whether pair consisting of a node and an atom are separated by more
     //! than the maximum interaction distance
     //!
     //! \param node1 Node value
     //! \param atom Atom value
     //! \param translation Translation between reference frames of node1 and atom coordinate
     //! \param rotation Rotation between reference frames of node1 and atom coordinate
     //! \return if distance exceeds cutoff distance
     bool distance_within_limit(Node<Residue, BvType> *node1, 
                                Atom &atom, 
                                Vector_3D &translation, 
                                Matrix_3D &rotation) {
          distance = node1->bv->computeDistance(atom.position, translation, rotation);

          return !(distance > settings.max_interaction_distance_bv);
     }
          
     //! Determine whether pair consisting of a node and a bounding volume are separated by more
     //! than the maximum interaction distance
     //!
     //! \param node1 Node value
     //! \param bv Bounding volume
     //! \param translation Translation between reference frames of node1 and atom coordinate
     //! \param rotation Rotation between reference frames of node1 and atom coordinate
     //! \return if distance exceeds cutoff distance
     bool distance_within_limit(Node<Residue, BvType> *node1, 
                                BvType &bv, 
                                Vector_3D &translation, 
                                Matrix_3D &rotation) {
          distance = node1->bv->computeDistance(bv, translation, rotation);
          return !(distance > settings.max_interaction_distance_bv);
     }


     //! Find the next pair of nodes
     //! This is the core iteration functionality
     //!
     //! \return True of a new valid node pair was found
     bool find_next_node_pair() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // Find new pair of bounding volumes.
          while (this->node_stack.size()) {
               Node<Residue, BvType> *node = this->node_stack.top(); this->node_stack.pop();

               // A NULL pointer on the node stack indicates that a pair should
               // be popped from the nodePair stack
               if (node==NULL) { // Check pairwise interaction
                    this->node1 = this->node_pair_stack.top(); this->node_pair_stack.pop();
                    this->node2 = this->node_pair_stack.top(); this->node_pair_stack.pop();
                    rotation = this->rotation_stack.top(); this->rotation_stack.pop();
                    translation = this->translation_stack.top(); this->translation_stack.pop();
                    separation = this->separation_stack.top();this->separation_stack.pop();

                    // // Skip the pair if the distance between them is sufficiently large
                    // if (!distance_within_limit(this->node1, this->node2, translation, rotation)) {
                    //      ((DERIVED_CLASS*)this)->on_node_pair_distance_exceeds_cutoff(this->node1,this->node2);
                    //      continue;
                    // }

                    // Skip the pair if the distance between them is sufficiently large
                    // Skip first levels without checking (gives slight performance boost)
                    if ((node1->level < (ct->nodes[ct->nodes.size()-1]->level - 5)) &&
                        !distance_within_limit(this->node1, this->node2, translation, rotation)) {
                         ((DERIVED_CLASS*)this)->on_node_pair_distance_exceeds_cutoff(this->node1,this->node2);
                         continue;
                    }

                    // Whether we are at the end level of iteration
                    bool terminal_level = this->node1->is_leaf || (settings.omit_lowest_level && this->node1->child1->is_leaf);

                    // Check whether nodes are internal
                    // if (!this->node1->is_leaf && !this->node2->is_leaf) {
                    if (!terminal_level) {

                         ((DERIVED_CLASS*)this)->on_inner_node_pair_begin();

                         if (this->node2->child2) {

                              // Translation and rotation between node1:child1 and this->node2:child2
                              Vector_3D translation_11_22 = (translation +
                                                             rotation*this->node2->child1->translation);
                              Matrix_3D rotation_11_22 = rotation * this->node2->child1->rotation;
                              if (this->node1->child2) {

                                   // Determine whether (this->node1:child2,this->node2:child2)
                                   // should be visited
                                   bool separation_12_22 = (separation ||
                                                            this->node2->child1->is_separator(this->ct->time));
                                   bool visit_12_22 = ((!settings.only_modified_pairs && 
                                                        !((DERIVED_CLASS*)this)->extra_exclusion_criterion(this->node1->child2,this->node2->child2)) ||
                                                       (this->ct->last_move_type != LOCAL && separation_12_22) ||
                                                       this->node1->child2->is_separator(this->ct->time) ||
                                                       this->node1->child2->is_modified(this->ct->time) ||
                                                       this->node2->child2->is_modified(this->ct->time));
                                   if (visit_12_22) {

                                        // Add (this->node1:child2, this->node2:child2) to stack
                                        this->node_pair_stack.push(this->node2->child2);
                                        this->node_pair_stack.push(this->node1->child2);
                                        this->node_stack.push(NULL);

                                        // Translation and rotation between this->node1:child2 and this->node2:child2
                                        Vector_3D translation_12_22 = (transpose(this->node1->child1->rotation)*
                                                                       (translation_11_22 - this->node1->child1->translation));
                                        Matrix_3D rotation_12_22 = (transpose(this->node1->child1->rotation) *
                                                                    rotation_11_22);
					
                                        // Add rotation and translation to stack
                                        this->rotation_stack.push(rotation_12_22);
                                        this->translation_stack.push(translation_12_22);
                                        this->separation_stack.push(separation_12_22);
                                   }
                                   ((DERIVED_CLASS*)this)->on_potential_child_nodes(this->node1->child2,
                                                                                    this->node2->child2);

                              }

                              // Determine whether (this->node1:child1,this->node2:child2)
                              // should be visited
                              bool separation_11_22 = (separation ||
                                                       (this->node1->child2 &&
                                                        this->node1->child2->is_separator(this->ct->time)) ||
                                                       this->node2->child1->is_separator(this->ct->time));
                              bool visit_11_22 = ((!settings.only_modified_pairs && 
                                                   !((DERIVED_CLASS*)this)->extra_exclusion_criterion(this->node1->child1, this->node2->child2)) ||
                                                  (this->ct->last_move_type != LOCAL && separation_11_22) ||
                                                  this->node1->child1->is_separator(this->ct->time) ||
                                                  this->node1->child1->is_modified(this->ct->time) ||
                                                  this->node2->child2->is_modified(this->ct->time));
                              if (visit_11_22) {
                                   
                                   // Add (this->node1:child1, this->node2:child2) to stack
                                   this->node_pair_stack.push(this->node2->child2);
                                   this->node_pair_stack.push(this->node1->child1);
                                   this->node_stack.push(NULL);

                                   // Add rotation and translation to stack
                                   this->rotation_stack.push(rotation_11_22);
                                   this->translation_stack.push(translation_11_22);
                                   this->separation_stack.push(separation_11_22);
                              } 
                              ((DERIVED_CLASS*)this)->on_potential_child_nodes(this->node1->child1,
                                                                               this->node2->child2);
                                   
                         }

                         if (this->node1->child2) {

                              // Determine whether (this->node1:child2,this->node2:child1)
                              // should be visited
                              bool separation_12_21 = separation;
                              bool visit_12_21 = ((!settings.only_modified_pairs && 
                                                   !((DERIVED_CLASS*)this)->extra_exclusion_criterion(this->node1->child2, this->node2->child1)) ||
                                                  (this->ct->last_move_type != LOCAL && separation_12_21) ||
                                                  this->node1->child2->is_separator(this->ct->time) ||
                                                  this->node1->child2->is_modified(this->ct->time) ||
                                                  this->node2->child1->is_modified(this->ct->time));
                              if (visit_12_21) {
                                   
                                   // Add (this->node1:child2, this->node2:child1) to stack
                                   this->node_pair_stack.push(this->node2->child1);
                                   this->node_pair_stack.push(this->node1->child2);
                                   this->node_stack.push(NULL);
                                        
                                   Vector_3D translation_12_21 = (transpose(this->node1->child1->rotation)
                                                                  *(translation -
                                                                    this->node1->child1->translation));
                                   Matrix_3D rotation_12_21 = (transpose(this->node1->child1->rotation) *
                                                               rotation);

                                   // Add rotation and translation to stack
                                   this->rotation_stack.push(rotation_12_21);
                                   this->translation_stack.push(translation_12_21);
                                   this->separation_stack.push(separation_12_21);
                              } 
                              ((DERIVED_CLASS*)this)->on_potential_child_nodes(this->node1->child2,
                                                                               this->node2->child1);
                         }
			 
                         // Determine whether (this->node1:child1,this->node2:child1)
                         // should be visited
                         bool separation_11_21 = (separation ||
                                                  (this->node1->child2 &&
                                                   this->node1->child2->is_separator(this->ct->time)));
                         bool visit_11_21 = ((!settings.only_modified_pairs && 
                                              !((DERIVED_CLASS*)this)->extra_exclusion_criterion(this->node1->child1,this->node2->child1)) ||
                                             (this->ct->last_move_type != LOCAL && separation_11_21) ||
                                             this->node1->child1->is_separator(this->ct->time) ||
                                             this->node1->child1->is_modified(this->ct->time) ||
                                             this->node2->child1->is_modified(this->ct->time));
                         if (visit_11_21) {
                              
                              // Add (this->node1:child1, this->node2:child1) to stack
                              this->node_pair_stack.push(this->node2->child1);
                              this->node_pair_stack.push(this->node1->child1);
                              this->node_stack.push(NULL);

                              // The transformation between the reference systems of the two nodes
                              // is given by the current rotation and translation
                              this->rotation_stack.push(rotation);
                              this->translation_stack.push(translation);
                              this->separation_stack.push(separation_11_21);
                         } 
                         ((DERIVED_CLASS*)this)->on_potential_child_nodes(this->node1->child1,
                                                                          this->node2->child1);

                    // Pair of leaf nodes 
                    } else {

                         // Ignore pair if nodes belong to the same residue,
                         // and intraResidueInteractions are turned off
                         if (std::fabs(this->node1->frame->res->index - this->node2->frame->res->index)>=
                             settings.minimum_residue_distance) {

                              bool complete = ((DERIVED_CLASS*)this)->on_leaf_node_pair_begin();
                              if (complete)
                                   return true;                   

                         } 
                    }

               // Single node popped off node_stack (inner-node interactions)
               // leaf node
               } else if (node->is_leaf || (settings.omit_lowest_level && node->child1->is_leaf)) {

                    if (settings.minimum_residue_distance == 0) {
                         this->node1 = this->node2 = node;

                         bool complete = ((DERIVED_CLASS*)this)->on_leaf_node_single_begin();
                         if (complete)
                              return true;
                    }

               // Single node popped off node_stack (inner-node interactions) 
               // inner-node
               } else {
                    this->node1 = this->node2 = node;

                    ((DERIVED_CLASS*)this)->on_inner_node_pair_begin();

                    if (!settings.only_modified_pairs || node->is_modified(this->ct->time)) {
                         if (node->child2) {

                              // Push the two children on the node_stack (intra-node interactions)
                              this->node_stack.push(node->child2);
                              this->node_stack.push(NULL);
                              this->node_stack.push(node->child1);

                              // Push the two children on the node-pair stack (inter-node interactions)
                              this->node_pair_stack.push(node->child2);
                              this->node_pair_stack.push(node->child1);
                              this->rotation_stack.push(node->child1->rotation);
                              this->translation_stack.push(node->child1->translation);
                              this->separation_stack.push(false);
			      
                         } else {

                              // If only one child was present, push it on the node_stack
                              this->node_stack.push(node->child1);
                         }
                    }                    
               }
          }

          return false;
     }
     
     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorBase& operator=(const PairIteratorBase& other) {
          this->settings = other.settings;
          this->node1 = other.node1;
          this->node2 = other.node2;
          this->entity_pair = other.entity_pair;
          this->node_stack = other.node_stack;
          this->node_pair_stack = other.node_pair_stack;
          this->rotation_stack = other.rotation_stack;
          this->translation_stack = other.translation_stack;
          this->separation_stack = other.separation_stack;
          this->rotation = other.rotation;
          this->translation = other.translation;
          this->separation = other.separation;
          return *this;
     }

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(const DERIVED_CLASS& other) const {
          return(this->node1->index == other.node1->index &&
                 this->node2->index == other.node2->index);
     }

     //! Greater than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is larger than other
     bool operator>(const DERIVED_CLASS& other) const {
          if (this->node1->index == other.node1->index) {
               return (this->node2->index > other.node2->index);
          } else {
               return (this->node1->index > other.node1->index);
          }
          return true;
     }

     //! Smaller than operator
     //!
     //! \param other Object to compare with
     //! \return True if current object is smaller than other
     bool operator<(const DERIVED_CLASS& other) const {
          if (this->node1->index == other.node1->index) {
               return (this->node2->index < other.node2->index);
          } else {
               return (this->node1->index < other.node1->index);
          }
          return true;
     }

     //! Dereference operator (*)
     //!
     //! \return underlying entity type
     const PairWithDistance<ENTITY_TYPE1*,ENTITY_TYPE2*> &operator*() const {
          return entity_pair;
     }
     
     //! Test for end of iteration
     //!
     //! \return True if iteration is at its end point
     bool end() const {
          return entity1 == NULL;
     }

     //! Code executed when potential child nodes are considered
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void on_potential_child_nodes(NodeType *node1, NodeType *node2) {
          // Nothing to be done here (only used for caching purposes)
     }

     //! Code executed when node pair is excluded due to distance cutoff
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void on_node_pair_distance_exceeds_cutoff(NodeType *node1, NodeType *node2) {
          // Nothing to be done here (only used for caching purposes)
     }

     //! Code executed when all pairs within leaf pair have been iterated over
     void on_leaf_node_pair_end() {
          // Nothing to be done here (only used for caching purposes)
     }

     //! Code executed when new inner node pair is iterated over
     void on_inner_node_pair_begin() {
          // Nothing to be done here (only used for caching purposes)
     }

     //! Makes it possible to specify an addition criterion for
     //! excluding child nodes from iteration. 
     //! Defaults to false, meaning NO EXCLUSION 
     //!
     //! \param node1 First node
     //! \param node2 Second node
     //! \return Whether interaction should be excluded
     bool extra_exclusion_criterion(NodeType *node1, NodeType *node2) {
          return false;
     }
};



//! Base class for all chaintree pair iterators for which EntityType2=Atom
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename ENTITY_TYPE1>
class PairIteratorAtomBase: public PairIteratorBase<DERIVED_CLASS,
                                                    CHAIN_TYPE,ENTITY_TYPE1,Atom> {
public:
     //! Local settings class for PairIterator<CHAIN_TYPE,Atom,Atom>.
     //! It is useful to declare your settings once and for all, instead
     //! of on every construction of the iterator, because some computations
     //! are involved
     class Settings: public PairIteratorBase<DERIVED_CLASS,CHAIN_TYPE,ENTITY_TYPE1,Atom>::Settings {

          //! Local definition of base class settings class
          typedef typename PairIteratorBase<DERIVED_CLASS,CHAIN_TYPE,ENTITY_TYPE1,Atom>::Settings PairIteratorBaseSettings;

     protected:
          //@{
          //! Indices at which wildcards are found in the selected_atoms vector
          int index_XS;
          int index_XO;
          int index_XN;
          int index_XC;
          int index_XH;
          int index_XX;
          //@}
     public:

          //! Two different ways of iterating over the atoms in a residue:
          //! ATOMS_IN_RESIDUE: iterate over all residues - checking which ones are relevant
          //! SELECTED_ATOMS: iterating over a vector of selected atoms
          enum IterationMode {ATOMS_IN_RESIDUE, SELECTED_ATOMS};

          //! Mode of iteration
          IterationMode iteration_mode;
          
          //! Selection of atoms which are of relevance to the energy term
          //! only pairs of such atoms will be returned by the iterator
          std::vector<definitions::AtomEnum> selected_atoms;

          //! For each atom, this vector contains an index to where it can
          //! be found in selected_atoms (or -1 if not present)
          //! Only used when defaultMaxInteractionDistance is set
          std::vector<int> selected_atoms_indices;

          //! Internal variable
          //! Normal mode: = selected_atoms
          //! When default_max_interaction_distance is set: = NULL
          std::vector<definitions::AtomEnum> *selected_atoms_inner_iterator;

          //! Constructor
          //!
          //! \param selected_atoms Selection of atoms which are of relevance to the energy term
          //!                       only pairs of such atoms will be returned by the iterator
          //! \param default_max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration
          Settings(const std::vector<definitions::AtomEnum> &selected_atoms,
                   double default_max_interaction_distance=UNINITIALIZED,
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false)
               : PairIteratorBaseSettings(default_max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 selected_atoms(selected_atoms) {


               // Import protein definitions (such as residue names)
               using namespace definitions;

               // If selected_atoms is set, use these to do the iteration
               if (this->selected_atoms.size() > 0)
                    // Unless default_max_interaction_distance is specified
                    if (is_initialized(default_max_interaction_distance)) {
                         this->iteration_mode = ATOMS_IN_RESIDUE;
                    } else {
                         this->iteration_mode = SELECTED_ATOMS;
                    }
               else {
                    this->iteration_mode = ATOMS_IN_RESIDUE;
               }
                              
               // Create a std::vector<int>(ATOM_ENUM_SIZE) vector, which contains the
               // index of where specific atoms can be found in maxInteractionDistances.
               // This way, an all-atom-iterator can quickly lookup these values.
               selected_atoms_indices = std::vector<int>(ATOM_ENUM_SIZE, -1);

               index_XS = -1;
               index_XO = -1;
               index_XN = -1;
               index_XC = -1;
               index_XH = -1;
               index_XX = -1;
               for (unsigned int i=0; i<selected_atoms.size(); i++) {

                    // All atoms in selected_atoms link to themselves
                    selected_atoms_indices[(int)(selected_atoms[i])] = i;

                    // Keep track of any wildcards
                    if (selected_atoms[i] == XS)
                         index_XS = i;

                    if (selected_atoms[i] == XO)
                         index_XO = i;
                         
                    if (selected_atoms[i] == XN)
                         index_XN = i;
                         
                    if (selected_atoms[i] == XC)
                         index_XC = i;
                         
                    if (selected_atoms[i] == XH)
                         index_XH = i;

                    if (selected_atoms[i] == XX)
                         index_XX = i;                         
               }

               // Iterate over all atom types and link them together with the
               // present wildcards
               for (unsigned int i=0; i<ATOM_ENUM_SIZE; i++) {
                    
                    // only consider if atom hasn't been assigned an index yet
                    if (selected_atoms_indices[i] == -1) {
                         
                         if (is_atom_XS(AtomEnum(i)) && index_XS != -1) {
                              selected_atoms_indices[i] = index_XS;
                         } else if (is_atom_XO(AtomEnum(i)) && index_XO != -1) {
                              selected_atoms_indices[i] = index_XO;
                         } else if (is_atom_XN(AtomEnum(i)) && index_XN != -1) {
                              selected_atoms_indices[i] = index_XN;
                         } else if (is_atom_XC(AtomEnum(i)) && index_XC != -1) {
                              selected_atoms_indices[i] = index_XC;
                         } else if (is_atom_XH(AtomEnum(i)) && index_XH != -1) {
                              selected_atoms_indices[i] = index_XH;
                         } else if (index_XX != -1) {
                              selected_atoms_indices[i] = index_XX;
                         }
                    }
               }


               // If no default maximum interaction distance has been set
               // and if XX is present among the selected_atoms, change to
               // ATOMS_IN_RESIDUE mode. This means that we will iterate
               // over all atoms in a resiue, rather rather than just the 
               // atoms specified in selected_atoms
               // NOTE: Derived classes should make sure that
               // default_max_interaction_distance is set. 
               if (!is_initialized(default_max_interaction_distance) && index_XX != -1) {
                    // Change iteration mode
                    this->iteration_mode = ATOMS_IN_RESIDUE;
               }

               // If we are iterating over selected atoms, expand wildcards in selected_atoms vector
               if (this->iteration_mode == SELECTED_ATOMS) {
                    std::vector<AtomEnum> new_atom_types;
                    for (std::vector<AtomEnum>::iterator it=this->selected_atoms.begin();
                         it!=this->selected_atoms.end();) {

                         if (is_atom_wildcard(*it)) {

                              // The most general wildcard should have been
                              // filtered out (if XX is present, iteration_mode
                              // should have been changed to ATOMS_IN_RESIDUE)
                              assert(*it != XX);

                              // Expand wildcard
                              AtomEnum atom_type = *it;
                              for (int j=0; j<atom_type_wildcards_size; j++) {
                                   if (atom_type_wildcards[atom_type-XS][j] == atom_type)
                                        break;
                                   new_atom_types.push_back(atom_type_wildcards[atom_type-XS][j]);
                              }

                              // // Erase the wildcard
                              it = this->selected_atoms.erase(it);
                         } else  {
                              ++it;
                         }
                    }
                    // Add the new types
                    for (unsigned int j=0; j<new_atom_types.size(); j++) {
                         this->selected_atoms.push_back(new_atom_types[j]);
                    }
               }
          }

          //! Constructor
          //!
          //! \param max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration
          Settings(double max_interaction_distance=std::numeric_limits<double>::infinity(),
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false)
               : PairIteratorBaseSettings(max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 iteration_mode(ATOMS_IN_RESIDUE),
                 selected_atoms(std::vector<definitions::AtomEnum>()),
                 selected_atoms_indices(0),
                 selected_atoms_inner_iterator(NULL) {
          }

     } settings;    //!< Local settings object

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIteratorAtomBase(CHAIN_TYPE &chain,
                          const Settings &settings)
          : PairIteratorBase<DERIVED_CLASS,
                             CHAIN_TYPE,ENTITY_TYPE1,Atom>(chain, settings),
            settings(settings) {
     }     

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorAtomBase& operator=(const PairIteratorAtomBase& other) {
          PairIteratorBase<DERIVED_CLASS,CHAIN_TYPE,ENTITY_TYPE1,Atom>::operator=(other);

          this->settings = other.settings;

          return *this;
     }
};


//! The general PairIterator implementation class. 
//! This is an empty class - the different specializations contain the functionality
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename ENTITY_TYPE1=Atom, typename ENTITY_TYPE2=Atom>
class PairIteratorImpl: public PairIteratorBase<PairIteratorImpl<DERIVED_CLASS,CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>,
                                                CHAIN_TYPE,ENTITY_TYPE1, ENTITY_TYPE2> {

protected:

     //! Settings are identical to PairIteratorBase
     typedef typename PairIteratorBase<PairIteratorImpl,
                                       CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::Settings Settings;

     //! Local settings object
     Settings settings;

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorImpl& operator=(const PairIteratorImpl& other) {
          PairIteratorBase<PairIteratorImpl,CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::operator=(other);

          this->settings = other.settings;

          return *this;
     }     
};


//! Pair-iterator implementation specialization for Atom-Atom pairs 
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class PairIteratorImpl<DERIVED_CLASS,CHAIN_TYPE,Atom,Atom>: public PairIteratorAtomBase<DERIVED_CLASS,
                                                                                        CHAIN_TYPE,Atom> {

public:
     //! Local settings class.
     //! Settings for PairIterator<CHAIN_TYPE,Atom,Atom>
     //! It is useful to declare your settings once and for all, instead
     //! of on every construction of the iterator, because some computations
     //! are involved
     class Settings: public PairIteratorAtomBase<DERIVED_CLASS,
                                                 CHAIN_TYPE,Atom>::Settings {

          //! Local definition of base class settings class
          typedef typename PairIteratorAtomBase<DERIVED_CLASS,
                                                CHAIN_TYPE,Atom>::Settings PairIteratorBaseSettings;

     public:

          //! For all pairs of selected atoms, specifies the maxInteractionDistance
          std::vector<std::vector<double> > max_interaction_distance_matrix;

          //! Constructor
          //!
          //! \param selected_atoms Selection of atoms which are of relevance to the energy term
          //!                       only pairs of such atoms will be returned by the iterator
          //! \param max_interaction_distance_matrix Matrix of maximum distances for all atoms in selected atoms vector
          //! \param default_max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration          
          Settings(const std::vector<definitions::AtomEnum> &selected_atoms,
                   const std::vector<std::vector<double> > &max_interaction_distance_matrix,
                   double default_max_interaction_distance=UNINITIALIZED,
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level = false)
               : PairIteratorBaseSettings(selected_atoms,
                                          default_max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 max_interaction_distance_matrix(max_interaction_distance_matrix) {

               // Find max_interaction_distance_bv: max(max_interaction_distance_matrix)
               this->max_interaction_distance_bv = 0.0;
               if (is_initialized(default_max_interaction_distance))
                    this->max_interaction_distance_bv = default_max_interaction_distance;                    
               for (unsigned int i=0; i<max_interaction_distance_matrix.size(); i++) {
                    for (unsigned int j=0; j<max_interaction_distance_matrix[i].size(); j++) {
                         if (max_interaction_distance_matrix[i][j] > this->max_interaction_distance_bv) {
                              this->max_interaction_distance_bv = max_interaction_distance_matrix[i][j];
                         }
                    }
               }

               // If no default maximum interaction distance has been set
               // and if XX is present among the selected_atoms, set the
               // default_max_interaction_distance. This means that we will
               // iterate over all atoms in a residue, rather than just the atoms
               // specified in selected_atoms               
               if (!is_initialized(default_max_interaction_distance) && this->index_XX != -1) {

                    // Change iteration mode (just to make sure - this should already have been done
                    // by base class)
                    this->iteration_mode = PairIteratorBaseSettings::ATOMS_IN_RESIDUE;

                    // Find default interaction distance (for distances not in distance matrix)
                    this->default_max_interaction_distance = 0.0;
                    for (unsigned int i=0; i<definitions::ATOM_ENUM_SIZE; i++) {
                         if (max_interaction_distance_matrix[i][this->index_XX] > this->default_max_interaction_distance) {
                              this->default_max_interaction_distance = max_interaction_distance_matrix[i][this->index_XX];
                         }
                    }
               }
          }


          //! Constructor
          //!
          //! \param max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration
          Settings(double max_interaction_distance=std::numeric_limits<double>::infinity(),
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level = false)
               : PairIteratorBaseSettings(max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 max_interaction_distance_matrix(0) {
          }

     } settings;    //!< Local settings object

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     // Friends
     friend class PairIteratorOneVsAll<CHAIN_TYPE,Atom,Atom>;
     friend class PairIteratorOneVsAllBase<PairIteratorOneVsAll<CHAIN_TYPE,Atom,Atom>,CHAIN_TYPE,Atom,Atom>;

     //! Iterator for node 1
     typename Node<Residue, BvType>::Iterator node_iterator1;

     //! Iterator for node 2
     typename Node<Residue, BvType>::Iterator node_iterator2;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIteratorImpl(CHAIN_TYPE &chain,
                      const Settings &settings=Settings())
          : PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,Atom>(chain, settings),
            settings(settings) {
     }     

     //! Initializer
     void init() {

          // If default_max_interaction_distance is set, the node_iterator
          // should ignore the selected_atoms vector, since we will
          // have to consider all atoms
          if (is_initialized(settings.default_max_interaction_distance)) {
               settings.selected_atoms_inner_iterator = NULL;
          }

          // Call base class init member function
          PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,Atom>::init();                    
     }

     //! Test whether atoms are within cutoff distance
     bool atom_pair_within_distance() {
          int selected_atom_index1 = (settings.selected_atoms_indices)[this->entity1->atom_type];
          int selected_atom_index2 = (settings.selected_atoms_indices)[this->entity2->atom_type];
          double max_distance = settings.max_interaction_distance_matrix[selected_atom_index1][selected_atom_index2];
          return (this->entity1 && ((this->distance=(this->entity1->position - this->entity2->position).norm()) <= max_distance));
     }

     //! Test whether atoms are within cutoff distance - falling back on default interaction distance
     bool atom_pair_within_distance_or_default_distance() {
          int selected_atom_index1 = (settings.selected_atoms_indices)[this->entity1->atom_type];
          int selected_atom_index2 = (settings.selected_atoms_indices)[this->entity2->atom_type];
          double max_distance = settings.default_max_interaction_distance;
          if (selected_atom_index1 >=0 && selected_atom_index2 >= 0) {
               max_distance = settings.max_interaction_distance_matrix[selected_atom_index1][selected_atom_index2];
          }
          return (this->entity1 && ((this->distance=(this->entity1->position - this->entity2->position).norm()) <= max_distance));
     }

     //! Test whether atoms are within cutoff distance - using default interaction distance
     bool atom_pair_within_default_distance() {
          double max_distance = settings.default_max_interaction_distance;
          return (this->entity1 && ((this->distance=(this->entity1->position - this->entity2->position).norm()) <= max_distance));
     }

     //! Code executed when base class iteration hits a pair of leaves
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_pair_begin() {
          node_iterator1 = typename Node<Residue, BvType>::
               Iterator(this->node1, 
                        settings.selected_atoms_inner_iterator,
                        &settings.selected_atoms_indices);
          node_iterator2 = typename Node<Residue, BvType>::
               Iterator(this->node2, 
                        settings.selected_atoms_inner_iterator,
                        &settings.selected_atoms_indices);

          bool complete = false;
          if (*node_iterator1 && *node_iterator2) {
               this->entity1 = *node_iterator1;
               this->entity2 = *node_iterator2;
               complete = true;
          }
          return complete;
     }

     //! Code executed when all pairs within leaf pair have been iterated over
     void on_leaf_node_pair_end() {
          // Nothing to be done here (only used for caching purposes)
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {

          node_iterator1 = typename Node<Residue, BvType>::
               Iterator(this->node1, settings.selected_atoms_inner_iterator, &settings.selected_atoms_indices);
          // Due to symmetry, only half of the pairs
          // need to be considered
          node_iterator2 = typename Node<Residue, BvType>::
               Iterator(node_iterator1);
          ++node_iterator2;

          bool complete = false;
          if (*node_iterator1 && *node_iterator2) {
               this->entity1 = *node_iterator1;
               this->entity2 = *node_iterator2;
               complete = true;
          }
          return complete;

     }

     //! Find next pair of atoms
     void find_next_pair() {

          // First attempt iteration over current node pair
          ++node_iterator2;
          if (!*node_iterator2) {
               ++node_iterator1;
               if (*node_iterator1) {

                    // within node iteration
                    if (this->node1 == this->node2) {

                         // Because of symmetry, only half of the pairs
                         // need to be considered
                         node_iterator2 = typename Node<Residue, BvType>::Iterator(node_iterator1);
                         ++node_iterator2;
                              
                         // pair iteration 
                    } else {
                         node_iterator2 = typename Node<Residue, BvType>::
                              Iterator(this->node2,
                                       settings.selected_atoms_inner_iterator,
                                       &settings.selected_atoms_indices);
                    } 
               }
          }

          if (*node_iterator1 && *node_iterator2) {
               this->entity1 = *node_iterator1;
               this->entity2 = *node_iterator2;
               return;
          }

          // Execute code when all pairs within leaf pair have been iterated over
          ((DERIVED_CLASS*)this)->on_leaf_node_pair_end();

          // If no atom-pairs left in current node pair, find new pair.
          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }

     //! Increment to find next relevant pair
     void increment() {

          // Find next pair of atoms
          find_next_pair();

          // If selected_atoms is set, only return atoms within specified distance
          if (settings.selected_atoms.size()>0) {

               // Iterate until an atom pair within max distance is found
               if (is_initialized(settings.default_max_interaction_distance)) {
                    while (!this->end()) {
                         if (!this->atom_pair_within_distance_or_default_distance())
                              find_next_pair();
                         else
                              break;
                    }
               } else {
                    while (!this->end()) {
                         if (!this->atom_pair_within_distance())
                              find_next_pair();
                         else
                              break;
                    }
               }

          // Otherwise, use default distance 
          } else {
               while (!this->end()) {
                    if (!this->atom_pair_within_default_distance())
                         find_next_pair();
                    else
                         break;
               }               
          }
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     DERIVED_CLASS& operator++() {
          this->increment();
          return (*((DERIVED_CLASS*)this));          
     }
     
     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DERIVED_CLASS& operator=(const DERIVED_CLASS& other) {
          PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,Atom>::operator=(other);
          this->settings = other.settings;
          this->node_iterator1 = other.node_iterator1;
          this->node_iterator2 = other.node_iterator2;
          return (*((DERIVED_CLASS*)this));          
     }     
};


//! Pair-iterator specialization for node-node pairs 
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class PairIteratorImpl<DERIVED_CLASS, 
                       CHAIN_TYPE,
                       typename CHAIN_TYPE::ChainTree::NodeType,
                       typename CHAIN_TYPE::ChainTree::NodeType>
     : public PairIteratorBase<DERIVED_CLASS,
                               CHAIN_TYPE,
                               typename CHAIN_TYPE::ChainTree::NodeType,
                               typename CHAIN_TYPE::ChainTree::NodeType> {

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     // Friends
     friend class PairIteratorOneVsAll<CHAIN_TYPE,RSS,RSS>;

public:

     //! Local settings class - identical to PairIteratorBase
     typedef typename PairIteratorBase<DERIVED_CLASS,
                                       CHAIN_TYPE,NodeType,NodeType>::Settings Settings;
     //! Local settings object
     Settings settings;

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIteratorImpl(CHAIN_TYPE &chain,
                      const Settings &settings=Settings())
          : PairIteratorBase<DERIVED_CLASS,
                             CHAIN_TYPE,
                             NodeType,
                             NodeType > (chain, settings),
            settings(settings) {
     }     


     //! Code executed when base class iteration hits a pair of leaves
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_pair_begin() {
          this->entity1 = this->node1;
          this->entity2 = this->node2;
          return true;
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {
          this->entity1 = this->node1;
          this->entity2 = this->node2;
          return true;
     }

     //! Find next node,node pair
     void find_next_pair() {

          // Execute code when all pairs within leaf pair have been iterated over
          ((DERIVED_CLASS*)this)->on_leaf_node_pair_end();

          // Since we are only interested in pairs of nodes
          // simply call find_next_node_pair
          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }

     //! Increment to find next relevant pair
     void increment() {
          this->find_next_pair();
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     DERIVED_CLASS& operator++() {
          increment();
          return (*((DERIVED_CLASS*)this));          
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DERIVED_CLASS& operator=(const DERIVED_CLASS& other) {
          PairIteratorBase<DERIVED_CLASS,CHAIN_TYPE,NodeType,NodeType>::operator=(other);
          this->settings = other.settings;
          return (*((DERIVED_CLASS*)this));          
     }
};


//! Pair-iterator specialization for node-atom pairs 
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class PairIteratorImpl<DERIVED_CLASS,
                       CHAIN_TYPE,
                       typename CHAIN_TYPE::ChainTree::NodeType,
                       Atom>
     : public PairIteratorAtomBase<DERIVED_CLASS,
                                   CHAIN_TYPE,
                                   typename CHAIN_TYPE::ChainTree::NodeType> {

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Iterator for node 2
     typename Node<Residue, BvType>::Iterator node_iterator2;

public:

     //! Local Settings class.
     class Settings: public PairIteratorAtomBase<DERIVED_CLASS,
                                                 CHAIN_TYPE,NodeType>::Settings {

          //! Local definition of base class settings class
          typedef typename PairIteratorAtomBase<DERIVED_CLASS,
                                                CHAIN_TYPE,NodeType>::Settings PairIteratorBaseSettings;

     public:

          //! For all pairs of selected atoms, specifies the max_interaction_distance
          std::vector<double> max_interaction_distance_vector;

          //! Constructor
          //!
          //! \param selected_atoms Selection of atoms which are of relevance to the energy term
          //!                       only pairs of such atoms will be returned by the iterator
          //! \param max_interaction_distance_vector Vector of maximum distances for all atoms in selected atoms vector
          //! \param default_max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration                    
          Settings(const std::vector<definitions::AtomEnum> &selected_atoms,
                   const std::vector<double> &max_interaction_distance_vector,
                   double default_max_interaction_distance=UNINITIALIZED,
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false)
               : PairIteratorBaseSettings(selected_atoms,
                                          default_max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 max_interaction_distance_vector(max_interaction_distance_vector) {

               // Find max_interaction_distance_bv: max(max_interaction_distance_vector)
               this->max_interaction_distance_bv = 0.0;
               if (is_initialized(default_max_interaction_distance))
                    this->max_interaction_distance_bv = default_max_interaction_distance;                    
               for (unsigned int i=0; i<max_interaction_distance_vector.size(); i++) {
                    if (max_interaction_distance_vector[i] > this->max_interaction_distance_bv) {
                         this->max_interaction_distance_bv = max_interaction_distance_vector[i];
                    }
               }

               // If no default maximum interaction distance has been set
               // and if XX is present among the selected_atoms, set the
               // default_max_interaction_distance. This means that we will
               // iterate over all atoms in a residue, rather than just the atoms
               // specified in selected_atoms               
               if (!is_initialized(default_max_interaction_distance) && this->index_XX != -1) {

                    // Change iteration mode (just to make sure - this should already have been done
                    // by base class)
                    this->iteration_mode = PairIteratorBaseSettings::ATOMS_IN_RESIDUE;

                    // Set default interaction distance (for distances not in distance matrix)
                    this->default_max_interaction_distance = max_interaction_distance_vector[this->index_XX];
               }
          }

          //! Constructor
          //!
          //! \param max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration
          Settings(double max_interaction_distance=std::numeric_limits<double>::infinity(),
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false)
               : PairIteratorBaseSettings(max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 max_interaction_distance_vector(0) {
          }

     } settings;    //!< Local settings object

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIteratorImpl(CHAIN_TYPE &chain,
                      const Settings &settings=Settings())
          : PairIteratorAtomBase<DERIVED_CLASS,
                                 CHAIN_TYPE,
                                 NodeType> (chain, settings),
            settings(settings){
     }     

     //! Initializer
     void init() {

          // If default_max_interaction_distance is set, the node_iterator
          // should ignore the selected_atoms vector, since we will
          // have to consider all atoms
          if (is_initialized(settings.default_max_interaction_distance)) {
               settings.selected_atoms_inner_iterator = NULL;
          }

          // Call base class init member function
          PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,NodeType>::init();                    
     }

     //! Test whether atoms are within cutoff distance
     bool atom_node_pair_within_distance() {
          // For this iterator, self interactions are always within distance
          if (this->node1 == this->node2) {
               this->distance = 0.0;
               return true;
          }
          int selected_atom_index2 = (settings.selected_atoms_indices)[this->entity2->atom_type];
          double max_distance = (settings.max_interaction_distance_vector)[selected_atom_index2];
          if (this->entity2) {
               Vector_3D local_pos = this->node2->local_positions[this->node2->atomIndex[this->entity2->atom_type]];
               this->distance = this->entity1->bv->computeDistance(local_pos, this->translation, this->rotation);
               return this->distance <= max_distance;
          } else {
               return false;
          }
     }

     //! Test whether atoms are within cutoff distance - falling back on default interaction distance
     bool atom_node_pair_within_distance_or_default_distance() {

          // For this iterator, self interactions are always within distance
          if (this->node1 == this->node2) {
               this->distance = 0.0;
               return true;
          }
          double max_distance = settings.default_max_interaction_distance;
          int selected_atom_index2 = (settings.selected_atoms_indices)[this->entity2->atom_type];
          if (selected_atom_index2 >= 0) {
               max_distance = settings.max_interaction_distance_vector[selected_atom_index2];
          }
          if (this->entity2) {
               Vector_3D local_pos = this->node2->local_positions[this->node2->atomIndex[this->entity2->atom_type]];
               this->distance = this->entity1->bv->computeDistance(local_pos, this->translation, this->rotation);
               return this->distance <= max_distance;
          } else {
               return false;
          }
     }

     //! Test whether atoms are within cutoff distance - using default interaction distance
     bool atom_node_pair_within_default_distance() {
          // For this iterator, self interactions are always within distance
          if (this->node1 == this->node2) {
               this->distance = 0.0;
               return true;
          }
          double max_distance = settings.default_max_interaction_distance;
          if (this->entity2) {
               Vector_3D local_pos = this->node2->local_positions[this->node2->atomIndex[this->entity2->atom_type]];
               this->distance = this->entity1->bv->computeDistance(local_pos, this->translation, this->rotation);
               return this->distance <= max_distance;
          } else {
               return false;
          }
     }

     //! Code executed when base class iteration hits a pair of leaves
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_pair_begin() {
          node_iterator2 = typename Node<Residue, BvType>::
               Iterator(this->node2, 
                        settings.selected_atoms_inner_iterator,
                        &settings.selected_atoms_indices);

          bool complete = false;
          if (*node_iterator2) {
               this->entity1 = this->node1;
               this->entity2 = *node_iterator2;
               complete = true;
          }
          return complete;
     }

     //! Code executed when all pairs within leaf pair have been iterated over
     void on_leaf_node_pair_end() {
          // Nothing to be done here (only used for caching purposes)
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {

          node_iterator2 = typename Node<Residue, BvType>::
               Iterator(this->node2,
                        settings.selected_atoms_inner_iterator,
                        &settings.selected_atoms_indices);

          bool complete = false;
          if (*node_iterator2) {
               this->entity1 = this->node1;
               this->entity2 = *node_iterator2;
               complete = true;
          }
          return complete;
     }

     //! Find next node,atom pair
     void find_next_pair() {

          // First attempt iteration over current node pair
          ++node_iterator2;
          if (*node_iterator2) {
               this->entity1 = this->node1;
               this->entity2 = *node_iterator2;
               return;
          }

          // Execute code when all pairs within leaf pair have been iterated over
          ((DERIVED_CLASS*)this)->on_leaf_node_pair_end();

          // If no pairs left in current node pair, find new pair.
          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }



     //! Increment to find next relevant pair
     void increment() {

          // Find next node,atom pair
          this->find_next_pair();

          // If selected_atoms is set, only return atoms within specified distance
          if (settings.selected_atoms.size()>0) {

               // Iterate until an atom pair within max distance is found
               if (is_initialized(settings.default_max_interaction_distance)) {
                    while (!this->end()) {
                         if (!this->atom_node_pair_within_distance_or_default_distance())
                              this->find_next_pair();
                         else
                              break;
                    }
               } else {
                    while (!this->end()) {
                         if (!this->atom_node_pair_within_distance())
                              this->find_next_pair();
                         else
                              break;
                    }
               }
          // Otherwise, use default distance 
          } else {
               while (!this->end()) {
                    if (!this->atom_node_pair_within_default_distance())
                         this->find_next_pair();
                    else
                         break;
               }               
          }
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     DERIVED_CLASS& operator++() {
          increment();
          return (*((DERIVED_CLASS*)this));          
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DERIVED_CLASS& operator=(const DERIVED_CLASS& other) {
          PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,NodeType>::operator=(other);
          this->settings = other.settings;
          this->node_iterator2 = other.node_iterator2;
          return (*((DERIVED_CLASS*)this));          
     }
};



//! Pair-iterator specialization for boundingvolume-atom pairs 
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename DERIVED_CLASS, typename CHAIN_TYPE>
class PairIteratorImpl<DERIVED_CLASS,
                       CHAIN_TYPE,
                       typename CHAIN_TYPE::ChainTree::BvType,
                       Atom>
     : public PairIteratorAtomBase<DERIVED_CLASS,
                                   CHAIN_TYPE,
                                   typename CHAIN_TYPE::ChainTree::BvType> {

protected:

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Iterator for node 2
     typename Node<Residue, BvType>::Iterator node_iterator2;

     //! Local Settings class.
     class Settings: public PairIteratorAtomBase<DERIVED_CLASS,
                                                 CHAIN_TYPE,BvType>::Settings {

          //! Local definition of base class settings class
          typedef typename PairIteratorAtomBase<DERIVED_CLASS,
                                                CHAIN_TYPE,BvType>::Settings PairIteratorBaseSettings;

     public:

          //! For all pairs of selected atoms, specifies the max_interaction_distance
          std::vector<double> max_interaction_distance_vector;

          //! Constructor
          //!
          //! \param selected_atoms Selection of atoms which are of relevance to the energy term
          //!                       only pairs of such atoms will be returned by the iterator
          //! \param max_interaction_distance_vector Matrix of maximum distances for all atoms in selected atoms vector
          //! \param default_max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration          
          Settings(const std::vector<definitions::AtomEnum> &selected_atoms,
                   const std::vector<double> &max_interaction_distance_vector,
                   double default_max_interaction_distance=UNINITIALIZED,
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false)
               : PairIteratorBaseSettings(selected_atoms,
                                          default_max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 max_interaction_distance_vector(max_interaction_distance_vector) {
          }

          //! Constructor
          //!
          //! \param max_interaction_distance An overall value for the maximum interaction distance between elements in a pair
          //! \param only_modified_pairs Whether only to report pairs that were modified in the last move
          //! \param minimum_residue_distance Minimum distance along chain between elements in a pair (measured in residues)
          //! \param omit_lowest_level Whether to omit the most fine grained level of the iteration
          Settings(double max_interaction_distance=std::numeric_limits<double>::infinity(),
                   bool only_modified_pairs=false,
                   int minimum_residue_distance=0,
                   bool omit_lowest_level=false)
               : PairIteratorBaseSettings(max_interaction_distance,
                                          only_modified_pairs,
                                          minimum_residue_distance,
                                          omit_lowest_level),
                 max_interaction_distance_vector(0) {
          }

     } settings;    //!< Local settings object

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIteratorImpl(CHAIN_TYPE &chain,
                      const Settings &settings=Settings())
          : PairIteratorAtomBase<DERIVED_CLASS,
                                 CHAIN_TYPE,
                                 BvType> (chain, settings),
            settings(settings){
     }     

     //! Find next bv-atom pair
     void find_next_pair() {

          // First attempt iteration over current node pair
          ++node_iterator2;
          if (*node_iterator2) {
               this->entity1 = this->node1->bv;
               this->entity2 = *node_iterator2;
               return;
          }

          // Execute code when all pairs within leaf pair have been iterated over
          ((DERIVED_CLASS*)this)->on_leaf_node_pair_end();

          // If no pairs left in current node pair, find new pair.
          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }



     //! Code executed when base class iteration hits a pair of leaves
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_pair_begin() {
          node_iterator2 = typename Node<Residue, BvType>::
               Iterator(this->node2, 
                        settings.selected_atoms_inner_iterator,
                        &settings.selected_atoms_indices);

          bool complete = false;
          if (*node_iterator2) {
               this->entity1 = this->node1->bv;
               this->entity2 = *node_iterator2;
               complete = true;
          }
          return complete;
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {

          node_iterator2 = typename Node<Residue, BvType>::
               Iterator(this->node2,
                        settings.selected_atoms_inner_iterator,
                        &settings.selected_atoms_indices);

          bool complete = false;
          if (*node_iterator2) {
               this->entity1 = this->node1->bv;
               this->entity2 = *node_iterator2;
               complete = true;
          }
          return complete;
     }

     //! Test whether atoms are within cutoff distance
     bool atom_node_pair_within_distance() {
          // For this iterator, self interactions are always within distance
          if (this->node1 == this->node2) {
               this->distance = 0.0;
               return true;
          }
          int selected_atom_index2 = (settings.selected_atoms_indices)[this->entity2->atom_type];
          double max_distance = (settings.max_interaction_distance_vector)[selected_atom_index2];
          if (this->entity2) {
               Vector_3D local_pos = this->node2->local_positions[this->node2->atomIndex[this->entity2->atom_type]];
               this->distance = this->entity1->computeDistance(local_pos, this->translation, this->rotation);
               return this->distance <= max_distance;
          } else {
               return false;
          }
     }

     //! Test whether atoms are within cutoff distance - falling back on default interaction distance
     bool atom_node_pair_within_distance_or_default_distance() {

          // For this iterator, self interactions are always within distance
          if (this->node1 == this->node2) {
               this->distance = 0.0;
               return true;
          }
          int selected_atom_index2 = (settings.selected_atoms_indices)[this->entity2->atom_type];
          double max_distance = settings.default_max_interaction_distance;
          if (selected_atom_index2 >= 0) {
               max_distance = settings.max_interaction_distance_vector[selected_atom_index2];
          }
          if (this->entity2) {
               Vector_3D local_pos = this->node2->local_positions[this->node2->atomIndex[this->entity2->atom_type]];
               this->distance = this->entity1->computeDistance(local_pos, this->translation, this->rotation);
               return this->distance <= max_distance;
          } else {
               return false;
          }
     }

     //! Test whether atoms are within cutoff distance - using default interaction distance
     bool atom_node_pair_within_default_distance() {

          // For this iterator, self interactions are always within distance
          if (this->node1 == this->node2) {
               this->distance = 0.0;
               return true;
          }
          double max_distance = settings.default_max_interaction_distance;
          if (this->entity2) {
               Vector_3D local_pos = this->node2->local_positions[this->node2->atomIndex[this->entity2->atom_type]];
               this->distance = this->entity1->computeDistance(local_pos, this->translation, this->rotation);
               return this->distance <= max_distance;
          } else {
               return false;
          }
     }


     //! Increment to find next relevant pair
     void increment() {

          // Find next bv-atom pair
          this->find_next_pair();

          // If selected_atoms is set, only return atoms within specified distance
          if (settings.selected_atoms.size()>0) {

               // Iterate until an atom pair within max distance is found
               if (is_initialized(settings.default_max_interaction_distance)) {
                    while (!this->end()) {
                         if (!this->atom_node_pair_within_distance_or_default_distance())
                              this->find_next_pair();
                         else
                              break;
                    }
               } else {
                    while (!this->end()) {
                         if (!this->atom_node_pair_within_distance())
                              this->find_next_pair();
                         else
                              break;
                    }
               }
          // Otherwise, use default distance 
          } else {
               while (!this->end()) {
                    if (!this->atom_node_pair_within_default_distance())
                         this->find_next_pair();
                    else
                         break;
               }               
          }
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     DERIVED_CLASS& operator++() {
          increment();
          return (*((DERIVED_CLASS*)this));          
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DERIVED_CLASS& operator=(const DERIVED_CLASS& other) {
          PairIteratorAtomBase<DERIVED_CLASS,CHAIN_TYPE,BvType>::operator=(other);
          this->settings = other.settings;
          this->node_iterator2 = other.node_iterator2;
          return (*((DERIVED_CLASS*)this));          
     }
};





//! The PairIterator class.
//! This is a simple front end to the PairIteratorImpl classes 
//!
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
template <typename CHAIN_TYPE, typename ENTITY_TYPE1=Atom, typename ENTITY_TYPE2=Atom>
class PairIterator: public PairIteratorImpl<PairIterator<CHAIN_TYPE, ENTITY_TYPE1, ENTITY_TYPE2>,
                                            CHAIN_TYPE,ENTITY_TYPE1, ENTITY_TYPE2> {

public:

     //! Trick to simulate a templated typedef
     //! Allows Cached iterator to access the implementation class of the PairIterator that is passed
     //! along as template parameter
     template <typename DERIVED_CLASS>
     struct Implementation {
          //! Implementation type definition
          typedef PairIteratorImpl<DERIVED_CLASS,
                                   CHAIN_TYPE,ENTITY_TYPE1, ENTITY_TYPE2> Type;
     };

     //! Local settings class - identical to PairIteratorImpl     
     typedef typename PairIteratorImpl<PairIterator,CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::Settings Settings;

     //! Local settings object     
     Settings settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local settings object
     PairIterator(CHAIN_TYPE &chain,
                  const Settings &settings=Settings())
          : PairIteratorImpl<PairIterator,CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>(chain, settings) {
          this->init();
     }     

     //! Overload assignment operator.
     //! This method is inherited automatically by GCC compilers, but needed explicitly by ICC.
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIterator &operator=(const PairIterator &other) {
          PairIteratorImpl<PairIterator<CHAIN_TYPE, ENTITY_TYPE1, ENTITY_TYPE2>,
                           CHAIN_TYPE,ENTITY_TYPE1, ENTITY_TYPE2>::operator=(other);
          return *this;
     }
};








//! Base class for all one-vs-all chaintree pair iterators
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2>
class PairIteratorOneVsAllBase: public PairIteratorImpl<DERIVED_CLASS, 
							CHAIN_TYPE,
							ENTITY_TYPE1,
							ENTITY_TYPE2> {

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     


protected:
     //! Settings are identical to PairIteratorImpl
     typedef typename PairIteratorImpl<DERIVED_CLASS, CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::Settings Settings;

     //! Local settings object
     Settings settings;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param entity The fixed element that we are comparing to
     //! \param settings Local settings object
     PairIteratorOneVsAllBase(CHAIN_TYPE &chain,
                              ENTITY_TYPE2 &entity,
                              const Settings &settings=Settings())
          : PairIteratorImpl<DERIVED_CLASS, CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>(chain, settings),
            settings(settings) {
     }     

     //! Initializer
     void init() {

          // Initialize node stack at root node
          // In this iteration mode, only single nodes are pushed/popped
          // on the stack
          Node<Residue, BvType> *root = this->ct->nodes[this->ct->nodes.size()-1];
          this->node_stack.push(root);

          // The frame of the root node is not necessarily
          // consistent with the current positions -- only frames
          // of residues in which dofs are modified are updated. We
          // therefore update the frame information in the root node
          root->frame->update();

          // Specify rotation and translation between the root node
          // and the global reference system
          this->rotation_stack.push(transpose(*root->frame->orientation));
          this->translation_stack.push(root->frame->transform(Vector_3D(0,0,0)));

          // Start at first pair
          ++(*((DERIVED_CLASS*)this));
          
     }

     //! Find the next pair of nodes
     //! This is the core iteration functionality
     //!
     //! \return True of a new valid node pair was found
     bool find_next_node_pair() {

          // Find new pair of bounding volumes.
          while (this->node_stack.size()) {

               Node<Residue, BvType> *node = this->node_stack.top(); this->node_stack.pop();

               this->rotation = this->rotation_stack.top(); this->rotation_stack.pop();
               this->translation = this->translation_stack.top(); this->translation_stack.pop();

               // Check whether distance is within cutoff
               // Node: This works because the rotation and translation at the root 
               // node have been set appropriately (see init)
               // if (!distance_within_limit(node, this->entity1->position, this->translation, this->rotation)) {
               if (!distance_within_limit(node, *this->entity2, this->translation, this->rotation)) {
                    continue;
               }

               if (node->is_leaf) {

                    // Ignore node the minimum_residue_check is fulfilled
                    if (!(((DERIVED_CLASS*)this)->within_minimum_residue_distance(node,this->entity2))) {
                         
                         this->node1 = node;

                         bool complete = ((DERIVED_CLASS*)this)->on_leaf_node_single_begin();
                         if (complete)
                              return true;
                    }
               } else {
                    if (node->child2) {
                         // Push child2 on stack
                         this->node_stack.push(node->child2);

                         // Push rotation and translation on stack
                         this->rotation_stack.push(transpose(node->child1->rotation)*this->rotation);
                         this->translation_stack.push(transpose(node->child1->rotation)*
                                                      (this->translation - node->child1->translation));
                    }

                    // Push child1 on stack
                    this->node_stack.push(node->child1);
                    
                    // Child1 shares frame with parent
                    this->rotation_stack.push(this->rotation);
                    this->translation_stack.push(this->translation);
                    
                    // Separation detection not implemented for OneVsAll iteration
                    // separationStack.push(????);
               }
          }
          return false;
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DERIVED_CLASS& operator=(const DERIVED_CLASS& other) {
          PairIteratorImpl<DERIVED_CLASS, CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::operator=(other);
          this->settings = other.settings;
          return *((DERIVED_CLASS*)this);
     }
};



//! The general PairIteratorOneVsAll class. 
//! This is an empty class - the different specializations contain the functionality
//!
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
template <typename CHAIN_TYPE, typename ENTITY_TYPE1=Atom, typename ENTITY_TYPE2=Atom>
class PairIteratorOneVsAll: public PairIteratorOneVsAllBase<PairIteratorOneVsAll<CHAIN_TYPE,
                                                                                 ENTITY_TYPE1, ENTITY_TYPE2>,
                                                            CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2> {     
};



//! Pair-iterator specialization for atom-vs-all - returning atoms
//!
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename CHAIN_TYPE>
class PairIteratorOneVsAll<CHAIN_TYPE,Atom,Atom>: public PairIteratorOneVsAllBase<PairIteratorOneVsAll<CHAIN_TYPE,
                                                                                                       Atom,Atom>,
                                                                                  CHAIN_TYPE,Atom,Atom> {

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

public:

     //! Settings are identical to PairIterator
     typedef typename PairIterator<CHAIN_TYPE,Atom,Atom>::Settings Settings;

     //! Local settings object
     Settings settings;

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param entity The fixed element that we are comparing to
     //! \param settings Local settings object
     PairIteratorOneVsAll(CHAIN_TYPE &chain,
                          Atom &entity,
                          const Settings &settings=Settings())
          : PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,Atom,Atom>(chain, entity, settings),
            settings(settings) {

          this->entity2 = &entity;
          this->init();
     }     

     //! Initializer
     void init() {
          
          // If defaultMaxInteractionDistance is set, the nodeIterator
          // should ignore the selected_atoms vector, since we will
          // have to consider all atoms
          if (is_initialized(settings.default_max_interaction_distance)) {
               settings.selected_atoms_inner_iterator = NULL;
          }

          PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,Atom,Atom>::init();          
     }

     //! Check whether node and entity are within residue distance
     //! Called by base class
     //!
     //! \param node Current node value
     //! \param entity The fixed element that we are comparing to
     //! \return True if the residues of node and entity are not separated by more than the minimum distance
     bool within_minimum_residue_distance(NodeType *node, Atom *entity) {
          return (fabs(node->frame->res->index - entity->residue->index)<
                  settings.minimum_residue_distance);
     }


     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {

          this->node_iterator1 = 
               typename Node<Residue, BvType>::Iterator(this->node1, 
                                                         settings.selected_atoms_inner_iterator,
                                                         &settings.selected_atoms_indices);
          bool complete = false;
          if (*this->node_iterator1) {
               this->entity1 = *this->node_iterator1;
               complete = true;
          }
          return complete;

     }

     //! Find next atom-atom pair
     void find_next_pair() {

          // First attempt iteration over atoms in current node
          ++this->node_iterator1;
          if (*this->node_iterator1) {
               this->entity1 = *this->node_iterator1;
               return;
          }

          // If no atom-pairs left in current node pair, find new pair.
          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }


     //! Increment operator
     //!
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator++() {

          // Find next atom-atom pair
          find_next_pair();

          // If selected_atoms is set, only return atoms within specified distance
          if (settings.selected_atoms.size()>0) {

               // Iterate until an atom pair within max distance is found
               if (is_initialized(settings.default_max_interaction_distance)) {
                    while (!this->end()) {
                         if (!this->atom_pair_within_distance_or_default_distance())
                              find_next_pair();
                         else
                              break;
                    }
               } else {
                    while (!this->end()) {
                         if (!this->atom_pair_within_distance())
                              find_next_pair();
                         else
                              break;
                    }
               }
          // Otherwise, use default distance 
          } else {
               while (!this->end()) {
                    if (!this->atom_pair_within_default_distance())
                         find_next_pair();
                    else
                         break;
               }               
          } 
          return (*this);
     }


     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator=(const PairIteratorOneVsAll& other) {
          PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,Atom,Atom>::operator=(other);
          this->settings = other.settings;
          return *this;
     }
};


//! Pair-iterator specialization for Bv-vs-all - returning atoms
//!
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename CHAIN_TYPE>
class PairIteratorOneVsAll<CHAIN_TYPE,
                           typename CHAIN_TYPE::ChainTree::BvType,
                           Atom>
     : public PairIteratorOneVsAllBase<PairIteratorOneVsAll<CHAIN_TYPE,
                                                            typename CHAIN_TYPE::ChainTree::BvType,
                                                            Atom>,
                                       CHAIN_TYPE,
                                       typename CHAIN_TYPE::ChainTree::BvType,
                                       Atom> {

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Settings are identical to PairIterator
     typedef typename PairIterator<CHAIN_TYPE,BvType,Atom>::Settings Settings;

     //! Local settings object
     Settings settings;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param entity The fixed element that we are comparing to
     //! \param settings Local settings object
     PairIteratorOneVsAll(CHAIN_TYPE &chain,
                          BvType &entity,
                          const Settings &settings=Settings())
          : PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,BvType,Atom>(chain, entity, settings),
            settings(settings) {

          this->entity2 = &entity;
          this->init();
     }     


     //! Initializer
     void init() {
          
          // If defaultMaxInteractionDistance is set, the nodeIterator
          // should ignore the selected_atoms vector, since we will
          // have to consider all atoms
          if (is_initialized(settings.default_max_interaction_distance)) {
               settings.selected_atoms_inner_iterator = NULL;
          }

          // Call base class initializer
          PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,BvType,Atom>::init();
     }

     //! Check whether node and entity are within residue distance
     //! Called by base class
     //! Since a BV object has no residue information, we always return false
     //!
     //! \param node Current node value
     //! \param entity The fixed element that we are comparing to
     //! \return True if the residues of node and entity are not separated by more than the minimum distance
     bool within_minimum_residue_distance(NodeType *node, BvType *entity) {
          return false;
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {

          this->node_iterator1 = 
               typename Node<Residue, BvType>::Iterator(this->node1, 
                                                         settings.selected_atoms_inner_iterator,
                                                         &settings.selected_atoms_indices);
          bool complete = false;
          if (*this->node_iterator1) {
               this->entity1 = *this->node_iterator1;
               complete = true;
          }
          return complete;

     }

     //! Find next bv-atom pair
     void find_next_pair() {

          // First attempt iteration over atoms in current node
          ++this->node_iterator1;
          if (*this->node_iterator1) {
               this->entity1 = *this->node_iterator1;
               return;
          }

          // If no atom-pairs left in current node pair, find new pair.
          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }


     //! Increment operator
     //!
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator++() {

          // Find next bv-atom pair
          find_next_pair();

          return (*this);
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator=(const PairIteratorOneVsAll& other) {
          PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,BvType,Atom>::operator=(other);
          this->settings = other.settings;
          return *this;
     }
};


//! Pair-iterator specialization for atom-vs-all - returning nodes
//!
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename CHAIN_TYPE>
class PairIteratorOneVsAll<CHAIN_TYPE,
                           typename CHAIN_TYPE::ChainTree::NodeType,
			   Atom>
     : public PairIteratorOneVsAllBase<PairIteratorOneVsAll<CHAIN_TYPE,
                                                            typename CHAIN_TYPE::ChainTree::NodeType,
							    Atom>,
                                       CHAIN_TYPE,
                                       typename CHAIN_TYPE::ChainTree::NodeType,
				       Atom> {

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

public:

     //! Settings are identical to PairIteratorOneVsAllBase
     typedef typename PairIteratorOneVsAllBase<PairIteratorOneVsAll,
                                               CHAIN_TYPE,
                                               NodeType,
					       Atom>::Settings Settings;

     //! Local settings object
     Settings settings;

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param entity The fixed element that we are comparing to
     //! \param settings Local settings object
     PairIteratorOneVsAll(CHAIN_TYPE &chain,
                          Atom &entity,
                          const Settings &settings=Settings())
          : PairIteratorOneVsAllBase<PairIteratorOneVsAll,
                                     CHAIN_TYPE,
                                     NodeType,
				     Atom> (chain, entity, settings),
            settings(settings) {

          this->entity2 = &entity;
          this->init();
     }     

     //! Check whether node and entity are within residue distance
     //! Called by base class
     //!
     //! \param node Current node value
     //! \param entity The fixed element that we are comparing to
     //! \return True if the residues of node and entity are not separated by more than the minimum distance
     bool within_minimum_residue_distance(NodeType *node, Atom *entity) {
          return (fabs(node->frame->res->index - entity->residue->index)<
                  settings.minimum_residue_distance);
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {
          this->entity1 = this->node1;
          return true;
     }

     //! Find next atom-node pair
     void find_next_pair() {

          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }


     //! Increment operator
     //!
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator++() {

          // Find next atom-node pair
          find_next_pair();

          return (*this);
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator=(const PairIteratorOneVsAll& other) {
          PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,NodeType,Atom>::operator=(other);
          this->settings = other.settings;
          return *this;
     }
};


//! Pair-iterator specialization for bv-vs-all - returning nodes
//!
//! \tparam CHAIN_TYPE Type of molecule chain
template <typename CHAIN_TYPE>
class PairIteratorOneVsAll<CHAIN_TYPE,
                           typename CHAIN_TYPE::ChainTree::BvType,
                           typename CHAIN_TYPE::ChainTree::NodeType>
     : public PairIteratorOneVsAllBase<PairIteratorOneVsAll<CHAIN_TYPE,
                                                            typename CHAIN_TYPE::ChainTree::BvType,
                                                            typename CHAIN_TYPE::ChainTree::NodeType>,
                                       CHAIN_TYPE,
                                       typename CHAIN_TYPE::ChainTree::BvType,
                                       typename CHAIN_TYPE::ChainTree::NodeType> {

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Settings are identical to PairIteratorOneVsAllBase::Settings
     typedef typename PairIteratorOneVsAllBase<PairIteratorOneVsAll,
                                               CHAIN_TYPE,
                                               BvType,
                                               NodeType>::Settings Settings;

     //! Local settings object
     Settings settings;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     //! \param entity The fixed element that we are comparing to
     //! \param settings Local settings object
     PairIteratorOneVsAll(CHAIN_TYPE &chain,
                          BvType &entity,
                          const Settings &settings=Settings())
          : PairIteratorOneVsAllBase<PairIteratorOneVsAll,
                                     CHAIN_TYPE,
                                     BvType,
                                     NodeType > (chain, entity, settings),
            settings(settings) {

          this->entity2 = &entity;
          this->init();
     }     

     //! Check whether node and entity are within residue distance
     //! Called by base class
     //! Since a BV object has no residue information, we always return false
     //!
     //! \param node Current node value
     //! \param entity The fixed element that we are comparing to
     //! \return True if the residues of node and entity are not separated by more than the minimum distance
     bool within_minimum_residue_distance(NodeType *node, BvType *entity) {
          return false;
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {
          this->entity1 = this->node1;
          return true;
     }

     //! Find next bv-node pair
     void find_next_pair() {

          bool complete = this->find_next_node_pair();
          if (!complete) {
               this->entity1 = NULL;
               this->entity2 = NULL;
          }
     }


     //! Increment operator
     //!
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator++() {

          // Find next bv-node pair
          find_next_pair();

          return (*this);
     }

     //! Overload assignment operator
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     PairIteratorOneVsAll& operator=(const PairIteratorOneVsAll& other) {
          PairIteratorOneVsAllBase<PairIteratorOneVsAll,CHAIN_TYPE,BvType,NodeType>::operator=(other);
          this->settings = other.settings;
          return *this;
     }
};




//! Base class for Cached version of PairIterator
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
//! \tparam CONTRIBUTION_TYPE Type of value registered to cache
//! \tparam RETURN_VALUE_TYPE Type of value returned as result
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2, 
          typename CONTRIBUTION_TYPE, typename RETURN_VALUE_TYPE>
class CachedIteratorBase
     : public ::phaistos::CachedIteratorBase<CONTRIBUTION_TYPE, RETURN_VALUE_TYPE>, 
       public chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type { 

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;


     //! Stack of inner node-pairs that have been processed (for postprocessing)
     std::stack<chaintree::Node<typename CHAIN_TYPE::Residue, BvType> *> node_pair_stack_post;

     //! Cache datastructure
     //! For each level in the tree, a matrix of interactions between nodes is maintained
     //! each entry is a <value,bool> pair, where the bool indicates if the node is the 
     //! root of a fully initialized subtree
     std::vector<std::vector<std::vector<std::pair<RETURN_VALUE_TYPE,bool> > > > cache;

     //! This is a vector of backup values needed in case of a reject
     //! Each time a cache value is changed, the corresponding node pointers and a <value,bool> 
     //! pair is stored
     std::vector<std::pair<std::pair<NodeType*,NodeType*>,std::pair<RETURN_VALUE_TYPE,bool> > > cache_backup;
     
     //! Sum for current node pair
     RETURN_VALUE_TYPE node_pair_sum;

     //! Total sum
     RETURN_VALUE_TYPE sum;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     CachedIteratorBase(CHAIN_TYPE &chain)
          // : chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>(chain) {
          : chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type(chain) {

          // Allocate cache
          cache.resize(chain.chain_tree->get_height());
          for (unsigned int i=0; i<cache.size(); i++) {
               cache[i].resize(chain.chain_tree->get_nodes_at_level(i));
               for (unsigned int j=0; j<cache[i].size(); j++) {
                    // cache[i][j].resize(j+1, 0);
                    cache[i][j].resize(j+1, std::make_pair(0,true));
               }
          }

          // Initialize sums
          node_pair_sum = 0;
          sum = 0.0;

          // Reset chaintree
          chain.get_chain_tree()->reset();
     }


     //! Overload () operator to act as a constructor (has the same interface as the PairIterator constructor)
     //!
     //! \param chain Molecule chain
     //! \param settings PairIterator settings object
     void operator()(CHAIN_TYPE &chain, 
                     typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::Settings &settings) {

          // Ensure that cache is in sync with chain;
          if (this->enforce_cache_sync(chain))
               // Only resetchaintree if it hasn't just been reset
               // this is to avoid that several energy terms reset the chaintree
               if (this->ct->time != -1)
                    this->ct->reset();

          // Clean backup
          cache_backup.clear();

          // Here we explicitly cast the pair iterator settings object to the cached iterator settings
          // These have the same type, but cannot be assigned normally since the parent class
          // is parameterized with the derived class
          typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::Settings &settings_cache = 
               (typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::Settings&)settings;

          // Assign pair iterator to cached pair iterator
          typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type it(chain,
                                                                                                                                  settings_cache);
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::operator=(it);

          this->init();

     }

     //! Register a contribution to the cached iterator
     //!
     //! \param contribution Value to register
     //! \return reference to cache entry in which value is stored
     RETURN_VALUE_TYPE &register_contribution(const CONTRIBUTION_TYPE &contribution) {
          node_pair_sum += contribution;
          return get_cache_entry(this->node1, this->node2).first;
     }

     //! Compute total value
     //!
     //! \return Total sum
     RETURN_VALUE_TYPE &compute_total() {
          // Update cache for ancestors of nodes that were updated
          // in this iteration
          while (this->node_pair_stack_post.size() > 0) {

               // Pop nodes from post-iteration stack
               NodeType *node1 = this->node_pair_stack_post.top();
               this->node_pair_stack_post.pop();
               NodeType *node2 = this->node_pair_stack_post.top();
               this->node_pair_stack_post.pop();

               RETURN_VALUE_TYPE child_sum(0.0);
               // Intra-node interactions
               if (node1 == node2) {
                    child_sum += get_cache_entry(node1->child1, node1->child1).first;

                    if (node1->child2) {
                         child_sum += get_cache_entry(node1->child1, node1->child2).first;
                         child_sum += get_cache_entry(node1->child2, node1->child2).first;
                    }

               // Inter-node interactions
               } else {
                    child_sum += get_cache_entry(node1->child1, node2->child1).first;

                    if (node1->child2) {
                         child_sum += get_cache_entry(node1->child2, node2->child1).first;

                         if (node2->child2) {
                              child_sum += get_cache_entry(node1->child1, node2->child2).first;
                              child_sum += get_cache_entry(node1->child2, node2->child2).first;
                         }
                    } else if (node2->child2) {
                         child_sum += get_cache_entry(node1->child1, node2->child2).first;
                    }
               }
               set_cache_entry(node1, node2, child_sum);
          }

          // The total sum is given as the root node's interaction with itself
          NodeType *root = this->ct->nodes[this->ct->nodes.size()-1];               

          sum = cache[root->level][root->index][root->index].first;
          
          return sum;
     }     

     //! Accept last evaluation
     void accept() {
          // Clear backup
          cache_backup.clear();
          this->time_stamp++;
     }

     //! Reject last evaluation
     void reject() {
          for (unsigned int i=0; i<cache_backup.size(); ++i) {
               NodeType *node1 = cache_backup[i].first.first;
               NodeType *node2 = cache_backup[i].first.second;
               RETURN_VALUE_TYPE &value = cache_backup[i].second.first;
               bool fully_initialized_subtree = cache_backup[i].second.second;
               cache[node1->level][node2->index][node1->index].first = value;
               cache[node1->level][node2->index][node1->index].second = fully_initialized_subtree;
          }
          // Clear backup
          cache_backup.clear();
     }


     //! Set an entry in the cache
     //!
     //! \param node1 First node
     //! \param node2 Second node
     //! \param value Value to set in cache
     //! \param fully_initialized_subtree Whether subtree with this node-pair as root is fully initialized
     void set_cache_entry(NodeType *node1,
                          NodeType *node2,
                          const RETURN_VALUE_TYPE &value,
                          bool fully_initialized_subtree=true) {
          int level = node1->level;
          assert(node1->index <= node2->index);
          cache_backup.push_back(std::make_pair(std::make_pair(node1,node2), 
                                                std::make_pair(cache[level][node2->index][node1->index].first,
                                                               cache[level][node2->index][node1->index].second)));
          cache[level][node2->index][node1->index].first = value;
          cache[level][node2->index][node1->index].second = fully_initialized_subtree;
     }

     //! Get an entry from the cache
     //!
     //! \param node1 First node
     //! \param node2 Second node
     //! \return (value,bool) pair, where the bool indicates if the node is the 
     //! root of a fully initialized subtree
     std::pair<RETURN_VALUE_TYPE,bool> &get_cache_entry(NodeType *node1,
                                                        NodeType *node2) {
          int level = node1->level;

          assert(node1->index <= node2->index);

          return cache[level][node2->index][node1->index];
     }

     //! For debuggin purposes: print the current contents of the cache
     void print_cache_contents() {

          NodeType *root = this->ct->nodes[this->ct->nodes.size()-1];               

          int start_index = 0;
          int end_index = 0;
          for (int level =0; level < root->level; ++level) {
               start_index = end_index;
               end_index = start_index + this->ct->get_nodes_at_level(level);

               std::cout << "LEVEL: " << level << "\n";

               for (int j1=start_index; j1<end_index; ++j1) {
                    NodeType *node1 = this->ct->nodes[j1];
                    for (int j2=start_index; j2<=j1; ++j2) {
                         NodeType *node2 = this->ct->nodes[j2];

                         if (node1->index > node2->index) {
                              if (cache[level][node1->index][node2->index].first > 0)
                                   std::cout << "cache entry: " << node1->id << " " << node2->id << ": " << cache[level][node1->index][node2->index].first << "(" << cache[level][node1->index][node2->index].second << ")\n";
                         } else {

                              if (cache[level][node2->index][node1->index].first > 0)
                                   std::cout << "cache entry: " << node1->id << " " << node2->id << ": " << cache[level][node2->index][node1->index].first << "(" << cache[level][node2->index][node1->index].second << ")\n";

                         }
                    }
               }
          }
     }     


     //// Call backs ////

     //! Code executed when node pair is excluded due to distance cutoff
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void on_node_pair_distance_exceeds_cutoff(NodeType *node1, NodeType *node2) {

          // When a node pair is excluded, the whole subtree is pruned
          // we do not actually update the cache if the entire subtree (to zero)
          // but instead just mark the subtree as not fully initialized
          bool fully_initialized_subtree = false;
          set_cache_entry(node1, node2, 0, fully_initialized_subtree);
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_node_pair_distance_exceeds_cutoff(node1,node2);
     }
     
     //! Code executed when potential child nodes are considered
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void on_potential_child_nodes(NodeType *node1, NodeType *node2) {

          // Ensure that cache value is set to zero
          // The chaintree cache uses a lazy initialization
          // when subtrees are pruned: only the root node 
          // is marked in the pruned subtree.
          // When traversing the tree, if a node-pair below the
          // root is used, we must ensure that the zero is passed down.
          NodeType *parent1 = node1->parent;
          NodeType *parent2 = node2->parent;
          if (!cache[parent1->level][parent2->index][parent1->index].second) {
               bool fully_initialized_subtree = false;
               set_cache_entry(node1, node2, 0, fully_initialized_subtree);
          }
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_potential_child_nodes(node1,node2);
     }

     //! Code executed when all pairs within leaf pair have been iterated over
     void on_leaf_node_pair_end() {
          if (this->node1!=NULL) {
               // Update cache with node_pair_sum
               set_cache_entry(this->node1, this->node2, node_pair_sum);
          }
          node_pair_sum = 0.0;
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_leaf_node_pair_end();
     }

     //! Code executed when new inner node pair is iterated over.
     //! pushes pair on stack for post iteration traversal
     void on_inner_node_pair_begin() {
          node_pair_stack_post.push(this->node2);
          node_pair_stack_post.push(this->node1);
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_inner_node_pair_begin();
     }

};

//! Base class for Cached version of PairIterator - with vector as each cache entry
//!
//! \tparam DERIVED_CLASS Type of the class inheriting from this base class
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
//! \tparam CONTRIBUTION_TYPE Type of value registered to cache
//! \tparam VECTOR_RETURN_VALUE_TYPE Type of value returned as result
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2, 
          typename CONTRIBUTION_TYPE, typename VECTOR_RETURN_VALUE_TYPE>
class CachedIteratorVectorBase
     : public ::phaistos::CachedIteratorBase<std::pair<CONTRIBUTION_TYPE,CONTRIBUTION_TYPE>,
                                             VECTOR_RETURN_VALUE_TYPE>, 
       public chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type { 

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

protected:

     //! Stack of inner node-pairs that have been processed (for postprocessing)
     std::stack<chaintree::Node<typename CHAIN_TYPE::Residue, BvType> *> node_pair_stack_post;

     //! Cache datastructure
     //! For each level in the tree, a matrix of interactions between nodes is maintained
     //! each entry is a <value,bool> pair, where the bool indicates if the node is the 
     //! root of a fully initialized subtree
     std::vector<std::vector<std::vector<std::pair<VECTOR_RETURN_VALUE_TYPE,bool> > > > cache;

     //! Datastructure keeping track of where the nodes in a child pair are located in a parent pair
     std::vector<std::vector<std::vector<std::vector<int> > > > node_pair_parent_indices;

     //! This is a vector of backup values needed in case of a reject
     //! Each time a cache value is changed, the corresponding node pointers and a <value,bool> 
     //! pair is stored
     std::vector<std::pair<std::pair<NodeType*,NodeType*>,std::pair<VECTOR_RETURN_VALUE_TYPE,bool> > > cache_backup;
     
     //! Sum for current node pair
     VECTOR_RETURN_VALUE_TYPE node_pair_sum;

     //! Initializer
     //!
     //! \param chain Molecule chain
     void initialize(CHAIN_TYPE &chain) {
          // Allocate cache
          cache.resize(chain.chain_tree->get_height());
          node_pair_parent_indices.resize(chain.chain_tree->get_height());
          unsigned int offset = 0;
          for (unsigned int i=0; i<cache.size(); ++i) {
               cache[i].resize(chain.chain_tree->get_nodes_at_level(i));
               node_pair_parent_indices[i].resize(chain.chain_tree->get_nodes_at_level(i));
               for (unsigned int j=0; j<cache[i].size(); ++j) {
                    cache[i][j].resize(j+1);
                    node_pair_parent_indices[i][j].resize(j+1);
               }
          }

          // Initialize node_pair_parent_indices - first phase
          // Initially the indices are initialized to global indices
          for (unsigned int i=0; i<cache.size(); ++i) {
               for (unsigned int j=0; j<cache[i].size(); ++j) {
                    for (unsigned int k=0; k<cache[i][j].size(); ++k) {

                         NodeType *node1 = chain.chain_tree->nodes[offset+j];
                         NodeType *node2 = chain.chain_tree->nodes[offset+k];
                         if (node1->index < node2->index) {
                              node2 = chain.chain_tree->nodes[offset+j];
                              node1 = chain.chain_tree->nodes[offset+k];
                         }
                         unsigned int index1 = node1->index;
                         unsigned int index2 = node2->index;

                         if (i==0) { // leaf-node

                              int leafnodes_in_subtree = 2;
                              node_pair_parent_indices[i][index1][index2].resize(leafnodes_in_subtree);
                              node_pair_parent_indices[i][index1][index2][0] = index1;
                              node_pair_parent_indices[i][index1][index2][1] = index2;

                         } else {

                              VECTOR_RETURN_VALUE_TYPE atoms_in_leafnodes;
                              if (j==k) {

                                   NodeType *child1 = node1->child1;
                                   NodeType *child2 = node1->child2;

                                   std::vector<int> leaf_indices_child1 = node_pair_parent_indices[i-1][child1->index][child1->index];

                                   node_pair_parent_indices[i][index1][index2]
                                        .insert(node_pair_parent_indices[i][index1][index2].end(),
                                                leaf_indices_child1.begin(),
                                                leaf_indices_child1.end());   
                                   if (child2) {
                                        std::vector<int> leaf_indices_child2 = node_pair_parent_indices[i-1][child2->index][child2->index];
                                        node_pair_parent_indices[i][index1][index2]
                                             .insert(node_pair_parent_indices[i][index1][index2].end(),
                                                     leaf_indices_child2.begin(),
                                                     leaf_indices_child2.end());   
                                   }

                              } else {
                                   NodeType *child11 = node1->child1;
                                   NodeType *child12 = node1->child2;
                                   NodeType *child21 = node2->child1;
                                   NodeType *child22 = node2->child2;

                                   std::vector<int> leaf_indices_child11 = node_pair_parent_indices[i-1][child11->index][child11->index];
                                   node_pair_parent_indices[i][index1][index2]
                                        .insert(node_pair_parent_indices[i][index1][index2].end(),
                                                leaf_indices_child11.begin(),
                                                leaf_indices_child11.end());                                      
                                   std::vector<int> leaf_indices_child21 = node_pair_parent_indices[i-1][child21->index][child21->index];
                                   node_pair_parent_indices[i][index1][index2]
                                        .insert(node_pair_parent_indices[i][index1][index2].end(),
                                                leaf_indices_child21.begin(),   
                                                leaf_indices_child21.end());                                      
                                   if (child12) {
                                        std::vector<int> leaf_indices_child12 = node_pair_parent_indices[i-1][child12->index][child12->index];
                                        node_pair_parent_indices[i][index1][index2]
                                             .insert(node_pair_parent_indices[i][index1][index2].end(),
                                                     leaf_indices_child12.begin(),
                                                     leaf_indices_child12.end());                                      
                                   }
                                   if (child22) {
                                        std::vector<int> leaf_indices_child22 = node_pair_parent_indices[i-1][child22->index][child22->index];
                                        node_pair_parent_indices[i][index1][index2]
                                             .insert(node_pair_parent_indices[i][index1][index2].end(),
                                                     leaf_indices_child22.begin(),   
                                                     leaf_indices_child22.end());                                      
                                   }                                   
                              }
                         }

                         // Sort leaf indices
                         std::sort(node_pair_parent_indices[i][index1][index2].begin(), 
                                   node_pair_parent_indices[i][index1][index2].end());
                                   
                         // // Remove duplicates
                         node_pair_parent_indices[i][index1][index2].erase(std::unique(node_pair_parent_indices[i][index1][index2].begin(),
                                                                                     node_pair_parent_indices[i][index1][index2].end()),
                                                                         node_pair_parent_indices[i][index1][index2].end());

                         // std::vector<unsigned int> sizes(node_pair_parent_indices[i][index1][index2].size());
                         // for (unsigned l=0; l<sizes.size(); ++l) {
                         //      unsigned int leaf_index = node_pair_parent_indices[i][index1][index2][l];
                         //      sizes[l] = chain.chain_tree->nodes[leaf_index]->size();
                         // }
                         // ((DERIVED_CLASS*)this)->allocate_cache_entry(cache[i][index1][index2].first, sizes);

                         // VECTOR_RETURN_VALUE_TYPE cache_vector(node_pair_parent_indices[i][index1][index2].size());
                         // for (unsigned l=0; l<cache_vector.size(); ++l) {
                         //      unsigned int leaf_index = node_pair_parent_indices[i][index1][index2][l];
                         //      cache_vector[l] = 

                         cache[i][index1][index2].first.resize(node_pair_parent_indices[i][index1][index2].size());
                         for (unsigned l=0; l<node_pair_parent_indices[i][index1][index2].size(); ++l) {
                              unsigned int leaf_index = node_pair_parent_indices[i][index1][index2][l];
                              cache[i][index1][index2].first[l].initialize(chain.chain_tree->nodes[leaf_index]);
                         }
                         //      // cache[i][index1][index2].first[l].resize(chain.chain_tree->nodes[leaf_index]->atoms.size());
                         //      // RETURN_VALUE_TYPE tmp(chain.chain_tree->nodes[leaf_index]);

                         //      ((DERIVED_CLASS*)this)->initialize_cache_entry(cache[i][index1][index2].first[l], chain.chain_tree->nodes[leaf_index]);
                         //      // cache[i][index1][index2].first[l] = RETURN_VALUE_TYPE(chain.chain_tree->nodes[leaf_index]);

                         //      // cache[i][index1][index2].first[l] = RETURN_VALUE_TYPE(tmp);
                         //      // cache[i][index1][index2].first[l].initialize(*chain.chain_tree->nodes[leaf_index]);
                         // }
                    }
               }
               offset += chain.chain_tree->get_nodes_at_level(i);
          }

          // Initialize node_pair_parent_indices - second phase
          // The indices are now reinitialized so that children now at which indices there nodes reside in their parents
          offset = 0;
          for (unsigned int i=0; i<cache.size()-1; ++i) {

               if (i==cache.size()-1) {
                    continue;
               }

               for (unsigned int j=0; j<cache[i].size(); ++j) {
                    for (unsigned int k=0; k<cache[i][j].size(); ++k) {

                         NodeType *node1 = chain.chain_tree->nodes[offset+j];
                         NodeType *node2 = chain.chain_tree->nodes[offset+k];
                         if (node1->index < node2->index) {
                              node2 = chain.chain_tree->nodes[offset+j];
                              node1 = chain.chain_tree->nodes[offset+k];
                         }
                         unsigned int index1 = node1->index;
                         unsigned int index2 = node2->index;

                         unsigned int parent_index1 = node1->parent->index;
                         unsigned int parent_index2 = node2->parent->index;

                         for (unsigned int l=0; l<node_pair_parent_indices[i][index1][index2].size(); ++l) {
                              int &leaf_node_index = node_pair_parent_indices[i][index1][index2][l];
                              bool found = false;
                              for (unsigned int m=0; m<node_pair_parent_indices[i+1][parent_index1][parent_index2].size(); ++m) {
                                   if (node_pair_parent_indices[i+1][parent_index1][parent_index2][m] == leaf_node_index) {
                                        leaf_node_index = m;
                                        found=true;
                                        break;
                                   }
                              }
                              assert(found);
                         }
                    }
               }
               offset += chain.chain_tree->get_nodes_at_level(i);
          }

          // Reset chaintree
          chain.get_chain_tree()->reset();          
     }

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     CachedIteratorVectorBase(CHAIN_TYPE &chain)
          : chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type(chain) {

          initialize(chain);
     }


     //! Overload assignment operator.
     //! This method is inherited automatically by GCC compilers, but needed explicitly by ICC.
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     DERIVED_CLASS& operator=(const DERIVED_CLASS& other) {
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::operator=(other);
          return (*((DERIVED_CLASS*)this));
     }

     //! Register a contribution to the cached iterator
     //!
     //! \param contribution Value to register
     //! \return reference to cache entry in which value is stored
     VECTOR_RETURN_VALUE_TYPE &register_contribution(const std::pair<CONTRIBUTION_TYPE,
                                                                     CONTRIBUTION_TYPE> &contribution) {
          bool potential_duplicate = false;
          if (this->node1->index == this->node2->index) {
               node_pair_sum[0].register_contribution(this->node1, (*this)->first, (*this)->second, 
                                                      contribution.first, potential_duplicate);

               potential_duplicate = true;
               node_pair_sum[0].register_contribution(this->node2, (*this)->second, (*this)->first, 
                                                      contribution.second, potential_duplicate);
          } else {
               node_pair_sum[0].register_contribution(this->node1, (*this)->first, (*this)->second, 
                                                      contribution.first, potential_duplicate);
               node_pair_sum[1].register_contribution(this->node2, (*this)->second, (*this)->first, 
                                                      contribution.second, potential_duplicate);
          }

          return get_cache_entry(this->node1, this->node2).first;
     }


     //! Copy child value entries in vector to corresponding parent values
     //!
     //! \param child1 Child node 1
     //! \param child2 Child node 2
     //! \param parent_pair_vector Destination vector
     void copy_child_pair_to_parent_pair(NodeType *child1, NodeType *child2,
                                         VECTOR_RETURN_VALUE_TYPE &parent_pair_vector) {
          // Check if child pair is flagged as a zero contribution
          if (!get_cache_entry(child1, child2).second) {
               return;
          }
          std::vector<int> &child_parent_node_indices = 
               node_pair_parent_indices[child1->level][child2->index][child1->index];
          VECTOR_RETURN_VALUE_TYPE &child_pair_vector = 
               get_cache_entry(child1, child2).first;
          for (unsigned int i=0; i<child_parent_node_indices.size(); ++i) {
               unsigned int parent_index = child_parent_node_indices[i];
               parent_pair_vector[parent_index] += child_pair_vector[i];
               // for (unsigned int j=0; j<child_pair_vector[i].size(); ++j) {
               //      parent_pair_vector[parent_index][j] += child_pair_vector[i][j];
               // }
          }
     }

     //! Compute total value
     //!
     //! \return Total sum
     VECTOR_RETURN_VALUE_TYPE &compute_total() {

          // Update cache for ancestors of nodes that were updated
          // in this iteration
          while (this->node_pair_stack_post.size() > 0) {

               // Pop nodes from post-iteration stack
               NodeType *node1 = this->node_pair_stack_post.top();
               this->node_pair_stack_post.pop();
               NodeType *node2 = this->node_pair_stack_post.top();
               this->node_pair_stack_post.pop();

               initialize_node_pair_sum(node1, node2);

               // Intra-node interactions
               if (node1 == node2) {

                    // Copy from child1-child1
                    ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child1, node1->child1, node_pair_sum);

                    if (node1->child2) {
                         // Copy from child1-child2
                         ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child1, node1->child2, node_pair_sum);

                         // Copy from child2-child2
                         ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child2, node1->child2, node_pair_sum);
                    }
               // Inter-node interactions
               } else {

                    // Copy from child1-child1
                    ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child1, node2->child1, node_pair_sum);

                    if (node1->child2) {
                         // Copy from child2-child1
                         ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child2, node2->child1, node_pair_sum);

                         if (node2->child2) {
                              // Copy from child1-child2
                              ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child1, node2->child2, node_pair_sum);

                              // Copy from child-child2
                              ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child2, node2->child2, node_pair_sum);
                         } 
                    } else if (node2->child2) {
                         // Copy from child1-child2
                         ((DERIVED_CLASS*)this)->copy_child_pair_to_parent_pair(node1->child1, node2->child2, node_pair_sum);
                    }
               }

               ((DERIVED_CLASS*)this)->set_cache_entry(node1, node2, node_pair_sum);
          }

          // The total sum is given as the root node's interaction with itself
          NodeType *root = this->ct->nodes[this->ct->nodes.size()-1];               

          // sum = cache[root->level][root->index][root->index].first;
          
          // return sum;
          return cache[root->level][root->index][root->index].first;
     }


     //! Accept last evaluation
     void accept() {
          // Clear backup
          cache_backup.clear();
          this->time_stamp++;
     }

     //! Reject last evaluation
     void reject() {
          for (unsigned int i=0; i<cache_backup.size(); ++i) {
               NodeType *node1 = cache_backup[i].first.first;
               NodeType *node2 = cache_backup[i].first.second;
               VECTOR_RETURN_VALUE_TYPE &value = cache_backup[i].second.first;
               bool fully_initialized_subtree = cache_backup[i].second.second;
               if (value.size() > 0) {
                    cache[node1->level][node2->index][node1->index].first = value;
               }
               cache[node1->level][node2->index][node1->index].second = fully_initialized_subtree;
          }
          // Clear backup
          cache_backup.clear();
     }

     //! Overload () operator to act as a constructor (has the same interface as the PairIterator constructor)
     //!
     //! \param chain Molecule chain
     //! \param settings PairIterator settings object
     void operator()(CHAIN_TYPE &chain, 
                     typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::Settings &settings) {

          // Ensure that cache is in sync with chain;
          if (this->enforce_cache_sync(chain))
               this->ct->reset();

          // Clean backup
          cache_backup.clear();

          // Here we explicitly cast the pair iterator settings object to the cached iterator settings
          // These have the same type, but cannot be assigned normally since the parent class
          // is parameterized with the derived class
          typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template 
               Implementation<DERIVED_CLASS>::Type::Settings &settings_cache = 
               (typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template 
                Implementation<DERIVED_CLASS>::Type::Settings&)settings;

          // Assign pair iterator to cached pair iterator
          typename chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template 
               Implementation<DERIVED_CLASS>::Type it(chain, settings_cache);
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template 
               Implementation<DERIVED_CLASS>::Type::operator=(it);

          this->init();
     }

     //! Get an entry from the cache
     //!
     //! \param node1 First node
     //! \param node2 Second node
     //! \return (value,bool) pair, where the bool indicates if the node is the 
     //! root of a fully initialized subtree
     std::pair<VECTOR_RETURN_VALUE_TYPE,bool> &get_cache_entry(NodeType *node1,
                                                               NodeType *node2) {
          int level = node1->level;

          assert(node1->index <= node2->index);

          return cache[level][node2->index][node1->index];
     }

     //! Set an entry in the cache
     //!
     //! \param node1 First node
     //! \param node2 Second node
     //! \param value Value to set in cache
     //! \param fully_initialized_subtree Whether subtree with this node-pair as root is fully initialized
     void set_cache_entry(NodeType *node1,
                          NodeType *node2,
                          const VECTOR_RETURN_VALUE_TYPE &value,
                          bool fully_initialized_subtree=true) {
          int level = node1->level;
          assert(node1->index <= node2->index);
          cache_backup.push_back(std::make_pair(std::make_pair(node1,node2), 
                                                std::make_pair(cache[level][node2->index][node1->index].first,
                                                               cache[level][node2->index][node1->index].second)));
          cache[level][node2->index][node1->index].first = value;
          cache[level][node2->index][node1->index].second = fully_initialized_subtree;
     }

     //! Initialize the sum vector for the current node pair
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void initialize_node_pair_sum(NodeType *node1, NodeType *node2) {
          int local_cache_size = get_cache_entry(node1, node2).first.size();
          node_pair_sum = VECTOR_RETURN_VALUE_TYPE(local_cache_size);
          for (unsigned int i=0; i<node_pair_sum.size(); ++i) {
               node_pair_sum[i].initialize_from_other(get_cache_entry(node1, node2).first[i]);
               // node_pair_sum[i].resize(get_cache_entry(node1, node2).first[i].size());
          }
     }


     //// Call-backs ////

     //! Code executed when node pair is excluded due to distance cutoff
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void on_node_pair_distance_exceeds_cutoff(NodeType *node1, NodeType *node2) {

          // When a node pair is excluded, the whole subtree is pruned
          // we do not actually update the cache if the entire subtree (to zero)
          // but instead just mark the subtree as not fully initialized
          bool fully_initialized_subtree = false;
          cache_backup.push_back(std::make_pair(std::make_pair(node1,node2), 
                                                std::make_pair(VECTOR_RETURN_VALUE_TYPE(),
                                                               cache[node1->level][node2->index][node1->index].second)));
          get_cache_entry(node1, node2).second = fully_initialized_subtree;
     }
     
     //! Code executed when potential child nodes are considered
     //!
     //! \param node1 First node
     //! \param node2 Second node
     void on_potential_child_nodes(NodeType *node1, NodeType *node2) {

          // Ensure that cache value is set to zero
          // The chaintree cache uses a lazy initialization
          // when subtrees are pruned: only the root node 
          // is marked in the pruned subtree.
          // When traversing the tree, if a node-pair below the
          // root is used, we must ensure that the zero is passed down.
          NodeType *parent1 = node1->parent;
          NodeType *parent2 = node2->parent;
          if (!cache[parent1->level][parent2->index][parent1->index].second) {
               bool fully_initialized_subtree = false;
               get_cache_entry(node1, node2).second = fully_initialized_subtree;
          }
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_potential_child_nodes(node1,node2);
     }

     //! Code executed when all pairs within leaf pair have been iterated over
     void on_leaf_node_pair_end() {
          if (this->node1!=NULL) {
               ((DERIVED_CLASS*)this)->set_cache_entry(this->node1, this->node2, node_pair_sum);
          }

          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_leaf_node_pair_end();
     }

     //! Code executed when base class iteration hits a pair of leaves
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_pair_begin() {
          initialize_node_pair_sum(this->node1, this->node2);
          return chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_leaf_node_pair_begin();
     }

     //! Code executed when base class iteration hits a single leaf
     //!
     //! \return False if leaf iteration is finished
     bool on_leaf_node_single_begin() {
          initialize_node_pair_sum(this->node1, this->node2);
          return chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_leaf_node_single_begin();
     }

     //! Code executed when new inner node pair is iterated over
     //! push pair on stack for post iteration traversal
     void on_inner_node_pair_begin() {
          node_pair_stack_post.push(this->node2);
          node_pair_stack_post.push(this->node1);
          chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>::template Implementation<DERIVED_CLASS>::Type::on_inner_node_pair_begin();
     }


     //! For debugging purposes: print the current contents of the cache
     void print_cache_contents() {

          NodeType *root = this->ct->nodes[this->ct->nodes.size()-1];               

          int start_index = 0;
          int end_index = 0;
          for (int level =0; level <= root->level; ++level) {
               start_index = end_index;
               end_index = start_index + this->ct->get_nodes_at_level(level);

               std::cout << "LEVEL: " << level << "\n";

               for (int j1=start_index; j1<end_index; ++j1) {
                    NodeType *node1 = this->ct->nodes[j1];
                    for (int j2=start_index; j2<=j1; ++j2) {
                         NodeType *node2 = this->ct->nodes[j2];

                         if (node1->index > node2->index) {
                              if (cache[level][node1->index][node2->index].first.size() > 0)
                                   std::cout << "cache entry: " << node1->id << " " << node2->id << ": " << cache[level][node1->index][node2->index].first << "(" << cache[level][node1->index][node2->index].second << ")\n";
                         } else {

                              if (cache[level][node2->index][node1->index].first.size() > 0)
                                   std::cout << "cache entry: " << node1->id << " " << node2->id << ": " << cache[level][node2->index][node1->index].first << "(" << cache[level][node2->index][node1->index].second << ")\n";

                         }
                    }
               }
          }          
     }     
};


}


//! Specialization of CachedIterator for PairIterator
//!
//! \tparam CHAIN_TYPE Type of molecule chain
//! \tparam ENTITY_TYPE1 Type of first element in pairs
//! \tparam ENTITY_TYPE2 Type of second element in pairs
//! \tparam CONTRIBUTION_TYPE Type of value registered to cache
//! \tparam RETURN_VALUE_TYPE Type of value returned as result
template <typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2, typename CONTRIBUTION_TYPE, typename RETURN_VALUE_TYPE>
class CachedIterator<chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>,CONTRIBUTION_TYPE, RETURN_VALUE_TYPE>
     : public chaintree::CachedIteratorBase<CachedIterator<chaintree::PairIterator<CHAIN_TYPE,
                                                                                   ENTITY_TYPE1,
                                                                                   ENTITY_TYPE2>,
                                                         CONTRIBUTION_TYPE, RETURN_VALUE_TYPE>,
                                          CHAIN_TYPE,
                                          ENTITY_TYPE1,
                                          ENTITY_TYPE2,
                                          CONTRIBUTION_TYPE,
                                          RETURN_VALUE_TYPE> { 

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     CachedIterator(CHAIN_TYPE &chain)
          : chaintree::CachedIteratorBase<CachedIterator,
                                          CHAIN_TYPE,
                                          ENTITY_TYPE1,
                                          ENTITY_TYPE2,
                                          CONTRIBUTION_TYPE,
                                          RETURN_VALUE_TYPE>(chain) {
     }     

     //! Overload assignment operator.
     //! This method is inherited automatically by GCC compilers, but needed explicitly by ICC.
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     CachedIterator& operator=(const CachedIterator& other) {
          chaintree::CachedIteratorBase<CachedIterator<chaintree::PairIterator<CHAIN_TYPE,
                                                                               ENTITY_TYPE1,
                                                                               ENTITY_TYPE2>,
                                                       CONTRIBUTION_TYPE, RETURN_VALUE_TYPE>,
                                        CHAIN_TYPE,
                                        ENTITY_TYPE1,
                                        ENTITY_TYPE2,
                                        CONTRIBUTION_TYPE,
                                        RETURN_VALUE_TYPE>::operator=(other);
          return *this;
     }

};


//! Specialization of CachedIterator for Vector-cache PairIterator
template <typename CHAIN_TYPE, typename ENTITY_TYPE1, typename ENTITY_TYPE2, 
          typename CONTRIBUTION_TYPE, typename RETURN_VALUE_TYPE>
class CachedIterator<chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>,
                     std::pair<CONTRIBUTION_TYPE,CONTRIBUTION_TYPE>,
                     std::vector<RETURN_VALUE_TYPE> >
     : public chaintree::CachedIteratorVectorBase<CachedIterator<chaintree::PairIterator<CHAIN_TYPE,ENTITY_TYPE1,ENTITY_TYPE2>,
                                                      std::pair<CONTRIBUTION_TYPE,CONTRIBUTION_TYPE>,
                                                      std::vector<RETURN_VALUE_TYPE> >,
                                                CHAIN_TYPE, ENTITY_TYPE1, ENTITY_TYPE2, 
                                                CONTRIBUTION_TYPE, std::vector<RETURN_VALUE_TYPE> > { 

     //! Define BvType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::BvType BvType;     

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

public:

     //! Constructor
     //!
     //! \param chain Molecule chain
     CachedIterator(CHAIN_TYPE &chain)
          : chaintree::CachedIteratorVectorBase<CachedIterator,
                                                CHAIN_TYPE,
                                                ENTITY_TYPE1,
                                                ENTITY_TYPE2,
                                                CONTRIBUTION_TYPE, std::vector<RETURN_VALUE_TYPE> >(chain) {
     }

     //! Increment operator
     //!
     //! \return Current iterator (this)
     CachedIterator &operator++() {
          // Call increment of PairIterator, but with ourselves as 
          // call-back class
          this->increment();
          return *this;
     }


     //! Overload assignment operator.
     //! This method is inherited automatically by GCC compilers, but needed explicitly by ICC.
     //!
     //! \param other Source object from which assignment is made.
     //! \return Current iterator (this)
     CachedIterator& operator=(const CachedIterator& other) {
          chaintree::CachedIteratorVectorBase<CachedIterator,
                                              CHAIN_TYPE,
                                              ENTITY_TYPE1,
                                              ENTITY_TYPE2,
                                              CONTRIBUTION_TYPE, std::vector<RETURN_VALUE_TYPE> >::operator=(other);
          return *this;
     }
};


//! Data structure which can be used as RETURN_TYPE in 
//! vector version of CachedIterator
template <typename CHAIN_TYPE, typename CONTRIBUTION_TYPE>
class VectorReturnType {

protected:

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

     //! Internal data container
     std::vector<CONTRIBUTION_TYPE> data;

public:

     //! Return size of vector
     //!
     //! \return Size of internal vector
     unsigned size() const {
          return data.size();
     }

     //! Define how data should be initialized given another value
     //! When this method is called, information about which leaf node this value corresponds to
     //! is not known. However, we do know which entry in the cache this value will be written
     //! to (which has the right size), and we thus initialize based on this value
     //!
     //! \param other Source value
     void initialize_from_other(const VectorReturnType &other) {
          data.resize(other.data.size());
     }

     //! Overload += operator
     //! Specifies how new values are added to an VectorReturnType
     //! \param other Object with elements to be added to current object
     VectorReturnType &operator+=(const VectorReturnType &other) {
          for (unsigned int i=0; i<other.data.size(); ++i) {
               data[i] += other.data[i];
          }
          return (*this);
     }


     // These functions are not used by pair iterator, but provided for convenience

     //! Overload [] operator
     //!
     //! \param index index into vector
     //! \return vector element
     CONTRIBUTION_TYPE operator[](const unsigned int index) const {
          return data[index];
     }

     // Overload [] operator
     //!
     //! \param index index into vector
     //! \return vector element
     CONTRIBUTION_TYPE &operator[](const unsigned int index) {
          return data[index];
     }

     //! Overload output operator
     friend std::ostream & operator<<(std::ostream &o, const VectorReturnType &av) {
          o << av.data;
          return o;
     }
};


//! Data structure which can be used as RETURN_TYPE in 
//! vector version of CachedIterator.
//! Keeps track of an atom vector for each node in the tree
template <typename CHAIN_TYPE, typename CONTRIBUTION_TYPE>
class AtomVectorReturnType: public VectorReturnType<CHAIN_TYPE,
                                                    CONTRIBUTION_TYPE> {

     //! Define NodeType locally for ease of reference
     typedef typename CHAIN_TYPE::ChainTree::NodeType NodeType;

public:

     //! Initialize based on node
     //!
     //! \param node Node value
     void initialize(NodeType *node) {
          this->data.resize(node->size());
     }

     //! Specifies how a contribution is stored in this type - given
     //! the current node in the chaintree and the entity type
     //!
     //! \param node Node value
     //! \param atom_self Atom for which contribution is stored
     //! \param atom_other Not used here
     //! \param contribution Value to register in cache
     //! \param potential_duplicate Not used here
     void register_contribution(const NodeType *node, 
                                const Atom *atom_self, 
                                const Atom *atom_other, 
                                const CONTRIBUTION_TYPE &contribution,
                                bool potential_duplicate) {
          this->data[node->get_atom_index(atom_self->atom_type)] += contribution;
     }
};

}

#endif
