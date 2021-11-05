// chain_tree.h --- Data structure to efficiently detect pairwise interactions
//                 Based on: Lotan, Schwarzer, Halperin, Latombe,
//                           Journal of Computational Biology, 2004 
// Copyright (C) 2008 Wouter Boomsma
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


#ifndef CHAINTREE_H
#define CHAINTREE_H

#include <stack>
#include <set>
#include "utils/matrix.h"
#include "moves/move_info.h"
#include "iterators/covalent_bond_iterator.h"
#include "utils/eigen_system_3x3.h"

namespace phaistos {
namespace chaintree {

// Forward declarations
bool rect_distance_within_limit(Vector_3D &translation_AB, Matrix_3D &rotation_AB,
                                double *side_lengths_A, double *side_lengths_B, double cutoff_distance);
double rect_distance(Vector_3D &translation_AB, Matrix_3D &rotation_AB,
                     double *side_lengths_A, double *side_lengths_B);
double rect_distance(Vector_3D &translation_AB, Matrix_3D &rotation_AB, double *side_lengths_A);

bool edge_pair_check(Vector_3D &shortestVector,
                     double edge_A_ref_B_vertex1, double edge_A_ref_B_vertex2,
                     double edge_B_ref_A_vertex1, double edge_B_ref_A_vertex2,
                     double edge_A_ref_A, int edge_A_direction, double edge_A_length,
                     double edge_B_ref_B, int edge_B_direction, double edge_B_length,
                     int edge_A_dot_X,

                     double edge_A_dot_edge_B, double edge_A_dot_edge_Bnorm,
                     double edge_A_norm_dot_edge_B, double edge_A_norm_dot_edge_Bnorm,

                     double translation_AB_dot_edge_A, double translation_AB_dot_edge_A_norm,
                     double translation_BA_dot_edge_B, double translation_BA_dot_edge_B_norm,
                     const Vector_3D &translation_AB, const Vector_3D &edge_B_ref_A,
                     const Vector_3D &edge_B_norm_ref_A);


//! Frame object
//! Each node in the tree is associated with a frame object
//! which defines its local reference frame.
template <typename RES_TYPE>
class Frame {
public:

     //! The frame is defined from backbone atoms in a specified residue
     RES_TYPE *res;

     //! Orientation relative to global reference
     Matrix_3D *orientation;

     //! Origin relative to global reference
     Vector_3D *origin;

     //! Backup: Orientation relative to global reference
     Matrix_3D *orientation_backup;

     //! Backup: Origin relative to global reference
     Vector_3D *origin_backup;

     //! Constructor
     //!
     //! \param res Residue pointer
     Frame(RES_TYPE *res)
          : res(res),
            orientation(NULL), origin(NULL),
            orientation_backup(NULL), origin_backup(NULL) {
          update();
     }

     //! Destructor
     ~Frame() {
          delete orientation;
          delete origin;
          delete orientation_backup;
          delete origin_backup;          
     }

     //! Update frame based on information in residue
     void update() {
          if (orientation_backup)
               delete orientation_backup;
          if (origin_backup)
               delete origin_backup;
          orientation_backup = orientation;
          origin_backup = origin;
          
          Atom *origin_atom = (*res)[definitions::CA];
          this->origin = new Vector_3D(origin_atom->position);
          Atom *atom2 = origin_atom->get_neighbour<-1>(definitions::BACKBONE);
          if (!atom2)
               atom2 = origin_atom->get_neighbour<+2>(definitions::BACKBONE);
          Atom *atom3 = origin_atom->get_neighbour<+1>(definitions::BACKBONE);
          if (!atom3)
               atom3 = origin_atom->get_neighbour<-2>(definitions::BACKBONE);
          Vector_3D e1 = (*origin - atom2->position).normalize();
          Vector_3D e2 = cross_product(atom3->position - atom2->position, e1).normalize();
          Vector_3D e3 = cross_product(e1, e2);

          // Use first atom residue to define basis
          // this->origin = new Vector_3D((*res)[N]->position);
          // Vector_3D e1 = ((*res)[CA]->position - *origin).normalize();
          // Vector_3D e2 = cross_product((*res)[C]->position - *origin, e1).normalize();
          // Vector_3D e3 = cross_product(e1, e2);

          bool as_columns=true;
          this->orientation = new Matrix_3D(e3, e2, e1, as_columns);

     }

     //! Undo last update
     void undo() {
          if (orientation)
               delete orientation;
          if (origin)
               delete origin;
          orientation = orientation_backup;
          origin = origin_backup;
          
          orientation_backup = NULL;
          origin_backup = NULL;

          assert(orientation != NULL);
     }

     //! Transform coordinate from global frame to local frame
     //!
     //! \param v 3D-coordinate
     //! \return Transformed 3D-coordinate
     Vector_3D transform(Vector_3D v) {
          return transpose(*orientation)*(v-*origin);
     }

     //! Transform vector of coordinates from global frame to local frame
     //!
     //! \param vec vector of 3D-coordinates
     //! \return Vector of transformed 3D-coordinates
     std::vector<Vector_3D> transform(std::vector<Vector_3D> vec) {
          std::vector<Vector_3D> ret_vec;
          for (int i=0; i<vec.size(); i++) {
               ret_vec.push_back(this->transform(vec[i]));
          }
          return ret_vec;
     }

     //! Transform vector of atom-positions from global frame to local frame
     //!
     //! \param vec vector of atoms
     //! \return Vector of transformed 3D-coordinates
     std::vector<Vector_3D> transform(std::vector<Atom *> vec) {
          std::vector<Vector_3D> ret_vec;
          for (unsigned int i=0; i<vec.size(); i++) {
               ret_vec.push_back(this->transform(vec[i]->position));
          }
          return ret_vec;
     }

     //! Transform coordinate from local frame to global frame
     //!
     //! \param v 3D-coordinate
     //! \return Transformed 3D-coordinate
     Vector_3D transform_back(Vector_3D v) {
          return *origin+*orientation*(v);
     }

     //! Transform vector of coordinates from local frame to global frame
     //!
     //! \param vec vector of 3D-coordinates
     //! \return Vector of transformed 3D-coordinates
     std::vector<Vector_3D> transform_back(std::vector<Vector_3D> vec) {
          std::vector<Vector_3D> retVec;
          for (int i=0; i<vec.size(); i++) {
               retVec.push_back(this->transform_back(vec[i]));
          }
          return retVec;
     }

     //! Transform vector of atom-positions from local frame to global frame
     //!
     //! \param vec vector of atoms
     //! \return Vector of transformed 3D-coordinates
     std::vector<Vector_3D> transform_back(std::vector<Atom *> vec) {
          std::vector<Vector_3D> ret_vec;
          for (unsigned int i=0; i<vec.size(); i++) {
               ret_vec.push_back(this->transform_back(vec[i]->position));
          }
          return ret_vec;
     }     
};



//! Bounding volume: Rectangular-Swept-Sphere.
//! Each node in the chaintree contains a bounding volume (e.g. a RSS)
//! containing all volumes of child nodes - or atom positions if the children
//! are leaf nodes.
//!
//! An RSS also maintains its own frame of reference, fascilitating
//! the calculation of distances between two volumes
//! Note that the vertices are NOT specified in this reference
//! frame, but in the reference frame in which the positions
//! are given upon construction
//!
//! For efficiency reasons, we also maintain the center and radius of a
//! circle in the plane of the RSS, used for a preliminary distance cutoff check.
class RSS {
public:

     //! Radius of swept sphere
     double radius;

     //! Sidelengths of volume     
     double side_lengths[2];

     //! Radius of the circle enclosing the rectangle
     double circle_radius;

     //! Center of the circle enclosing the rectangle;
     double circle_center[2];

     //! basis for local reference frame
     Matrix_3D basis;

     //! basis for local reference frame - transpose version (cached for efficiency)
     Matrix_3D basis_transpose;

     //! Position of origin of local frame in global reference system
     Vector_3D reference_pos;

     //! List of 8 vertices defining the box:
     std::vector<Vector_3D> vertices;

     //! A slightly tighter box can be used when considering
     //! the bounding volume as a simple box rather than a swept sphere
     //! This cannot be used for distance measurements, but is useful
     //! when parents base their bounding volumes of the vertices of their children
     std::vector<Vector_3D> vertices_tight;
     

     //! Specifies whether tight vertices should be maintained
     //! This allows parent boxes to be slightly tighter than otherwise,
     //! (giving more efficient checking), but there is a slight
     //! price to pay in updating time
     const static bool use_tight_vertices = false;
     
     //! Constructor (from positions).
     //!
     //! \param positions Vector of 3D coordinates
     RSS(std::vector<Vector_3D> *positions)
          : vertices(8) {
          compute_bounding_volume(*positions);          
     }

     //! Constructor (from 2 child volumes).
     //!
     //! \param child1 Child1 RSS volume
     //! \param child2 Child2 RSS volume
     //! \param translation Translation vector relating local reference frame of child1 and child2
     //! \param rotation Rotation matrix relating local reference frame of child1 and child2
     RSS(RSS *child1, RSS *child2, Vector_3D *translation, Matrix_3D *rotation)
          : vertices(8) {

          std::vector<Vector_3D> combined_vertices(16);
          get_vertices(combined_vertices, child1, child2, *translation, *rotation);
          compute_bounding_volume(combined_vertices);
     }

     //! Constructor (from 4 child volumes).
     //! This is an attempt at making the bounding volumes tighter
     //! by constructing boxes based on the vertices of grand-children
     //! rather than the direct children.
     //! NOTE: there is a payoff here, since every update requires
     //! 32 positions rather than just 16.
     //!
     //! \param child1 Grandchild1 RSS volume
     //! \param child2 Grandchild2 RSS volume
     //! \param child3 Grandchild3 RSS volume
     //! \param child4 Grandchild4 RSS volume
     //! \param translation1 Translation vector relating local reference frame of child1 and child2
     //! \param rotation1 Rotation matrix relating local reference frame of child1 and child2
     //! \param translation2 Translation vector relating local reference frame of child2 and child3
     //! \param rotation2 Rotation matrix relating local reference frame of child2 and child3
     //! \param translation3 Translation vector relating local reference frame of child3 and child4
     //! \param rotation3 Rotation matrix relating local reference frame of child3 and child4
     RSS(RSS *child1, RSS *child2, RSS *child3, RSS *child4,
         Vector_3D *translation1, Matrix_3D *rotation1,
         Vector_3D *translation2, Matrix_3D *rotation2,
         Vector_3D *translation3, Matrix_3D *rotation3)
          : vertices(8) {
          
          std::vector<Vector_3D> combined_vertices(16);
          get_vertices(combined_vertices, child1, child2, child3, child4,
                      *translation1, *rotation1,
                      *translation2, *rotation2,
                      *translation3, *rotation3);
          compute_bounding_volume(combined_vertices);
          
     }

     //! Destructor.
     ~RSS(){}

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(RSS &other) {
          bool equal = true;

          if (this->vertices.size() != other.vertices.size()) {
               equal = false;
          } else {
               for (unsigned int i=0; i<this->vertices.size(); i++) {
                    if (sum(this->vertices[i] - other.vertices[i]) > 0.1) {
                         equal = false;
                    }
               }
               for (unsigned int i=0; i<this->vertices_tight.size(); i++) {
                    if (sum(this->vertices_tight[i] - other.vertices_tight[i]) > 0.1) {
                         equal = false;
                    }
               }
          }

          return equal;
     }

     //! Finds the vertices of the current RSS, given the vertices of two
     //! child volumes.
     //!
     //! \param combined_vertices destination vector for vertex coordinates
     //! \param child1 Child1 RSS volume
     //! \param child2 Child2 RSS volume
     //! \param translation Translation vector relating local reference frame of child1 and child2
     //! \param rotation Rotation matrix relating local reference frame of child1 and child2
     void get_vertices(std::vector<Vector_3D> &combined_vertices, RSS *child1, RSS *child2,
                       Vector_3D &translation, Matrix_3D &rotation) {

          if (use_tight_vertices) {

               // A parent has the same reference system as its left child
               unsigned int i;
               for (i=0; i<child1->vertices_tight.size(); i++) {
                    combined_vertices[i] = child1->vertices_tight[i];
               }
               
               // The vertex positions of the right child are transformed to current frame
               for (unsigned int j=0; j<child2->vertices_tight.size(); j++) {
                    combined_vertices[i+j] = rotation*child2->vertices_tight[i]+translation;
               }

          } else {

               // combinedVertices.resize(child1->vertices.size() + child2->vertices.size());

               // A parent has the same reference system as its left child
               unsigned int i;
               for (i=0; i<child1->vertices.size(); i++) {
                    combined_vertices[i] = child1->vertices[i];
               }
               
               // The vertex positions of the right child are transformed to current frame
               for (unsigned int j=0; j<child2->vertices.size(); j++) {
                    combined_vertices[i+j] = rotation*child2->vertices[j]+translation;
               }

          }
     }

     //! Finds the vertices of the current RSS, given the vertices of four
     //! grandchild volumes
     //!
     //! \param combined_vertices destination vector for vertex coordinates
     //! \param child1 Grandchild1 RSS volume
     //! \param child2 Grandchild2 RSS volume
     //! \param child3 Grandchild3 RSS volume
     //! \param child4 Grandchild4 RSS volume
     //! \param translation1 Translation vector relating local reference frame of child1 and child2
     //! \param rotation1 Rotation matrix relating local reference frame of child1 and child2
     //! \param translation2 Translation vector relating local reference frame of child2 and child3
     //! \param rotation2 Rotation matrix relating local reference frame of child2 and child3
     //! \param translation3 Translation vector relating local reference frame of child3 and child4
     //! \param rotation3 Rotation matrix relating local reference frame of child3 and child4
     void get_vertices(std::vector<Vector_3D> &combined_vertices,
                       RSS *child1, RSS *child2, RSS *child3, RSS *child4,
                       Vector_3D &translation1, Matrix_3D &rotation1,
                       Vector_3D &translation2, Matrix_3D &rotation2,
                       Vector_3D &translation3, Matrix_3D &rotation3) {

          if (use_tight_vertices) {
                    
               // A parent has the same reference system as its left child
               unsigned int i;
               for (i=0; i<child1->vertices_tight.size(); i++) {
                    combined_vertices[i] = child1->vertices_tight[i];
               }

               // The vertex positions of child2 are transformed to current frame
               unsigned int j;
               for (j=0; j<child2->vertices_tight.size(); j++) {
                    combined_vertices[i+j] = rotation1*child2->vertices_tight[j]+translation1;
               }

               if (child3) {

                    Matrix_3D rotation12 = rotation1*rotation2;
                    Vector_3D translation12 = rotation1*translation2 + translation1;
               
                    // The vertex positions of child3 are transformed to current frame
                    unsigned int k;
                    for (k=0; k<child3->vertices_tight.size(); k++) {
                         combined_vertices[k+j] = rotation12*child3->vertices_tight[k]+translation12;
                    }

                    if (child4) {
                         Matrix_3D rotation13 = rotation12*rotation3;
                         Vector_3D translation13 = rotation12*translation3 + translation12;
                    
                         // The vertex positions of child4 are transformed to current frame
                         unsigned int l;
                         for (l=0; l<child4->vertices_tight.size(); l++) {
                              combined_vertices[l+k] = rotation13*child4->vertices_tight[l]+translation13;
                         }
                    } else {
                         combined_vertices.resize(12);
                    }
               } else {
                    combined_vertices.resize(8);
               }

          } else {

               // A parent has the same reference system as its left child
               unsigned int i;
               for (i=0; i<child1->vertices.size(); i++) {
                    combined_vertices.push_back(child1->vertices[i]);
               }

               // The vertex positions of child2 are transformed to current frame
               unsigned int j;
               for (j=0; j<child2->vertices.size(); j++) {
                    combined_vertices[i+j] = rotation1*child2->vertices[j]+translation1;
               }

               if (child3) {

                    Matrix_3D rotation12 = rotation1*rotation2;
                    Vector_3D translation12 = rotation1*translation2 + translation1;
               
                    // The vertex positions of child3 are transformed to current frame
                    unsigned int k;
                    for (k=0; k<child3->vertices.size(); k++) {
                         combined_vertices[k+j] = rotation12*child3->vertices[k]+translation12;
                    }

                    if (child4) {
                         Matrix_3D rotation13 = rotation12*rotation3;
                         Vector_3D translation13 = rotation12*translation3 + translation12;
                    
                         // The vertex positions of child4 are transformed to current frame
                         unsigned int l;
                         for (l=0; l<child4->vertices.size(); l++) {
                              combined_vertices[l+k] = rotation13*child4->vertices[l]+translation13;
                         }
                    } else {
                         combined_vertices.resize(12);
                    }
               } else {
                    combined_vertices.resize(8);
               }
          }
     }


     //! Calculate the bounding volume from a vector of positions
     //!
     //! \param positions Vector of 3D coordinates
     void compute_bounding_volume(std::vector<Vector_3D> &positions) {

          // Special case for one position
          if (positions.size() == 1) {
               radius = 0;
               reference_pos = positions[0];
               this->side_lengths[0] = 0;
               this->side_lengths[1] = 0;
               this->circle_radius = 0;
               this->circle_center[0] = 0;
               this->circle_center[1] = 0;
               basis = identity_matrix();
               basis_transpose = transpose(basis);
               if (use_tight_vertices)
                    compute_vertices_tight();
               compute_vertices();
               return;

          // Special case for two positions
          } else if (positions.size() ==2) {
               radius = 0;
               reference_pos = positions[0];
               Vector_3D eigvec1 = (positions[1] - positions[0]).normalize();
               Vector_3D eigvec2 = (Vector_3D(1,0,0) - (Vector_3D(1,0,0)*eigvec1)*eigvec1).normalize();
               Vector_3D eigvec3 = eigvec1%eigvec2;
               this->side_lengths[0] = (positions[1] - positions[0]).norm();
               this->side_lengths[1] = 0;
               circle_radius = this->side_lengths[0]/2.0;
               this->circle_center[0] = circle_radius;
               this->circle_center[1] = 0;
               bool asColumns = true;
               basis = Matrix_3D(eigvec1, eigvec2, eigvec3, asColumns);
               basis_transpose = transpose(basis);
               if (use_tight_vertices)
                    compute_vertices_tight();
               compute_vertices();
               return;
          }

          Vector_3D center_of_mass(0,0,0);
          for (unsigned int i=0; i<positions.size(); ++i) {
               center_of_mass += positions[i];
          }
          center_of_mass /= positions.size();

          Matrix_3D cov(null_matrix());
          for (unsigned int i=0; i<positions.size(); ++i) {
               Vector_3D v = (positions[i]-center_of_mass);

               // Outer product: cov += v^v;
               cov(0,0) += v[0]*v[0];
               cov(0,1) += v[0]*v[1];
               cov(0,2) += v[0]*v[2];
               cov(1,1) += v[1]*v[1];
               cov(1,2) += v[1]*v[2];
               cov(2,2) += v[2]*v[2];
          }          
          cov(1,0) = cov(0,1);
          cov(2,1) = cov(1,2);
          cov(2,0) = cov(0,2);

          // // Get eigenvectors
          // Matrix eigenVectors(3,3);
          // Vector_nD *eigenValues = NULL;
          // Matrix(cov).eigen_sym(eigenValues, &eigenVectors);

          // // Eigenvectors are sorted in ascending order
          // Vector_3D eigvec1(eigenVectors.get_column(2).get_array());
          // Vector_3D eigvec2(eigenVectors.get_column(1).get_array());
          // Vector_3D eigvec3 = eigvec1%eigvec2; // cross-product

          // // Create local reference system
          bool as_columns = true;
          // basis = Matrix_3D(eigvec1, eigvec2, eigvec3, asColumns);

          EigenSystem3x3 eigen_system(cov);
          basis = Matrix_3D(eigen_system.eigen_vectors[2], 
                            eigen_system.eigen_vectors[1], 
                            eigen_system.eigen_vectors[0],
                            as_columns);
          basis_transpose = transpose(basis);

          // Transform positions to local reference frame
          std::vector<Vector_3D> local_positions(positions.size());
          for (unsigned int i=0; i<positions.size(); i++) {
               local_positions[i] = basis_transpose*positions[i];
          }

          // Minimum and maximum values in x,y,z dimensions
          double min_x, max_x;
          double min_y, max_y;
          double min_z, max_z;

          // Minimum and maximum indices in x,y,z dimensions
          int min_index_x, max_index_x;
          int min_index_y, max_index_y;
          int min_index_z, max_index_z;

          // Find boundaries for z
          find_bounding_points(local_positions, 2, &min_index_z, &max_index_z, &min_z, &max_z);
          
          // maxz-minz is used to define the radius
          this->radius = 0.5*(max_z - min_z);

          // Set radius to zero if points are coplanar
          if (radius < 1e-10) {
               radius = 0.0;
          }

          // Middle of z-range
          double center_z = 0.5*(max_z + min_z);

          // Find boundaries for x
          find_bounding_points(local_positions, 0, &min_index_x, &max_index_x, &min_x, &max_x);

          // Find boundaries for y
          find_bounding_points(local_positions, 1, &min_index_y, &max_index_y, &min_y, &max_y);

          Vector_3D c;

          ////// TIGHT BOUNDING BOX //////

          if (use_tight_vertices) {
               
               // reference point of tight bounding box
               c = Vector_3D(min_x, min_y, center_z);

               // Transform back to original coordinates
               this->reference_pos = basis*c;
          
               // box side lengths
               this->side_lengths[0] = std::max(max_x-min_x, 0.0);
               this->side_lengths[1] = std::max(max_y-min_y, 0.0);

               // Compute tight bounding box
               compute_vertices_tight();
          }

          ////// RECTANGULAR SWEPT SPHERE BOUNDING BOX ////
          
          // Adjust boundaries (given round edges)
          adjust_boundaries(local_positions, center_z, radius, 0, min_index_x, max_index_x, &min_x, &max_x);
          adjust_boundaries(local_positions, center_z, radius, 1, min_index_y, max_index_y, &min_y, &max_y);

          // Adjust corners of bounding box
          adjust_corners(local_positions, center_z, radius, &min_x, &max_x, &min_y, &max_y);

          // reference point of bounding box
          c = Vector_3D(min_x, min_y, center_z);

          // Transform back to original coordinates
          this->reference_pos = basis*c;

          // box side lengths
          this->side_lengths[0] = std::max(max_x-min_x, 0.0);
          this->side_lengths[1] = std::max(max_y-min_y, 0.0);

          this->circle_center[0] = this->side_lengths[0]/2.0;
          this->circle_center[1] = this->side_lengths[1]/2.0;
          this->circle_radius = sqrt(circle_center[0]*circle_center[0]+circle_center[1]*circle_center[1]);

          // Calculate vertices of bounding volume
          compute_vertices();
     }


     //! Finds positions at the boundaries of the volume.
     //!
     //! \param positions Vector of 3D coordinates
     //! \param dimension x,y or z-dimension (x=0,y=1,z=2)
     //! \param min_index Index of minimum element in positions vector
     //! \param max_index Index of maximum element in positions vector
     //! \param min_val Minimum value
     //! \param max_val Maximum value
     void find_bounding_points(std::vector<Vector_3D> &positions, int dimension,
                               int *min_index, int *max_index, double *min_val, double *max_val) {

          *min_index=0, *max_index=0;
          for (unsigned int i=0; i<positions.size(); i++) {
               if (positions[i][dimension] < positions[*min_index][dimension]) {
                    *min_index = i;
               } else if (positions[i][dimension] > positions[*max_index][dimension]) {
                    *max_index = i;
               }
          }
          *min_val = positions[*min_index][dimension];
          *max_val = positions[*max_index][dimension];
     }

     //! Adjust boundaries of RSS. The bounding box can be reduced slightly since it is swept by a sphere.
     //!
     //! \param positions Vector of 3D coordinates
     //! \param center_z Center in z-direction
     //! \param radius of swept sphere
     //! \param dimension x,y or z-dimension (x=0,y=1,z=2)
     //! \param min_index Index of minimum element in positions vector
     //! \param max_index Index of maximum element in positions vector
     //! \param min_val Minimum value
     //! \param max_val Maximum value
     void adjust_boundaries(std::vector<Vector_3D> &positions, double center_z, double radius, int dimension,
                             int min_index, int max_index, double *min_val, double *max_val) {

          // Calculate z-distance from rectangle center
          double delta_z_min_index = positions[min_index][2] - center_z;
          double delta_z_max_Index = positions[max_index][2] - center_z;

          // Adjust boundaries:
          //    o
          // dz |\ rad
          //    | \          .
          //    o--o---------
          //    minx -> minx+sqrt(rad^2-dz^2)
          double r2 = radius*radius;
          *min_val += sqrt(std::max(r2 - delta_z_min_index*delta_z_min_index,0.0));
          *max_val -= sqrt(std::max(r2 - delta_z_max_Index*delta_z_max_Index,0.0));

          // Some of the other positions may now have fallen outside the range
          // Check for this and adjust boundaries if necessary
          for (unsigned int i=0; i<positions.size(); i++) {

               double dz = positions[i][2] - center_z;
               
               // Minimum
               if (positions[i][dimension] < *min_val) {
                    double val = positions[i][dimension] + sqrt(std::max(r2 - dz*dz,0.0));
                    if (val < *min_val) {
                         *min_val = val;
                    }
               }

               // Maximum
               if (positions[i][dimension] > *max_val) {
                    double val = positions[i][dimension] - sqrt(std::max(r2 - dz*dz,0.0));
                    if (val > *max_val) {
                         *max_val = val;
                    }
               }
          }
     }


     //! Adjust corners of RSS. The corners of the bounding box can be reduced slightly since it is swept by a sphere
     //!
     //! \param positions Vector of 3D coordinates
     //! \param center_z Center in z-direction
     //! \param radius Radius of swept sphere
     //! \param min_x Minimum value in x direction
     //! \param max_x Maximum value in x direction
     //! \param min_y Minimum value in y direction
     //! \param max_y Maximum value in y direction
     void adjust_corners(std::vector<Vector_3D> &positions, double center_z, double radius, 
                         double *min_x, double *max_x, double *min_y, double *max_y) {

          for (unsigned int i=0; i<positions.size(); i++) {

               Vector_3D unit_vector;
               double *boundary_x=NULL, *boundary_y=NULL;
               if (positions[i][0] > *max_x) {
                    if (positions[i][1] > *max_y) {
                         unit_vector = Vector_3D(Math<double>::inv_sqrt_2, Math<double>::inv_sqrt_2, 0);
                         boundary_x = max_x;
                         boundary_y = max_y;
                    } else if (positions[i][1] < *min_y) {
                         unit_vector = Vector_3D(Math<double>::inv_sqrt_2, -Math<double>::inv_sqrt_2, 0); 
                         boundary_x = max_x;
                         boundary_y = min_y;
                    }
               } else if (positions[i][0] < *min_x) {
                    if (positions[i][1] > *max_y) {
                         unit_vector = Vector_3D(-Math<double>::inv_sqrt_2, Math<double>::inv_sqrt_2, 0); 
                         boundary_x = min_x;
                         boundary_y = max_y;
                    } else if (positions[i][1] < *min_y) {
                         unit_vector = Vector_3D(-Math<double>::inv_sqrt_2, -Math<double>::inv_sqrt_2, 0); 
                         boundary_x = min_x;
                         boundary_y = min_y;
                    }
               }

               if (boundary_x) {

                    Vector_3D dv(positions[i][0] - *boundary_x,
                                 positions[i][1] - *boundary_y,
                                 positions[i][2] - center_z);

                    // Project in 45 degree direction
                    double u_len = (dv*unit_vector);
                    Vector_3D u = u_len*unit_vector;

                    // Calculate offset vector from 45 degree direction
                    Vector_3D t = dv-u;

                    // Reduce as much as is supported by radius (see adjust_boundaries)
                    u_len -= sqrt(std::max(radius*radius - (t[0]*t[0]+t[1]*t[1]+t[2]*t[2]), 0.0));

                    // if u_len is larger than zero, the position was not covered by the RSS
                    // and the boundaries should be adjusted
                    if (u_len > 0) {
                         *boundary_x += u_len*unit_vector[0];
                         *boundary_y += u_len*unit_vector[1];
                    }
               }
          }
          
     }

     //! Calculate the vertices of the bounding volume
     void compute_vertices() {

          // First vertex is the reference point offset by (radius,radius,radius)
          // vertices.push_back(reference_pos - basis*Vector_3D(radius, radius, radius));
          vertices[0] = reference_pos - basis*Vector_3D(radius, radius, radius);

          // Second vertex placed along z-axis
          // vertices.push_back(vertices[0] + basis.col_vector(2)*(2*radius));
          vertices[1] = vertices[0] + basis.col_vector(2)*(2*radius);

          // Third vertex placed along x-axis
          vertices[2] = vertices[0] + basis.col_vector(0)*(2*radius+side_lengths[0]);

          // Fourth vertex placed along z-axis from vertex 2
          // vertices.push_back(vertices[2] + basis.col_vector(2)*(2*radius));
          vertices[3] = vertices[2] + basis.col_vector(2)*(2*radius);

          // Fifth vertex placed along y-axis from vertex 0
          // vertices.push_back(vertices[0] + basis.col_vector(1)*(2*radius+side_lengths[1]));
          vertices[4] = vertices[0] + basis.col_vector(1)*(2*radius+side_lengths[1]);

          // Sixth vertex placed along z-axis from vertex 4
          // vertices.push_back(vertices[4] + basis.col_vector(2)*(2*radius));
          vertices[5] = vertices[4] + basis.col_vector(2)*(2*radius);

          // Seventh vertex placed along x-axis from vertex 4
          // vertices.push_back(vertices[4] + basis.col_vector(0)*(2*radius+side_lengths[0]));
          vertices[6] = vertices[4] + basis.col_vector(0)*(2*radius+side_lengths[0]);

          // Eighth vertex placed along z-axis from vertex 6
          // vertices.push_back(vertices[6] + basis.col_vector(2)*(2*radius));
          vertices[7] = vertices[6] + basis.col_vector(2)*(2*radius);
     }

     
     //! Calculate the vertices of the bounding volume - tighter version.
     //! A slightly tighter box without the swept sphere
     //! used only to update parent values in tree.
     void compute_vertices_tight() {

          vertices_tight.clear();
          
          // First vertex is the reference point offset by (radius,radius,radius)
          vertices_tight.push_back(reference_pos - basis*Vector_3D(0.0, 0.0, radius));

          // Second vertex placed along z-axis
          vertices_tight.push_back(vertices_tight[0] + basis.col_vector(2)*(2*radius));

          // Third vertex placed along x-axis
          vertices_tight.push_back(vertices_tight[0] + basis.col_vector(0)*(side_lengths[0]));

          // Fourth vertex placed along z-axis from vertex 2
          vertices_tight.push_back(vertices_tight[2] + basis.col_vector(2)*(2*radius));

          // Fifth vertex placed along y-axis from vertex 0
          vertices_tight.push_back(vertices_tight[0] + basis.col_vector(1)*(side_lengths[1]));

          // Sixth vertex placed along z-axis from vertex 4
          vertices_tight.push_back(vertices_tight[4] + basis.col_vector(2)*(2*radius));

          // Seventh vertex placed along x-axis from vertex 4
          vertices_tight.push_back(vertices_tight[4] + basis.col_vector(0)*(side_lengths[0]));

          // Eighth vertex placed along z-axis from vertex 6
          vertices_tight.push_back(vertices_tight[6] + basis.col_vector(2)*(2*radius));
     }


     //! Calculate distance between two bounding volumes
     //! translation and rotation denote how the frame in which the positions of rss2 are defined
     //! should be transformed into the frame of the positions of the current rss
     //! (note that this is not the transformation between the basis frames of the bounding volumes
     //! themselves, but of the frame used to specify the vertices (typically the local frames
     //! of the chaintree nodes). 
     //!
     //! \param rss2 object to which distance is measured
     //! \param translation Translation vector specifying how frame of rss2 
     //!                    should be translated to current local frame
     //! \param rotation Rotation matrix specifying how frame of rss2 should be rotated 
     //!                 to current local frame
     //! \param cutoff_distance distance cutoff
     //! \return Whether distance between current rss and rss2 is within cutoff distance.
     bool distance_within_limit(RSS &rss2, const Vector_3D &translation, const Matrix_3D &rotation, double cutoff_distance) {

          // Boxes are placed so first box is at the origin with its sides along the axes 
          // Place box with corner at (0,0,0) by subtracting referencePoint from translation vector
          Vector_3D translation_new = basis_transpose*((rotation*rss2.reference_pos) +
                                                       (translation-reference_pos));

          // Calculate rotation of second box relative to first box
          Matrix_3D rotation_new = basis_transpose*(rotation*rss2.basis);

          // The distance between two rectangular swept sphere volumes
          // is the distance between the inner rectangles substracted
          // by the distance of the radii of the swept spheres.
          // We therefore add the radii to the cutoff
          // double dist = rect_distance_within_limit(translation_new, rotation_new, this->side_lengths, rss2.side_lengths, cutoff_distance) - (this->radius+rss2.radius);
          bool within_limit = rect_distance_within_limit(translation_new, rotation_new, rss2, cutoff_distance + (this->radius+rss2.radius));

          return within_limit;
     }

     //! Calculate distance between two bounding volumes
     //! translation and rotation denote how the frame in which the positions of rss2 are defined
     //! should be transformed into the frame of the positions of the current rss
     //! (note that this is not the transformation between the basis frames of the bounding volumes
     //! themselves, but of the frame used to specify the vertices (typically the local frames
     //! of the chaintree nodes). 
     //!
     //! \param rss2 object to which distance is measured
     //! \param translation Translation vector specifying how frame of rss2 
     //!                    should be translated to current local frame
     //! \param rotation Rotation matrix specifying how frame of rss2 should be rotated 
     //!                 to current local frame
     //! \return Distance between current rss and rss2
     double compute_distance(RSS &rss2, const Vector_3D &translation, const Matrix_3D &rotation) {

          // Boxes are placed so first box is at the origin with its sides along the axes 
          // Place box with corner at (0,0,0) by subtracting referencePoint from translation vector
          Vector_3D translation_new = basis_transpose*((rotation*rss2.reference_pos) +
                                                       (translation-reference_pos));

          // Calculate rotation of second box relative to first box
          Matrix_3D rotation_new = basis_transpose*(rotation*rss2.basis);

          // The distance between two rectangular swept sphere volumes
          // is the distance between the inner rectangles substracted
          // by the distance of the radii of the swept spheres
          double dist = rect_distance(translation_new, rotation_new, this->side_lengths, rss2.side_lengths) - (this->radius+rss2.radius);

          return dist;
     }


     //! Calculate distance between bounding volume and a single position
     //!
     //! \param pos2 3D-coordinate to which distance is measured
     //! \param translation Translation vector specifying how pos2
     //!                    should be translated to current local frame
     //! \param rotation Rotation matrix specifying how pos2 should be rotated 
     //!                 to current local frame
     //! \return Distance between current rss and pos2
     double compute_distance(Vector_3D &pos2, const Vector_3D &translation, const Matrix_3D &rotation) {

          // Boxes are placed so first box is at the origin with its sides along the axes 
          
          // Place box with corner at (0,0,0) by subtracting referencePoint from translation vector
          Vector_3D translation_new = basis_transpose*((rotation*pos2) + (translation-reference_pos));

          // Calculate rotation of second box relative to first box
          // Matrix_3D rotation_new = transpose(basis)*(rotation*rss2.basis);
          Matrix_3D rotation_new = basis_transpose*rotation;

          // The distance between a rectangular swept sphere volume and
          // a position is the distance between the inner rectangle of the RSS and the
          // position substracted by the distance of the radii of the swept sphere
          double dist = rect_distance(translation_new, rotation_new, this->side_lengths) - this->radius;

          return dist;
     }
     
     //! Calculate distance to another bounding volume
     //! in same reference frame
     //!
     //! \param rss2 object to which distance is measured     
     //! \return Distance between current rss and rss2
     double compute_distance(RSS &rss2) {
          Vector_3D translation = null_vector();
          Matrix_3D rotation = identity_matrix();
          return compute_distance(rss2, translation, rotation);
     }

     //! Calculate distance to a position
     //! in same reference frame
     //!
     //! \param pos2 3D-coordinate to which distance is measured
     //! \return Distance between current rss and pos2
     double compute_distance(Vector_3D &pos2) {
          Vector_3D translation = null_vector();
          Matrix_3D rotation = identity_matrix();
          return compute_distance(pos2, translation, rotation);
     }     

     
     //! Determine whether the distance between two rectangles in 3D-space is within a provided cutoff.
     //!
     //! \param translation translation vector relating reference frame A (this) and B
     //! \param rotation rotation matrix relating reference frame A(this) and B
     //! \param rss2 object to which distance is measured     
     //! \param cutoff_distance distance cutoff
     //! \return Whether distance between current rss and rss2 is within cutoff distance.
     inline bool rect_distance_within_limit(Vector_3D &translation, Matrix_3D &rotation,
                                            RSS &rss2, double cutoff_distance) {

          double *side_lengths_A = this->side_lengths;
          double *side_lengths_B = rss2.side_lengths;

          double cutoff_distance_squared = cutoff_distance*cutoff_distance;

          // Compute translation vector in B's frame of reference
          Vector_3D translation_BA = transpose(rotation)*translation;

          // Column vectors in rotation matrix
          Vector_3D rotation0 = rotation.col_vector(0);
          Vector_3D rotation1 = rotation.col_vector(1);
          Vector_3D rotation2 = rotation.col_vector(2);
     
          // Calculate dot products between unit vectors in the two reference frames
          double AX_dot_BX = rotation(0,0);
          double AX_dot_BY = rotation(0,1);
          double AY_dot_BX = rotation(1,0);
          double AY_dot_BY = rotation(1,1);

          // Projection of A's vertices into in B's reference frame: x-coordinates
          double A_LL_X = -translation_BA[0];
          double A_LU_X = A_LL_X + side_lengths_A[1]*AY_dot_BX;
          double A_UL_X = A_LL_X + side_lengths_A[0]*AX_dot_BX;
          double A_UU_X = A_LU_X + side_lengths_A[0]*AX_dot_BX;

          // Projection of B's vertices into in A's reference frame: x-coordinates
          double B_LL_X = translation[0];
          double B_LU_X = B_LL_X + side_lengths_B[1]*AX_dot_BY;
          double B_UL_X = B_LL_X + side_lengths_B[0]*AX_dot_BX;
          double B_UU_X = B_LU_X + side_lengths_B[0]*AX_dot_BX;

          // Initial sphere-sphere check
          double AZ_dot_BX = rotation(2,0);
          double AZ_dot_BY = rotation(2,1);
          double b_center_x = translation[0] + rss2.circle_center[1]*AX_dot_BY + rss2.circle_center[0]*AX_dot_BX;
          double b_center_y = translation[1] + rss2.circle_center[1]*AY_dot_BY + rss2.circle_center[0]*AY_dot_BX;
          double b_center_z = translation[2] + rss2.circle_center[1]*AZ_dot_BY + rss2.circle_center[0]*AZ_dot_BX;
          double center_center_x = b_center_x - circle_center[0];
          double center_center_y = b_center_y - circle_center[1];
          double cutoff_circle = cutoff_distance+circle_radius+rss2.circle_radius;
          if ((center_center_x*center_center_x + center_center_y*center_center_y + b_center_z*b_center_z) > cutoff_circle*cutoff_circle) {
               return false;
          }


          Vector_3D shortest_vector;

          // Edges in rectangle:
          //   ^
          // LU|____ UU
          //   |    |
          //   |____|____>
          // LL      UL
     
          // 1. check edge A: UL-UU vs B: UL-UU
          if (edge_pair_check(shortest_vector,
                              A_UL_X, A_UU_X, B_UL_X, B_UU_X,
                              side_lengths_A[0], +1, side_lengths_A[1],
                              side_lengths_B[0], +1, side_lengths_B[1],
                              0,
                              AY_dot_BY, AY_dot_BX, AX_dot_BY, AX_dot_BX,
                              translation[1], translation[0],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 1\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
         
          // 2. check edge A: UL-UU vs B: LL-LU
          if (edge_pair_check(shortest_vector,
                              A_UL_X, A_UU_X, B_LL_X, B_LU_X,
                              side_lengths_A[0], +1, side_lengths_A[1],
                              0,               -1, side_lengths_B[1],
                              0,
                              AY_dot_BY, AY_dot_BX, AX_dot_BY, AX_dot_BX,
                              translation[1], translation[0],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 2\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
         
          // 3. check edge A: LL-LU vs B: UL-UU
          if (edge_pair_check(shortest_vector,
                              A_LL_X, A_LU_X, B_UL_X, B_UU_X,
                              0,               -1, side_lengths_A[1],
                              side_lengths_B[0], +1, side_lengths_B[1],
                              0,
                              AY_dot_BY, AY_dot_BX, AX_dot_BY, AX_dot_BX,
                              translation[1], translation[0],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 3\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }

          // 4. check edge A: LL-LU vs B: LL-LU
          if (edge_pair_check(shortest_vector,
                              A_LL_X, A_LU_X, B_LL_X, B_LU_X,
                              0,               -1, side_lengths_A[1],
                              0,               -1, side_lengths_B[1],
                              0,
                              AY_dot_BY, AY_dot_BX, AX_dot_BY, AX_dot_BX,
                              // AX_dot_BX, AX_dot_BY, AY_dot_BX, AY_dot_BY,
                              translation[1], translation[0],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 4   " << sqrt(shortestVector*shortestVector) << "\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
     
          // Projection of A's vertices into in B's reference frame: y-coordinates
          double A_LL_Y = -translation_BA[1];
          double A_LU_Y = A_LL_Y + side_lengths_A[1]*AY_dot_BY;
          double A_UL_Y = A_LL_Y + side_lengths_A[0]*AX_dot_BY;
          double A_UU_Y = A_LU_Y + side_lengths_A[0]*AX_dot_BY;


          // 5. check edge A: UL-UU vs B: LU-UU
          if (edge_pair_check(shortest_vector,
                              A_UL_Y, A_UU_Y, B_LU_X, B_UU_X,
                              side_lengths_A[0], +1, side_lengths_A[1],
                              side_lengths_B[1], +1, side_lengths_B[0],
                              0,
                              AY_dot_BX, AY_dot_BY, AX_dot_BX, AX_dot_BY,
                              translation[1], translation[0],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 5\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }


          // 6. check edge A: UL-UU vs B: LL-UL
          if (edge_pair_check(shortest_vector,
                              A_UL_Y, A_UU_Y, B_LL_X, B_UL_X,
                              side_lengths_A[0], +1, side_lengths_A[1],
                              0              , -1, side_lengths_B[0],
                              0,
                              AY_dot_BX, AY_dot_BY, AX_dot_BX, AX_dot_BY,
                              translation[1], translation[0],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 6\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }

          // 7. check edge A: LL-LU vs B: LU-UU
          if (edge_pair_check(shortest_vector,
                              A_LL_Y, A_LU_Y, B_LU_X, B_UU_X,
                              0,               -1, side_lengths_A[1],
                              side_lengths_B[1], +1, side_lengths_B[0],
                              0,
                              AY_dot_BX, AY_dot_BY, AX_dot_BX, AX_dot_BY,
                              translation[1], translation[0],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 7\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }

          // 8. check edge A: LL-LU vs B: LL-UL
          if (edge_pair_check(shortest_vector,
                              A_LL_Y, A_LU_Y, B_LL_X, B_UL_X,
                              0,               -1, side_lengths_A[1],
                              0,               -1, side_lengths_B[0],
                              0,
                              AY_dot_BX, AY_dot_BY, AX_dot_BX, AX_dot_BY,
                              translation[1], translation[0],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 8\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
     

          // Projection of B's vertices into in A's reference frame: y-coordinates
          double B_LL_Y = translation[1];
          double B_LU_Y = B_LL_Y + side_lengths_B[1]*AY_dot_BY;
          double B_UL_Y = B_LL_Y + side_lengths_B[0]*AY_dot_BX;
          double B_UU_Y = B_LU_Y + side_lengths_B[0]*AY_dot_BX;

     
          // 9. check edge A: LU-UU vs B: UL-UU
          if (edge_pair_check(shortest_vector,
                              A_LU_X, A_UU_X, B_UL_Y, B_UU_Y,
                              side_lengths_A[1], +1, side_lengths_A[0],
                              side_lengths_B[0], +1, side_lengths_B[1],
                              1,
                              AX_dot_BY, AX_dot_BX, AY_dot_BY, AY_dot_BX,
                              translation[0], translation[1],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 9   " << sqrt(shortestVector*shortestVector) << "\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }

          // 10. check edge A: LU-UU vs B: LL-LU
          if (edge_pair_check(shortest_vector,
                              A_LU_X, A_UU_X, B_LL_Y, B_LU_Y,
                              side_lengths_A[1], +1, side_lengths_A[0],
                              0,               -1, side_lengths_B[1],
                              1,
                              AX_dot_BY, AX_dot_BX, AY_dot_BY, AY_dot_BX,
                              translation[0], translation[1],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 10\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
     
          // 11. check edge A: LL-UL vs B: UL-UU
          if (edge_pair_check(shortest_vector,
                              A_LL_X, A_UL_X, B_UL_Y, B_UU_Y,
                              0,               -1, side_lengths_A[0],
                              side_lengths_B[0], +1, side_lengths_B[1],
                              1,
                              AX_dot_BY, AX_dot_BX, AY_dot_BY, AY_dot_BX,
                              translation[0], translation[1],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 11\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }

          // 12. check edge A: LL-UL vs B: LL-LU
          if (edge_pair_check(shortest_vector,
                              A_LL_X, A_UL_X, B_LL_Y, B_LU_Y,
                              0,               -1, side_lengths_A[0],
                              0,               -1, side_lengths_B[1],
                              1,
                              AX_dot_BY, AX_dot_BX, AY_dot_BY, AY_dot_BX,
                              translation[0], translation[1],
                              translation_BA[1], translation_BA[0],
                              translation, rotation1, rotation0)) {
               // std::cout << "rect_distance: 12\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }


          // 13. check edge A: LU-UU vs B: LU-UU
          if (edge_pair_check(shortest_vector,
                              A_LU_Y, A_UU_Y, B_LU_Y, B_UU_Y,
                              side_lengths_A[1], +1, side_lengths_A[0],
                              side_lengths_B[1], +1, side_lengths_B[0],
                              1,
                              AX_dot_BX, AX_dot_BY, AY_dot_BX, AY_dot_BY,
                              translation[0], translation[1],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 13\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }


          // 14. check edge A: LU-UU vs B: LL-UL
          if (edge_pair_check(shortest_vector,
                              A_LU_Y, A_UU_Y, B_LL_Y, B_UL_Y,
                              side_lengths_A[1], +1, side_lengths_A[0],
                              0,               -1, side_lengths_B[0],
                              1,
                              AX_dot_BX, AX_dot_BY, AY_dot_BX, AY_dot_BY,
                              translation[0], translation[1],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 14\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
     
     
          // 15. check edge A: LL-UL vs B: LU-UU
          if (edge_pair_check(shortest_vector,
                              A_LL_Y, A_UL_Y, B_LU_Y, B_UU_Y,
                              0,               -1, side_lengths_A[0],
                              side_lengths_B[1], +1, side_lengths_B[0],
                              1,
                              AX_dot_BX, AX_dot_BY, AY_dot_BX, AY_dot_BY,
                              translation[0], translation[1],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 15\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }
     

          // 16. check edge A: LL-UL vs B: LL-UL
          if (edge_pair_check(shortest_vector,
                              A_LL_Y, A_UL_Y, B_LL_Y, B_UL_Y,
                              0,               -1, side_lengths_A[0],
                              0,               -1, side_lengths_B[0],
                              1,
                              AX_dot_BX, AX_dot_BY, AY_dot_BX, AY_dot_BY,
                              translation[0], translation[1],
                              translation_BA[0], translation_BA[1],
                              translation, rotation0, rotation1)) {
               // std::cout << "rect_distance: 16\n";
               return (shortest_vector*shortest_vector < cutoff_distance_squared);
               // return sqrt(shortestVector*shortestVector);
          }

          // If none of the above tests were true, calculate maximum
          // distance along rectangle face normal vectors

          int sign1;
          if (translation[2] > 0) {     // B is "above" A
               sign1 = +1;
          } else {                              // B is "below" A
               sign1 = -1;
          }
          double separation1 = translation[2]*sign1;
          
          // Test if x-unit vector in B points towards A
          // if so, subtract B's x-edge projected onto normal
          if (rotation[2][0]*sign1 < 0) {
               separation1 += rotation[2][0]*side_lengths_B[0]*sign1;
          }
          // Test if y-unit vector in B points towards A
          // if so, subtract B's y-edge projected onto normal
          if (rotation[2][1]*sign1 < 0) {
               separation1 += rotation[2][1]*side_lengths_B[1]*sign1;
          }


          int sign2;
          if (translation_BA[2] < 0) {  // A is "above" B
               sign2 = +1;
          } else {                              // A is "below" B
               sign2 = -1;
          }
          double separation2 = -translation_BA[2]*sign2;

          // Test if x-unit vector in A points towards B
          // if so, subtract A's x-edge projected onto normal
          if (rotation[0][2]*sign2 < 0) {
               separation2 += rotation[0][2]*side_lengths_A[0]*sign2;
          }
          // Test if y-unit vector in A points towards B
          // if so, subtract A's y-edge projected onto normal
          if (rotation[1][2]*sign2 < 0) {
               separation2 += rotation[1][2]*side_lengths_A[1]*sign2;
          }

          // Return maximum separation if larger than 0
          double dist = max(0.0, separation1, separation2);
     
          return (dist < cutoff_distance);
     }

     //! Overload output operator - RSS reference
     friend std::ostream& operator<<(std::ostream& out, const RSS &r) {
          for (unsigned int i=0; i<r.vertices.size(); i++) {
               out << r.vertices[i] << " ";
          }
          return out;
     }

     //! Overload output operator - RSS pointer
     friend std::ostream& operator<<(std::ostream& out, const RSS *r) {
          for (unsigned int i=0; i<r->vertices.size(); i++) {
               out << r->vertices[i] << " ";
          }
          return out;
     }
     
};


//! ChainTree Node class
//! \tparam RES_TYPE Residue type
//! \tparam BV_TYPE Bounding volume type
template <typename RES_TYPE, typename BV_TYPE>
class Node {

public:
     
     //! Node identification string
     std::string id;

     //! Depth in tree
     int level;

     //! index in node list
     int index;
     
     //! Pointer to bounding volume object 
     BV_TYPE *bv;

     //! Backup: Pointer to bounding volume object      
     BV_TYPE *bv_backup;

     //! Pointer to node that owns the bounding volume pointer     
     Node *bv_owner;

     //! Flag determining whether node is a leaf node
     bool is_leaf;

     //! Pointer to frame of reference
     Frame<RES_TYPE> *frame;

     //! Specifies if node is the owner of the frame pointer
     bool frame_owner;

     //! Atom vector (for leaf nodes)
     std::vector<Atom *> atoms;

     //! Specifies which atoms are at which index     
     int *atom_index;

     //! Vector of 3D-coordinates (local reference frame)
     std::vector<Vector_3D> local_positions;

     //! Backup: Vector of 3D-coordinates (local reference frame)
     std::vector<Vector_3D> local_positions_backup;
     
     //! Pointers to next node at same level (depth) in tree
     Node *next;

     //! Pointers to previous node at same level (depth) in tree
     Node *prev;

     //! Pointer to child node 1
     Node *child1;

     //! Pointer to child node 2
     Node *child2;

     //! Pointer to parent
     Node *parent;
     
     //! Rotation relative to neighbour in same level
     Matrix_3D rotation;

     //! Translation relative to neighbour in same level
     Vector_3D translation;

     //! Backup: Rotation relative to neighbour in same level     
     Matrix_3D rotation_backup;

     //! Backup: Translation relative to neighbour in same level
     Vector_3D translation_backup;

     //! Time stamp of last update of node
     long int time_stamp;

     //! Backup: Time stamp of last update of node     
     long int time_stamp_backup;

     //! Specifies how this node was updated at last update
     enum UpdateType {NONE=0, FRAME=1, BV=2, TRANSLATION=4, ROTATION=8, TRANSFORM=12};

     //! Overload + operator for UpdateType
     friend UpdateType operator+(UpdateType v1, UpdateType v2) {return UpdateType((int)v1 | (int)v2);}

     //! Specifies how this node was updated at last update
     UpdateType update_type; 

     //! Specifies whether the vertices of grandchildren should be used
     //! instead of the vertices of children to obtain a tighter bound.
     //! Using this option, the update time will be slightly increased,
     //! but the tighter bonds might give rise to better pruning during
     //! when finding pairs within a given cutoff distance.
     const static bool use_grandchildren = false;

     //! Constructor (leaf node)
     //!
     //! \param frame Reference frame
     //! \param frame_owner Whether the node should clean up the reference frame
     //! \param atoms Vector of atom pointers
     //! \param id Node id
     //! \param index Index of this node in enclosing list
     Node(Frame<RES_TYPE> *frame, bool frame_owner, std::vector<Atom *> &atoms,
          std::string id="LEAF", int index=-1)
          : id(id), level(0), index(index),
            bv(NULL), bv_backup(NULL), bv_owner(this), is_leaf(true),
            frame(frame), frame_owner(frame_owner), atoms(atoms),
            next(NULL), prev(NULL), child1(NULL), child2(NULL),
            parent(NULL),
            rotation(identity_matrix()), translation(null_vector()),
            time_stamp(0) {
          
          // Register atoms for quick lookup
          atom_index = new int[definitions::ATOM_ENUM_SIZE];       
          for (int i=0; i<definitions::ATOM_ENUM_SIZE; i++) {
               atom_index[i] = -1;
          }
          for (unsigned int i=0; i<atoms.size(); i++) {
               atom_index[atoms[i]->atom_type] = i;
          }

          // Initialize bounding volume
          this->update_bounding_volume();

          // Initialize update_type
          update_type = NONE;
     }

     //! Constructor (internal node)
     //!
     //! \param child1 Child node 1
     //! \param child2 Child node 2
     //! \param id Node id
     //! \param index Index of this node in enclosing list
     Node(Node *child1, Node *child2,
          std::string id="", int index=-1)
          : id(id), index(index),
            bv(NULL), bv_backup(NULL), bv_owner(this), is_leaf(false),
            frame(NULL), frame_owner(false),
            atom_index(NULL), next(NULL), prev(NULL),
            parent(NULL),
            rotation(identity_matrix()), translation(null_vector()),
            time_stamp(0) {

          // Make sure that there is always a child1, and optionally a child2
          if (child1 && child2) {
               this->child1 = child1;
               this->child1->parent = this;

               this->child2 = child2;
               this->child2->parent = this;

          } else if (child1) {
               this->child1 = child1;
               this->child1->parent = this;
               this->child2 = NULL;
          } else if (child2) {
               this->child1 = child2;
               this->child1->parent = this;
               this->child2 = NULL;
          }

          // Initialize bounding volume
          this->update_bounding_volume();     
          
          this->frame = this->child1->frame;
          this->level = this->child1->level+1;

          // Initialize update_type
          update_type = NONE;          
     }

     //! Destructor
     ~Node() {
          if (frame_owner)
               delete frame;
          if (bv_owner==this) {
               delete bv;
               if (bv_backup)
                    delete bv_backup;
          }
          delete[] atom_index;
     }

     //! Equality operator
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(Node &other) {
          bool equal = (*this->bv == *other.bv) && (time_stamp == other.time_stamp);
          return equal;
     }

     //! Return ID string
     //! \return identification string
     std::string get_id() {
          std::string output = this->id;
          if (child1 && child2) {
               output += " (C1=" + this->child1->id + " C2=" + this->child2->id + ")";
          } else if (child1) {
               output += " (C1=" + this->child1->id + ")";
          } 
          return output;
     }

     //! Overload [] indexing operator - by atom type
     //! \return Atom pointer
     Atom *operator[](const definitions::AtomEnum atom_type) const {
          return atoms.at(atom_index[atom_type]);
     }

     //! Overload [] indexing operator - by index
     //! \return Atom pointer
     Atom *operator[](int index) const {
          return atoms.at(index);
     }

     //! Return size of node (number of atoms)
     //! \return number of atoms
     unsigned int size() const {
          return atoms.size();
     }

     //! Determines whether atom is present
     //!
     //! \param atom_type Type of atom (N,CA,C,...)
     //! \return True if node has specified atom
     bool has_atom(const definitions::AtomEnum atom_type) const {
          return (atom_index[atom_type]!=-1);
     }

     //! Return index of specified atom type
     //!
     //! \param atom_type Type of atom (N,CA,C,...)
     //! \return index in internal atom vector     
     int get_atom_index(const definitions::AtomEnum atom_type) const {
          return atom_index[atom_type];
     }

     //! Update bounding volume
     void update_bounding_volume() {          

          // Only update if the node owns the bv pointer
          if (bv_owner!=this) {
               bv = this->bv_owner->bv;

          } else {

               // Delete bounding volume backup
               if (bv_backup)
                    delete bv_backup;
               bv_backup = bv;
          
               if (is_leaf) {
                    
                    // Transform to local reference frame
                    local_positions_backup = local_positions;
                    local_positions = frame->transform(atoms);

                    // Create bounding volume
                    bv = new BV_TYPE(&local_positions);

               } else {

                    if (child1 && child2) {
                    
                         // If there are grandchildren, use their bv information
                         // instead of child bv information
                         if (this->use_grandchildren &&
                             this->child1->child1) {
                         
                              BV_TYPE *bv1 = this->child1->child1->bv;
                              BV_TYPE *bv2=NULL;
                              BV_TYPE *bv3=NULL;
                              BV_TYPE *bv4=NULL;
                              Matrix_3D *rotation1 = &this->child1->child1->rotation;
                              Matrix_3D *rotation2 = NULL;
                              Matrix_3D *rotation3 = NULL;
                              Vector_3D *translation1 = &this->child1->child1->translation;
                              Vector_3D *translation2 = NULL;
                              Vector_3D *translation3 = NULL;
                              if (this->child1->child2) {
                                   bv2 = this->child1->child2->bv;
                                   bv3 = this->child2->child1->bv;
                                   if (this->child2->child2) {
                                        bv4 = this->child2->child2->bv;
                                   }
                                   rotation2 = &this->child1->child2->rotation;
                                   rotation3 = &this->child2->child1->rotation;
                                   translation2 = &this->child1->child2->translation;
                                   translation3 = &this->child2->child1->translation;
                              } else {
                                   bv2 = this->child2->child1->bv;
                                   if (this->child2->child2) {
                                        bv3 = this->child2->child2->bv;
                                   }
                                   rotation2 = &this->child2->child1->rotation;
                                   translation2 = &this->child2->child1->translation;
                              }
                              bv = new BV_TYPE(bv1, bv2, bv3, bv4,
                                               translation1, rotation1,
                                               translation2, rotation2,
                                               translation3, rotation3);

                         } else {
                              // Otherwise just construct a bv based on the two
                              // children
                              bv = new BV_TYPE(this->child1->bv, this->child2->bv,
                                               &this->child1->translation, &this->child1->rotation);
                         }

                    } else {

                         // If only one child is present, reuse the bounding volume of
                         // that child
                         bv = this->child1->bv;

                         // This bv pointer should not be freed by this node
                         bv_owner = this->child1;
                    }
               }
          }
          
          // Register that BV has been updated
          update_type = update_type + BV;
     }
     

     //! Update reference frame.
     //! NOTE: the reference frame of a node is provided for convience
     //! only, and is used only when a leaf node is updated. Therefore
     //! only frames of modified leaf nodes are updated. Since parent
     //! nodes use the frame of one of their children, the frame of a
     //! node is not always consistent with its current position. For
     //! instance, the root node has the frame of the first residue in
     //! the chain (and the first leaf node). If a move is made
     //! effecting the middle of the chain, the frame of the first leaf
     //! node will not be updated and the frame of the root node will
     //! be inconsistent.
     void update_frame() {

          if (frame_owner) {
               frame->update();

               // Register that frame has been updated
               update_type = update_type + FRAME;
          }
     }

     
     //! Update rotation matrix and translation vector
     //! relating the current node to its neighbour
     void update_transform() {

          // Don't update if there is no neighbour or
          // neighbour has the same frame
          if (!this->next || (this->next->frame == this->frame)) {
               return;
          }
          
          update_translation();
          update_rotation();

          // Register that transform has been updated
          update_type = update_type + TRANSFORM;
     }

     //! Set translation vector relating current node to its neighbour
     void update_translation() {

          // Don't update if there is no neighbour or
          // neighbour has the same frame
          if (!this->next || (this->next->frame == this->frame)) {
               return;
          }

          // Save backup (in case of undo)
          translation_backup = translation;
          
          // In the case of leaf nodes, we need to determine the relative
          // rotation and translation between the current frame and the next.
          if (is_leaf) {
               Frame<RES_TYPE> *current_frame = this->frame;
               Frame<RES_TYPE> *next_frame = this->next->frame;
               translation = transpose(*current_frame->orientation)*
                    (*next_frame->origin - *current_frame->origin);

          // In the case of internal nodes, the translation is given
          // by the rotation and translation of the child nodes
          } else if (child2) {
               translation = (child1->rotation * child2->translation) + child1->translation;
          } else {
               translation = child1->translation;
          }

          // Register that transform has been updated
          update_type = update_type + TRANSLATION;          
     }
     
     //! Set rotation matrix relating current node to its neighbour
     void update_rotation() {

          // Don't update if there is no neighbour or
          // neighbour has the same frame
          if (!this->next || (this->next->frame == this->frame)) {
               return;
          }
          
          // Save backup (in case of undo)
          rotation_backup = rotation;

          // In the case of leaf nodes, we need to determine the relative
          // rotation and translation between the current frame and the next.
          if (is_leaf) {
               
               Frame<RES_TYPE> *currentFrame = this->frame;
               Frame<RES_TYPE> *nextFrame = this->next->frame;

               rotation = transpose(*currentFrame->orientation) * *nextFrame->orientation;        

          // In the case of internal nodes, the rotation is given
          // by the rotation and translation of the child nodes
          } else if (child2) {
               rotation = child1->rotation * child2->rotation;
          } else {
               rotation = child1->rotation;
          }

          // Register that transform has been updated
          update_type = update_type + ROTATION;          
     }


     //! Undo last update
     void undo() {

          // Reset time stamp
          time_stamp = time_stamp_backup;

          // Frame update undo
          if (update_type & FRAME) {
               frame->undo();
          }

          // Translation update undo
          if (update_type & TRANSLATION) {
               translation = translation_backup;
          }

          // Rotation update undo
          if (update_type & ROTATION) {
               rotation = rotation_backup;
          }

          // Bounding volume update undo
          if (update_type & BV) {
               local_positions = local_positions_backup;
               if (bv_owner==this) {
                    if (bv)
                         delete bv;
                    bv = bv_backup;
                    bv_backup = NULL;
               } else {
                    bv = this->bv_owner->bv;
               }
          }
     }

     //! Evaluate whether the node was modified in the last update
     //!
     //! \param global_time Current global time stamp
     //! \return True if node was modified in the last update
     bool is_modified(long int global_time) {
          // std::cout << (time_stamp >= globalTime) << " " << time_stamp << " " << globalTime << " " << id << "\n";
          return time_stamp >= global_time;
     }

     //! Evaluate whether the node serves as separator.
     //! A node is a separator if it was modified in the last update
     //! and contains a transform.
     //!
     //! \param global_time Current global time stamp
     //! \return True if node is a separator
     bool is_separator(long int global_time) {
          return ((time_stamp >= global_time) &&
                  (update_type & TRANSFORM));
     }

     //! Sets a timestamp. The timestamp is used to determine whether
     //! a node was recently modified.
     //!
     //! \param global_time Current global time stamp
     void set_time_stamp(long int global_time) {
          time_stamp_backup = time_stamp;
          time_stamp = global_time;

          // Reset update status of node
          update_type = NONE;
     }
     
     //! Get neighbouring node at same level
     //!
     //! \param offset Neighbour distance (nth neighbour)
     //! \return Neighbour node pointer
     Node *get_neighbour(int offset) {
          Node *curr_node = this;
          if (offset>0) {
               while (offset>0 && curr_node) {
                    curr_node = curr_node->next;
                    offset-=1;
               }
          } else {
               while (offset<0 && curr_node) {
                    curr_node = curr_node->prev;
                    offset+=1;
               }
          }
          return curr_node;
     }

     //! Check whether node is self-consistent 
     void check_consistency() {

          // Don't check for rotation/translation consistency if there is no neighbour
          if (next) {

               // Create new frames for current and next nodes
               // (corresponds to frame->update but without
               // effecting the current frames)
               Frame<RES_TYPE> current_frame(this->frame->res);
               Frame<RES_TYPE> next_frame(this->next->frame->res);

               // Calculate the transform between the new frames
               Vector_3D translation_real = transpose(*current_frame.orientation)*
                    (*next_frame.origin - *current_frame.origin);
               Matrix_3D rotation_real = transpose(*current_frame.orientation) * *next_frame.orientation;   

               // Compate the calculated transforms with the information
               // stored in the node
               if (!(sum(fabs(translation_real - translation)) < 0.01))
                    std::cout << "\nNOT IDENTICAL at node " << id << "!!!\n\n\n";

               if (!(sum(fabs(rotation_real - rotation)) < 0.01))
                    std::cout << "\nNOT IDENTICAL at node " << id << "!!!\n\n\n";
               
               assert(sum(fabs(translation_real - translation)) < 0.01);
               assert(sum(fabs(rotation_real - rotation)) < 0.01);
          }          

          // Check BV
          if (this->is_leaf) {

               // Create new frame for current node
               // (corresponds to frame->update but without
               // effecting the current frame)
               Frame<RES_TYPE> new_frame(this->frame->res);

               // Transform to local reference frame
               std::vector<Vector_3D> local_positions = new_frame.transform(atoms);

               // Check that all positions are within the
               // current bounding volume
               for (uint i=0; i<local_positions.size(); i++) {
                    // std::vector<Vector_3D> position;
                    // position.push_back(local_positions[i]);
                    // // Create bounding volume
                    // BV_TYPE newBV(&position);

                    // Make sure that all positions are inside BV
                    double dist = bv->compute_distance(local_positions[i]);
                    assert(dist < 0.01);
               }
               
          } else {

               // In case of internal nodes, check whether
               // parents are compatible with their children (or grandchildren)
               if (child1 && child2) {
                    
                    // If there are grandchildren, use their bv information
                    // instead of child bv information
                    if (this->use_grandchildren &&
                        this->child1->child1) {
                         
                         BV_TYPE *bv1 = this->child1->child1->bv;
                         BV_TYPE *bv2=NULL;
                         BV_TYPE *bv3=NULL;
                         BV_TYPE *bv4=NULL;
                         Matrix_3D *rotation1 = &this->child1->child1->rotation;
                         Matrix_3D *rotation2 = NULL;
                         Matrix_3D *rotation3 = NULL;
                         Vector_3D *translation1 = &this->child1->child1->translation;
                         Vector_3D *translation2 = NULL;
                         Vector_3D *translation3 = NULL;
                         if (this->child1->child2) {
                              bv2 = this->child1->child2->bv;
                              bv3 = this->child2->child1->bv;
                              if (this->child2->child2) {
                                   bv4 = this->child2->child2->bv;
                              }
                              rotation2 = &this->child1->child2->rotation;
                              rotation3 = &this->child2->child1->rotation;
                              translation2 = &this->child1->child2->translation;
                              translation3 = &this->child2->child1->translation;
                         } else {
                              bv2 = this->child2->child1->bv;
                              if (this->child2->child2) {
                                   bv3 = this->child2->child2->bv;
                              }
                              rotation2 = &this->child2->child1->rotation;
                              translation2 = &this->child2->child1->translation;
                         }
                         BV_TYPE new_bv(bv1, bv2, bv3, bv4,
                                       translation1, rotation1,
                                       translation2, rotation2,
                                       translation3, rotation3);
                         
                         // Check that it matches the node's BV
                         assert(*bv == new_bv);

                    } else {

                         // Otherwise, just use the BV information of the direct children
                         BV_TYPE new_bv(this->child1->bv, this->child2->bv,
                                       &this->child1->translation, &this->child1->rotation);

                         // Check that it matches the node's BV
                         assert(*bv == new_bv);                         
                    }

               } else {

                    // If only one child was present, it should share
                    // its bounding volume with the parent
                    assert(bv == this->child1->bv);
               }               
          }
     }


     //! Overload output operator - Node reference
     friend std::ostream& operator<<(std::ostream& out, const Node &n) {
          out << *n.bv;
          return out;
     }

     //! Overload output operator - Node pointer
     friend std::ostream& operator<<(std::ostream& out, const Node *n) {
          out << *n->bv;
          return out;
     }

     
     //! Iterator for iterating over atoms in a node
     class Iterator: public std::iterator<std::random_access_iterator_tag, Atom > {

          //! Specifies how to iterate
          enum IterationMode {ATOMS_IN_RESIDUE, SELECTED_ATOMS};

          //! Specifies how to iterate
          IterationMode iteration_mode;
          
          //! Currently active atom
          Atom *current_atom;

          //! Associated node
          Node <RES_TYPE, BV_TYPE> *node;

          //! Vector of atom types to include in the iteration.
          //! Empty vector denotes all atom types.
          std::vector<definitions::AtomEnum> *selected_atoms;

          //! Vector of indexes denoting where atoms can be found in a
          //! distance matrix. Here, the vector is used only to specify
          //! which atoms are included in the selected_atoms vector.
          std::vector<int> *selected_atoms_indices;
          
          //! Choose atom based on current_index
          void assign_atom() {
               // If there is no node, return NULL
               if (!node) {
                    current_atom = NULL;

               } else if (selected_atoms) {

                    if (iteration_mode == ATOMS_IN_RESIDUE) {

                         // Iterate over all atoms, but only return selected ones
                         current_atom = NULL;
                         while (current_index < node->atoms.size()) {
                              if ((*selected_atoms_indices)[(node->atoms)[current_index]->atom_type] != -1) {
                                   current_atom = (node->atoms)[current_index];
                                   break;
                              } else {
                                   current_index++;
                              }
                         }
                         
                    } else {

                         // Iterate directly over selected_atoms vector
                         current_atom = NULL;
                         while (current_index < selected_atoms->size()) {
                              if (node->has_atom((*selected_atoms)[current_index])) {
                                   current_atom = (*node)[(*selected_atoms)[current_index]];
                                   break;
                              } else {
                                   current_index++;
                              }
                         }
                    }
               } else {

                    // If no specific atoms are selected
                    // iterate over everything
                    if (current_index < node->atoms.size()) {
                         current_atom = (node->atoms)[current_index];
                    } else {
                         current_atom = NULL;
                    }
               }
          }
          
     public:

          //! Current iteration index
          unsigned int current_index;

          //! Constructor (default)
          Iterator(): node(NULL),selected_atoms(NULL), selected_atoms_indices(NULL) {
               current_index = 0;
               current_atom = NULL;
          }

          //! Constructor
          //!
          //! \param node Associated node
          //! \param selected_atoms Vector of atom types to include in the iteration
          //! \param selected_atoms_indices Vector of indexes denoting where atoms 
          //!        can be found in a distance matrix
          Iterator(Node<RES_TYPE, BV_TYPE> *node,
                   std::vector<definitions::AtomEnum> *selected_atoms,
                   std::vector<int> *selected_atoms_indices):
               node(node),selected_atoms(selected_atoms),
               selected_atoms_indices(selected_atoms_indices) {

               current_index = 0;

               // if (selected_atoms)
               //      std::cout << *selected_atoms << "\n";
               
               if (selected_atoms) {
                    // if the selected_atoms vector is zero
                    // ignore it
                    if (selected_atoms->size()==0) {
                         this->selected_atoms=NULL;

                    // Set iteration style 
                    } else {
                         // If the selected_atoms vector is larger
                         // than the number of atoms in the node, iterate
                         // directly over the atoms
                         if (node->atoms.size() < selected_atoms->size()) {
                              iteration_mode = ATOMS_IN_RESIDUE;

                         // otherwise, iterate over the vector of selected atoms
                         } else {
                              iteration_mode = SELECTED_ATOMS;
                         }
                    }
               }
               
               // Pick an atom based on current_index
               assign_atom();
          }

          //! Copy constructor.
          //!
          //! \param other Source object from which copy is made.
          Iterator(const Iterator& other) {
               current_index = other.current_index;
               current_atom = other.current_atom;
               node = other.node;
               iteration_mode = other.iteration_mode;
               selected_atoms = other.selected_atoms;
               selected_atoms_indices = other.selected_atoms_indices;
          }
          
          //! Assignment operator
          //!
          //! \param other Source object from which assignment is made.
          //! \return Current iterator (this)
          Iterator& operator=(const Iterator& other) {
               current_index = other.current_index;
               current_atom = other.current_atom;
               node = other.node;
               iteration_mode = other.iteration_mode;
               selected_atoms = other.selected_atoms;
               selected_atoms_indices = other.selected_atoms_indices;
               return(*this);
          }
          
          //! Overload dereference operator (*)
          //!
          //! \return Current atom
          Atom *operator*() {
               return current_atom;
          }

          //! Overload increment operator (*)
          //!
          //! \return Current iterator (this)
          Iterator &operator++() {

               current_index++;
               assign_atom();

               return (*this);
          }

     };
};

//! Chain Tree datastructure for rapid detection
//! of interacting pairs in a molecule.
//!
//! \tparam CHAIN_TYPE Molecule chain type
//! \tparam BV_TYPE Bounding volume type
template <typename CHAIN_TYPE, typename BV_TYPE>
class ChainTree {

     //! Define Residue locally for ease of reference
     typedef typename CHAIN_TYPE::Residue Residue;

public:
     //! Local BV_TYPE
     typedef BV_TYPE BvType;

     //! Local Node type
     typedef Node<Residue, RSS> NodeType;
private:

     //! Protein chain
     CHAIN_TYPE *chain;

public:
     //! Vector of nodes (containing bounding volumes)
     std::vector<NodeType *> nodes;

     //! Keep track of last move type
     definitions::MoveTypeEnum last_move_type;

private:
     //! Keep track of which residues are in which nodes
     std::vector<std::vector<int> > residue_to_node_indices;

     //! Nodes updated in last update
     std::vector<Node<Residue, BV_TYPE> *> updated_nodes;

     //! Keep track of which nodes belong to which level
     std::vector<int> level_end_indices;
     
     //! Helper function to add atom to an atom vector
     //!
     //! \param atoms Destination vector of atoms
     //! \param chain Molecule chain containing atoms
     //! \param residue_index Selected residue index
     //! \param atom_type Selected atom type
     void add_atom(std::vector<Atom *> &atoms, CHAIN_TYPE *chain, 
                  int residue_index, definitions::AtomEnum atom_type) {
          Atom *atom = (*chain)(residue_index, atom_type);
          if (atom) {
               atoms.push_back(atom);
          }
     }
     
public:

     //! Reference time to keep track of updates in tree
     long int time;

     //! Whether chaintree is currently consistent with chain.
     //! The chain is sometimes temporarily out of sync (on purpose).
     bool consistent_with_chain;

     //! Constructor
     //!
     //! \param chain Molecule chain
     ChainTree(CHAIN_TYPE *chain)
          : chain(chain),
            last_move_type(definitions::NON_LOCAL),
            time(0),
            consistent_with_chain(true) {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          for (int i=0; i<chain->size(); i++) {

               // Select residue
               Residue &res = (*chain)[i];

               // Compute frame of reference for residue
               Frame<Residue> *frame = new Frame<Residue>(&res);

               // Keep track of node-residue correspondance
               residue_to_node_indices.push_back(std::vector<int>());

               // Vector of atoms for backbone node
               // note that nodes do not strictly follow
               // the residues, but are split in peptide
               // plane units.
               std::vector<Atom *> node1_atoms;
               add_atom(node1_atoms, chain, i-1, C);
               add_atom(node1_atoms, chain, i-1, O);
               add_atom(node1_atoms, chain, i  , N);
               add_atom(node1_atoms, chain, i  , CA);
               add_atom(node1_atoms, chain, i  , H);
               add_atom(node1_atoms, chain, i  , H1);
               add_atom(node1_atoms, chain, i  , H2);
               add_atom(node1_atoms, chain, i  , H3);
               add_atom(node1_atoms, chain, i  , HA);
               add_atom(node1_atoms, chain, i  , HA2);
               add_atom(node1_atoms, chain, i  , HA3);

               // This node maintains the frame
               bool frame_owner = true;

               // Create node and add it to node vector
               nodes.push_back(new Node<Residue, BV_TYPE>(frame, frame_owner, 
                                                          node1_atoms,
                                                          "leaf_"+stringify(i)+"_bb1",
                                                          nodes.size()));
               // Register node is res->node map
               residue_to_node_indices[i].push_back(nodes.size()-1);


               // vector of atoms for sidechain node
               std::vector<Atom *> node2_atoms;
               if (res.has_atom(CB)) {

                    // Iterate over sidechain atoms
                    std::vector<Atom *> forbidden_neighbours;
                    if (res.has_atom(CA))
                        forbidden_neighbours.push_back(res[CA]);
                    if (res.has_atom(N))
                        forbidden_neighbours.push_back(res[N]);
                    CovalentBondIterator<CHAIN_TYPE> it(res[CB], forbidden_neighbours);
                    for (; !it.end(); ++it) {
                         node2_atoms.push_back(&*it);
                    }
               }
               if (res.has_atom(PS)) {
                    node2_atoms.push_back(res[PS]);
               }
               if (node2_atoms.size()>0) {

                    // This node uses the frame of the backbone node
                    bool frame_owner = false;
                    
                    // Create node and add it to node vector
                    nodes.push_back(new Node<Residue, BV_TYPE>(frame, frame_owner, 
                                                               node2_atoms,
                                                               "leaf_"+stringify(i)+"_sc",
                                                               nodes.size()));
                    // Register node is res->node map
                    residue_to_node_indices[i].push_back(nodes.size()-1);
               }

               // The C-terminus has an additional node containing terminal atoms
               if (res.terminal_status == CTERM) {

                    // Vector of atoms for terminal node
                    std::vector<Atom *> node3_atoms;
                    add_atom(node3_atoms, chain, i, C);
                    add_atom(node3_atoms, chain, i, O);
                    add_atom(node3_atoms, chain, i, OXT);

                    // This node uses the frame of the backbone node
                    bool frame_owner = false;

                    if (!node3_atoms.empty()) {
                         // Create node and add it to node vector
                         nodes.push_back(new Node<Residue, BV_TYPE>(frame, frame_owner, 
                                                                    node3_atoms,
                                                                    "leaf_"+stringify(i)+"_bb2",
                                                                    nodes.size()));
                         residue_to_node_indices[i].push_back(nodes.size()-1);
                    }
               }
          }

          
          // Set connectivity in leaf nodes:
          for (unsigned int i=0; i<nodes.size(); i++) {
               if (i>0) {
                    nodes[i]->prev = nodes[i-1];
               } else {
                    nodes[i]->prev = NULL;
               }
               if (i<nodes.size()-1) {
                    nodes[i]->next = nodes[i+1];
               } else {
                    nodes[i]->next = NULL;
               }
          }


          // Set orientation and translation values
          for (unsigned int i=0; i<nodes.size(); i++) {
               if (i<nodes.size()-1) {    
                    nodes[i]->update_transform();
               } else {
                    nodes[i]->translation = Vector_3D(0,0,0);
                    nodes[i]->rotation = identity_matrix();
               }
          }

          // Register number of leaf nodes
          level_end_indices.push_back(nodes.size());
          
          // Setup internal nodes in tree
          build_hierarchy();

          // // Set time=-1 to indicate that the chaintree should be reset at 
          // // the next update
          // reset();
     }

     //! Destructor
     ~ChainTree() {
          for (uint i=0; i<nodes.size(); i++) {
               delete nodes[i];
          }
     }

     //! Equality operator.
     //!
     //! \param other Object to compare with
     //! \return True if objects are identical
     bool operator==(ChainTree &other) {
          bool equal = true;
          
          if (nodes.size() != other.nodes.size()) {
               equal = false;
          } else {
               for (unsigned int i=0; i<nodes.size(); i++) {
                    if (!(*nodes[i] == *other.nodes[i])) {
                         equal = false;
                         break;
                    }
               }
          }
          return equal;
     }

     //! Create internal nodes in tree
     void build_hierarchy() {

          int counter = 0;
          int direction = 1;
          
          Node<Residue, BV_TYPE> *curr_level = nodes[0];
          Node<Residue, BV_TYPE> *prev_node = NULL;
          Node<Residue, BV_TYPE> *curr_node = NULL;
          Node<Residue, BV_TYPE> *child_node = curr_level;
          while (curr_level->get_neighbour(direction) != NULL) {

               while(child_node != NULL) {

                    // Create new internal node
                    if (direction>0) {
                         curr_node = new Node<Residue, BV_TYPE>(child_node,
                                                               child_node->get_neighbour(direction),
                                                               "internal_"+stringify(counter),
                                                               (nodes.size()-level_end_indices.back())*direction);
                    } else {
                         curr_node = new Node<Residue, BV_TYPE>(child_node->get_neighbour(direction),
                                                               child_node,
                                                               "internal_"+stringify(counter),
                                                               (nodes.size()-level_end_indices.back())*direction);
                    }
                    nodes.push_back(curr_node);
                    counter++;

                    // // Calculate transformation to next node
                    // currNode->update_transform();

                    // Set prev/next attributes
                    if (direction>0) {
                         curr_node->prev = prev_node;
                         if (prev_node) {
                              prev_node->next = curr_node;
                              prev_node->update_transform();

                         }
                    } else {
                         curr_node->next = prev_node;
                         if (prev_node) {
                              prev_node->prev = curr_node;
                              curr_node->update_transform();
                         }
                    }

                    // Set new prevNode to currNode
                    prev_node = curr_node;

                    // Jump two steps ahead at child level
                    child_node = child_node->get_neighbour(direction*2);
               }

               // Choose the last node as starting point for next iteration
               curr_level = curr_node;

               // Reset prevNode for this level
               prev_node = NULL;

               // Set first childNode
               child_node = curr_level;

               // If negative direction, adjust indices so that all pairs
               // in the iteration have node1->index < node2->index
               // this makes it easier to look up values in the cache
               if (direction < 0) {
                    int nodes_in_level = nodes.size() - level_end_indices.back() - 1;
                    for (unsigned int i=level_end_indices.back(); i<nodes.size(); i++) {
                         nodes[i]->index += nodes_in_level;
                    }
               }
               
               // Switch direction (to get balanced trees
               direction *= -1;

               // Register number of nodes at current level
               level_end_indices.push_back(nodes.size());
          }
     }

     //! Do a full reset (update)
     //!
     //! \param time Time stamp to reset to
     void reset(int time=-1) {
          std::cout << "Chaintree: Full reset\n";

          // After a reset, no undo can be done. Set time=-1 to indicate this
          this->time = time;
          
          // Leaf nodes
          for (int i=0; i<=chain->size()-1; i++) {
               // Get indices of node containing residue atoms
               std::vector<int> node_indices = residue_to_node_indices[i];

               // Iterate over leaf nodes
               for (unsigned int j=0; j<node_indices.size(); j++) {
                    nodes[node_indices[j]]->update_frame();
               }
          }

          // All nodes
          for (unsigned int i=0; i<nodes.size(); i++) {
               nodes[i]->update_bounding_volume();
               nodes[i]->update_transform();

               // Time stamp of the nodes is set to time
               // so that all nodes will are flagged as 
               // updated after the next update. Any iteration
               // done right after a reset followed by a move
               // will therefore evaluate all pairs
               nodes[i]->set_time_stamp(time);
          }

          consistent_with_chain = true;
     }
     
     //! Update the chaintree.
     //!
     //! \param move_info Information about last move.
     void update(MoveInfo &move_info) {

          // If modified_angles is empty, do nothing
          if (move_info.modified_angles.size() == 0) {
               return;
          }

          // Type of move made
          last_move_type = move_info.move_type;
          
          // Clear vector of previously updated nodes
          updated_nodes.clear();
          
          // Increment time
          time++;

          // Reset if we are inconsistent with chain
          // This is a delayed reset of the chaintree. If a chaintree
          // is reset during an energy evaluation, it might become necessary
          // to reject the reset. The Undo method does not do the reset
          // itself, since the chain might not have been reverted to its
          // old state at that point. Instead we reset the chain at the 
          // following update.
          if (!consistent_with_chain) {
               // We now know that the chain is initialized, so if accepted, an update
               // after this update whould not be reset.
               int time = 0;
               reset(time);

               return;
          }
          
          // Node containers for first-in-first-out iteration
          std::queue<Node<Residue, BV_TYPE> *> node_queue_leaves;
          std::queue<Node<Residue, BV_TYPE> *> node_queue_internal;

          // Make sure no duplicates are registered;
          std::set<int> node_set;
          std::set<int> node_set_leaves;

          for (unsigned int range_index=0; range_index < move_info.modified_angles.size(); ++range_index) {

               // Determine range of update
               int start_index = std::max(0,               move_info.modified_angles[range_index].first  - 1);
               int end_index   = std::min(chain->size()-1, move_info.modified_angles[range_index].second + 1);

               for (int i=start_index; i<=end_index; i++) {

                    // Get indices of node containing residue atoms
                    std::vector<int> node_indices = residue_to_node_indices[i];

                    // Iterate over leaf nodes
                    for (unsigned int j=0; j<node_indices.size(); j++) {

                         Node<Residue, BV_TYPE> *node = nodes[node_indices[j]];

                         // Skip if not node has already been registered (overlapping start-end ranges)
                         if (node_set.count(node_indices[j]) > 0) {

                              // In case of overlapping regions, check whether the current node
                              // was at the very end of its region when it was registered previously
                              // if so, we now add it to nodeQueueLeaves (since we skipped it before)
                              if ((node_set_leaves.count(node_indices[j]) == 0) &&
                                  (i < end_index || j<node_indices.size()-1 || i==chain->size()-1)) {
                                   node_queue_leaves.push(node);
                                   node_set_leaves.insert(node_indices[j]);
                              }
                              continue;
                         }

                         // Register index in set to keep track of duplicates
                         node_set.insert(node_indices[j]);

                         // Keep track of which nodes were updated
                         updated_nodes.push_back(node);

                         // Set new time stamp for node
                         // NOTE: this clears the update_type value so
                         // this has to be done BEFORE any updates
                         node->set_time_stamp(time);

                         // Start by updating frames of all nodes
                         node->update_frame();

                         // Register nodes for updates (last node is not updated further - only frame)
                         if (i<end_index || j<node_indices.size()-1 || i==chain->size()-1) {
                              node_queue_leaves.push(node);
                              node_set_leaves.insert(node_indices[j]);
                         }
                    }
               }
          }

          // Iterate over leaf nodes again (frames have now been updated)
          while (node_queue_leaves.size() > 0) {
               Node<Residue, BV_TYPE> *node = node_queue_leaves.front();
               node_queue_leaves.pop();
          
               // Update both transforms and boxes
               // unless only a sidechains was modified
               node->update_bounding_volume();
               if (move_info.move_type != definitions::SIDECHAIN) {
                    node->update_transform();
               }

               // Add parent to queue
               if (node->parent) {
                    node_queue_internal.push(node->parent);
               }
          }

          // Iterate over all pushed nodes
          while (node_queue_internal.size() > 0) {

               // Pop node of queue
               Node<Residue, BV_TYPE> *node = node_queue_internal.front();
               node_queue_internal.pop();

               // Make sure that node has not alreay been updated
               if (node->time_stamp < time) {

                    updated_nodes.push_back(node);
                    
                    // Set modification time
                    node->set_time_stamp(time);

                    // No need to update transforms if move was a sidechain move
                    if (move_info.move_type == definitions::SIDECHAIN) {
                         node->update_bounding_volume();

                    // } else if (moveInfo.moveType == LOCAL) {

                    //      // THIS OPTIMIZATION DOES NOT WORK!
                    //      // No need to update transforms if the current node incorporates
                    //      // all dofs modified by move meaning only one node remains
                    //      // on the current level of the tree
                    //      // (and this one node has just been popped)
                    //      if (nodeQueueInternal.size()==0) {
                    //           node->update_bounding_volume();                         
                    //           node->update_transform();
                    //           // std::cout << "diff: " << (sum(fabs(node->translation_backup - node->translation))) << "\n";
                    //           // assert((sum(fabs(node->translation_backup - node->translation)) < 0.01));
                    //      } else {
                    //           node->update_bounding_volume();                         
                    //           node->update_transform();
                    //      }
                              
                    // Otherwise, update both boxes and transforms
                    } else {
                         node->update_bounding_volume();                         
                         node->update_transform();
                    }
                    
                    // Add parent to queue
                    if (node->parent)
                         node_queue_internal.push(node->parent);
               }
          }          
          // check_consistency();
     }

     //! Undo last update
     void undo() {

          if (time<=0) {
               // We do not have information to do an undo on the first iteration
               // so we should do a reset at this point. However, the chain might
               // not have been updated to its old state. Therefore, we do nothing
               // and leave it up to the next call to update() to reset the chain.
               time = -1;
               consistent_with_chain=false;
               return;
          }
          
          for (unsigned int i=0; i<updated_nodes.size(); i++) {
               updated_nodes[i]->undo();               
          }
     }

     //! Check self-consistency of chain tree
     void check_consistency() {

          // This output should remain so that it is clear
          // to the user if the consistency is constantly checked (very slow)
          std::cout << "Checking ChainTree consistency";

          // The chaintree can be temporarily out of sync with the chain,
          // until the next call to update()
          // If we are in this state, do not check consistency
          if (!consistent_with_chain) {
               std::cout << " - skipped (known temporary inconsistency)\n";
               return;
          } else {
               std::cout << "\n";
          }

          for (unsigned int i=0; i<nodes.size(); i++) {
               nodes[i]->check_consistency();
          }
     }
     
     //! Return height of tree
     //!
     //! \return Number of levels in tree
     int get_height() {
          return level_end_indices.size();
     }

     //! Return number of nodes at a specific level (height) in the tree.
     //!
     //! \param level Specified height (level)
     //! \return Number of nodes
     int get_nodes_at_level(int level) {
          if (level==0) {
               return level_end_indices[0];
          } else {
               return (level_end_indices[level]-level_end_indices[level-1]);
          }
     }

     //! Return start index at a specific height in the tree
     //!
     //! \param level Specified height (level)
     //! \return Start index into global node vector
     int get_level_start_index(int level) {
          if (level==0) {
               return 0;
          } else {
               return level_end_indices[level-1];
          }
     }

     //! Output to dot-format (tree-graph visualization)
     //!
     //! \param filename Output filename
     void export_to_dot_file(std::string filename) {
          std::ofstream file;
          file.open(filename.c_str());
          file << "# DOT FORMAT (viewable in graphviz)\n";
          file << "digraph G {\n";
          file << "overlap=scale;\n";
          file << "mode=major;\n";
          file << "model=shortpath;\n";
          for (unsigned int i=0; i<nodes.size(); i++) {
               if (nodes[i]->child1)
                    file << nodes[i]->id << "->" << nodes[i]->child1->id << ";\n";
               if (nodes[i]->child2)
                    file << nodes[i]->id << "->" << nodes[i]->child2->id << ";\n";
          }          
          file << "}";          
          file.close();
     }

     //! Overload output operator. 
     //! Output chaintree as string in a format plotable by Pymol.
     friend std::ostream& operator<<(std::ostream& out, const ChainTree &c) {

          // Define Pymol functions
          out << "from pymol.cgo import *\nfrom pymol import cmd\n\ndef box(id, v1, v2, v3, v4, v5, v6, v7, v8):\n\n    box = [\n       LINEWIDTH, 1.0,\n       BEGIN, LINES,\n       COLOR, 0.8, 0.8, 0.8,\n       VERTEX,  v1[0], v1[1], v1[2],\n       VERTEX,  v2[0], v2[1], v2[2],\n\n       VERTEX,  v1[0], v1[1], v1[2],\n       VERTEX,  v3[0], v3[1], v3[2],\n\n       VERTEX,  v2[0], v2[1], v2[2],\n       VERTEX,  v4[0], v4[1], v4[2],\n\n       VERTEX,  v3[0], v3[1], v3[2],\n       VERTEX,  v4[0], v4[1], v4[2],\n\n\n       VERTEX,  v5[0], v5[1], v5[2],\n       VERTEX,  v6[0], v6[1], v6[2],\n\n       VERTEX,  v5[0], v5[1], v5[2],\n       VERTEX,  v7[0], v7[1], v7[2],\n\n       VERTEX,  v6[0], v6[1], v6[2],\n       VERTEX,  v8[0], v8[1], v8[2],\n\n       VERTEX,  v7[0], v7[1], v7[2],\n       VERTEX,  v8[0], v8[1], v8[2],\n\n\n       VERTEX,  v1[0], v1[1], v1[2],\n       VERTEX,  v5[0], v5[1], v5[2],\n\n       VERTEX,  v2[0], v2[1], v2[2],\n       VERTEX,  v6[0], v6[1], v6[2],\n\n       VERTEX,  v3[0], v3[1], v3[2],\n       VERTEX,  v7[0], v7[1], v7[2],\n\n       VERTEX,  v4[0], v4[1], v4[2],\n       VERTEX,  v8[0], v8[1], v8[2],\n\n       END\n       ]\n       \n    cmd.load_cgo(box, id, 1)\n    \n\ndef sphere(id, pos, radius):\n    sphere = [\n       LINEWIDTH, 1.0,\n       BEGIN, LINES,\n       COLOR, 0.8, 0.8, 0.8,\n       SPHERE,  pos[0], pos[1], pos[2], radius,\n\n       END\n       ]\n       \n    cmd.load_cgo(sphere, id, 1)\n";

          for (unsigned int i=0; i<c.nodes.size(); i++) {

               out << "box(\"" << c.nodes[i]->get_id() << "\"";
               for (unsigned int j=0; j<c.nodes[i]->bv->vertices.size(); j++) {
                    out << ", " << c.nodes[i]->frame->transform_back(c.nodes[i]->bv->vertices[j]);
               }
               out << ");\n";

               out << "sphere(\"ref1_" << c.nodes[i]->get_id() << "\"";
               out << ", " << c.nodes[i]->frame->transform_back(c.nodes[i]->bv->reference_pos);
               out << ", 0.1);\n";

          }          
          return out;
     }
     
};


///////////////////////////////////////////////////////////////////////////////////////
//// Functionality for determining the distance between two rectangles in 3D-space //// 
///////////////////////////////////////////////////////////////////////////////////////

//! Truncate value to be within range
//!
//! \param var pointer to variable containing value
//! \param min Range minimum
//! \param max Range maximum
inline void clip_to_range(double *var, double min, double max) {
     if (*var < min) {
          *var = min;
     } else if (*var > max) {
          *var = max;
     }
}

//! Calculates the closest points on two line segments A and B. The
//! segments are defined as:
//!    Pa + a*t,   0 <= t <= lengthA
//!    Pb + b*u,   0 <= u <= lengthB
//! where Pa and Pb are one of the endpoints of edge_A and edge_B, respectively, and
//! a and b are unit vectors. 
//! The details of the calculation can be found in:
//!     Vladimir J. Lumelsky, On fast computation of distance between line segments,
//!     Information Processing Letters 21, 1985
//! Note that u and t in this formulation are bound by
//! the lengths of A and B rather than being in the range [0,1] as in the original
//! description.
inline void line_segment_distance(double *t, double *u,
                                  double length_A, double length_B,
                                  double a_dot_b, double a_dot_T,
                                  double b_dot_T) {

     // The expression for t is rewritten compared to the paper:
     //   t = ((S1D1)/(D1D2) - (S2R)/(D1D2)) / ((1-R^2)/(D1D2))
     // In our case:
     // R=A*B, S1=A*t, S2=B*t, D1=A*A=|A|^2, D2=B*B=|B|^2
     //   t = ((a*t)/|A| - (a*b)(b*t)/|A|) / (1-(a*b)(a*b))
     // here we express t in units of |A|
     
     double denominator = 1 - (a_dot_b)*(a_dot_b);

     // Calculate point t on A
     if (denominator < 1e-7) {
          *t = 0;
     } else {
          *t =  (a_dot_T - b_dot_T*a_dot_b)/denominator;
          clip_to_range(t, 0, length_A);
     }

     // Calculate point u on B based on t
     *u = *t*a_dot_b - b_dot_T;

     if (*u < 0) {
          *u = 0;

          // Set final value of t
          *t = a_dot_T;
          clip_to_range(t, 0, length_A);
          
     } else if (*u > length_B) {
          *u = length_B;

          // set final value of t
          *t = (*u)*a_dot_b + a_dot_T;
          clip_to_range(t, 0, length_A);
     }
}


//! Determines whether two edges are in each others Voronoi regions
//!
//! \param edge_A_length Length of active edge of box A
//! \param edge_B_length Length of active edge of box B
//! \param edge_A_norm_dot_edge_B The norm of the edge of box A dotted with the edge of box B
//! \param edge_A_norm_dot_T The norm of the edge of box A dotted with the translation
//! \param edge_A_dot_edge_B The edge of box A dotted with the edge of B
//! \param edge_A_dot_T The edge of box A dotted with the translation
//! \param edge_B_dot_T The edge of box B dotted with the translation
//! \return Whether two edges are in each others Voronoi regions
inline bool in_voronoi(double edge_A_length, double edge_B_length,
                       double edge_A_norm_dot_edge_B, double edge_A_norm_dot_T,
                       double edge_A_dot_edge_B, double edge_A_dot_T, double edge_B_dot_T) {

     // Return false if edges are parallel
     if (std::fabs(edge_A_norm_dot_edge_B) < 1e-7) {
          return false;
     }

     // Calculate u, the length of the hypotenuse in the triangle
     // edge_A_norm_dot_T, with cos(theta)=-edge_A_norm_dot_edge_B
     // Measure where edge_B intersects with the plane defined by edge_A_norm
     //     
     //     '.  u
     //     . '.
     //     .   '  b
     //     ...../;.
     //  ___'   / /
     // |  a|  / /
     // |___|   '
     double u = -edge_A_norm_dot_T / edge_A_norm_dot_edge_B;

     // Ensure that u is within [0,b]
     clip_to_range(&u, 0, edge_B_length);

     // Based on the position u, calculate position t on edge_A with
     // shortest distance to the point defined by u. Uses Lumelsky's
     // algorithm.
     double t = u*edge_A_dot_edge_B + edge_A_dot_T;

     // Ensure that t is within [0,a]
     clip_to_range(&t, 0, edge_A_length);

     // Calculate new position v on edge_B, this time finding the point
     // with shortest distance to the position defined by t
     double v = t*edge_A_dot_edge_B - edge_B_dot_T;

     // Check that the new position v is "smaller" than u
     if (edge_A_norm_dot_edge_B > 0) {
          if (v > (u + 1e-7)) {
               return true;
          }
     } else {
          if (v < (u - 1e-7)) {
               return true;
          }
     }
     return false;
}

//! Checks whether two edges are in eachothers external Voronoi regions - first by doing
//! some simple tests, and doing a full in_voronoi check if the simple tests are
//! inconclusive. Finds the shortest vector between edges if the test is positive. 
//!
//! \param shortest_vector The resulting vector between edges
//! \param edge_A_ref_B_vertex1 position of endpoint 1 of an edge of A in B's reference
//!                             frame
//! \param edge_A_ref_B_vertex2 position of endpoint 2 of an edge of A in B's reference
//!                             frame
//! \param edge_B_ref_A_vertex1 position of endpoint 1 of an edge of B in A's reference
//!                             frame
//! \param edge_B_ref_A_vertex2 position of endpoint 2 of an edge of B in A's reference
//!                             frame
//! \param edge_A_ref_A The positions of the edges in their own reference frame (where both
//!                     endpoints by construction will have the same value)
//! \param edge_A_direction Keep track of the outward direction of edge_A using this sign (-1|1)
//! \param edge_A_length length of edge_A
//! \param edge_B_ref_B the positions of the edges in their own reference frame (where both
//!                      endpoints by construction will have the same value)
//! \param edge_B_direction Keep track of the outward direction of edge_B using this sign (-1|1)
//! \param edge_B_length  length of edge_B
//! \param edge_A_dot_X whether A is along X axis
//! \param edge_A_dot_edge_B Dot product between edge_A vector and edge_B vector
//! \param edge_A_dot_edge_B_norm Dot product between edge_A and the norm of edge_B
//! \param edge_A_norm_dot_edge_B Dot product between the norm of edge_A and edge_B
//! \param edge_A_norm_dot_edge_B_norm Dot product between the two norms
//! \param translation_AB_dot_edge_A Dot product between translation_AB vector and edge_A
//! \param translation_AB_dot_edge_A_norm Dot product between translation_AB vector and norm of edge_A
//! \param translation_BA_dot_edge_B Dot product between translation_BA vector and edge_B
//! \param translation_BA_dot_edge_B_norm Dot product between translation_BA vector and norm of edge_B
//! \param translation_AB Complete translation vector
//! \param edge_B_ref_A edge_B vector in A's reference system
//! \param edge_B_norm_ref_A Norm of edge_B vector in A's reference system
inline bool PHAISTOS_ALWAYS_INLINE edge_pair_check(Vector_3D &shortest_vector,
                                                   double edge_A_ref_B_vertex1, double edge_A_ref_B_vertex2,
                                                   double edge_B_ref_A_vertex1, double edge_B_ref_A_vertex2,
                                                   double edge_A_ref_A, int edge_A_direction, double edge_A_length,
                                                   double edge_B_ref_B, int edge_B_direction, double edge_B_length,
                                                   int edge_A_dot_X,

                                                   double edge_A_dot_edge_B, double edge_A_dot_edge_B_norm,
                                                   double edge_A_norm_dot_edge_B, double edge_A_norm_dot_edge_B_norm,

                                                   double translation_AB_dot_edge_A, double translation_AB_dot_edge_A_norm,
                                                   double translation_BA_dot_edge_B, double translation_BA_dot_edge_B_norm,
                                                   const Vector_3D &translation_AB, const Vector_3D &edge_B_ref_A,
                                                   const Vector_3D &edge_B_norm_ref_A) {
     
     // Find lower and upper vertex for A
     double edge_A_ref_B_lower = edge_A_ref_B_vertex1;
     double edge_A_ref_B_upper = edge_A_ref_B_vertex2;
     if (edge_A_ref_B_vertex1 > edge_A_ref_B_vertex2) {
          edge_A_ref_B_lower = edge_A_ref_B_vertex2;
          edge_A_ref_B_upper = edge_A_ref_B_vertex1;
     }

     // Find lower and upper vertex for B
     double edge_B_ref_A_lower = edge_B_ref_A_vertex1;
     double edge_B_ref_A_upper = edge_B_ref_A_vertex2;
     if (edge_B_ref_A_vertex1 > edge_B_ref_A_vertex2) {
          edge_B_ref_A_lower = edge_B_ref_A_vertex2;
          edge_B_ref_A_upper = edge_B_ref_A_vertex1;
     }

     // if edge_B's direction is negative, switch lower and upper on A
     if (edge_B_direction<0) {
          double tmp = edge_A_ref_B_lower;
          edge_A_ref_B_lower = edge_A_ref_B_upper;
          edge_A_ref_B_upper = tmp;
     }
     
     // if edge_A's direction is negative, switch lower and upper on B
     if (edge_A_direction<0) {
          double tmp = edge_B_ref_A_lower;
          edge_B_ref_A_lower = edge_B_ref_A_upper;
          edge_B_ref_A_upper = tmp;
     }
     
     // Check if upper vertices are "outside" edge in both reference systems
     if ((((edge_A_ref_B_upper - edge_B_ref_B)*edge_B_direction) > 0) &&
         (((edge_B_ref_A_upper - edge_A_ref_A)*edge_A_direction) > 0)) {

          // Check if lower vertices are "outside" edge or do full in_voronoi check 
          if (((((edge_A_ref_B_lower - edge_B_ref_B)*edge_B_direction) > 0) ||
               in_voronoi(edge_B_length, edge_A_length,
                         edge_A_dot_edge_B_norm*edge_B_direction,
                         (edge_A_norm_dot_edge_B_norm*edge_A_ref_A -
                          edge_B_ref_B - translation_BA_dot_edge_B_norm)*edge_B_direction,
                         edge_A_dot_edge_B, edge_A_norm_dot_edge_B*edge_A_ref_A - translation_BA_dot_edge_B,
                         -translation_AB_dot_edge_A - edge_A_dot_edge_B_norm*edge_B_ref_B))
              &&
              ((((edge_B_ref_A_lower - edge_A_ref_A)*edge_A_direction) > 0) ||
               in_voronoi(edge_A_length, edge_B_length,
                         edge_A_norm_dot_edge_B*edge_A_direction,
                         (edge_A_norm_dot_edge_B_norm*edge_B_ref_B -
                          edge_A_ref_A + translation_AB_dot_edge_A_norm)*edge_A_direction,
                         edge_A_dot_edge_B, edge_A_dot_edge_B_norm*edge_B_ref_B + translation_AB_dot_edge_A,
                         translation_BA_dot_edge_B - edge_A_norm_dot_edge_B*edge_A_ref_A))) {

               // If edge is in external Voronoi region, compute distance between segments
               double t, u;
               line_segment_distance(&t, &u, edge_A_length, edge_B_length,
                                   edge_A_dot_edge_B,
                                   edge_A_dot_edge_B_norm*edge_B_ref_B + translation_AB_dot_edge_A,
                                   translation_BA_dot_edge_B - edge_A_norm_dot_edge_B*edge_A_ref_A);

               // Edge A is along Y axis
               if (edge_A_dot_X == 0) {
                    shortest_vector[0] = translation_AB[0] + edge_B_norm_ref_A[0]*edge_B_ref_B + edge_B_ref_A[0]*u -
                         edge_A_ref_A;
                    shortest_vector[1] = translation_AB[1] + edge_B_norm_ref_A[1]*edge_B_ref_B + edge_B_ref_A[1]*u - t;
                    shortest_vector[2] = translation_AB[2] + edge_B_norm_ref_A[2]*edge_B_ref_B + edge_B_ref_A[2]*u;

               // Edge A is along X axis 
               } else {
                    shortest_vector[0] = translation_AB[0] + edge_B_norm_ref_A[0]*edge_B_ref_B + edge_B_ref_A[0]*u - t;
                    shortest_vector[1] = translation_AB[1] + edge_B_norm_ref_A[1]*edge_B_ref_B + edge_B_ref_A[1]*u -
                         edge_A_ref_A;
                    shortest_vector[2] = translation_AB[2] + edge_B_norm_ref_A[2]*edge_B_ref_B + edge_B_ref_A[2]*u;
               }
               
               return true;
          }
     }
     return false;
}



////////////////////////////////////////////////////////////
//// Simpler functionality for point-vs-rectangle tests ////
////////////////////////////////////////////////////////////


//! Determine distance between a line segment and a point
inline void line_segment_point_distance(double *t, 
                                        double lengthA,
                                        double a_dot_T) {
     *t = a_dot_T;
     clip_to_range(t, 0, lengthA);     
}


//! Checks whether a point is in the external Voronoi region of an edge.
//! Finds the shortest vector between edges if the test is positive. 
//!
//! \param shortest_vector The resulting vector between edges
//! \param pos_B_ref_A Point's position in A's reference frame
//! \param edge_A_ref_A The positions of edge_A in its own reference frame
//! \param edge_A_length Length of edge_A
//! \param edge_A_direction Keep track of the outward direction of edge_A using this sign (-1|1)
//! \param edge_A_dot_X Whether A is along X axis
//! \param translation_AB_dot_edge_A Dot product between translation_AB vector and edge_A
//! \param translation_AB Complete translation vector
//! \return Whether point is in external Voronoi region of edge
inline bool edge_point_check(Vector_3D &shortest_vector,
                             double pos_B_ref_A,
                             double edge_A_ref_A, int edge_A_direction, double edge_A_length,
                             int edge_A_dot_X,
                             double translation_AB_dot_edge_A,
                             const Vector_3D &translation_AB) {

     if (((pos_B_ref_A - edge_A_ref_A)*edge_A_direction) > 0) {

          double t;
          line_segment_point_distance(&t, edge_A_length,
                                   translation_AB_dot_edge_A);

          // Edge A is along Y axis
          if (edge_A_dot_X == 0) {
               shortest_vector[0] = translation_AB[0] - edge_A_ref_A;
               shortest_vector[1] = translation_AB[1] - t;
               shortest_vector[2] = translation_AB[2];

               // A is along X axis 
          } else {
               shortest_vector[0] = translation_AB[0] - t;
               shortest_vector[1] = translation_AB[1] - edge_A_ref_A;
               shortest_vector[2] = translation_AB[2];
          }

          return true;
     }
     return false;
}


//! Determine the distance between a rectangles and a point in 3D-space
//!
//! \param translation_AB Translation vector relating reference frame A and B
//! \param rotation_AB Rotation matrix relatinv reference frame A and B
//! \param side_lengths_A (A's sidelength in x-direction, A's sidelength in y-direction)
//! \return Distance between rectangle and point
inline double rect_distance(Vector_3D &translation_AB, Matrix_3D &rotation_AB, double *side_lengths_A) {

     // Compute translation vector in B's frame of reference
     Vector_3D translation_BA = transpose(rotation_AB)*translation_AB;

     // column vectors in rotation matrix
     Vector_3D rotation_AB0 = rotation_AB.col_vector(0);
     Vector_3D rotation_AB1 = rotation_AB.col_vector(1);
     Vector_3D rotation_AB2 = rotation_AB.col_vector(2);
     
     // Projection of B's vertices into in A's reference frame: x-coordinates
     double B_X = translation_AB[0];
     double B_Y = translation_AB[1];
     
     Vector_3D shortest_vector;

     // Edges in rectangle:
     //   ^
     // LU|____ UU
     //   |    |
     //   |____|____>
     // LL      UL
     
     // 1. check edge A: UL-UU vs B: UL-UU
     if (edge_point_check(shortest_vector,
                        B_X,
                        side_lengths_A[0], +1, side_lengths_A[1],
                        0,
                        translation_AB[1],
                        translation_AB)) {
          // std::cout << "rect_distance: 1\n";
          return sqrt(shortest_vector*shortest_vector);
     }
         
     // 2. check edge A: LL-LU vs B: UL-UU
     if (edge_point_check(shortest_vector,
                        B_X,
                        0,               -1, side_lengths_A[1],
                        0,
                        translation_AB[1],
                        translation_AB)) {
          // std::cout << "rect_distance: 2\n";
          return sqrt(shortest_vector*shortest_vector);
     }
     
     // 3. check edge A: LU-UU vs B: UL-UU
     if (edge_point_check(shortest_vector,
                        B_Y,
                        side_lengths_A[1], +1, side_lengths_A[0],
                        1,
                        translation_AB[0],
                        translation_AB)) {
          // std::cout << "rect_distance: 3\n";
          return sqrt(shortest_vector*shortest_vector);
     }

     // 4. check edge A: LL-UL vs B: UL-UU
     if (edge_point_check(shortest_vector,
                        B_Y,
                        0,               -1, side_lengths_A[0],
                        1,
                        translation_AB[0],
                        translation_AB)) {
          // std::cout << "rect_distance: 4\n";
          return sqrt(shortest_vector*shortest_vector);
     }

     // If none of the above tests were true, calculate maximum
     // distance along rectangle face normal vectors

     int sign1;
     if (translation_AB[2] > 0) {       // B is "above" A
          sign1 = +1;
     } else {                           // B is "below" A
          sign1 = -1;
     }
     double separation1 = translation_AB[2]*sign1;
          
     int sign2;
     if (translation_BA[2] < 0) {       // A is "above" B
          sign2 = +1;
     } else {                           // A is "below" B
          sign2 = -1;
     }
     double separation2 = -translation_BA[2]*sign2;

     // Test if x-unit vector in A points towards B
     // if so, subtract A's x-edge projected onto normal
     if (rotation_AB[0][2]*sign2 < 0) {
          separation2 += rotation_AB[0][2]*side_lengths_A[0]*sign2;
     }
     // Test if y-unit vector in A points towards B
     // if so, subtract A's y-edge projected onto normal
     if (rotation_AB[1][2]*sign2 < 0) {
          separation2 += rotation_AB[1][2]*side_lengths_A[1]*sign2;
     }
     
     // Return maximum separation if larger than 0
     return max(0.0, separation1, separation2);
}

}
}

#endif
