// sphere.h --- Sphere class for calculating visible volume derived measures
// Copyright (C) 2011-2013 Kristoffer Enøe Johansson
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



// calc_res/atom_contacts: consider passing a vector pointer to avoid copying the whole thing on return
// window: better variable names

#ifndef VISIBLE_VOLUME_SPHERE_H
#define VISIBLE_VOLUME_SPHERE_H

#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>

// for debug
#include <iostream>

//! Class to contain all properties of a single angle point of the sphere
class AnglePoint {
     
private:

     //! Unit vectior of the angle point
     std::vector<double> unit_vector;
     
     //! Distance from the center to the surface element before any neighbours are added
     double initial_distance;
     
     //! Is angle point cover by a neighbour
     bool covered;
     
     //! Distance to nearest neighbour
     double distance;
     
     //! Residue index and type of nearest neighbour
     int res_index, res_type;
     
     //! Atom index and type of nearest neighbour
     int atom_index, atom_type;
     
public:
     
     //! Window through which this angle point is seen
     int window;
     
     //! Default Constructor
     AnglePoint() {
          unit_vector = std::vector<double>(0);
          initial_distance = -1.0;
          reset();
     };
     
     //! Main Constructor - always pass unit vector and initial distance
     AnglePoint(std::vector<double> unit_vector, double initial_distance): unit_vector(unit_vector),
                                                                           initial_distance(initial_distance) {
          assert(unit_vector.size() == 3);
          reset();
     };
     
     // //! Copy constructor - kill state upon copy
     // AnglePoint(const AnglePoint &other): unit_vector(other.unit_vector),
     //                                      initial_distance(initial_distance) {
     // };
     
     //! Destructor - just default destroy everything
     ~AnglePoint() {};
     
     void reset() {
          this->covered = false;
          this->res_index = 0;
          this->atom_index = 0;
          this->distance = initial_distance;
          this->res_type = 0;
          this->atom_type = 0;
          this->window = 0;
     };
     
     inline void cover(int residue_index, int atom_index, double distance) {
          this->covered = true;
          this->res_index = residue_index;
          this->atom_index = atom_index;
          // this should be checked elseware
          assert(distance < this->distance);
          this->distance = distance;
     };
     
     // getters
     inline bool is_covered() { return this->covered; }
     inline std::vector<double> *get_unit_vec() { return &(this->unit_vector); }
     inline double get_distance() { return this->distance; }
     inline int get_res_index() { return this->res_index; }
     inline int get_res_type() { return this->res_type; }
     inline int get_atom_index() { return this->atom_index; }
     inline int get_atom_type() { return this->atom_type; }
     
     // setters
     // atom_type and res_type is non-negative by definition
     inline void set_atom_type(int atom_type) {
          assert(atom_type>=0);
          this->atom_type = atom_type;
     }
     
     inline void set_res_type(int res_type) {
          assert(res_type>=0);
          this->res_type = res_type;
     }
};

class Sphere {
     
public:
          
     //!Cached quantity
     double two_pi, cos_four_pi_third;
     
     //! Radius of the sphere
     double radius;

     //! Actual number of points on the sphere surface
     int n_points;

     //! The number of different windows that the angle points can be assigned to
     int n_windows;

     //! Count how many shadows are added to the object
     int shade_counter;

     // Surface points are distributed approximately uniformly by a construction based
     // on the spherical coordinates, theta (azimuth angle from positive z-axis in [0,pi])
     // and phi (polar angle from positive x-axis in [0,2pi]):
     // *d_theta is constant
     // *theta boundaries are in 0, 1*dt, 2*dt, ..., pi
     // *point theta values are at 1/2*dt, 3/2*dt, ..., pi-1/2*dt
     // *d_phi scales with sin(d_theta) to get points distributed as uniformly as possible
     //
     // 0                   1 x d_phi[0]        2 x d_phi[0]        3 x d_phi[0] = 2pi
     // |-------------------|-------------------|-------------------| 0
     // |     phi[0][0]     |         p         |         p         | theta[0] = 0.5 x d_theta
     // |--------------|----|---------|---------|----|--------------| 1 x d_theta
     // |  phi[1][0]   |       p      |      p       |      p       | theta[1] = 1.5 x d_theta
     // |-----------|--|--------|-----|-----|--------|--|-----------| 2 x d_theta
     // |     p     |     p     |     p     | phi[2][3] |     p     | theta[2] = 2.5 x d_theta
     // |-----------|--|--------|-----|-----|--------|--|-----------| 3 x d_theta
     // |       p      |       p      |      p       |       p      | theta[3] = 3.5 x d_theta
     // |--------------|----|---------|---------|----|--------------| 4 x d_theta
     // |         p         |         p         |         p         | theta[4] = 4.5 x d_theta
     // |-------------------|-------------------|-------------------| 5 x d_theta = pi
     // 
     
     //! List of theta values
     std::vector<double> theta_list;
     
     //! List of lists containing phi values for each theta value
     std::vector<std::vector<double> > phi_list_list;
     
     //! Constant theta variation of each point
     double d_theta;
     
     //! Theta dependent phi vaiation
     std::vector<double> d_phi_list;

     //! The angle points
     std::vector<std::vector<AnglePoint> > angle_point;
     
     //! Constructor
     Sphere(double radius, int min_points) {

          assert(radius > 0.0);
          assert(min_points > 3);
          
          // init cached values
          this->two_pi = 2*M_PI;
          this->cos_four_pi_third = cos(4.0*M_PI/3.0);
          
          // init other
          this->radius = radius;
          this->n_points = 0;
          this->n_windows = 1;
          this->shade_counter = 0;
          
          // find point coordinates
          double unit_sphere_surface = 4.0*M_PI;
          // ideal area pr point on unit sphere
          double surface_element = unit_sphere_surface/min_points;
          // angle element for square surface element at equator
          this->d_theta = 2*asin(sqrt(surface_element)/2.0);
          // integer number of intervals
          int theta_intervals = int(ceil(M_PI/this->d_theta));
          
          // there's not really any reason for this but I think it gives a better distribution
          // even number of intervals => halfsphere symmetry
          // if ( theta_intervals%2 != 0 )
          //      theta_intervals += 1;
          
          // reset angle element to ensure boundaries in 0 and pi
          this->d_theta = M_PI/theta_intervals;

          // square elements at equator
          double equatorial_d_phi = two_pi/this->d_theta;
          for (int t=0; t<theta_intervals; t++) {
               // put points in the middle of each interval
               double theta = (t+0.5)*this->d_theta;
               this->theta_list.push_back(theta);
               // int n_phi = int(round(equatorial_d_phi*sin(theta)));
               int n_phi = int(ceil(equatorial_d_phi*sin(theta)));
               this->d_phi_list.push_back(two_pi/n_phi);
               std::vector<double> phi_list(0);
               std::vector<AnglePoint> angle_point_buffer(0);
               for (int p=0; p<n_phi; p++) {
                    double phi = (p+0.5)*this->d_phi_list[t];
                    phi_list.push_back(phi);
                    std::vector<double> point_unit_vec = sphere_to_cart(1.0, theta, phi);
                    // distance from center to pyramid base and not sphere surface
                    double corner_theta = 0.0;
                    if (theta < M_PI/2.0)
                         corner_theta = theta - this->d_theta/2.0;
                    else
                         corner_theta = theta + this->d_theta/2.0;
                    double corner_phi = phi - this->d_phi_list[t]/2.0;
                    std::vector<double> corner_unit_vec = sphere_to_cart(1.0, corner_theta, corner_phi);
                    double cos_ang = this->dot(point_unit_vec, corner_unit_vec);
                    double initial_distance = this->radius*cos_ang;
                    angle_point_buffer.push_back(AnglePoint(point_unit_vec, initial_distance));
               }
               this->n_points += phi_list.size();
               this->phi_list_list.push_back(phi_list);
               this->angle_point.push_back(angle_point_buffer);
          }
     }
     
     //! Destructor
     ~Sphere() {};
     
     //! Shade the sphere surface (looking from the origin, i.e contact center) with a 
     //! disk of radius 'radius' at position 'position' relative to the contact center.
     //! The shade is marked with atom type and sequence index.
     void shade(double *position, double disk_radius, int res_index, int atom_index,
                int res_type = -1, int atom_type = -1) {

          // check disk radius
          assert( disk_radius>0.0 );

          // std::cout<<"resi "<<res_index<<" type "<<res_type<<" and atomi "<<atom_index<<" type "<<atom_type<<std::endl;
          // std::cout<<position[0]<<", "<<position[1]<<", "<<position[2]<<std::endl;
          assert(position);
          // std::cout<<"assertion ok"<<std::endl;

          // transform disk position to spherical coordinates
          std::vector<double> disk_vector = this->cart_to_sphere(position[0], position[1], position[2]);

          // check disk distance to contact center
          if (disk_vector[0] > this->radius)
               return;
          
          // registre shade
          this->shade_counter += 1;
          
          // shadow cone angle
          double cone_angle = atan(disk_radius/disk_vector[0]);
          double cos_cone_angle = cos(cone_angle);
          // printf("Cone angle: %8.3f\n",cone_angle*180.0/M_PI);
          
          // unit vector of disk center - faster using sphere_to_cart(1.0,t,p)?
          std::vector<double> disk_unit_vec(position,position+3);
          normalize(disk_unit_vec);

          // find theta index range to search for shaded points
          int start_theta_index, end_theta_index;
          double start_theta = disk_vector[1] - cone_angle;
          if (start_theta < 0.0)
               start_theta_index = 0;
          else
               start_theta_index = int( floor(start_theta/this->d_theta) );
          double end_theta = disk_vector[1] + cone_angle;
          if (end_theta > M_PI)
               end_theta_index = this->theta_list.size();
          else
               end_theta_index = int(ceil( end_theta/this->d_theta ));
          
          // search through phi points in positive and negative direction
          for (int theta_index=start_theta_index; theta_index<end_theta_index; theta_index++) {
               // search phi in both directions until out of shadow
               // int start_phi_index = int(ceil( disk_vector[2]/this->d_phi_list[theta_index]));
               int start_phi_index = int(floor( disk_vector[2]/this->d_phi_list[theta_index]));
               int phi_indices = this->phi_list_list[theta_index].size();
               if (start_phi_index == phi_indices)
                    start_phi_index -= 1; // in case disk_vector[2] is _exactly_ 2pi
               int covered_points = 0;

               // search verbose
               // printf("Search t=%2d from p=%2d+ and p=%2d-, lpi %d\n",
               //        theta_index, start_phi_index, start_phi_index-1, phi_indices);
               
               bool in_shadow = true;
               int phi_index = start_phi_index;
               // while (in_shadow  &&  phi_index < int(this->phi_list_list[theta_index].size())) {
               while (in_shadow  &&  covered_points < phi_indices) {
                    std::vector<double> *point_unit_vec = this->angle_point[theta_index][phi_index].get_unit_vec();
                    // angle from disk center vector to point vector
                    double cos_point_angle = dot(*point_unit_vec, disk_unit_vec);
                    if (cos_point_angle > cos_cone_angle && cos_point_angle > 0.0) {
                         covered_points += 1;
                         double dist = disk_vector[0]/cos_point_angle;
                         
                         // search verbose
                         // printf("t=%2d  p=%2d cover %3d by res index %3d atom index %2d\n",
                         //        theta_index, phi_index, covered_points, res_index, atom_index);

                         if (this->angle_point[theta_index][phi_index].get_distance() > dist) {
                              this->angle_point[theta_index][phi_index].cover(res_index, atom_index, dist);
                              if (atom_type >= 0)
                                   this->angle_point[theta_index][phi_index].set_atom_type(atom_type);
                              if (res_type >= 0)
                                   this->angle_point[theta_index][phi_index].set_res_type(res_type);
                         }
                    } else {
                         in_shadow = false;
                    }
                    // continue from p=0 if last phi index i reached
                    if (phi_index < (phi_indices-1))
                         phi_index += 1;
                    else
                         phi_index = 0;
               }

               // search negative phi direction
               in_shadow = true;
               
               // start from last phi index if start phi index is zero
               if (start_phi_index > 0)
                    phi_index = start_phi_index-1;
               else
                    phi_index = phi_indices-1;
               // while (in_shadow && phi_index>=0) {
               while (in_shadow && covered_points < phi_indices) {
                    std::vector<double> *point_unit_vec = this->angle_point[theta_index][phi_index].get_unit_vec();
                    double cos_point_angle = dot(*point_unit_vec, disk_unit_vec);
                    if (cos_point_angle > cos_cone_angle && cos_point_angle > 0.0) {
                         covered_points += 1;
                         double dist = disk_vector[0]/cos_point_angle;

                         // search verbose
                         // printf("t=%2d  p=%2d cover %3d by res index %3d atom index %2d\n",
                         //        theta_index, phi_index, covered_points, res_index, atom_index);

                         if (this->angle_point[theta_index][phi_index].get_distance() > dist) {
                              this->angle_point[theta_index][phi_index].cover(res_index, atom_index, dist);
                              if (atom_type >= 0)
                                   this->angle_point[theta_index][phi_index].set_atom_type(atom_type);
                              if (res_type >= 0)
                                   this->angle_point[theta_index][phi_index].set_res_type(res_type);
                         }
                    } else {
                         in_shadow = false;
                    }
                    // continue from last phi index if index zero is reached
                    if (phi_index > 0)
                         phi_index -= 1;
                    else
                         phi_index = phi_indices-1;
               }
          }
     }
     
     //! Assign all angle points to the same window. Default.
     void assign_single_window() {
          this->n_windows = 1;

          // assign all angle points
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    this->angle_point[t][p].window = 0;
               }
          }
     }

     //! Assign window according to the "half sphere exposure" model described in
     //! T. Hamelryck: "An Amino Acid Has Two Sides"
     //! PROTEINS: Structure, Function, and Bioinformatics 59:38-48 (2005)
     // /param ca_cc_vec pointer to vector from Ca atom to ContactCenter
     void assign_halfsphere_window(std::vector<double> *ca_cc_vec) {
          this->n_windows = 2;

          // assign all angle points
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    std::vector<double> *point_unit_vec = angle_point[t][p].get_unit_vec();
                    double projection = dot(*ca_cc_vec, *point_unit_vec);
                    if (projection > 0)
                         this->angle_point[t][p].window = 0;
                    else
                         this->angle_point[t][p].window = 1;
               }
          }
     }

     //! Assign tetrahedral window model
     // /param ca_cc_vec pointer to vector from Ca atom to ContactCenter
     // /param ca_c_vec pointer to vector from Ca to backbone C atom
     void assign_tetrahedral_window(std::vector<double> ca_cc_vec, std::vector<double> ca_c_vec) {

          this->n_windows = 4;

          // Normalize vectors
          normalize(ca_cc_vec);
          normalize(ca_c_vec);
          
          // ca_c vector in plane perpendicular on ca_cc_vec (window plane, wp)
          std::vector<double> wp_zero = cross(ca_cc_vec,ca_c_vec);
          assert( (wp_zero[0]*wp_zero[0]+wp_zero[1]*wp_zero[1]+wp_zero[2]*wp_zero[2]) > 1E-6);
          normalize(wp_zero);


          // assign all angle points
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    std::vector<double> *point_unit_vec = angle_point[t][p].get_unit_vec();
                    // check for backbone window
                    double cos = dot(ca_cc_vec, *point_unit_vec);
                    if (cos < this->cos_four_pi_third) {
                         this->angle_point[t][p].window = 3;
                         continue;
                    }
                    // check for H-N window
                    std::vector<double> wp_point = cross(ca_cc_vec, *point_unit_vec); // length sin(angle)
                    normalize(wp_point);
                    // unsigned angle
                    cos = dot(wp_point, wp_zero);
                    if (cos < this->cos_four_pi_third) {
                         // H-N window
                         this->angle_point[t][p].window = 1;
                         continue;
                    }
                    // check for C-H or N-C window
                    std::vector<double> wp_zero_point_cross = cross(wp_zero, wp_point);
                    // parallel or anti-parallel to CA-CC vector?
                    double direction_vector = dot(wp_zero_point_cross, ca_cc_vec);
                    if (direction_vector > 0) {
                         // C-H window
                         this->angle_point[t][p].window = 0;
                    } else {
                         // N-C window
                         this->angle_point[t][p].window = 2;
                    }
               }
          }
     }
          
     //! Calculate the total visible volume of the current state of the object
     std::vector<double> calc_visible_volume() {
          // return vector
          std::vector<double> vv(this->n_windows, 0.0);
          
          // sum pyramid volumes for each angle points
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    double r = this->angle_point[t][p].get_distance();
                    int w = this->angle_point[t][p].window;

                    // pyramid volume for each point
                    double theta_base = 2.0*r*tan(this->d_theta/2.0);
                    double phi_base = 2.0*r*tan(this->d_phi_list[t]/2.0)*sin(this->theta_list[t]);
                    double base_area = theta_base*phi_base;
                    double vol = base_area*r/3.0;

                    // assign to window volume of that angle point
                    vv[w] += vol;
               }
          }
          return vv;
     }
     
     //! The total visible volume back calculated to an 'effective' sphere radius
     double calc_effective_radius() {
          double vv = 0.0;
          std::vector<double> wv = calc_visible_volume();
          for (int w=0; w<n_windows; w++)
               vv += wv[w];
          double r = pow(vv*3.0/4.0/M_PI, 1.0/3.0);
          return r;
     }

     //! Calculate a list of neighbours in contact. A contact is defined as a shading disk that 
     //! covers at least one angle points in the current state of the object. 
     //! A contact is returned as a [res_index, atom_index, covered points] vector
     //! Note that the same residue,atom contact can show up in different windows. These contributions should be added.
     std::vector<std::vector<std::vector<int> > > get_contacts() {
          int max_res_index=0, min_res_index=0;
          int max_atom_index=0, min_atom_index=0;
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    int ri = angle_point[t][p].get_res_index();
                    int ai = angle_point[t][p].get_atom_index();
                    if (ri > max_res_index)
                         max_res_index = ri;
                    if (ri < min_res_index)
                         min_res_index = ri;
                    if (ai > max_atom_index)
                         max_atom_index = ai;
                    if (ai < min_atom_index)
                         min_atom_index = ai;
               }
          }

          // check for no shades
          if (max_res_index == min_res_index) {
               // zero contacts in all windows
               return std::vector<std::vector<std::vector<int> > >(this->n_windows, std::vector<std::vector<int> >(0));
          }
               
          // check for negative indeices
          int res_offset = 0;
          if (min_res_index < 0) {
               res_offset = -min_res_index;
               max_res_index -= min_res_index;
          }
          int atom_offset = 0;
          if (min_atom_index < 0) {
               atom_offset = -min_atom_index;
               max_atom_index -= min_atom_index;
          }

          // // 3D vector array counter
          // int ***shade_counts;
          const int res_size = max_res_index+1;
          const int atom_size = max_atom_index+1;
          // shade_counts = new int**[this->n_windows];
          // for (int i=0; i<(this->n_windows); i++) {
          //      shade_counts[i] = new int*[max_res_index+1];
          //      for (int j=0; j<(max_res_index+1); j++) {
          //           shade_counts[i][j] = new int[max_atom_index+1];
          //           for (int k=0; k<(max_atom_index+1); k++)
          //                shade_counts[i][j][k] = 0;
          //      }
          // }

          std::vector<std::vector<std::vector<short int> > > shade_counts(
               this->n_windows, std::vector<std::vector<short int> >(
                    res_size, std::vector<short int>(
                         atom_size, 0)));
          
          // loop over angle points
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    // check for no shade
                    if (! this->angle_point[t][p].is_covered())
                         continue;
                    int res_index = this->angle_point[t][p].get_res_index() + res_offset;
                    int atom_index = this->angle_point[t][p].get_atom_index() + atom_offset;
                    int w = this->angle_point[t][p].window;
                    shade_counts[w][res_index][atom_index] += 1;
               }
          }

          // put contacts in list
          std::vector<std::vector<std::vector<int> > > contacts(0);
          
          // for all windows
          for (int w=0; w<this->n_windows; w++) {
               std::vector<std::vector<int> > win_contacts(0);
               // find all contacts
               for (int r=0; r<res_size; r++) {
                    for (int a=0; a<atom_size; a++) {
                         // and put them in a list
                         if (shade_counts[w][r][a] > 0) {
                              std::vector<int> single_contact(3);
                              single_contact[0] = r-res_offset;
                              single_contact[1] = a-atom_offset;
                              single_contact[2] = (int) shade_counts[w][r][a];
                              win_contacts.push_back(single_contact);
                         }
                    }
               }
               contacts.push_back(win_contacts);
          }

          // // release memory
          // for (int i=0; i<(this->n_windows); i++) {
          //      for (int j=0; j<(max_res_index+1); j++) {
          //           delete [] shade_counts[i][j];
          //      }
          //      delete [] shade_counts[i];
          // }
          // delete [] shade_counts;

          return contacts;
     }

     //! Calculate a histogram for each window of how many angle points is covered by each residue
     //! type. After the first 'n_types' elements an extra element gives the number of exposed
     //! (not covered) angle points. Types with index higher than n_types are ignored.
     //! /param n_types Highest type index included in the histogram
     std::vector<std::vector<int> > residue_type_histogram(int n_types) {
          assert(n_types >= 0);

          // return vector
          std::vector<std::vector<int> > counts =
               std::vector<std::vector<int> >(n_windows, std::vector<int>(n_types+1,0));

          // loop over sphere points
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    int w = this->angle_point[t][p].window;
                    // check for exposure
                    if (! this->angle_point[t][p].is_covered()) {
                         counts[w][n_types] += 1;
                    } else {
                         int res_type = this->angle_point[t][p].get_res_type();
                         assert(res_type >= 0);
                         if (res_type<n_types)
                              counts[w][res_type] += 1;
                    }
               }
          }
          return counts;
     }

     //! Calculate a histogram for each window of how many angle points is covered by each atom
     //! type. After the first 'n_types' elements an extra element gives the number of exposed
     //! (not covered) angle points. Types with index higher than n_types are ignored.
     //! /param n_types Highest type index included in the histogram
     std::vector<std::vector<int> > atom_type_histogram(int n_types) {          
          assert(n_types >= 0);

          // return vector
          std::vector<std::vector<int> > counts =
               std::vector<std::vector<int> >(n_windows, std::vector<int>(n_types+1,0));

          // loop over sphere points
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    int w = this->angle_point[t][p].window;
                    // check for exposure
                    if (! this->angle_point[t][p].is_covered()) {
                         counts[w][n_types] += 1;
                    } else {
                         int atom_type = this->angle_point[t][p].get_atom_type();
                         assert(atom_type >= 0);
                         if (atom_type<n_types)
                              counts[w][atom_type] += 1;
                    }
               }
          }
          return counts;
     }
     
     //! Remove all shades from the object
     void reset() {
          this->shade_counter = 0;
          this->n_windows = 1;
          int t_size = this->theta_list.size();
          for (int t=0; t<t_size; t++) {
               int p_size = phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    this->angle_point[t][p].reset();
               }
          }
     }

     //! Output operator - dump a list with status of each point in the current state of the object.
     friend std::ostream &operator<<(std::ostream &o, Sphere &sphere) {
          // const char *aa_name[22] = {"ala", "cys", "asp", "glu", "phe", "gly", "his", "ile", "lys", "leu",
          //                            "met", "asn", "pro", "gln", "arg", "ser", "thr", "val", "trp", "tyr"};
          
          // char buffer 
          char buffer[512]; //ugly but I'm oldschool
          sprintf(buffer,"Status after %d shades. Sphere has %d points in %d windows\n",
                  sphere.shade_counter, sphere.n_points, sphere.n_windows);
          o << buffer;
          
          int counter = 1;
          std::vector<int> w_counter(sphere.n_windows, 0);
          std::vector<double> ba_list(0);
          double to_deg = 180.0/M_PI;
          int t_size = sphere.theta_list.size();
          for (int t=0; t<t_size; t++) {
               double st = t*sphere.d_theta;
               sprintf(buffer, "theta = %5.1f covers from %5.1f to %5.1f\n", sphere.theta_list[t]*to_deg, st*to_deg,
                       (st+sphere.d_theta)*to_deg);
               o << buffer;
               int p_size = sphere.phi_list_list[t].size();
               for (int p=0; p<p_size; p++) {
                    // printf("theta = %5.1f covers from %5.1f to %5.1f\n", sphere.theta_list[t]*to_deg, st*to_deg,
                    //        (st+sphere.d_theta)*to_deg);
                    // pyramid base area for angle point
                    double dist = sphere.angle_point[t][p].get_distance();
                    double theta_base = 2.0*dist*tan(sphere.d_theta/2.0);
                    double phi_base = 2.0*dist*tan(sphere.d_phi_list[t]/2.0)*sin(sphere.theta_list[t]);
                    double base_area = theta_base*phi_base;
                    ba_list.push_back(base_area);
                    int w = sphere.angle_point[t][p].window;
                    sprintf(buffer, "%4d:  t=%5.1f  p=%5.1f  r=%6.3f, pyramid base area: %5.2f, "
                            "shadow from res %3d(type %2d) atom %2d (type %2d), through window %d\n", 
                            counter, sphere.theta_list[t]*to_deg, sphere.phi_list_list[t][p]*to_deg, dist, base_area,
                            sphere.angle_point[t][p].get_res_index(), sphere.angle_point[t][p].get_res_type(), 
                            sphere.angle_point[t][p].get_atom_index(), sphere.angle_point[t][p].get_atom_type(), w);
                    o << buffer;
                    counter += 1;
                    w_counter[w] += 1;
               }
          }
          double surface = 4*M_PI*sphere.radius*sphere.radius;
          double sum_ba = 0.0, std_ba = 0.0, max_ba = 0.0, min_ba = 10000.0;
          int ba_length = ba_list.size();
          for (int i=0; i<ba_length; i++) {
               sum_ba += ba_list[i];
               if (ba_list[i] > max_ba)
                    max_ba = ba_list[i];
               if (ba_list[i] < min_ba)
                    min_ba = ba_list[i];
          }
          double mean_ba = sum_ba/ba_length;
          for (int i=0; i<ba_length; i++) {
               double diff = ba_list[i]-mean_ba;
               std_ba += diff*diff;
          }
          std_ba /= (ba_length-1);
          sprintf(buffer, "Analytical sphere surface area: %.2f.  Sum of base area: %.2f\n", surface, sum_ba);
          o << buffer;
          sprintf(buffer, "Base area: %.2f +- %.2e,   max = %.2f   min = %.2f\n", mean_ba, std_ba, max_ba, min_ba);
          o << buffer;
          
          o << "Window points: ";
          int w_sum = 0;
          for (int w=0; w<sphere.n_windows; w++) {
               sprintf(buffer, "%2d ", w_counter[w]);
               o << buffer;
               if (w < (sphere.n_windows-1))
                    o << "+ ";
               w_sum += w_counter[w];
          }
          sprintf(buffer, "= %2d\n", w_sum);
          o << buffer;
          return o;
     }

private:
     
     //! Length of vector
     inline double norm(std::vector<double> &v) {
          assert(v.size() == 3);
          return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
     }
     
     //! Normalize vector
     inline void normalize(std::vector<double> &v) {
          assert(v.size() == 3);
          double length = norm(v);
          for (int i=0; i<3; i++)
               v[i] /= length;
     }
     
     //! Vector dot product
     inline double dot(std::vector<double> &v1, std::vector<double> &v2) {
          assert(v1.size() == 3);
          assert(v2.size() == 3);
          return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
     }
     
     //! Vector cross product
     inline std::vector<double> cross(std::vector<double> &v1, std::vector<double> &v2) {
          assert(v1.size() == 3);
          assert(v2.size() == 3);
          std::vector<double> ret(3);
          ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
          ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
          ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
          return ret;
     }
     
     //! Transformes (r, theta, phi) in radians to vector in cartesian corrdinates (x,y,z)
     inline std::vector<double> sphere_to_cart(double radius, double theta, double phi) {
          std::vector<double> ret(3);
          ret[0] = radius*sin(theta)*cos(phi);
          ret[1] = radius*sin(theta)*sin(phi);
          ret[2] = radius*cos(theta);
          return ret;
     }
     
     //! Transformes Vector with (x,y,z) to (r,theta,phi) in radians
     inline std::vector<double> cart_to_sphere(double x, double y, double z) {
          std::vector<double> ret(3);
          double r = sqrt(x*x + y*y + z*z);
          ret[0] = r;
          // theta defined from positive z-axis
          ret[1] = acos(z/r);
          // phi defined from positive x-axis towards positive y-axis in [0,2pi[
          ret[2] = atan2(y,x);
          if (ret[2] < 0.0)
               ret[2] += 2*M_PI;
          return ret;
     }
};

#endif
