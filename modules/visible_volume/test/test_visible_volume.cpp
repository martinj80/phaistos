
#include <cmath>
#include <iostream>
#include <cstdio>

#include "models/visible_volume/sphere.h"
#include "energy/visible_volume.h"

#include "energy/energy.h"

#include "protein/chain_fb.h"

#include "utils/vector_matrix_3d.h"

using namespace std;
using namespace phaistos;
using namespace phaistos::definitions;

int main (int argc, char *argv[]) {

     if (argc < 2) {
          printf("\nusage: ./test_visible_volume <PDB file> <res index, opt>\n");
          return 1;
     }

     // Sphere settings
     double sphere_radius=11.0;
     int min_sphere_points=1000;
     double atom_radius = 1.8;
     
     // Construct chain
     ChainFB *chain = new ChainFB(argv[1], ALL_ATOMS);
     // Add atoms missing in the pdb structure
     chain->add_atoms(ALL_PHYSICAL_ATOMS);
     int chain_length = chain->size();
     
     // Parse argument - find residue array index
     int index_begin=-1, index_end=-1;
     if (argc == 3) {
          int seq_index = atoi(argv[2]);
          for (int i=0; i<chain_length; i++) {
               if ((*chain)[i].index_res_seq == seq_index) {
                    index_begin = i;
                    break;
               }
          }
          if (index_begin < 0) {
               printf("\nCould not find index \"%s\" - quitting\n",argv[2]);
               return 1;
          }
          index_end = index_begin+1;
     } else {
          index_begin = 0;
          index_end = chain_length;
     }

     // energy collection
     Energy<ChainFB> *energy;
     energy = new Energy<ChainFB>(chain);

     // energy term
     TermVisibleVolume<ChainFB>::Settings settings;
     settings.angexp = true;
     settings.min_sphere_points = min_sphere_points;
     settings.sphere_radius = sphere_radius;
     settings.atom_radius = atom_radius;
     TermVisibleVolume<ChainFB> *energy_term;
     energy_term = new TermVisibleVolume<ChainFB> (chain, settings);
     energy->add_term(energy_term);

     // for ( int i=0; i<999; i++)
     //      energy->evaluate();
     
     cout << energy->evaluate() << endl;
     
     // Energy<ChainFB> e2 = *energy;
     // delete energy;
     // cout << e2.evaluate() << endl;

     // // construct sphere object
     // Sphere *sphere;
     // sphere = new Sphere(sphere_radius, min_sphere_points);

     // // test copy constructor
     // Sphere s2 = *sphere;
     // delete sphere;
     // cout << s2.angle_point[0].size() << endl;
     
     // // calc visible volume for choosen indices
     // for (int i=index_begin; i<index_end; i++) {
     //      ResidueFB *resi = &(*chain)[i];
     //      int resi_type = resi->residue_type;
          
     //      // find contact center
     //      Vector_3D cc;
     //      if (resi_type == definitions::ALA)
     //           cc = TermVisibleVolume<ChainFB>::get_cb_pos(resi);
     //      else if (resi_type == definitions::GLY)
     //           cc = TermVisibleVolume<ChainFB>::calc_cb_pos(resi);
     //      else
     //           cc = TermVisibleVolume<ChainFB>::calc_geo_center(resi);
          
     //      // add neighbours
     //      int chain_size = chain->size();
     //      for (int j=0; j<chain_size; j++) {
               
     //           // no shadow from own side chain
     //           if (i==j)
     //                continue;
               
     //           ResidueFB *resj = &(*chain)[j];
     //           int resj_type = (*chain)[j].residue_type;
               
     //           int atom_size = resj->atoms.size();
     //           for (int a=0; a<atom_size; a++) {
     //                Atom *atom = resj->atoms[a];
     //                int atom_type = atom->atom_type;
                    
     //                // atom position relative to contact center
     //                Vector_3D pos = atom->position - cc;
                    
     //                // put in sphere object
     //                sphere.shade(pos.get_array(), atom_radius, j, a, resj_type, atom_type);
     //           }

     //           // calculate C alpha - contact center vector
     //           Vector_3D ca_cc_vec3d = cc - (*resi)[CA]->position;
     //           double *ca_cc_array = ca_cc_vec3d.get_array();
     //           std::vector<double> ca_cc_vec(ca_cc_array, ca_cc_array+3);
               
     //           // calculate C alpha - C vector
     //           Vector_3D ca_c_vec3d = (*resi)[C]->position - (*resi)[CA]->position;
     //           double *ca_c_array = ca_c_vec3d.get_array();
     //           std::vector<double> ca_c_vec(ca_c_array, ca_c_array+3);
               
     //           // Orient sphere
     //           sphere.assign_tetrahedral_window(ca_cc_vec, ca_c_vec);
               
     //           // Calculate visible volume seen through each window
     //           vector<double> wv = sphere.calc_visible_volume();
               
     //           // Calculate contacts multinomial
     //           vector<vector<int> > contacts = sphere.residue_type_contacts(20);
               
     //           // Dump sphere info for debug
     //           // cout << sphere << endl;
               
     //           printf("%s-%s%04d:  %8.3f  %8.3f  %8.3f  %8.3f\n",
     //                  "xxxxx", residue_name[resi->residue_type], resi->index_res_seq, wv[0], wv[1], wv[2], wv[3]);
     //           int c_size = contacts[0].size();
     //           for (int w=0; w<sphere.n_windows; w++) {
     //                printf("%1d: ",w);
     //                int c_sum = 0;
     //                for (int c=0; c<c_size; c++) {
     //                     printf("%3d ",contacts[w][c]);
     //                     c_sum += contacts[w][c];
     //                }
     //           printf("= %3d\n",c_sum);
     //           }
     //           int c_sum = 0;
     //           printf("T: ");
     //           for (int c=0; c<c_size; c++) {
     //                int w_sum = 0;
     //                for (int w=0; w< sphere.n_windows; w++)
     //                     w_sum += contacts[w][c];
     //                printf("%3d ",w_sum);
     //                c_sum += w_sum;
     //           }
     //           printf("= %3d\n\n",c_sum);
               
     //           // Reset sphere 
     //           sphere.reset();
     //      }
          
     // }
     
     return 0;
};
