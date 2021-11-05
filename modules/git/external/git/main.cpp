// main.cpp --- Simple example using the Gauss Integrals Tuned (GIT) library code
//              Original C-code by Peter Roegen, Copyright (C) 2011
//              Please cite: 
//              P. R{\o}gen, Evaluating protein structure descriptors and tuning Gauss integral
//              based descriptors: Journal of Physics Condensed Matter, vol: 17, pages: 1523-1538, 2005.
//              
//              C++ version by Tim Harder and Wouter Boomsma, Copyright (C) 2011
//
// This file is part of Git
//
// Git is free software: you can redistribute it and/or modify
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
// along with Git.  If not, see <http://www.gnu.org/licenses/>.
//

#include "git.h"

#include <dirent.h>
#include <cstring>

using namespace git;

//! Output git vector
//!
//! \param git_vector Gauss Integral vector
//! \param n_res Number of residues
//! \param smoothen_backbone_mode Whether to run in smoothen_backbone mode
//! \param output_file File pointer to which output is written
void print_git_vector(const std::vector<double> &git_vector, int n_res, bool smoothen_backbone_mode, FILE *output_file) {


     unsigned int active_vector_size = git_vector.size();
     if (smoothen_backbone_mode) {
          fprintf(output_file ,"%i %5.4e",n_res, 19.11*pow(n_res,1.0/3.0));
     } else {
          fprintf(output_file, "%i %5.4e", n_res, git_vector.back());
          active_vector_size--;
     }
     for (unsigned int i=0; i<active_vector_size; ++i) {
          fprintf(output_file," %5.4e", git_vector[i]);
     }
     fprintf(output_file,"\n");     
}

//! Parse atom line
//!
//! \param fi File pointer
//! \param line Data container to write to
int get_atom(FILE *fi, char *line) {
     char *cp;
     while (cp=fgets(line,999,fi),strstr(line,"ATOM")!=line && cp!=NULL);
     if (cp==NULL) return 0;
     return 1;
}


//! Parse input from PDB file and calculate Gaus integrals
//!
//! \param git Git object used to calculate Gauss integrals
//! \param dir Directory name
//! \param filename Filename in directory
//! \param output_file File pointer to write output to
//! \param error_file File pointer to write errors to
//! \param dipeptide_mode Whether to run in dipeptide mode
//! \return Number of chains parsed
int parse_pdb_calc_gauss_integrals(Git &git, std::string dir, std::string filename, FILE *output_file, FILE *error_file, bool dipeptide_mode) {

     FILE *fi;
     double x, y, z;
     int nr, no=-1500, link_count=0, start_chain,i;
     char chain_name, new_chain_name;
     char ThreeLeterAmno[3];
     std::string dir_filename;
     char *cp;
     std::vector<Vec3> polyl;
     char line[1000];
     int n_proteins = 0;

     dir_filename = dir;
     dir_filename += filename;

     fi=fopen(dir_filename.c_str(), "r");
     chain_name='a';

     while (true) {
          if (get_atom(fi, line)==0) { 
               if(link_count<git.get_max_length() && link_count>9) {
                    if(dipeptide_mode==0){
                         fprintf(output_file,"%s %c %i ", filename.c_str(), chain_name, no-start_chain+1-link_count);
                    }
                    if(dipeptide_mode==1){
                         fprintf(output_file,"%s ",filename.c_str());
                    }
                    printf("  number of residues in chain %c is %i\n",chain_name,link_count);
                    std::vector<double> git_vector = git.generate_gauss_integrals(polyl);
                    print_git_vector(git_vector, polyl.size(), git.get_smoothen_backbone_mode(), output_file);
                    fclose(fi);
                    n_proteins++;
                    return n_proteins;
               }
               fclose(fi);  
               return n_proteins;
          }
          sscanf(line+22,"%4i",&nr);
          sscanf(line+21,"%1c",&new_chain_name);
          if(new_chain_name==' ') new_chain_name='0';
  
          if (nr!=no && line[13]=='C' && line[14]=='A') {
               if(new_chain_name!=chain_name && chain_name!='a') {
                    if(link_count<git.get_max_length() && link_count>9) {
                         /* printing pdb entry  fprintf(output_file,"%c%c%c%c:%c  %i %i",line[72],line[73],line[74],line[75],chain_name,link_count,no-start_chain+1-link_count);            */
                         if(dipeptide_mode==0){
                              fprintf(output_file,"%s %c %i ", filename.c_str(), chain_name,no-start_chain+1-link_count);
                         }
                         if(dipeptide_mode==1){
                              fprintf(output_file,"%s ",filename.c_str());
                         }
	     
                         printf("  number of residues in chain %c is %i\n",chain_name,link_count);
                         std::vector<double> git_vector = git.generate_gauss_integrals(polyl);
                         print_git_vector(git_vector, polyl.size(), git.get_smoothen_backbone_mode(), output_file);
                         n_proteins++;
                    }  
                    link_count=0;
               }
               link_count++;
               if (link_count>git.get_max_length()-1) {
                    printf(" >MAXLENGTH=%i residues\n",git.get_max_length());
                    fprintf(error_file, "WARNING!  Chain %c of <%s> is longer than MAXLENGTH=%i. \n          No further attempt to calculate Gauss Integrals for chains of\n          <%s> has been made.\n", chain_name, filename.c_str(), git.get_max_length(), filename.c_str());fclose(fi); return n_proteins;}
               no=nr;
               if(link_count==1) start_chain=no;
               sscanf(line+30,"%lf %lf %lf",&x,&y,&z);
               ThreeLeterAmno[0]=line[17];
               ThreeLeterAmno[1]=line[18];
               ThreeLeterAmno[2]=line[19];
               // SequinceInteger[link_count]=Three2IntegerProteinId(ThreeLeterAmno);
               polyl.push_back(Vec3(x,y,z));
               // for (unsigned int i=0; i<polyl.size(); ++i) {
               //      std::cout << polyl[i] << " ";
               // }
               // std::cout << "\n";
               // polyl[link_count][0]=x;
               // polyl[link_count][1]=y;
               // polyl[link_count][2]=z;
               chain_name=new_chain_name;
      
          }
     }
     fclose(fi);
     return n_proteins;
}


//! Iterate over directory
//!
//! \param git Git object used to calculate Gauss integrals
//! \param dir Directory name
//! \param output_file File pointer to write output to
//! \param error_file File pointer to write errors to
//! \param dipeptide_mode Whether to run in dipeptide mode
//! \return Number of chains parsed
int dir_walk(Git &git, std::string dir, FILE *output_file, FILE *error_file, bool dipeptide_mode=0) {

     char name[1024];
     struct dirent *dp;
     DIR *dfd;

     int n_proteins = 0;

     if ((dfd=opendir(dir.c_str()))==NULL) {
          printf("Sorry - could not open directory <%s>\n",dir.c_str());
          return n_proteins; 
     }

     while ((dp=readdir(dfd))!=NULL) {
          if (strcmp(dp->d_name,".")==0 || strcmp(dp->d_name,"..")==0) 
               continue;
          sprintf(name,"%s/%s",dir.c_str(),dp->d_name);
          if (dp->d_name[strlen(dp->d_name)-4]=='.' && 
              dp->d_name[strlen(dp->d_name)-3]=='p' && 
              dp->d_name[strlen(dp->d_name)-2]=='d' && 
              dp->d_name[strlen(dp->d_name)-1]=='b') {

               printf("%i It's %s\n",n_proteins,dp->d_name);
               
               n_proteins += parse_pdb_calc_gauss_integrals(git, dir ,dp->d_name, output_file, error_file, dipeptide_mode);
          }
     }
     closedir(dfd);

     return n_proteins;
}



//! Main function
int main(int argc, char *argv[]) {

     // Parse command line input
     if (argc<2) {
          printf("Please specify input directory (with trailing \"/\") followed by \n Averge_gauss_integral_file (optional),\n output file (optional) and error file name (optional).\n");
          return 1;
     }

     bool smoothen_backbone_mode = 1;
     std::string directory_name = argv[1];
     std::string average_gauss_integral_file = "";
     FILE *output_file = stdout;
     FILE *error_file = stderr;

     if (argc>2)
          average_gauss_integral_file = argv[2];

     if (argc>3)
          output_file = fopen(argv[3], "w");

     if (argc>4)
          error_file = fopen(argv[4], "w");

     // Create git object
     Git git(average_gauss_integral_file, smoothen_backbone_mode);

     // Calculate GIT vectors for the given input directory
     int n_proteins = dir_walk(git, directory_name, output_file, error_file);

     
     printf("Calculated 29 Tuned Gauss Integrals for %i protein strands.\n", n_proteins);

     // Cleanup
     if (output_file != stdout)
          fclose(output_file);
     if (error_file != stderr)
          fclose(error_file);
          
     return 0;
}
