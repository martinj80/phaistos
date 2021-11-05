// saxs.h --- Small Angle X-ray scattering energy function
// Copyright (C) 2008-2009 Christian Andreetta, Kasper Stovgaard, Wouter Boomsma
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

// SAXS energy function: Debye O(N**2) implementation

#ifndef SAXS_H
#define SAXS_H

#include <fstream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "energy/energy_term.h"


namespace phaistos {

//! Numerical constants
namespace saxs_definitions {
const static double two_pi = 2*M_PI;
const static double sqrt_two_pi = sqrt(two_pi);
const static int sine_lut_size = 4096;
const static int sine_lut_bitmask = 4095;
}

//! Probability density function for a Normal distribution
//!
//! \param x value to compute pdf for
//! \param mu mean of the Normal distribution
//! \param sigma standard deviation of the Normal distribution
//! \return probability density value
double normal_pdf(double x, double mu=0.0, double sigma=1.0) {
     // gaussian pdf= 1/(sigma*sqrt(2*PI)) * exp( -(x-mu)**2/(2*sigma**2) )
     return 1/( sigma*saxs_definitions::sqrt_two_pi ) * std::exp( -1 * pow(x-mu,2)/(2*pow(sigma,2)) );
}

//! Probability density function for a Normal distribution (log version)
//!
//! \param x value to compute pdf for
//! \param mu mean of the Normal distribution
//! \param sigma standard deviation of the Normal distribution
//! \return log(probability density value)
double ln_normal_pdf(double x, double mu=0.0, double sigma=1.0) {
     // gaussian pdf= 1/(sigma*sqrt(2*PI)) * exp( -(x-mu)**2/(2*sigma**2) )
     return -std::log(sigma*saxs_definitions::sqrt_two_pi) -1*(x-mu)*(x-mu)/(2*sigma*sigma);
}

// order follows residueEnum
/*const double PS_WEIGHT[20] = {
     0.0,32.065,44.0095,56.0202,72.0642,0.0,64.0455,36.0321,
     50.0388,36.0321,56.0864,42.0168,24.0214,54.0275,78.0522,  
     15.9994,28.0101,24.0214,110.0923,88.0636
     };
*/

//! Returns the char* "true" or "false" according to a boolean value (helper function)
//!
//! \param b boolean value
//! \return "true" or "false" according to the input boolean value
inline const char * const bool_to_str(bool b) {
     // returns char* describing the value of a boolean
     return b ? "true" : "false";
}

//! Debye SAXS energy term
//!
//! Pairwise Debye computation
template <typename CHAIN_TYPE>
class TermSaxsDebye: public EnergyTermCommon<TermSaxsDebye<CHAIN_TYPE>, CHAIN_TYPE> {
private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermSaxsDebye<CHAIN_TYPE>, CHAIN_TYPE> EnergyTermCommon;         

     //! Form factor enum - two body
     enum TWO_BODIES_FF_NAMES {
          BACKBONE_NORM=0, // normal backbone *with* C_beta
          ALA_COMPLETE, ARG_SC, ASN_SC,       ASP_SC, CYS_SC,  // GLY and ALA don't have a side chain
          GLN_SC,       GLU_SC, GLY_COMPLETE, HIS_SC, ILE_SC,
          LEU_SC,       LYS_SC, MET_SC,       PHE_SC, PRO_SC,
          SER_SC,       THR_SC, TRP_SC,       TYR_SC, VAL_SC
          };

     //! Form factors enum - one body
     enum ONE_BODY_FF_NAMES {
          ALA_C=0, ARG_C, ASN_C, ASP_C, CYS_C, GLN_C, GLU_C, GLY_C, HIS_C, ILE_C,
          LEU_C, LYS_C, MET_C, PHE_C, PRO_C, SER_C, THR_C, TRP_C, TYR_C, VAL_C
          };

     //! Form factors database
     std::vector< std::vector<double> > form_factors_db; // size: Q_NUM * FF_NUM

     //@{
     //! saxs profile
     std::vector<double> q_vec;                          // size: Q_NUM
     std::vector<double> saxs_profile;                   // size: Q_NUM
     std::vector<double> saxs_profile_accepted;          // size: Q_NUM
     std::vector<double> saxs_profile_reference;         // size: Q_NUM
     std::vector<double> saxs_profile_reference_errors;  // size: Q_NUM
     //@}

     //! form factors containers for SAXS computation - two body
     std::vector<TWO_BODIES_FF_NAMES> ff_ids;

     //! form factors containers for SAXS computation - one body
     std::vector<ONE_BODY_FF_NAMES> ff_one_body_ids;

     //! Whether a dummy center changed location
     std::vector<bool> ff_centers_changed;

     //! Vector of dummy centers     
     std::vector<Vector_3D> ff_centers;

     //! Center distances
     std::vector< std::vector<double> > ff_centers_dist;

     //! sinc(x) matrix
     std::vector< std::vector< std::vector<double> > > ff_sinc;

     //! cache: ACCEPTED container
     std::vector<Vector_3D> ff_centers_accepted;

     //! cache: PROPOSAL container
     std::vector< std::vector<double> > ff_centers_dist_proposal;

     //! sinc(x) matrix proposal
     std::vector< std::vector< std::vector<double> > > ff_sinc_proposal;

     //! Table of sin(x) values
     std::vector<double> sine_table;

     //! Number of dymmy bodies
     unsigned int dummy_bodies_num;

     //! Whether cache is outdated
     bool cache_is_outdated;

     //! Current energy value (class variables)
     double current_energy;

public:

     //! Local settings class
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! Path to file containing SAXS intensities
          std::string saxs_intensities_filename;

          //! Path to file containing SAXS form factors
          std::string saxs_form_factors_filename;

          //! Alpha parameter for experimental error term
          double exp_errors_alpha;

          //! Beta parameter for experimental error term
          double exp_errors_beta;

          //! Number of q-bins used in the curve calculation
          unsigned int q_bins;

          //! Energy evaluation starts at this q-bin (default: 0). the first bin is used for normalization
          unsigned int q_bins_first;

          //! Use one-body instead of two-body form factor model
          bool one_body_model;

          //! Use sin() lookup table instead of the full trigonometric evaluation
          bool sine_lookup_table;

          //! Constructor
          Settings(std::string saxs_intensities_filename="",
                   std::string saxs_form_factors_filename="",
                   double exp_errors_alpha=0.15,
                   double exp_errors_beta=0.30,
                   unsigned int q_bins=51,
                   unsigned int q_bins_first=0,
                   bool one_body_model=false,
                   bool sine_lookup_table=true)
               : saxs_intensities_filename(saxs_intensities_filename),
                 saxs_form_factors_filename(saxs_form_factors_filename),
                 exp_errors_alpha(exp_errors_alpha),
                 exp_errors_beta(exp_errors_beta),
                 q_bins(q_bins),
                 q_bins_first(q_bins_first),
                 one_body_model(one_body_model),
                 sine_lookup_table(sine_lookup_table) {
          }

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "saxs-intensities-filename:" << settings.saxs_intensities_filename << "\n";
               o << "saxs-form-factors-filename:" << settings.saxs_form_factors_filename << "\n";
               o << "exp-errors-alpha:" << settings.exp_errors_alpha << "\n";
               o << "exp-errors-beta:" << settings.exp_errors_beta << "\n";
               o << "q-bins:" << settings.q_bins << "\n";
               o << "q-bins-first:" << settings.q_bins_first << "\n";
               o << "one-body-model:" << settings.one_body_model << "\n";
               o << "use-sine-lookup-table:" << settings.sine_lookup_table << "\n";
               o << static_cast<typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }
     } settings;    //!< Local settings object 


     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermSaxsDebye(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                   RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain,"saxs-debye", settings, random_number_engine),
            settings(settings) {

          this->current_energy=0;

          if (settings.saxs_intensities_filename!="") { 
               initialize(); 
          }
     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermSaxsDebye(const TermSaxsDebye &other, 
                   RandomNumberEngine *random_number_engine,
                   int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),

            form_factors_db(other.form_factors_db),

            q_vec(other.q_vec),
            saxs_profile(other.saxs_profile),
            saxs_profile_accepted(other.saxs_profile_accepted),
            saxs_profile_reference(other.saxs_profile_reference),
            
            saxs_profile_reference_errors(other.saxs_profile_reference_errors),
            
            ff_ids(other.ff_ids),
            ff_one_body_ids(other.ff_one_body_ids),
            ff_centers_changed(other.ff_centers_changed),
            ff_centers(other.ff_centers),
            ff_centers_dist(other.ff_centers_dist),
            ff_sinc(other.ff_sinc),
            ff_centers_accepted(other.ff_centers_accepted),
            ff_centers_dist_proposal(other.ff_centers_dist_proposal),
            ff_sinc_proposal(other.ff_sinc_proposal),
            sine_table(other.sine_table),

            dummy_bodies_num(other.dummy_bodies_num),
            cache_is_outdated(other.cache_is_outdated),

            current_energy(other.current_energy),
            settings(other.settings) {
     }


     //! move acceptance function: updates cache
     void accept() {
          // std::vector: operator= is a copy assignment
          if (settings.debug>=100) { printf("# SAXS: MOVE: accepted\n"); }
          if (cache_is_outdated) {
               saxs_profile_accepted=saxs_profile;
               ff_centers_accepted=ff_centers;
               // updates centers positions
               for (unsigned int j=0; j<dummy_bodies_num; j++) {
                    if (ff_centers_changed[j]) {
                         // upper part of the lower triangle, row-wise
                         for (unsigned int k=0; k<j; k++) {
                              ff_centers_dist[j][k]=ff_centers_dist_proposal[j][k];
                         }
                         // lower part of the lower triangle, column-wise
                         for (unsigned int r=j+1; r<dummy_bodies_num; r++) {
                              ff_centers_dist[r][j]=ff_centers_dist_proposal[r][j];
                         }
                    }
               }
               // updates sinc() matrices
               for (unsigned int qj=0; qj<settings.q_bins; qj++) {
                    for (unsigned int ffj=0; ffj<dummy_bodies_num; ffj++) {
                         if (ff_centers_changed[ffj]) {
                              // upper part of the lower triangle, row-wise
                              for (unsigned int ffk=0; ffk<=ffj; ffk++) {
                                   ff_sinc[qj][ffj][ffk] = ff_sinc_proposal[qj][ffj][ffk];
                              }
                              // lower part of the lower triangle, column-wise
                              for (unsigned int ffr=ffj+1; ffr<dummy_bodies_num; ffr++) {
                                   ff_sinc[qj][ffr][ffj] = ff_sinc_proposal[qj][ffr][ffj];
                              }
                         }
                    }
               }
               // cache is now updated
               cache_is_outdated=false;
          }
     }


     //! cache reject: restore proposal from cache
     void reject() {
          // std::vector: operator= is a copy assignment
          if (settings.debug>=100) { printf("# SAXS: MOVE: rejected\n"); }
          if (cache_is_outdated) {
               saxs_profile=saxs_profile_accepted;
               ff_centers=ff_centers_accepted;
               cache_is_outdated=false;
               if (settings.debug>=100) {
                    // sanity checks
                    assert( &(ff_centers[0]) != &(ff_centers_accepted[0]) );
                    assert(  saxs_profile == saxs_profile_accepted );
                    assert( &saxs_profile != &saxs_profile_accepted );
               }
          }
     }

     //! energy initialization
     //! - sin() lookup table definition
     //! - sets the number of form factors according to database size
     //! - allocation of the cache and database memory
     //! - calls: read_ff(), reading of form factors table
     //! - calls: read_reference(), reading of SAXS experimental reference
     //! - calls: generate_experimental_errors(), generation of the parameters of the Normal distribution around each data point
     //! - calls: saxs_debye_ff_ids(), computation of the bodies centers
     //! - calls: saxs_debye_vectors_initialize(), computation of cache vectors and matrices
     void initialize() {
          // sine look-up table
          if (settings.sine_lookup_table==true) {
               sine_table.resize(saxs_definitions::sine_lut_size);
               for (int i = 0; i < saxs_definitions::sine_lut_size; i++) {
                    sine_table[i] = sin( (saxs_definitions::two_pi*i)/saxs_definitions::sine_lut_size );
               }
          }
          // FF DB
          int number_of_ffs = settings.one_body_model ? 20 : 21;
          if (settings.debug >= 5) {
               printf("# DEBUG: SAXS: number of form factors used = %d\n",number_of_ffs);
          }
          form_factors_db.resize(settings.q_bins);
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               form_factors_db[qj].resize(number_of_ffs);
          }
          // DB space allocation
          q_vec.resize(settings.q_bins);
          saxs_profile.resize(settings.q_bins);
          saxs_profile_reference.resize(settings.q_bins);
          saxs_profile_reference_errors.resize(settings.q_bins);
          // FF initialization
          read_ff();              // ff vectors
          read_reference();       // crysol intensities, q vector
          generate_experimental_errors();  // "experimental error" definition
          // CHAIN initialization
          cache_is_outdated=false;
          saxs_debye_ff_ids();
          saxs_debye_vectors_initialize();
          // reference: total energy
          test_energy_reference();
     }

     //! reads the form factor database, assumed as a python list of lists
     void read_ff() {
          /*
               reads python generated form factors, stores in internal database
          */
          unsigned int ffj=0;
          int qj=-1;
          std::string elem_str="";
          double elem_val=0;
          std::ifstream f_h;

          f_h.open( (settings.saxs_form_factors_filename).c_str() );
          if (!f_h) {
               printf("#ERROR: saxs input: unable to open file '%s'\n", (settings.saxs_form_factors_filename).c_str() );
               exit(1); // terminate with error
          }

          if (this->settings.debug>=100) { printf("# SAXS: reading FF file\n"); }
          while(!f_h.eof() and qj<int(settings.q_bins)) {
               f_h >> elem_str;
               if (settings.debug>=100) { printf("ff: elem_str: '%s'",elem_str.c_str()); }
               if ( elem_str.size()<=1 ) { continue; }
               if ( elem_str.compare(0,1,"[")==0 ) {
                    elem_str.erase(0,1);
                    ffj=0; qj++;
               }
               // have we finished parsing the form factors datafile?
               if (qj==int(settings.q_bins)) { break; }
               // store form factor in database
               elem_val=strtod( elem_str.c_str(),NULL );
               this->form_factors_db[qj][ffj]=elem_val;

               if (settings.debug>=100) { printf("    ff: elem_val: %f (qj: %d, ffj: %d)\n",elem_val,qj,ffj); }
               ffj++;
          }
          if (this->settings.debug>=100) { printf("\n"); }
          f_h.close();
     }

     //! reads the SAXS experimental reference, and the scattering angle vector from data file
     void read_reference() {
          /*
          reads the crysol profile and the q vector from datafile
          */
          double q_elem=0, profile_elem=0;
          unsigned int rec_num=0;
          std::ifstream profile_file_h;

          profile_file_h.open( (settings.saxs_intensities_filename).c_str() );
          if (!profile_file_h) {
               printf("#ERROR: saxs input: unable to open file '%s'\n", (settings.saxs_intensities_filename).c_str() );
               exit(1); // terminate with error
          }

          if (this->settings.debug>=100) { printf("# SAXS: reading reference file\n"); }
          while(!profile_file_h.eof() and rec_num<settings.q_bins) {
               profile_file_h >> q_elem >> profile_elem;
               this->q_vec[rec_num]=q_elem;
               this->saxs_profile_reference[rec_num]=profile_elem;

               if (settings.debug>=100) { printf("%2d -> %f -> %f\n",rec_num,q_elem,profile_elem); }
               rec_num++;
          }
          profile_file_h.close();
     }

     //! Given SAXS reference and the parameters alpha and beta, generate the paramters of the Normal distribution centered on the data file
     //! - if alpha >=0: sigma[j] = reference[j] * ( q_vec[j] + alpha ) * beta
     //! - if alpha <0:  sigma[j] = reference[j] * beta
     void generate_experimental_errors() {
          if (this->settings.debug>=10) { printf("# DEBUG: alpha: %.2f, beta: %.2f\n",settings.exp_errors_alpha,settings.exp_errors_beta); }
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               this->saxs_profile_reference_errors[qj]=this->saxs_profile_reference[qj]*settings.exp_errors_beta;
               if (settings.exp_errors_alpha >= 0) {
                    this->saxs_profile_reference_errors[qj]*=(this->q_vec[qj]+settings.exp_errors_alpha);
               }
          }
     }

     //! computation of the energy term, normalizing the SAXS proposal
     double evaluate_do( std::vector<double> *saxs_vec ) {
          double energy = 0.0;
          double int_scaled = 0.0;
          double scale_factor = this->saxs_profile_reference[settings.q_bins_first]/(*saxs_vec)[settings.q_bins_first];

          for (unsigned int qj=settings.q_bins_first; qj<settings.q_bins; qj++) {
               int_scaled = scale_factor*(*saxs_vec)[qj];
               //energy+= log( normal_pdf( int_scaled,this->saxs_profile_reference[qj],this->saxs_profile_reference_errors[qj] ) );
               energy+= ln_normal_pdf( int_scaled,this->saxs_profile_reference[qj],this->saxs_profile_reference_errors[qj] );
               if (this->settings.debug>=100) {
                    printf("%d -> q: %.3f -> saxs: %.3f, saxs_scaled: %.3f, ref: %.3f, sigma: %.3f -> energy (partial): %.3f\n",
                         qj,  this->q_vec[qj],int_scaled,int_scaled,this->saxs_profile_reference[qj],
                         //this->saxs_profile_reference_errors[qj],log( normal_pdf( int_scaled,
                         this->saxs_profile_reference_errors[qj],ln_normal_pdf( int_scaled,
                         this->saxs_profile_reference[qj],this->saxs_profile_reference_errors[qj] )
                    );
               }
          }
          /*
               monte carlo method expects energies to proceed from +infinity to zero
          */
          energy=-energy;
          return energy;
     }

     //! prints the reference energy
     void test_energy_reference() {
          // test: prints the reference energy
          printf("# DEBUG: reference energy: %.3f\n", evaluate_do(&saxs_profile_reference));
          this->chain->output_as_pdb_file("reference.pdb");
     }

     //! caller to identification of the centers, and computation of pairwise distances
     void saxs_debye_prepare() {
          saxs_debye_ff_centers();
          saxs_debye_ff_pairwise_distances();
     }

     //! computation of SAXS profile from scattering angle vector, bodies centers and form factors
     //!
     //! writes on internal vector this->saxs_profile
     void saxs_debye_compute() {
          //! - reads from this->ff_centers and this->ff_values_vec
          //! - compute saxs profile
          //! - intensity=sum over qj,ffi,ffj of ffi*ffj*sinc( qj * distance(ffi,ffj) )
          //! - remove outdated contributions: from ACCEPTED matrices
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               // loops over bodies
               for (unsigned int ffj=0; ffj<dummy_bodies_num; ffj++) {
                    if (ff_centers_changed[ffj]) {
                         // upper part of the lower triangle, row-wise
                         for (unsigned int ffk=0; ffk<=ffj; ffk++) {
                              saxs_profile[qj] -= ff_sinc[qj][ffj][ffk];
                         }
                         // lower part of the lower triangle, column-wise
                         for (unsigned int ffr=ffj+1; ffr<dummy_bodies_num; ffr++) {
                              // avoid to process twice
                              if (ff_centers_changed[ffr]) {continue;}
                              saxs_profile[qj] -= ff_sinc[qj][ffr][ffj];
                         }
                    }
               }
          }
          // update sinc() matrices
          if (!settings.one_body_model) {
               saxs_debye_ff_sinc_matrices();
          }
          else {
               saxs_debye_one_body_ff_sinc_matrices();
          }
          // adds the updated cached contributions: from PROPOSAL matrices
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               // loops over bodies
               for (unsigned int ffj=0; ffj<dummy_bodies_num; ffj++) {
                    if (ff_centers_changed[ffj]) {
                         // upper part of the lower triangle, row-wise
                         for (unsigned int ffk=0; ffk<=ffj; ffk++) {
                              saxs_profile[qj] += ff_sinc_proposal[qj][ffj][ffk];
                         }
                         // lower part of the lower triangle, column-wise
                         for (unsigned int ffr=ffj+1; ffr<dummy_bodies_num; ffr++) {
                              // avoid to process twice
                              if (ff_centers_changed[ffr]) {continue;}
                              saxs_profile[qj] += ff_sinc_proposal[qj][ffr][ffj];
                         }
                    }
               }
          }
          cache_is_outdated=true;
     }

     //! - fills the form factors id vector: it is done only once (the residues order and number doesn't change)
     //! - matches the AA names with the form factors IDs
     void saxs_debye_ff_ids() {
          // Import protein definitions (such as residue names)
          using namespace definitions;

          CHAIN_TYPE *chain=this->chain;
          unsigned int chain_len=(*chain).size();
          // two-body model
          if (!(this->settings.one_body_model)) {
               if ( this->ff_ids.empty() ) {
                    for (unsigned int j=0; j<chain_len; j++) {
                         if (this->settings.debug>=100) { std::cout << "# residue: "<<j<<"-> "<< (*chain)[j].residue_type << "\n"; }
                         if      ( (*chain)[j].residue_type == ALA ) { this->ff_ids.push_back(ALA_COMPLETE); }
                         else if ( (*chain)[j].residue_type == ARG ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(ARG_SC); }
                         else if ( (*chain)[j].residue_type == ASN ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(ASN_SC); }
                         else if ( (*chain)[j].residue_type == ASP ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(ASP_SC); }
                         else if ( (*chain)[j].residue_type == CYS ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(CYS_SC); }
                         else if ( (*chain)[j].residue_type == GLN ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(GLN_SC); }
                         else if ( (*chain)[j].residue_type == GLU ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(GLU_SC); }
                         else if ( (*chain)[j].residue_type == GLY ) { this->ff_ids.push_back(GLY_COMPLETE); }
                         else if ( (*chain)[j].residue_type == HIS ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(HIS_SC); }
                         else if ( (*chain)[j].residue_type == ILE ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(ILE_SC); }
                         else if ( (*chain)[j].residue_type == LEU ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(LEU_SC); }
                         else if ( (*chain)[j].residue_type == LYS ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(LYS_SC); }
                         else if ( (*chain)[j].residue_type == MET ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(MET_SC); }
                         else if ( (*chain)[j].residue_type == PHE ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(PHE_SC); }
                         else if ( (*chain)[j].residue_type == PRO ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(PRO_SC); }
                         else if ( (*chain)[j].residue_type == SER ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(SER_SC); }
                         else if ( (*chain)[j].residue_type == THR ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(THR_SC); }
                         else if ( (*chain)[j].residue_type == TRP ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(TRP_SC); }
                         else if ( (*chain)[j].residue_type == TYR ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(TYR_SC); }
                         else if ( (*chain)[j].residue_type == VAL ) { this->ff_ids.push_back(BACKBONE_NORM); this->ff_ids.push_back(VAL_SC); }
                         else {
                              std::cerr << "#ERROR: unsupported residue: '" << (*chain)[j].residue_type << "'\n";
                              exit(1);
                         }
                    }
               } 
          } //two-body
          //one-body case
          else {
               if ( this->ff_one_body_ids.empty() ) {
                    for (unsigned int j=0; j<chain_len; j++) {
                         if (this->settings.debug>=100) { std::cout << "# residue: "<<j<<"-> "<< (*chain)[j].residue_type << "\n"; }

                         if      ( (*chain)[j].residue_type == ALA ) { this->ff_one_body_ids.push_back(ALA_C); }
                         else if ( (*chain)[j].residue_type == ARG ) { this->ff_one_body_ids.push_back(ARG_C); }
                         else if ( (*chain)[j].residue_type == ASN ) { this->ff_one_body_ids.push_back(ASN_C); }
                         else if ( (*chain)[j].residue_type == ASP ) { this->ff_one_body_ids.push_back(ASP_C); }
                         else if ( (*chain)[j].residue_type == CYS ) { this->ff_one_body_ids.push_back(CYS_C); }
                         else if ( (*chain)[j].residue_type == GLN ) { this->ff_one_body_ids.push_back(GLN_C); }
                         else if ( (*chain)[j].residue_type == GLU ) { this->ff_one_body_ids.push_back(GLU_C); }
                         else if ( (*chain)[j].residue_type == GLY ) { this->ff_one_body_ids.push_back(GLY_C); }
                         else if ( (*chain)[j].residue_type == HIS ) { this->ff_one_body_ids.push_back(HIS_C); }
                         else if ( (*chain)[j].residue_type == ILE ) { this->ff_one_body_ids.push_back(ILE_C); }
                         else if ( (*chain)[j].residue_type == LEU ) { this->ff_one_body_ids.push_back(LEU_C); }
                         else if ( (*chain)[j].residue_type == LYS ) { this->ff_one_body_ids.push_back(LYS_C); }
                         else if ( (*chain)[j].residue_type == MET ) { this->ff_one_body_ids.push_back(MET_C); }
                         else if ( (*chain)[j].residue_type == PHE ) { this->ff_one_body_ids.push_back(PHE_C); }
                         else if ( (*chain)[j].residue_type == PRO ) { this->ff_one_body_ids.push_back(PRO_C); }
                         else if ( (*chain)[j].residue_type == SER ) { this->ff_one_body_ids.push_back(SER_C); }
                         else if ( (*chain)[j].residue_type == THR ) { this->ff_one_body_ids.push_back(THR_C); }
                         else if ( (*chain)[j].residue_type == TRP ) { this->ff_one_body_ids.push_back(TRP_C); }
                         else if ( (*chain)[j].residue_type == TYR ) { this->ff_one_body_ids.push_back(TYR_C); }
                         else if ( (*chain)[j].residue_type == VAL ) { this->ff_one_body_ids.push_back(VAL_C); }
                         else {
                              std::cerr << "#ERROR: unsupported residue: '" << (*chain)[j].residue_type << "'\n";
                              exit(1);
                         }
                    }
               }
          } //end one-body

          if (!settings.one_body_model) {
               this->dummy_bodies_num = this->ff_ids.size();
          }
          else {
               this->dummy_bodies_num = this->ff_one_body_ids.size();
          }

          // output
          if (this->settings.debug>=50) {
               printf("# ff_ids vector: ");
               if (!settings.one_body_model) {
                    for (unsigned int j=0; j<this->dummy_bodies_num; j++) {
                         printf("%d: %d, ",j,ff_ids[j]);
                    }
               }
               else {
                    for (unsigned int j=0; j<this->dummy_bodies_num; j++) {
                         printf("%d: %d, ",j,ff_one_body_ids[j]);
                    }
               }
               printf("\n");
          }
     }

     //! Initialization of internal vectors and matrices: centers, pairwise distances, sinc() matrix, caches
     void saxs_debye_vectors_initialize() {
          // centers
          ff_centers.resize(dummy_bodies_num);
          ff_centers_changed.resize(dummy_bodies_num);
          // dummy initializations, to avoid valgrind warnings (can be skipped)
          for (unsigned int j=0; j<dummy_bodies_num; j++) {
               ff_centers[j]=Vector_3D(0,0,0);
               ff_centers_changed[j]=false;
          }
          // centers distances
          ff_centers_dist.resize(dummy_bodies_num);
          for (unsigned int j=0; j<dummy_bodies_num; j++) {
               ff_centers_dist[j].resize(j+1);
               ff_centers_dist[j][j]=0.0;
          }

          ff_centers_dist_proposal=ff_centers_dist;

          // sinc() multiarray
          ff_sinc.resize(settings.q_bins);
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               ff_sinc[qj].resize(dummy_bodies_num);
               for (unsigned int dbj=0; dbj<dummy_bodies_num; dbj++) {
                    ff_sinc[qj][dbj].resize(dbj+1);
               }
          }
          ff_sinc_proposal=ff_sinc;
          cache_is_outdated=false;
     }

     //! - reads ff_centers
     //! - computes distances
     //! - stores in triangular matrix (std::vector of std::vectors)
     void saxs_debye_ff_pairwise_distances() {
          for (unsigned int j=0; j<dummy_bodies_num; j++) {
               // did this dummy body change position?
               if (ff_centers_changed[j]) {
                    // upper part of the lower triangle, row-wise
                    for (unsigned int k=0; k<j; k++) {
                         ff_centers_dist_proposal[j][k]=( ff_centers[j]-ff_centers[k] ).norm();
                    }
                    // lower part of the lower triangle, column-wise
                    for (unsigned int r=j+1; r<dummy_bodies_num; r++) {
                         // avoid to process twice
                         if (ff_centers_changed[r]) {continue;}
                         ff_centers_dist_proposal[r][j]=( ff_centers[r]-ff_centers[j] ).norm();
                    }
               }
          }
     }

     //! computes the sinc() function: uses the sin() lookup table, and divides for the argument
     //! \param x argument of the sinc() function
     //! \return value of sinc(x)
     double lookup_sinc(double x) {
          double res=1.0;
          if (std::fabs(x) <= 1e-8) {return 1.0;}
          // lookup table?
          if (settings.sine_lookup_table==true) {
               int pos = (int) (x * 651.8986469f) & saxs_definitions::sine_lut_bitmask;
               res =  sine_table[pos]/x;
               // sanity check
               if (settings.debug >= 50) {
                    double sinc_val=sin(x)/x;
                    double diff = res-sinc_val;
                    if (std::fabs(diff) > 1e-3) {
                         printf("# DEBUG: WARN: sine: lookup for x: %.3e (table: pos: %d, val: %.3e) -> sinc: diff: %.3e (table: %.3e, full sinc: %.3e)\n",
                              x,pos,sine_table[pos],diff,res,sinc_val
                         );
                    }
               }
          }
          else {
               // trigonometric function
               res=sin(x)/x;
          }
          return res;
     }

     //! two body version: using scattering angle vector and pairwise distances, computes sinc() for all the terms
     //!
     //! storing the values in settings.q_bins * triangular matrices
     void saxs_debye_ff_sinc_matrices() {
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               // loops over bodies
               for (unsigned int ffj=0; ffj<dummy_bodies_num; ffj++) {
                    // cache is valid?
                    if (ff_centers_changed[ffj]) {
                         // same body?
                         ff_sinc_proposal[qj][ffj][ffj] =
                              form_factors_db[qj][ ff_ids[ffj] ] *
                              form_factors_db[qj][ ff_ids[ffj] ];
                         // upper part of the lower triangle, row-wise
                         for (unsigned int ffk=0; ffk<ffj; ffk++) {
                              ff_sinc_proposal[qj][ffj][ffk] = 2.0 *
                                   form_factors_db[qj][ ff_ids[ffj] ] *
                                   lookup_sinc( q_vec[qj] * ff_centers_dist_proposal[ffj][ffk] ) *
                                   form_factors_db[qj][ ff_ids[ffk] ];
                         }
                         // lower part of the lower triangle, column-wise
                         for (unsigned int ffr=ffj+1; ffr<dummy_bodies_num; ffr++) {
                              // avoid to process twice
                              if (ff_centers_changed[ffr]) {continue;}
                              ff_sinc_proposal[qj][ffr][ffj] = 2.0 *
                                   form_factors_db[qj][ ff_ids[ffr] ] *
                                   lookup_sinc( q_vec[qj] * ff_centers_dist_proposal[ffr][ffj] ) *
                                   form_factors_db[qj][ ff_ids[ffj] ];
                         }
                    }
               }
          }
     }

     //! one body version: using scattering angle vector and pairwise distances, computes sinc() for all the terms
     //!
     //! storing the values in settings.q_bins * triangular matrices
     void saxs_debye_one_body_ff_sinc_matrices() {
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               // loops over bodies
               for (unsigned int ffj=0; ffj<dummy_bodies_num; ffj++) {
                    // cache is valid?
                    if (ff_centers_changed[ffj]) {
                         // same body?
                         ff_sinc_proposal[qj][ffj][ffj] =
                              form_factors_db[qj][ ff_one_body_ids[ffj] ] *
                              form_factors_db[qj][ ff_one_body_ids[ffj] ];
                         // upper part of the lower triangle, row-wise
                         for (unsigned int ffk=0; ffk<ffj; ffk++) {
                              ff_sinc_proposal[qj][ffj][ffk] = 2.0 *
                                   form_factors_db[qj][ ff_one_body_ids[ffj] ] *
                                   lookup_sinc( q_vec[qj] * ff_centers_dist_proposal[ffj][ffk] ) *
                                   form_factors_db[qj][ ff_one_body_ids[ffk] ];
                         }
                         // lower part of the lower triangle, column-wise
                         for (unsigned int ffr=ffj+1; ffr<dummy_bodies_num; ffr++) {
                              // avoid to process twice
                              if (ff_centers_changed[ffr]) {continue;}
                              ff_sinc_proposal[qj][ffr][ffj] = 2.0 *
                                   form_factors_db[qj][ ff_one_body_ids[ffr] ] *
                                   lookup_sinc( q_vec[qj] * ff_centers_dist_proposal[ffr][ffj] ) *
                                   form_factors_db[qj][ ff_one_body_ids[ffj] ];
                         }
                    }
               }
          }
     }

     //! reads the chain object, iterates over residues
     //!   - writes to this->ff_centers (vector of coordinates)
     //!     - centers of backbone: get non-SC atoms -> center of residue -> C,N,CA,O,CB (see below for exact weight of each atom. hydrogens are excluded from calculation)
     //!     - center of SC: atom[PS] (check that C_beta and hydrogens are not present!!)
     //!   - from python code: mass_weights= { 'C':12.0107, 'N':14.0067, 'O':15.9994, 'S':32.065 }
     void saxs_debye_ff_centers() {

          // Import protein definitions (such as residue names)
          using namespace definitions;

          CHAIN_TYPE *chain=this->chain;
          unsigned int chain_len=(*chain).size();
          //
          unsigned int ff_centers_j=0;

          if (settings.one_body_model) {
               for (unsigned int j=0; j<chain_len; j++) {
                    Vector_3D CM =  (*chain)[j].calc_center_of_mass();

                    // positions changed?
                    if ( ff_centers[ff_centers_j]==CM ) {
                         ff_centers_changed[ff_centers_j]=false;
                    }
                    else {
                         ff_centers_changed[ff_centers_j]=true;
                         ff_centers[ff_centers_j]=CM;
                    }

                    ff_centers_j++;

                    if (this->settings.debug>=100) {
                         std::cout<< "id: " << ff_centers_j << ", " <<
                              (*chain)[j].residue_type << " -> " << (*chain)[j].atoms << "\n" <<
                              "    CM: "<< CM << "(changed: " << bool_to_str(ff_centers_changed[j]) << ")\n";
                    }
               }
          }
          // two-body model
          else {
               Vector_3D CM;
               Vector_3D PS_CM; 

               for (unsigned int j=0; j<chain_len; j++) {

                    // weighted by atomic masses from python code, 'C':12.0107, 'N':14.0067, 'O':15.9994
                    CM = (*chain)[j].calc_backbone_center_of_mass();

                    // position: update?
                    if ( ff_centers[ff_centers_j]==CM ) {
                         ff_centers_changed[ff_centers_j]=false;
                    }
                    else {
                         ff_centers_changed[ff_centers_j]=true;
                         ff_centers[ff_centers_j]=CM;
                    }
                    //
                    ff_centers_j++; 
                    //
                    if ((*chain)[j].residue_type != GLY && (*chain)[j].residue_type != ALA) {
                         // position: update?
                         PS_CM = (*chain)[j].calc_sidechain_center_of_mass();
                         if ( ff_centers[ff_centers_j] == PS_CM ) {
                             ff_centers_changed[ff_centers_j]=false;
                         }
                         else {
                             ff_centers_changed[ff_centers_j]=true;
                             ff_centers[ff_centers_j]=PS_CM;
                         }
                         //
                         ff_centers_j++;
                    }

                    if (this->settings.debug>=100) {
                         std::cout<< "id: " << ff_centers_j << ", " <<
                              (*chain)[j].residue_type << " -> " << (*chain)[j].atoms << "\n" <<
                              "    backbone: "<< CM << "(changed: " << bool_to_str(ff_centers_changed[j]) <<
                              "), side chain: " << PS_CM <<
                              " (changed: " << bool_to_str(ff_centers_changed[j]) << ")\n";
                    }
               }
          } 
     }

     //! prints all structures
     //! - scattering vector
     //! - reference SAXS profile
     //! - parameters of Normal distribution describing the error
     //! - form factors
     //! - bodies centers
     //! - cache of sinc() elementary pairwise contributions
     //! - proposed SAXS profile
     //! - total monte carlo energy
     void print_all() {
          printf("########################################################################\n");
          printf("# SAXS: q_vec\n");
          for (unsigned int j=0; j<settings.q_bins; j++) { printf("%.3f, ",q_vec[j]); }
          printf("\n");
          //
          printf("# SAXS: reference\n");
          for (unsigned int j=0; j<settings.q_bins; j++) { printf("%.3f, ",saxs_profile_reference[j]); }
          printf("\n");
          //
          printf("# SAXS: 'exp' errors\n");
          for (unsigned int j=0; j<settings.q_bins; j++) { printf("%.3f, ",saxs_profile_reference_errors[j]); }
          printf("\n");
          //
          printf("# SAXS: ff\n");
          int nof_ffs = settings.one_body_model ? 20 : 21; 
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               printf("q: %.3f -> ",q_vec[qj]);
               for (int ffj=0; ffj<nof_ffs; ffj++) {
                    printf("%.2f, ",form_factors_db[qj][ffj]);
               }
               printf("\n");
          }
          printf("\n");
          //
          printf("# SAXS: ff centers\n");
          if (!settings.one_body_model) {
               for (unsigned int j=0; j<dummy_bodies_num; j++) {
                    printf("%3d -> %2d -> %.3f %.3f %.3f (changed: %s)\n",j,ff_ids[j],
                         ff_centers[j][0],ff_centers[j][1],ff_centers[j][2],
                         bool_to_str(ff_centers_changed[j])
                    );
               }
          }
          else {
               for (unsigned int j=0; j<dummy_bodies_num; j++) {
                    printf("%3d -> %2d -> %.3f %.3f %.3f (changed: %s)\n",j,ff_one_body_ids[j],
                        ff_centers[j][0],ff_centers[j][1],ff_centers[j][2],
                        bool_to_str(ff_centers_changed[j])
                    );
               }
          }
          //
          unsigned int qj=1;
          printf("# SAXS: ff_sinc: scattering momentum q[%d]=%.4f\n",qj,q_vec[qj]);
          for (unsigned int ffj=0; ffj<dummy_bodies_num; ffj++) {
               printf("# ffj: %3d: ",ffj);
               for (unsigned int ffk=0; ffk<=ffj; ffk++) {
                    printf("%.2f, ",ff_sinc[qj][ffj][ffk]);
               }
               printf("\n");
          }
          //
          printf("# SAXS: profile\n");
          printf("# index -> q -> saxs , crysol\n");
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               printf("%d -> %.3f  -> %.3f %.3f\n",qj,q_vec[qj],saxs_profile[qj],saxs_profile_reference[qj]);
          }
          //
          print_energy();
     }

     //! prints monte carlo energy
     void print_energy() {
          printf("# SAXS energy: %.3f\n",current_energy);
     }

     //! saves a file with the computed SAXS profile
     void dump_internals(const char *filename_pre) {
          // dumps q_vec, SAXS intensity for all the q bins
          char filename[2000];
          sprintf(filename,"%s.%s",filename_pre,"int");

          FILE *f_h;
          f_h = fopen(filename,"w");

          if (!f_h) {
               printf("#ERROR: saxs output: unable to open file '%s'\n", filename );
               return;
          }

          // write scattering curve to file
          fprintf(f_h,"#  q         i \n");
          for (unsigned int qj=0; qj<settings.q_bins; qj++) {
               fprintf(f_h,"%.3f  %.4f\n",q_vec[qj],saxs_profile[qj]);
          }
          fclose(f_h);
     }

     //! evaluation of the monte carlo energy: overwrites virtual function of base class
     double evaluate(MoveInfo *moveInfo=NULL) {
          // SAXS: Debye formula
          saxs_debye_prepare();
          saxs_debye_compute();

          this->current_energy=evaluate_do( &saxs_profile );

          if (this->settings.debug>=100) { print_all(); }
          else if (this->settings.debug>=10) { print_energy(); }

          return this->current_energy;
     }
};

}

#endif
