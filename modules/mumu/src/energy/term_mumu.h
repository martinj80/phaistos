// term_mumu.h --- MUMU Multibody energy term
// Copyright (C) 2012-2013 Kristoffer Enøe Johansson
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

#ifndef TERM_MUMU_H
#define TERM_MUMU_H

#include "mocapy.h"
#include "models/mumu/mumu_data_miner.h"
#include "models/mumu/mumu_inf.h"

#include "energy/energy_term.h"

namespace phaistos {

//! MUMU multibody energy term - base class
template <typename DERIVED_CLASS, typename CHAIN_TYPE, typename MUMU_TYPE>
class TermMumuBase:public EnergyTermCommon<DERIVED_CLASS, CHAIN_TYPE> {
private:

     // For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS,CHAIN_TYPE> EnergyTermCommon;

     //! Index of hidden node in mismask
     int hidden_index;
     
     //! Initialize mumu term
     void init(CHAIN_TYPE *chain) {

          // Save chain
          this->chain = chain;
          chain_length = chain->size();
          cache.resize(chain_length);
          
          // Init data miner
          data_miner = new MUMU_TYPE;
          data_miner->init(chain);

          // hidden node mismask index
          this->hidden_index = data_miner->hidden_index;
          
          // Get missing masks pointers from data_miner
          mismask_joint = &data_miner->mismask_joint;
          mismask_margi = &data_miner->mismask_margi;
          
          // Init model
          model = new mocapy::DBN;
          std::string filename;
          if ( settings.model_filename == "default" ) {
               filename = settings.mocapy_dbn_dir+"/"+data_miner->default_dbn_filename;
          } else {
               filename = settings.mocapy_dbn_dir+"/"+settings.model_filename;
          }
//          std::cout<<"# MuMu energy init: Loading dbn file "<<filename<<std::endl;
          model->load(filename.c_str());
          if ( settings.use_ratio ) {
                if ( settings.reference_model_filename == "default" ) {
                     filename = settings.mocapy_dbn_dir+"/"+data_miner->default_reference_dbn_filename;
                } else {
                     filename = settings.mocapy_dbn_dir+"/"+settings.reference_model_filename;
                }
                ref_model = new mocapy::DBN;
                // std::cout<<"# MuMu energy init: Loading reference dbn file "<<filename<<std::endl;
                ref_model->load(filename.c_str());
          } else {
               ref_model = NULL;
          }
          
          // Init inference engine
          inf_engine = new MumuTest14Inf(model);
     };
     
protected:

     mocapy::DBN *model;
     mocapy::DBN *ref_model;
     MumuTest14Inf *inf_engine;

     MUMU_TYPE *data_miner;

     mocapy::MDArray<mocapy::eMISMASK> *mismask_joint, *mismask_margi;

     int chain_length;
     std::vector<double> cache;
     
public:

     //! Local settings class.     
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

          //! Whether to use the ratio method
          bool use_ratio;

          //! Path to model file
          std::string model_filename;

          //! Path to reference model file 
          std::string reference_model_filename;

          //! Directory in which to search for model files
          std::string mocapy_dbn_dir;

          //! Constructor
          Settings(bool use_ratio=false,
                   std::string model_filename="default",
                   std::string reference_model_filename="default",
                   std::string mocapy_dbn_dir="../data/mocapy_dbns")
               : use_ratio(use_ratio),
                 model_filename(model_filename),
                 reference_model_filename(reference_model_filename),
                 mocapy_dbn_dir(mocapy_dbn_dir) {
          }

          //! Output operator          
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "use-ratio:" << settings.use_ratio << "\n";
               o << "model-filename:" << settings.model_filename << "\n";
               o << "reference-model-filename:" << settings.reference_model_filename << "\n";
               o << "mocapy-dbn-dir:" << settings.mocapy_dbn_dir << "\n";
               o << static_cast<const typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
               return o;
          }                    
     } settings;  //!< Local settings object 


     //! Constructor
     //! \param chain Molecule chain
     //! \param name Name of term
     //! \param settings Local Settings object     
     TermMumuBase(CHAIN_TYPE *chain, std::string name, const Settings &settings=Settings(),
                  RandomNumberEngine *random_number_engine = &random_global) 
          : EnergyTermCommon(chain, name, settings, random_number_engine),
            settings(settings) {
          
          // construct model(s)
          init(chain);
     };
     
     //! Copy Constructor
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermMumuBase(const TermMumuBase &other, RandomNumberEngine *random_number_engine,
                  int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            settings(other.settings) {
          
          // Copy everything and make new data_miner, model and inference engine
          this->init(chain);
     };

     //! Destructor
     ~TermMumuBase() {
          if (model) {
               delete model;
               model = NULL;
          }
          if (inf_engine) {
               delete inf_engine;
               inf_engine = NULL;
          }
          if (ref_model) {
               delete ref_model;
               ref_model = NULL;
          }
          if (data_miner) {
               delete data_miner;
               data_miner = NULL;
          }
     };

     //! Evaluate energy term.
     //! \param move_info Object containing information about the last executed move     
     double evaluate(MoveInfo *move_info=NULL) {
          double ll = 0.0;
          mocapy::Sequence *data;

          data_miner->update(move_info);

          for (int i=0; i<chain_length; i++) {
               // get data and set mismasks
               data = data_miner->make_featurevector(i);

               // calc conditional log likelihood for mumu model
               double res_ll = inf_engine->calc_aa_ll(data->get_view(0).get_values());
               cache[i] = -res_ll;
               ll += res_ll;
          }

          return -ll;
     }
};


//! MUMU multibody energy term 
template <typename CHAIN_TYPE>
class TermMumu:public TermMumuBase<TermMumu<CHAIN_TYPE>, CHAIN_TYPE, MumuDataMinerTest14> {
public:

     // For convenience, define local EnergyTermCommon
     typedef phaistos::TermMumuBase<TermMumu<CHAIN_TYPE>,CHAIN_TYPE,MumuDataMinerTest14> TermMumuBase;

     // Use settings from base class
     typedef typename TermMumuBase::Settings Settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object     
     TermMumu(CHAIN_TYPE *chain, const Settings &settings=Settings(),
              RandomNumberEngine *random_number_engine = &random_global) 
          : TermMumuBase(chain, "mumu", settings, random_number_engine) {}
     
     //! Copy Constructor
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermMumu(const TermMumu &other,
              RandomNumberEngine *random_number_engine,
              int thread_index, CHAIN_TYPE *chain)
          : TermMumuBase(other, random_number_engine, thread_index, chain) {};
};

//! Observable specialization for TermMumu
template <typename CHAIN_TYPE>
class Observable<TermMumu<CHAIN_TYPE> >: public TermMumu<CHAIN_TYPE>, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermMumu<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename TermMumu<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<const ObservableBase::Settings>(settings);
               return o;
          }          
     } settings; //!< Local settings object 
     
     //! Constructor.
     //! \param energy_term VisibleVolume energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermMumu<CHAIN_TYPE> &energy_term,
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermMumu<CHAIN_TYPE>(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index,
                typename TermMumu<CHAIN_TYPE>::ChainType *chain)
          : TermMumu<CHAIN_TYPE>(other, thread_index, chain),
            settings(other.settings) {
     }     


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermMumu<CHAIN_TYPE> *clone(int thread_index=0,
                                 typename TermMumu<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermMumu<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL,
                                 PHAISTOS_LONG_LONG current_iteration=0,
                                 bool register_only=false) {

          this->evaluate(move_info);

          std::string output = "";

          for (ResidueIterator<CHAIN_TYPE> it(*this->chain); !it.end(); ++it) {
               std::string output_entry = ObservableBase::vector_output_tag(*it);
               
               output_entry += boost::lexical_cast<std::string>(this->cache[it->index]);
               
               if (output != "")
                    output += ",";
               output += output_entry;
          }

          return output;
     }
     
};

}

#endif
