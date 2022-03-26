#include "energy/term_opls_charge.h"
#include "energy/term_opls_vdw.h"
#include "energy/term_opls_angle_bend.h"
#include "energy/term_opls_torsion.h"
#include "energy/term_opls_imptor.h"
#include "energy/term_opls_bond_stretch.h"
#include "energy/term_opls_non_bonded.h"
#include "energy/term_gbsa.h"

namespace phaistos {

// This definition needs to appear right after the includes and is therefore included here rather than global.cpp
namespace module_opls {

// Module: mode definition
struct ModeDefinitionInitialization {

     // General case - do nothing
     template <typename CHAIN_TYPE>
     ModeDefinitionInitialization(ProgramOptionParser &target, CHAIN_TYPE *chain, 
                                  std::string &atom_types_default, bool &implicit_energies_allowed) {
     }

     // ChainFB case
     ModeDefinitionInitialization(ProgramOptionParser &target, ChainFB *chain,
                                  std::string &atom_types_default, bool &implicit_energies_allowed) {
          if (target.has_key("mode") && target["mode"].as<std::string>() == "opls-mc-dynamics") {
               implicit_energies_allowed = false;
               atom_types_default  = "ALL_PHYSICAL_ATOMS";
               target.super_group_default("procedure") = "--procedure fold";
               target.super_group_default("monte carlo") = "--monte-carlo metropolis-hastings[declash-on-reinitialize:false,reinitialization-interval:0]";
               target.super_group_default("move")   = "--move crisp[weight:0.2,std-dev-bond-angle:0.5,std-dev-phi-psi:4,std-dev-omega:0.5] pivot-local[weight:0.05,std-dev-bond-angle:0.8,std-dev-phi-psi:1,std-dev-omega:0.8] sidechain-rotamer[weight:0.5,rotamer-state-resample-frequency:0.33] sidechain-local[weight:0.125,sample-minor-dofs:1,sigma-major-dofs:0.1,mode:constrain-one-endpoint,lagrange-multiplier:200] sidechain-local[weight:0.125,sample-minor-dofs:1,sigma-major-dofs:0.1,mode:constrain-one-endpoint,lagrange-multiplier:500000]";
               target.super_group_default("energy") = "--energy opls-cached";
          }
     }
};

}
}
