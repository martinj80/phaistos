#include "energy/term_constrain_distances.h"

namespace phaistos {

// This definition needs to appear right after the includes and is therefore included here rather than global.cpp
namespace module_typhon {

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
          if (target.has_key("mode") && target["mode"].as<std::string>() == "typhon") {
               implicit_energies_allowed = true;
               atom_types_default  = "ALL_PHYSICAL_ATOMS";
               target.super_group_default("procedure") = "--procedure fold";
               target.super_group_default("monte carlo") = "--monte-carlo metropolis-hastings";
               target.super_group_default("move")   = "--move crisp sc-basilisk";
               target.super_group_default("energy") = "--energy clash-fast local-dbn opls-angle-bend-cached[omit-sidechains:1] constrain-distances[verbose:1]";
               // target.super_group_default("move")   = "--move crisp-dbn-eh sc-basilisk";
               // target.super_group_default("energy") = "--energy clash-fast opls-angle-bend-cached constrain-distances[verbose:1]";
          }
     }
};

}
}

