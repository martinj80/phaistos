
// Insert if not already present
if (std::find(phaistos_modes.begin(), phaistos_modes.end(), "typhon") == phaistos_modes.end())
     phaistos_modes.push_back("typhon");

// Module: Initialize mode definition (defined in includes.cpp)
module_typhon::ModeDefinitionInitialization(target, chain, atom_types_default, implicit_energies_allowed);
