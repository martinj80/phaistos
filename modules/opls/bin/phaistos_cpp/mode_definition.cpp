
// Insert if not already present
if (std::find(phaistos_modes.begin(), phaistos_modes.end(), "opls-mc-dynamics") == phaistos_modes.end())
     phaistos_modes.push_back("opls-mc-dynamics");

// Module: Initialize mode definition (defined in includes.cpp)
module_opls::ModeDefinitionInitialization(target, chain, atom_types_default, implicit_energies_allowed);
