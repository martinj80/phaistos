# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(libphaistos
  protein/atom.cpp 
  protein/atom.h 
  protein/residue.cpp 
  protein/residue.h 
  protein/chain.h 
  protein/chain_ca.h 
  protein/chain_fb.h 
  protein/iterators/atom_iterator.h 
  protein/iterators/atom_pair_iterator.h
  protein/iterators/covalent_bond_iterator.h
  protein/iterators/dof_iterator.h
  protein/iterators/iterator_base.h
  protein/iterators/pair_iterator_chaintree.h
  protein/iterators/residue_iterator.h 
  protein/pdb_input.cpp
  protein/protein_ideal_values.h
  protein/protein_data.h
  utils/eigen_system_3x3.h
  utils/eigen_system_3x3.cpp
  utils/math.h
  utils/math.cpp
  moves/move.h
  moves/move_collection.h 
  moves/move_fixed_structure.h 
  energy/energy.h 
  energy/energy_term.h 
  energy/term_contact_map.h
  energy/term_angle_histogram.h
  energy/term_backbone_dbn.h
  energy/term_clash_fast.h
  energy/term_contact_map.h
  energy/term_energy_sum.h
  energy/term_energy_term_wrapper.h
  energy/term_helix_content.h
  energy/term_q_factor.h
  energy/term_rg.h 
  energy/term_rmsd.h 
  utils/matrix.h
  utils/matrix_vec.h
  utils/program_option_parser.cpp
  utils/svd.cpp
  utils/debug_levels.h
  # utils/math_utils.h
  utils/optimize.h
  utils/utils.h
  utils/vector_matrix_3d.h
  utils/vector_nd.h
  ${phaistos_MODULE_SRC}
)

# Tell cmake not to add additional "lib" prefix to library (it's already there)
include_directories(${phaistos_MODULE_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

set_target_properties (libphaistos PROPERTIES PREFIX "")

add_subdirectory(utils)
add_subdirectory(models)
