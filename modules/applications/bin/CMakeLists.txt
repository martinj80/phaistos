require_module("trajectory")

# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/bin ${phaistos_SOURCE_DIR}/src/models/backbone_dbn) 

# Make sure the linker can find the libphaistos library once it is built. 
link_directories (${phaistos_BINARY_DIR}/libs)

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

if(Boost_FOUND)

  include_directories(${Boost_INCLUDE_DIRS} ${MUNINN_INCLUDE_DIRS} ${MOCAPY_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS})

  # Create executable: evaluate_observable
  add_executable(evaluate_observable evaluate_observable.cpp)

  # Link executable with libraries
  target_link_libraries(evaluate_observable libphaistos dbn integrate ${phaistos_MODULE_LIBS} ${MUNINN_LIBRARIES} ${MOCAPY_LIBRARIES} ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${PTHREAD} ${CMAKE_FLIB} ${Boost_LIBRARIES})


  # Create executable: calc_rmsd
  add_executable(calc_rmsd calc_rmsd.cpp)

  # Link executable with libraries
  target_link_libraries(calc_rmsd libphaistos ${LAPACK_LIBRARY})


  # Create executable: calc_rg
  add_executable(calc_rg calc_rg.cpp)

  # Link executable with libraries
  target_link_libraries(calc_rg libphaistos ${LAPACK_LIBRARY})


  # Create executable: calc_q_score
  add_executable(calc_q_score calc_q_score.cpp)

  # Link executable with libraries
  target_link_libraries(calc_q_score libphaistos ${LAPACK_LIBRARY})

  INSTALL(TARGETS evaluate_observable calc_rmsd calc_rg calc_q_score DESTINATION bin)
endif()
