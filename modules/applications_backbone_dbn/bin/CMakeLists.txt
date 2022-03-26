# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/bin ${phaistos_SOURCE_DIR}/src/models/backbone_dbn) 

# Make sure the linker can find the libphaistos library once it is built. 
link_directories (${phaistos_BINARY_DIR}/libs)

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

if(Boost_FOUND)

  include_directories(${Boost_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS})

  # Create executable: check_mutations
  add_executable(check_mutations check_mutations.cpp)

  # Link executable with libraries
  target_link_libraries(check_mutations libphaistos dbn integrate ${Boost_LIBRARIES})


  # Create executable: entropy
  add_executable(entropy entropy.cpp)

  # Link executable with libraries
  target_link_libraries(entropy libphaistos dbn integrate ${Boost_LIBRARIES})


  # Create executable: likelihood
  add_executable(likelihood likelihood.cpp)

  # Link executable with libraries
  target_link_libraries(likelihood libphaistos dbn integrate ${Boost_LIBRARIES})


  # Create executable: predictor
  add_executable(predictor predictor.cpp)

  # Link executable with libraries
  target_link_libraries(predictor libphaistos dbn integrate ${Boost_LIBRARIES})


  # Create executable: sampler
  add_executable(sampler sampler.cpp)

  # Link executable with libraries
  target_link_libraries(sampler libphaistos dbn integrate ${Boost_LIBRARIES})

  INSTALL( TARGETS check_mutations entropy likelihood predictor sampler DESTINATION bin )
endif()
