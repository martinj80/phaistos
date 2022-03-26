# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/bin ${phaistos_SOURCE_DIR}/src/models/backbone_dbn ${Boost_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS} ${MUNINN_INCLUDE_DIRS}) 

# Make sure the linker can find the libphaistos library once it is built. 
link_directories (${phaistos_BINARY_DIR}/libs)

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Create executable
add_executable(trajectory_subset trajectory_subset.cpp)

# Link executable with libraries
target_link_libraries(trajectory_subset libphaistos ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${Boost_LIBRARIES} ${phaistos_MODULE_LIBS} ${MUNINN_LIBRARIES} )

