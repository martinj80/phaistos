# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/src/models/backbone_dbn ${CMAKE_CURRENT_SOURCE_DIR}/../src ${MOCAPY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS} ${MUNINN_INCLUDE_DIRS}) 

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Create executable: pdb2git
add_executable(pdb2git pdb2git.cpp)
add_dependencies(pdb2git svnrev)

# Link executable with library
target_link_libraries(pdb2git libphaistos integrate ${LAPACK_LIBRARY} ${BLAS_LIBRARY} git ${MUNINN_LIBRARIES} ${CMAKE_FLIB} ${Boost_LIBRARIES})

INSTALL(TARGETS pdb2git DESTINATION bin)
