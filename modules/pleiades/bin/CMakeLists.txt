# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/src/models/backbone_dbn ${CMAKE_CURRENT_SOURCE_DIR}/../src ${MOCAPY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS} ${MUNINN_INCLUDE_DIRS}) 

#
# Create executable: pleiades 
add_executable(pleiades pleiades.cpp pleiades_options.cpp pleiades_options.h)
add_dependencies(pleiades svnrev)

# Link executable with library
target_link_libraries(pleiades libpleiades git libphaistos integrate ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${MUNINN_LIBRARIES} ${CMAKE_FLIB} ${Boost_LIBRARIES})

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Create executable: pleiades RMSD
add_executable(pleiades_rmsd pleiades_rmsd.cpp pleiades_rmsd_options.cpp pleiades_rmsd_options.h)
add_dependencies(pleiades_rmsd svnrev)

# Link executable with library
target_link_libraries(pleiades_rmsd libpleiades libphaistos dbn integrate ${LAPACK_LIBRARY} ${BLAS_LIBRARY} git ${MUNINN_LIBRARIES} ${CMAKE_FLIB} ${Boost_LIBRARIES})

INSTALL(TARGETS pleiades pleiades_rmsd DESTINATION bin )
