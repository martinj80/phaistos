# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/bin ${phaistos_SOURCE_DIR}/src/models/backbone_dbn) 

# Make sure the linker can find the libphaistos library once it is built. 
link_directories (${phaistos_BINARY_DIR}/libs)

# Find pthreads
find_library(PTHREAD pthread)

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

if(Boost_FOUND)

  include_directories(${Boost_INCLUDE_DIRS} ${MUNINN_INCLUDE_DIRS} ${MOCAPY_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS})

  # Create executable: phaistos
  add_executable(phaistos phaistos.cpp)
  add_dependencies(phaistos svnrev)

  # Link executable with libraries
  target_link_libraries(phaistos libphaistos dbn integrate ${phaistos_MODULE_LIBS} ${MUNINN_LIBRARIES} ${MOCAPY_LIBRARIES} ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${PTHREAD} ${CMAKE_FLIB} ${Boost_LIBRARIES})

  INSTALL(TARGETS phaistos DESTINATION bin)
endif()


