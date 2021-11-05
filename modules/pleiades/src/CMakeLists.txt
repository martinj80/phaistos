require_module("git")

# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR} ${MUNINN_INCLUDE_DIRS} ${MOCAPY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Define basilisk library
add_library(libpleiades
  ${CMAKE_CURRENT_SOURCE_DIR}/cluster_git.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/cluster_git.cpp 
  ${CMAKE_CURRENT_SOURCE_DIR}/cluster_rmsd.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/cluster_rmsd.cpp 
  ${CMAKE_CURRENT_SOURCE_DIR}/git_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/git_element.h
)

# Tell cmake not to add additional "lib" prefix to library (it's already there)
set_target_properties (libpleiades PROPERTIES PREFIX "")

# Add libraries to phaistos_MODULE_LIBS variable
list(APPEND 
  phaistos_MODULE_LIBS
  libpleiades
)
