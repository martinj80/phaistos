# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR} ${MUNINN_INCLUDE_DIRS} ${MOCAPY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Define basilisk library
add_library(git
  ${CMAKE_CURRENT_SOURCE_DIR}/git.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/energy/observable_git.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/git/git.cpp 
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/git/git.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/git/git_parameters.h 
)

# Add libraries to phaistos_MODULE_LIBS variable
list(APPEND 
  phaistos_MODULE_LIBS
  git
)
