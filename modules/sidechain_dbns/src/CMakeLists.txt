require_module("mocapy")

# Add source files
list(APPEND 
  phaistos_MODULE_SRC
  energy/term_sidechain_dbn.h
)

# Define basilisk library
add_library(libbasilisk
  ${CMAKE_CURRENT_SOURCE_DIR}/models/basilisk/basilisk_dbn.h
  ${CMAKE_CURRENT_SOURCE_DIR}/models/basilisk/basilisk_dbn.cpp 
)
# Remove "lib" prefix
set_target_properties (libbasilisk PROPERTIES PREFIX "")

# Copy data directory to build directory
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/models/basilisk/data
                                                           ${CMAKE_BINARY_DIR}/src/models/basilisk/data)


# Define basilisk library
add_library(libcompas
  ${CMAKE_CURRENT_SOURCE_DIR}/models/compas/compas_dbn.h 
  ${CMAKE_CURRENT_SOURCE_DIR}/models/compas/compas_dbn.cpp
)
# Remove "lib" prefix
set_target_properties (libcompas PROPERTIES PREFIX "")

# Copy data directory to build directory
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/models/compas/data
                                                           ${CMAKE_BINARY_DIR}/src/models/compas/data)


# Set include directories
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR} ${MOCAPY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Add libraries to phaistos_MODULE_LIBS variable
list(APPEND 
  phaistos_MODULE_LIBS
  libbasilisk
  libcompas
)

