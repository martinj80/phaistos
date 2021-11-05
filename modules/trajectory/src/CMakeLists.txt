add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../external)

list(APPEND 
  phaistos_MODULE_SRC
  energy/observable_xtc_trajectory.h
)

# Add libraries to phaistos_MODULE_LIBS variable
list(APPEND 
  phaistos_MODULE_LIBS
  xdrfile
)

