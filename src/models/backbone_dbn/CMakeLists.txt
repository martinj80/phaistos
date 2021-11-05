# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add src/ directory to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src) 

# Add boost include directories
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()


add_library(dbn
  parameters.h
  backbone_dbn.h
  nodes/node.h
  nodes/fb5.cpp 
  nodes/fb5.h
  nodes/torus.cpp
  nodes/torus.h
  nodes/discrete.h
  nodes/discrete.cpp
  nodes/hidden.h
  nodes/vonmises.h
  nodes/vonmises.cpp
  nodes/gaussian.h
  nodes/gaussian.cpp)
