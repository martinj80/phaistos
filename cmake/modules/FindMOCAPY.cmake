# Try to find Mocapy library
# This will define:
#
#  MOCAPY_FOUND - system has Mocapy
#  MOCAPY_INCLUDE_DIRS - the Mocapy include directories
#  MOCAPY_LIBRARIES - link these to use Mocapy


# Macros with tools for finding libraries
# include(LibFindMacros)
set(_MOCAPY_IN_CACHE FALSE)
if(MOCAPY_INCLUDE_DIR AND MOCAPY_LIBRARY)
  set(_MOCAPY_IN_CACHE TRUE)
endif(MOCAPY_INCLUDE_DIR AND MOCAPY_LIBRARY)

# Include dir
find_path(MOCAPY_INCLUDE_DIR
  NAMES mocapy.h
  PATHS ${MOCAPY_ROOT}/src
)
if (MOCAPY_INCLUDE_DIR)
  set(MOCAPY_INCLUDE_DIRS ${MOCAPY_INCLUDE_DIR})
  set(MOCAPY_INCLUDE_DIRS ${MOCAPY_INCLUDE_DIRS} PARENT_SCOPE)
endif (MOCAPY_INCLUDE_DIR)

# The library
find_library(MOCAPY_LIBRARY
  NAMES Mocapy
  PATHS ${MOCAPY_ROOT}/libs ${MOCAPY_ROOT}/src
  PATH_SUFFIXES . .libs
  NO_DEFAULT_PATH
)
if (MOCAPY_LIBRARY)
  set(MOCAPY_LIBRARIES ${MOCAPY_LIBRARY})
  set(MOCAPY_LIBRARIES ${MOCAPY_LIBRARIES} PARENT_SCOPE)
endif (MOCAPY_LIBRARY)


if (MOCAPY_INCLUDE_DIR AND MOCAPY_LIBRARY)

  # Report if values were not found previously
  if (NOT _MOCAPY_IN_CACHE) 
    if (NOT MOCAPY_FIND_QUIETLY)
      message(STATUS "Found Mocapy: ${MOCAPY_LIBRARY}")
    endif (NOT MOCAPY_FIND_QUIETLY)  
  endif (NOT _MOCAPY_IN_CACHE)

  set(MOCAPY_FOUND TRUE)
  add_definitions(-DHAVE_MOCAPYLIB)
  set(_MOCAPY_IN_CACHE TRUE)
endif (MOCAPY_INCLUDE_DIR AND MOCAPY_LIBRARY)

if (NOT MOCAPY_FOUND)
   if (MOCAPY_FIND_REQUIRED)
     message("Could not find Mocapy (MOCAPY_ROOT = ${MOCAPY_ROOT})")
   # else (MOCAPY_FIND_REQUIRED)     
   #   message (STATUS "Mocapy not found (MOCAPY_ROOT = ${MOCAPY_ROOT})")
   endif (MOCAPY_FIND_REQUIRED)
endif(NOT MOCAPY_FOUND)