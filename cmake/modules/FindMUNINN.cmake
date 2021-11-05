# Try to find Muninn library
# This will define:
#
#  MUNINN_FOUND - system has Muninn
#  MUNINN_INCLUDE_DIRS - the Muninn include directories
#  MUNINN_LIBRARIES - link these to use Muninn
#  MUNINN_DEFINITIONS - definitions set by Muninn (HAVE_MUNINNLIB). You should
#                       add_definitions(${MUNINN_DEFINITIONS}) before compiling
#                       code that includes Muninn. 

set(_MUNINN_IN_CACHE FALSE)
if(MUNINN_INCLUDE_DIR AND MUNINN_LIBRARY)
  set(_MUNINN_IN_CACHE TRUE)
endif(MUNINN_INCLUDE_DIR AND MUNINN_LIBRARY)

# Include dir
find_path(MUNINN_INCLUDE_DIR
  NAMES muninn/GE.h
  PATHS ${MUNINN_ROOT}
)
if (MUNINN_INCLUDE_DIR)
  set(MUNINN_INCLUDE_DIRS ${MUNINN_INCLUDE_DIR})
  set(MUNINN_INCLUDE_DIRS ${MUNINN_INCLUDE_DIR} PARENT_SCOPE)

  set(MUNINN_FOUND TRUE)  
  set(MUNINN_DEFINITIONS "HAVE_MUNINNLIB") 

endif (MUNINN_INCLUDE_DIR)

# The library
find_library(MUNINN_LIBRARY
  NAMES muninn
  PATHS ${MUNINN_ROOT}/muninn
  PATH_SUFFIXES . .libs
  NO_DEFAULT_PATH
)
if (MUNINN_LIBRARY)
  set(MUNINN_LIBRARIES ${MUNINN_LIBRARY})
  set(MUNINN_LIBRARIES ${MUNINN_LIBRARY} PARENT_SCOPE)
endif (MUNINN_LIBRARY)

if (MUNINN_INCLUDE_DIR AND MUNINN_LIBRARY)

  # Report if values were not found previously
  if (NOT _MUNINN_IN_CACHE) 
    if (NOT MUNINN_FIND_QUIETLY)
      message(STATUS "Found Muninn: ${MUNINN_LIBRARY}")
    endif (NOT MUNINN_FIND_QUIETLY)  
  endif (NOT _MUNINN_IN_CACHE)

  set(_MUNINN_IN_CACHE TRUE)
endif (MUNINN_INCLUDE_DIR AND MUNINN_LIBRARY)

if (NOT MUNINN_FOUND)
   if (MUNINN_FIND_REQUIRED)
     message(FATAL_ERROR "Could not find Muninn (MUNINN_ROOT = ${MUNINN_ROOT})")
   # else (MUNINN_FIND_REQUIRED)
   #   message (STATUS "Muninn not found (MUNINN_ROOT = ${MUNINN_ROOT})")
   endif (MUNINN_FIND_REQUIRED)
endif(NOT MUNINN_FOUND)