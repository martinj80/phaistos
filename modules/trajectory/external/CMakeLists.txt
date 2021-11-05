# set root directory
set(LIBXDRFILE_DIR "xdrfile-1.1b")

set(phaistos_MODULE_INCLUDE_DIRS ${phaistos_MODULE_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/${LIBXDRFILE_DIR}/include")
set(phaistos_MODULE_INCLUDE_DIRS ${phaistos_MODULE_INCLUDE_DIRS} PARENT_SCOPE)

# Add root directory to the list of includes
include_directories (${LIBXDRFILE_DIR}/include) 

add_library(xdrfile
  ${LIBXDRFILE_DIR}/src/xdrfile.c
  ${LIBXDRFILE_DIR}/src/xdrfile_xtc.c
)