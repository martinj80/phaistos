# The name of our project is "MOCAPY". CMakeLists files in this project can
# refer to the root source directory of the project as ${MOCAPY_SOURCE_DIR} and
# to the root binary directory of the project as ${MOCAPY_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6.4)
project (MOCAPY C CXX Fortran)
# version number
set(MOCAPY_VERSION 0.1)

# Find Boost package
SET(Boost_ADDITIONAL_VERSIONS "1.45" "1.45.0" "1.44" "1.44.0" "1.43" "1.43.0" "1.42" "1.42.0" "1.41" "1.41.0" "1.40" "1.40.0" "1.39" "1.39.0")
SET(Boost_USE_MULTITHREAD OFF)
FIND_PACKAGE(Boost COMPONENTS serialization program_options thread REQUIRED)

# find lapack
FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack lapack-3 )
IF(NOT LAPACK_LIBRARY)
message("error: could not find lapack library")
ENDIF(NOT LAPACK_LIBRARY)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_FLIB "gfortran")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_FLIB "g2c")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
#  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#  message ("No optimized Fortran link flags are known, we just try -lgfortran")
  set (CMAKE_FLIB "gfortran")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

# Recurse into the "src" and "examples" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)
add_subdirectory (examples)
add_subdirectory (tests)

# settings for making distribution package
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Mocapy is a Dynamic Bayesian Network toolkit")
set(CPACK_PACKAGE_VENDOR "Thomas Hamelryck and Martin Paluszewski")
#set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README)
set(CPACK_GENERATOR TGZ)
set(
        CPACK_SOURCE_PACKAGE_FILE_NAME
        "Mocapy-${MOCAPY_VERSION}"
        CACHE INTERNAL "tarball basename"
        )
set(CPACK_SOURCE_GENERATOR TGZ)
include(CPack)

