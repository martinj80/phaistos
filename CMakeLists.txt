# Project name
project (phaistos)

# Set version number
set (VERSION 1.0)

set( LICENSE_FILE "COPYING" )
enable_language (Fortran)

# Minimum required version of cmake
cmake_minimum_required (VERSION 2.6.2)

# Only allow out-of-source build 
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a build directory (e.g. build) and run CMake from there.")
endif()

# Set default make type
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" )
endif()

# Set parameters for different build types
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb")

# Add local cmake module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Require proper c++ code
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ansi -pedantic -Wno-long-long")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Fortran specifics
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  set (CMAKE_FLIB "gfortran")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  set (CMAKE_FLIB "g2c")
endif()

# Boost initialization
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.41 COMPONENTS serialization thread program_options regex unit_test_framework system filesystem REQUIRED)

# Find lapack
find_library(LAPACK_LIBRARY NAMES lapack )
if(NOT LAPACK_LIBRARY)
  message("error: could not find lapack library")
endif()

find_library(BLAS_LIBRARY NAMES blas)
if(NOT BLAS_LIBRARY)
  message("error: could not find blas library")
endif()

# Create all libraries in lib
set(LIBRARY_OUTPUT_PATH ${phaistos_BINARY_DIR}/libs)

# Setup regular expression patterns for enabling and disabling modules
# Only if not already set
if (NOT PHAISTOS_MODULE_ENABLE)
  # Check environment variable
  if (NOT $ENV{PHAISTOS_MODULE_ENABLE} STREQUAL "")
    set (PHAISTOS_MODULE_ENABLE $ENV{PHAISTOS_MODULE_ENABLE})
  # Set to default value
  else (NOT $ENV{PHAISTOS_MODULE_ENABLE} STREQUAL "")
    set (PHAISTOS_MODULE_ENABLE "")
  endif (NOT $ENV{PHAISTOS_MODULE_ENABLE} STREQUAL "")
endif (NOT PHAISTOS_MODULE_ENABLE)
if (NOT PHAISTOS_MODULE_DISABLE)
  # Check environment variable
  if (NOT $ENV{PHAISTOS_MODULE_DISABLE} STREQUAL "")
    set (PHAISTOS_MODULE_DISABLE $ENV{PHAISTOS_MODULE_DISABLE})
  # Set to default value
  else (NOT $ENV{PHAISTOS_MODULE_DISABLE} STREQUAL "")
    set (PHAISTOS_MODULE_DISABLE "^$")
  endif (NOT $ENV{PHAISTOS_MODULE_DISABLE} STREQUAL "")
endif (NOT PHAISTOS_MODULE_DISABLE)

# Setup regular expression patterns for enabling and disabling module binaries
# Only if not already set
if (NOT PHAISTOS_MODULE_BINARIES_ENABLE)
  # Check environment variable
  if (NOT $ENV{PHAISTOS_MODULE_BINARIES_ENABLE} STREQUAL "")
    set (PHAISTOS_MODULE_BINARIES_ENABLE $ENV{PHAISTOS_MODULE_BINARIES_ENABLE})
  # Set to default value
  else (NOT $ENV{PHAISTOS_MODULE_BINARIES_ENABLE} STREQUAL "")
    set (PHAISTOS_MODULE_BINARIES_ENABLE ".*")
  endif (NOT $ENV{PHAISTOS_MODULE_BINARIES_ENABLE} STREQUAL "")
endif (NOT PHAISTOS_MODULE_BINARIES_ENABLE)
if (NOT PHAISTOS_MODULE_BINARIES_DISABLE)
  # Check environment variable
  if (NOT $ENV{PHAISTOS_MODULE_BINARIES_DISABLE} STREQUAL "")
    set (PHAISTOS_MODULE_BINARIES_DISABLE $ENV{PHAISTOS_MODULE_BINARIES_DISABLE})
  # Set to default value
  else (NOT $ENV{PHAISTOS_MODULE_BINARIES_DISABLE} STREQUAL "")
    set (PHAISTOS_MODULE_BINARIES_DISABLE "^$")
  endif (NOT $ENV{PHAISTOS_MODULE_BINARIES_DISABLE} STREQUAL "")
endif (NOT PHAISTOS_MODULE_BINARIES_DISABLE)

# Copy data directory to build directory
# Reset phaistos_cpp directory
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)
INSTALL( DIRECTORY data DESTINATION . REGEX ".svn" EXCLUDE )

# Copy scripts directory to build directory
# Reset phaistos_cpp directory
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/scripts ${CMAKE_BINARY_DIR}/scripts)
INSTALL( DIRECTORY scripts DESTINATION . REGEX ".svn" EXCLUDE )

# Add current source directory and build directories to list of include directories
set(phaistos_MODULE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Construct svn revision header file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/revision.h.in ${CMAKE_CURRENT_BINARY_DIR}/revision.h.in @ONLY)

# for revision header
include(${phaistos_SOURCE_DIR}/cmake/modules/svnrev.cmake)

# Options
option(PHAISTOS_INCLUDE_CHAIN_CA "Whether to include support for CA chains" OFF)
if (PHAISTOS_INCLUDE_CHAIN_CA)
  add_definitions(-DPHAISTOS_INCLUDE_CHAIN_CA)
endif (PHAISTOS_INCLUDE_CHAIN_CA)
option(PHAISTOS_INCLUDE_OMEGA_DBNS "Whether to include DBN models with omega node" OFF)
if (PHAISTOS_INCLUDE_OMEGA_DBNS)
  add_definitions(-DPHAISTOS_INCLUDE_OMEGA_DBNS)
endif (PHAISTOS_INCLUDE_OMEGA_DBNS)

# Optionally include doc directory
option (PHAISTOS_INCLUDE_DOCS "Whether to exclude doc directory" OFF)

# Add subdirectories
add_subdirectory(modules ${phaistos_BINARY_DIR}/modules/.phase1)
add_subdirectory(src)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(bin)
if (PHAISTOS_INCLUDE_DOCS)
  add_subdirectory(doc)
endif (PHAISTOS_INCLUDE_DOCS)
# The module subdirectory is added twice - parsed in two phases
add_subdirectory(modules)


# Packaging support - generate package using "make package_source"
file(GLOB phaistos_BUILD_DIRS "build*")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${VERSION}${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
  "${phaistos_BUILD_DIRS};${CPACK_SOURCE_IGNORE_FILES};\\\\.svn")
install( FILES "${LICENSE_FILE}" ChangeLog DESTINATION doc )

##########################
# Binary packaging support
##########################
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Molecular modelling software package for simulation of proteins.")
set(CPACK_PACKAGE_VENDOR "University of Copenhagen")
set( CMAKE_INSTALL_PREFIX "/opt/${CMAKE_PROJECT_NAME}" )
set(CPACK_SET_DESTDIR  ON)
find_program(DPKGEXE "dpkg")
find_program(RPMEXE "rpmbuild")

if( APPLE )
  set( CPACK_GENERATOR "DragNDrop" )
  set( CPACK_DMG_FORMAT "UDBZ" )
  set( CPACK_DMG_VOLUME_NAME "Phaistos" )
  set( CPACK_SYSTEM_NAME "OSX" )
elseif(EXISTS ${DPKGEXE})
  message( STATUS "debian system detected")
  set(CPACK_GENERATOR "DEB")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Mikael Borg <mikael.borg@gmail.com>") #required
  set(CPACK_DEBIAN_PACKAGE_DEBUG "ON")
  set(CPACK_GENERATOR "DEB")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, libgcc1, libquadmath0, libgfortran3, libblas-dev, liblapack3, libstdc++6")
elseif(EXISTS ${RPMEXE})
  # note: rpmbuild includes dependencies automagically
  message( STATUS "redhat-based system detected")
  set(CPACK_RPM_PACKAGE_RELOCATABLE OFF)
  set(CPACK_GENERATOR "RPM")
  set(CPACK_RPM_PACKAGE_LICENSE "GPLv3")
  set(CPACK_RPM_PACKAGE_DESCRIPTION "Phaistos is a molecular modelling software
package for simulation of proteins.The distinctive features of Phaistos are
the use of the Markov chain Monte Carlo (MCMC) approach to molecular
simulation, and the use of sophisticated probabilistic models of protein
structure in continuous space. Phaistos can be used for molecular simulation
using empirical energy functions such as OPLS-AA or Profasi, or for
probabilistic, Bayesian inference of protein structure based on a likelihood
and conformational priors.")
else()
  set(CPACK_GENERATOR "TGZ")
endif()

set (CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${VERSION}-${CMAKE_SYSTEM_PROCESSOR}")
message( STATUS "Package name: ${CPACK_PACKAGE_FILE_NAME}")

include(CPack)
# Display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------- Phaistos Cmake settings ----------------------------" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}" )
MESSAGE( STATUS "LAPACK_LIBRARY = ${LAPACK_LIBRARY}" )
MESSAGE( STATUS "BLAS_LIBRARY = ${BLAS_LIBRARY}" )
MESSAGE( STATUS "PHAISTOS_INCLUDE_CHAIN_CA = ${PHAISTOS_INCLUDE_CHAIN_CA}" )
MESSAGE( STATUS "PHAISTOS_INCLUDE_OMEGA_DBNS = ${PHAISTOS_INCLUDE_OMEGA_DBNS}" )
MESSAGE( STATUS "PHAISTOS_INCLUDE_DOCS = ${PHAISTOS_INCLUDE_DOCS}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "For a list of all variables, use cmake -LA" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )
