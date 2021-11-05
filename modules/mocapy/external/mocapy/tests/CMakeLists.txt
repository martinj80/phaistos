SET(BUILD_RELEASE ON)
SET(BUILD_DEBUG OFF)
SET(BUILD_PROFILE OFF)
SET(BUILD_DUMMYRANDOM OFF)

# Make sure the compiler can find include files from mocapy.
include_directories (${MOCAPY_SOURCE_DIR}/src)


INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ${MOCAPY_SOURCE_DIR}/src)

IF(BUILD_RELEASE)
SET(CMAKE_CXX_FLAGS "-O3 -Wall -fmessage-length=0")
SET(MOCAPYLIB Mocapy CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_RELEASE)

IF(BUILD_DEBUG)
SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fmessage-length=0")
SET(MOCAPYLIB MocapyD CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_DEBUG)

IF(BUILD_PROFILE)
SET(CMAKE_CXX_FLAGS "-pg -O3 -Wall -fmessage-length=0")
SET(MOCAPYLIB MocapyP CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_PROFILE)

IF(BUILD_DUMMYRANDOM)
SET(CMAKE_CXX_FLAGS "-DRANDDUMMY -O3 -Wall -fmessage-length=0")
SET(MOCAPYLIB MocapyDR CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_DUMMYRANDOM)

SET(PROGS infenginehmm_test infenginemm_test test_mdarray test_mdarray_permute)

FOREACH(p ${PROGS})
        add_executable(${p} ${p}.cpp)
        target_link_libraries (${p} ${MOCAPYLIB} ${Boost_SERIALIZATION_LIBRARY} ${LAPACK_LIBRARY} ${CMAKE_FLIB})
ENDFOREACH(p)
