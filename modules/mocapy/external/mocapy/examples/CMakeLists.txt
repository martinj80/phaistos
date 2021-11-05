SET(BUILD_RELEASE ON)
SET(BUILD_DEBUG OFF)
SET(BUILD_PROFILE OFF)
SET(BUILD_DUMMYRANDOM OFF)
SET(CMAKE_VERBOSE_MAKEFILE OFF) 

SET(Boost_USE_STATIC_LIBS ON)


FIND_LIBRARY(PTHREAD pthread)

# Make sure the compiler can find include files from mocapy and boost.
INCLUDE_DIRECTORIES(${MOCAPY_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})


IF(BUILD_RELEASE)
SET(CMAKE_CXX_FLAGS "-O3 -fPIC -Wall -fmessage-length=0")
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


SET(PROGS mdarray hmm_bippo hmm_factorial hmm_simple hmm_discrete discrete_hmm_with_prior hmm_gauss hmm_gauss_1d hmm_multinomial hmm_dirichlet hmm_poisson hmm_vonmises infenginehmm_example infenginemm_example hmm_torus hmm_kent print_dbn train_example)
SET(PROGS_BOOST_OPTIONS_THREADS mlr-uni)

FOREACH(p ${PROGS})
        add_executable(${p} ${p}.cpp)
        target_link_libraries (${p} ${MOCAPYLIB} ${Boost_SERIALIZATION_LIBRARY} ${LAPACK_LIBRARY} ${CMAKE_FLIB})
ENDFOREACH(p)

FOREACH(p ${PROGS_BOOST_OPTIONS_THREADS})
        add_executable(${p} ${p}.cpp)
        target_link_libraries (${p} ${MOCAPYLIB} ${Boost_SERIALIZATION_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_THREAD_LIBRARY} ${PTHREAD} ${LAPACK_LIBRARY} ${CMAKE_FLIB})
ENDFOREACH(p)

# copy data directory
exec_program(
        ${CMAKE_COMMAND}
        ARGS
        "-E create_symlink \"${CMAKE_CURRENT_SOURCE_DIR}/data\" \"${CMAKE_CURRENT_BINARY_DIR}/data\""
        )
#exec_program(
#        ${CMAKE_COMMAND}
#        ARGS
#        "-E copy_directory \"${CMAKE_CURRENT_SOURCE_DIR}/data\" \"${CMAKE_CURRENT_BINARY_DIR}\""
#        )
