require_module("sidechain_dbns")
require_module("crisp")
require_module("git")
require_module("opls")

# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/src/models/backbone_dbn ${CMAKE_CURRENT_SOURCE_DIR}/../src ${MOCAPY_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${phaistos_MODULE_INCLUDE_DIRS} ${MUNINN_INCLUDE_DIRS}) 

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Create executable: typhon
add_executable(typhon typhon.cpp typhon_options.cpp)
add_dependencies(typhon svnrev)

# Link executable with library
target_link_libraries(typhon libphaistos dbn integrate ${MUNINN_LIBRARIES} ${LAPACK_LIBRARY} ${BLAS_LIBRARY} libbasilisk libcompas git ${MOCAPY_LIBRARIES} ${CMAKE_FLIB} ${Boost_LIBRARIES})

INSTALL(TARGETS typhon DESTINATION bin )
