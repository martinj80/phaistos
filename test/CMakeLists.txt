# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/src/models/backbone_dbn) 

# Make sure the linker can find the libphaistos library once it is built. 
link_directories (${phaistos_BINARY_DIR}/libs)

# Find pthreads
find_library(PTHREAD pthread)

# Create executable
add_executable(test_optimization test_optimization.cpp)

# Link executable with library
target_link_libraries(test_optimization libphaistos)


