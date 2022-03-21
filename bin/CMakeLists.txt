# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${phaistos_SOURCE_DIR}/src/models/backbone_dbn) 

# Make sure the linker can find the libphaistos library once it is built. 
link_directories (${phaistos_BINARY_DIR}/libs)

# Find pthreads
find_library(PTHREAD pthread)
