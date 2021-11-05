# Use internal copy of Mocapy if MUNINN_ROOT is not set
if (NOT MUNINN_ROOT OR MUNINN_ROOT STREQUAL "")

  set(MUNINN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../external/muninn/")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../external)

  # Find Muninn
  find_package(MUNINN)
  set(MUNINN_LIBRARIES "muninn" PARENT_SCOPE)

else (NOT MUNINN_ROOT OR MUNINN_ROOT STREQUAL "")

  # Find Muninn
  find_package(MUNINN REQUIRED)

endif (NOT MUNINN_ROOT OR MUNINN_ROOT STREQUAL "")

set_property(DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY COMPILE_DEFINITIONS ${MUNINN_DEFINITIONS})

# Add source files
list(APPEND 
  phaistos_MODULE_SRC
  monte_carlo/monte_carlo_muninn.h
)
