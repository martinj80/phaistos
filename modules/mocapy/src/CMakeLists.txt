# Use internal copy of Mocapy if MOCAPY_ROOT is not set
if (NOT MOCAPY_ROOT OR MOCAPY_ROOT STREQUAL "")

  set(MOCAPY_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../external/mocapy/")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../external)

  # Find Mocapy
  find_package(MOCAPY)
  set(MOCAPY_LIBRARIES "Mocapy" PARENT_SCOPE)

  
else (NOT MOCAPY_ROOT OR MOCAPY_ROOT STREQUAL "")

  # Find Mocapy
  find_package(MOCAPY REQUIRED)

endif (NOT MOCAPY_ROOT OR MOCAPY_ROOT STREQUAL "")





