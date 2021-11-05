
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mocapy/CMakeLists.txt")
  # Output compiled library in Phaistos build directory
  set(MOCAPY_SOURCE_DIR "${CMAKE_BINARY_DIR}")

  add_subdirectory("mocapy/src")

endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mocapy/CMakeLists.txt")