
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/muninn/CMakeLists.txt")

  # Output compiled library in Phaistos build directory
  set(muninn_SOURCE_DIR "${CMAKE_BINARY_DIR}")

  # Add root directories to the list of includes
  include_directories (${CMAKE_CURRENT_SOURCE_DIR}/muninn ${CMAKE_CURRENT_SOURCE_DIR}/muninn/external) 

  add_subdirectory("muninn/muninn")

endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/muninn/CMakeLists.txt")

# Add scripts to scripts directory
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/muninn/scripts/ ${CMAKE_BINARY_DIR}/scripts/muninn/)
execute_process(COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_BINARY_DIR}/scripts/muninn/__init__.py)
