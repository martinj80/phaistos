
function(require_module module_name)
  list(FIND phaistos_MODULE_NAMES_ACTIVE ${module_name} found_index)
  if (found_index LESS 0)
    message(SEND_ERROR "ERROR (require_module): Required module " ${module_name} " not found. Please enable or include module.")
  endif(found_index LESS 0)
endfunction(require_module)

# Determine whether to run phase 1 or 2
if (NOT phaistos_MODULE_PHASE )
  set(phaistos_MODULE_PHASE 1)
else (NOT phaistos_MODULE_PHASE)
  set(phaistos_MODULE_PHASE 2)
endif (NOT phaistos_MODULE_PHASE)
set(phaistos_MODULE_PHASE ${phaistos_MODULE_PHASE} PARENT_SCOPE)

# Report phase
message(STATUS "Modules - phase: " ${phaistos_MODULE_PHASE})

if (${phaistos_MODULE_PHASE} EQUAL 1)

  # Find all module directories
  file(GLOB phaistos_MODULE_FILES "${phaistos_SOURCE_DIR}/modules/phaistos_cpp")
  file(GLOB phaistos_MODULE_FILES_tmp "${phaistos_SOURCE_DIR}/modules/[a-zA-Z0-9]*")
  list(SORT phaistos_MODULE_FILES_tmp)
  list(APPEND phaistos_MODULE_FILES ${phaistos_MODULE_FILES_tmp})
  list(REMOVE_DUPLICATES phaistos_MODULE_FILES)
  foreach(file ${phaistos_MODULE_FILES})
    if (IS_DIRECTORY ${file})
      # Check if the directory is another modules directory
      string(REGEX MATCH "${phaistos_SOURCE_DIR}/modules/.*modules.*" submodules_dir ${file})
      if (submodules_dir)
        list(REMOVE_ITEM phaistos_MODULE_FILES ${file})
        file(GLOB phaistos_SUBMODULE_FILES "${file}/[a-zA-Z0-9]*")
        list(APPEND phaistos_MODULE_FILES ${phaistos_SUBMODULE_FILES})
        # message(STATUS "${phaistos_MODULE_FILES}")
      endif()
    endif()
  endforeach()  
  foreach(file ${phaistos_MODULE_FILES})
    if (IS_DIRECTORY ${file})
      file(RELATIVE_PATH relative_name ${phaistos_SOURCE_DIR}/modules ${file})
      list(APPEND phaistos_MODULE_DIRS ${relative_name})
    endif()
  endforeach()  

  message(STATUS "  Found modules:")
  foreach(module_dir ${phaistos_MODULE_DIRS})

    get_filename_component(module_name ${module_dir} NAME)

    # Check if module specifies any settings
    set(PHAISTOS_MODULE_DEFAULT_ENABLED 1)
    if (EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/settings.cmake)
      include(${phaistos_SOURCE_DIR}/modules/${module_dir}/settings.cmake)
    endif(EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/settings.cmake)

    # Define a local version of the PHAISTOS_MODULE_ENABLE regex. In case
    # the global version is not explicitly set by the user (or set to ""), modules can
    # be configured to be turned off by default
    set(PHAISTOS_CURRENT_MODULE_ENABLE ${PHAISTOS_MODULE_ENABLE})
    if ("${PHAISTOS_CURRENT_MODULE_ENABLE}" STREQUAL "")
      if (${PHAISTOS_MODULE_DEFAULT_ENABLED})
        set(PHAISTOS_CURRENT_MODULE_ENABLE ".*")
      endif (${PHAISTOS_MODULE_DEFAULT_ENABLED})
    endif ("${PHAISTOS_CURRENT_MODULE_ENABLE}" STREQUAL "")

    if (NOT ${PHAISTOS_CURRENT_MODULE_ENABLE} STREQUAL "")
      string(REGEX MATCH ${PHAISTOS_CURRENT_MODULE_ENABLE}  match_enabled   ${module_name})
      string(REGEX MATCH ${PHAISTOS_MODULE_DISABLE} match_disabled  ${module_name})
      if (match_enabled)
        if (${match_enabled} STREQUAL ${module_name})
          set(enabled true)
        else ()
          set(enabled false)
        endif()
      else()
        set(enabled false)
      endif()
      if (match_disabled)
        if (${match_disabled} STREQUAL ${module_name})
          set(disabled true)
        else ()
          set(disabled false)
        endif()
      else()
        set(disabled false)
      endif()

      if (enabled AND NOT disabled)
        message(STATUS "    ${module_name}")
        list(APPEND phaistos_MODULE_DIRS_ACTIVE ${module_dir})
        list(APPEND phaistos_MODULE_NAMES_ACTIVE ${module_name})
      else()
        message(STATUS "    ${module_name} (disabled)")
      endif()
    else (NOT ${PHAISTOS_CURRENT_MODULE_ENABLE} STREQUAL "")
        message(STATUS "    ${module_name} (disabled by default)")      
    endif (NOT ${PHAISTOS_CURRENT_MODULE_ENABLE} STREQUAL "")

  endforeach(module_dir)

  # Set variable in parent scope (so that it is available for phase 2)
  set(phaistos_MODULE_DIRS_ACTIVE ${phaistos_MODULE_DIRS_ACTIVE} PARENT_SCOPE)
  set(phaistos_MODULE_NAMES_ACTIVE ${phaistos_MODULE_NAMES_ACTIVE} PARENT_SCOPE)

  message(STATUS "  Adding sources to libphaistos...")
  foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})

    # Iterate over all directories containing sources for libphaistos (currently only src)
    foreach(module_subdir src)

      if (EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/${module_subdir})

        # subdirectory
        set(_CMakeLists_file ${phaistos_SOURCE_DIR}/modules/${module_dir}/${module_subdir}/CMakeLists.txt)
        if (EXISTS ${_CMakeLists_file})

          # Save current SRC variable
          set(_phaistos_MODULE_SRC_TMP ${phaistos_MODULE_SRC})
          set(phaistos_MODULE_SRC "")

          set(CMAKE_CURRENT_SOURCE_DIR_BACKUP ${CMAKE_CURRENT_SOURCE_DIR})
          set(CMAKE_CURRENT_SOURCE_DIR ${phaistos_SOURCE_DIR}/modules/${module_dir}/${module_subdir}/)
          set(CMAKE_CURRENT_BINARY_DIR_BACKUP ${CMAKE_CURRENT_BINARY_DIR})
          set(CMAKE_CURRENT_BINARY_DIR ${phaistos_BINARY_DIR}/modules/${module_dir}/${module_subdir}/)
          include(${_CMakeLists_file})
          set(CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR_BACKUP})
          set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR_BACKUP})

          # Save new SRC variable
          set(_phaistos_MODULE_SRC_NEW ${phaistos_MODULE_SRC})
          set(phaistos_MODULE_SRC ${_phaistos_MODULE_SRC_TMP})

          # Prepend directory to all source files
          foreach(src_file ${_phaistos_MODULE_SRC_NEW})
            list(APPEND phaistos_MODULE_SRC "${phaistos_SOURCE_DIR}/modules/${module_dir}/${module_subdir}/${src_file}")
          endforeach()

          # Restore backed up version
          set(phaistos_MODULE_SRC ${_phaistos_MODULE_SRC_tmp} ${phaistos_MODULE_SRC})

        endif(EXISTS ${_CMakeLists_file})

        # Add current source directory to list of include directories
        set(phaistos_MODULE_INCLUDE_DIRS ${phaistos_MODULE_INCLUDE_DIRS} ${phaistos_SOURCE_DIR}/modules/${module_dir}/${module_subdir}/)
        set(phaistos_MODULE_INCLUDE_DIRS ${phaistos_MODULE_INCLUDE_DIRS} PARENT_SCOPE)
        
        # Add module libraries to list
        set(phaistos_MODULE_LIBS ${phaistos_MODULE_LIBS} PARENT_SCOPE)

      endif(EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/${module_subdir})

    endforeach(module_subdir)
  endforeach(module_dir)



  message(STATUS "  Preparing phaistos.cpp includes...")

  # Reset phaistos_cpp directory
  file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/modules/phaistos_cpp)

  # Extract all possible module includes from bin directory
  set(phaistos_cpp_matches "")
  set(phaistos_BIN_FILES "")
  file(GLOB phaistos_BIN_FILES_tmp "${phaistos_SOURCE_DIR}/bin/*.cpp")
  list(APPEND phaistos_BIN_FILES ${phaistos_BIN_FILES_tmp})
  file(GLOB phaistos_BIN_FILES_tmp "${phaistos_SOURCE_DIR}/bin/*.h")
  list(APPEND phaistos_BIN_FILES ${phaistos_BIN_FILES_tmp})
  file(GLOB phaistos_BIN_FILES_tmp "${phaistos_SOURCE_DIR}/modules/phaistos_cpp/bin/*.cpp")
  list(APPEND phaistos_BIN_FILES ${phaistos_BIN_FILES_tmp})
  file(GLOB phaistos_BIN_FILES_tmp "${phaistos_SOURCE_DIR}/modules/phaistos_cpp/bin/*.h")
  list(APPEND phaistos_BIN_FILES ${phaistos_BIN_FILES_tmp})
  foreach(file ${phaistos_BIN_FILES})

    # Extract all possible includes from phaistos.cpp
    file(READ ${file} phaistos_bin_file)
    string(REGEX MATCHALL "\"modules/phaistos_cpp/[^\"]*\"" phaistos_cpp_matches_tmp ${phaistos_bin_file})

    list(APPEND phaistos_cpp_matches ${phaistos_cpp_matches_tmp})
  endforeach()

  # Generate empty includes
  foreach (phaistos_cpp_include ${phaistos_cpp_matches})
    string(LENGTH ${phaistos_cpp_include} phaistos_cpp_include_length)
    math(EXPR phaistos_cpp_include_new_length "${phaistos_cpp_include_length}-2")
    string(SUBSTRING ${phaistos_cpp_include} 1 ${phaistos_cpp_include_new_length} phaistos_cpp_include_short)
    file(WRITE ${CMAKE_BINARY_DIR}/${phaistos_cpp_include_short} "")
  endforeach()

  foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
    # phaistos_cpp directory
    set(_phaistos_cpp_directory ${phaistos_SOURCE_DIR}/modules/${module_dir}/bin/phaistos_cpp)
    if (EXISTS ${_phaistos_cpp_directory})

      file(GLOB _phaistos_cpp_files ${_phaistos_cpp_directory}/*.cpp)
      foreach(_phaistos_cpp_file ${_phaistos_cpp_files})
        file(RELATIVE_PATH relative_name ${_phaistos_cpp_directory} ${_phaistos_cpp_file})
        set(_phaistos_cpp_file_merged ${CMAKE_BINARY_DIR}/modules/phaistos_cpp/${relative_name})
        file(APPEND ${_phaistos_cpp_file_merged} "#include \"${_phaistos_cpp_file}\"\n") 
      endforeach(_phaistos_cpp_file ${_phaistos_cpp_files})
      
    endif (EXISTS ${_phaistos_cpp_directory})
  endforeach(module_dir)

  # Communicate variable to parent scope
  set (phaistos_MODULE_SRC ${phaistos_MODULE_SRC} PARENT_SCOPE)


  message(STATUS "  Copying data directories...")
  foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
    # data directory
    set(_data_dir ${phaistos_SOURCE_DIR}/modules/${module_dir}/data)
    if (EXISTS ${_data_dir})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${_data_dir} ${CMAKE_BINARY_DIR}/data)
      message(STATUS "    ${_data_dir}")
    endif (EXISTS ${_data_dir})
  endforeach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})

  message(STATUS "  Copying script directories...")
  foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
    # scripts directory
    set(_scripts_dir ${phaistos_SOURCE_DIR}/modules/${module_dir}/scripts)
    if (EXISTS ${_scripts_dir})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${_scripts_dir} ${CMAKE_BINARY_DIR}/scripts)
      message(STATUS "    ${_scripts_dir}")
    endif (EXISTS ${_scripts_dir})
  endforeach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})



  message(STATUS "  Preparing doc includes...")

  # Define authors file name
  set(_authors_file_merged ${CMAKE_BINARY_DIR}/modules/doc/authors.txt)

  # Reset authors file
  file(WRITE ${_authors_file_merged} "")


  # Extract all possible includes from main manual.tex
  file(READ ${CMAKE_SOURCE_DIR}/doc/manual.tex manual_tex_source)
  string(REGEX MATCHALL "\\.\\./modules/doc/[.0-9a-zA-Z_-]*" manual_tex_matches ${manual_tex_source})

  # Generate empty includes
  foreach (manual_tex_include ${manual_tex_matches})
    string(REPLACE "../modules/" "" short_name ${manual_tex_include})
    file(WRITE ${CMAKE_BINARY_DIR}/modules/${short_name} "")
  endforeach (manual_tex_include ${manual_tex_matches})



  # # Extract all possible includes from main manual.texi
  # file(READ ${CMAKE_SOURCE_DIR}/doc/manual.texi manual_texi_source)
  # string(REGEX MATCHALL "\\.\\./modules/doc/[.0-9a-zA-Z_-]*" manual_texi_matches ${manual_texi_source})

  # # Generate empty includes
  # foreach (manual_texi_include ${manual_texi_matches})
  #   string(REPLACE "../modules/" "" short_name ${manual_texi_include})
  #   file(WRITE ${CMAKE_BINARY_DIR}/modules/${short_name} "")
  # endforeach (manual_texi_include ${manual_texi_matches})

  # Include tex files in modules
  set (phaistos_TEX_FILES "")
  set (phaistos_BIB_FILES "")
  foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
    # doc directory
    set(_doc_directory ${phaistos_SOURCE_DIR}/modules/${module_dir}/doc)
    if (EXISTS ${_doc_directory})

      # Find tex files in module
      file(GLOB _manual_files ${_doc_directory}/*.tex)
      foreach(_tex_file ${_manual_files})
        file(RELATIVE_PATH relative_name ${_doc_directory} ${_tex_file})
        file(RELATIVE_PATH relative_name_docs "${CMAKE_SOURCE_DIR}/doc" ${_tex_file})
        if (EXISTS "${CMAKE_BINARY_DIR}/modules/doc/${relative_name}")
          message(STATUS "    Found: ${_tex_file}")
          file(APPEND ${CMAKE_BINARY_DIR}/modules/doc/${relative_name} "\\input{${_tex_file}}\n") 
          set(phaistos_TEX_FILES "${phaistos_TEX_FILES};${relative_name_docs}")
        endif (EXISTS "${CMAKE_BINARY_DIR}/modules/doc/${relative_name}")          
      endforeach(_tex_file ${_manual_files})

      # Find bib files in module
      file(GLOB _manual_bib_files ${_doc_directory}/*.bib)
      foreach(_bib_file ${_manual_bib_files})
        file(RELATIVE_PATH relative_name_docs "${CMAKE_SOURCE_DIR}/doc" ${_bib_file})
        set(phaistos_BIB_FILES "${phaistos_BIB_FILES};${relative_name_docs}")
      endforeach(_bib_file ${_manual_bib_files})

    endif(EXISTS ${_doc_directory})

    # Find authors file
    set(_authors_file ${phaistos_SOURCE_DIR}/modules/${module_dir}/doc/authors.txt)
    if (EXISTS ${_authors_file})

      file(READ ${_authors_file} _authors_file_contents)
      file(APPEND ${_authors_file_merged} ${_authors_file_contents})

    endif (EXISTS ${_authors_file})
    
  endforeach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})

  # Set variable in parent scope (so that it is available for phase 2)
  set(phaistos_TEX_FILES ${phaistos_TEX_FILES} PARENT_SCOPE)
  set(phaistos_BIB_FILES ${phaistos_BIB_FILES} PARENT_SCOPE)



  # # Extract all possible includes from main manual.texi
  # file(READ ${CMAKE_SOURCE_DIR}/doc/manual.texi manual_texi_source)
  # string(REGEX MATCHALL "\\.\\./modules/doc/[.0-9a-zA-Z_-]*" manual_texi_matches ${manual_texi_source})

  # # Generate empty includes
  # foreach (manual_texi_include ${manual_texi_matches})
  #   string(REPLACE "../modules/" "" short_name ${manual_texi_include})
  #   file(WRITE ${CMAKE_BINARY_DIR}/modules/${short_name} "")
  # endforeach (manual_texi_include ${manual_texi_matches})

  # # Include texi files in modules
  # foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
  #   # doc directory
  #   set(_doc_directory ${phaistos_SOURCE_DIR}/modules/${module_dir}/doc)
  #   if (EXISTS ${_doc_directory})

  #     # Find texi files in module
  #     file(GLOB _manual_files ${_doc_directory}/*.texi)
  #     foreach(_texi_file ${_manual_files})
  #       file(RELATIVE_PATH relative_name ${_doc_directory} ${_texi_file})
  #       if (EXISTS "${CMAKE_BINARY_DIR}/doc/${relative_name}")
  #         message(STATUS "    Found: ${_texi_file}")
  #         file(APPEND ${CMAKE_BINARY_DIR}/modules/doc/${relative_name} "@include ${_texi_file}\n") 
  #       endif (EXISTS "${CMAKE_BINARY_DIR}/doc/${relative_name}")          
  #     endforeach(_texi_file ${_manual_files})
  #   endif(EXISTS ${_doc_directory})

  #   # Find authors file
  #   set(_authors_file ${phaistos_SOURCE_DIR}/modules/${module_dir}/doc/authors.txt)
  #   if (EXISTS ${_authors_file})

  #     file(READ ${_authors_file} _authors_file_contents)
  #     file(APPEND ${_authors_file_merged} ${_authors_file_contents})

  #   endif (EXISTS ${_authors_file})
    
  # endforeach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})

  # set(_doc_file_merged ${CMAKE_BINARY_DIR}/modules/doc/manual.texi)
  # set(_authors_file_merged ${CMAKE_BINARY_DIR}/modules/doc/authors.txt)

  # # Reset authors file
  # file(WRITE ${_authors_file_merged} "")
  # foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
  #   # doc directory
  #   set(_doc_file ${phaistos_SOURCE_DIR}/modules/${module_dir}/doc/manual.texi)
  #   if (EXISTS ${_doc_file})

  #     message(STATUS "    Found: ${_doc_file}")
  #     file(APPEND ${_doc_file_merged} "@include ${_doc_file}\n") 

  #   endif (EXISTS ${_doc_file})

  #   set(_authors_file ${phaistos_SOURCE_DIR}/modules/${module_dir}/doc/authors.txt)
  #   if (EXISTS ${_authors_file})

  #     file(READ ${_authors_file} _authors_file_contents)
  #     file(APPEND ${_authors_file_merged} ${_authors_file_contents})

  #   endif (EXISTS ${_authors_file})

  # endforeach(module_dir)

  # Sorty entries in author file
  file(STRINGS ${_authors_file_merged} _authors)
  set(_authors_names_extended "")
  foreach(_author ${_authors})
    string(REPLACE " " ";" _author ${_author})
    list(LENGTH _author _author_names_length)
    math(EXPR _last_name_index "${_author_names_length}-1")
    list(GET _author ${_last_name_index} _last_name)
    set(_extended_name ${_last_name} ${_author})
    string(REPLACE ";" " " _extended_name "${_extended_name}")
    list(APPEND _authors_names_extended ${_extended_name})
  endforeach()
  list(SORT _authors_names_extended)
  set(_${_authors_names_sorted} "")
  foreach(_author ${_authors_names_extended})
    string(REPLACE " " ";" _author ${_author})
    list(REMOVE_AT _author 0) 
    string(REPLACE ";" " " _author "${_author}")
    list(APPEND _authors_names_sorted ${_author})
  endforeach()
  string(REPLACE ";" "\n" _authors_names_sorted "${_authors_names_sorted}")
  file(WRITE ${_authors_file_merged} ${_authors_names_sorted})
  # message(STATUS "${_authors_names_sorted} ${_authors_file_merged}")


# Phase 2: add modules as subdirectories. The CMakeLists.txt in the root directory of 
# the module is responsible for setting up the correct building environment 
elseif (${phaistos_MODULE_PHASE} EQUAL 2)

  # Add all subdirectories
  foreach(module_dir ${phaistos_MODULE_DIRS_ACTIVE})
    if (EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/CMakeLists.txt)
      
      string(REGEX MATCH ${PHAISTOS_MODULE_BINARIES_ENABLE}  match_enabled   ${module_dir})
      string(REGEX MATCH ${PHAISTOS_MODULE_BINARIES_DISABLE} match_disabled  ${module_dir})
      if (match_enabled)
        if (${match_enabled} STREQUAL ${module_dir})
          set(enabled true)
        else ()
          set(enabled false)
        endif()
      else()
        set(enabled false)
      endif()
      if (match_disabled)
        if (${match_disabled} STREQUAL ${module_dir})
          set(disabled true)
        else ()
          set(disabled false)
        endif()
      else()
        set(disabled false)
      endif()

      if (enabled AND NOT disabled)
        message(STATUS "  Preparing binaries in ${module_dir}...")
        add_subdirectory(${module_dir})
      else (enabled AND NOT disabled)
        message(STATUS "  Preparing binaries in ${module_dir}...SKIPPED (binaries disabled)")
      endif (enabled AND NOT disabled)
    else (EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/CMakeLists.txt)
      message(STATUS "  Preparing binaries in ${module_dir}...SKIPPED (no CMakeLists.txt file)")
    endif (EXISTS ${phaistos_SOURCE_DIR}/modules/${module_dir}/CMakeLists.txt)
  endforeach()

endif (${phaistos_MODULE_PHASE} EQUAL 1)



