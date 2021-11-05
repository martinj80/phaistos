# Always update subversion revision
add_custom_target(svnrev ALL
                  COMMENT "Updating revision.h")
add_custom_command(TARGET svnrev
COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${phaistos_SOURCE_DIR}
                         -DINPUT=${phaistos_BINARY_DIR}/revision.h.in
                         -DOUTPUT=${phaistos_BINARY_DIR}/revision.h
                         -P ${phaistos_SOURCE_DIR}/cmake/modules/updateRevision.cmake)

set_source_files_properties(${phaistos_BINARY_DIR}/revision.h
    PROPERTIES GENERATED TRUE
    HEADER_FILE_ONLY TRUE)


