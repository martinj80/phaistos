# This is some hoopla in order to always have updated subversion revision
# number. See:
# http://stackoverflow.com/questions/3780667/use-cmake-to-get-build-time-svn-revision
# http://www.cmake.org/pipermail/cmake/2010-July/038015.html

find_package(Subversion)
if(Subversion_FOUND)
  if(EXISTS "${SOURCE_DIR}/.svn/")
    Subversion_WC_INFO(${SOURCE_DIR} PHAISTOS)
    set(SUBVERSION_REVISION ${PHAISTOS_WC_REVISION})
  else(EXISTS "${SOURCE_DIR}/.svn/")
    set(SUBVERSION_REVISION "N/A")
  endif(EXISTS "${SOURCE_DIR}/.svn/")
endif(Subversion_FOUND)
CONFIGURE_FILE(${INPUT} ${OUTPUT})




