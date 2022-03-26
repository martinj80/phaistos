# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/home/martinjuhas/CLionProjects/phaistos;/home/martinjuhas/CLionProjects/phaistos/build-debug")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEBIAN_PACKAGE_DEBUG "ON")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, libgcc1, libquadmath0, libgfortran3, libblas-dev, liblapack3, libstdc++6")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Mikael Borg <mikael.borg@gmail.com>")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/snap/clion/184/bin/cmake/linux/share/cmake-3.21/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "phaistos built using CMake")
set(CPACK_GENERATOR "DEB")
set(CPACK_INSTALL_CMAKE_PROJECTS "/home/martinjuhas/CLionProjects/phaistos/build-debug;phaistos;ALL;/")
set(CPACK_INSTALL_PREFIX "/opt/phaistos")
set(CPACK_MODULE_PATH "/home/martinjuhas/CLionProjects/phaistos/cmake/modules/")
set(CPACK_NSIS_DISPLAY_NAME "phaistos 0.1.1")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "phaistos 0.1.1")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OUTPUT_CONFIG_FILE "/home/martinjuhas/CLionProjects/phaistos/build-debug/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/snap/clion/184/bin/cmake/linux/share/cmake-3.21/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Molecular modelling software package for simulation of proteins.")
set(CPACK_PACKAGE_FILE_NAME "phaistos-1.0-x86_64")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "phaistos 0.1.1")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "phaistos 0.1.1")
set(CPACK_PACKAGE_NAME "phaistos")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "University of Copenhagen")
set(CPACK_PACKAGE_VERSION "0.1.1")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_RESOURCE_FILE_LICENSE "/snap/clion/184/bin/cmake/linux/share/cmake-3.21/Templates/CPack.GenericLicense.txt")
set(CPACK_RESOURCE_FILE_README "/snap/clion/184/bin/cmake/linux/share/cmake-3.21/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/snap/clion/184/bin/cmake/linux/share/cmake-3.21/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "ON")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/home/martinjuhas/CLionProjects/phaistos/build;/home/martinjuhas/CLionProjects/phaistos/build-debug;/home/martinjuhas/CLionProjects/phaistos/build-release;;\\.svn")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/martinjuhas/CLionProjects/phaistos/build-debug/CPackSourceConfig.cmake")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "phaistos-1.0")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Linux")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/martinjuhas/CLionProjects/phaistos/build-debug/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
