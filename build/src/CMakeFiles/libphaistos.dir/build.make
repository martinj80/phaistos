# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/juhasm/Documents/phaistos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/juhasm/Documents/phaistos/build

# Include any dependencies generated for this target.
include src/CMakeFiles/libphaistos.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/libphaistos.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/libphaistos.dir/flags.make

src/CMakeFiles/libphaistos.dir/protein/atom.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/protein/atom.cpp.o: ../src/protein/atom.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/libphaistos.dir/protein/atom.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/protein/atom.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/protein/atom.cpp

src/CMakeFiles/libphaistos.dir/protein/atom.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/protein/atom.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/protein/atom.cpp > CMakeFiles/libphaistos.dir/protein/atom.cpp.i

src/CMakeFiles/libphaistos.dir/protein/atom.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/protein/atom.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/protein/atom.cpp -o CMakeFiles/libphaistos.dir/protein/atom.cpp.s

src/CMakeFiles/libphaistos.dir/protein/residue.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/protein/residue.cpp.o: ../src/protein/residue.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/libphaistos.dir/protein/residue.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/protein/residue.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/protein/residue.cpp

src/CMakeFiles/libphaistos.dir/protein/residue.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/protein/residue.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/protein/residue.cpp > CMakeFiles/libphaistos.dir/protein/residue.cpp.i

src/CMakeFiles/libphaistos.dir/protein/residue.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/protein/residue.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/protein/residue.cpp -o CMakeFiles/libphaistos.dir/protein/residue.cpp.s

src/CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.o: ../src/protein/pdb_input.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/protein/pdb_input.cpp

src/CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/protein/pdb_input.cpp > CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.i

src/CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/protein/pdb_input.cpp -o CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.s

src/CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.o: ../src/utils/eigen_system_3x3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/utils/eigen_system_3x3.cpp

src/CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/utils/eigen_system_3x3.cpp > CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.i

src/CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/utils/eigen_system_3x3.cpp -o CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.s

src/CMakeFiles/libphaistos.dir/utils/math.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/utils/math.cpp.o: ../src/utils/math.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/libphaistos.dir/utils/math.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/utils/math.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/utils/math.cpp

src/CMakeFiles/libphaistos.dir/utils/math.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/utils/math.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/utils/math.cpp > CMakeFiles/libphaistos.dir/utils/math.cpp.i

src/CMakeFiles/libphaistos.dir/utils/math.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/utils/math.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/utils/math.cpp -o CMakeFiles/libphaistos.dir/utils/math.cpp.s

src/CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.o: ../src/utils/program_option_parser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/utils/program_option_parser.cpp

src/CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/utils/program_option_parser.cpp > CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.i

src/CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/utils/program_option_parser.cpp -o CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.s

src/CMakeFiles/libphaistos.dir/utils/svd.cpp.o: src/CMakeFiles/libphaistos.dir/flags.make
src/CMakeFiles/libphaistos.dir/utils/svd.cpp.o: ../src/utils/svd.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/CMakeFiles/libphaistos.dir/utils/svd.cpp.o"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libphaistos.dir/utils/svd.cpp.o -c /mnt/c/Users/juhasm/Documents/phaistos/src/utils/svd.cpp

src/CMakeFiles/libphaistos.dir/utils/svd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libphaistos.dir/utils/svd.cpp.i"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/juhasm/Documents/phaistos/src/utils/svd.cpp > CMakeFiles/libphaistos.dir/utils/svd.cpp.i

src/CMakeFiles/libphaistos.dir/utils/svd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libphaistos.dir/utils/svd.cpp.s"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && /usr/bin/g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/juhasm/Documents/phaistos/src/utils/svd.cpp -o CMakeFiles/libphaistos.dir/utils/svd.cpp.s

# Object files for target libphaistos
libphaistos_OBJECTS = \
"CMakeFiles/libphaistos.dir/protein/atom.cpp.o" \
"CMakeFiles/libphaistos.dir/protein/residue.cpp.o" \
"CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.o" \
"CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.o" \
"CMakeFiles/libphaistos.dir/utils/math.cpp.o" \
"CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.o" \
"CMakeFiles/libphaistos.dir/utils/svd.cpp.o"

# External object files for target libphaistos
libphaistos_EXTERNAL_OBJECTS =

libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/protein/atom.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/protein/residue.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/protein/pdb_input.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/utils/eigen_system_3x3.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/utils/math.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/utils/program_option_parser.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/utils/svd.cpp.o
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/build.make
libs/libphaistos.a: src/CMakeFiles/libphaistos.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/juhasm/Documents/phaistos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX static library ../libs/libphaistos.a"
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && $(CMAKE_COMMAND) -P CMakeFiles/libphaistos.dir/cmake_clean_target.cmake
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libphaistos.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/libphaistos.dir/build: libs/libphaistos.a

.PHONY : src/CMakeFiles/libphaistos.dir/build

src/CMakeFiles/libphaistos.dir/clean:
	cd /mnt/c/Users/juhasm/Documents/phaistos/build/src && $(CMAKE_COMMAND) -P CMakeFiles/libphaistos.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/libphaistos.dir/clean

src/CMakeFiles/libphaistos.dir/depend:
	cd /mnt/c/Users/juhasm/Documents/phaistos/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/juhasm/Documents/phaistos /mnt/c/Users/juhasm/Documents/phaistos/src /mnt/c/Users/juhasm/Documents/phaistos/build /mnt/c/Users/juhasm/Documents/phaistos/build/src /mnt/c/Users/juhasm/Documents/phaistos/build/src/CMakeFiles/libphaistos.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/libphaistos.dir/depend
