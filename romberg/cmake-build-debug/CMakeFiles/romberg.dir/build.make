# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /opt/clion-2018.3.4/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2018.3.4/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/artfin/Desktop/repos/NumerovProject/romberg

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/romberg.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/romberg.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/romberg.dir/flags.make

CMakeFiles/romberg.dir/main.cpp.o: CMakeFiles/romberg.dir/flags.make
CMakeFiles/romberg.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/romberg.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/romberg.dir/main.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/romberg/main.cpp

CMakeFiles/romberg.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/romberg.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/romberg/main.cpp > CMakeFiles/romberg.dir/main.cpp.i

CMakeFiles/romberg.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/romberg.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/romberg/main.cpp -o CMakeFiles/romberg.dir/main.cpp.s

# Object files for target romberg
romberg_OBJECTS = \
"CMakeFiles/romberg.dir/main.cpp.o"

# External object files for target romberg
romberg_EXTERNAL_OBJECTS =

romberg: CMakeFiles/romberg.dir/main.cpp.o
romberg: CMakeFiles/romberg.dir/build.make
romberg: CMakeFiles/romberg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable romberg"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/romberg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/romberg.dir/build: romberg

.PHONY : CMakeFiles/romberg.dir/build

CMakeFiles/romberg.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/romberg.dir/cmake_clean.cmake
.PHONY : CMakeFiles/romberg.dir/clean

CMakeFiles/romberg.dir/depend:
	cd /home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/artfin/Desktop/repos/NumerovProject/romberg /home/artfin/Desktop/repos/NumerovProject/romberg /home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug /home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug /home/artfin/Desktop/repos/NumerovProject/romberg/cmake-build-debug/CMakeFiles/romberg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/romberg.dir/depend

