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
CMAKE_SOURCE_DIR = /home/artfin/Desktop/repos/NumerovProject

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/artfin/Desktop/repos/NumerovProject/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/matrix_numerov.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/matrix_numerov.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/matrix_numerov.dir/flags.make

CMakeFiles/matrix_numerov.dir/main.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/matrix_numerov.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/main.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/main.cpp

CMakeFiles/matrix_numerov.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/main.cpp > CMakeFiles/matrix_numerov.dir/main.cpp.i

CMakeFiles/matrix_numerov.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/main.cpp -o CMakeFiles/matrix_numerov.dir/main.cpp.s

CMakeFiles/matrix_numerov.dir/filereader.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/filereader.cpp.o: ../filereader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/matrix_numerov.dir/filereader.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/filereader.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/filereader.cpp

CMakeFiles/matrix_numerov.dir/filereader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/filereader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/filereader.cpp > CMakeFiles/matrix_numerov.dir/filereader.cpp.i

CMakeFiles/matrix_numerov.dir/filereader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/filereader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/filereader.cpp -o CMakeFiles/matrix_numerov.dir/filereader.cpp.s

CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.o: ../generalizedmatrixnumerov.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/generalizedmatrixnumerov.cpp

CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/generalizedmatrixnumerov.cpp > CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.i

CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/generalizedmatrixnumerov.cpp -o CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.s

CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.o: ../matrixnumerov.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/matrixnumerov.cpp

CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/matrixnumerov.cpp > CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.i

CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/matrixnumerov.cpp -o CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.s

CMakeFiles/matrix_numerov.dir/matrixreader.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/matrixreader.cpp.o: ../matrixreader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/matrix_numerov.dir/matrixreader.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/matrixreader.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/matrixreader.cpp

CMakeFiles/matrix_numerov.dir/matrixreader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/matrixreader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/matrixreader.cpp > CMakeFiles/matrix_numerov.dir/matrixreader.cpp.i

CMakeFiles/matrix_numerov.dir/matrixreader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/matrixreader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/matrixreader.cpp -o CMakeFiles/matrix_numerov.dir/matrixreader.cpp.s

CMakeFiles/matrix_numerov.dir/parameters.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/parameters.cpp.o: ../parameters.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/matrix_numerov.dir/parameters.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/parameters.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/parameters.cpp

CMakeFiles/matrix_numerov.dir/parameters.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/parameters.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/parameters.cpp > CMakeFiles/matrix_numerov.dir/parameters.cpp.i

CMakeFiles/matrix_numerov.dir/parameters.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/parameters.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/parameters.cpp -o CMakeFiles/matrix_numerov.dir/parameters.cpp.s

CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.o: CMakeFiles/matrix_numerov.dir/flags.make
CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.o: ../eigenvalue.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.o -c /home/artfin/Desktop/repos/NumerovProject/eigenvalue.cpp

CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/artfin/Desktop/repos/NumerovProject/eigenvalue.cpp > CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.i

CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/artfin/Desktop/repos/NumerovProject/eigenvalue.cpp -o CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.s

# Object files for target matrix_numerov
matrix_numerov_OBJECTS = \
"CMakeFiles/matrix_numerov.dir/main.cpp.o" \
"CMakeFiles/matrix_numerov.dir/filereader.cpp.o" \
"CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.o" \
"CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.o" \
"CMakeFiles/matrix_numerov.dir/matrixreader.cpp.o" \
"CMakeFiles/matrix_numerov.dir/parameters.cpp.o" \
"CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.o"

# External object files for target matrix_numerov
matrix_numerov_EXTERNAL_OBJECTS =

matrix_numerov: CMakeFiles/matrix_numerov.dir/main.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/filereader.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/generalizedmatrixnumerov.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/matrixnumerov.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/matrixreader.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/parameters.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/eigenvalue.cpp.o
matrix_numerov: CMakeFiles/matrix_numerov.dir/build.make
matrix_numerov: /usr/lib/x86_64-linux-gnu/libgsl.so
matrix_numerov: /usr/lib/x86_64-linux-gnu/libgslcblas.so
matrix_numerov: CMakeFiles/matrix_numerov.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable matrix_numerov"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix_numerov.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/matrix_numerov.dir/build: matrix_numerov

.PHONY : CMakeFiles/matrix_numerov.dir/build

CMakeFiles/matrix_numerov.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/matrix_numerov.dir/cmake_clean.cmake
.PHONY : CMakeFiles/matrix_numerov.dir/clean

CMakeFiles/matrix_numerov.dir/depend:
	cd /home/artfin/Desktop/repos/NumerovProject/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/artfin/Desktop/repos/NumerovProject /home/artfin/Desktop/repos/NumerovProject /home/artfin/Desktop/repos/NumerovProject/cmake-build-debug /home/artfin/Desktop/repos/NumerovProject/cmake-build-debug /home/artfin/Desktop/repos/NumerovProject/cmake-build-debug/CMakeFiles/matrix_numerov.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/matrix_numerov.dir/depend
