# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.15

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2019.2.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2019.2.2\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MNUM_23_09.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MNUM_23_09.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MNUM_23_09.dir/flags.make

CMakeFiles/MNUM_23_09.dir/main.cpp.obj: CMakeFiles/MNUM_23_09.dir/flags.make
CMakeFiles/MNUM_23_09.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MNUM_23_09.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\MNUM_23_09.dir\main.cpp.obj -c C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\main.cpp

CMakeFiles/MNUM_23_09.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MNUM_23_09.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\main.cpp > CMakeFiles\MNUM_23_09.dir\main.cpp.i

CMakeFiles/MNUM_23_09.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MNUM_23_09.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\main.cpp -o CMakeFiles\MNUM_23_09.dir\main.cpp.s

# Object files for target MNUM_23_09
MNUM_23_09_OBJECTS = \
"CMakeFiles/MNUM_23_09.dir/main.cpp.obj"

# External object files for target MNUM_23_09
MNUM_23_09_EXTERNAL_OBJECTS =

MNUM_23_09.exe: CMakeFiles/MNUM_23_09.dir/main.cpp.obj
MNUM_23_09.exe: CMakeFiles/MNUM_23_09.dir/build.make
MNUM_23_09.exe: CMakeFiles/MNUM_23_09.dir/linklibs.rsp
MNUM_23_09.exe: CMakeFiles/MNUM_23_09.dir/objects1.rsp
MNUM_23_09.exe: CMakeFiles/MNUM_23_09.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MNUM_23_09.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\MNUM_23_09.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MNUM_23_09.dir/build: MNUM_23_09.exe

.PHONY : CMakeFiles/MNUM_23_09.dir/build

CMakeFiles/MNUM_23_09.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\MNUM_23_09.dir\cmake_clean.cmake
.PHONY : CMakeFiles/MNUM_23_09.dir/clean

CMakeFiles/MNUM_23_09.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09 C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09 C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\cmake-build-debug C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\cmake-build-debug C:\Users\andre_4mxog39\CLionProjects\MNUM_23_09\cmake-build-debug\CMakeFiles\MNUM_23_09.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MNUM_23_09.dir/depend

