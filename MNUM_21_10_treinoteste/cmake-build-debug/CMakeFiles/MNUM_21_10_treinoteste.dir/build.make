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
CMAKE_SOURCE_DIR = "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/MNUM_21_10_treinoteste.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MNUM_21_10_treinoteste.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MNUM_21_10_treinoteste.dir/flags.make

CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.obj: CMakeFiles/MNUM_21_10_treinoteste.dir/flags.make
CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\MNUM_21_10_treinoteste.dir\main.cpp.obj -c "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\main.cpp"

CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\main.cpp" > CMakeFiles\MNUM_21_10_treinoteste.dir\main.cpp.i

CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\main.cpp" -o CMakeFiles\MNUM_21_10_treinoteste.dir\main.cpp.s

# Object files for target MNUM_21_10_treinoteste
MNUM_21_10_treinoteste_OBJECTS = \
"CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.obj"

# External object files for target MNUM_21_10_treinoteste
MNUM_21_10_treinoteste_EXTERNAL_OBJECTS =

MNUM_21_10_treinoteste.exe: CMakeFiles/MNUM_21_10_treinoteste.dir/main.cpp.obj
MNUM_21_10_treinoteste.exe: CMakeFiles/MNUM_21_10_treinoteste.dir/build.make
MNUM_21_10_treinoteste.exe: CMakeFiles/MNUM_21_10_treinoteste.dir/linklibs.rsp
MNUM_21_10_treinoteste.exe: CMakeFiles/MNUM_21_10_treinoteste.dir/objects1.rsp
MNUM_21_10_treinoteste.exe: CMakeFiles/MNUM_21_10_treinoteste.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MNUM_21_10_treinoteste.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\MNUM_21_10_treinoteste.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MNUM_21_10_treinoteste.dir/build: MNUM_21_10_treinoteste.exe

.PHONY : CMakeFiles/MNUM_21_10_treinoteste.dir/build

CMakeFiles/MNUM_21_10_treinoteste.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\MNUM_21_10_treinoteste.dir\cmake_clean.cmake
.PHONY : CMakeFiles/MNUM_21_10_treinoteste.dir/clean

CMakeFiles/MNUM_21_10_treinoteste.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\cmake-build-debug" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\cmake-build-debug" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\MNUM_21_10_treinoteste\cmake-build-debug\CMakeFiles\MNUM_21_10_treinoteste.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/MNUM_21_10_treinoteste.dir/depend
