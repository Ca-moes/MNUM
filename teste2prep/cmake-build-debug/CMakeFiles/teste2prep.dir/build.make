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
CMAKE_SOURCE_DIR = "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/teste2prep.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/teste2prep.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/teste2prep.dir/flags.make

CMakeFiles/teste2prep.dir/main.cpp.obj: CMakeFiles/teste2prep.dir/flags.make
CMakeFiles/teste2prep.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/teste2prep.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\teste2prep.dir\main.cpp.obj -c "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\main.cpp"

CMakeFiles/teste2prep.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/teste2prep.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\main.cpp" > CMakeFiles\teste2prep.dir\main.cpp.i

CMakeFiles/teste2prep.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/teste2prep.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\main.cpp" -o CMakeFiles\teste2prep.dir\main.cpp.s

# Object files for target teste2prep
teste2prep_OBJECTS = \
"CMakeFiles/teste2prep.dir/main.cpp.obj"

# External object files for target teste2prep
teste2prep_EXTERNAL_OBJECTS =

teste2prep.exe: CMakeFiles/teste2prep.dir/main.cpp.obj
teste2prep.exe: CMakeFiles/teste2prep.dir/build.make
teste2prep.exe: CMakeFiles/teste2prep.dir/linklibs.rsp
teste2prep.exe: CMakeFiles/teste2prep.dir/objects1.rsp
teste2prep.exe: CMakeFiles/teste2prep.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable teste2prep.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\teste2prep.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/teste2prep.dir/build: teste2prep.exe

.PHONY : CMakeFiles/teste2prep.dir/build

CMakeFiles/teste2prep.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\teste2prep.dir\cmake_clean.cmake
.PHONY : CMakeFiles/teste2prep.dir/clean

CMakeFiles/teste2prep.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\cmake-build-debug" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\cmake-build-debug" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\teste2prep\cmake-build-debug\CMakeFiles\teste2prep.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/teste2prep.dir/depend

