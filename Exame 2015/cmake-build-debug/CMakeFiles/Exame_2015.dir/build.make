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
CMAKE_SOURCE_DIR = "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/Exame_2015.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Exame_2015.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Exame_2015.dir/flags.make

CMakeFiles/Exame_2015.dir/main.cpp.obj: CMakeFiles/Exame_2015.dir/flags.make
CMakeFiles/Exame_2015.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Exame_2015.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Exame_2015.dir\main.cpp.obj -c "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\main.cpp"

CMakeFiles/Exame_2015.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Exame_2015.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\main.cpp" > CMakeFiles\Exame_2015.dir\main.cpp.i

CMakeFiles/Exame_2015.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Exame_2015.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\main.cpp" -o CMakeFiles\Exame_2015.dir\main.cpp.s

# Object files for target Exame_2015
Exame_2015_OBJECTS = \
"CMakeFiles/Exame_2015.dir/main.cpp.obj"

# External object files for target Exame_2015
Exame_2015_EXTERNAL_OBJECTS =

Exame_2015.exe: CMakeFiles/Exame_2015.dir/main.cpp.obj
Exame_2015.exe: CMakeFiles/Exame_2015.dir/build.make
Exame_2015.exe: CMakeFiles/Exame_2015.dir/linklibs.rsp
Exame_2015.exe: CMakeFiles/Exame_2015.dir/objects1.rsp
Exame_2015.exe: CMakeFiles/Exame_2015.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Exame_2015.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Exame_2015.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Exame_2015.dir/build: Exame_2015.exe

.PHONY : CMakeFiles/Exame_2015.dir/build

CMakeFiles/Exame_2015.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Exame_2015.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Exame_2015.dir/clean

CMakeFiles/Exame_2015.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\cmake-build-debug" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\cmake-build-debug" "C:\Data\Andre\Work\MIEIC - ANO 2 - SEMESTRE 1\MNUM\Exame 2015\cmake-build-debug\CMakeFiles\Exame_2015.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/Exame_2015.dir/depend
