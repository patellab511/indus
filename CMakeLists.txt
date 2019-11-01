cmake_minimum_required(VERSION 3.15)

# Include directory with custom CMake modules
# - These contain useful helper functions
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(setIfUnset)

project(indus)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Default to 'Release' build type
setIfUnset(CMAKE_BUILD_TYPE Release)

# Verbose output
# - TODO: restrict to 'debug'?
set(CMAKE_VERBOSE_MAKEFILE ON)

# Put all binaries and libraries into subdirectories, instead of the main 'build' directory
# - Adapted from CMake tutorial (v3.15)
setIfUnset(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
setIfUnset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
setIfUnset(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")


#######################################
### Important variables and options ###
#######################################

set(XDRFILE_DIR "${HOME}/programs/xdrfile" CACHE STRING "Location of xdrfile installation")

# GPTL (for benchmarking)
option(GPTL_ENABLED "Enable GPTL" ON)
set(GPTL_DIR "${HOME}/programs/gptl" CACHE STRING "Location of GPTL installation")

# Parallelization
option(MPI_ENABLED "Enable MPI" ON)
option(OPENMP_ENABLED "Enable OpenMP" ON)

# Output names
set(LIB_INDUS    "libindus"        CACHE STRING "Name of the library to be created")
set(INDUS_DRIVER "${PROJECT_NAME}" CACHE STRING "Name of the executable")


#################
### Libraries ###
#################

# TODO: check here that libraries will work?

### xdrfile ###

message("-- xdrfile location: ${XDRFILE_DIR}")

# Look for static library
set(LIB_XDRFILE "xdrfile")
find_library(XDRFILE_LIB_STATIC
             NAMES lib${LIB_XDRFILE}.a ${LIB_XDRFILE}.a
             HINTS "${XDRFILE_DIR}/lib"
             )
if(XDRFILE_LIB_STATIC)
	message("-- Found static library for ${LIB_XDRFILE}: ${XDRFILE_LIB_STATIC}")
endif()


### GPTL ###

if(GPTL_ENABLED)
	# TODO: make target-wise? put in config file, that also gets included in the plumed patch?
	add_definitions(-DHAVE_GPTL)

	message("-- gptl enabled")
	message("-- gptl location: ${GPTL_DIR}")

	# Look for static library
	set(LIB_GPTL "gptl")
	find_library(LIB_GPTL_STATIC
							 NAMES lib${LIB_GPTL}.a ${LIB_GPTL}.a
							 HINTS "${GPTL_DIR}/lib"
							 )
	if(LIB_GPTL_STATIC)
		message("-- Found static library for ${LIB_GPTL}: ${LIB_GPTL_STATIC}")
	endif()
endif()


### MPI ###

if(MPI_ENABLED)
	message("-- MPI enabled")
	# Add MPI
	find_package(MPI)
	if(MPI_FOUND)
		# TODO Move to config.h?
		add_definitions(-DMPI_ENABLED)
		add_definitions(-DHAVE_LIBMPI)
	else()
    message(FATAL_ERROR "MPI could not be located")
	endif()
endif()


### OpenMP ###

if(OPENMP_ENABLED)
	message("-- OpenMP enabled")
	# Add OpenMP
	find_package(OpenMP)
	if(NOT OpenMP_FOUND)
		if(NOT APPLE)
			message(FATAL_ERROR "OpenMP could not be located")
		else()
			# Mac is notoriously difficult: try to find Homebrew OpenMP semi-manually
			message("-- Will try to find Homebrew OpenMP")
			include(findHomebrewOpenMP)
			findHomebrewOpenMP()
		endif()
	endif()

	# Only set a linker flag if the user didn't set them manually
	if(NOT DEFINED OpenMP_LINKER_FLAGS)
			set(OpenMP_LINKER_FLAGS "${OpenMP_CXX_FLAGS}")
	endif()
	if(NOT DEFINED OpenMP_SHARED_LINKER_FLAGS)
			set(OpenMP_SHARED_LINKER_FLAGS "${OpenMP_CXX_FLAGS}")
	endif()
else()
	# Unset all OpenMP flags in case OpenMP was disabled either by the user
	# or because it was only partially detected (e.g. only for C but not C++ compiler)
	unset(OpenMP_C_FLAGS CACHE)
	unset(OpenMP_CXX_FLAGS CACHE)
endif()


######################
### Compiler Flags ###
######################

# Basic flags
add_compile_options(-g -Wall)
#target_compile_options(${PROJECT_NAME} PUBLIC -g -Wall)

# Optimizations
set(INDUS_OPT_DEBUG_FLAGS   -O0)
set(INDUS_OPT_RELEASE_FLAGS -O3 -ffast-math -march=native)
if(GNU)
	# FIXME this does not seem to be working
	string(APPEND INDUS_OPT_RELEASE_FLAGS -ftree-vectorize)
endif()
add_compile_options("$<$<CONFIG:DEBUG>:${INDUS_OPT_DEBUG_FLAGS}>")
add_compile_options("$<$<CONFIG:RELEASE>:${INDUS_OPT_RELEASE_FLAGS}>")


#############
### Build ###
#############

add_definitions(-DINDUS_STANDALONE_MODE)
add_subdirectory(src)


############
### Test ###
############

enable_testing()
add_subdirectory(test)
