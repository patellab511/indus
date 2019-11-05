# Testing
# - All tests run the program for a set of sample input, and compare the resulting output files to 
#   "gold standard" output files stored in the repo
#   - Simplistic, but easy to implement and modify
# - TODO: Convert from bash-driven testing to cmake-driven testing
#   - Would make a number of things much cleaner
#   - Could possibly leverage more cmake/ctest features

#############
### Setup ###
#############

# Find bash
find_program (BASH_PROGRAM bash)
if(NOT BASH_PROGRAM)
	message(FATAL_ERROR "Unable to find bash for testing")
endif()

# Testing driver script (absolute path)
set(TEST_DRIVER_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/run_test.sh)

set(ECHO_FAILED_DIFFS "-diff" CACHE STRING "tests: whether to print diffs for tests that fail")

# Path to program (absolute)
set(PROGRAM_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${INDUS_DRIVER}")

# Directory with testing files
# - This is inside <repo>/test, *not* the build directory
set(TEST_DIR_ROOT ${CMAKE_CURRENT_LIST_DIR})


#############
### Tests ###
#############

### Builtin Tests ###

# Root directory where testing binaries will be built
set(TEST_BIN_ROOT ${CMAKE_BINARY_DIR}/test/bin)

# Name of the test driver program
# - Tests are differentiated by their directories
set(TEST_PROGRAM test_driver)

# CellGrid
# - TODO: move to function: build test driver with custom 'main' (link to library)
set(TEST_NAME     cell_grid/basic)                 # directory under <repo>/test for this test
set(TEST_DIR      ${TEST_DIR_ROOT}/${TEST_NAME})   # directory with source files for this test
set(TEST_BIN      ${TEST_BIN_ROOT}/${TEST_NAME})   # where to build this testing binary
add_executable(${TEST_PROGRAM} ${TEST_DIR}/main.cpp)
set_target_properties(${TEST_PROGRAM}
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${TEST_BIN}
                      )
target_link_libraries(${TEST_PROGRAM} PUBLIC ${LIB_INDUS})
target_include_directories(${TEST_PROGRAM} PUBLIC ${CMAKE_SOURCE_DIR}/src)
# - Force target 'check' (i.e. 'make check') to depend on building the test executable;
#   otherwise, it won't get built
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TEST_PROGRAM})
add_test(NAME ${TEST_NAME}
         WORKING_DIRECTORY ${TEST_BIN}
         COMMAND ${TEST_PROGRAM}
         )


### Serial Tests ###

# TODO List from searching file names?
set(SERIAL_TESTS "")
list(APPEND SERIAL_TESTS "indus/box/bulk_water")
list(APPEND SERIAL_TESTS "indus/cylinder/bulk_water")
list(APPEND SERIAL_TESTS "indus/sphere_near_box_edge/bulk_water")

foreach(TEST_NAME IN LISTS SERIAL_TESTS)
	set(TEST_DIR ${TEST_DIR_ROOT}/${TEST_NAME})
	add_test(NAME ${TEST_NAME}
					 COMMAND ${BASH_PROGRAM} ${TEST_DRIVER_SCRIPT} -e ${PROGRAM_PATH} -d ${TEST_DIR} ${ECHO_FAILED_DIFFS}
					 )
endforeach()


### Parallel Tests ###

# TODO Testings for MPI *or* OpenMP, not just both at the same time?
if(MPI_ENABLED AND OPENMP_ENABLED)
	# TODO allow this to be specified on the command line
	set(MPI_NUM_RANKS   4)
	set(OMP_NUM_THREADS 2)

	# TODO List from file/searching directories?
	set(PARALLEL_TESTS "")
	list(APPEND PARALLEL_TESTS "indus_mpi/sphere_near_box_edge/bulk_water")

	# TODO set WORKING_DIRECTORY?
	foreach(TEST_NAME IN LISTS PARALLEL_TESTS)
		set(TEST_DIR ${TEST_DIR_ROOT}/${TEST_NAME})
		add_test(NAME ${TEST_NAME}
						 COMMAND ${BASH_PROGRAM} ${TEST_DRIVER_SCRIPT} -e ${PROGRAM_PATH} -d ${TEST_DIR} ${ECHO_FAILED_DIFFS}
										 -np ${MPI_NUM_RANKS} -nt ${OMP_NUM_THREADS}
						 )
		         #WORKING_DIRECTORY ${TEST_ROOT_DIR}
	endforeach()
endif()
