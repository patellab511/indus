function(findHomebrewOpenMP)
	# Look for Homebrew
	find_program(HOMEBREW_PROGRAM brew)
	if(NOT HOMEBREW_PROGRAM)
		message(FATAL_ERROR "Unable to find Homebrew")
	endif()

	# Get the location of Homebrew OpenMP
	execute_process(COMMAND ${HOMEBREW_PROGRAM} --prefix libomp
									OUTPUT_VARIABLE HOMEBREW_OPENMP_PREFIX
									OUTPUT_STRIP_TRAILING_WHITESPACE
									)

	# Try again to find OpenMP
	# - Compiler flags
	set(OpenMP_C_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_OPENMP_PREFIX}/include")
	set(OpenMP_C_FLAGS     ${OpenMP_C_CXX_FLAGS})
	set(OpenMP_CXX_FLAGS   ${OpenMP_C_CXX_FLAGS})
	# - Library
	set(OpenMP_C_LIB_NAMES   "omp")
	set(OpenMP_CXX_LIB_NAMES "omp")
	set(OpenMP_omp_LIBRARY   "${HOMEBREW_OPENMP_PREFIX}/lib/libomp.dylib")
	# - Try again
	find_package(OpenMP)
	if(OpenMP_FOUND)
		message("-- Found Homebrew OpenMP at: ${HOMEBREW_OPENMP_PREFIX}")
	else()
		message(FATAL_ERROR "Homebrew OpenMP not found")
	endif()
endfunction()
