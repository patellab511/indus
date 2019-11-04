// OpenMP - static wrappers for OpenMP functions

#pragma once
#ifndef OPEN_MP_H
#define OPEN_MP_H

#include <cstddef>  // defines std::size_t

#ifdef _OPENMP
#include <omp.h>
#endif // ifdef _OPENMP

#include "System.h"

class OpenMP {
 public:
	// Returns whether the program has been compiled with OpenMP
	static bool is_enabled();

	static int get_num_threads();

	static int get_max_threads();

	// Returns the 'number' (ID) of the thread, indexed from 0
	static int get_thread_num();

	static int get_num_procs();

	// Returns 'true' is the program is in an OpenMP 'parallel' region
	static bool in_parallel();

	// TODO Remaining functions

	// TODO omp_get_wtime / omp_get_wtick

	// Returns the cache line size (in bytes)
	static std::size_t get_cache_line_size();
};

// TODO locks
// class OpenMP_Lock

#endif // ifndef OPEN_MP_H
