#include "OpenMP.h"

bool OpenMP::is_enabled()
{
#ifdef _OPENMP
	return true;
#else
	return false;
#endif // ifdef _OPENMP
}

int OpenMP::get_num_threads()
{
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif // ifdef _OPENMP
}


int OpenMP::get_max_threads()
{
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif // ifdef _OPENMP
}


int OpenMP::get_thread_num()
{
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif // ifdef _OPENMP
}


int OpenMP::get_num_procs()
{
#ifdef _OPENMP
	return omp_get_num_procs();
#else
	return 1;
#endif // ifdef _OPENMP
}


bool OpenMP::in_parallel()
{
#ifdef _OPENMP
	return static_cast<bool>( omp_in_parallel() );
#else
	return false;
#endif // ifdef _OPENMP
}


std::size_t OpenMP::get_cache_line_size()
{
	return System::get_cache_line_size();
}

/*
std::size_t OpenMP::get_cache_line_size()
{
	// Note: PLUMED uses 512?!
	// - Cache lines are typically 64 bytes --> 64*8 = 512 bits (could be the explanation)
	// - For a page size of 4096 bytes, note: 4096/512 = 8
	unsigned num_bytes = 64;
	return num_bytes;
}
*/

