#include "MpiEnvironment.h"


namespace MpiEnvironment {


bool is_enabled()
{
#ifdef MPI_ENABLED
	return true;
#else
	return false;
#endif /* MPI_ENABLED */
}


bool is_initialized()
{
	int is_init = false;
#ifdef MPI_ENABLED
	MPI_Initialized(&is_init);
#endif /* MPI_ENABLED */
	if ( is_init ) { return true;  }
	else           { return false; }
}


void initialize()
{
#ifdef MPI_ENABLED
	// Start MPI
	int    dummy_argc = 0;
	char** dummy_argv = nullptr;
	MPI_Init(&dummy_argc, &dummy_argv);
#endif // ifdef MPI_ENABLED
}


void finalize()
{
#ifdef MPI_ENABLED
	MPI_Finalize();
#endif // ifdef MPI_ENABLED
}


void init()
{
	initialize();
}


} // end namespace MpiEnvironment
