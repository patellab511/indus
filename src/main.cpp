/*	main.cpp
 *
 *	ABOUT:
 *	AUTHOR: Sean M. Marks
 *	NOTES:
 *	DEVELOPMENT:
 */

#include "main.h"


/*
	Input:
		argv[1] = indus_cpp.input options file
		argv[2] = gro file (for single-frame calculations/debugging)
*/
int main(int argc, char* argv[]) 
{
	try {
		// Input
		if ( argc < 2 ) {
			std::cerr << "main: ERROR: Must submit an input file!" << std::endl;
			return 1;
		}
		std::string opInputFile(argv[1]);

	#ifdef MPI_ENABLED
		// Start MPI
		int    dummy_argc = 0;
		char** dummy_argv = nullptr;
		MPI_Init(&dummy_argc, &dummy_argv);
	#endif /* MPI_ENABLED */

		// Place the following in a different scope so that the MPI communicator
		// in Indus in cleaned up before MPI_Finalize() is called
		{
			Indus indus(opInputFile);	
			indus.do_indus_standalone();
		}

	#ifdef MPI_ENABLED
		// Finish MPI
		MPI_Finalize();
	#endif /* MPI_ENABLED */
	}
	catch (const std::exception& ex) {
		std::cout << "main caught an exception\n"
		          << "  what() = " << ex.what() << "\n"
		          << "\n"
		          << "The program will now exit.\n";
		exit(1);
	}

	return 0;
} // main
