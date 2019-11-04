#include "main.h"

int main(int argc, char* argv[]) 
{
	try {
		// Input checking
		if ( argc < 2 ) {
			std::cerr << "indus: ERROR: missing command-line input file" << std::endl;
			return 1;
		}
#ifdef MPI_ENABLED
		// Start MPI
		MpiEnvironment::init();
#endif /* MPI_ENABLED */

		// Place the following in a different scope so that the MPI communicator
		// in OrderParametersDriver in cleaned up before MPI_Finalize() is called
		{
			// GPTL: start
			GPTL::GlobalOptions options;  // sets default options
			GPTL::initialize();           // Initialize GPTL
			GPTL::Timer timer("main");
			timer.start();

			// Input file path
			std::string op_input_file(argv[1]);

			// Initialize communicator
			MpiCommunicator mpi_communicator( MpiCommunicator::get_mpi_comm_world() );

			// Run calculation
			OrderParametersDriver driver(op_input_file, mpi_communicator);
			driver.run();

			// GPTL: done
			timer.stop();
			driver.finalize();  // prints output files
		}

#ifdef MPI_ENABLED
		// Finish MPI
		MpiEnvironment::finalize();
#endif /* MPI_ENABLED */
	} // end try
	catch (const std::exception& ex) {
		std::cout << "main caught an exception\n"
		          << "  what() = " << ex.what() << "\n"
		          << "\n"
		          << "The program will now exit.\n";
		exit(1);
	}

	return 0;
} // end main
