#include "GptlWrappers.h"

// Parallelization
#include "MpiCommunicator.h"
#include "OpenMP.h"

namespace GPTL {

//----- Loose Functions -----//

// This flag helps prevent multiple calls to GPTLinitialize
// - Please don't mess with it. Please.
// - TODO Move to a class for encapsulation?
static bool is_gptl_initialized_ = false;

void initialize()
{
	#ifdef HAVE_GPTL
	// If OpenMP is enabled, only the master thread should initialize GPTL
	// - Multiple initializations are an error
	// - It seems to be okay if MPI_Init has already been called, and multiple
	//   independent ranks call GPTLinitialize()
	//   - TODO: Verify with devs
	// - If the program is not in a parallel region when initialize() is called,
	//   the pragmas are ignored
	// - Note that 'master' does not impose an implicit barrier, so an explicit one is needed
	// - TODO Is it possible to make this threadsafe for use with other kinds of threads?
	//   - pthreads? C++ threads?
	#pragma omp master
	{
		if ( not is_gptl_initialized_ ) {
			GPTLinitialize();
			is_gptl_initialized_ = true;
		}
	}
	#pragma omp barrier
	#endif // ifdef HAVE_GPTL
}


void print(const std::string& file) 
{
	#ifdef HAVE_GPTL
	GPTLpr_file( file.c_str() );
	#endif // ifdef HAVE_GPTL
}


void print(const std::string& file, const MpiCommunicator& mpi_communicator)
{
	#ifdef HAVE_GPTL
	int my_rank = mpi_communicator.get_rank();
	std::string rank_file = file + "." + std::to_string(my_rank);
	GPTLpr_file( rank_file.c_str() );
	#endif // ifdef HAVE_GPTL
}


void printSummary(const std::string& file, const MpiCommunicator& mpi_communicator)
{
	#ifdef HAVE_GPTL
	GPTLpr_summary_file( mpi_communicator.get_communicator(), file.c_str() );
	#endif // ifdef HAVE_GPTL
}


//----- Timer -----//

Timer::Timer(const std::string& name):
	name_(name)
{
	#ifdef HAVE_GPTL
	GPTLinit_handle( name_.c_str(), &hash_index_ );
	#endif
}


void Timer::start()
{
	#ifdef HAVE_GPTL
	GPTLstart_handle( name_.c_str(), &hash_index_ );
	//GPTLstart( name_.c_str() );
	#endif
}


void Timer::stop()
{
	#ifdef HAVE_GPTL
	GPTLstop_handle( name_.c_str(), &hash_index_ );
	//GPTLstop( name_.c_str() );
	#endif // ifdef HAVE_GPTL
}


//----- GlobalOptions -----//

GlobalOptions::GlobalOptions()
{	
	setDefaults();
}


void GlobalOptions::setDefaults()
{
	#ifdef HAVE_GPTL
	GPTLsetoption(GPTLpercent, 1);  // Print a column with % of first timer
	//GPTLsetoption(GPTLoverhead, 0);  // Don't print overhead estimate
	#endif // ifdef HAVE_GPTL
}

} // end namespace GPTL
