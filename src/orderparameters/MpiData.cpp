#include "MpiData.h"
//#include "MpiCommunicator.h" // FIXME redundant

// Default constructor
MpiData::MpiData()
 : data_ptr(nullptr), max_size(0), size(0), 
#ifdef MPI_ENABLED
	// FIXME
	//data_type(MPI_DATATYPE_NULL) 
	data_type() 
#else
	data_type()
#endif /* MPI_ENABLED */
{}
