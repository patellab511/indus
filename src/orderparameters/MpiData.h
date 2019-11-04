#ifndef MPI_DATA_H
#define MPI_DATA_H

#include <array>
#include <complex>
#include <vector>

// Library headers
#ifdef MPI_ENABLED
#include <mpi.h>
#endif /* MPI_ENABLED */

#include "MpiDatatype.h"
#include "utils.h"

// Forward declaration to prevent circular dependency
class MpiCommunicator;

// Wrapper around several variables that define an MPI buffer
// - A "data element" is a single piece of data of type "data_type"
// - The difference between "size" and "max_size" has the same semantics
//   as the difference between std::vector's size and capacity
// - Note: need to pass an MpiCommunicator to constructors because it contains
//   the mapping from C++ types to MPI_Datatypes
class MpiData {
 public:
	// Pointer to contiguous block of memory where data may be stored,
	// and the maximum number of data elements (of type "data_type") 
	// that can be stored at that location
	void* data_ptr;
	int   max_size;

	// Number of elements of MPI_Datatype "data_type" actually stored 
	// in the buffer
	int size;

	// Type of data; used to determine the size of the buffer in bytes
	MPI_Datatype data_type;
	//MpiDatatype data_type; 
	// FIXME: need to figure out how to handle comparing types which are equivalent

	//----- Constructors -----//

	// Default: empty
	MpiData(); 

	// From a reference to a single element
	template<typename T> explicit MpiData(T& data, MpiCommunicator& comm);

	// From a pointer to a single data element or array of primitive type
	// - Default: pointer to single element
	template<typename T> explicit MpiData(T* data_ptr_in, const int num_elements_in, MpiCommunicator& comm);

	// From a std::vector of a registered MPI type
	template<typename T> MpiData(std::vector<T>& vec, MpiCommunicator& comm);
};


//----- Templates -----//

#include "MpiCommunicator.h"

// From a reference to a supported contiguous type
template<typename T>
MpiData::MpiData(T& data, MpiCommunicator& comm)
 : data_ptr( static_cast<void*>(&data) ), 
   max_size(1),
   size(1), 
   data_type( comm.map_mpi_datatype<T>().get_MPI_Datatype() )
{}


// From a raw pointer to a supported contiguous type
template<typename T>
MpiData::MpiData(T* data_ptr_in, const int num_elements_in, MpiCommunicator& comm)
 : data_ptr( static_cast<void*>(data_ptr_in) ), 
   max_size(num_elements_in),
   size(num_elements_in), 
   data_type( comm.map_mpi_datatype<T>().get_MPI_Datatype() )
{
	if ( data_ptr == nullptr ) {
		size = 0;
	}
}

// From a std::vector of a primitive type
template<typename T> 
MpiData::MpiData(std::vector<T>& vec, MpiCommunicator& comm)
 : data_ptr( static_cast<void*>(vec.data()) ),
   max_size(vec.capacity()), 
   size(vec.size()), 
   data_type( comm.map_mpi_datatype<T>().get_MPI_Datatype() )
{}

#endif /* MPI_DATA_H */
