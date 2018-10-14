// MpiCommunicator
//
// ABOUT: Wrapper around MPI communicator and other library functions
//	- Based heavily on the Communicator implemented in PLUMED 2.4.0 (see http://plumed.org)
//  - Contains lots of template code for commonly-used STL containers that store
//    primitive types (e.g. int, double) in contiguous memory
//  - Shields the user from many of the tedious aspects of MPI library calls
// 
// FOR THE USER:
//  - Make sure std::vector buffers are prepared using resize(), not just reserve()
//
// NOTES:
//	- Most functions (e.g. get/set, send/recv) throw if MPI is not initialized
//  - MPI routines here almost always *assume* that the receive buffers have enough
//		memory reserved to handle the maximum amount of data that may be sent at a time.
//    They do not pre-reserve space for you, but will auto-resize based on the 
//    amount of data received.
//    - They are not designed to receive data over multiple calls
//    - Exceptions:
//       - allgatherv(): calls gather() to determine how much data each rank will send,
//                       and resizes all output arrays accordingly
//
//  - send functions:
//    - Infer the number of data elements to send based on MpiData::size
//
//  - recv functions:
//    - Set MpiData::size based on how much data was received 
//      (except if MPI_STATUS_IGNORE is requested)
//       - If MPI_STATUS_IGNORE is requested, the functions assume that
//         all available space was filled
//          - e.g. vectors are filled to capacity, irrespective of starting size
//       - If you want a vector to receive an amount of data other
//         than its current size, you *MUST* pass it an MpiStatus object, even if
//         you do nothing with it after the call
//
// DEVELOPMENT: TODO
//  - Isend routines
//  - Include PLMD::Vector overloads?
//  - Pointer overloads
//  - Generalize template code for std::vectors of contiguous data arrays?
//    - Such a setup would cover std::vector<T>, PLMD::Vector, std::array<T,dim>, ...
//    - Could make use of MPI's own custom type system with MPI_Type_contiguous
//      - Reference:  https://tech.io/playgrounds/349/introduction-to-mpi/custom-types

#ifndef MPI_COMMUNICATOR_H
#define MPI_COMMUNICATOR_H

// Turn on MPI_ENABLED when compiling for PLUMED
#if ( defined(__PLUMED_HAS_MPI) && ! defined(MPI_ENABLED) )
#define MPI_ENABLED
#endif


// Standard headers
#include <array>
#include <complex>
#include <cstdlib>
#include <exception>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Library headers
#ifdef MPI_ENABLED
#include <mpi.h>
#endif /* MPI_ENABLED */

// Project headers TODO


#ifndef MPI_ENABLED
// Dummy classes for when MPI is not available
class MPI_Comm {};
class MPI_Datatype {};
class MPI_Op {};
class MPI_Status {};
class MPI_Request {};

/*
// Dummy MPI_Init
static void MPI_Init(int* argc, char** argv[]);

// Dummy MPI_Finalize
static void MPI_Finalize();
*/
#endif /* MPI_ENABLED */


class MpiCommunicator
{
 public:
	//-------------------------------------------//
	//----- Constructors and Initialization -----//
	//-------------------------------------------//

	// Default constructor
	// - Sets the communicator to MPI_COMM_SELF
	// - Will do MPI_Init if necessary
	MpiCommunicator();

	// Initialize from raw MPI communicator
	// - This is the "primary" constructor
	// - Duplicates the provided communicator, creating a new context
	// - Will do MPI_Init if necessary
	MpiCommunicator(const MPI_Comm& communicator);

	// Copy constructor 
	MpiCommunicator(const MpiCommunicator& communicator);

	// Destructor
	~MpiCommunicator();

	// Checks whether MPI is available (whether the code was compiled with MPI)
	static bool is_mpi_enabled();

	// Checks whether the MPI library has been initialized
	// - Returns "false" if MPI is not enabled
	// - Note: MPI mandates that MPI_Init may only be called *once* by each process
	static bool is_mpi_initialized();

	// Returns MPI_COMM_WORLD if MPI is enabled, else returns a dummy MPI_Comm
	// - Use this in order to avoid any direct references to MPI_COMM_WORLD in your code,
	//   so that in can compile if MPI is not enabled
	static MPI_Comm get_mpi_comm_world();

	// Returns MPI_COMM_SELF if MPI is enabled, else returns a dummy MPI_Comm
	// - Use this in order to avoid any direct references to MPI_COMM_WORLD in your code,
	//   so that in can compile if MPI is not enabled
	static MPI_Comm get_mpi_comm_self();

	// Assuming a dim-dimensional system in Cartesian space, this function
	// uses MPI_Dims_create() to determine the number of domain decomposition
	// cells that should be placed along each axis
	template<std::size_t dim>
	std::array<int,dim> calculateGridDimensions() const;


	//-----------------------------//
	//----- Get/Set Functions -----//
	//-----------------------------//

	// Duplicates the communicator provided
	void set_communicator(const MPI_Comm& communicator);
	void set_communicator(void* communicator_ptr);

	// Returns the rank of this process in the commmunicator
	// - Returns 0 if MPI is not enabled
	int get_rank() const;

	int get_master_rank() const { return master_rank_; }

	// Returns the number of ranks in this communicator
	// - Returns 1 if MPI is not enabled
	int get_size() const;

	// Static mapping from primitive C++ types to MPI types
	template<typename T>
	static MPI_Datatype get_MPI_Datatype();

	bool is_serial() const {
		return ( this->get_size() == 1 );
	}

	//-------------------------------//
	//----- MPI Struct Wrappers -----//
	//-------------------------------//

	// MPI_Status wrapper
	class MpiStatus {
	 public:
		MpiStatus() {}
		MpiStatus(const MPI_Status& status) 
		 : status_(status) {}

		// Get the number of elements transferred
		// - Wrapper for private method that deals with actual MPI_Datatype
		template<typename T>
		int get_count() const {
			return get_count( get_MPI_Datatype<T>() );
		}

		const MPI_Status& get_status() const { return status_; };
		MPI_Status&       access_status()    { return status_; };

	 private:
		// Underlying MPI Status
		MPI_Status status_;

		int get_count(const MPI_Datatype& data_type) const;
	};


	// Pass as MpiStatus to obtain MPI_STATUS_IGNORE behavior
	static MpiStatus mpi_status_ignore;

	// MPI_Request wrapper
	class MpiRequest {
		MpiRequest() {}
		MpiRequest(const MPI_Request& request) 
		 : request_(request) {}

		// Wrapper for MPI_Wait
		// - Default: MPI_STATUS_IGNORE
		void wait(MpiStatus& mpi_status = mpi_status_ignore);

		const MPI_Request& get_request() const { return request_; };
		MPI_Request&       access_request()    { return request_; };

	 private:
		// Underlying MPI Request
		MPI_Request request_;
	};

	//------------------------------//
	//----- Misc. MPI Routines -----//
	//------------------------------//

	// TODO abort


	//-----------------------------//
	//----- MPI Data Wrappers -----//
	//-----------------------------//

	/* TODO 
	   Contiguous data types like std::array<T,dim>,
	   std::array<std::complex<T>,dim>, std::array<std::array<dim2>,dim1>,
	   and so on can all be handled within the same generalized framework.
     Then you only need one extra template for a std::vector of a 
	   ContiguousData/GeneralizedArray object

	template typename<ContainerType, PrimitiveType>
	struct GeneralizedArray {
		ContainerType<T> array;

		void* data_ptr,
		int 
		int num_primitive_elements;
	}
	*/


	// Wrapper around several variables that define an MPI buffer
	// - A "data element" is a single piece of data of type "data_type"
	// - The difference between "size" and "max_size" has the same semantics
	//   as the difference between std::vector's size and capacity
	// - TODO: Make class with public data? Would allow inheritance
	struct MpiData {
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


		//----- Constructors -----//

		// Default: empty
		MpiData(); 

		// From a reference to a primitive type
		template<typename T> explicit MpiData(T& data);

		// From a pointer to a single data element or array of primitive type
		// - Default: pointer to single element
		template<typename T> explicit MpiData(T* data_ptr_in, const int num_elements_in = 1);

		// From a std::vector of a primitive type
		template<typename T> MpiData(std::vector<T>& vec);

		// From a std::vector of std::complex of a primitive type
		template<typename T> MpiData(std::vector<std::complex<T>>& vec);

		// From a std::array of a primitive type 
		template<typename T, std::size_t dim> MpiData(std::array<T,dim>& arr);

		// From a std::vector of std::arrays of a primitive type
		template<typename T, std::size_t dim> MpiData(std::vector<std::array<T,dim>>& vec);

		// From a std::array of std::arrays of a primitive type
		// - dim1 = num_rows, dim2 = num_cols
		template<typename T, std::size_t dim1, std::size_t dim2> MpiData(
				std::array<std::array<T,dim2>, dim1>& matrix
		);
	}; // end struct MpiData

	// Pass as MpiData to perform MPI_IN_PLACE
	static MpiData mpi_in_place;


	//------------------------------------------------//
	//----- Point-to-Point Communication Methods -----//
	//------------------------------------------------//

	//----- Send -----//

	// Wrapper around underlying MPI_Send
	void send(
		const MpiData& data, 
		const int destination,  // rank of receiving process
		const int tag           // unique identifier (alternative to communicators)
	);

	// primitive type (by reference)
	template<typename T>
	void send(const T& value, const int destination, const int tag);

	// std::vector of a primitive type
	template<typename T>
	void send(const std::vector<T>& vec, const int destination, const int tag);

	// std::vector of std::arrays of a primitive type
	template<typename T, std::size_t dim>
	void send(const std::vector<std::array<T,dim>>& vec, const int destination, const int tag);


	//----- Recv -----//

	// Wrapper around underlying MPI_Send
	void recv(
		MpiData& data, 
		const int source, // rank of source
		const int tag, 
		MpiStatus& status = mpi_status_ignore
	);

	// primitive type (by reference)
	template<typename T>
	void recv(T& value, const int source, const int tag, 
	          MpiStatus& status = mpi_status_ignore);

	// std::vector of a primitive type
	template<typename T>
	void recv(std::vector<T>& vec, const int source, const int tag, 
	          MpiStatus& status = mpi_status_ignore);

	// std::vector of std::arrays of a primitive type
	template<typename T, std::size_t dim>
	void recv(std::vector<std::array<T,dim>>& vec, const int source, const int tag, 
	          MpiStatus& status = mpi_status_ignore);


	//--------------------------------------------//
	//----- Collective Communication Methods -----//
	//--------------------------------------------//

	// Wrapper for MPI_Barrier()
	void barrier();

	// Supported MPI_Op operations for allreduce
	// - To add a new MPI_Op, include it here and in the map_mpi_op() function below
	enum class MpiOp {
		null, max, min, sum, product
	};


	//----- Bcast -----//

	// Wrapper for MPI_Bcast
	// - "data" is copied from rank "root" to all ranks (including root!)
	void bcast(
		MpiData&  data, 
		const int root
	);

	// reference to a block of contiguous memory of a primitive type 
	template<typename T>
	void bcast(T& data, const int root);


	//----- Allreduce -----//

	// Wrapper around underlying MPI_Allreduce
	// - Pass MpiCommunicator::mpi_in_place as "data_in" to perform MPI_IN_PLACE
	void allreduce(
		const MpiData&  data_in, 
		MpiData&        data_out,
		const MpiOp& op
	);

	// reference to a block of contiguous memory of a primitive type
	// - This covers many relevant cases of interest (P = primitive type):
	//    - Single element, P&
	//    - std::vector<P>
	//    - std::vector<std::array<P>>
	// - The only requirement is that the data type has an MpiData constructor
	template<typename T>
	void allreduce(const T& data_in, T& data_out, const MpiOp& op);

	// Same as above function, but operation is performed in place
	template<typename T>
	void allreduce_in_place(T& data, const MpiOp& op);

	// Same as above function, but operation is summation (and it is performed in place)
	// - Defined for maximum laziness
	template<typename T>
	void allreduce_sum_in_place(T& data);

	// Same as above function, but operation is multiplication (and it is performed in place)
	// - Defined for maximum laziness
	template<typename T>
	void allreduce_multiply_in_place(T& data);


	//----- Allgather -----//

	// Wrapper around underlying MPI_Allgather
	// - Amount of data to send/receive from each rank is inferred from the size
	//   of data_to_send/received_data
	void allgather(
		const MpiData& data_to_send,
		MpiData&       received_data
	);
		
	// Receive a single value of primitive type from each rank and put the results
	// into a std::vector
	template<typename T>
	void allgather(
		const T&        value,
		std::vector<T>& received_values
	);


	//----- Allgatherv -----//

	// Wrapper around underlying MPI_Allgatherv  TODO
	// - Assumes that appropriate blocks have already been prepared, that that 
	//   "received_data" is the correct size to receive it all
	// - Variables
	//     received_data: buffer for all data received
	//     block_offsets: index in "received_data.data_ptr" at which each block of data starts
	//     block_sizes: number of elements received
	void allgatherv(
		const MpiData&          data,
		const std::vector<int>& block_offsets,
		const std::vector<int>& block_sizes,
		MpiData&                received_data
	);

	// std::vectors of a primitive type 
	// - Calls this->allgather() first to determine the number of data elements
	//   to receive from each rank
	template<typename T>
	void allgatherv(
		const std::vector<T>& data,
		std::vector<T>&       received_data,
		std::vector<int>&     block_offsets,
		std::vector<int>&     block_sizes
	);

	// std::vectors of std::arrays of a primitive type 
	template<typename T, std::size_t dim>
	void allgatherv(
		const std::vector<std::array<T,dim>>& data,
		std::vector<std::array<T,dim>>&       received_data,
		std::vector<int>&                     block_offsets,
		std::vector<int>&                     num_arrays_per_block
	);


	//----------------------//
	//----- Exceptions -----//
	//----------------------//

	class MpiNotEnabledException : public std::exception {
	 public:
		MpiNotEnabledException() {}
		const char* what() const noexcept override {
			return "MPI is not enabled";
		};
	};

	class MpiNotInitializedException : public std::exception {
	 public:
		MpiNotInitializedException() {}
		const char* what() const noexcept override {
			return "MPI has not been initialized";
		};
	};

	class MpiOpNotFoundException : public std::exception {
	 public:
		MpiOpNotFoundException() {}
		const char* what() const noexcept override {
			return "Given MpiOp not found";
		};
	};


 private:
	// Underlying MPI communicator
	bool is_communicator_initialized_ = false;
	MPI_Comm communicator_;

	// Rank (index) of the master process
	const int master_rank_ = 0;

	//

	// Manually generate the map from C++ enum to MPI_Op since std::map
	// would require a custom comparator function, and the underlying 
	// MPI_Op type is not portable
	// - Throws an exception if the desired operator is not available
	MPI_Op map_mpi_op(const MpiOp op) const;
}; // end class MpiCommunicator


//--------------------------------//
//--------------------------------//
//----- TEMPLATE DEFINITIONS -----//
//--------------------------------//
//--------------------------------//

// Templates... templates everywhere...


// Wrapper for MPI_Dims_create
template<std::size_t dim>
std::array<int,dim> MpiCommunicator::calculateGridDimensions() const
{
	if ( dim == 0 ) {
		throw std::runtime_error("MPI: Can't parition ranks along 0 dimensions.\n");
	}

	std::array<int,dim> grid_dimensions;

	// Use an int to avoid compiler warnings about comparing signed and
	// unsigned integers
	int dim_int = static_cast<int>(dim);

#ifdef MPI_ENABLED
	// Use MPI_Dimms_create to auto-select grid dimensions
	int num_ranks = this->get_size();
	for ( int d=0; d<dim_int; ++d ) {
		// MPI_Dimms create will only change grid_dimensions[d] if it's 0; any
		// positive values will be unchanged!
		grid_dimensions[d] = 0;
	}
	MPI_Dims_create( num_ranks, dim_int, grid_dimensions.data() );

#else
	for ( int d=0; d<dim_int; ++d ) {
		grid_dimensions[d] = 1;
	}
#endif /* MPI_ENABLED */

	return grid_dimensions;
}



//-----------------------------//
//----- MPI Data Wrappers -----//
//-----------------------------//

//----- Constructors -----//

// From a reference to a primitive type
template<typename T>
MpiCommunicator::MpiData::MpiData(T& data)
 : data_ptr( static_cast<void*>(&data) ), 
   max_size(1),
   size(1), 
   data_type( get_MPI_Datatype<T>() )
{}


// From a raw pointer to a primitive type
template<typename T>
MpiCommunicator::MpiData::MpiData(T* data_ptr_in, const int num_elements_in)
 : data_ptr( static_cast<void*>(data_ptr_in) ), 
   max_size(num_elements_in),
   size(num_elements_in), 
   data_type( get_MPI_Datatype<T>() )
{
	if ( data_ptr == nullptr ) {
		size = 0;
	}
}


// From a std::vector of a primitive type
template<typename T> 
MpiCommunicator::MpiData::MpiData(std::vector<T>& vec)
 : data_ptr( static_cast<void*>(vec.data()) ),
   max_size(vec.capacity()), 
   size(vec.size()), 
   data_type( get_MPI_Datatype<T>() )
{}


// From a std::vector of std::complex of a primitive type
template<typename T> 
MpiCommunicator::MpiData::MpiData(std::vector<std::complex<T>>& vec)
 : data_ptr( static_cast<void*>(vec.data()) ),
   max_size( 2*vec.capacity() ), 
   size( 2*vec.size() ), 
   data_type( get_MPI_Datatype<T>() )
{}


// From a std::array of a primitive type 
template<typename T, std::size_t dim> 
MpiCommunicator::MpiData::MpiData(std::array<T,dim>& arr)
 : data_ptr( static_cast<void*>(arr.data()) ),
   max_size(dim),
   size(dim),
   data_type( get_MPI_Datatype<T>() )
{
	if ( dim == 0 ) {
		throw std::runtime_error("Can't construct MpiData from an array with dim = 0");
	}

	// I don't think this is an issue, but just to be safe ...
	if ( data_ptr == nullptr ) {
		throw std::runtime_error("pointer to data (from std::vector) is null");
	}
}


// From a std::vector of std::arrays of a primitive type
template<typename T, std::size_t dim> 
MpiCommunicator::MpiData::MpiData(std::vector<std::array<T,dim>>& vec)
 : data_ptr( static_cast<void*>(vec.data()) ),
   max_size( dim*vec.capacity() ), 
   size( dim*vec.size() ), 
   data_type( get_MPI_Datatype<T>() )
{
	if ( dim == 0 ) {
		throw std::runtime_error("Can't construct MpiData from a vector of arrays with dim = 0");
	}
}


// From a std::array of std::arrays of a primitive type
// - dim1 = num_rows, dim2 = num_cols
template<typename T, std::size_t dim1, std::size_t dim2> 
MpiCommunicator::MpiData::MpiData(std::array<std::array<T,dim2>, dim1>& matrix)
 : max_size( dim1*dim2 ),
   size( dim1*dim2 ), 
   data_type( get_MPI_Datatype<T>() )
{
	if ( dim1 == 0 or dim2 == 0 ) {
		throw std::runtime_error("Can't construct MpiData from a matrix with dim1 = 0 or dim2 = 0");
	}

	// Safe to take the address of the first element
	data_ptr = static_cast<void*>( &matrix[0][0] );
}


//------------------------------------------------//
//----- Point-to-Point Communication Methods -----//
//------------------------------------------------//


//----- Send -----//

// primitive type by reference
template<typename T>
void MpiCommunicator::send(const T& value, const int destination, const int tag)
{
	// MpiData constructor will not actually change the data
	const MpiData data_to_send( const_cast<T&>(value) );
	this->send(data_to_send, destination, tag);
}


// std::vector of a primitive type
template<typename T>
void MpiCommunicator::send(
		const std::vector<T>& vec, const int destination, const int tag)
{
	// MpiData constructor will not actually change the data
	const MpiData data_to_send( const_cast<std::vector<T>&>(vec) );
	this->send(data_to_send, destination, tag);
}


// std::vector of std::arrays of a primitive type
template<typename T, std::size_t dim>
void MpiCommunicator::send(
		const std::vector<std::array<T,dim>>& vec, const int destination, const int tag)
{
	// MpiData constructor will not actually change the data
	auto non_const_vec_ref = const_cast< std::vector<std::array<T,dim>>& >(vec);
	const MpiData data_to_send( non_const_vec_ref );
	this->send(data_to_send, destination, tag);
}


//----- Recv -----//

// primitive type by reference
template<typename T>
void MpiCommunicator::recv(
		T& value, const int source, const int tag, 
		MpiCommunicator::MpiStatus& status)
{
	MpiData data_to_recv(value);
	this->recv(data_to_recv, source, tag, status);
}


// std::vector of a primitive type
template<typename T>
void MpiCommunicator::recv(
		std::vector<T>& vec, const int source, const int tag, 
		MpiCommunicator::MpiStatus& status)
{
	// For safety, extend the buffer to its maximum capacity
	vec.resize( vec.capacity() );
	MpiData data_to_recv(vec);
	this->recv(data_to_recv, source, tag, status);

	// If the MPI_Status is not being ignored, set the vector's final size based on 
	// how much data was actually received
	// - Else it must be assumed that the entire vector was filled
	if ( &status != &mpi_status_ignore ) {
		data_to_recv.size = status.get_count<T>();
		vec.resize( data_to_recv.size );
	}
}


// std::vector of std::arrays of a primitive type
template<typename T, std::size_t dim>
void MpiCommunicator::recv(
		std::vector<std::array<T,dim>>& vec, const int source, const int tag, 
		MpiCommunicator::MpiStatus& status)
{
	// For safety, extend the buffer to its maximum capacity
	vec.resize( vec.capacity() );
	MpiData data_to_recv(vec);
	this->recv(data_to_recv, source, tag, status);

	// If the MPI_Status is not being ignored, set the vector's final size based on 
	// how much data was actually received
	// - Else it must be assumed that the entire vector was filled
	if ( &status != &mpi_status_ignore ) {
		data_to_recv.size = status.get_count<T>();
		if ( data_to_recv.size % dim == 0 ) {
			vec.resize( data_to_recv.size/dim );
		}
		else {
			throw std::runtime_error("MPI recv obtained a partial array");
		}
	}
}


//------------------------------------//
//----- Collective Communication -----//
//------------------------------------//

//----- Bcast -----//

// reference to a block of contiguous memory of a primitive type 
template<typename T>
void MpiCommunicator::bcast(T& data, const int root)
{
	MpiData mpi_data( data );
	this->bcast(mpi_data, root);
}


//----- Allreduce -----//

// reference to a block of contiguous memory of a primitive type 
template<typename T>
void MpiCommunicator::allreduce(
		const T& data_in, T& data_out, const MpiOp& op)
{
	MpiData mpi_data_in( const_cast<T&>(data_in) );
	MpiData mpi_data_out( data_out );
	this->allreduce(mpi_data_in, mpi_data_out, op);
}


// reference to a block of contiguous memory of a primitive type 
template<typename T>
void MpiCommunicator::allreduce_in_place(
		T& data, const MpiOp& op)
{
	MpiData mpi_data(data);
	this->allreduce(mpi_in_place, mpi_data, op);
}


// reference to a block of contiguous memory of a primitive type 
template<typename T>
void MpiCommunicator::allreduce_sum_in_place(T& data)
{
	MpiData mpi_data(data);
	this->allreduce(mpi_in_place, mpi_data, MpiOp::sum);
}


// reference to a block of contiguous memory of a primitive type 
template<typename T>
void MpiCommunicator::allreduce_multiply_in_place(T& data)
{
	MpiData mpi_data(data);
	this->allreduce(mpi_in_place, mpi_data, MpiOp::product);
}


//----- Allgather -----//


// Receive a single value of primitive type from each rank and put the results
// into a std::vector
template<typename T>
void MpiCommunicator::allgather(const T& value, std::vector<T>& received_values)
{
	// Value to send will not actually be modified
	MpiData data_to_send( const_cast<T&>(value) );

	// Allocate space in the receive buffer
	int num_ranks = this->get_size();
	received_values.resize(num_ranks);
	MpiData received_data(received_values);

	this->allgather(data_to_send, received_data);
}


//----- Allgatherv -----//


// From std::vectors of a primitive type 
template<typename T>
void MpiCommunicator::allgatherv(
	const std::vector<T>& data,
	std::vector<T>& received_data, std::vector<int>& block_offsets, std::vector<int>& block_sizes)
{
	// First, gather the number of elements that each rank will send.
	// This is required by MPI_Allgatherv.
	int num_ranks = this->get_size();
	block_sizes.assign(num_ranks, 0);
	int local_block_size = data.size();
	this->allgather(local_block_size, block_sizes);

	// Set up local receive buffer
	int buffer_size = 0;
	block_offsets.resize(num_ranks);
	for ( int r=0; r<num_ranks; ++r ) {
		// Start of each block
		if ( r == 0 )
			block_offsets[r] = 0;
		else {
			block_offsets[r] = block_offsets[r-1] + block_sizes[r-1];
		}

		// Total number of elements of primitive type, across all ranks
		buffer_size += block_sizes[r];
	}

	// Prepare receive buffer
	received_data.resize(buffer_size);

	// Gather the data
	MpiData mpi_data_to_send( const_cast<std::vector<T>&>(data) );
	MpiData mpi_received_data( received_data );
	this->allgatherv(mpi_data_to_send, block_offsets, block_sizes,
	                 mpi_received_data);
}


// From std::vectors of std::arrays of a primitive type TODO
template<typename T, std::size_t dim>
void MpiCommunicator::allgatherv(
	const std::vector<std::array<T,dim>>& data,
	std::vector<std::array<T,dim>>& received_data, 
	std::vector<int>& block_offsets, std::vector<int>& num_arrays_per_block)
{
	// First, gather the number of arrays that each rank will send.
	// This is required by MPI_Allgatherv.
	int num_ranks = this->get_size();
	num_arrays_per_block.assign(num_ranks, 0);
	int local_num_arrays = data.size();
	this->allgather(local_num_arrays, num_arrays_per_block);

	// Set up the output arrays
	int num_arrays_total = 0;
	block_offsets.resize(num_ranks);
	for ( int r=0; r<num_ranks; ++r ) {
		// Start of each block
		if ( r == 0 )
			block_offsets[r] = 0;
		else {
			block_offsets[r] = block_offsets[r-1] + num_arrays_per_block[r-1];
		}

		// Total number of elements of primitive type, across all ranks
		num_arrays_total += num_arrays_per_block[r];
	}

	// MPI_Allgatherv expects non-const input (data will not actually be changed)
	auto data_to_send_ref = const_cast< std::vector<std::array<T,dim>>& >(data);
	MpiData mpi_data_to_send( data_to_send_ref );

	// Prepare receive buffer
	// - Resizing the output array ensures that the MpiData's size is correctly set
	//   based on the array's dimension
	received_data.resize(num_arrays_total);
	MpiData mpi_received_data( received_data );

	// MPI expects a buffer of contiguous data elements of a primitive type
	std::vector<int> mpi_block_offsets(num_ranks), mpi_block_sizes(num_ranks);
	int dim_int = static_cast<int>(dim);
	for ( int r=0; r<num_ranks; ++r ) {
		mpi_block_offsets[r] = dim_int*block_offsets[r];
		mpi_block_sizes[r]   = dim_int*num_arrays_per_block[r];
	}

	// Gather the data
	this->allgatherv(mpi_data_to_send, mpi_block_offsets, mpi_block_sizes,
	                 mpi_received_data);
}


#endif /* MPI_COMMUNICATOR_H */
