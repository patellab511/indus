// MpiCommunicator.cpp
// - Mainly implementations of non-templated wrappers around MPI functions

#include "MpiCommunicator.h"
//#include "MpiDatatype.h"

// Useful typedefs
using MpiStatus = MpiCommunicator::MpiStatus;

#ifndef MPI_ENABLED
/*
// Dummy MPI_Init and MPI_Finalize for when MPI isn't available
void MPI_Init(int* argc, char** argv[]) {
	// Suppress "unused arguments" warnings
	(void) argc; (void) argv;
};

// Dummy MPI_Finalize
void MPI_Finalize() {};
*/
#endif /* MPI_ENABLED */


//-------------------------------------------//
//----- Constructors and Initialization -----//
//-------------------------------------------//

namespace MpiEnvironment {

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

} // end namespace MpiEnvironment


MpiCommunicator::MpiCommunicator()
{
#ifdef MPI_ENABLED
	//do_MPI_Init_if_not_initialized();
	set_communicator(MPI_COMM_SELF);
	register_default_MpiDatatypes();
#endif /* MPI_ENABLED */
}


MpiCommunicator::MpiCommunicator(const MPI_Comm& communicator)
{
#ifdef MPI_ENABLED
	//do_MPI_Init_if_not_initialized();
	set_communicator(communicator);
	register_default_MpiDatatypes();
#else
	// Suppress "unused arguments" warnings
	(void) communicator;
#endif /* MPI_ENABLED */
}


void MpiCommunicator::register_default_MpiDatatypes()
{
	// First, register MPI primitives. These will be used to build derived MPI_Datatypes.
	register_primitive_MPI_Datatypes();

	// Register a few commonly-used derived types
	register_type< std::array<int,1> >();
	register_type< std::array<int,2> >();
	register_type< std::array<int,3> >();
	register_type< std::array<float,1> >();
	register_type< std::array<float,2> >();
	register_type< std::array<float,3> >();
	register_type< std::array<double,1> >();
	register_type< std::array<double,2> >();
	register_type< std::array<double,3> >();
	register_type< std::complex<float> >();
	register_type< std::complex<double> >();

	register_type< std::array<std::array<double,1>,1> >();
	register_type< std::array<std::array<double,2>,2> >();
	register_type< std::array<std::array<double,3>,3> >();
}


// Registers all predefined MPI_Datatypes in 'mpi_datatype_map_'
void MpiCommunicator::register_primitive_MPI_Datatypes()
{
	// TODO mappings for MPI_Byte and MPI_Packed
	register_type<float>();
	register_type<double>();
	register_type<long double>();

	register_type<int>();
	register_type<long int>();
	register_type<short int>();
	register_type<unsigned int>();
	register_type<unsigned long int>();
	register_type<unsigned short int>();

	register_type<char>();
	register_type<unsigned char>();
}


// TODO Copy other member variables
MpiCommunicator::MpiCommunicator(const MpiCommunicator& communicator)
{
#ifdef MPI_ENABLED
	set_communicator(communicator.communicator_);
#else
	// Suppress "unused arguments" warnings
	(void) communicator;
#endif /* MPI_ENABLED */
}


bool MpiCommunicator::is_mpi_enabled()
{
	return MpiEnvironment::is_enabled();
}


bool MpiCommunicator::is_mpi_initialized()
{
	return MpiEnvironment::is_initialized();
}


MpiCommunicator::~MpiCommunicator() 
{
#ifdef MPI_ENABLED
	// Deallocate the duplicated communicator
	MPI_Comm_free(&communicator_);
#endif /* MPI_ENABLED */
}

// Returns MPI_COMM_WORLD if MPI is enabled, else returns a dummy MPI_Comm
MPI_Comm MpiCommunicator::get_mpi_comm_world()
{
#ifdef MPI_ENABLED
	return MPI_COMM_WORLD;
#else
	return MPI_Comm();
#endif
}

// Returns MPI_COMM_SELF if MPI is enabled, else returns a dummy MPI_Comm
MPI_Comm MpiCommunicator::get_mpi_comm_self()
{
#ifdef MPI_ENABLED
	return MPI_COMM_SELF;
#else
	return MPI_Comm();
#endif
}


//-------------------------------------------------//
//----- Get/Set Functions for MpiCommunicator -----//
//-------------------------------------------------//


void MpiCommunicator::set_communicator(const MPI_Comm& communicator)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		if ( is_communicator_initialized_ ) {
			// Free the old communicator
			MPI_Comm_free(&communicator_);
		}

		MPI_Comm_dup(communicator, &communicator_);
		is_communicator_initialized_ = true;
	}
	else {
		throw MpiNotInitializedException();
	}

#else
	// Set dummy communicator
	communicator_ = communicator;

	/*
	// Suppress "unused function" warnings
	int argc_dummy = 0;
	char** argv_dummy = nullptr;
	MPI_Init(&argc_dummy, &argv_dummy);
	MPI_Finalize();
	*/
#endif /* MPI_ENABLED */
}


void MpiCommunicator::set_communicator(void* communicator_ptr)
{
	if ( communicator_ptr != nullptr ) {
		set_communicator( *(static_cast<MPI_Comm*>(communicator_ptr)) );
	}
	else {
		throw std::runtime_error("communicator pointer is null");
	}
}


int MpiCommunicator::get_rank() const
{
	int rank = 0;
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) { 
		MPI_Comm_rank(communicator_, &rank);
	}
	else {
		throw MpiNotInitializedException(); 
	}
#endif /* MPI_ENABLED */
	return rank;
}


int MpiCommunicator::get_size() const
{
	int size = 1;
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) { 
		MPI_Comm_size(communicator_, &size);
	}
	else {
		throw MpiNotInitializedException(); 
	}
#endif /* MPI_ENABLED */
	return size;
}


//-------------------------------------------//
//----- MPI Status and Request Wrappers -----//
//-------------------------------------------//


int MpiStatus::get_count(const MPI_Datatype& data_type) const
{
	int count = 0;
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) { 
		// Perform const_cast since MPI_Get_count() expects a non-const MPI_Status
		MPI_Get_count(const_cast<MPI_Status*>(&status_), data_type, &count);
	}
	else {
		throw MpiNotInitializedException();  
	}
#else
	(void) data_type;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
	return count;
}


// Definition of static variable
MpiStatus MpiCommunicator::mpi_status_ignore;


void MpiCommunicator::MpiRequest::wait(MpiStatus& mpi_status)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		if ( &mpi_status == &(MpiCommunicator::mpi_status_ignore) ) {
			// Using "alias" for MPI_STATUS_IGNORE
			MPI_Wait(&request_, MPI_STATUS_IGNORE);
		}
		else {
			MPI_Wait(&request_, &(mpi_status.access_MPI_Status()));
		}
	}
	else {
		throw MpiNotInitializedException(); 
	}
#else
	(void) mpi_status;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


//-----------------------------//
//----- MPI Data Wrappers -----//
//-----------------------------//

// Definition of static variable
const MpiData MpiCommunicator::mpi_in_place;


//------------------------------------------------//
//----- Point-to-Point Communication Methods -----//
//------------------------------------------------//


// Wrapper for MPI_Send
void MpiCommunicator::send(
	const MpiData& data, const int destination, const int tag)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		// MPI_Send expects non-const arguments
		auto mpi_data = const_cast<MpiData&>(data);
		MPI_Send(mpi_data.data_ptr, mpi_data.size, mpi_data.data_type, //.get_MPI_Datatype(),
		         destination, tag, this->communicator_);
	}	
	else { 
		throw MpiNotInitializedException(); 
	}
#else
	(void) data; (void) destination; (void) tag;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper for MPI_Recv
void MpiCommunicator::recv(
	MpiData& data, const int source, const int tag, MpiStatus& status)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		if ( &status == &mpi_status_ignore ) {
			// Ignore status
			MPI_Recv(data.data_ptr, data.size, data.data_type, //.get_MPI_Datatype(),
							 source, tag, this->communicator_, MPI_STATUS_IGNORE);
		}
		else {
			MPI_Recv(data.data_ptr, data.max_size, data.data_type, //.get_MPI_Datatype(),
							 source, tag, this->communicator_, &status.access_MPI_Status());
		}
	}
	else {
		throw MpiNotInitializedException(); 
	}	
#else
	(void) data; (void) source; (void) tag; (void) status;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper for MPI_Isend
void MpiCommunicator::isend(
	const MpiData& data, const int destination, const int tag, MpiRequest& request)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		MPI_Isend(data.data_ptr, data.size, data.data_type, //.get_MPI_Datatype(),
		          destination, tag, this->communicator_, &request.access_MPI_Request());
	}
	else {
		throw MpiNotInitializedException(); 
	}	
#else
	(void) data; (void) destination; (void) tag; (void) request;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper for MPI_Irecv
void MpiCommunicator::irecv(
	MpiData& data, const int source, const int tag, MpiRequest& request)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		MPI_Irecv(data.data_ptr, data.size, data.data_type, //.get_MPI_Datatype(),
		          source, tag, this->communicator_, &request.access_MPI_Request());
	}
	else {
		throw MpiNotInitializedException(); 
	}	
#else
	(void) data; (void) source; (void) tag; (void) request;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


//------------------------------------//
//----- Collective Communication -----//
//------------------------------------//


// Wrapper for MPI_Barrier
void MpiCommunicator::barrier() 
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		MPI_Barrier(this->communicator_);
	}
	else { 
		throw MpiNotInitializedException(); 
	}
#else
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper for MPI_Bcast
void MpiCommunicator::bcast(MpiData& data, const int root)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		MPI_Bcast(data.data_ptr, data.size, data.data_type, //.get_MPI_Datatype(),
		          root, this->communicator_);
	}
	else { 
		throw MpiNotInitializedException(); 
	}
#else
	(void) data; (void) root;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper for MPI_Allreduce
void MpiCommunicator::allreduce(
	const MpiData& data_in, MpiData& data_out, const MpiOp& op)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		if ( &data_in == &(MpiCommunicator::mpi_in_place) ) {
			// Perform allreduce in place
			MPI_Allreduce(MPI_IN_PLACE, data_out.data_ptr, data_out.size, data_out.data_type, //.get_MPI_Datatype(),
			              op.mpi_op_, this->communicator_);
		}
		else {
			// Perform allreduce and place result in output buffer
			// - MPI_Allreduce expects non-constant inputs
			auto data_in_ref = const_cast<MpiData&>(data_in);
			if ( data_in_ref.data_type == data_out.data_type and
			     data_in_ref.size      == data_out.size 
			) {
				MPI_Allreduce(data_in_ref.data_ptr, data_out.data_ptr, data_out.size, data_out.data_type, //.get_MPI_Datatype(),
											op.mpi_op_, this->communicator_);
			}
			else {
				throw std::runtime_error("allreduce: data type and/or size mismatch");
			}
		}
	}
	else { 
		throw MpiNotInitializedException(); 
	}
#else
	(void) data_in; (void) data_out; (void) op;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper around underlying MPI_Allgather
void MpiCommunicator::allgather(
	const MpiData& data_to_send,
	MpiData& received_data)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		// Check input
		if ( data_to_send.data_type != received_data.data_type ) {
			throw std::runtime_error("Error in MPI allgather: data type mismatch");
		}

		// MPI_Allgather expects non-const inputs
		auto data_to_send_ref = const_cast<MpiData&>(data_to_send);

		MPI_Allgather(data_to_send_ref.data_ptr, data_to_send_ref.size, data_to_send_ref.data_type, //.get_MPI_Datatype(),
		              received_data.data_ptr,    data_to_send_ref.size, data_to_send_ref.data_type, //.get_MPI_Datatype(),
		              this->communicator_);
	}
	else {
		throw MpiNotInitializedException(); 
	}
#else
	(void) data_to_send; (void) received_data;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}


// Wrapper around underlying MPI_Allgatherv 
void MpiCommunicator::allgatherv(
		const MpiData& data, const std::vector<int>& block_offsets, const std::vector<int>& block_sizes,
		MpiData& received_data)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		// MPI_Allgatherv expects non-constant send buffer
		auto data_ref = const_cast<MpiData&>(data);
		MPI_Allgatherv(data_ref.data_ptr, data_ref.size, data_ref.data_type, //.get_MPI_Datatype(),
		               received_data.data_ptr, block_sizes.data(), block_offsets.data(),
		               received_data.data_type /*.get_MPI_Datatype()*/, this->communicator_);
	}
	else { 
		throw MpiNotInitializedException(); 
	}
#else
	(void) data; (void) block_offsets; (void) block_sizes; (void) received_data;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}
