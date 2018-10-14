// MpiCommunicator.cpp
// - Mainly implementations of non-templated wrappers around MPI functions

#include "MpiCommunicator.h"


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


MpiCommunicator::MpiCommunicator()
{
#ifdef MPI_ENABLED
	//do_MPI_Init_if_not_initialized();
	set_communicator(MPI_COMM_SELF);
#endif /* MPI_ENABLED */
}


MpiCommunicator::MpiCommunicator(const MPI_Comm& communicator)
{
#ifdef MPI_ENABLED
	//do_MPI_Init_if_not_initialized();
	set_communicator(communicator);
#else
	// Suppress "unused arguments" warnings
	(void) communicator;
#endif /* MPI_ENABLED */
}


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
#ifdef MPI_ENABLED
	return true;
#else
	return false;
#endif /* MPI_ENABLED */
}


bool MpiCommunicator::is_mpi_initialized()
{
	int is_init = false;
#ifdef MPI_ENABLED
	MPI_Initialized(&is_init);
#endif /* MPI_ENABLED */
	if ( is_init ) { return true;  }
	else           { return false; }
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


//-----------------------------//
//----- Get/Set Functions -----//
//-----------------------------//


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


#ifdef MPI_ENABLED
// Mapping from C++ primitives to MPI_Datatypes (as template specializations; clever, PLUMED team!)
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<float>()         { return MPI_FLOAT;         }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<double>()        { return MPI_DOUBLE;        }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<int>()           { return MPI_INT;           }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<char>()          { return MPI_CHAR;          }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<unsigned>()      { return MPI_UNSIGNED;      }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<long unsigned>() { return MPI_UNSIGNED_LONG; }
#else
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<float>()         { return MPI_Datatype(); }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<double>()        { return MPI_Datatype(); }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<int>()           { return MPI_Datatype(); }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<char>()          { return MPI_Datatype(); }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<unsigned>()      { return MPI_Datatype(); }
template<> MPI_Datatype MpiCommunicator::get_MPI_Datatype<long unsigned>() { return MPI_Datatype(); }
#endif /* MPI_ENABLED */


//------------------------//
//----- MPI Wrappers -----//
//------------------------//


int MpiCommunicator::MpiStatus::get_count(const MPI_Datatype& data_type) const
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
MpiCommunicator::MpiStatus MpiCommunicator::mpi_status_ignore;


void MpiCommunicator::MpiRequest::wait(MpiStatus& mpi_status)
{
#ifdef MPI_ENABLED
	if ( is_mpi_initialized() ) {
		if ( &mpi_status == &(MpiCommunicator::mpi_status_ignore) ) {
			// Using "alias" for MPI_STATUS_IGNORE
			MPI_Wait(&request_, MPI_STATUS_IGNORE);
		}
		else {
			MPI_Wait(&request_, &(mpi_status.access_status()));
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

MpiCommunicator::MpiData::MpiData()
 : data_ptr(nullptr), max_size(0), size(0), 
#ifdef MPI_ENABLED
	data_type(MPI_DATATYPE_NULL) 
#else
	data_type()
#endif /* MPI_ENABLED */
{}


// Definition of static variable
MpiCommunicator::MpiData MpiCommunicator::mpi_in_place;


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
		MPI_Send(mpi_data.data_ptr, mpi_data.size, mpi_data.data_type, destination, 
		         tag, this->communicator_);
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
			MPI_Recv(data.data_ptr, data.max_size, data.data_type, source,
							 tag, this->communicator_, MPI_STATUS_IGNORE);
		}
		else {
			MPI_Recv(data.data_ptr, data.max_size, data.data_type, source,
							 tag, this->communicator_, &status.access_status());
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


//------------------------------------//
//----- Collective Communication -----//
//------------------------------------//


MPI_Op MpiCommunicator::map_mpi_op(const MpiOp op) const
{
#ifdef MPI_ENABLED
	if      ( op == MpiOp::null     ) { return MPI_OP_NULL; }
	else if ( op == MpiOp::max      ) { return MPI_MAX;  }
	else if ( op == MpiOp::min      ) { return MPI_MIN;  }
	else if ( op == MpiOp::sum      ) { return MPI_SUM;  }
	else if ( op == MpiOp::product  ) { return MPI_PROD; }
	else {
		throw MpiOpNotFoundException();
	}
#else
	(void) op;
	throw MpiNotEnabledException();
	return MPI_Op(); // prevent compiler from complaining
#endif /* MPI_ENABLED */
};


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
		MPI_Bcast(data.data_ptr, data.size, data.data_type, root,
		          this->communicator_);
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

		// Get underlying MPI_Op
		MPI_Op mpi_op = map_mpi_op(op);

		if ( &data_in == &(MpiCommunicator::mpi_in_place) ) {
			// Perform allreduce in place
			MPI_Allreduce(MPI_IN_PLACE, data_out.data_ptr, data_out.size, data_out.data_type, 
			              mpi_op, this->communicator_);
		}
		else {
			// Perform allreduce and place result in output buffer
			// - MPI_Allreduce expects non-constant inputs
			auto data_in_ref = const_cast<MpiData&>(data_in);
			if ( data_in_ref.data_type == data_out.data_type and
			     data_in_ref.size      == data_out.size 
			) {
				MPI_Allreduce(data_in_ref.data_ptr, data_out.data_ptr, data_out.size, data_out.data_type, 
											mpi_op, this->communicator_);
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

		MPI_Allgather(data_to_send_ref.data_ptr, data_to_send_ref.size, data_to_send_ref.data_type,
		              received_data.data_ptr,    data_to_send_ref.size, data_to_send_ref.data_type,
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
		MPI_Allgatherv(data_ref.data_ptr, data_ref.size, data_ref.data_type,
		               received_data.data_ptr, block_sizes.data(), block_offsets.data(),
		               received_data.data_type, this->communicator_);
	}
	else { 
		throw MpiNotInitializedException(); 
	}
#else
	(void) data; (void) received_data; (void) received_counts;
	throw MpiNotEnabledException();
#endif /* MPI_ENABLED */
}
