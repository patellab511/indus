#ifndef MPI_DATATYPE_REGISTRAR_H
#define MPI_DATATYPE_REGISTRAR_H

// Registar for primitive types, and ones that support standard arithmetic operations
//, typename std::enable_if< std::is_arithmetic<T>::value >::type* = nullptr  // FIXME
//, class = typename std::enable_if< std::is_arithmetic<T>::value >
//, typename std::enable_if< std::is_arithmetic<T>::value, T >::type = 0
//, typename std::enable_if< std::is_arithmetic<T>::value, T >::type* = nullptr
//typename std::enable_if< std::is_integral<T>::value,   T >::type* = nullptr  // example
template<typename T>
class MpiDatatypeRegistrar 
{
 public:
	MpiDatatypeRegistrar() = delete;  // This class is only used for its static templated methods
	static MpiDatatype makeMpiDatatype(MpiCommunicator& comm);
	static void registerMpiOps(MpiOp::StandardOpsMap& map);

	static const bool mpi_sum_commutes = true;
	static void mpi_sum(T* invec, T* inoutvec, int* len, MPI_Datatype* data_type_ptr)
	{
		for ( int i=0; i<(*len); ++i ) {
			*inoutvec += *invec;
			++inoutvec;  ++invec;
		}
	}
};


// Helper class for registering std::arrays
//, typename std::enable_if< std::is_arithmetic<T>::value >::type* = nullptr  // FIXME
//, class = typename std::enable_if< std::is_arithmetic<T>::value >
//, typename std::enable_if< std::is_arithmetic<T>::value, T >::type = 0
//, class = typename std::enable_if< std::is_arithmetic<T>::value, T >::type
template<typename T, std::size_t dim>
class MpiDatatypeRegistrar<std::array<T,dim>> 
{
 public:
	MpiDatatypeRegistrar() = delete;  // This class is only used for its static templated methods

	static MpiDatatype makeMpiDatatype(MpiCommunicator& comm);
	static void registerMpiOps(MpiOp::StandardOpsMap& map);

	static const bool mpi_sum_commutes = true;
	static void mpi_sum(std::array<T,dim>* invec, std::array<T,dim>* inoutvec, 
											int* len, MPI_Datatype* data_type_ptr)
	{
		for ( int i=0; i<(*len); ++i ) {
			for ( unsigned a=0; a<dim; ++a ) {
				(*inoutvec)[a] += (*invec)[a];
			}
			++inoutvec;  ++invec;
		}
	}

	// TODO What should mpi_prod be for array types? element-wise?
};


// Helper class for registering std::array<std::array>
template<typename T, std::size_t dim1, std::size_t dim2>
class MpiDatatypeRegistrar<std::array<std::array<T, dim2>, dim1>>
{
 public:
	using Matrix = std::array<std::array<T, dim2>, dim1>;

	MpiDatatypeRegistrar() = delete;  // This class is only used for its static templated methods

	static MpiDatatype makeMpiDatatype(MpiCommunicator& comm);
	static void registerMpiOps(MpiOp::StandardOpsMap& map);

	static const bool mpi_sum_commutes = true;
	static void mpi_sum(Matrix* invec, Matrix* inoutvec, 
											int* len, MPI_Datatype* data_type_ptr)
	{
		for ( int i=0; i<(*len); ++i ) {
			for ( unsigned a=0; a<dim1; ++a ) {
				for ( unsigned b=0; b<dim2; ++b ) {
					(*inoutvec)[a][b] += (*invec)[a][b];
				}
			}
			++inoutvec;  ++invec;
		}
	}

	// TODO What should mpi_prod be for matrix types? element-wise?
};


// Helper class for registering std::complex
template<typename T>
class MpiDatatypeRegistrar<std::complex<T>>
{
 public:
	MpiDatatypeRegistrar() = delete;  // This class is only used for its static templated methods

	static MpiDatatype makeMpiDatatype(MpiCommunicator& comm);
	static void registerMpiOps(MpiOp::StandardOpsMap& map);

	static const bool mpi_sum_commutes = true;
	static void mpi_sum(std::complex<T>* invec, std::complex<T>* inoutvec, 
											int* len, MPI_Datatype* data_type_ptr)
	{
		for ( int i=0; i<(*len); ++i ) {
			*inoutvec += *invec;
			++inoutvec;  ++invec;
		}
	}

	static const bool mpi_prod_commutes = false;  // complex multiplication is not comutative
	static void mpi_prod(std::complex<T>* invec, std::complex<T>* inoutvec, 
											int* len, MPI_Datatype* data_type_ptr)
	{
		for ( int i=0; i<(*len); ++i ) {
			*inoutvec *= *invec;
			++inoutvec;  ++invec;
		}
	}
};

#include "MpiCommunicator.h"

// Registar for primitive types
template<typename T>
MpiDatatype MpiDatatypeRegistrar<T>::makeMpiDatatype(MpiCommunicator& comm)
{
	// Find the corresponding primitive
	std::type_index ti(typeid(T));
	const auto primitive_it = MpiDatatype::primitive_MPI_Datatype_map_.find(ti);
	if ( primitive_it == MpiDatatype::primitive_MPI_Datatype_map_.end() ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
					 << "  type \"" << ti.name() << "\" is not a recognized MPI primitive type \n"
					 << "  (NOTE: the type name printed above is implementation-dependent)\n";
		throw std::runtime_error( err_ss.str() );
	}

	return MpiDatatype( primitive_it->second );
}

template<typename T>
void MpiDatatypeRegistrar<T>::registerMpiOps(MpiOp::StandardOpsMap& map)
{
	using UserFunction = MpiOp::MpiUserFunction<T>;
	if ( std::is_arithmetic<T>::value ) {
		map.insert( std::make_pair( MpiOp::StandardOp::Sum, MpiOp(UserFunction(mpi_sum), true)));
	}
}


// Helper class for registering std::arrays
template<typename T, std::size_t dim>
MpiDatatype MpiDatatypeRegistrar<std::array<T,dim>>::makeMpiDatatype(MpiCommunicator& comm)
{
	// First, get the MpiDatatype corresponding to T
	const MpiDatatype& primitive_type = comm.map_mpi_datatype<T>();

	// Make a contiguous type
	MPI_Datatype mpi_array_type;
	auto count = static_cast<int>(dim);
#ifdef MPI_ENABLED
	MPI_Type_contiguous(count, primitive_type.get_MPI_Datatype(), &mpi_array_type);
#endif /* MPI_ENABLED */

	// Wrap it up in an MpiDatatype and return it
	return MpiDatatype( mpi_array_type );
}


template<typename T, std::size_t dim>
void MpiDatatypeRegistrar<std::array<T,dim>>::registerMpiOps(MpiOp::StandardOpsMap& map)
{
	using Array        = std::array<T,dim>;
	using UserFunction = MpiOp::MpiUserFunction<Array>;
	if ( std::is_arithmetic<T>::value ) {
		map.insert( std::make_pair( MpiOp::StandardOp::Sum, MpiOp(UserFunction(mpi_sum), true)));
	}
}


// Helper class for registering std::array<std::array>
template<typename T, std::size_t dim1, std::size_t dim2>
MpiDatatype MpiDatatypeRegistrar<std::array<std::array<T, dim2>, dim1>>::makeMpiDatatype(MpiCommunicator& comm)
{
	// First, get the MpiDatatype corresponding to T
	const MpiDatatype& primitive_type = comm.map_mpi_datatype<T>();

	// Make a contiguous type
	MPI_Datatype mpi_array_type;
	auto count = static_cast<int>(dim1*dim2);
#ifdef MPI_ENABLED
	MPI_Type_contiguous(count, primitive_type.get_MPI_Datatype(), &mpi_array_type);
#endif /* MPI_ENABLED */

	// Wrap it up in an MpiDatatype and return it
	return MpiDatatype( mpi_array_type );
}


template<typename T, std::size_t dim1, std::size_t dim2>
void MpiDatatypeRegistrar<std::array<std::array<T, dim2>, dim1>>
	::registerMpiOps(MpiOp::StandardOpsMap& map)
{
	using Matrix       = std::array<std::array<T, dim2>, dim1>;
	using UserFunction = MpiOp::MpiUserFunction<Matrix>;
	if ( std::is_arithmetic<T>::value ) {
		map.insert( std::make_pair( MpiOp::StandardOp::Sum, MpiOp(UserFunction(mpi_sum), true)));
	}
}


// Helper class for registering std::complex
template<typename T>
MpiDatatype MpiDatatypeRegistrar<std::complex<T>>::makeMpiDatatype(MpiCommunicator& comm)
{
	// First, get the MpiDatatype corresponding to T
	const MpiDatatype& primitive_type = comm.map_mpi_datatype<T>();

	// Make a contiguous type with 2 elements
	MPI_Datatype mpi_complex_type;
#ifdef MPI_ENABLED
	MPI_Type_contiguous(2, primitive_type.get_MPI_Datatype(), &mpi_complex_type);
#endif /* MPI_ENABLED */

	// Wrap it up in an MpiDatatype and return it
	return MpiDatatype( mpi_complex_type );
}

template<typename T>
void MpiDatatypeRegistrar<std::complex<T>>::registerMpiOps(MpiOp::StandardOpsMap& map)
{
	using Complex      = std::complex<T>;
	using UserFunction = MpiOp::MpiUserFunction<Complex>;
	if ( std::is_arithmetic<T>::value ) {
		map.insert( std::make_pair( MpiOp::StandardOp::Sum, MpiOp(UserFunction(mpi_sum), true)));
	}
}

#endif // #ifndef MPI_DATATYPE_REGISTRAR_H
