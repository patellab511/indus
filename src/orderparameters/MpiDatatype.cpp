#include "MpiDatatype.h"


MpiDatatype::MpiDatatype():
	is_allocated_(false) 
{}


// Explicit copy constructor implementation is required in order to duplicate 
// and commit the underlying MPI_Datatype
MpiDatatype::MpiDatatype(const MpiDatatype& mpi_datatype)
{
#ifdef MPI_ENABLED
	if ( mpi_datatype.is_allocated_ ) {
		// Duplicate the type (using const_cast is safe since the existing type will not be modified)
		MPI_Type_dup(const_cast<MPI_Datatype&>(mpi_datatype.mpi_datatype_), &(this->mpi_datatype_));
		MPI_Type_commit(&(this->mpi_datatype_));  // TODO need to commit? is this done automatically?
		this->is_allocated_ = true;
	}
	else {
		this->is_allocated_ = false;
	}
#else
	this->mpi_datatype_ = mpi_datatype.mpi_datatype_;
	this->is_allocated_ = mpi_datatype.is_allocated_;
#endif /* MPI_ENABLED */
}


// Copy assignment operator
MpiDatatype& MpiDatatype::operator=(const MpiDatatype& mpi_datatype)
{
	if ( &mpi_datatype != this ) {  // check for self-assignment
		if ( mpi_datatype.is_allocated_ ) {
			// Free stored type (if applicable)
			free();

			// Duplicate the type (using const_cast is safe since the existing type will not be modified)
#ifdef MPI_ENABLED
			MPI_Type_dup(const_cast<MPI_Datatype&>(mpi_datatype.mpi_datatype_), &(this->mpi_datatype_));
#endif /* MPI_ENABLED */
			commit();

			// TODO need to commit? is this done automatically? if so, need to free later.
			// - What does "same committed state" ultimately imply?
		}
		else {
			this->is_allocated_ = false;
		}
	}
	return *this;
}


// Create and register a contiguous type derived from an existing MpiDatatype
MpiDatatype::MpiDatatype(const MpiDatatype& type, const std::size_t dim):
	is_allocated_(false)
{
#ifdef MPI_ENABLED
	auto count = static_cast<int>(dim);
	MPI_Type_contiguous(count, type.get_MPI_Datatype(), &(this->mpi_datatype_));
	commit();
#endif /* MPI_ENABLED */
}


// Copy from existing MPI datatype directly
MpiDatatype::MpiDatatype(const MPI_Datatype& mpi_datatype):
	is_allocated_(false)
{
#ifdef MPI_ENABLED
	MPI_Type_dup(const_cast<MPI_Datatype&>(mpi_datatype), &(this->mpi_datatype_));
	commit();
#else
	(void) mpi_datatype;
#endif /* MPI_ENABLED */
}


// Commit the type to MPI
void MpiDatatype::commit() {
	if ( not is_allocated_ ) {
#ifdef MPI_ENABLED
		MPI_Type_commit(&mpi_datatype_);
		is_allocated_ = true;
#endif /* MPI_ENABLED */
	}
}


// Free (de-commit) the datatype
void MpiDatatype::free() {
	if ( is_allocated_ ) {
#ifdef MPI_ENABLED
		MPI_Type_free(&mpi_datatype_);
		is_allocated_ = false;
#endif /* MPI_ENABLED */
	}
}


#ifdef MPI_ENABLED
// Mapping from C++ primitives to MPI_Datatype primitives
// TODO add remaining primitives
const std::unordered_map<std::type_index, MPI_Datatype> 
		MpiDatatype::primitive_MPI_Datatype_map_ = {
	// Floating-point types
	{ std::type_index(typeid(float)),               MPI_FLOAT          },
	{ std::type_index(typeid(double)),              MPI_DOUBLE         },
	{ std::type_index(typeid(long double)),         MPI_LONG_DOUBLE    },
	// Integers
	{ std::type_index(typeid(int)),                 MPI_INT            },
	{ std::type_index(typeid(long int)),            MPI_LONG           },
	{ std::type_index(typeid(short int)),           MPI_SHORT          },
	{ std::type_index(typeid(unsigned int)),        MPI_UNSIGNED       },
	{ std::type_index(typeid(unsigned long int)),   MPI_UNSIGNED_LONG  },
	{ std::type_index(typeid(unsigned short int)),  MPI_UNSIGNED_SHORT },
	// Characters
	{ std::type_index(typeid(char)),                MPI_CHAR           },
	{ std::type_index(typeid(unsigned char)),       MPI_UNSIGNED_CHAR  }
};
#else
// TODO populate dummy version of map? (maybe by declaring MPI_Datatype dummy classes with names MPI_FLOAT, etc.)
const std::unordered_map<std::type_index, MPI_Datatype> 
		MpiDatatype::primitive_MPI_Datatype_map_;
#endif /* MPI_ENABLED */
