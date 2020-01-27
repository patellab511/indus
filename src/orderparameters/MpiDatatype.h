// MpiDatatype: Wrapper around MPI_Datatype

#pragma once
#ifndef MPI_DATATYPE_H
#define MPI_DATATYPE_H

#include "MpiEnvironment.h"

// Standard headers
#include <array>
#include <complex>
#include <cstdlib>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include "MpiOp.h"
#include "utils.h"


#ifndef MPI_ENABLED
// Dummy classes for when MPI is not available
class MPI_Datatype {};
#endif /* MPI_ENABLED */


class MpiDatatype 
{
	friend class MpiCommunicator;

	template<typename T>
	friend class MpiDatatypeRegistrar;

 public:
	MpiDatatype();

	// Copy constructor
	// - Note that the copied MPI_Datatype has the same 'commitment state' as the
	//   MPI_Datatype from which it was copied
	//   - i.e. if the original MPI_Datatype was committed (see MPI_Type_commit), then
	//     the resulting MPI_Datatype is also committed
	MpiDatatype(const MpiDatatype& mpi_datatype);

	// Copy assignment operator
	MpiDatatype& operator=(const MpiDatatype& mpi_datatype);

	// Destructor must clean up underlying MPI_Datatype
	// - Note that MPI_Datatypes copied from existing MPI_Datatypes are unaffected
	//   when the original is freed
	~MpiDatatype() { free(); }

	// Create and register a contiguous type derived from an existing MpiDatatype
	MpiDatatype(const MpiDatatype& type, const std::size_t dim);

	// Copy from existing MPI datatype directly
	MpiDatatype(const MPI_Datatype& mpi_datatype);

	// Access underlying MPI_Datatype (for calls to MPI routines)
	const MPI_Datatype& get_MPI_Datatype() const { return mpi_datatype_; }

	// Map from C++ type to primitives such as MPI_DOUBLE, MPI_INT, etc.
	template<typename T>
	static const MPI_Datatype& map_primitive_MPI_Datatype();

#ifdef MPI_ENABLED
	// FIXME: need to figure out how to handle comparing types which are equivalent
	//        in order to use MpiDatatype (instead of MPI_Datatype) in MpiData
	// Checks whether two MpiDatatypes correspond to the same MPI_Datatype
	// - Does *not* check whether the type has been committed!
	bool operator==(const MpiDatatype& rhs) const {
		return (this->mpi_datatype_ == rhs.mpi_datatype_);
	}
	bool operator!=(const MpiDatatype& rhs) const {
		return (not (*this == rhs));
	}
#endif // #ifdef MPI_ENABLED

 private:
	// Underlying MPI struct
	MPI_Datatype mpi_datatype_;

	// Whether the MPI_Datatype was committed
	bool is_allocated_ = false;

	// Commit/decommit the type with the MPI core
	void commit();
	void free();

	// Map from C++ types to predefined MPI_Datatypes like int, char, ...
	static const std::unordered_map<std::type_index, MPI_Datatype> primitive_MPI_Datatype_map_;
};


//--------------------------------//
//----- Template Definitions -----//
//--------------------------------//


// Map from C++ type to primitives such as MPI_DOUBLE, MPI_INT, etc.
template<typename T>
const MPI_Datatype& MpiDatatype::map_primitive_MPI_Datatype() 
{
	std::type_index ti(typeid(T));
	const auto it = primitive_MPI_Datatype_map_.find(ti);
	if ( it != primitive_MPI_Datatype_map_.end() ) {
		return it->second;
	}
	else {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
					 << "  type \"" << ti.name() << "\" is not a registered primitive type \n"
					 << "  (NOTE: the type name printed above is implementation-dependent)\n";
		throw std::runtime_error( err_ss.str() );
	}
}


#endif /* MPI_DATATYPE_H */
