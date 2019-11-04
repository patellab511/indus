
#pragma once
#ifndef ORDER_PARAMETER_H
#define ORDER_PARAMETER_H

#include <array>
#include <vector>

#include "AtomGroup.h"
#include "CommonTypes.h"
#include "DomainDecomposition.h"
#include "GptlWrappers.h"
#include "MpiCommunicator.h"
#include "SimulationState.h"
#include "utils.h"

class OrderParameter
{
 public:
	//----- Basic Typedefs and Data Structures -----//

	static constexpr int DIM_  = CommonTypes::DIM_;   // Dimensionality of simulation	
  static constexpr int X_DIM = CommonTypes::X_DIM; // x-axis index
  static constexpr int Y_DIM = CommonTypes::Y_DIM; // y-axis index
  static constexpr int Z_DIM = CommonTypes::Z_DIM; // z-axis index

	// Arrays with 3 entries
	using Rvec  = CommonTypes::Rvec;
	using Real3 = CommonTypes::Real3;

	// Vectors (indefinite length)
	using RvecArray   = CommonTypes::RvecArray;
	using Vector      = CommonTypes::Vector;
	using VectorReal3 = CommonTypes::VectorReal3;

	// Matrices/tensors
	using Matrix = CommonTypes::Matrix; 
	using Box    = CommonTypes::Box; 

  using Range = CommonTypes::Range;

	struct InputPack {
		const SimulationState&     simulation_state;
		const DomainDecomposition& domain_decomposition;
		MpiCommunicator&           mpi_communicator;
		bool                       need_derivatives;
	};


	//----- Core Methods -----//

	OrderParameter(
		// std::string name,  // TODO e.g. 'q', 'ntilde', ...
		InputPack& input_pack
	);

	virtual ~OrderParameter();

	// Adds the AtomGroup with the given name to the list of dependencies
	// - Returns the index of the added AtomGroup in this OP's arrays
	//   (*not* the global index of the AtomGroup)
	int addAtomGroupDependency(const std::string& group_name);

	// Updates the OrderParameter at each step
	// - Called after all state variables and AtomGroups have been updated,
	//   but before calculate()
	virtual void update() {
		return;
	}

	virtual void calculate() = 0;

	// Prefer to put all global communication in this routine for performance reasons
	virtual void synchronize() {
		return;
	}

	// Wrap up (e.g. close output files managed directly by the OP)
	virtual void finalize() { return; }


	//----- Get Functions -----//

	double get_value() const {
		return value_; 
	}
	const std::vector<int>& get_atom_group_indices() const {
		return atom_group_indices_; 
	}
	const std::vector<std::vector<Real3>>& get_local_derivatives() const { 
		return local_derivatives_;
	}
	const std::vector<std::vector<int>>& get_local_group_indices() const {
		return local_group_indices_; 
	}
	const Matrix& get_local_sum_r_cross_derivatives() const {
		return local_sum_r_cross_derivatives_; 
	}
	const Matrix& get_sum_r_cross_derivatives() const {
		return sum_r_cross_derivatives_; 
	}

	// Returns a string with complete, formatted information about the order parameter,
	// which can be included in output files. 
	// - Default: Each line begins with a "#" which signifies that it's a comment line.
	virtual std::string getInputSummary(const std::string& prepend_string = "# ") const {
		return "";
	}

 protected:
	//std::string                   name_;  // TODO
	const SimulationState&        simulation_state_;
	const SimulationBox&          simulation_box_;
	const std::vector<AtomGroup>& atom_groups_;

	const DomainDecomposition& domain_decomposition_;

	MpiCommunicator& mpi_communicator_;

	bool need_derivatives_ = false;

	// The value of the order parameter
	double value_;

	// List of indices of AtomGroups used by this OrderParameter
	// - TODO better name?
	std::vector<int> atom_group_indices_;

	// Derivatives and indices for local atoms, organized by relevant AtomGroup
	std::vector<std::vector<Real3>> local_derivatives_;
	std::vector<std::vector<int>> local_group_indices_;

	Matrix local_sum_r_cross_derivatives_;
	Matrix sum_r_cross_derivatives_;

	void calculate_local_sum_r_cross_derivatives();
	void calculate_sum_r_cross_derivatives();

	mutable GPTL::Timer allreduce_timer_ = GPTL::Timer("OrderParameter_allreduce_r_cross_derivs");
};

#endif /* ORDER_PARAMETER_H */
