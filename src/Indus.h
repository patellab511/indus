/* Indus.h
 *
 * ABOUT: Driver for doing INDUS
 *   
 */

#pragma once
#ifndef INDUS_H
#define INDUS_H

// Check whether PLUMED is defined using one of its preprocessor variables
#ifndef INDUS_STANDALONE_MODE
#define INDUS_PLUMED_MODE
#endif

// Standard headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

// Project headers
#include "AtomGroup.h"
#include "Bias.h"
#include "CellGrid.h"
#include "CommonTypes.h"
#include "GenericFactory.h"
#include "InputParser.h"
#include "ProbeVolume.h"
#include "SimulationBox.h"
#include "SimulationState.h"
#include "StringTools.h"
#include "Topology.h"

// Non-PLUMED headers
#ifndef INDUS_PLUMED_MODE
#include "XdrFileTools.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment" /* ignore useless warnings from this header */
#include "xdrfile_xtc.h" /* read XTC files */
#pragma GCC diagnostic pop
#endif


class Indus 
{
 public:
	Indus(
		const std::string& input_file
	);


	//----- Basic Typedefs and Data Structures -----//

	static const int DIM_  = CommonTypes::DIM_;   // Dimensionality of simulation	
  static const int X_DIM = CommonTypes::X_DIM; // x-axis index
  static const int Y_DIM = CommonTypes::Y_DIM; // y-axis index
  static const int Z_DIM = CommonTypes::Z_DIM; // z-axis index

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

	// Arrays of zeros
	//using zero_3 = CommonTypes::zero_3;

	//---------------------------------//
	//----- Core Driver Functions -----//
	//---------------------------------//

#ifndef INDUS_PLUMED_MODE
	// Compute N_v and Ntilde_v for frames from an XTC file
	void do_indus_standalone();

#endif /* ifndef INDUS_PLUMED_MODE */

	// Updates simulation state variables for each frame (coordinates + box), and
	// re-does the domain decomposition if running with MPI
	// - Standalone MODE
	//    - RvecArray = std::SimpleVector<rvec> where rvec = float[3] from xdrfile
	//    - atom_positions: all atom positions
	// - PLUMED mode
	//    - RvecArray = std::vector<PLMD::Vector>
	//    - atom_positions: requested atoms [obtained via getPositions()]
	void updateSimulationState(
		const RvecArray& atom_positions,
		const Box&   box_matrix,
		const int    step,
		const double time
	);

	// Runs steinhardt_ptr->calculate() followed by this->calculateTotalBias()
	void calculate();

	// Calculates the total bias as well as its derivatives wrt. atom positions
	// (which are the negatives of the biasing forces)
	void calculateTotalBias();


	//---------------------------//
	//----- Misc. Functions -----//
	//---------------------------//

	// Compute P_v(N) from a series of numbers and bin it over the range [0,max_N],
	// where each bin "b" spans the range [b,b+1)
	void make_PvN(
		const std::vector<double>& n_values,
		const int max_N,
		// Output
		std::vector<double>& bins,
		std::vector<double>& pvn,
		double& avg_N,
		double& var_N
	) const;


	//-------------------------//
	//----- Get Functions -----//
	//-------------------------//

	// Global indices of atoms involved in calculations
	const std::vector<int>& get_target_atoms()           const { return target_atoms_; }
	const std::vector<int>& get_indus_global_indices() const { return indus_atom_global_indices_; }

	// Check object state
	bool areDerivativesReady() const { return derivatives_are_ready_; }

	// Access to ProbeVolume (INDUS) variables
	int get_n_v() const {
		return probe_volume_ptr_->get_n_v();
	}
	double get_ntilde_v() const {
		return probe_volume_ptr_->get_ntilde_v();
	}
	int get_num_local_atoms_in_vtilde() const {
		return probe_volume_ptr_->get_num_local_atoms_in_vtilde();
	}
	const std::vector<Real3>& get_derivatives_ntilde_v() const {
		 return derivatives_ntilde_v_;
	}
	const Matrix& get_sum_r_cross_dntilde_v_dr() const { 
		return sum_r_cross_dntilde_v_dr_; 
	}

	// Access biasing potential variables
	double get_u_bias_ntilde_v() const { return u_bias_ntilde_v_; }
	double get_u_bias_total()    const { return u_bias_total_; }
	const std::vector<Real3>& get_derivatives_u_bias_total() const {
		 return derivatives_u_bias_total_;
	}
	const Matrix& get_virial_u_bias_total() const { 
		return virial_u_bias_total_; 
	}
	const std::vector<Real3>& get_numerical_derivatives_u_bias_total() const {
		 return numerical_derivatives_u_bias_total_;
	}
	const Matrix& get_numerical_virial_u_bias_total() const { 
		return numerical_virial_u_bias_total_; 
	}

	bool printing_numerical_forces() const { return print_numerical_forces_; }

	// Returns a big string with a whole bunch of input variables.
	std::string getInputSummary() const;

	// Access to internal objects (for debugging purposes)
	const ProbeVolume* getProbeVolumePtr() const {
		return probe_volume_ptr_.get();
	}

	// Print the indicated derivatives to the indicated file stream
	void printForcesForFrame(
		std::ofstream& forces_ofs,
		const double time,
		const Matrix& box_matrix,
		const std::vector<Real3>& derivatives_u_bias_total,
		const Matrix& virial_u_bias_total
	) const;


	//-------------------------//
	//----- Set Functions -----//
	//-------------------------//

	void setTargetAtomsByGlobalIndices(
		const std::vector<int>& target_atoms
	);


 private:
	//----- Input -----//

	// Indus input files
	std::string input_file_;

	ParameterPack input_parameter_pack_;

	// Reference/data files (for e.g. topology, esp. atom types)
	std::string xtc_file_, gro_file_, top_file_;

	// Probe volume geometry
	std::string probe_volume_type_;

	// Object with handy routines for dealing with std strings
	StringTools string_tools_;


	//----- Simulation State -----//

	// Bundles together a number of useful simulation state variables
	// - Simulation topology
	// - Simulation box
	// - Atom groups
	SimulationState simulation_state_;


	//----- Domain Decomposition and MPI -----//

	// MPI communicator (or dummy, if MPI is not enabled)
	// - Default: MPI_COMM_WORLD
	MpiCommunicator mpi_communicator_;

	// Domain decomposition handler
	CellGrid domain_decomposition_;

	// MPI ranks
	const int my_rank_;
	const int master_rank_;


	//----- Indus Atoms -----//

	// Target atoms: the "primary" set of atoms upon which the biases act
	// - e.g. in simulations of water, these are typically the water oxygens
	std::string           target_atoms_string_; // "Target" line from input file
	std::vector<int>      target_atoms_;        //  global indices
	Topology::TargetFlags target_flags_;

	// Global indices of *all* atoms that are involved in computing the desired 
	// order parameter(s)
	// - Target atoms are listed first, followed by ProbeVolume atoms
	// - Created by concatenating the target atom and ProbeVolume atom vectors:
	//   if there is overlap between the two lists, it is preserved here
	std::vector<int> indus_atom_global_indices_;

	// Mapping from target/ProbeVolume atoms to Indus atoms
	std::vector<int> target_atom_indus_indices_;


	//----- Order Parameters -----//

	// Probe volume manages n_v and ntilde_v
	std::shared_ptr<ProbeVolume> probe_volume_ptr_;


	//----- Biases and Order Parameter Derivatives -----//

	// Controls whether derivatives are computed
	bool need_derivatives_;

	// Whether derivatives have been computed for the current set of coordinates
	bool derivatives_are_ready_;

	// Derivatives of ntilde_v_ wrt. the positions of all Indus atoms
	std::vector<Real3> derivatives_ntilde_v_;

	// Buffers for communicating derivatives
	template<typename T>
	struct AllgathervBuffer {
		std::vector<T>   buffer;
		std::vector<int> block_offsets;
		std::vector<int> block_sizes;
	};

	AllgathervBuffer<Real3> derivatives_ntilde_v_buffer_;
	AllgathervBuffer<int>   group_indices_for_ntilde_v_buffer_;


	// Matrices needed to compute the contribution of biases on ntilde_v
	// to the virial
	// - "cross" refers to the vecter outer product (vector direct product)
	// - This is the negative of PLUMED's "boxDerivatives"
	// - Ex. pseudo-LaTeX of the matrix for Ntilde_v:
	//     sum_r_cross_dntilde_v_dr = sum_{i=1}^{num_indus_atoms} r_i \cross dntilde_v/dr_i
	Box sum_r_cross_dntilde_v_dr_;

	// Input to Bias object
	std::string bias_ntilde_v_string_; 

	std::unique_ptr<Bias> bias_ntilde_v_ptr_;

	// Bias values (and derivatives)
	double      u_bias_ntilde_v_;           // Bias on ntilde_v
	double      u_bias_total_;              // Total bias
	VectorReal3 derivatives_u_bias_total_;  // derivatives wrt. Indus atoms' positions
	Matrix      virial_u_bias_total_;       // total virial

	// Numerical derivatives
	double delta_x_;  // how much to perturb each coordinate
	std::vector<Real3> numerical_derivatives_ntilde_v_;
	std::vector<Real3> numerical_derivatives_u_bias_total_;
	Matrix numerical_virial_u_bias_total_;
	

	//----- Output Handling -----//

	// Output control
	bool want_output_;
	bool print_forces_;
	bool print_numerical_forces_;
	std::string output_file_suffix_;  // Appended to output file names (before .out extension)

	// Production phase
	double t_min_, t_max_;  // (ps)


	//----- Derivatives and Biasing -----//

	// Gather local derivatives of ntilde_v from each MPI rank
	// on all processes
	void gatherDerivatives();

	// Takes derivatives obtained via a call to allgatherv() and computes
	// the total derivatives and the matrix for the virial
	void sumGatheredDerivatives(
		const AllgathervBuffer<Real3>& derivatives_x_buffer,
		const AllgathervBuffer<int>&   group_indices_for_x_buffer,
		const AtomGroup&               atom_group,
		const std::vector<int>&        op_indices_of_group_atoms,
		// Output
		std::vector<Real3>& derivatives_x,
		Box&                sum_r_cross_dx_dr
	);

	// Compute numerical derivatives and store the results in member variables
	// *** !!! This functions calls this->calculate() and overwrites member vars !!! ***
	void computeNumericalDerivatives();

	//----- Output Functions -----//

	// If "number" is finite, sends "number" to "ofs";
	// Else sends "nan" to "ofs"
	// - Standardizes what is printed when a number becomes non-finite
	//   across different systems (otherwise, regtests with NaNs fail)
	template<typename T>
	void print_nan_if_not_finite(std::ofstream& ofs, const T number) const {
		if ( std::isfinite(number) ) { ofs << number; }
		else                         { ofs << "nan";  }
	};
};

#endif /* INDUS_H */
