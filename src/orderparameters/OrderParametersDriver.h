/* OrderParametersDriver.h
 *
 * ABOUT: Driver object
 */

#pragma once
#ifndef ORDER_PARAMETERS_DRIVER_H
#define ORDER_PARAMETERS_DRIVER_H

#ifndef INDUS_STANDALONE_MODE
#define ORDER_PARAMETERS_DRIVER_PLUMED_MODE
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
#include "Bins.h"
#include "CellGrid_DomainDecomposition.h"
#include "DomainDecomposition.h"
#include "CommonTypes.h"
#include "GenericFactory.h"
#include "GenericFactoryInitialization.h"
#include "Indus.h"
#include "InputParser.h"
#include "OpenMP.h"
#include "OrderParameter.h"
#include "SimulationBox.h"
#include "SimulationState.h"
#include "StringTools.h"
#include "Topology.h"
#include "utils.h"

#include "GptlWrappers.h"

// Non-PLUMED headers
#ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE
#include "Statistics.h"
#include "XdrFileTools.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment" /* ignore useless warnings from this header */
#include "xdrfile/xdrfile_xtc.h" /* read XTC files */
#pragma GCC diagnostic pop
#endif


class OrderParametersDriver 
{
 public:
	OrderParametersDriver(
		const std::string& input_file,
		MpiCommunicator& mpi_communicator
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

#ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE
	// Main driver function
	void run();
#endif /* ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE */

	// Updates simulation state variables for each frame (coordinates + box), and
	// re-does the domain decomposition if running with MPI
	// - Standalone MODE
	//    - RvecArray = std::SimpleVector<rvec> where rvec = float[3] from xdrfile
	//    - atom_positions: all atom positions
	// - PLUMED mode
	//    - RvecArray = std::vector<PLMD::Vector>
	//    - atom_positions: requested atoms [obtained via getPositions()]
	void update(
		const RvecArray& atom_positions,
		const Box&       box_matrix,
		const int        step,
		const double     time
	);

	void updateAtomGroups(
		const RvecArray& global_atom_positions
	);

	void updateAtomGroupsPlumed(
		const RvecArray& op_atom_positions
	);

	void calculate();

	// Calculates the total bias as well as its derivatives wrt. atom positions
	// (which are the negatives of the biasing forces)
	void calculateTotalBias();

	// Call this once calculations are done in order to wrap up and
	// finish any last-minute tasks (e.g. write profiling info to file)
	void finalize();

	void printProfilingInfo();


	//-------------------------//
	//----- Get Functions -----//
	//-------------------------//

	// Global indices of atoms involved in calculations
	const std::vector<int>& get_op_atom_global_indices() const { return op_atom_global_indices_; }

	// Check object state
	bool areDerivativesReady() const { return derivatives_are_ready_; }

	double get_value() const { 
		return this->get_ntilde_v();
	}

	// Access to INDUS variables
	int get_n_v() const {
		return indus_ptr_->get_n_v();
	}
	double get_ntilde_v() const {
		return indus_ptr_->get_ntilde_v();
	}
	int get_num_local_atoms_in_vtilde() const {
		return indus_ptr_->get_num_local_atoms_in_vtilde();
	}
	const std::vector<Real3>& get_derivatives_ntilde_v() const {
		 return derivatives_ops_[index_indus_];
	}
	const Matrix& get_sum_r_cross_dntilde_v_dr() const { 
		return sum_r_cross_op_derivatives_[index_indus_]; 
	}

	// Access biasing potential variables
	double get_u_bias_ntilde_v() const { return u_bias_ops_[index_indus_]; }
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
	const std::vector<AtomGroup>& getAtomGroups() const {
		return simulation_state_.get_atom_groups();
	}
	const ProbeVolume* getProbeVolumePtr() const {
		return probe_volume_ptr_.get();
	}

	MpiCommunicator& access_mpi_communicator() {
		return mpi_communicator_;
	}

	const SimulationBox& get_simulation_box() const {
		return simulation_state_.get_simulation_box();
	}


	//---------------------------//
	//----- Output Handling -----//
	//---------------------------//

	// Print the indicated derivatives to the indicated file stream
	void printForcesForFrame(
		std::ofstream& forces_ofs,
		const double time,
		const Matrix& box_matrix,
		const std::vector<Real3>& derivatives_u_bias_total,
		const Matrix& virial_u_bias_total
	) const;

	void set_share_all_derivatives(const bool share) {
		share_all_op_derivatives_ = share;
	}
	bool share_all_derivatives() const {
		return share_all_op_derivatives_;
	}

	// For debugging
	friend void calculateForFrame(
		const std::string& op_input_file,
		const std::string& coord_file,
		MpiCommunicator& mpi_communicator
	);


	//---------------------------//
	//----- Setup Functions -----//
	//---------------------------//

	// Initialize all AtomGroups
	void initializeAtomGroups();


 private:
	//----- Input -----//

	// OrderParametersDriver input files
	std::string input_file_;

	ParameterPack input_parameter_pack_;

	// Reference/data files (for e.g. topology, esp. atom types)
	std::string xtc_file_, gro_file_, top_file_;

	// Used to initialize probe volume
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

	MpiCommunicator& mpi_communicator_;  // owned by the driver's owner

	DomainDecomposition domain_decomposition_;

	CellGrid_DomainDecomposition& dd_cell_grid_;

	// MPI ranks
	const int my_rank_;
	const int master_rank_;


	//----- OrderParametersDriver Atoms -----//

	// Global indices of *all* atoms that are involved in computing the desired 
	// order parameter(s)
	// - Target atoms are listed first, followed by ProbeVolume atoms
	// - Created by concatenating the target atom and ProbeVolume atom vectors:
	//   if there is overlap between the two lists, it is preserved here
	std::vector<int> op_atom_global_indices_;

	// Mapping from AtomGroup indices to OrderParametersDriver indices
	// - Size: [ num_groups x num_atoms_per_group ]
	std::vector<std::vector<int>> atom_group_to_op_indices_;


	//----- Order Parameters -----//

	// Probe volume
	bool only_do_indus_;  // If true, only compute n and ntilde
	std::shared_ptr<ProbeVolume> probe_volume_ptr_;

	std::unique_ptr<Indus> indus_ptr_;

	// Organize OPs for conveniene (e.g. ease of handling derivatives)
	// - TODO: Make these the owning ptrs?
	std::vector<OrderParameter*> op_ptrs_;
	int index_indus_ = -1;


	//----- Biases and Order Parameter Derivatives -----//

	// Controls whether derivatives are computed
	bool need_derivatives_;

	// Whether derivatives have been computed for the current set of coordinates
	bool derivatives_are_ready_;

	// Buffers for communicating derivatives
	template<typename T>
	struct AllgathervBuffer {
		std::vector<T>   buffer;
		std::vector<int> block_offsets;
		std::vector<int> block_sizes;
	};
	mutable AllgathervBuffer<Real3> derivatives_buffer_;
	mutable AllgathervBuffer<int>   group_indices_buffer_;

	// Input to Bias objects
	std::string bias_ntilde_v_string_; 
	std::unique_ptr<Bias> bias_ntilde_v_ptr_;

	//
	std::vector<Bias*>  bias_ptrs_;  // Non-owning pointers to above biases
	std::vector<std::vector<Real3>> local_derivatives_ops_;  // derivatives of each OP, across all relevant AtomGroups
	std::vector<std::vector<int>>   local_op_indices_for_derivatives_;  // ^corresponding "OP atom" indices

	// Matrices needed to compute the contribution of biases to the virial
	// - "cross" refers to the vecter outer product (vector direct product)
	// - This is the negative of PLUMED's "boxDerivatives"
	// - Ex. pseudo-LaTeX the matrix for 'x':
	//     sum_r_cross_dx_dr = sum_{i=1}^{num_op_atoms} r_i \cross dx/dr_i
	std::vector<Matrix> sum_r_cross_op_derivatives_;

	// Biases on individual OPs
	// - TODO: Move to separate class?
	std::vector<double>             u_bias_ops_;
	std::vector<std::vector<Real3>> local_derivatives_u_bias_ops_;
	std::vector<Matrix>             virial_u_bias_ops_;

	// Total bias, across all OPs (and its derivatives)
	double             u_bias_total_;                    // Total bias
	std::vector<Real3> local_derivatives_u_bias_total_;  // derivatives computed locally
	std::vector<int>   local_op_indices_u_bias_total_;   // ^corresponding "OP atom" indices
	std::vector<Real3> derivatives_u_bias_total_;  // derivatives wrt. *all* OrderParametersDriver atoms' positions
	Matrix             virial_u_bias_total_;       // total virial

	// Derivatives of each OP wrt. the positions of *all* OrderParameters atoms
	// - Only shared across MPI ranks if the flag is set; otherwise, only derivatives of
	//   u_bias_total are shared across all ranks
	bool share_all_op_derivatives_ = false;
	std::vector<std::vector<Real3>> derivatives_ops_;

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

	// Gather local derivatives of all OPs
	// - For each OP, collect derivatives etc. across all relevant AtomGroups
	void collectLocalOrderParameterDerivatives();

	// Share all derivatives of the OPs across all ranks (in addition to derivatives of u_bias_total)
	// - This is only useful if you're going to do something *else* with the derivatives,
	//   like pass them to another PLUMED function
	// - Only use this feature if you know you need it: it's expensive
	void shareAllOrderParameterDerivatives();

	// Compute numerical derivatives and store the results in member variables
	// *** !!! This functions calls this->calculate() and overwrites member vars !!! ***
	void computeNumericalDerivatives();


	//----- Output Functions -----//

	// TODO move to a 'Distribution' object for handling free energy distributions, F(x)

	// When computed from a histogram, F(x) is finite only when num_samples > 0
	template<typename T, typename I>
	static bool is_f_x_finite(const T f_x, const I num_samples) {
		return ( num_samples > 0 ? true : false );
	}

	// Get the minimum of F(x), given that not all bins have samples
	template<typename T, typename I>
	static T get_min_f_x(const std::vector<T>& f_x, const std::vector<I>& num_samples) {
		T min      = std::numeric_limits<T>::max();
		I num_bins = num_samples.size();
		for ( int b=0; b<num_bins; ++b ) {
			if ( is_f_x_finite(f_x[b], num_samples[b]) ) {
				min = std::min( min, f_x[b] );
			}
		}
		return min;
	}

	template<typename T, typename I>
	static std::vector<T> shift_f_x_to_zero(const std::vector<T>& f_x, const std::vector<I>& num_samples) {
		T min_f_x = get_min_f_x(f_x, num_samples);

		// Subtract minimum
		std::vector<T> f_x_shifted(f_x);
		std::for_each( f_x_shifted.begin(), f_x_shifted.end(), [=](T& f) { f -= min_f_x; } );
		return f_x_shifted;
	}

	// If 'f_x' is finite, print it
	// Else print 'nan'
	// - Standardizes what is printed when a number becomes non-finite
	//   across different systems (otherwise, regtests with NaNs fail)
	// TODO Make an object and overload operator<< to make usage less clunky
	template<typename T, typename I>
	static void print_free_energy(std::ofstream& ofs, const T f_x, const I num_samples) {
		if ( is_f_x_finite(f_x, num_samples) ) { ofs << f_x;   }
		else                                   { ofs << "nan"; }    
	}


	//----- GPTL -----//

	using Timer = GPTL::Timer;

	// Initialization
	mutable Timer setup_timer_ = Timer("OP_Driver_setup");

	mutable Timer read_xtc_timer_       = Timer("OP_Driver_read_xtc");
	mutable Timer calculate_timer_      = Timer("OP_Driver_calculate");
	mutable Timer op_calculate_timer_   = Timer("OP_Driver_OP_calculate()");
	mutable Timer op_synchronize_timer_ = Timer("OP_Driver_OP_synchronize()");

	// Handling derivatives and biases
	mutable Timer bias_timer_        = Timer("OP_Driver_calculateTotalBias");
	mutable Timer bias_gather_timer_ = Timer("calculateTotalBias_allgatherv");
	mutable Timer bias_sum_timer_    = Timer("calculateTotalBias_sum");

	// Handling derivatives
	mutable Timer collect_local_derivs_timer_ = 
		Timer("OP_Driver_collectLocalOrderParameterDerivatives");
	mutable Timer share_all_derivs_timer_ = 
		Timer("OP_Driver_shareAllOrderParameterDerivatives");

	// Updating
	mutable Timer update_timer_             = Timer("OP_Driver_update");
	mutable Timer update_atom_groups_timer_ = Timer("OP_Driver_update_atom_groups");

	// Output
	mutable Timer print_forces_timer_ = Timer("OP_Driver_printForcesForFrame");
};

#endif /* ORDER_PARAMETERS_DRIVER_H */
