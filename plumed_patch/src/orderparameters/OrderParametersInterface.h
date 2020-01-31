// PLUMED headers
#ifndef ORDER_PARAMETERS_INTERFACE
#define ORDER_PARAMETERS_INTERFACE

// Standard headers
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// PLUMED headers
#include "../core/ActionRegister.h"
#include "../core/Colvar.h"

// OrderParameters headers
#include "CommonTypes.h"
#include "GptlWrappers.h"
#include "MpiCommunicator.h"
#include "OpenMP.h"
#include "OrderParametersDriver.h"
#include "SimulationBox.h"
#include "System.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR NAME
/*
TODO documentation
*/
//+ENDPLUMEDOC

// *** Notes on box derivatives and biasing ***
// - The "box derivatives" are proportional to a colvar's contribution to othe virial
//   - In general, the box derivatives have a bit of an odd definition. Let ARG=x be the colvar. 
//     Then:
//            boxDerivatives(ARG) = -1 * sum_{i=1}^{N} r_i (x) d(ARG)/dr_i
//
//   - Let "Xi_bias" be the contribution to the total virial from of biasing ARG. Then:
//            Xi_bias = 1/2 * [-dU_bias/d(ARG)] * boxDerivatives(ARG)
//     - In OrderParamtersInterface, ARG=ubias
//     
// - The desired biasing forces are applied to the atoms by passing ubias to a RESTRAINT as:
//                RESTRAINT ARG=op.ubias SLOPE=1.0 AT=0.0 KAPPA=0.0
//  
//  - Choosing ARG=ubias means that dU_bias/d(ARG) = dU_bias/dU_bias = 1
//    - PLUMED will subsequently multiply both the "atomsDerivatives" and "boxDerivatives"
//      by the "force" on the colvar (ubias): f_ubias = -dU_bias/d(ARG) = -1
//    - Hence the "atomsDerivatives" are turned into biasing forces and the "extra" sign in
//      the "boxDerivatives" cancelled out, turning it into 2*Xi_bias
//  - NOTE: PLUMED requires the prefactor of 1/2 be left out in setBoxDerivatives, and does not 
//    apply it in RESTRAINT. Instead, it is included later later when PLUMED interfaces
//    with GROMACS
//    - The prefactor is added in "at the last second" because it's a hotfix: the original code
//      forgot to include it, and code was written around that error, so the most 
//      straightforward way to fix everything was to resolve the error on the back end


class OrderParametersInterface : public Colvar 
{
 public:
	static constexpr int DIM_ = 3;
	static constexpr int X_DIM = 0;
	static constexpr int Y_DIM = 1;
	static constexpr int Z_DIM = 2;
	using Real3  = CommonTypes::Real3;
	using Box    = CommonTypes::Box;
	using Matrix = CommonTypes::Matrix;

	// Creates the descriptions of all the keywords used by the CV
	static void registerKeywords( Keywords& keys ); 

	OrderParametersInterface(const ActionOptions& ao);

	virtual void calculate() override;

	// This function is run at the very end, as PLUMED wraps up
	// FIXME: It should be called at this stage (according to PLUMED documentation),
	// but in practice, it doesn't seem to ever get called
	virtual void runFinalJobs() override;

	//----- Debugging Routines -----//

	enum class FrameFileType {
		gro,  // standard gro file, but with derivatives of q after the velocities
		xyz   // custom xyz file
	};

	void printFrameWithDerivatives(
		const std::string& file_name, 
		const std::string& message,
		const FrameFileType& file_type
	) const;


 private:
	// Input file with OrderParameters parameters 
	std::string input_file_;

	// Underlying MPI_Comm is copied from from PLUMED's Communicator (this->comm)
	MpiCommunicator mpi_communicator_;

	//----- PLUMED variables -----//

	// MPI rank of this processor
	int my_rank_;

	// Used to request positions from the PLUMED core
	std::vector<PLMD::AtomNumber> op_atom_numbers_;

	// Reference to PLUMED's simulation box object
	const PLMD::Tensor& plumed_simulation_box_;

	//----- Interface with OrderParameters -----//

	// Pointers to external objects
	std::unique_ptr<OrderParametersDriver> driver_ptr_;

	// Copy of PLUMED's simulation box
	CommonTypes::Box box_matrix_;

	// INDUS variables
	int    n_v_;
	double ntilde_v_;
	PLMD::Tensor box_derivatives_ntilde_v_;  // proportional to contribution of bias to virial

	// Biasing potential
	double u_bias_;
	double u_bias_ntilde_v_;
	PLMD::Tensor box_derivatives_u_bias_;  // proportional to contribution of bias to virial

	// Global index of the atom with the largest norm(derivative) and the norm itself
	// - Useful for identifying abnormally large forces
	int    global_index_max_norm_deriv_;
	double max_norm_deriv_;
	double max_norm_deriv_tol_; // tolerance
	bool   check_max_norm_deriv_;

	// Arrays based on no. atoms: indices, positions, and derivatives
	PLMD::Vector total_force_;  // TODO total force on OP


	//----- Debugging -----//

	void calculateMaxNormDeriv(
		const std::vector<Real3>& derivatives,
		// Output
		double& max_norm_deriv,
		int&    global_index_max_norm_deriv
	) const;


	//----- GPTL -----//

	// TODO More timers
	GPTL::Timer update_timer_           = GPTL::Timer("OP_Interface_update");
	GPTL::Timer calculate_timer_        = GPTL::Timer("OP_Interface_calculate");
	GPTL::Timer component_output_timer_ = GPTL::Timer("OP_Interface_component_output");


	//----- Pressure due to Bias -----//

	static const std::map<int, std::string> axis_index_to_label_map_;
	PLMD::Tensor p_bias_;
};

} // namespace colvar
} // namespace PLMD

#endif // ORDER_PARAMETERS_INTERFACE
