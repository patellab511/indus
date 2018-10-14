// Patch between INDUS and PLUMED

#ifndef INDUS_INTERFACE
#define INDUS_INTERFACE

// PLUMED headers
#include "ActionRegister.h"
#include "Colvar.h"

// INDUS headers
#include "CommonTypes.h"
#include "MpiCommunicator.h"
#include "Indus.h"
#include "SimulationBox.h"

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
#include <vector>

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
//     - In IndusInterface, ARG=ubias
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


/* We begin by declaring a class for your colvar.  This class inherits everything from the Colvar class.
   This ensures it has a label, a place to store its value, places to the store the values of the derivatives
   and that it can access the various atoms it will employ. */

class IndusInterface : public Colvar 
{
 public:
	static const int DIM_ = 3;
	static const int X_DIM = 0;
	static const int Y_DIM = 1;
	static const int Z_DIM = 2;
	using Real3  = CommonTypes::Real3;
	using Box    = CommonTypes::Box;
	using Matrix = CommonTypes::Matrix;

	// This routine is used to create the descriptions of all the keywords used by your CV
	static void registerKeywords( Keywords& keys ); 

	/* This is the constructor for your colvar.  It is this routine that will do all the reading.
     Hence it takes as input a line from the input file. */
	IndusInterface(const ActionOptions& ao);

	/* This is the routine that will be used to calculate the value of the colvar, whenever its calculation is required.
     This routine and the constructor above must be present - if either of them are not the code will not compile. */
	void calculate() override;


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
	// Input file with MDAnalysis++ parameters 
	std::string input_file_;

	//----- PLUMED variables -----//

	// MPI rank of this processor
	int my_rank_;

	// Used to request positions from the PLUMED core
	std::vector<PLMD::AtomNumber> indus_numbers_;

	// Reference to PLUMED's simulation box object
	const PLMD::Tensor& plumed_simulation_box_;

	//----- Interface with MDAnalysis++ -----//

	// Pointers to external objects
	std::unique_ptr<Indus> indus_ptr_;

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

	// Arrays based on no. atoms: indices, positions, and derivatives
	PLMD::Vector total_force_;  // TODO total force on OP


	//----- Debugging -----//

	void calculateMaxNormDeriv(
		const std::vector<Real3>& derivatives,
		// Output
		double& max_norm_deriv,
		int&    global_index_max_norm_deriv
	) const;
};

} // namespace colvar
} // namespace PLMD

#endif // INDUS_INTERFACE
