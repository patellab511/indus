/* ProbeVolume.h
 *
 * ABOUT: Abstract base class for implementing different probe volume geometries
 *  - Probe volumes which require a neighbor sphere radius are designed for usage
 *    with, e.g., order parameters which need to know the local neighbor list. They
 *    should still be usable for order parameters that don't care about the distinction:
 *    just set r_neighborSphere to some small value > 0 nm. (TODO: verify)
 */

#pragma once

#ifndef PROBE_VOLUME_H
#define PROBE_VOLUME_H

// Standard headers
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>    // unique_ptr, shared_ptr
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

// Project headers
#include "AtomGroup.h"
#include "BoundingBox.h"
#include "CommonTypes.h"
#include "GenericFactory.h"
#include "InputParser.h"
#include "MpiCommunicator.h"
#include "OpenMP.h"
#include "Region.h"
#include "SimulationBox.h"
#include "SimulationState.h"
#include "StringTools.h"
#include "SwitchingFunction_Indus_r.h"
#include "SwitchingFunction_Indus_x.h"
#include "SwitchingFunction_WithRegions_r.h"
#include "SwitchingFunction_WithRegions_x.h"
#include "Topology.h"

class ProbeVolume
{
 public:
	//----- Data Types -----//

  static const int DIM_  = CommonTypes::DIM_;  // Dimensionality of simulation 
	static const int X_DIM = CommonTypes::X_DIM;
	static const int Y_DIM = CommonTypes::Y_DIM;
	static const int Z_DIM = CommonTypes::Z_DIM;
  using Real3 = CommonTypes::Real3;
	using Box   = CommonTypes::Box;
  using Rvec  = CommonTypes::Rvec;
	using RvecArray  = CommonTypes::RvecArray;
	using RegionEnum = Region::RegionEnum;
  using Range      = CommonTypes::Range;

	// ProbeVolumeInputPack: packages input to derived class probe volumes
	struct ProbeVolumeInputPack {
		const ParameterPack& input_parameter_pack;
		const SimulationState& simulation_state;  // owned by the driver
		MpiCommunicator& mpi_communicator;  // owned by the driver
		const bool need_derivatives;  // whether derivatives are needed
	};

	ProbeVolume(ProbeVolumeInputPack& input_pack);


	//----- Virtual Functions -----//

	// Virtual destructor, so that inheriting classes may implement destruction
	virtual ~ProbeVolume() {}

	// Sets the geometry of the probe volume, using only member variables
	// (such as sigma_, alpha_c_, width_shell_1_, etc. from the base class, 
	//  and any other variables needed by a derived class)
	// - This is used to update the probe volume if/when parameters change
	virtual void setGeometry() = 0;

	// Update probe volume parameters for each time step using state variables 
	// from SimulationState (e.g. box matrix, atom positions, etc.)
	// - Examples of usage
	//   - Dynamic probe volumes based on atom positions must be updated at each time
	//     step as the atoms positions change
	//   - ProbeVolume's implementation often needs to know the box lengths in order to
	//     properly implement probe volumes which can cross the periodic boundaries
	virtual void update() {
		setGeometry();
		setBoundingBox();
		return;
	}

	// Checks whether the particle is in the probe volume
	// - Returns true if the particle is in the probe volume (no coarse-graining)
	// - If the particle is in the coarse-grained probe volume, then htilde_v > 0
	virtual void calculateIndicator(
		const Real3& x,                  // Particle position {x,y,z} [nm]
		// Output
		double&      h_v,                // Indicator without coarse-graining, h_v(x)
		double&      htilde_v,           // Indicator with coarse-graining, htilde_v(x)
		Real3&       dhtilde_v_dx,       // Derivatives of htilde_v w.r.t. position x
		RegionEnum&  region
	) const = 0;

	// Set the bounding box using the ProbeVolume geometry
	// - Default: the entire SimulationBox
	void setBoundingBox();

	// Returns a string with complete, formatted information about the probe volume's
	// location and geometry which can be included in output files. 
	// - Default: Each line begins with a "#" which signifies that it's a comment line.
	virtual std::string getInputSummary(
		const std::string& prepend_string = "# "
	) const;


	//----- Get functions -----//

	const BoundingBox& getBoundingBox() const { return bounding_box_; }

	// Access to parameters
	double get_sigma() const { return sigma_; }
	double get_alpha_c() const { return alpha_c_; }
	void   getShellWidths(double& width_shell_1, double& width_shell_2) const {
		width_shell_1 = width_shell_1_;
		width_shell_2 = width_shell_2_;
	}

	bool has_probe_volume_atoms() {
		return has_probe_volume_atoms_;
	}

	// Get summary of attributes shared by all probe volumes
	// ( e.g. shell widths, coarse-graining parameters, etc.)
	std::string getSharedAttributesSummary(const std::string& prepend_string) const;

	/*
	// TODO probe volume group
	int get_target_atom_group_index() const {
		return atom_group_indices_[target_group_];
	}
	*/


	//----- Misc. Set Functions -----//

	// Updates sigma_, alpha_c_, and any other member variables that depend upon them
	void setCoarseGrainingParameters(const double sigma, const double alpha_c) {
		sigma_   = sigma;
		alpha_c_ = alpha_c;
		setGeometry();
	}

	// Sets the size of the shells around the probe volume which are searched when 
	// looking for atoms NEAR (but not in) the probe volume
	void setShellWidths(const double width_shell_1, const double width_shell_2) {
		width_shell_1_ = width_shell_1;
		width_shell_2_ = width_shell_2;
		setGeometry();
	}

 protected:
	const SimulationState& simulation_state_;
	const SimulationBox&   simulation_box_;
	MpiCommunicator&       mpi_communicator_;
	bool need_derivatives_ = false;

	// TODO
	bool has_probe_volume_atoms_ = false;
	int pv_group_ = -1;

	// Safety factor for setting bounding boxes
	double bounding_box_tol_ = 1.0e-3;


	//----- Bounding Box -----//

	// Triclinic region that encompasses the probe volume (for DD)
	BoundingBox bounding_box_;
	bool use_default_bounding_box_ = false;
	std::function<BoundingBox()> bounding_box_fxn_;

	// Constructs the box that bounds the probe volume, including any coarse-graining
	// - Note: this does *not* include the shells
	// - Default (base class implementation): entire simulation box
	virtual BoundingBox constructBoundingBox() const;


	//----- General Parameters -----//

	double sigma_   = 0.01;  // Coarse-graining length scale [nm]
	double alpha_c_ = 0.02;  // Size of Gaussian buffer [nm]

	// Size of shells around probe volume
	// - This *includes* any coarse-graining that other objects may need
	double width_shell_1_ = 0.0, width_shell_2_ = 0.0;

	double getTotalWidthOfShells() const {
		return width_shell_1_ + width_shell_2_;
	}

	// Flags
	bool is_dynamic_ = false;


	//----- Constants -----//

	// Map x,y,z axes to indices (TODO make static constexpr?)
	const int X_ = 0, Y_ = 1, Z_ = 2;
	const std::map<int, std::string> axis_index_to_label_map_ = {
		{ X_, "x" },
		{ Y_, "y" },
		{ Z_, "z" }
	};

 private:
};


// Register ProbeVolumes using a GenericFactory
// - Use these typedefs for brevity
namespace ProbeVolumeRegistry {
template<typename P>
using Register = RegisterInFactory<
	ProbeVolume, P,                      // base class and derived class
	std::string,                         // generating key
  ProbeVolume::ProbeVolumeInputPack&   // input type(s)
>;

using Factory = GenericFactory< ProbeVolume, std::string, ProbeVolume::ProbeVolumeInputPack& >;
} // end namespace ProbeVolumeRegistry


#endif // PROBE_VOLUME_H	
