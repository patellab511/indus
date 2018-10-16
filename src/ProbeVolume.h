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
#include <iostream>
#include <memory>    // unique_ptr, shared_ptr
#include <sstream>
#include <stdexcept>
#include <string>

// Project headers
#include "AtomGroup.h"
#include "CommonTypes.h"
#include "InputParser.h"
#include "MpiCommunicator.h"
#include "SimulationBox.h"
#include "SimulationState.h"
#include "StringTools.h"
#include "SwitchingFunction_r.h"
#include "SwitchingFunction_x.h"
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
	using RvecArray = CommonTypes::RvecArray;

	// ProbeVolumeInputPack: packages input to derived class probe volumes
	struct ProbeVolumeInputPack {
		const ParameterPack& input_parameter_pack;
		const SimulationState& simulation_state;  // owned by the driver
		MpiCommunicator& mpi_communicator;  // owned by the driver
		const bool need_derivatives;  // whether derivatives are needed
	};

	ProbeVolume(ProbeVolumeInputPack& input_pack);

	// Computes INDUS Variables (see list below, in private data)
	void doIndusWithShells();


	//----- Virtual Functions -----//

	// Virtual destructor, so that inheriting classes may implement destruction
	virtual ~ProbeVolume() {}

	// Updates sigma_, alpha_c_, and any other member variables that depend upon them
	virtual void setCoarseGrainingParameters(
		const double sigma, 
		const double alpha_c
	) = 0;

	// Sets the size of the shells around the probe volume which are searched when 
	// looking for atoms NEAR (but not in) the probe volume
	virtual void setShellWidths(
		const double width_shell_1, 
		const double width_shell_2
	) = 0;

	// Update probe volume parameters for each time step using state variables 
	// from SimulationState (e.g. box matrix, atom positions, etc.)
	virtual void updateUsingSimulationState() = 0;

	// Checks whether the particle is in the probe volume
	// - Returns true if the particle is in the probe volume (no coarse-graining)
	// - If the particle is in the coarse-grained probe volume, then htilde_v > 0
	virtual bool isInProbeVolume(
		const Real3& x,                  // Particle position {x,y,z} [nm]
		// Output
		double&      h_v,                // Indicator without coarse-graining, h_v(x)
		double&      htilde_v,           // Indicator with coarse-graining, htilde_v(x)
		Real3&       dhtilde_v_dx,       // Derivatives of htilde_v w.r.t. position x
		bool&        is_in_shell_1,
		bool&        is_in_shell_2
	) const = 0;

	// Set the bounding box using the ProbeVolume geometry
	// - Default: the entire SimulationBox
	// - Must encompass the furthest extent of the shells (including coarse-
	//   graining width)
	virtual void setBoundingBox();

	// Returns a string with complete, formatted information about the probe volume's
	// location and geometry which can be included in output files. Each line begins
	// with a "#" which signifies that it's a comment line.
	virtual std::string getInputSummary() const;


	//----- Get functions -----//

	void getBoundingBox(Real3& box_offset, Box& box_matrix) const {
		box_offset = bounding_box_offset_;
		box_matrix = bounding_box_matrix_;
	}

	// Local atoms in vtilde
	const std::vector<int>& get_group_indices_of_local_indus_atoms() const {
		return local_atom_group_indices_;
	}
	const std::vector<int>& get_global_indices_of_local_indus_atoms() const {
		return local_atom_global_indices_;
	}
	const std::vector<double>& getIndicatorFunctions() const {
		return htilde_v_;
	}
	const std::vector<Real3>&  getIndicatorFunctionDerivatives() const {
		return derivatives_htilde_v_;
	}

	// Nonlocal atoms in vtilde
	const std::vector<int>& get_group_indices_of_nonlocal_indus_atoms() const {
		return nonlocal_atom_group_indices_;
	}
	const std::vector<int>& get_global_indices_of_nonlocal_indus_atoms() const {
		return nonlocal_atom_global_indices_;
	}

	// Shell 1 atoms
	const std::vector<int>& get_group_indices_of_nearby_shell_1_atoms() const {
		return nearby_shell_1_atom_group_indices_;
	}
	const std::vector<int>& get_global_indices_of_nearby_shell_1_atoms() const {
		return nearby_shell_1_atom_global_indices_;
	}

	// Shell 2 atoms
	const std::vector<int>& get_group_indices_of_nearby_shell_2_atoms() const {
		return nearby_shell_2_atom_group_indices_;
	}
	const std::vector<int>& get_global_indices_of_nearby_shell_2_atoms() const {
		return nearby_shell_2_atom_global_indices_;
	}

	// Global INDUS variables
	int get_n_v()         const { return n_v_;      }
	double get_ntilde_v() const { return ntilde_v_; }

	// Quick access to atom counts
	int get_num_local_atoms_in_vtilde()   const { return local_atom_group_indices_.size();    }
	int get_num_nearby_atoms_in_vtilde()  const { 
		return local_atom_group_indices_.size() + 
		       nonlocal_atom_group_indices_.size();
	}
	int get_num_nearby_atoms_in_shell_1() const { return nearby_shell_1_atom_group_indices_.size(); }
	int get_num_nearby_atoms_in_shell_2() const { return nearby_shell_2_atom_group_indices_.size(); }

	// Access to parameters
	double get_sigma() const { return sigma_; }
	double get_alpha_c() const { return alpha_c_; }
	void   getShellWidths(double& width_shell_1, double& width_shell_2) const {
		width_shell_1 = width_shell_1_;
		width_shell_2 = width_shell_2_;
	}

	// Get summary of attributes shared by all probe volumes
	// ( e.g. shell widths, coarse-graining parameters, etc.)
	std::string getSharedAttributesSummary() const;


 protected:
	// Objects owned by the driver
	const SimulationState& simulation_state_;
	const SimulationBox& simulation_box_;  // alias for SimulationState's SimulationBox

	// Objects shared with the driver
	MpiCommunicator& mpi_communicator_;


	//----- INDUS Variables -----//

	// Integer and coarse-grained number of particles in the probe volume
	// - doIndusWithShells() makes these values available on all ranks
	int    n_v_;
	double ntilde_v_;

	// Triclinic region that encompasses the probe volume
	Real3 bounding_box_offset_;  // position of lower-left corner
	Box   bounding_box_matrix_;  // shape of box

	// Vectors for local INDUS atoms
	// - These are atoms in vtilde which are in the local DD cell
	// - Derivatives of htilde_v are with respect to the corresponding atom position
	std::vector<int>    local_atom_global_indices_;
	std::vector<int>    local_atom_group_indices_;
	std::vector<double> htilde_v_;
	std::vector<Real3>  derivatives_htilde_v_;

	// Atoms in vtilde which are *not* local
	std::vector<int> nonlocal_atom_group_indices_;
	std::vector<int> nonlocal_atom_global_indices_;

	// Nearby atoms in the ProbeVolume's shell 1 and shell 2
	// - In the local DD cell or one of its shells
	std::vector<int> nearby_shell_1_atom_group_indices_;
	std::vector<int> nearby_shell_1_atom_global_indices_;

	// Nearby atoms in the ProbeVolume's shell 2
	// - In the local DD cell or one of its shells
	std::vector<int> nearby_shell_2_atom_group_indices_;
	std::vector<int> nearby_shell_2_atom_global_indices_; 


	//----- INDUS Parameters -----//

	double sigma_;    // Coarse-graining length scale [nm]
	double alpha_c_;  // Size of Gaussian buffer [nm]

	// Size of shells around probe volume
	// - This *includes* any coarse-graining that other objects may need
	double width_shell_1_, width_shell_2_;

	// Flags
	bool need_derivatives_;


	//----- Constants -----//

	// Map x,y,z axes to indices (TODO make static constexpr?)
	const int X_ = 0, Y_ = 1, Z_ = 2;
	const std::map<int, std::string> axis_index_to_label_map_ = {
		{ X_, "x" },
		{ Y_, "y" },
		{ Z_, "z" }
	};
};


#endif // PROBE_VOLUME_H	
