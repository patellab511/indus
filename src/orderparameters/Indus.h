/* Indus.h
 *
 * ABOUT: Abstract base class for implementing different probe volume geometries
 *  - Probe volumes which require a neighbor sphere radius are designed for usage
 *    with, e.g., order parameters which need to know the local neighbor list. They
 *    should still be usable for order parameters that don't care about the distinction:
 *    just set r_neighborSphere to some small value > 0 nm. (TODO: verify)
 */

#pragma once

#ifndef INDUS_H
#define INDUS_H

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
#include "CommonTypes.h"
#include "GenericFactory.h"
#include "InputParser.h"
#include "MpiCommunicator.h"
#include "OpenMP.h"
#include "OrderParameter.h"
#include "ProbeVolume.h"
#include "Region.h"
#include "SimulationBox.h"
#include "SimulationState.h"
#include "StringTools.h"
#include "SwitchingFunction_Indus_r.h"
#include "SwitchingFunction_Indus_x.h"
#include "SwitchingFunction_WithRegions_r.h"
#include "SwitchingFunction_WithRegions_x.h"
#include "Topology.h"

class Indus : public OrderParameter
{
 public:
	//----- Data Types -----//

	using RegionEnum = Region::RegionEnum;
	using ProbeVolumeInputPack = ProbeVolume::ProbeVolumeInputPack;

	// TODO use own InputPack object? what about a generic OrderParameterInputPack object?
	Indus(
		std::shared_ptr<ProbeVolume>& probe_volume_ptr,
		ProbeVolumeInputPack& input_pack,
		OrderParameter::InputPack& op_input_pack
	);

	// Virtual destructor, so that inheriting classes may implement destruction
	virtual ~Indus() {}

	// Computes INDUS Variables, such as N_v and Ntilde_v
	// - See list below, in private data
	virtual void calculate() override;

	virtual void synchronize() override;


	//----- Get functions -----//

	// Local atoms in vtilde
	const std::vector<int>& get_group_indices_of_local_indus_atoms() const {
		return local_group_indices_[target_group_];
	}
	const std::vector<double>& getIndicatorFunctions() const {
		return htilde_v_;
	}
	const std::vector<Real3>&  getIndicatorFunctionDerivatives() const {
		return local_derivatives_[target_group_];
	}

	// Nonlocal atoms in vtilde
	const std::vector<int>& get_group_indices_of_nonlocal_indus_atoms() const {
		return nonlocal_atom_group_indices_;
	}

	// Shell 1 atoms
	const std::vector<int>& get_group_indices_of_nearby_shell_1_atoms() const {
		return nearby_shell_1_atom_group_indices_;
	}

	// Shell 2 atoms
	const std::vector<int>& get_group_indices_of_nearby_shell_2_atoms() const {
		return nearby_shell_2_atom_group_indices_;
	}

	// Core INDUS output
	int get_n_v()         const { return n_v_;      }
	double get_ntilde_v() const { return this->get_value(); }

	// Quick access to atom counts
	int get_num_local_atoms_in_vtilde()   const { return local_group_indices_[target_group_].size();    }
	int get_num_nearby_atoms_in_vtilde()  const { 
		return local_group_indices_[target_group_].size() + 
		       nonlocal_atom_group_indices_.size();
	}
	int get_num_nearby_atoms_in_shell_1() const { return nearby_shell_1_atom_group_indices_.size(); }
	int get_num_nearby_atoms_in_shell_2() const { return nearby_shell_2_atom_group_indices_.size(); }

	// Get summary of attributes shared by all probe volumes
	// ( e.g. shell widths, coarse-graining parameters, etc.)
	std::string getSharedAttributesSummary() const;

	int get_target_atom_group_index() const {
		return atom_group_indices_[target_group_];
	}


 protected:
	int target_group_ = -1;

	//----- INDUS Variables -----//

	// Integer number of particles in the probe volume
	// - 'value_'(base class variable) = ntilde_v
	int n_v_;

	// Coarse-grained indicator functions for local INDUS atoms
	// - These are atoms in vtilde which are in the local DD cell
	std::vector<double> htilde_v_;

	// Atoms in vtilde which are *not* local
	std::vector<int> nonlocal_atom_group_indices_;

	// Indices of nearby atoms in the ProbeVolume's shell 1 and shell 2
	// - Nearby --> in the local DD cell or one of its shells
	std::vector<int> nearby_shell_1_atom_group_indices_;
	std::vector<int> nearby_shell_2_atom_group_indices_;


	//----- Constants -----//

	const std::map<int, std::string> axis_index_to_label_map_ = {
		{ X_DIM, "x" },
		{ Y_DIM, "y" },
		{ Z_DIM, "z" }
	};

 private:
	std::shared_ptr<ProbeVolume> probe_volume_ptr_;

	// Timers
	mutable GPTL::Timer calculate_timer_              = GPTL::Timer("Indus_calculate");
	mutable GPTL::Timer calculate_omp_timer_          = GPTL::Timer("Indus_omp");
	mutable GPTL::Timer calculate_omp_critical_timer_ = GPTL::Timer("Indus_omp_critical");
	mutable GPTL::Timer allreduce_timer_              = GPTL::Timer("Indus_allreduce");
};

#endif // ifndef INDUS_H	
