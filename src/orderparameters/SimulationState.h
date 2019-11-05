// SimulationState
//  - Unified place for simulation state variables to live
//    (e.g. atom_positions, box, etc.)
//
//  - Simplifies sending all this information to objects that need it 
//    (probe volumes, Steinhardt OPs, ...)
//     - No need to explicitly pass around the latest set of atom_positions
//     - Get functions control access: an object that receives a const SimulationState&
//       can see the current values of the state vars, but not modify them
//
// FIXME Make a new object, AtomGroupRegistry, to manage AtomGroups?

#pragma once
#ifndef SIMULATION_STATE_H
#define SIMULATION_STATE_H

// Standard headers
#include <array>
#include <vector>

// Project headers
#include "AtomGroup.h"
#include "CommonTypes.h"
#include "InputParser.h"
#include "SimulationBox.h"
#include "Topology.h"

class SimulationState 
{
 public:
	using RvecArray      = CommonTypes::RvecArray;
	using SelectionFlags = Topology::TargetFlags;

	SimulationState();

	//----------------------------//
	//----- Accessor Methods -----//
	//----------------------------//

	// Atom_positions list
	RvecArray& access_atom_positions() { return atom_positions_; }
	const RvecArray& get_atom_positions() const { return atom_positions_; }
	
	// SimulationBox
	SimulationBox& access_simulation_box() { return simulation_box_; }
	const SimulationBox& get_simulation_box() const { return simulation_box_; }

	// Topology
	Topology& access_topology() { return topology_; }
	const Topology& get_topology() const { return topology_; }

	// TODO Add an AtomGroup by passing an AtomGroup ParameterPack

	// Create a new AtomGroup and return its index in the array of groups
	int addAtomGroup(
		const ParameterPack& atom_group_pack
	);
	int addAtomGroup(
		const std::string& name,
		const std::vector<std::string>& selection_tokens,
		const std::vector<int>& global_indices,
		const SelectionFlags& selection_flags
	);
	int addAtomGroup();  // make an empty group with a default-format name

	// Returns 'true' if a group with the given name exists
	bool atomGroupExists(const std::string& name) const;

	// Generic AtomGroup access
	// - By index (throw an exception if the group cannot be found)
	AtomGroup& access_atom_group(const int index);
	const AtomGroup& get_atom_group(const int index) const;
	// - Get all groups
	std::vector<AtomGroup>& access_atom_groups() { return atom_groups_; }
	const std::vector<AtomGroup>& get_atom_groups() const { return atom_groups_; }

	// Simulation step and time (TODO: use longer int for step?)
	int    get_step() const { return step_; }
	double get_time() const { return time_; }

	void set_step(const int step)    { step_ = step; }
	void set_time(const double time) { time_ = time; }

	// Whether to print debugging messages
	bool debug_mode() const {
		return is_debug_mode_;
	}
	void set_debug_mode(const bool is_debug_mode) { 
		is_debug_mode_ = is_debug_mode; 
	}

 private:
	// Current step
	int    step_;
	double time_;  // [ps]

	// List of *all* atom_positions (FIXME MPI)
	// - Standalone: only available on master rank 
	// - PLUMED mode: TODO
	RvecArray atom_positions_;

	// SimulationBox: holds box matrix and manages periodic boundary conditions (PBCs)
	// - Use to compute minimum image distances
	SimulationBox simulation_box_;

	Topology topology_;

	// Organize atoms into groups
	std::vector<AtomGroup> atom_groups_;

	bool is_debug_mode_ = false;
};

#endif /* SIMULATION_STATE_H */
