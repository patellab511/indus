// SimulationState
//  - Unified place for simulation state variables to live
//    (e.g. atom_positions, box, etc.)
//
//  - Simplifies sending all this information to objects that need it 
//     - No need to explicitly pass around the latest set of atom_positions
//     - Get functions control access: an object that receives a const SimulationState&
//       can see the current values of the state vars, but not modify them

#pragma once
#ifndef SIMULATION_STATE_H
#define SIMULATION_STATE_H

// Standard headers
#include <array>
#include <vector>

// Project headers
#include "AtomGroup.h"
#include "CommonTypes.h"
#include "MpiCommunicator.h"
#include "SimulationBox.h"
#include "Topology.h"

class SimulationState 
{
 public:
	using RvecArray = CommonTypes::RvecArray;

	SimulationState();

	//----------------------------//
	//----- Accessor Methods -----//
	//----------------------------//

	// Atom_positions list
	RvecArray& access_atom_positions();
	const RvecArray& get_atom_positions() const;
	
	// SimulationBox
	SimulationBox& access_simulation_box(); 
	const SimulationBox& get_simulation_box() const;

	// Topology
	Topology& access_topology();
	const Topology& get_topology() const;

	// Target AtomGroup
	AtomGroup& access_target_atom_group();
	const AtomGroup& get_target_atom_group() const;

	int get_step() const { return step_; }
	void set_step(const int step) { step_ = step; }

	double get_time() const { return time_; }
	void set_time(const double time) { time_ = time; }

 private:
	// Current step
	int    step_;
	double time_;  // [ps]

	// List of "all" atom_positions (FIXME MPI)
	// - Standalone MDAnalysis++: only available on master rank 
	// - PLUMED mode:
	RvecArray atom_positions_;

	// SimulationBox: holds box matrix and manages periodic boundary conditions (PBCs)
	// - Use to compute minimum image distances
	SimulationBox simulation_box_;

	Topology topology_;

	// Atoms which are the targeted by a bias
	AtomGroup target_atom_group_;
};

#endif /* SIMULATION_STATE_H */
