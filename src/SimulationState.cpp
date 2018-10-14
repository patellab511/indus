#include "SimulationState.h"


SimulationState::SimulationState()
 : step_(0), time_(0.0),
   simulation_box_(),
   topology_()
{}

// Global atom_positions list
CommonTypes::RvecArray& SimulationState::access_atom_positions() { 
	return atom_positions_; 
}

const CommonTypes::RvecArray& SimulationState::get_atom_positions() const { 
	return atom_positions_; 
}

// SimulationBox
SimulationBox& SimulationState::access_simulation_box() { 
	return simulation_box_; 
}
const SimulationBox& SimulationState::get_simulation_box() const { 
	return simulation_box_; 
}

// Topology
Topology& SimulationState::access_topology() { 
	return topology_;
}

const Topology& SimulationState::get_topology() const { 
	return topology_; 
}

// Target AtomGroup
AtomGroup& SimulationState::access_target_atom_group() {
	return target_atom_group_;
}

const AtomGroup& SimulationState::get_target_atom_group() const {
	return target_atom_group_;
}
