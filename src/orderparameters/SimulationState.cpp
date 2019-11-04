#include "SimulationState.h"


SimulationState::SimulationState():
	step_(0), time_(0.0),
	simulation_box_(),
	topology_()
{}


// Create a new AtomGroup and return its index in the array of groups
int SimulationState::addAtomGroup(const ParameterPack& atom_group_pack)
{

	using KeyType = ParameterPack::KeyType;

	std::string name;
	atom_group_pack.readString("name", KeyType::Required, name);

	std::vector<std::string> selection_tokens;
	atom_group_pack.readVector("selection", KeyType::Required, selection_tokens);

	// Use Topology to parse selection
	// FIXME Pass flags to AtomGroup along with indices
	// - Have AtomGroup parse its own input?!
	std::vector<int> global_indices;
	Topology::TargetFlags selection_flags;
	topology_.getTargets(selection_tokens, global_indices, selection_flags);

	// Register the new group
	return addAtomGroup( name, selection_tokens, global_indices, selection_flags );
}


int SimulationState::addAtomGroup(
	const std::string& name, const std::vector<std::string>& selection_tokens,
	const std::vector<int>& global_indices, const SelectionFlags& selection_flags )
{
	// Make sure a group with this name does not already exist
	if ( atomGroupExists(name) ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "More than one AtomGroup was given the name \"" << name << "\"\n";
		throw std::runtime_error( err_ss.str() );
	}

	atom_groups_.push_back( AtomGroup(name, selection_tokens, global_indices, selection_flags) );
	int new_group_index = static_cast<int>(atom_groups_.size()) - 1;

	if ( debug_mode() ) {
		std::cout << "ADDED ATOMGROUP WITH NAME \"" << name << "\" [" << FANCY_FUNCTION << "]\n";
	}

	return new_group_index;
}


int SimulationState::addAtomGroup()
{
	std::stringstream default_name_ss;
	default_name_ss << "atom_group_" << atom_groups_.size();
	
	std::vector<std::string> dummy_tokens = {{ "n/a" }};
	std::vector<int> empty_indices(0);
	Topology::TargetFlags dummy_flags;

	return addAtomGroup( default_name_ss.str(), dummy_tokens, empty_indices, dummy_flags );
}


bool SimulationState::atomGroupExists(const std::string& name) const
{
	int num_groups = atom_groups_.size();
	for ( int n=0; n<num_groups; ++n ) {
		if ( atom_groups_[n].get_name() == name ) {
			return true;
		}
	}

	return false;
}


// AtomGroup management and access
// - Throws an exception if the group cannot be found
AtomGroup& SimulationState::access_atom_group(const int index)
{
	if ( index >= static_cast<int>(atom_groups_.size()) ) {
		std::stringstream err_ss;
		err_ss << "Error in SimulationState::access_atom_group: "
		       << "Atom group index " << index << " is out of bounds";
		throw std::runtime_error( err_ss.str() );
	}

	return atom_groups_[index];
}


const AtomGroup& SimulationState::get_atom_group(const int index) const
{
	if ( index >= static_cast<int>(atom_groups_.size()) ) {
		std::stringstream err_ss;
		err_ss << "Error in SimulationState::get_atom_group: "
		       << "Atom group index " << index << " is out of bounds";
		throw std::runtime_error( err_ss.str() );
	}

	return atom_groups_[index];
}
