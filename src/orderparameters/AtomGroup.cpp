#include "AtomGroup.h"

//----- Constructors -----//

AtomGroup::AtomGroup():
	name_("n/a"), selection_string_("n/a")
{}


AtomGroup::AtomGroup(
	const std::string& name, 
	const std::vector<std::string>& selection_tokens,
	const std::vector<int>& global_atom_indices,
	const SelectionFlags& selection_flags
):
	name_(name),
	selection_tokens_(selection_tokens),
	selection_flags_(selection_flags)
{
	selection_string_ = "";
	int num_tokens = selection_tokens_.size();
	for ( int i=0; i<num_tokens; ++i ) {
		if ( i > 0 ) { selection_string_ += " "; }
		selection_string_ += selection_tokens_[i];
	}

	set_global_atom_indices(global_atom_indices);
}


//----- Set Functions -----//


void AtomGroup::set_global_atom_indices(
	const std::vector<int>& global_atom_indices)
{
	global_atom_indices_ = global_atom_indices;

	// Preallocate memory
	int num_atoms_in_group = global_atom_indices_.size();
	atom_positions_.resize(num_atoms_in_group);
}


void AtomGroup::set_atom_positions(
	const std::vector<CommonTypes::Real3>& atom_positions)
{
	if ( atom_positions.size() == global_atom_indices_.size() ) {
		atom_positions_ = atom_positions;
	}
	else {
		std::stringstream err_ss;
		err_ss << "Group " << name_ << ": size of input positions array (" << atom_positions.size()
					 << ") does not match the number of atoms in the group ("
					 << global_atom_indices_.size() << ")\n";
		throw std::runtime_error( err_ss.str() );
	}
}


/*
void AtomGroup::set_atom_positions_using_global_positions_array(
	const CommonTypes::RvecArray& global_positions_array) 
{
	int num_atoms_total    = global_positions_array.size();
	int num_atoms_in_group = global_atom_indices_.size();

	// Ensure there's enough memory
	atom_positions_.resize(num_atoms_in_group);

	int global_index;
	for ( int k=0; k<num_atoms_in_group; ++k ) {
		global_index = global_atom_indices_[k];
		if ( global_index < num_atoms_total ) {
			for ( int d=0; d<DIM_; ++d ) {
				atom_positions_[k][d] = global_positions_array[global_index][d];
			}
		}
		else {
			std::stringstream err_ss;
			err_ss << "global index of atom in group (" << global_index 
						 << ") exceeds the global number of atoms ("
						 << num_atoms_total << ")\n";
			throw std::runtime_error( err_ss.str() );
		}
	}
}
*/
