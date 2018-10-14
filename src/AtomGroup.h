// AtomGroup - Holds global *and* local information about a group of atoms

// NOTES:
// - A "local" atom is typically defined as an atom that is either:
//    (1) in the processor's DD cell (processor "owns" this atom); or
//    (2) in shell 1 of the processor's DD cell; or
//    (3) in shell 2 of the processor's DD cell
//   - "Shell" semantics are the same as for ProbeVolumes
//
// TODO instead of copying over all positions, could instead have a wrapper
//      around a ptr to the positions array owned by the driver

#ifndef ATOM_GROUP_H
#define ATOM_GROUP_H

// Standard headers
#include <exception>
#include <sstream>
#include <string>
#include <vector>

// Project headers
#include "CommonTypes.h"
#include "MpiCommunicator.h"

class AtomGroup
{
 public:
	const int DIM_  = CommonTypes::DIM_;
	using Real3     = CommonTypes::Real3;
	using RvecArray = CommonTypes::RvecArray;

	// Group atom which is "nearby" relative to this rank
	// - Three cases (see NOTES)
	struct NearbyAtom {
		int  group_index;  // index in the AtomGroup
		bool is_in_local_cell;
		bool is_in_local_cell_shell_1;
		bool is_in_local_cell_shell_2;
	};


	//----- Constructors -----//

	AtomGroup();

	AtomGroup(
		const std::vector<int>& global_atom_indices
	);


	//----- Get Functions -----//

	int get_num_atoms() const {
		return static_cast<int>( global_atom_indices_.size() );
	}
	int get_num_nearby_atoms() const {
		return static_cast<int>( nearby_atoms_.size() );
	}


	const std::vector<int>& get_global_atom_indices() const {
		return global_atom_indices_;
	}
	std::vector<int>& access_global_atom_indices() {
		return global_atom_indices_;
	}

	int getGlobalIndexUsingGroupIndex(const int group_index) const {
		return global_atom_indices_[group_index];
	}


	const std::vector<Real3>& get_atom_positions() const {
		return atom_positions_;
	}
	std::vector<Real3>& access_atom_positions() {
		return atom_positions_;
	}

	
	const std::vector<NearbyAtom>& get_nearby_atoms() const {
		return nearby_atoms_;
	}
	std::vector<NearbyAtom>& access_nearby_atoms() {
		return nearby_atoms_;
	}


	//----- Set Functions -----//

	void set_global_atom_indices(
		const std::vector<int>& global_atom_indices
	);

	// These functions check that the input is consistent with the number of atoms
	// in the group [determined by global_atom_indices_.size()]
	void set_atom_positions(const std::vector<Real3> atom_positions);
	//void set_atom_positions_using_global_positions_array(const RvecArray& global_positions_array);


 private:
	//----- All Atoms in Group -----//

	// Global indices of all atoms in the group
	// - Size: [num_atoms_in_group]
	std::vector<int> global_atom_indices_;

	// Vector with positions of atoms in the group
	// - Size: [num_atoms_in_group]
	std::vector<Real3> atom_positions_;


	//----- Local Group Atoms -----//

	// CellGrid cell index of local atoms for this rank
	// - Same as MPI rank
	int local_cell_index_;

	std::vector<NearbyAtom> nearby_atoms_;
};

#endif /* ATOM_GROUP_H */
