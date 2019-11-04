// AtomGroup - Holds global *and* local information about a group of atoms

// NOTES:
// - A "local" atom is typically defined as an atom that is either:
//    (1) in the processor's DD cell (processor "owns" this atom); or
//    (2) in shell 1 of the processor's DD cell; or
//    (3) in shell 2 of the processor's DD cell
//   - "Shell" semantics are the same as for ProbeVolumes
// - Domain decomposition 
//   - Currently, all MPI ranks have copies all atom positions in the group
//   - Each rank uses a simple algorithm to figure out which atoms it's responsible
//     for using a simple spatial decomposition scheme (see CellGrid)
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
#include "Region.h"
#include "Topology.h"
#include "utils.h"

class AtomGroup
{
 public:
	const int DIM_ = CommonTypes::DIM_;
	using Real3          = CommonTypes::Real3;
	using RvecArray      = CommonTypes::RvecArray;
	using SelectionFlags = Topology::TargetFlags;
	using RegionEnum     = Region::RegionEnum;

	// Group atom which is "nearby" relative to this rank
	// - Three cases (see NOTES)
	struct NearbyAtom {
		int        group_index;       // index in the AtomGroup
		bool       is_in_local_cell;
		RegionEnum region;            // note: 'vtilde' => 'v' here
	};


	//----- Constructors -----//

	AtomGroup();

	AtomGroup(
		const std::string& name,
		const std::vector<std::string>& selection_tokens,
		const std::vector<int>& global_atom_indices,
		const SelectionFlags& selection_flags
	);


	//----- Get Functions -----//

	const std::string& get_name() const { return name_; }
	const std::string& get_selection_string() const { return selection_string_; }

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

	// Access single-atom properties by group index
	int getGlobalIndexUsingGroupIndex(const int group_index) const {
		return global_atom_indices_[group_index];
	}
	const Real3& getPositionUsingGroupIndex(const int group_index) const {
		return atom_positions_[group_index];
	}
	int get_global_index(const int group_index) const {
		return global_atom_indices_[group_index];
	}
	const Real3& get_position(const int group_index) const {
		return atom_positions_[group_index];
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
	void set_atom_positions(const std::vector<Real3>& atom_positions);


 private:
	// Selection
	std::string name_;
	std::vector<std::string> selection_tokens_;
	SelectionFlags selection_flags_;

	// Concatenated form of selection_tokens_
	std::string selection_string_;


	//----- All Atoms in Group -----//

	// Global indices of all atoms in the group
	// - Size: [num_atoms_in_group]
	std::vector<int> global_atom_indices_;

	// Vector with positions of all atoms in the group
	// - Size: [num_atoms_in_group]
	std::vector<Real3> atom_positions_;



	//----- Local Atoms -----//

	// CellGrid cell index of local atoms for this rank
	// - Same as MPI rank
	//int local_cell_index_;

	std::vector<NearbyAtom> nearby_atoms_;
};

#endif /* ATOM_GROUP_H */
