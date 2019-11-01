// CellGrid_DomainDecomposition
// - Partitions an orthorhombic domain into cells of equal size
//   and shape
//   - The total number of cells is equal to the number of MPI ranks
// - Contains methods for assigning particles to domains, as well as identifying
//   particles which are close to (but not inside) these domains
//   - This is useful for domain decomposition with Steinhardt parameters
// - Implementation Notes
//   - The template cell is determined automatically using the specified 
//     grid dimensions and the bounding box

#pragma once
#ifndef CELL_GRID_DOMAIN_DECOMPOSITION_H
#define CELL_GRID_DOMAIN_DECOMPOSITION_H

#include "CellGrid.h"
#include "BoundingBox.h"
#include "ProbeVolume_Box.h"   // for finding atoms

class CellGrid_DomainDecomposition : public CellGrid
{
 public:
	using RegionEnum = ProbeVolume::RegionEnum;

	CellGrid_DomainDecomposition(
		const Int3&  grid_dimensions, 
		const double width_shell_1,
		const double width_shell_2, 
		const SimulationState& simulation_state,
		MpiCommunicator& mpi_communicator
	);

	// Updates the template cell using stored grid dimensions
	// - Assumes that the bounding box has been set
	virtual void updateGrid() override;

	// Takes a local AtomGroup (which contains copies of the positions of all its atoms) and
	// determines which atoms are "nearby" the DD cell with the given index
	// - Updates the local atoms list in the given AtomGroup
	void findNearbyAtoms(
		const int  cell_index,  // linear index of cell to search
		AtomGroup& atom_group   // atoms to check
	) const;


	void setShellWidths(const double width_shell_1, const double width_shell_2)
	{
		width_shell_1_ = width_shell_1;
		width_shell_2_ = width_shell_2;

		if ( local_cell_probe_volume_ptr_ != nullptr ) {
			// Set shell widths for probe box representing local DD cell
			local_cell_probe_volume_ptr_->setShellWidths(width_shell_1_, width_shell_2_);
		}
		else {
			throw std::runtime_error("CellGrid_DomainDecomposition: missing local cell probe volume");
		}
	}

	void getShellWidths(double& width_shell_1, double& width_shell_2) const
	{
		width_shell_1 = width_shell_1_;
		width_shell_2 = width_shell_2_;
	}

	const BoundingBox& getLocalCellBoundingBox() const {
		return local_cell_probe_volume_ptr_->getBoundingBox();
	}


	//----- Debugging -----//

	// Check that the given AtomGroup was correctly partitioned
	void printAtomGroupDomainDecompositionCheck(
		std::ostream&    os, 
		const AtomGroup& atom_group
	) const;


 private:
	// Shells widths (for working with probe volume shells)
	double width_shell_1_ = 0.0;
	double width_shell_2_ = 0.0;

	// Probe box which represents the local cell
	std::unique_ptr<ProbeVolume_Box> local_cell_probe_volume_ptr_ = nullptr;

	// Buffers for identifying NearbyAtoms with multiple threads
	mutable std::vector<std::vector<AtomGroup::NearbyAtom>> nearby_atom_buffers_;

	// Timers
	mutable GPTL::Timer find_nearby_atoms_timer_ = GPTL::Timer("CellGrid_DD::findNearbyAtoms");
	mutable GPTL::Timer dd_omp_timer_            = GPTL::Timer("DD_findNearbyAtoms_omp");
	mutable GPTL::Timer dd_omp_combine_timer_    = GPTL::Timer("DD_findNearbyAtoms_omp_combine");
};

#endif /* CELL_GRID_DOMAIN_DECOMPOSITION_H */
