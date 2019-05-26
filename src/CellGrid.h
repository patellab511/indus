/* CellGrid.h
 *
 * ABOUT: Divides a specified region of the simulation into a grid of cells
 *  - The region of space to be partitioned is defined by member variables
 *    "bounding_box_matrix_" (size and shape) and "bounding_box_offset_"
 *    (location of bottom-left corner)
 *    - The vectors making up bounding_box_matrix_ are parallel to the 
 *      corresponding simulation box vectors
 *  - Uses: 
 *    - Partition atoms into regions of space for domain decomposition
 *    - Efficiently generate a neighbor list for a specified cutoff
 *
 * PROGRAMMING NOTES
 * DEVELOPMENT: (TODO)
 *  - Currently, it is assumed that all box matrices are orthorhombic
 *    - This includes the bounding box, individual cells, *and* the simulation box itself
 *    
 *
 * CITIATIONS 
 *
 *   
 */

#pragma once
#ifndef CELL_GRID_H
#define CELL_GRID_H

// Check whether PLUMED is defined using one of its preprocessor variables
#ifndef INDUS_STANDALONE_MODE
#define CELL_GRID_PLUMED_MODE
#endif

// Standard headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

// Project headers
#include "AtomGroup.h"
#include "CommonTypes.h"
#include "MpiCommunicator.h"
#include "ProbeVolume_Box.h" // for finding local atoms
#include "SimulationBox.h"
#include "SimulationState.h"

class CellGrid
{
 public:

	//----- Typedefs -----//

	static constexpr int MAX_NUM_NEIGHBORING_CELLS = 26;

  static const int DIM_  = CommonTypes::DIM_;  // Dimensionality of simulation 
  static const int X_DIM = CommonTypes::X_DIM; // x-axis index
  static const int Y_DIM = CommonTypes::Y_DIM; // y-axis index
  static const int Z_DIM = CommonTypes::Z_DIM; // z-axis index
  using rvec      = CommonTypes::Rvec;
	using RvecArray = CommonTypes::RvecArray;
  using Real3     = CommonTypes::Real3;
	using Int3      = CommonTypes::Int3;
	using Bounds3D  = CommonTypes::Bounds3D;
	using Box       = CommonTypes::Box;

	enum class PartitionType {
		LocalNeighborList,
		SimulationBoxDomainDecomposition
	};

	// A "cell" is a (triclinic) subvolume within a larger (triclinic) region of interest
	// - Each cell has up to 26 neighboring cells.
	struct Cell {
		Int3  indices;              // indices in the 3D grid
		Box   box_matrix;           // box matrix for this cell
		Real3 x_lower;              // position of lower-left corner
		Real3 x_upper;              // position of upper-right corner
		Real3 x_center;             // center point
		std::vector<int> neighbors; // indices of neighboring cells in "cells_" array
	};


	//------ Public Functions -----//

	// Constructor for PartitionType::LocalNeighborList
	// - Cells will (mostly) be cubes with the indicated side length 
	CellGrid(
		const double r_cutoff,
		const double cell_side_length,
		const SimulationState& simulation_state,
		MpiCommunicator& mpi_communicator
	);

	// Constructor for PartitionType::SimulationBoxDomainDecomposition
	// - Cells will be created according to the given grid dimensions, which
	//   *should* create a number of cells equal to the number of MPI ranks
	// - TODO check this vs. MpiCommunicator?
	CellGrid(
		const Int3&  grid_dimensions, 
		const double width_shell_1,
		const double width_shell_2, 
		const SimulationState& simulation_state,
		MpiCommunicator& mpi_communicator
	);

	~CellGrid();

	// Returns true if "box" is orthorhombic
	bool is_orthorhombic_box(const Box& box) const;

	// Returns true if "box" has no non-negative components
	// - For now, it is also assumed that a valid box must be orthorhombic as well
	bool is_box_valid(const Box& box) const;

	// Partitions the region defined by "bounding_box_matrix_" and
	// "bounding_box_offset_" into cells
	void makeCells();

	// Create a cell list for target atoms
	// - "cell_list" records which atoms are in each cell by TARGET atom index (not global)
	//   - In other words, it records the atom's index j in "target_atoms", not target_atoms[j]
	//   - If "target_atoms" is sorted, then the entries of cell_list will also be sorted
	// - Dimensions
	//   - cell_list = numCells x [num atoms in each cell (varies by cell)]
	void sortTargetAtomsIntoCells(
		const std::vector<Real3>& coords,
		const std::vector<int>&   target_atoms, // indices in coords array
		// Output
		std::vector<std::vector<int>>& cell_list,    // indices in TARGET atoms array
		std::vector<int>&              cell_indices  // cell in which each target atom resides
	);

	// Uses a cell search algorithm to construct a neighbor list based on r_cutoff_
	// and the internally stored cell partitions
	void makeNeighborList(
		const std::vector<Real3>& coords,
		const std::vector<int>&   target_atoms, // indices in coords array 
		// Output
		std::vector<std::vector<int>>& neighbor_list
	);

	// Takes a local AtomGroup (which contains copies of the positions of all its atoms) and
	// determines which atoms are "nearby" the DD cell with the given index
	// - Updates the local atoms list in the given AtomGroup
	void findNearbyAtoms(
		const int  cell_index,  // linear index of cell to search
		AtomGroup& atom_group   // atoms to check
	) const;


	//----- Get/Set Functions -----//

	// Sets the region covered by the cell grid to be the
	// simulation box
	void setCellGridRegion();

	// Sets the cell region region to be the given box
	void setCellGridRegion(
		const Real3& box_offset,  // position of lower-left corner
		const Box&   box_matrix   // shape of box
	);

	// Sets the region covered by the cell grid to encompass all
	// of the given positions
	void setCellGridRegion(
		const std::vector<Real3>& atom_positions, 
		const std::vector<int>&   target_atoms   // indices of atom positions to check
	);

	// Cube with side length cell_side_length
	// - Sets: template_cell_box_matrix_
	void setTemplateCell(const double cell_side_length);

	// Sets the template cell using the stored grid dimensions
	// - Assumes that the CellGrid region of interest has been set
	// - Sets: template_cell_box_matrix_
	void setTemplateCellUsingGridDimensions();

	// Sets the grid dimensions using the stored template cell
	// - Assumes that the CellGrid region of interest has been set
	// - Sets: grid_dimensions_
	void setGridDimensionsUsingTemplateCell();

	// Directly set the grid dimensions
	void set_grid_dimensions(const Int3& grid_dimensions) {
		grid_dimensions_ = grid_dimensions;
	}

	void setShellWidths(const double width_shell_1, const double width_shell_2)
	{
		width_shell_1_  = width_shell_1;
		width_shell_2_  = width_shell_2;

		if ( local_cell_probe_volume_ptr_ != nullptr ) {
			// Set shell widths for probe box representing local DD cell
			local_cell_probe_volume_ptr_->setShellWidths(width_shell_1_, width_shell_2_);
		}

	}

	double get_r_cutoff() const { return r_cutoff_; }

	// Get linear cell index from indices of cell in the grid
	int get_cell_index(const Int3& grid_indices) const;

	const std::vector<Cell>& get_cells() const { return cells_; }

	//----- Debugging -----//

	// Printing cells
	friend std::ostream& operator<<(std::ostream& os, const CellGrid& cell_grid);
	void printCells(std::ostream& os) const;

	// Check that the given AtomGroup was correctly partitioned
	void printAtomGroupDomainDecompositionCheck(
		std::ostream&    os, 
		const AtomGroup& atom_group
	);


 private:
	// Indicates the region of interest and how it is partitioned
	PartitionType partition_type_;

	// Cutoff for neighbor list creation 
	double r_cutoff_ = 0.0, r_cutoff_sq_ = 0.0;

	// Shells widths (for domain decomposition)
	double width_shell_1_  = 0.0;
	double width_shell_2_  = 0.0;

	// Probe box which represents the local cell (for domain decomposition)
	std::unique_ptr<ProbeVolume_Box> local_cell_probe_volume_ptr_ = nullptr;

	// Box matrix describing the "model" individual cell
	// - Most cells created from this template will have the same size and shape. 
	//   However, cells near the boundaries of the region of interest may be
	//   made larger to make sure the grid of cells covers the entire region
	Box template_cell_box_matrix_;

	// Region covered by cell grid
	Box   bounding_box_matrix_; // Size and shape of bounding box
	Real3 bounding_box_offset_; // Location of lower-left corner of bounding box

	// Cell grid
	std::vector<Cell> cells_;
	Int3 grid_dimensions_;     // number of cells along x, y, and z

	// Upper and lower bounds along each linear axis
	std::array<std::vector<double>, DIM_> cell_lower_bounds_;
	std::array<std::vector<double>, DIM_> cell_upper_bounds_;

	// Store these cell list variable lists internally to minimize reallocations
	std::vector<std::vector<int>> cell_list_;
	std::vector<int> cell_indices_; // cell in which each target atom resides

	// Handle to simulation state
	const SimulationState& simulation_state_;

	// Need box to compute minimum image distances
	const SimulationBox& simulation_box_;

	// MPI communicator
	MpiCommunicator& mpi_communicator_;
};

#endif /* CELL_GRID_H */
