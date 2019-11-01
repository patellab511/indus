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
#include "utils.h"

// Project headers
#include "AtomGroup.h"
#include "BoundingBox.h"
#include "CommonTypes.h"
#include "GptlWrappers.h"
#include "MpiCommunicator.h"
#include "OpenMP.h"
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

	// A "cell" is a (triclinic) subvolume within a larger (triclinic) region of interest
	// - Each cell has up to 26 neighboring cells.
	struct Cell {
		Int3  indices;       // indices in the 3D grid
		Box   box_matrix;    // box matrix for this cell
		Real3 x_lower;       // position of lower-left corner
		Real3 x_upper;       // position of upper-right corner
		Real3 x_center;      // center point

		// Indices of neighboring cells in "cells_" array
		std::vector<int> neighbors;

		// For fast grid search with 3D PBCs: when a particle is in this cell,
		// search the cells with these indices for neighbors
		std::vector<int> pbc_grid_search_neighbors;


		//----- Managing neighbors -----//

		void reset() {
			neighbors.resize(0);                  neighbors.reserve(MAX_NUM_NEIGHBORING_CELLS);
			pbc_grid_search_neighbors.resize(0);  pbc_grid_search_neighbors.reserve(MAX_NUM_NEIGHBORING_CELLS);
		};

		// TODO Move to utilities header?
		template<typename T>
		static void removeDuplicates(std::vector<T>& vec) {
			std::sort( vec.begin(), vec.end() );
			vec.erase( std::unique(vec.begin(), vec.end()), vec.end() );
		}

		void removeDuplicates() {
			removeDuplicates( neighbors );
			removeDuplicates( pbc_grid_search_neighbors );
		}
	};


	//------ Core Functions -----//

	CellGrid(
		const SimulationState& simulation_state,
		MpiCommunicator& mpi_communicator
	);

	virtual ~CellGrid();

	// Returns true if "box" is orthorhombic (off-diagonal elements all zero)
	bool is_orthorhombic_box(const Box& box) const;

	// Returns true if "box" has no non-negative components
	// - For now, it is also assumed that a valid box must be orthorhombic as well
	bool is_box_valid(const Box& box) const;

	// Partitions the region defined by "bounding_box_matrix_" and
	// "bounding_box_offset_" into cells
	void makeCells();

	// Find the cell that contains position 'x'
	// - If a cell for 'x' cannot be found, throws an exception
	// THREADSAFE
	void findCell(
		const Real3& x,
		// Output
		int&  cell_index,   // linear cell index
		Int3& grid_indices  // indices of the cell in the grid
	) const;


	//----- Get/Set Functions -----//

	// Sets the CellGrid region of interest to be the entire simulation box
	void setCellGridRegion();

	// Sets the cell region region to be the given box
	void setCellGridRegion(const BoundingBox& bounding_box);

	// Directly set the grid dimensions
	// TODO delete/move?
	void set_grid_dimensions(const Int3& grid_dimensions) {
		grid_dimensions_ = grid_dimensions;

		this->makeCells();
	}

	// Called in makeCells() to update template cell and/or grid dimensions,
	// if one of them is inferred from the other in a derived class
	virtual void updateGrid() {};


	//----- Get Functions -----//

	// Get linear cell index from indices of cell in the grid
	int get_cell_index(const Int3& grid_indices) const;

	// Access stored cells
	const std::vector<Cell>& get_cells() const { return cells_; }

	// Access stored cells
	const Int3& get_grid_dimensions() const { return grid_dimensions_; }


	//----- Debugging -----//

	void checkGridDimensions() const {
		for ( int d=0; d<DIM_; ++d ) {
			if ( grid_dimensions_[d] < 1 ) {
				std::stringstream err_ss;
				err_ss << "Error in " << FANCY_FUNCTION << "\n"
							 << "  Number of cells along " << CommonTypes::axis_names[d] << " is " 
								 << grid_dimensions_[d] << ".\n";
				throw std::runtime_error( err_ss.str() ); 
			}
		}
	}

	// Printing cells
	friend std::ostream& operator<<(std::ostream& os, const CellGrid& cell_grid);
	void printCells(std::ostream& os) const;

	// For running tests
	// - TODO: Find a more elegant way to do this
	friend int main(int argc, char* argv[]);

 protected:
	// Handle to simulation state
	const SimulationState& simulation_state_;

	// Need box to compute minimum image distances
	const SimulationBox& simulation_box_;

	// MPI communicator
	MpiCommunicator& mpi_communicator_;

	// Box matrix describing the "model" individual cell
	// - Most cells created from this template will have the same size and shape. 
	//   However, cells near the boundaries of the region of interest may be
	//   made larger to make sure the grid of cells covers the entire region
	Box template_cell_box_matrix_;

	// Number of cells along x, y, and z
	Int3 grid_dimensions_;     

	// Region covered by cell grid
	BoundingBox bounding_box_;

	// Cell grid
	std::vector<Cell> cells_;

	// Upper and lower bounds along each linear axis
	std::array<std::vector<double>, DIM_> cell_lower_bounds_, cell_upper_bounds_;

	// Default bounding box: entire simulation box
	BoundingBox constructDefaultBoundingBox() const {
		return BoundingBox(simulation_box_);
	}

	// For debugging: thoroughly checks the cells
	void checkCells() const;

	// Tests the construction of cells for grids of a number of different sizes
	// - Only use on a dummy object in debug mode!
	void test_grid_dimensions();

 private:
	mutable GPTL::Timer make_cells_timer_     = GPTL::Timer("CellGrid::makeCells");
	mutable GPTL::Timer make_cells_omp_timer_ = GPTL::Timer("CellGrid::makeCells_omp");
};

#endif /* CELL_GRID_H */
