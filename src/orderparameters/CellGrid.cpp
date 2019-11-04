/* CellGrid.cpp
 *   
 */


#include "CellGrid.h"


CellGrid::CellGrid(const SimulationState& simulation_state, MpiCommunicator& mpi_communicator):
	simulation_state_( simulation_state ),
	simulation_box_( simulation_state.get_simulation_box() ),
	mpi_communicator_( mpi_communicator ),
	bounding_box_( constructDefaultBoundingBox() )
{
	// Default partition: one cell encompassing the whole simulation box
	grid_dimensions_.fill(1);
	template_cell_box_matrix_ = bounding_box_.get_box_matrix();

	this->makeCells();
}


CellGrid::~CellGrid() {}


bool CellGrid::is_orthorhombic_box(const Box& box) const
{
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( a != b and box[a][b] != 0.0 ) {
				return false;
			}
		}
	}

	return true;
}


// Imposes an orthorhombic box
bool CellGrid::is_box_valid(const Box& box) const
{
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( box[a][b] < 0.0 ) {
				return false;
			}
		}
	}

	if ( not is_orthorhombic_box(box) ) {
		return false;
	}

	return true;
}


// Sets the CellGrid region of interest to be the entire simulation box
void CellGrid::setCellGridRegion() 
{
	setCellGridRegion( constructDefaultBoundingBox() );
}

// Sets the cell region region to be the given box
void CellGrid::setCellGridRegion(const BoundingBox& bounding_box)
{
	bounding_box_ = bounding_box;
	this->makeCells();  // update cells
}


// Get linear cell index from indices of cell in the grid
int CellGrid::get_cell_index(const Int3& grid_indices) const
{
	return grid_indices[X_DIM]*grid_dimensions_[Y_DIM]*grid_dimensions_[Z_DIM] +
	       grid_indices[Y_DIM]*grid_dimensions_[Z_DIM] + 
         grid_indices[Z_DIM];
}


void CellGrid::makeCells()
{
	make_cells_timer_.start();

	// Currently this code assumes orthorhombic boxes and subvolumes
	const Box& bounding_box_matrix = bounding_box_.get_box_matrix();
	if ( not is_box_valid( bounding_box_matrix ) ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  Bounding box for region of interest is invalid.\n"
		       << "\n"
		       << "bounding_box_matrix = [";
		
		for ( int a=0; a<DIM_; ++a ) {
			err_ss << " [";
			for ( int b=0; b<DIM_; ++b ) { 
				err_ss << " " << bounding_box_matrix[a][b];
			}
			err_ss << " ]\n";
		}
		err_ss << " ]\n";

		throw std::runtime_error( err_ss.str() );
	}

	// Update the grid dimensions and/or template cell
	// TODO change order of things?
	updateGrid();

	// Determine the number of cells
	checkGridDimensions();
	int num_cells = 1;
	for ( int d=0; d<DIM_; ++d ) {
		num_cells *= grid_dimensions_[d];
	}

	// If there are fewer than 3 cells along any direction, the algorithm below will
	// produce a list which includes duplicates: remove them as needed
	bool duplicates_possible = false;
	for ( int d=0; d<DIM_; ++d ) {
		if ( grid_dimensions_[d] < 3 ) {
			duplicates_possible = true;
			break;
		}
	}

	// Reset array of cells
	cells_.resize(num_cells);
	for ( int k=0; k<num_cells; ++k ) {
		cells_[k].reset();
	}

	// Upper and lower edges of the cells along each dimension
	// - Compute them this way so that the resulting partition leaves no
	//   gaps due to roundoff error
	Real3 x_lower = bounding_box_.get_x_lower();
	Real3 x_upper = bounding_box_.get_x_upper();
	for ( int d=0; d<DIM_; ++d ) {
		cell_lower_bounds_[d].resize( grid_dimensions_[d] );
		cell_upper_bounds_[d].resize( grid_dimensions_[d] );

		// First cell is flush with the bounding box corner
		cell_lower_bounds_[d][0] = x_lower[d];

		for ( int i=1; i<grid_dimensions_[d]; ++i ) {
			double edge = cell_lower_bounds_[d][i-1] + template_cell_box_matrix_[d][d];
			cell_upper_bounds_[d][i-1] = edge;
			cell_lower_bounds_[d][i]   = edge;
		}

		// Uppermost bound extends to the boundary of the gridded region
		cell_upper_bounds_[d].back() = x_upper[d];
	}

	// Create each cell
	// - TODO: move to makeCell function?
	int   k, neighbor_index;   // indices in flattened cell array ("linear" indices)
	Int3  neighbor_indices, offset, grid_indices, actual_shift;
	const int num_cells_x = grid_dimensions_[X_DIM];
	const int num_cells_y = grid_dimensions_[Y_DIM];
	const int num_cells_z = grid_dimensions_[Z_DIM];
	#pragma omp parallel for collapse(3) \
		private(k, neighbor_index, neighbor_indices, offset, grid_indices, actual_shift) \
		shared(num_cells)
	for ( int kx=0; kx<num_cells_x; ++kx ) {
		for ( int ky=0; ky<num_cells_y; ++ky ) {
			for ( int kz=0; kz<num_cells_z; ++kz ) {
				make_cells_omp_timer_.start();
				grid_indices = {{ kx, ky, kz }};

				// Corresponding linear index in array of cells
				k = get_cell_index( grid_indices );
				cells_[k].indices = grid_indices;

				// Key points in cell (assumes orthorhombic box)
				for ( int d=0; d<DIM_; ++d ) {
					// Lower-left corner
					cells_[k].x_lower[d] = cell_lower_bounds_[d][ grid_indices[d] ];

					// Upper-right corner
					cells_[k].x_upper[d] = cell_upper_bounds_[d][ grid_indices[d] ];

					// Center
					cells_[k].x_center[d] = 0.5*( cells_[k].x_upper[d] + cells_[k].x_lower[d] );
				}

				// Cell box matrix (assumes orthorhombic box)
				for ( int a=0; a<DIM_; ++a ) {
					for ( int b=0; b<DIM_; ++b ) {
						if ( a == b ) {
							cells_[k].box_matrix[a][a] = cells_[k].x_upper[a] - 
																										cells_[k].x_lower[a];
						}
						else {
							cells_[k].box_matrix[a][b] = 0.0;
						}
					}
				}

				//----- Neighboring Cells -----//

				// a, b, and c index over the cells which are -1, 0, and +1 in each
				// direction relative to this cell
				for ( int a=-1; a<=1; ++a ) {
					offset[X_DIM] = a;
					for ( int b=-1; b<=1; ++b ) {
						offset[Y_DIM] = b;
						for ( int c=-1; c<=1; ++c ) {
							offset[Z_DIM] = c;

							// Apply PBCs
							// - TODO: for grids that aren't close to the periodic boundaries of
							//   the box, this will lead to unnecessary checking between cells
							//   at opposite edges of the grid. Maybe check the distance between 
							//   cell centers to eliminate this?
							for ( int d=0; d<DIM_; ++d ) {
								neighbor_indices[d] = cells_[k].indices[d] + offset[d];
								if ( neighbor_indices[d] >= grid_dimensions_[d] ) {
									neighbor_indices[d] = 0;
								}
								else if ( neighbor_indices[d] < 0 ) {
									neighbor_indices[d] = grid_dimensions_[d] - 1;
								}
								actual_shift[d] = neighbor_indices[d] - cells_[k].indices[d];
							}

							// Corresponding linear index
							neighbor_index = get_cell_index( neighbor_indices );

							// Add the following check to handle cases where the cell maps back to
							// itself under PBCs
							if ( neighbor_index != k ) {
								cells_[k].neighbors.push_back( neighbor_index );

								// In the typical case for grid search under PBCs, the central cell 'k' (*) searches
								// itself and half (13 out of 26) of its neighbors (o)
								//
								//               c = -1          c = 0:          c = +1  
								//            ___ ___ ___     ___ ___ ___     ___ ___ ___
								//        +1 |___|___|___|   |_o_|_o_|_o_|   |_o_|_o_|_o_|   
								//    b =  0 |___|___|___|   |___|_*_|_o_|   |_o_|_o_|_o_|   
								//        -1 |___|___|___|   |___|___|___|   |_o_|_o_|_o_|   
								//            -1   0  +1
								//                 a
								bool is_grid_search_neighbor = false;
								if ( (c == 1) or 
										 ((c == 0) and (b >= 0) and (not ((a == -1) and (b == 0))))
								) {
									is_grid_search_neighbor = true;
								}

								// Special cases: grid is small along an axis
								// TODO: Possible to make more concise?
								for ( int d=0; d<DIM_; ++d ) {
									if ( (grid_dimensions_[d] == 2) and
											 ( ((cells_[k].indices[d] == 0) and (offset[d] < 0)) or
												 ((cells_[k].indices[d] == 1) and (offset[d] > 0)) )
									) {
										is_grid_search_neighbor = false;
									}
									else if ( (grid_dimensions_[d] == 1) and
														( offset[d] != 0 )
									) {
										is_grid_search_neighbor = false;
									}
								}

								if ( is_grid_search_neighbor ) {
									cells_[k].pbc_grid_search_neighbors.push_back( neighbor_index );
								}
							}
						}
					}
				}  // end loops over neighboring cells

				// Sort indices (for convenience)
				//std::sort( cells_[k].neighbors.begin(), cells_[k].neighbors.end() );

				if ( duplicates_possible ) {
					cells_[k].removeDuplicates();
				}

				make_cells_omp_timer_.stop();
			}  // z
		}  // y
	}  // x  (end pragma omp parallel for)
	make_cells_timer_.stop();

	if ( simulation_state_.debug_mode() ) {
		checkCells();
	}

	return;
}


// TODO change this so that it doesn't necessarily throw if a position is outside the grid
void CellGrid::findCell(
	const Real3& x, int& cell_index, Int3& grid_indices
) const
{
	const Real3& x_lower = cells_.front().x_lower;  // bottom-left corner of the grid
	const Real3& x_upper = cells_.back().x_upper;   // top-rightcorner of the grid
	for ( int d=0; d<DIM_; ++d ) {
		if ( x[d] < x_lower[d] or x[d] > x_upper[d] ) {
			// Completely outside region of interest
			grid_indices.fill(-1);
			cell_index = -1;
			break;
		}

		// Most cells are boxes with size and shape given by "template_cell_box_matrix_"
		// - The ones at the upper edges along x, y, and z may be larger, but the following
		//   expression below will still work for these
		// - Assumes orthorhombic box and cells
		grid_indices[d] = static_cast<int>( 
				std::floor( (x[d] - x_lower[d])/template_cell_box_matrix_[d][d] ) 
		);

		// Edge cases (the above expression can suffer from roundoff error)
		if ( grid_indices[d] < 0 and x[d] >= x_lower[d]) {
			grid_indices[d] = 0;
		}
		else if ( grid_indices[d] >= 1 and x[d] < cell_lower_bounds_[d][ grid_indices[d] ] ) {
			// Atom is actually one cell down
			grid_indices[d] -= 1;
		}
		else if ( grid_indices[d] < (grid_dimensions_[d] - 1) and
		          x[d] > cell_upper_bounds_[d][ grid_indices[d] ]
		) {
			// Atom is actually one cell up
			grid_indices[d] += 1;
		}
		else if ( grid_indices[d] == grid_dimensions_[d] ) {
			// If the particle is at the upper bound along this direction, put it into the
			// cell right next to the edge
			grid_indices[d] = grid_dimensions_[d] - 1;
		}
	}

	// Store atom index in appropriate cell
	cell_index = get_cell_index(grid_indices);
	if ( cell_index < 0 or cell_index > static_cast<int>(cells_.size()) ) {
		// Invalid grid index (TODO Move to an exception)
		std::stringstream err_ss;
		err_ss << "(rank " << mpi_communicator_.get_rank() << ") "
		       << "CellGrid::findCell - could not find cell containing position {";
		for ( int d=0; d<DIM_; ++d ) { err_ss << " " << x[d]; }
		err_ss << " }\n";
		err_ss << "  You probably need to change the cell grid boundaries.\n";
		printCells(err_ss);
		
		throw std::runtime_error( err_ss.str() );
	}
}


std::ostream& operator<<(std::ostream& os, const CellGrid& cell_grid) 
{
	cell_grid.printCells(os);
	return os;
}


void CellGrid::printCells(std::ostream& os) const
{
	// Print simulation box for reference
	const Box& box_matrix = simulation_box_.get_box_matrix();
	os << "Simulation box";
	os << "  h = {\n";
	for ( int a=0; a<DIM_; ++a ) {
		os << "  ";
		for ( int b=0; b<DIM_; ++b ) {
			if ( b > 0 ) { os << ","; }
			os << " " << box_matrix[a][b];
		}
		os << "\n";
	}
	os << "  }\n";

	const Real3& x_lower = bounding_box_.get_x_lower();
	const Real3& x_upper = bounding_box_.get_x_upper();
	const Box& bounding_box_matrix = bounding_box_.get_box_matrix();

	// Region of interest
	os << "Region of interest\n";
	os << "  bounding_box_matrix = {\n";
	for ( int a=0; a<DIM_; ++a ) {
		os << "  ";
		for ( int b=0; b<DIM_; ++b ) {
			if ( b > 0 ) { os << ","; }
			os << " " << bounding_box_matrix[a][b];
		}
		os << "\n";
	}
	os << "  }\n";
	os << "  x_lower = {";
	for ( int d=0; d<DIM_; ++d ) {
		if ( d > 0 ) { os << ","; }
		os << " " << x_lower[d];
	}
	os << " }\n";
	os << "  x_upper = {";
	for ( int d=0; d<DIM_; ++d ) {
		if ( d > 0 ) { os << ","; }
		os << " " << x_upper[d];
	}
	os << " }\n";
	os << "  grid_dimensions = {";
	for ( int d=0; d<DIM_; ++d ) {
		if ( d > 0 ) { os << ","; }
		os << " " << grid_dimensions_[d];
	}
	os << " }\n";
	os << "  template_cell_box_matrix = {\n";
	for ( int a=0; a<DIM_; ++a ) {
		os << "  ";
		for ( int b=0; b<DIM_; ++b ) {
			if ( b > 0 ) { os << ","; }
			os << " " << template_cell_box_matrix_[a][b];
		}
		os << "\n";
	}
	os << "  }\n";

	// Description of each individual cell
	int num_cells = cells_.size();
	for ( int k=0; k<num_cells; ++k ) {
		os << "Cell " << k << "\n";

		os << "  x_lower =";
		for ( int d=0; d<DIM_; ++d ) {
			os << " " << cells_[k].x_lower[d];
		}
		os << "\n";

		os << "  x_upper =";
		for ( int d=0; d<DIM_; ++d ) {
			os << " " << cells_[k].x_upper[d];
		}
		os << "\n";

		os << "  x_center =";
		for ( int d=0; d<DIM_; ++d ) {
			os << " " << cells_[k].x_center[d];
		}
		os << "\n";

		os << "  box_matrix =\n";
		for ( int a=0; a<DIM_; ++a ) {
			os << "    ";
			for ( int b=0; b<DIM_; ++b ) {
				os << " " << cells_[k].box_matrix[a][b];
			}
			os << "\n";
		}

		os << "\n";
	}
}


void CellGrid::checkCells() const
{
	int num_cells = cells_.size();
	if ( num_cells == 0 ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  No cells are present\n";
		throw std::runtime_error( err_ss.str() );
	}

	std::vector<std::pair<int,int>> pairs_found;
	std::vector<int> pair_owners;
	for ( int k=0; k<num_cells; ++k ) {
		int num_neighbors = cells_[k].pbc_grid_search_neighbors.size();
		for ( int n=0; n<num_neighbors; ++n ) {
			for ( int m=n+1; m<num_neighbors; ++m ) {
				if ( cells_[k].pbc_grid_search_neighbors[n] == cells_[k].pbc_grid_search_neighbors[m] ) {
					std::stringstream err_ss;
					err_ss << "Error in " << FANCY_FUNCTION << "\n"
					       << "  Duplicate cells are present in the list of one cell's neighbors for grid search.\n";
					throw std::runtime_error( err_ss.str() );
				}
			}

			// Record each pair of neighboring cells, and ensure it only appears once
			int l = cells_[k].pbc_grid_search_neighbors[n];
			auto new_pair = std::make_pair( std::min(k,l), std::max(k,l) );  // enforce order for comparison
			auto find_it = std::find( pairs_found.begin(), pairs_found.end(), new_pair );
			if ( find_it == pairs_found.end() ) {
				pairs_found.push_back( new_pair );
				pair_owners.push_back( k );
			}
			else {
				int other_cell = pair_owners[ std::distance(pairs_found.begin(), find_it) ];

				std::stringstream err_ss;
				err_ss << "Error in " << FANCY_FUNCTION << "\n"
				       << "  duplicate neighboring cell pair\n";

				err_ss << "  k: indices =";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[k].indices[d]; }
				err_ss << "\n";

				err_ss << "  l: indices =";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[l].indices[d]; }
				err_ss << "\n";

				err_ss << "  grid dims:  ";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << grid_dimensions_[d]; }
				err_ss << "\n";

				err_ss << "  other cell: ";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[other_cell].indices[d]; }
				err_ss << "\n";

				std::cerr << err_ss.str() << std::endl;
				//throw std::runtime_error( err_ss.str() );
			}
		}
	}

	// Check that each pair expected is present
	for ( int k=0; k<num_cells; ++k ) {
		int num_neighbors = cells_[k].neighbors.size();
		for ( int n=0; n<num_neighbors; ++n ) {
			int l = cells_[k].neighbors[n];

			auto this_pair = std::make_pair( std::min(k,l), std::max(k,l) );  // enforce order for comparison
			if ( std::find( pairs_found.begin(), pairs_found.end(), this_pair ) == pairs_found.end() ) {
				std::stringstream err_ss;
				err_ss << "missing neighboring cell pair\n";

				err_ss << "  k: indices =";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[k].indices[d]; }
				err_ss << "\n";

				err_ss << "  l: indices =";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[l].indices[d]; }
				err_ss << "\n";

				err_ss << "  grid dims:  ";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << grid_dimensions_[d]; }
				err_ss << "\n";

				std::cerr << err_ss.str() << std::endl;
				//throw std::runtime_error( err_ss.str() );
			}
		}
	}
}


void CellGrid::test_grid_dimensions()
{
	if ( mpi_communicator_.get_size() != 1 
	     or ( not simulation_state_.debug_mode() )
	) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  Improper use of testing function";
		throw std::runtime_error( err_ss.str() );
	}

	// Test different grid sizes
	int max_num_cells = 5;
	std::array<int,DIM_> grid_dimensions;
	for ( grid_dimensions[X_DIM] = 1; grid_dimensions[X_DIM] < max_num_cells; ++grid_dimensions[X_DIM] ) {
		for ( grid_dimensions[Y_DIM] = 1; grid_dimensions[Y_DIM] < max_num_cells; ++grid_dimensions[Y_DIM] ) {
			for ( grid_dimensions[Z_DIM] = 1; grid_dimensions[Z_DIM] < max_num_cells; ++grid_dimensions[Z_DIM] ) {
				// This will invoke makeCells(), which determines the cell neighbor list
				set_grid_dimensions( grid_dimensions );
			}
		}
	}
}
