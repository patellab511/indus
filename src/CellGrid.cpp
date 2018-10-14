/* CellGrid.cpp
 *   
 */


#include "CellGrid.h"


CellGrid::CellGrid(const double r_cutoff, const double cell_side_length, 
                   const SimulationState& simulation_state, MpiCommunicator& mpi_communicator
)
 : partition_type_( PartitionType::LocalNeighborList ), 
   r_cutoff_(r_cutoff),
   r_cutoff_sq_(r_cutoff*r_cutoff),
   simulation_state_( simulation_state ),
   simulation_box_( simulation_state.get_simulation_box() ),
   mpi_communicator_( mpi_communicator )
{
	// Cubic template cell
	setTemplateCell(cell_side_length);
}


CellGrid::CellGrid(const Int3&  grid_dimensions, const double width_shell_1, const double width_shell_2,
                   const double alpha_c_shells, const SimulationState& simulation_state,
                   MpiCommunicator& mpi_communicator)
 : partition_type_( PartitionType::SimulationBoxDomainDecomposition ),
   width_shell_1_( width_shell_1 ),
   width_shell_2_( width_shell_2 ),
   alpha_c_shells_( alpha_c_shells ),
   grid_dimensions_ ( grid_dimensions ),
   simulation_state_( simulation_state ),
   simulation_box_( simulation_state.get_simulation_box() ),
   mpi_communicator_( mpi_communicator )
{
	// TODO anything here? input checks?
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
	simulation_box_.getBoxMatrix(bounding_box_matrix_);

	for ( int d=0; d<DIM_; ++d ) {
		// Origin
		bounding_box_offset_[d] = 0.0;
	}
}


// Sets the cell region region to be the given box
void CellGrid::setCellGridRegion(const Real3& box_offset, const Box& box_matrix)
{
	bounding_box_offset_ = box_offset;
	bounding_box_matrix_ = box_matrix;
}


void CellGrid::setCellGridRegion(
	const std::vector<Real3>& atom_positions, const std::vector<int>& target_atoms)
{
	// Setup 
	// - Assumes orthorhombic box
	Real3 box_lengths = simulation_box_.getLengths();
	for ( int d=0; d<DIM_; ++d ) {
		bounding_box_offset_[d]    = std::numeric_limits<double>::max();
		bounding_box_matrix_[d][d] = std::numeric_limits<double>::lowest();

		// Set off-diagonal elements to zero
		for ( int e=0; e<DIM_; ++e ) {
			if ( d != e ) { bounding_box_matrix_[d][e] = 0.0; }
		}
	}

	// Set bounding box
	int num_target_atoms = target_atoms.size();
	if ( num_target_atoms > 0 ) {
		// Search positions to put bounds on region
		for ( int u=0; u<num_target_atoms; ++u ) {
			int index = target_atoms[u];
			// Compare to current bounds
			for ( int d=0; d<DIM_; ++d ) {
				// Lower bound
				if ( atom_positions[index][d] < bounding_box_offset_[d] ) {
					bounding_box_offset_[d] = atom_positions[index][d];
				}
				// Upper bound
				if ( atom_positions[index][d] > bounding_box_matrix_[d][d] ) {
					bounding_box_matrix_[d][d] = atom_positions[index][d];
				}
			}
		}
	}
	else {
		// No target atoms: default is entire box
		for ( int d=0; d<DIM_; ++d ) {
			bounding_box_offset_[d]    = 0.0;
			bounding_box_matrix_[d][d] = box_lengths[d];
		}
	}
}


void CellGrid::setTemplateCell(const double cell_side_length)
{
	// Assumes orthorhombic box)
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( a == b ) {
				template_cell_box_matrix_[a][a] = cell_side_length;
			}
			else {
				template_cell_box_matrix_[a][b] = 0.0;
			}
		}
	}
}


void CellGrid::setTemplateCellUsingGridDimensions() 
{
	// Check grid dimensions
	for ( int d=0; d<DIM_; ++d ) {
		if ( grid_dimensions_[d] < 1 ) {
			std::stringstream err_ss;
			err_ss << "Error in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ":" << __LINE__ << ")\n"
						 << "  Number of cells along each axis must be positive.\n";

			err_ss << "  (Grid dimensiuons are currently {";
			for ( int a=0; a<DIM_; ++a ) {
				err_ss << " " << grid_dimensions_[a];
			}
			err_ss << " }.)\n";

			throw std::runtime_error( err_ss.str() );
		}
	}

	// Assumes orthorhombic box
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( a == b ) {
				// Main diagonal
				template_cell_box_matrix_[a][a] = (bounding_box_matrix_[a][a] - bounding_box_offset_[a])/
				                                  static_cast<double>(grid_dimensions_[a]);
			}
			else {
				// Off-diagonal entries are zero
				template_cell_box_matrix_[a][b] = 0.0;
			}
		}
	}
}


void CellGrid::setGridDimensionsUsingTemplateCell()
{
	// Assumes orthorhombic box
	double dx;
	for ( int d=0; d<DIM_; ++d ) {
		dx = bounding_box_matrix_[d][d] - bounding_box_offset_[d];
		// Require min > max: box cannot cross periodic boundaries or have zero volume
		if ( dx < 0.0 ) {
			std::stringstream err_ss;
			err_ss << "CellGrid:setGridDimensionsUsingTemplateCell \n"
			       << "  Invalid bounds along dimension " << d+1 << " of 3 (min < max).\n"
						 << "  (input: min=" << bounding_box_matrix_[d][d] << ", max=" << bounding_box_offset_[d] << ")\n";
			throw std::runtime_error( err_ss.str() );
		}

		// Round down so that all cells are at least as large as the template
		// (cells at the edges may be larger to ensure that the grid covers
		// the entire region of interest)
		grid_dimensions_[d] = static_cast<int>( std::floor(dx/template_cell_box_matrix_[d][d]) );
		if ( grid_dimensions_[d] < 1 ) {
			grid_dimensions_[d] = 1;
		}
	}
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
	// Currently this code assumes orthorhombic boxes and subvolumes
	if ( not is_box_valid(bounding_box_matrix_) ) {
		std::stringstream err_ss;
		err_ss << "Error in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ":" << __LINE__ << ")\n"
		       << "  Bounding box for region of interest is invalid.\n"
		       << "\n"
		       << "bounding_box_matrix = [";
		
		for ( int a=0; a<DIM_; ++a ) {
			err_ss << " [";
			for ( int b=0; b<DIM_; ++b ) { 
				err_ss << " " << bounding_box_matrix_[a][b];
			}
			err_ss << " ]\n";
		}
		err_ss << " ]\n";

		throw std::runtime_error( err_ss.str() );
	}

	// Setup
	if ( partition_type_ == PartitionType::LocalNeighborList ) {
		// Template cell already made: need to set grid dimensions
		// - Assumes the bounding box has been set
		setGridDimensionsUsingTemplateCell();
	}
	else if ( partition_type_ == PartitionType::SimulationBoxDomainDecomposition ) {
		// Grid dimensions are fixed: need to make template cell
		// - Assumes that the bounding box has been set
		setTemplateCellUsingGridDimensions();
	}
	else {
		throw std::runtime_error("partition type not recognized");
	}

	// Determine the number of cells
	int num_cells = 1;
	for ( int d=0; d<DIM_; ++d ) {
		num_cells *= grid_dimensions_[d];

		// TODO move check to function
		if ( grid_dimensions_[d] < 1 ) {
			std::stringstream err_ss;
			err_ss << "Error in CellGrid::makeCells() [";
			if ( partition_type_ == PartitionType::LocalNeighborList ) {
				err_ss << "LocalNeighborList";
			}
			else if ( partition_type_ == PartitionType::SimulationBoxDomainDecomposition ) {
				err_ss << "SimulationBoxDomainDecomposition";
			}
			err_ss << " mode]\n";

			err_ss << "  Number of cells along " << CommonTypes::axis_names[d] << " is " 
			       << grid_dimensions_[d] << ".\n";
			throw std::runtime_error( err_ss.str() ); 
		}
	}

	// Allocate memory
	cells_.resize(num_cells);
	for ( int k=0; k<num_cells; ++k ) {
		cells_[k].neighbors.resize(0);
		cells_[k].neighbors.reserve(MAX_NUM_NEIGHBORING_CELLS);
	}

	// Upper and lower edges of the cells along each dimension
	// - Compute them this way so that the resulting partition leaves no
	//   gaps due to roundoff error
	for ( int d=0; d<DIM_; ++d ) {
		cell_lower_bounds_[d].resize( grid_dimensions_[d] );
		cell_upper_bounds_[d].resize( grid_dimensions_[d] );

		// First cell is flush with the bounding box offset
		cell_lower_bounds_[d][0] = bounding_box_offset_[d];

		for ( int i=1; i<grid_dimensions_[d]; ++i ) {
			cell_upper_bounds_[d][i-1] = cell_lower_bounds_[d][i-1] + 
			                             template_cell_box_matrix_[d][d];

			cell_lower_bounds_[d][i] = cell_lower_bounds_[d][i-1] + 
			                           template_cell_box_matrix_[d][d];
		}

		// Uppermost bound extends to the boundary of the gridded region
		cell_upper_bounds_[d].back() = bounding_box_matrix_[d][d];
	}

	// Create each cell
	int   cell_index, neighbor_index;   // indices in flattened cell array ("linear" indices)
	Int3  neighbor_indices, offset, grid_indices;
	for ( int l=0; l<grid_dimensions_[X_DIM]; ++l ) { // cells along x
		for ( int m=0; m<grid_dimensions_[Y_DIM]; ++m ) { // cells along y
			for ( int n=0; n<grid_dimensions_[Z_DIM]; ++n ) { // cells along z
				// Indices of this cell in the grid
				grid_indices = {{ l, m, n }};
				cell_index = get_cell_index( grid_indices );  // linear index in cells_
				cells_[cell_index].indices = grid_indices;

				// Key points in cell (assumes orthorhombic box)
				for ( int d=0; d<DIM_; ++d ) {
					// Lower-left corner
					cells_[cell_index].x_lower[d] = cell_lower_bounds_[d][ grid_indices[d] ];

					// Upper-right corner
					cells_[cell_index].x_upper[d] = cell_upper_bounds_[d][ grid_indices[d] ];

					// Center
					cells_[cell_index].x_center[d] = 0.5*( cells_[cell_index].x_upper[d] - 
					                                       cells_[cell_index].x_lower[d] );
				}

				// Cell box matrix (assumes orthorhombic box)
				for ( int a=0; a<DIM_; ++a ) {
					for ( int b=0; b<DIM_; ++b ) {
						if ( a == b ) {
							cells_[cell_index].box_matrix[a][a] = cells_[cell_index].x_upper[a] - 
							                                      cells_[cell_index].x_lower[a];
						}
						else {
							cells_[cell_index].box_matrix[a][b] = 0.0;
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
							// - TODO for grids that aren't close to the periodic boundaries of
							//   the box, this will lead to unnecessary checking between cells
							//   at opposite edges of the grid. Maybe check the distance between 
							//   cell centers to eliminate this?
							for ( int d=0; d<DIM_; ++d ) {
								neighbor_indices[d] = cells_[cell_index].indices[d] + offset[d];
								if ( neighbor_indices[d] >= grid_dimensions_[d] ) {
									neighbor_indices[d] = 0;
								}
								else if ( neighbor_indices[d] < 0 ) {
									neighbor_indices[d] = grid_dimensions_[d] - 1;
								}
							}

							// Corresponding linear index
							neighbor_index = get_cell_index( neighbor_indices );

							// Add the following check to handle cases where a=b=c=0,
							// and when the number of cells along a given direction is less than 3
							if ( neighbor_index != cell_index ) {
								cells_[cell_index].neighbors.push_back( neighbor_index );
							}
						}
					}
				}

				// Sort indices (for convenience)
				std::sort( cells_[cell_index].neighbors.begin(), cells_[cell_index].neighbors.end() );
			}
		}
	}

	// If there are fewer than 3 cells along any direction, the algorithm above will
	// produce a list which includes duplicates: remove them
	bool duplicates_possible = false;
	for ( int d=0; d<DIM_; ++d ) {
		if ( grid_dimensions_[d] < 3 ) {
			duplicates_possible = true;
			break;
		}
	}
	if ( duplicates_possible ) {
		for ( int k=0; k<num_cells; ++k ) {
			cells_[k].neighbors.erase( 
					unique( cells_[k].neighbors.begin(), cells_[k].neighbors.end() ),
					cells_[k].neighbors.end() 
			);
		}
	}

	// DEBUG
	/*
	if ( partition_type_ == PartitionType::SimulationBoxDomainDecomposition ) {
		std::cout << "Simulation box:  L =";
		Real3 box_lengths = simulation_box_.getLengths();
		for ( int d=0; d<DIM_; ++d ) {
			std::cout << " " << box_lengths[d];
		}
		std::cout << "\n";

		int num_cells = cells_.size();
		for ( int k=0; k<num_cells; ++k ) {
			std::cout << "Cell " << k << "\n";

			std::cout << "  x_lower =";
			for ( int d=0; d<DIM_; ++d ) {
				std::cout << " " << cells_[k].x_lower[d];
			}
			std::cout << "\n";

			std::cout << "  x_upper =";
			for ( int d=0; d<DIM_; ++d ) {
				std::cout << " " << cells_[k].x_upper[d];
			}
			std::cout << "\n";

			std::cout << "  x_center =";
			for ( int d=0; d<DIM_; ++d ) {
				std::cout << " " << cells_[k].x_center[d];
			}
			std::cout << "\n";

			std::cout << "  box_matrix =\n";
			for ( int a=0; a<DIM_; ++a ) {
				std::cout << "    ";
				for ( int b=0; b<DIM_; ++b ) {
					std::cout << " " << cells_[k].box_matrix[a][b];
				}
				std::cout << "\n";
			}

			std::cout << "\n";
		}
	}
	*/

	return;
}


//
void CellGrid::sortTargetAtomsIntoCells(
	const std::vector<Real3>& coords, const std::vector<int>& target_atoms,
	std::vector<std::vector<int>>& cell_list, std::vector<int>& cell_indices)
{
	// Prepare output arrays
	int num_cells = cells_.size();
	if ( num_cells <= 0 ) {
		std::cerr << "CellGrid::sortTargetAtomsIntoCells - Function was called before "
		             << "cells were defined.\n" 
		          << "Use \"makeCells\" to specify and partition the target domain.\n";
	}
	cell_list.resize(num_cells);
	for ( int k=0; k<num_cells; ++k ) {
		cell_list[k].resize(0); // hopefully won't lead to a reallocation
	}
	int num_target_atoms = target_atoms.size();
	cell_indices.resize(num_target_atoms);

	// Sort atoms
	int    i, cell_index;
	Int3   grid_indices;
	Real3  x_offset = cells_[0].x_lower; // bottom-left corner of grid
	for ( int j=0; j<num_target_atoms; ++j ) {
		i = target_atoms[j]; // global index
		for ( int d=0; d<DIM_; ++d ) {
			// Most cells are boxes with size and shape given by "template_cell_box_matrix_"
			// - The ones at the upper edges along x, y, and z may be larger, but the following
			//   expression below will still work for these
			grid_indices[d] = static_cast<int>( std::floor( (coords[i][d] - x_offset[d])/template_cell_box_matrix_[d][d] ) );

			// Edge cases (the above expressions can suffer from roundoff error)
			if ( coords[i][d] < cell_lower_bounds_[d][ grid_indices[d] ] and
			     grid_indices[d] >= 1 
			) {
				grid_indices[d] -= 1;
			}
			else if ( coords[i][d] > cell_upper_bounds_[d][ grid_indices[d] ] and 
			          grid_indices[d] < grid_dimensions_[d] - 1 
			) {
				grid_indices[d] += 1;
			}

			if ( grid_indices[d] == grid_dimensions_[d] ) {
				// If the particle is at the upper bound along this direction, put it into the
				// cell right next to the edge
				grid_indices[d] = grid_dimensions_[d] - 1;
			}
		}

		// Store *target* atom index in appropriate cell
		cell_index = get_cell_index(grid_indices);
		if ( cell_index < 0 or cell_index > num_cells ) {
			// Invalid grid index
			std::stringstream err_ss;
			err_ss << "CellGrid::sortTargetAtomsIntoCells - could not place atom " << i+1
			             << " into a cell (position =";
			for ( int d=0; d<DIM_; ++d ) { err_ss << " " << coords[i][d]; }
			err_ss << ")\n"
			           << "You probably need to change the cell grid boundaries.\n";
			throw std::runtime_error( err_ss.str() );
		}
		cell_list[cell_index].push_back( j );
		cell_indices[j] = cell_index;

		// Check placement (TODO necessary check?)
		for ( int d=0; d<DIM_; ++d ) {
			if ( coords[i][d] < cells_[cell_index].x_lower[d] or
			     coords[i][d] > cells_[cell_index].x_upper[d] 
			) {
				std::stringstream err_ss;

				err_ss << "CellGrid::sortTargetAtomsIntoCells - "
				          << "Atom " << i+1 << " was placed in the wrong cell!\n";
				// Atom position
				err_ss << "Atom at {";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << coords[i][d]; }
				err_ss << " } nm\n";
				
				// Cell indices
				err_ss << "Cell " << cell_index << ", indices {";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[cell_index].indices[d]; }
				err_ss << " } of {";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << grid_dimensions_[d]-1; }
				err_ss << " }\n";

				// Cell: x_lower
				err_ss << "  Bottom left corner = {";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[cell_index].x_lower[d]; }
				err_ss << " } nm\n";

				// Cell: x_upper
				err_ss << "  Top right corner = {";
				for ( int d=0; d<DIM_; ++d ) { err_ss << " " << cells_[cell_index].x_upper[d]; }
				err_ss << " } nm\n";
	
				throw std::runtime_error( err_ss.str() );
			}
		}
	}

	return;
}


// Construct the neighbor list
void CellGrid::makeNeighborList(
	const std::vector<Real3>& coords, const std::vector<int>& target_atoms,
	// Output
	std::vector<std::vector<int>>& neighbor_list)
{
	// Sort atoms into cells based on partition defined by member variable cells_
	this->sortTargetAtomsIntoCells(coords, target_atoms, 
	                               cell_list_, cell_indices_);

	// Allocate memory
	int num_target_atoms = target_atoms.size();
	neighbor_list.resize(num_target_atoms);
	for ( int t=0; t<num_target_atoms; ++t ) {
		neighbor_list[t].resize(0);
	}

	// Perform neighbor search
	int i, j, k, u, cell_index, num_cells_to_search, num_atoms_in_cell;
	std::vector<int> cells_to_search;
	double dist2;
	Real3  x_ij;
	for ( int t=0; t<num_target_atoms; ++t ) {
		i = target_atoms[t];
		cell_index = cell_indices_[t];

		// Need to search this cell and all of its neighbors
		cells_to_search = cells_[cell_index].neighbors;
		cells_to_search.push_back(cell_index);
		num_cells_to_search = cells_to_search.size();
		for ( int l=0; l<num_cells_to_search; ++l ) {
			// Search for neighbors of particle 'i' in neighboring cell 'k'
			k = cells_to_search[l];
			num_atoms_in_cell = cell_list_[k].size();
			for ( int m=0; m<num_atoms_in_cell; ++m ) {
				u = cell_list_[k][m]; // cell list has target atom indices
				j = target_atoms[u];

				if ( j > i ) { // only need to consider each unique pair once
					simulation_box_.calculateDistance(coords[i], coords[j], x_ij, dist2);
					if ( dist2 < r_cutoff_sq_ ) {
						neighbor_list[t].push_back( u ); // neighbor of i
						neighbor_list[u].push_back( t ); // neighbor of j
					}
				}
			}
		}
	}
}


void CellGrid::findNearbyAtoms(
		const int cell_index, AtomGroup& atom_group
) const
{
	// Check input
	int num_cells = cells_.size();
	if ( num_cells < 1 ) {
		throw std::runtime_error("no cells are set: unable to find local atoms");
	} 
	else if ( cell_index >= num_cells or cell_index < 0 ) {
		throw std::runtime_error("bad cell index: unable to find local atoms");
	}

	//----- Create ProbeVolume for the Cell -----//

	// TODO Make a persistent member variable

	double sigma = 0.0, alpha_c = 0.0;  // No coarse-graining for cell
	bool   need_derivatives = false;

	ParameterPack probe_volume_input_pack("CellGridProbeVolume");

	ProbeVolume::ProbeVolumeInputPack input_pack = {
		probe_volume_input_pack,
		simulation_state_,
		mpi_communicator_,
		need_derivatives
	};

	const Cell& local_cell = cells_[cell_index];

	ProbeVolume_Box cell_probe_volume(input_pack, local_cell.x_lower, local_cell.box_matrix, 
	                                  sigma, alpha_c);
	cell_probe_volume.setShellParameters(width_shell_1_, width_shell_2_, alpha_c_shells_);


	//----- Search for Local Atoms -----//

	// Get handles to necessary AtomGroup member variables
	const std::vector<Real3>&          atom_positions = atom_group.get_atom_positions();
	std::vector<AtomGroup::NearbyAtom>& nearby_atoms    = atom_group.access_nearby_atoms();
	nearby_atoms.resize(0);

	// Working variables
	double h_v, dummy_htilde_v;
	Real3  dummy_derivatives, position;
	bool   is_in_local_cell, is_in_local_cell_shell_1, is_in_local_cell_shell_2;

	int num_atoms_in_group = atom_positions.size();
	for ( int i=0; i<num_atoms_in_group; ++i ) {
		// Ensure that the atom is inside the box
		// - In XTC files, atoms can be just outside the box
		position = atom_positions[i];
		simulation_box_.putInBox(position);

		// Check this cell
		is_in_local_cell = cell_probe_volume.isInProbeVolume(
					position, h_v, dummy_htilde_v, dummy_derivatives,
					is_in_local_cell_shell_1, is_in_local_cell_shell_2);

		if ( is_in_local_cell or is_in_local_cell_shell_1 or is_in_local_cell_shell_2 ) {
			// Check edge cases (assumes orthorhombic box)
			if ( is_in_local_cell ) {
				for ( int d=0; d<DIM_; ++d ) {
					if ( position[d] < local_cell.x_lower[d] or position[d] > local_cell.x_upper[d] ) {
						is_in_local_cell = false;
						is_in_local_cell_shell_1 = true;
					}
					else if ( position[d] == local_cell.x_upper[d] and grid_dimensions_[d] > 1 ) {
						// Cells don't own the atoms exactly at their upper edges 
						// (unless there's only one cell along that axis)
						is_in_local_cell = false;
						is_in_local_cell_shell_1 = true;
						break;
					}
				}
			}

			// Store the local atom's index and classification
			nearby_atoms.push_back( 
				AtomGroup::NearbyAtom{ i, is_in_local_cell, 
				                      is_in_local_cell_shell_1, is_in_local_cell_shell_2 }
			);
		}
	}

	
	//----- DEBUG -----//

	/*
	int num_nearby_atoms = nearby_atoms.size();
	std::cout << "\n" << cell_probe_volume.getInputSummary() << "\n"
	          << "  num_nearby_atoms = " << num_nearby_atoms << "\n";

	// Copy local atom indices
	std::vector<int> local_atom_indices(num_nearby_atoms);	
	for ( int j=0; j<num_nearby_atoms; ++j ) {
		local_atom_indices.push_back( nearby_atoms[j].group_index );
	}

	// Find target atoms which are *not* local
	const std::vector<int>& global_atom_indices = atom_group.get_global_atom_indices();
	std::cout << "Non-local target atoms\n";
 	for ( int i=0; i<num_atoms; ++i ) {
		
		auto local_atom_it = std::find( local_atom_indices.begin(), local_atom_indices.end(), i );
		if ( local_atom_it == local_atom_indices.end() ) {
			// If std::find fails, atom with group index i is not a local target atom
			std::cout << " global " << global_atom_indices[i]+1 << " (local " << i+1 << "):  x = {";

			// Position
			for ( int d=0; d<DIM_; ++d ) { std::cout << " " << atom_positions[i][d]; };
			std::cout << " },  ";

			// Minimum image distance to center of box (simulation box center = cell center for a single rank)
			double dist_sq, dist;
			Real3 x_center_i;
			simulation_box_.calculateDistance( local_cell.x_center, atom_positions[i], x_center_i, dist_sq );
			dist = sqrt(dist_sq);

			std::cout << "x_center_i = {";
			for ( int d=0; d<DIM_; ++d ) { std::cout << " " << x_center_i[d]; };
			std::cout << " },  ";
				
			std::cout << "dist = " << dist;

			std::cout << "\n";  // end the line
		}
	}
	*/
}

std::ostream& operator<<(std::ostream& os, const CellGrid& cell_grid) 
{
	cell_grid.printCells(os);
	return os;
}


void CellGrid::printCells(std::ostream& os) const
{
	// Print simulation box for reference (assumes orthorhombic box)
	os << "Simulation box:  L =";
	Real3 box_lengths = simulation_state_.get_simulation_box().getLengths();
	for ( int d=0; d<DIM_; ++d ) {
		os << " " << box_lengths[d];
	}
	os << "\n";

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


void CellGrid::printAtomGroupDomainDecompositionCheck(std::ostream& os, const AtomGroup& atom_group)
{
	if ( mpi_communicator_.is_mpi_initialized() and 
			partition_type_ == PartitionType::SimulationBoxDomainDecomposition 
	) {
		// Coordinates of atoms in the group
		const std::vector<Real3>& group_atom_positions = atom_group.get_atom_positions();

		// Global indices (left here for reference)
		//const std::vector<int>& group_global_indices = atom_group.get_global_atom_indices();

		// Group atoms which are local to this rank
		const std::vector<AtomGroup::NearbyAtom>& nearby_atoms = atom_group.get_nearby_atoms();
		int num_nearby_atoms = nearby_atoms.size();

		// MPI ranks
		int my_rank     = mpi_communicator_.get_rank();
		int master_rank = mpi_communicator_.get_master_rank();
		int num_ranks   = mpi_communicator_.get_size();

		// Find the group indices of the atoms owned by this rank
		std::vector<int> my_owned_atoms;
		int num_atoms_in_local_dd_cell = 0;
		for ( int k=0; k<num_nearby_atoms; ++k ) {
			if ( nearby_atoms[k].is_in_local_cell ) {
				// Atom belongs to this rank's DD cell
				++num_atoms_in_local_dd_cell;
				my_owned_atoms.push_back( nearby_atoms[k].group_index );
			}
		}

		// Header
		if ( my_rank == master_rank ) {
			os << "\n(master rank) DomainDecomposition check\n" << std::flush;
		}
		mpi_communicator_.barrier();

		// Print number of atoms "owned" by each processor
		for ( int r=0; r<num_ranks; ++r ) {
			if ( r == my_rank ) {
				os << "(rank " << my_rank << ") " << num_atoms_in_local_dd_cell 
					 << " atoms in local dd cell\n" << std::flush;
			}
			mpi_communicator_.barrier();
		}

		// Sum of local DD cell atoms over all ranks should be the same as the total
		// number of atoms in the group
		// - Each atom in the group is assigned a unique cell that "owns" it
		int num_atoms          = atom_group.get_num_atoms();
		int num_dd_atoms_total = num_atoms_in_local_dd_cell;
		mpi_communicator_.allreduce_sum_in_place(num_dd_atoms_total);
		if ( my_rank == master_rank and num_dd_atoms_total != num_atoms ) {
			os << "num_dd_atoms_total = " << num_dd_atoms_total 
			   << "  (should be " << num_atoms << ")\n" << std::flush;
		}
		mpi_communicator_.barrier();

		// Collect local atom group indices on master
		std::vector<std::vector<int>> owned_atoms_lists(num_ranks);
		if ( my_rank == master_rank ) {

			for ( int r=0; r<num_ranks; ++r ) {
				// Reserve space for recv
				owned_atoms_lists[r].resize( num_atoms );

				if ( r == master_rank ) {
					// Master rank's owned target atoms
					owned_atoms_lists[master_rank] = my_owned_atoms;
				}
				else {
					// Receive owned target atoms from rank 'r'
					int tag = r;
					MpiCommunicator::MpiStatus status;
					mpi_communicator_.recv( owned_atoms_lists[r], r, tag, status );
				}
				
				//os << "MASTER: rank " << r << " has " << owned_atoms_lists[r].size()
				//					<< " owned target atoms\n";
			}
		}
		else {
			// Workers send owned target atoms to the master
			int tag = my_rank;
			mpi_communicator_.send( my_owned_atoms, master_rank, tag );
		}
		mpi_communicator_.barrier();

		/*
		// Check numbers received by master vs. those sent by workers
		int max_i = 5;
		for ( int r=0; r<num_ranks; ++r ) {
			// Locally owned
			if ( my_rank == r ) {
				os << "(rank " << r << ") group indices =";
				for ( int i=0; i<max_i; ++i ) {
					os << " " << my_owned_atoms[i];
				}
				os << " ... \n";
			}
			mpi_communicator_.barrier();

			// Copied to masster
			if ( my_rank == master_rank ) {
				os << "(master) rank " << r << " group indices =";
				for ( int i=0; i<max_i; ++i ) {
					os << " " << owned_atoms_lists[r][i];
				}
				os << " ... \n";
			}
			mpi_communicator_.barrier();
		}
		mpi_communicator_.barrier();
		*/

		// Master rank searches for duplicates
		if ( my_rank == master_rank ) {
			for ( int i=0; i<num_atoms; ++i ) {
				// Find the ranks with this atom
				std::vector<int> ranks_with_atom;
				for ( int r=0; r<num_ranks; ++r ) {
					auto vec_begin = owned_atoms_lists[r].begin();
					auto vec_end   = owned_atoms_lists[r].end();
					auto it = std::find(vec_begin, vec_end, i);
					if ( it != vec_end ) {
						// Found atom
						ranks_with_atom.push_back(r);
					}
				}

				// Print atoms that are present on more than 1 rank (or none!)
				int num_ranks_with_atom = ranks_with_atom.size();
				if ( num_ranks_with_atom != 1 ) {

					// Serial number and position
					os << "DUPLICATE: atom " << i+1 << ": x = {";
					for ( int d=0; d<DIM_; ++d ) {
						os << " " << group_atom_positions[i][d];
					}
					os << " },  ranks =";

					// Print the ranks
					for ( int j=0; j<num_ranks_with_atom; ++j ) {
						os << " " << ranks_with_atom[j];
					}
					os << "\n" << std::flush;
				}
			}
		}
		mpi_communicator_.barrier();
	}
	else {
		os << "DomainDecomposition check: nothing to print\n";
		if ( not mpi_communicator_.is_mpi_initialized() ) {
			os << "  MPI is not initialized.\n";
		}
		if ( partition_type_ != PartitionType::SimulationBoxDomainDecomposition ) {
			os << "  This CellGrid is not being used for DomainDecomposition.\n";
		}
		os << std::flush;
	}
}
