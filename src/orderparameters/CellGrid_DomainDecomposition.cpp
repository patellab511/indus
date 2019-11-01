#include "CellGrid_DomainDecomposition.h"

CellGrid_DomainDecomposition::CellGrid_DomainDecomposition(
	const Int3& grid_dimensions, const double width_shell_1, const double width_shell_2,
	const SimulationState& simulation_state, MpiCommunicator& mpi_communicator
):
	CellGrid(simulation_state, mpi_communicator),
	width_shell_1_( width_shell_1 ),
	width_shell_2_( width_shell_2 )
{
	grid_dimensions_ = grid_dimensions;

	// Default: grid the entire simulation box
	setCellGridRegion();


	//-----  Create a probe box to represent the local DD cell -----//

	// No coarse-graining or derivatives required
	double sigma = 0.0, alpha_c = 0.0;
	bool   need_derivatives = false;

	// Probe volume base class requires a ParameterPack
	ParameterPack probe_volume_input_pack("CellGridProbeVolume");

	ProbeVolume::ProbeVolumeInputPack input_pack = {
		probe_volume_input_pack,
		simulation_state_,
		mpi_communicator_,
		need_derivatives
	};

	if ( simulation_state.debug_mode() ) {
		std::cout << "MAKING CELLGRID PROBE VOLUME [" << LOCATION_IN_SOURCE_STRING << "]\n";
	}

	// Initialize with dummy geometry that will be overwritten later
	Real3 x_lower_tmp;  x_lower_tmp.fill(0.0);
	const Box& box_matrix_tmp = simulation_box_.get_box_matrix();
	local_cell_probe_volume_ptr_ = std::unique_ptr<ProbeVolume_Box>(
		new ProbeVolume_Box{ input_pack, x_lower_tmp, box_matrix_tmp, sigma, alpha_c }
	);
	setShellWidths(width_shell_1_, width_shell_2_);

	if ( simulation_state.debug_mode() ) {
		std::cout << "DONE MAKING CELLGRID PROBE VOLUME [" << LOCATION_IN_SOURCE_STRING << "]\n";
	}

	this->makeCells();
}


void CellGrid_DomainDecomposition::updateGrid() 
{
	// Check grid dimensions
	for ( int d=0; d<DIM_; ++d ) {
		if ( grid_dimensions_[d] < 1 ) {
			std::stringstream err_ss;
			err_ss << "Error in " << FANCY_FUNCTION << "\n"
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
	const Box& bounding_box_matrix = bounding_box_.get_box_matrix();
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( a == b ) {
				// Main diagonal
				template_cell_box_matrix_[a][a] = bounding_box_matrix[a][a]/
				                                    static_cast<double>(grid_dimensions_[a]);
			}
			else {
				// Off-diagonal entries are zero
				template_cell_box_matrix_[a][b] = 0.0;
			}
		}
	}
}


void CellGrid_DomainDecomposition::findNearbyAtoms(
		const int cell_index, AtomGroup& atom_group
) const
{
	find_nearby_atoms_timer_.start();

	// Check input
	const int num_cells = cells_.size();
	if ( num_cells < 1 ) {
		throw std::runtime_error("no cells are set: unable to find local atoms");
	} 
	else if ( cell_index >= num_cells or cell_index < 0 ) {
		throw std::runtime_error("bad cell index: unable to find local atoms");
	}

	// Update local cell
	const Cell& local_cell = cells_[cell_index];
	local_cell_probe_volume_ptr_->setGeometry(local_cell.x_lower, local_cell.x_upper);
	local_cell_probe_volume_ptr_->update();

	//----- Search for Local Atoms -----//

	// Get handles to necessary AtomGroup member variables
	const std::vector<Real3>&           atom_positions = atom_group.get_atom_positions();
	std::vector<AtomGroup::NearbyAtom>& nearby_atoms   = atom_group.access_nearby_atoms();
	nearby_atoms.resize(0);
	const int num_atoms_in_group = atom_positions.size();

	// Set up buffers
	// TODO OPENMP Flexible number of threads? Move under 'omp parallel' with '#pragma omp master'?
	const int num_threads = OpenMP::get_max_threads();
	nearby_atom_buffers_.resize(num_threads);

	#pragma omp parallel num_threads(num_threads)
	{
		// Working variables
		double h_v, dummy_htilde_v;
		Real3  dummy_derivatives;
		bool   is_in_local_cell;  // is_in_local_cell_shell_1, is_in_local_cell_shell_2;
		bool   is_in_local_cell_check;  //is_in_simulation_box
		RegionEnum region;

		// Use the appropriate buffer for this thread
		const int thread_id = OpenMP::get_thread_num();
		std::vector<AtomGroup::NearbyAtom>* my_nearby_atoms_ptr = nullptr;
		if ( thread_id == 0 ) {
			my_nearby_atoms_ptr = &nearby_atoms;  // put result directly into output vector
		}
		else {
			my_nearby_atoms_ptr = &( nearby_atom_buffers_[thread_id] );
		}
		std::vector<AtomGroup::NearbyAtom>& my_nearby_atoms = *my_nearby_atoms_ptr;
		my_nearby_atoms.resize(0);

		#pragma omp for
		for ( int i=0; i<num_atoms_in_group; ++i ) {
			dd_omp_timer_.start();

			// Check whether the image available is inside the simulation box
			// - In XTC files and some DD schemes, atoms can be just outside the box
			// - One could try and shift the atom into the box, but roundoff error from the shift
			//   might cause the atom to be lost
			const Real3& position = atom_positions[i];
			//is_in_simulation_box = simulation_box_.is_in_box(position);
			//Real3 position = simulation_box_.putInBox(atom_positions[i]);

			// Check this cell
			local_cell_probe_volume_ptr_->calculateIndicator(
					position, h_v, dummy_htilde_v, dummy_derivatives, region );
			is_in_local_cell = ( h_v == 1.0 );

			// Manually check to catch edge cases when the atom is near the edges of the local cell
			is_in_local_cell_check = true;
			for ( int d=0; d<DIM_; ++d ) {
				if ( position[d] < local_cell.x_lower[d] or
						 position[d] > local_cell.x_upper[d] 
				) {
					is_in_local_cell_check = false;
					break;
				}
			}
			if ( is_in_local_cell_check ) {
				// Check one last edge case (assumes orthorhombic box)
				for ( int d=0; d<DIM_; ++d ) {
					if ( position[d] == local_cell.x_upper[d] and grid_dimensions_[d] > 1 ) {
						// Cells don't own the atoms exactly at their upper edges 
						// (unless there's only one cell along that axis)
						is_in_local_cell_check = false;
						break;
					}
				}
			}

			if ( is_in_local_cell_check ) {
				// Definitely in local cell, and not shells
				is_in_local_cell = true;
				region = RegionEnum::Vtilde;
			}
			else if ( is_in_local_cell ) {  // implicit: "and (not is_in_local_cell_check)"
				// Particle Was mistakenly placed placed in the local cell; it's actually in the shells
				is_in_local_cell = false;
				if ( width_shell_1_ > 0.0 ) {
					region = RegionEnum::Shell_1;
				}
				else if ( width_shell_2_ > 0.0 ) {
					region = RegionEnum::Shell_2;
				}
				else {
					region = RegionEnum::Unimportant;
				}

				/*
				// DEBUG
				if ( region != RegionEnum::Unimportant ) {
					std::cout << "(rank " << mpi_communicator_.get_rank() << ") MISPLACED BUT CAUGHT: x =";
					for ( int d=0; d<DIM_; ++d ) {
						std::cout << "  " << position[d];
					}
					std::cout << std::endl;
				}
				*/
			}

			if ( region != RegionEnum::Unimportant ) {
				// Store the local atom's index and classification
				my_nearby_atoms.push_back( 
					AtomGroup::NearbyAtom{ i, is_in_local_cell, region }
				);
			}

			dd_omp_timer_.stop();
		} // end loop over atoms in the group
	} // end pragma omp parallel

	if ( num_threads > 1 ) {
		// Merge the buffers
		dd_omp_combine_timer_.start();
		for ( int t=1; t<num_threads; ++t ) {
			nearby_atoms.insert( nearby_atoms.end(), nearby_atom_buffers_[t].begin(), 
			                     nearby_atom_buffers_[t].end() );
		}
		dd_omp_combine_timer_.stop();
	}

	find_nearby_atoms_timer_.stop();

	//----- DEBUG -----//

	/*
	int my_rank     = mpi_communicator_.get_rank();
	int master_rank = mpi_communicator_.get_master_rank();
	int num_ranks   = mpi_communicator_.get_size();

	if ( my_rank == master_rank ) {
		for ( int r=0; r<num_ranks; ++r ) {
			std::cout << "(rank " << r << ")\n";
			std::cout << "  x_lower =";
			for ( int d=0; d<DIM_; ++d ) {
				std::cout << "  " << cells_[r].x_lower[d];
			}
			std::cout << std::endl;
			std::cout << "  x_upper =";
			for ( int d=0; d<DIM_; ++d ) {
				std::cout << "  " << cells_[r].x_upper[d];
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl << std::flush;
	mpi_communicator_.barrier();

	printAtomGroupDomainDecompositionCheck(std::cout, atom_group);

	int num_nearby_atoms = nearby_atoms.size();
	std::cout << "\n" << local_cell_probe_volume_ptr_->getInputSummary("") << "\n"
	          << "  num_nearby_atoms = " << num_nearby_atoms << "\n";

	// Copy local atom indices
	std::vector<int> local_atom_indices(num_nearby_atoms);	
	for ( int j=0; j<num_nearby_atoms; ++j ) {
		local_atom_indices.push_back( nearby_atoms[j].group_index );
	}

	// Find target atoms which are *not* local
	const std::vector<int>& global_atom_indices = atom_group.get_global_atom_indices();
	std::cout << "Non-local target atoms\n";
 	for ( int i=0; i<num_nearby_atoms; ++i ) {
		
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


void CellGrid_DomainDecomposition::printAtomGroupDomainDecompositionCheck(
	std::ostream& os, const AtomGroup& atom_group) const
{
	if ( mpi_communicator_.is_mpi_initialized() ) {
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
			os << "\n"
			   << "(master rank) DomainDecomposition check\n"
			   << "  t = " << simulation_state_.get_time() << "\n"
			   << "  AtomGroup = " << atom_group.get_name() << "\n"
			   << "  width_shell_1 = " << width_shell_1_ << "\n"
			   << "  width_shell_2 = " << width_shell_2_ << "\n"
			   << std::flush;
		}
		mpi_communicator_.barrier();

		// Print number of atoms "owned" by each processor
		for ( int r=0; r<num_ranks; ++r ) {
			// Share count so that messages are printed in the right order
			int count_r;
			if ( r != master_rank ) {
				if ( my_rank == r ) {
					// This worker sends its count
					int tag = r;
					mpi_communicator_.send( num_atoms_in_local_dd_cell, master_rank, tag );
				}
				else if ( my_rank == master_rank ) {
					// Master rank receives the count
					int tag = r;
					MpiCommunicator::MpiStatus status;
					mpi_communicator_.recv( count_r, r, tag, status );
				}
			}
			else {
				count_r = num_atoms_in_local_dd_cell;  // master already has its own count
			}

			// Master rank prints the count
			if ( my_rank == master_rank ) {
				os << "(rank " << r << ") " << count_r
					 << " atoms in local dd cell\n" << std::flush;
			}
		}
		mpi_communicator_.barrier();

		// Sum of local DD cell atoms over all ranks should be the same as the total
		// number of atoms in the group
		// - Each atom in the group is assigned a unique cell that "owns" it
		int num_atoms          = atom_group.get_num_atoms();
		int num_dd_atoms_total = num_atoms_in_local_dd_cell;
		mpi_communicator_.allreduce_sum_in_place(num_dd_atoms_total);
		if ( my_rank == master_rank and num_dd_atoms_total != num_atoms ) {
			os << "ERROR: num_dd_atoms_total = " << num_dd_atoms_total 
			   << "  (number in group: " << num_atoms << ")\n" << std::flush;
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
			}
		}
		else {
			// Workers send owned target atoms to the master
			int tag = my_rank;
			mpi_communicator_.send( my_owned_atoms, master_rank, tag );
		}
		mpi_communicator_.barrier();

		// Master rank searches for missing atoms and duplicates
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
				if ( num_ranks_with_atom > 1 ) {
					// Serial number and position
					os << "DUPLICATE: atom " << i+1 << ": x = {";
					for ( int d=0; d<DIM_; ++d ) {
						os << " " << group_atom_positions[i][d];
					}
					os << " },  ranks =";
					for ( int j=0; j<num_ranks_with_atom; ++j ) {
						os << " " << ranks_with_atom[j];
					}
					os << "\n" << std::flush;
				}
				else if ( num_ranks_with_atom == 0 ) {
					// Serial number and position
					os << "MISSING: atom " << i+1 << ": x = {";
					for ( int d=0; d<DIM_; ++d ) {
						os << " " << group_atom_positions[i][d];
					}
					os << " }\n" << std::flush;
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
		os << std::flush;
	}
}

