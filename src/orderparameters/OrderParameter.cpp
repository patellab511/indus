#include "OrderParameter.h"


OrderParameter::OrderParameter(InputPack& input_pack):
	// Access to state variables
	simulation_state_( input_pack.simulation_state ),
	simulation_box_( simulation_state_.get_simulation_box() ),
	atom_groups_( simulation_state_.get_atom_groups() ),
	// Domain decomposition (DD)
	domain_decomposition_( input_pack.domain_decomposition ),
	// MPI (utilizied with DD)
	mpi_communicator_( input_pack.mpi_communicator ),
	need_derivatives_( input_pack.need_derivatives )
{
}


OrderParameter::~OrderParameter() {}


int OrderParameter::addAtomGroupDependency(const std::string& group_name)
{
	// Search for desired group
	int index = -1;
	bool found_group = false;
	int num_groups = atom_groups_.size();
	for ( int k=0; k<num_groups; ++k ) {
		if ( atom_groups_[k].get_name() == group_name ) {
			found_group = true;
			atom_group_indices_.push_back(k);
			index = static_cast<int>(atom_group_indices_.size()) - 1;

			// Allocate memory
			int num_groups = atom_group_indices_.size();
			local_derivatives_.resize( num_groups );
			local_group_indices_.resize( num_groups );

			break;
		}
	}

	if ( not found_group ) {
		throw std::runtime_error("Atom group \"" + group_name + "\" does not exist.");
	}

	return index;
}


void OrderParameter::calculate_local_sum_r_cross_derivatives()
{
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			local_sum_r_cross_derivatives_[a][b] = 0.0;
		}
	}

	int num_groups_for_op = atom_group_indices_.size();
	for ( int k=0; k<num_groups_for_op; ++k ) {
		// Get positions of atoms in this group
		const std::vector<Real3> group_positions = 
				atom_groups_[atom_group_indices_[k]].get_atom_positions();

		// Loop over local atoms
		int num_local_atoms = local_derivatives_[k].size();
		for ( int i=0; i<num_local_atoms; ++i ) { 
			for ( int a=0; a<DIM_; ++a ) {
				for ( int b=0; b<DIM_; ++b ) {
					local_sum_r_cross_derivatives_[a][b] 
						+= group_positions[ local_group_indices_[k][i] ][a] * local_derivatives_[k][i][b];
				}
			}
		}
	}
}


void OrderParameter::calculate_sum_r_cross_derivatives()
{
	// Contribution from this rank
	this->calculate_local_sum_r_cross_derivatives();

	sum_r_cross_derivatives_ = local_sum_r_cross_derivatives_;
	if ( mpi_communicator_.is_mpi_initialized() and mpi_communicator_.get_size() > 1 ) {
		allreduce_timer_.start();
		mpi_communicator_.allreduce_sum_in_place( sum_r_cross_derivatives_ );
		allreduce_timer_.stop();
	}
}
