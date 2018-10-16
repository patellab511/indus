/* ProbeVolume.cpp
 */

#include "ProbeVolume.h"

ProbeVolume::ProbeVolume(ProbeVolumeInputPack& input_pack)
 : simulation_state_( input_pack.simulation_state ),
   simulation_box_( simulation_state_.get_simulation_box() ),
   mpi_communicator_( input_pack.mpi_communicator ),
   // Coarse-graining
   sigma_(0.01),
   alpha_c_(0.02),
   // Shells
   width_shell_1_(0.0),
   width_shell_2_(0.0),
   // Flags
   is_dynamic_(false),
   need_derivatives_( input_pack.need_derivatives )
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// FIXME delete?
	// Coarse-graining parameters
	//input_parameter_pack.readNumber("sigma",   KeyType::Optional, sigma_);
	//input_parameter_pack.readNumber("alpha_c", KeyType::Optional, alpha_c_);
	//switching_function_.setParameters(sigma_, alpha_c_);

	// Flags
	input_parameter_pack.readFlag("is_dynamic", KeyType::Optional, is_dynamic_);
}


void ProbeVolume::doIndusWithShells()
{
	// Reset INDUS variables
	n_v_ = 0;
	ntilde_v_ = 0.0;
	local_atom_group_indices_.resize(0);
	local_atom_global_indices_.resize(0);

	htilde_v_.resize(0);
	derivatives_htilde_v_.resize(0);

	nonlocal_atom_group_indices_.resize(0);
	nonlocal_atom_global_indices_.resize(0);

	nearby_shell_1_atom_group_indices_.resize(0);
	nearby_shell_1_atom_global_indices_.resize(0);

	nearby_shell_2_atom_group_indices_.resize(0);
	nearby_shell_2_atom_global_indices_.resize(0);


	// Update the probe volume (as needed for particular geometries)
	this->updateUsingSimulationState();

	// Determine the location of each target atom relative to the probe volume
	// - Also compute the regular and coarse-grained number of atoms in the probe volume
	double h_v, htilde_v;
	bool   is_in_probe_volume, is_in_shell_1, is_in_shell_2;
	int    group_index;
	Real3  derivs_htilde_v; // derivatives of probe volume htilde_v for i, dh^[v]/dr

	// Group information
	const AtomGroup& target_atom_group = simulation_state_.get_target_atom_group();
	const std::vector<Real3>& group_atom_positions = target_atom_group.get_atom_positions();
	const std::vector<int>&   group_global_indices = target_atom_group.get_global_atom_indices();

	// Group atoms which are local to this rank
	const std::vector<AtomGroup::NearbyAtom>& nearby_target_atoms = target_atom_group.get_nearby_atoms();
	int num_nearby_target_atoms = nearby_target_atoms.size();
	for ( int k=0; k<num_nearby_target_atoms; ++k ) {
		// Atom's index in the group 
		group_index = nearby_target_atoms[k].group_index;

		// Result is "true" if the target is in the *non*-coarse-grained probe volume
		is_in_probe_volume = this->isInProbeVolume(
				group_atom_positions[group_index],
				h_v, htilde_v, derivs_htilde_v, is_in_shell_1, is_in_shell_2
		);

		if ( htilde_v > 0.0 ) {
			// Atom is in the coarse-grained probe volume
			if ( nearby_target_atoms[k].is_in_local_cell ) {
				// Atom belongs to this rank's cell
				ntilde_v_ += htilde_v;
				if ( is_in_probe_volume ) {
					// Atom is in the unsmoothed probe volume
					++n_v_;
				}

				// Save atom indices
				local_atom_group_indices_.push_back( group_index );
				local_atom_global_indices_.push_back( group_global_indices[group_index] );

				// Store information about indicator function (aka probe volume htilde_v)
				htilde_v_.push_back( htilde_v );
				if ( need_derivatives_ ) {
					derivatives_htilde_v_.push_back( derivs_htilde_v ); 
				}
			}
			else {
				// Atom is in vtilde, but doesn't belong to this rank
				nonlocal_atom_group_indices_.push_back( group_index );
				nonlocal_atom_global_indices_.push_back( group_global_indices[group_index] );
			}
		}
		else if ( is_in_shell_1 ) {
			// Shell 1
			nearby_shell_1_atom_group_indices_.push_back( group_index );
			nearby_shell_1_atom_global_indices_.push_back( group_global_indices[group_index] );
		}
		else if ( is_in_shell_2 ) {
			// Shell 2
			nearby_shell_2_atom_group_indices_.push_back( group_index );
			nearby_shell_2_atom_global_indices_.push_back( group_global_indices[group_index] );
		}
	}

	// MPI: Sum contributions to N_v and Ntilde_v over all ranks
	if ( mpi_communicator_.is_mpi_initialized() ) {
		mpi_communicator_.allreduce_sum_in_place(n_v_);
		mpi_communicator_.allreduce_sum_in_place(ntilde_v_);
	}

	return;
}


void ProbeVolume::setBoundingBox()
{
	simulation_box_.getBoxMatrix(bounding_box_matrix_);

	// No offset from the origin
	for ( int d=0; d<DIM_; ++d ) {
		bounding_box_offset_[d] = 0.0;
	}
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume::getInputSummary() const
{
	return "# Probe_volume output_summary_not_yet_implemented\n" + 
	       this->getSharedAttributesSummary();
}


// Get summary of attributes shared by all probe volumes
// ( e.g. shell widths, coarse-graining parameters, etc.)
std::string ProbeVolume::getSharedAttributesSummary() const
{
	std::stringstream ss;
	ss << "#   Coarse-graining\n"
		 << "#     sigma            = " << sigma_   << " [nm]\n"
	   << "#     alpha_c          = " << alpha_c_ << " [nm]\n"
	   << "#   Shells\n"
	   << "#     width_shell_1  = " << width_shell_1_   << " [nm]\n"
	   << "#     width_shell_2  = " << width_shell_2_   << " [nm]\n";

	return ss.str();
}
