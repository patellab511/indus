/* Indus.cpp
 */

#include "Indus.h"

Indus::Indus(
	std::shared_ptr<ProbeVolume>& probe_volume_ptr, ProbeVolumeInputPack& input_pack,
	OrderParameter::InputPack& op_input_pack
):
	OrderParameter(op_input_pack),
	probe_volume_ptr_(probe_volume_ptr)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Target atoms
	std::string target_atom_group_name("target_atoms");
	input_parameter_pack.readString("target_atom_group", KeyType::Optional, target_atom_group_name);
	target_group_ = this->addAtomGroupDependency( target_atom_group_name );
}


void Indus::calculate()
{
	calculate_timer_.start();

	if ( simulation_state_.debug_mode() ) {
		std::cout << "PROBE VOLUME -> calculate\n";
	}
	
	// Reset INDUS variables
	n_v_ = 0;
	value_ = 0.0;
	local_group_indices_[target_group_].resize(0);

	htilde_v_.resize(0);
	local_derivatives_[target_group_].resize(0);

	nonlocal_atom_group_indices_.resize(0);
	nearby_shell_1_atom_group_indices_.resize(0);
	nearby_shell_2_atom_group_indices_.resize(0);

	// Group information
	// FIXME This is a bit cludgy to write out each time
	const AtomGroup& target_atom_group = 
		simulation_state_.get_atom_group( atom_group_indices_[target_group_] );
	const std::vector<Real3>& group_atom_positions = target_atom_group.get_atom_positions();

	// Group atoms which are local to this rank
	const std::vector<AtomGroup::NearbyAtom>& nearby_target_atoms = target_atom_group.get_nearby_atoms();
	const int num_nearby_target_atoms = nearby_target_atoms.size();

	#pragma omp parallel
	{
		// Determine the location of each target atom relative to the probe volume
		// - Also compute the regular and coarse-grained number of atoms in the probe volume
		double     h_v, htilde_v;
		RegionEnum region;
		int        group_index;
		Real3      derivs_htilde_v; // derivatives of probe volume htilde_v for i, dh^[v]/dr

		// Some atoms will be quickly identified as far away from the probe volume, so use a
		// dynamic schedule to keep the threads busy
		#pragma omp for schedule(dynamic)
		for ( int k=0; k<num_nearby_target_atoms; ++k ) {
			calculate_omp_timer_.start();

			// Atom's index in the group 
			group_index = nearby_target_atoms[k].group_index;

			// Result is "true" if the target is in the *non*-coarse-grained probe volume
			// - TODO: Change signature for probe volumes with "probe volume atoms"
			probe_volume_ptr_->calculateIndicator(
					group_atom_positions[group_index], h_v, htilde_v, derivs_htilde_v, region );

			if ( region != RegionEnum::Unimportant ) {
				// Save this atom
				// TODO: Refactor to avoid a critical section?
				calculate_omp_critical_timer_.start();
				#pragma omp critical
				{
					if ( region == RegionEnum::Vtilde ) {
						// Atom is in the coarse-grained probe volume
						if ( nearby_target_atoms[k].is_in_local_cell ) {
							// Atom belongs to this rank's cell
							value_ += htilde_v;
							if ( h_v == 1.0 ) {
								// Atom is in the unsmoothed probe volume
								++n_v_;
							}

							// Save atom indices
							local_group_indices_[target_group_].push_back( group_index );

							// Store information about indicator function (aka probe volume htilde_v)
							htilde_v_.push_back( htilde_v );
							if ( need_derivatives_ ) {
								local_derivatives_[target_group_].push_back( derivs_htilde_v ); 
							}
						}
						else {
							// Atom is in vtilde, but doesn't belong to this rank
							nonlocal_atom_group_indices_.push_back( group_index );
						}
					}
					else if ( region == RegionEnum::Shell_1 ) {
						nearby_shell_1_atom_group_indices_.push_back( group_index );
					}
					else if ( region == RegionEnum::Shell_2 ) {
						nearby_shell_2_atom_group_indices_.push_back( group_index );
					}
				} // end pragma omp critical
				calculate_omp_critical_timer_.stop();
			} // end if not unimportant
			calculate_omp_timer_.stop();
		} // end loop over atoms
	} // end pragma omp parallel

	if ( simulation_state_.debug_mode() ) {
		std::cout << "PROBE VOLUME -> calculate DONE\n";
	}

	calculate_timer_.stop();

	return;
}


void Indus::synchronize()
{
	// MPI: Sum contributions to N_v and Ntilde_v over all ranks
	if ( mpi_communicator_.is_mpi_initialized() ) {
		allreduce_timer_.start();
		mpi_communicator_.allreduce_sum_in_place(n_v_);
		mpi_communicator_.allreduce_sum_in_place(value_);
		allreduce_timer_.stop();
	}

	if ( need_derivatives_ ) {
		// Calculate matrix required for virial
		this->calculate_sum_r_cross_derivatives();
	}
}
