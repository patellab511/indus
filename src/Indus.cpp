#include "Indus.h"


Indus::Indus(const std::string& input_file)
 : input_file_(input_file),
   // Simulation state
   simulation_state_(),
   // Domain decomposition/MPI
   mpi_communicator_( MpiCommunicator::get_mpi_comm_world() ),
   domain_decomposition_( std::array<int,DIM_>({{1, 1, 1}}), 0.0, 0.0,
                          simulation_state_, mpi_communicator_ ),
   my_rank_( mpi_communicator_.get_rank() ),
   master_rank_( mpi_communicator_.get_master_rank() ),
   // Indus atoms
	 target_atoms_string_("atom_type OW"), // default: water oxygens
   // Biasing
	 need_derivatives_(false),
   derivatives_are_ready_(false),
   delta_x_( sqrt(std::numeric_limits<double>::epsilon()) ),
   // Output control
   want_output_(true),
   print_forces_(false), print_numerical_forces_(false),
   output_file_suffix_(""),
   t_min_( std::numeric_limits<double>::lowest() ),
   t_max_( std::numeric_limits<double>::max() )
{
#ifdef INDUS_PLUMED_MODE
	// Generally need derivatives if running with PLUMED
	// (almost always biasing something)
	need_derivatives_ = true;
#endif

	xtc_file_.clear();
	gro_file_.clear();
	top_file_.clear();

	bias_ntilde_v_string_.clear();
	bias_ntilde_v_ptr_.reset(nullptr);

	//----- Process input -----//

	// Read input file
	InputParser input_parser;
	input_parser.parseFile(input_file, input_parameter_pack_);

	using KeyType = ParameterPack::KeyType;
	std::string line, token, lowercase_token;
	std::string* value_ptr;

	// Target atoms
	const std::vector<std::string>& target_tokens = 
			input_parameter_pack_.findRequiredVector("Target");
	target_atoms_string_ = "";
	int num_tokens = target_tokens.size();
	for ( int i=0; i<num_tokens; ++i ) {
		if ( i > 0 ) { 
			target_atoms_string_ += " "; 
		}
		target_atoms_string_ += target_tokens[i];
	}

	// Probe volume
	const ParameterPack* probe_volume_pack_ptr = 
			input_parameter_pack_.findParameterPack("ProbeVolume", KeyType::Required);

	// Biases 
	std::vector<const ParameterPack*> bias_pack_ptrs = 
			input_parameter_pack_.findParameterPacks("Bias", KeyType::Optional);

	// Output control
	input_parameter_pack_.readFlag("PrintOutput", KeyType::Optional, want_output_);
	input_parameter_pack_.readString("OutputFileSuffix", KeyType::Optional, output_file_suffix_);

	// Computing and printing forces
	input_parameter_pack_.readFlag("PrintForces", KeyType::Optional, print_forces_);
	input_parameter_pack_.readFlag("PrintNumericalForces", KeyType::Optional, print_numerical_forces_);
	input_parameter_pack_.readFlag("NeedForces", KeyType::Optional, need_derivatives_);
	if ( print_forces_ or print_numerical_forces_ ) {
		// Ensure derivatives are computed
		need_derivatives_ = true;
	}

	// Input files
	input_parameter_pack_.readString("XtcFile", KeyType::Optional, xtc_file_);
	input_parameter_pack_.readString("GroFile", KeyType::Optional, gro_file_);
	input_parameter_pack_.readString("TopFile", KeyType::Optional, top_file_);

	// Sampling phase (only relevant for post-processing)
	input_parameter_pack_.readNumber("t0", KeyType::Optional, t_min_);
	input_parameter_pack_.readNumber("tf", KeyType::Optional, t_max_);


	//----- Check input -----//

	// TODO extra checks here?


	//----- Initialize input-dependent objects -----//

	try {
		// Biases: initialize objects and infer whether derivatives are needed
		std::string token, lowercase_token;
		int num_bias_packs = bias_pack_ptrs.size();
		for ( int i=0; i<num_bias_packs; ++i ) {
			// Figure out which Bias block corresponds to which order parameter
			bias_pack_ptrs[i]->readString("order_parameter", KeyType::Required, token);
			lowercase_token = string_tools_.toLowercase(token);

			if ( lowercase_token == "ntilde" ) {
				// Bias on ntilde_v
				if ( bias_ntilde_v_ptr_ == nullptr ) {
					bias_ntilde_v_ptr_ = std::unique_ptr<Bias>( 
							new Bias( *bias_pack_ptrs[i], simulation_state_.get_simulation_box()) );
					need_derivatives_ = true;
				}
				else {
					throw std::runtime_error("Multiple Biases on ntilde_v are not allowed");
				}
			}
			else {
				std::string err_msg = "Error: unrecognized Bias order parameter \"" + token + "\".\n";
				throw std::runtime_error( err_msg );
			}
		}

		// Topology (for selecting target atoms)
		Topology& topology = simulation_state_.access_topology();
		if ( (not top_file_.empty()) and (not gro_file_.empty()) ) {
			topology.setFromFile(top_file_, gro_file_);
		}
		else if ( not top_file_.empty() ) {
			topology.setFromFile(top_file_);
		}
		else if ( not gro_file_.empty() ) {
			topology.setFromFile(gro_file_);
		}

		// Probe volume 
		probe_volume_pack_ptr->readString("type", KeyType::Required, probe_volume_type_);
		ProbeVolume::ProbeVolumeInputPack probe_volume_input_pack = { 
			*probe_volume_pack_ptr,
			simulation_state_,
			mpi_communicator_,
			need_derivatives_
		};
		probe_volume_ptr_ = std::shared_ptr<ProbeVolume>(
				GenericFactory<ProbeVolume, const std::string, ProbeVolume::ProbeVolumeInputPack
				               >::factory().create(probe_volume_type_, probe_volume_input_pack) );
	}
	catch (const std::exception& e) {
		std::cerr << "Indus constructor - Exception occurred (" << e.what() << ")\n";
		throw;
	}


	//----- Target Atoms -----//

	// Get target atom indices
	Topology& topology = simulation_state_.access_topology();
	topology.getTargets(target_tokens,
	                    target_atoms_, target_flags_);
	if ( target_atoms_.size() == 0 ) {
		std::stringstream err_ss;
		err_ss << "Error in Indus constructor - No target atoms were chosen.\n";
		throw std::runtime_error( err_ss.str() );
	}

	this->setTargetAtomsByGlobalIndices( target_atoms_ );


	//----- Domain Decomposition ------//

	// Set the number of cells along each axis
	std::array<int, DIM_> dd_grid_dimensions = mpi_communicator_.calculateGridDimensions<DIM_>();
	domain_decomposition_.set_grid_dimensions( dd_grid_dimensions );

	// Set shell parameters
	double width_shell_1 = 0.0, width_shell_2 = 0.0;
	domain_decomposition_.setShellWidths(width_shell_1, width_shell_2);


	//----- Derivatives Setup -----//

	if ( need_derivatives_ ) {
		// Allocate memory
		int num_indus_atoms = indus_atom_global_indices_.size();
		derivatives_ntilde_v_.assign( num_indus_atoms, CommonTypes::zero_3 );
		derivatives_u_bias_total_.assign( num_indus_atoms, CommonTypes::zero_3 );

		// Zero out all the virial arrays
		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				sum_r_cross_dntilde_v_dr_[a][b] = 0.0;
				virial_u_bias_total_[a][b] = 0.0;
			}
		}
	}

	// TODO extra checks? feedback to user?
}


void Indus::setTargetAtomsByGlobalIndices(const std::vector<int>& target_atoms)
{
	target_atoms_ = target_atoms;

	// Set target atom indices in group
	AtomGroup& target_atom_group = simulation_state_.access_target_atom_group();
	target_atom_group.set_global_atom_indices( target_atoms );

	//----- Indus atoms -----//

	indus_atom_global_indices_ = target_atoms;

	// Mapping from target atoms to Indus atoms
	int num_target_atoms = target_atoms.size();
	target_atom_indus_indices_.resize(num_target_atoms);
	for ( int i=0; i<num_target_atoms; ++i ) {
		target_atom_indus_indices_[i] = i;
	}
}


#ifndef INDUS_PLUMED_MODE
// Calculate P_v(N) for a generic order parameter q
void Indus::do_indus_standalone()
{

	//----- Basic input checks -----//

	if ( xtc_file_.empty() ) {
		throw std::runtime_error("error in do_indus_standalone: no XTC file provided.");
	}

	//----- Perform setup -----//

	// XTC variables
	int    step;			     // Simulation step counter
	float  xdrfile_time;   // Simulation time
	matrix xdrfile_box;    // Box coordinates in a 3x3 matrix
	float  xdrfile_prec;   // Precision of the xtc file
	RvecArray& xtc_atom_positions = simulation_state_.access_atom_positions();

	// C++ versions of XTC variables
	Box    box_matrix;
	double time;

	// Get total number of atoms from xtc file
	// - Create a copy of the file name to safely cast away the const qualifier from string::c_str
	//   without endagering the original string
	// TODO check numAtoms vs. Topology, if possible
	int num_atoms;
	std::string xtc_file_copy = xtc_file_;
	int flag = read_xtc_natoms(const_cast<char*>(xtc_file_copy.c_str()), &num_atoms);
	if ( flag != exdrOK ) {
		std::stringstream err_ss;
		err_ss << "Error in " << __PRETTY_FUNCTION__ << "(" << __FILE__ << ":" << __LINE__ << ")\n"
		       << "  read_xtc_natoms threw error code " << flag << "\n"
		       << "  XTC file:  " << xtc_file_ << "\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Open XTC file
	XDRFILE* xtc_file_ptr = nullptr;
	if ( my_rank_ == master_rank_ ) {
		// Only master rank reads the file 
		xtc_file_ptr = xdrfile_open(xtc_file_.c_str(), "r");

		if ( xtc_file_ptr == nullptr ) {
			std::stringstream err_ss;
			err_ss << "  Indus::do_indus_standalone: Unable to open xtc file." << "\n";
			throw std::runtime_error( err_ss.str() );
		}
	}

	// Allocate memory for xtc_atom_positions (and derivatives, if requested)
	try { 
		xtc_atom_positions.resize(num_atoms);
	}
	catch (const std::bad_alloc& ba) {
    std::cerr << "  Indus::do_indus_standalone: Unable to allocate memory." << "\n"
		          << "  .what() = " << ba.what() << "\n";

		// Close the xtc file (don't leave the file in a bad state!)
		xdrfile_close( xtc_file_ptr );
    throw;
	}

	// Time series
	std::vector<double> t_samples, ntilde_v_samples; 
	std::vector<int>    n_v_samples;

	// For printing biasing forces for analyzed frames
	std::ofstream forces_ofs, numerical_forces_ofs;
	if ( (need_derivatives_ or print_forces_ or print_numerical_forces_ ) and 
	     want_output_ and my_rank_ == master_rank_ 
	) {
		// Header
		std::stringstream header_ss;
		header_ss << "# Biasing forces on Indus atoms\n"
		          << "#   Units: t(ps), box_matrix(nm), virial(kJ/mol), force(kJ/mol/nm)\n"
		          << "#   Each line starts with the global serial number of the affected atom\n"
		          << "N_total= " << num_atoms << "  N_biased= " 
		             << indus_atom_global_indices_.size() << "\n";

		if ( print_forces_ ) {
			forces_ofs.open("forces_debug.out");
			if ( not forces_ofs.is_open() ) {
				throw std::runtime_error("couldn't open file for printing forces");
			}
			forces_ofs << header_ss.str();
		}
		
		if ( print_numerical_forces_ ) {
			numerical_forces_ofs.open("numerical_forces_debug.out");
			if ( not numerical_forces_ofs.is_open() ) {
				throw std::runtime_error("couldn't open file for printing forces");
			}
			numerical_forces_ofs << header_ss.str();
		}
	}


	//----- Calculate N_v and Ntilde_v -----//

	int num_frames = 0, num_samples_total_n_v = 0;
	double avg_n_v = 0.0, avg_n_v_sq = 0.0;

	// User feedback
	if ( my_rank_ == master_rank_ ) {
		std::cout << "  INDUS calculation initialized." << "\n"
		          << "  Beginning to read XTC file for atom positions." << "\n"
		          << "  Number of atoms: " << "\t" << num_atoms << "\n"
		          << "  Number of target atoms: " << "\t" << target_atoms_.size() << "\n";
	}

	int reading_xtc_file = 1;
	int analyze_frame    = 1;

	while (true) {
		// Master rank reads the XTC file
		if ( my_rank_ == master_rank_ ) {
			int xdr_return_code = read_xtc(xtc_file_ptr, num_atoms,
			                               &step, &xdrfile_time, xdrfile_box, 
			                               xtc_atom_positions.data(), &xdrfile_prec);
			// Convert to double precision
			time = static_cast<double>(xdrfile_time);

			// Convert box to C++ type
			for ( int a=0; a<DIM_; ++a ) {
				for ( int b=0; b<DIM_; ++b ) {
					box_matrix[a][b] = xdrfile_box[a][b];
				}
			}

			if ( time > t_max_ or xdr_return_code != exdrOK ) {
				// Done reading the XTC file
				reading_xtc_file = 0;
				analyze_frame    = 0;
			}
			else if ( time >= t_min_ ) {
				// In production phase
				analyze_frame = 1;
			}
			else {  // time < t_min_
				// Before production phase
				analyze_frame = 0;
			}
		}


		if ( mpi_communicator_.is_mpi_initialized() ) {
			// Share flags so that all ranks know:
			//  (1) whether to break (done reading file); and
			//  (2) whether to analyze the frame
			std::array<int,2> flags = {{ reading_xtc_file, analyze_frame }};
			mpi_communicator_.bcast(flags, master_rank_);

			reading_xtc_file = flags[0];
			analyze_frame    = flags[1];
		}

		// Determine whether to analyze the frame, continue, or break
		if ( reading_xtc_file == 0 ) {
			break; 
		}
		else if ( analyze_frame == 1 ) {
			num_frames++;

			// Update simulation state
			// - Box lengths
			// - Atomic positions
			// - Share data via MPI as necessary
			this->updateSimulationState(xtc_atom_positions, box_matrix, step, time);

			// Calculate q (and/or Ntilde_v)
			this->calculate();

			// Get INDUS order parameters from probe volume
			const int    n_v      = probe_volume_ptr_->get_n_v();
			const double ntilde_v = probe_volume_ptr_->get_ntilde_v();

			// User feedback
			if ( (num_frames % 500) == 0 and my_rank_ == master_rank_ ) {
				std::cout << "Time[ps] " << time << ":  N_v = "
				          << n_v << ",  ntilde_v = " << ntilde_v << "\n";
			}

			// Process results for N_v
			avg_n_v    += static_cast<double>(n_v);
			avg_n_v_sq += static_cast<double>(n_v*n_v);
			++num_samples_total_n_v;

			// Store time series of samples
			t_samples.push_back( time );
			n_v_samples.push_back( n_v );
			ntilde_v_samples.push_back( ntilde_v );


			//----- Printing Forces -----//

			// Print analytic forces
			if ( print_forces_ and my_rank_ == master_rank_ ) {
				printForcesForFrame( forces_ofs, time, box_matrix, derivatives_u_bias_total_,
														 virial_u_bias_total_ );
			}

			// Numerical forces
			if ( print_numerical_forces_ ) {
				this->computeNumericalDerivatives();

				if ( my_rank_ == master_rank_ ) {
					printForcesForFrame( numerical_forces_ofs, time, box_matrix, numerical_derivatives_u_bias_total_,
															 numerical_virial_u_bias_total_ );
				}
			}
		} // for frames in production phase
	} // loop over xtc frames

	// Close files
	if ( my_rank_ == master_rank_ ) {
		xdrfile_close( xtc_file_ptr );
		if ( forces_ofs.is_open() ) {
			forces_ofs.close();
		}
		if ( numerical_forces_ofs.is_open() ) {
			numerical_forces_ofs.close();
		}
	}


	//----- Calculate statistics -----//

	double num_samples_total_n_v_d = static_cast<double>(num_samples_total_n_v);
	avg_n_v    /= num_samples_total_n_v_d;
	avg_n_v_sq /= num_samples_total_n_v_d;
	double var_n_v = num_samples_total_n_v_d/(num_samples_total_n_v_d - 1.0)
	                  *( avg_n_v_sq - avg_n_v*avg_n_v );


	//----- Write output files -----//

	if ( want_output_ and my_rank_ == master_rank_ ) {
		// Suggested biasing parameters (kappa-values in kBT)
		double kappa_0_n_v    = 1.0/var_n_v;
		double kappa_n_v      = 2.0*kappa_0_n_v;
		double delta_n_v_star = 4.0/sqrt(3.0)*sqrt(var_n_v);

		// Generic header for output files associated with the calculation
		std::stringstream header_stream;
		header_stream << this->getInputSummary()
		             << "#\n"
		             << "# Output\n"
		             << "#   Number of particles\n"
		             << "#     Avg_n_v     " << avg_n_v       << "\n"
		             << "#     Var_n_v     " << var_n_v       << "\n"
		             << "#     Stdev_n_v   " << sqrt(var_n_v) << "\n"
		             << "#     kappa_0_n_v " << kappa_0_n_v   << " [kBT]\n"
		             << "#   Suggested biasing parameters (2*kappa_0 and 4kBT overlap rules of thumb)\n"
		             << "#     kappa_n_v      " << kappa_n_v      << " [kBT]\n"
		             << "#     delta_n_v_star " << delta_n_v_star << "\n"
		             << " " << "\n";

		//----- Write time series of N_v and Ntilde_v to file -----//

		std::ostringstream timeSeriesFileName;
		timeSeriesFileName << "time_samples_indus.out";
		std::ofstream ofs( timeSeriesFileName.str().c_str() );

		ofs << "# RESULTS: Time series of N_v and Ntilde_v" << "\n";
		ofs << header_stream.str();

		// Table header
		ofs << "# Time[ps]" << "\t" << "N_v" << "\t" << "Ntilde_v\n";

		// Table
		int num_times = n_v_samples.size();
		for ( int i=0; i<num_times; ++i ) {
			ofs << t_samples[i] << "\t" << n_v_samples[i] << "\t" << ntilde_v_samples[i] << "\n";
		}

		ofs.close();

	} // end if ( want_output )

	return;
}
#endif // #ifndef INDUS_PLUMED_MODE


void Indus::updateSimulationState(
		const RvecArray& atom_positions, const Box& box_matrix, const int step, const double time)
{
	//----- Update atom positions in target atom group(s) -----//

	AtomGroup& target_atom_group = simulation_state_.access_target_atom_group();
	std::vector<Real3>& target_atom_positions = target_atom_group.access_atom_positions();

	// Reserve space for target atom positions
	int num_target_atoms = target_atom_group.get_num_atoms();
	target_atom_positions.resize(num_target_atoms);

	// Create local copies that can be changed in standalone mode
	Box    local_box_matrix = box_matrix;
	int    local_step = step;
	double local_time = time;

#ifndef INDUS_PLUMED_MODE
	//----------------------------//
	//----- Standalone Mode ------//
	//----------------------------//

	if ( my_rank_ == master_rank_ ) {
		// Copy target atoms positions from XTC coordinates
		// - atom_positions contains global set of positions
		const std::vector<int>& global_atom_indices = target_atom_group.get_global_atom_indices();
		for ( int k=0; k<num_target_atoms; ++k ) {
			int global_index = global_atom_indices[k];
			for ( int d=0; d<DIM_; ++d ) {
				target_atom_positions[k][d] = atom_positions[global_index][d];
			}
		}
	}

	if ( mpi_communicator_.is_mpi_initialized() ) {
		//----- MPI: Share data -----//

		// Share target atom positions
		mpi_communicator_.bcast(target_atom_positions, master_rank_);

		// Share box
		mpi_communicator_.bcast(local_box_matrix, master_rank_);

		// Share step and time
		mpi_communicator_.bcast(local_step, master_rank_);
		mpi_communicator_.bcast(local_time, master_rank_);

	}

#else
	//------------------------//
	//----- PLUMED Mode ------//
	//------------------------//

	// - No further communication necessary (all ranks already have the necessary 
	//   atom positions and access to the simulation box)
	// - Notes on the input:
	//    - atom_positions contains atoms requested from PLMD core, which are
	//      the Indus atoms
	//    - Box (above) was copied from the PLUMED box

	for ( int k=0; k<num_target_atoms; ++k ) {
		for ( int d=0; d<DIM_; ++d ) {
			target_atom_positions[k][d] = atom_positions[k][d];
		}
	}

#endif // #ifndef INDUS_PLUMED_MODE

	//------ Other State Variables -----//

	// Simulation box (ssumes orthonormal box)
	SimulationBox& simulation_box = simulation_state_.access_simulation_box();
	simulation_box.setLengths( local_box_matrix[X_DIM][X_DIM], 
	                           local_box_matrix[Y_DIM][Y_DIM], 
	                           local_box_matrix[Z_DIM][Z_DIM] );

	simulation_state_.set_step( local_step );
	simulation_state_.set_time( local_time );


	//----- Update Domain Decomposition -----//

	// TODO update the number of cells in each direction based
	//      on box dimensions?
	// - e.g. Put the most cells along the longest box dimension

	// Update the ProbeVolume's bounding box
	probe_volume_ptr_->setBoundingBox();

	// Set DD region based on the ProbeVolume
	Real3 bounding_box_offset;
	Box   bounding_box_matrix;
	probe_volume_ptr_->getBoundingBox(bounding_box_offset, bounding_box_matrix);
	domain_decomposition_.setCellGridRegion(bounding_box_offset, bounding_box_matrix);

	// Set DD region to be equal to SimulationBox
	//domain_decomposition_.setCellGridRegion(); // OLD

	// Partition the simulation box into cells
	domain_decomposition_.makeCells();

	// Find the atoms which are local to this processor
	domain_decomposition_.findNearbyAtoms(my_rank_, target_atom_group);

	//----- DEBUG -----//

	/* 
	domain_decomposition_.printAtomGroupDomainDecompositionCheck(std::cout, target_atom_group);
	*/

	//----- END DEBUG -----//
}


void Indus::calculate()
{
	if ( need_derivatives_ ) {
		derivatives_are_ready_ = false;
	}

	// Perform INDUS
	probe_volume_ptr_->doIndusWithShells();

	if ( need_derivatives_ ) {
		// Gather all derivatives of q and/or ntilde_v across all ranks,
		// and make them available on all ranks
		this->gatherDerivatives();
	}

	// Compute biasing potentials and their derivatives 
	// (the negatives of the biasing forces)
	this->calculateTotalBias();

	if ( need_derivatives_ ) {
		derivatives_are_ready_ = true;
	}

	return;
}


void Indus::gatherDerivatives()
{
	const AtomGroup& target_atom_group = simulation_state_.get_target_atom_group();

	// Bias on ntilde_v
	if ( need_derivatives_ ) {
		// Derivatives of ntilde_v computed locally
		const std::vector<Real3>& local_derivatives_ntilde_v = 
				probe_volume_ptr_->getIndicatorFunctionDerivatives();
		const std::vector<int>& local_group_indices =
				probe_volume_ptr_->get_group_indices_of_local_indus_atoms();

		if ( mpi_communicator_.is_mpi_initialized() ) {
			// Share local derivatives and corresponding group indices across all processors
			mpi_communicator_.allgatherv( 
					local_derivatives_ntilde_v,
					// Output
					derivatives_ntilde_v_buffer_.buffer, derivatives_ntilde_v_buffer_.block_offsets,
					derivatives_ntilde_v_buffer_.block_sizes );
			mpi_communicator_.allgatherv( 
					local_group_indices,
					// Output
					group_indices_for_ntilde_v_buffer_.buffer, group_indices_for_ntilde_v_buffer_.block_offsets,
					group_indices_for_ntilde_v_buffer_.block_sizes );
		}
		else {
			int num_ranks = 1;
			derivatives_ntilde_v_buffer_.buffer        = local_derivatives_ntilde_v;
			derivatives_ntilde_v_buffer_.block_offsets = std::vector<int>(num_ranks, 0);
			derivatives_ntilde_v_buffer_.block_sizes   = std::vector<int>(num_ranks, 1);
			group_indices_for_ntilde_v_buffer_.buffer        = local_group_indices;
			group_indices_for_ntilde_v_buffer_.block_offsets = std::vector<int>(num_ranks, 0);
			group_indices_for_ntilde_v_buffer_.block_sizes   = std::vector<int>(num_ranks, 1);
		}

		this->sumGatheredDerivatives(
				derivatives_ntilde_v_buffer_, group_indices_for_ntilde_v_buffer_, 
				target_atom_group, target_atom_indus_indices_,
				// Output
				derivatives_ntilde_v_, sum_r_cross_dntilde_v_dr_ );
	}
}


void Indus::sumGatheredDerivatives(
		const AllgathervBuffer<Real3>& derivatives_x_buffer,
		const AllgathervBuffer<int>& group_indices_for_x_buffer,
		const AtomGroup& atom_group, const std::vector<int>& op_indices_of_group_atoms,
		// Output
		std::vector<Real3>& derivatives_x, Box& sum_r_cross_dx_dr)
{
	// TODO check size agreement between derivatives and indices?

	// Prepare output variables
	int num_indus_atoms = indus_atom_global_indices_.size();
	derivatives_x.assign( num_indus_atoms, CommonTypes::zero_3 );
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			sum_r_cross_dx_dr[a][b] = 0.0;
		}
	}

	// Sum over all contributions to the total derivatives
	int block_offset, block_size, group_index, buffer_offset, indus_index;
	int num_ranks = mpi_communicator_.get_size();
	for ( int r=0; r<num_ranks; ++r ) {
		block_offset = derivatives_x_buffer.block_offsets[r];
		block_size   = derivatives_x_buffer.block_sizes[r];

		for ( int i=0; i<block_size; ++i ) {
			buffer_offset = block_offset + i;
			group_index   = group_indices_for_x_buffer.buffer[buffer_offset];
			indus_index = op_indices_of_group_atoms[group_index];

			// Target atoms come first in IndusAtoms arrays
			for ( int d=0; d<DIM_; ++d ) {
				derivatives_x[indus_index][d] += derivatives_x_buffer.buffer[buffer_offset][d];
			}
		}
	}

	// Matrix for the virial
	const std::vector<Real3>& atom_positions = atom_group.get_atom_positions();
	for ( int k=0; k<num_indus_atoms; ++k ) {
		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				sum_r_cross_dx_dr[a][b] += atom_positions[k][a]*derivatives_x[k][b];
			}
		}
	}
}


void Indus::calculateTotalBias()
{
	// Prepare output variables
	u_bias_ntilde_v_ = 0.0;
	u_bias_total_    = 0.0;
	if ( need_derivatives_ ) {
		int num_indus_atoms = indus_atom_global_indices_.size();
		derivatives_u_bias_total_.assign(num_indus_atoms, CommonTypes::zero_3);
		virial_u_bias_total_ = CommonTypes::zero_matrix;
	}

	// Bias on Ntilde_v
	if ( bias_ntilde_v_ptr_ != nullptr ) {
		const double ntilde_v = probe_volume_ptr_->get_ntilde_v();

		bias_ntilde_v_ptr_->accumulateBias(
				ntilde_v, need_derivatives_, derivatives_ntilde_v_, sum_r_cross_dntilde_v_dr_,
				// Output
				u_bias_ntilde_v_, u_bias_total_, derivatives_u_bias_total_, virial_u_bias_total_
		);
	}
}


// Compute P_v(N) from a series of numbers with bins over the range [0,max_n],
// where each bin "b" spans the range [b,b+1)
void Indus::make_PvN(
	const std::vector<double>& n_values, const int max_n, 
	std::vector<double>& bins, std::vector<double>& p_v_n,
  double& avg_n, double& var_n
) const
{
	if ( max_n < 0 ) {
		std::stringstream err_ss;
		err_ss << "Indus::make_PvN: ERROR: max(N) can't be less than zero (input: max_n="
		          << max_n << ").\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Set up bins
	int num_bins = max_n + 1;
	bins.resize(num_bins);
	for ( int b=0; b<num_bins; ++b ) {
		bins[b] = static_cast<double>( b );
	}
	p_v_n.resize(num_bins, 0.0);

	avg_n = 0.0;
	var_n = 0.0;

	int numValues  = static_cast<int>( n_values.size() );
	int num_samples = 0;
	int    n;
	double n_d;

	for ( int i=0; i<numValues; ++i ) {
		n = n_values[i];
		if ( n >= 0 and n <= max_n ) {
			p_v_n[n] += 1.0;

			n_d    = static_cast<double>(n);
			avg_n += n_d;
			var_n += n_d*n_d;
			++num_samples;
		}
	}

	// Normalize bins (note that bin size = 1)
	double num_samples_d = static_cast<double>( num_samples );
	for ( int b=0; b<num_bins; ++b ) {
		p_v_n[b] /= num_samples_d;
	}

	// Finalize statistics
	avg_n /= num_samples_d; 
	var_n = num_samples_d/(num_samples_d - 1.0)*( var_n/num_samples_d - avg_n*avg_n );

	return;
}


// Returns a big string with a whole bunch of input variables.
std::string Indus::getInputSummary() const 
{
	std::stringstream ss;
	ss << "# Indus Input\n"
	   << "#\n";
	if ( mpi_communicator_.is_mpi_enabled() ) {
		ss << "# MPI enabled\n"
		   << "#   num_ranks = " << mpi_communicator_.get_size() << "\n";
	}
	else {
		ss << "# MPI disabled\n";
	}
	ss << "#\n"
	   << probe_volume_ptr_->getInputSummary()
	   << "#\n"
	   << "# Target_Atoms_String = " << target_atoms_string_ << "\n"
	   << "# Gro_File = " << gro_file_ << "\n"
		 << "# Xtc_File = " << xtc_file_ << "\n"
	   << "#   Production_phase\n"  
		 << "#     Start = " << t_min_ << " [ps]\n"
		 << "#     Stop  = " << t_max_ << " [ps]\n"
	   << "#\n"
	   << "# Biasing_potentials\n"
	   << "#   Need_Derivatives = " << need_derivatives_ << "\n";
	if ( bias_ntilde_v_ptr_ != nullptr ) {
		ss << "#   Bias_ntilde_v\n"
		   << bias_ntilde_v_ptr_->getInputSummary("#    ");
	}

	return ss.str();
}


void Indus::printForcesForFrame(
		std::ofstream& forces_ofs, const double time, const Matrix& box_matrix, 
		const std::vector<Real3>& derivatives_u_bias_total, const Matrix& virial_u_bias_total
) const 
{
	// Simulation state
	forces_ofs << "t= " << time << "\n";
	forces_ofs << "box= ";
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			forces_ofs << " " << std::fixed << std::setprecision(5) << box_matrix[a][b];
		}
	}
	forces_ofs << "\n";

	// Virial
	forces_ofs << "virial= ";
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			forces_ofs << " " << std::fixed << std::setprecision(5) << virial_u_bias_total[a][b];
		}
	}
	forces_ofs << "\n";

	// Biasing forces (negatives of total biasing potential)
	int num_indus_atoms = indus_atom_global_indices_.size();
	for ( int k=0; k<num_indus_atoms; ++k ) {
		forces_ofs << indus_atom_global_indices_[k] + 1;
		//for ( int d=0; d<DIM_; ++d ) { forces_ofs << "  " << xtc_atom_positions[k][d]; }
		for ( int d=0; d<DIM_; ++d ) { 
			forces_ofs << "  " << -derivatives_u_bias_total[k][d]; 
		}
		forces_ofs << "\n";
	}
}


void Indus::computeNumericalDerivatives()
{
	// Save initial values of key variables
	const double ntilde_v_0        = probe_volume_ptr_->get_ntilde_v();
	const double u_bias_ntilde_v_0 = u_bias_ntilde_v_;
	const double u_bias_total_0    = u_bias_total_;

	// Allocate memory and set everything to zero
	int num_indus_atoms = indus_atom_global_indices_.size();
	numerical_derivatives_ntilde_v_.assign( num_indus_atoms, CommonTypes::zero_3 );
	numerical_derivatives_u_bias_total_.assign( num_indus_atoms, CommonTypes::zero_3 );
	for ( int a=0; a<DIM_; ++a ) {
		numerical_virial_u_bias_total_[a] = CommonTypes::zero_3;
	}

	//----- Target atoms -----//

	// Access positions array
	AtomGroup& target_atom_group = simulation_state_.access_target_atom_group();
	std::vector<Real3>& target_atom_positions = target_atom_group.access_atom_positions();
	int num_target_atoms = target_atom_group.get_num_atoms();

	for ( int i=0; i<num_target_atoms; ++i ) {
		// Save actual position
		Real3 x_i_0 = target_atom_positions[i];

		try {
			int indus_index = target_atom_indus_indices_[i];

			for ( int d=0; d<DIM_; ++d ) {
				// Perturb position
				target_atom_positions[i][d] += delta_x_;

				// Recompute everything
				this->calculate();
				const double ntilde_v_new = probe_volume_ptr_->get_ntilde_v();

				// Accumulate numerical derivatives
				numerical_derivatives_ntilde_v_[indus_index][d] += 
						(ntilde_v_new - ntilde_v_0)/delta_x_;
				numerical_derivatives_u_bias_total_[indus_index][d] += 
						(u_bias_total_ - u_bias_total_0)/delta_x_;
				for ( int a=0; a<DIM_; ++a ) {
					numerical_virial_u_bias_total_[a][d] += 
							0.5*x_i_0[a]*numerical_derivatives_u_bias_total_[indus_index][d];
				}

				// Put atom back where it started
				target_atom_positions[i][d] = x_i_0[d];
			}
		}
		catch (const std::exception& ex) {
			// Replace original position
			target_atom_positions[i] = x_i_0;

			// Rethrow
			throw;
		}
	}
}
