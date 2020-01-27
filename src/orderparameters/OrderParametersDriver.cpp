#include "OrderParametersDriver.h"

OrderParametersDriver::OrderParametersDriver(const std::string& input_file, MpiCommunicator& mpi_communicator):
	input_file_(input_file),
	// Simulation state
	simulation_state_(),
	mpi_communicator_( mpi_communicator ),
	// Domain decomposition/MPI (DD options are overridden later)
	domain_decomposition_( simulation_state_, mpi_communicator ),
	dd_cell_grid_( domain_decomposition_.access_cell_grid() ),
	my_rank_( mpi_communicator_.get_rank() ),
	master_rank_( mpi_communicator_.get_master_rank() ),
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
	setup_timer_.start();

	// Checks whether the GenericFactory objects of interest have been populated
	checkGenericFactoryInitialization();

#ifdef ORDER_PARAMETERS_DRIVER_PLUMED_MODE
	// Generally need derivatives if using Steinardt with PLUMED
	// (almost always biasing something)
	need_derivatives_ = true;

	// By default, do not share all derivatives for best performance
	share_all_op_derivatives_ = false;
#endif

	xtc_file_.clear();
	gro_file_.clear();
	top_file_.clear();

	bias_ntilde_v_string_.clear();


	//----- Process input -----//

	// Read input file
	InputParser input_parser;
	input_parser.parseFile(input_file, input_parameter_pack_);

	using KeyType = ParameterPack::KeyType;
	std::string line, token, lowercase_token;

	// Check for debug mode
	bool is_debug_mode = false;
	bool found = input_parameter_pack_.readFlag("Debug", KeyType::Optional, is_debug_mode);
	if ( found ) {
		simulation_state_.set_debug_mode( is_debug_mode );
	}

	// Probe volume
	const ParameterPack* probe_volume_pack_ptr = 
			input_parameter_pack_.findParameterPack("ProbeVolume", KeyType::Required);

	// Check whether this calculation is INDUS only
	input_parameter_pack_.readFlag("OnlyDoIndus", KeyType::Optional, only_do_indus_);

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


	//----- TODO MULTIPLE OPS -----//

	int num_ops = 1;
	index_indus_ = 0;

	op_ptrs_.assign(num_ops, nullptr);

	// Derivatives handling
	local_derivatives_ops_.resize(num_ops);
	local_op_indices_for_derivatives_.resize(num_ops);
	sum_r_cross_op_derivatives_.resize(num_ops);

	// Biasing potentials
	bias_ptrs_.assign(num_ops, nullptr);
	u_bias_ops_.resize(num_ops);
	local_derivatives_u_bias_ops_.resize(num_ops);
	virial_u_bias_ops_.resize(num_ops);


	//----- Initialize input-dependent objects -----//

	// First, establish the Topology
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

	// Set up all AtomGroups, which order parameters may need below
	initializeAtomGroups();

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
							new Bias( *bias_pack_ptrs[i], simulation_state_ ) );
					need_derivatives_ = true;
					bias_ptrs_[index_indus_] = bias_ntilde_v_ptr_.get();
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

		// Probe volume (will be shared with various OrderParameters)
		probe_volume_pack_ptr->readString("type", KeyType::Required, probe_volume_type_);
		ProbeVolume::ProbeVolumeInputPack probe_volume_input_pack = { 
			*probe_volume_pack_ptr,
			simulation_state_,
			mpi_communicator_,
			need_derivatives_
		};
		probe_volume_ptr_ = std::shared_ptr<ProbeVolume>(
			ProbeVolumeRegistry::Factory::factory().create(probe_volume_type_, probe_volume_input_pack) 
		);

		// Input to each OrderParameter
		OrderParameter::InputPack op_input_pack {
			simulation_state_, domain_decomposition_,
			mpi_communicator_, need_derivatives_
		};
	
		indus_ptr_ = std::unique_ptr<Indus>(
			new Indus(probe_volume_ptr_, probe_volume_input_pack, op_input_pack)
		);
	}
	catch (const std::exception& e) {
		std::cerr << "OrderParametersDriver constructor - Exception occurred (" << e.what() << ")\n";
		throw;
	}

	// TODO MULTIPLE OPS
	op_ptrs_[index_indus_] = dynamic_cast<OrderParameter*>( indus_ptr_.get() );
	for ( int i=0; i<num_ops; ++i ) {
		if ( op_ptrs_[i] == nullptr ) {
			std::stringstream err_ss;
			err_ss << "Error in " << FANCY_FUNCTION << "\n"
						 << "  dynamic_cast to base class OrderParameter failed unexpectedly\n";
			throw std::runtime_error( err_ss.str() );
		}
	}


	//----- Domain Decomposition ------//

	// Set the number of cells along each axis
	// TODO Adjust which axis has the most cells based on the length of the (bounding) box
	//      in that direction?
	std::array<int, DIM_> dd_grid_dimensions = mpi_communicator_.calculateGridDimensions<DIM_>();
	dd_cell_grid_.set_grid_dimensions( dd_grid_dimensions );

	// Set shell parameters
	// - It's important to do this after setting up the Steinhardt parameter, which determines how
	//   large the shells must be
	double width_shell_1, width_shell_2;
	probe_volume_ptr_->getShellWidths(width_shell_1, width_shell_2);

	dd_cell_grid_.setShellWidths(width_shell_1, width_shell_2);

	setup_timer_.stop();
}


void OrderParametersDriver::initializeAtomGroups()
{
	if ( simulation_state_.debug_mode() ) {
		std::cout << "INITIALIZE ATOM GROUPS (" << LOCATION_IN_SOURCE_STRING << ")" << std::endl;
	}

	using KeyType = ParameterPack::KeyType;

	// Parse AtomGroup input packs
	std::vector<const ParameterPack*> atom_group_pack_ptrs = 
			input_parameter_pack_.findParameterPacks("AtomGroup", KeyType::Optional);
	int num_packs = atom_group_pack_ptrs.size();
	for ( int i=0; i<num_packs; ++i ) {
		simulation_state_.addAtomGroup( *atom_group_pack_ptrs[i] );
	}

	//----- Backwards compatibility (DEPRECATED SYNTAX) -----//

	// Target atoms
	std::string standard_name = "target_atoms";
	std::vector<std::string> tokens;
	bool found_target = input_parameter_pack_.readVector("Target", KeyType::Optional, tokens);
	if ( found_target and (not simulation_state_.atomGroupExists(standard_name))) {
		if ( simulation_state_.debug_mode() ) {
			std::cout << "MAKING DEFAULT TARGET ATOMGROUP (" << LOCATION_IN_SOURCE_STRING << ")" << std::endl;
		}

		// Create a group with the default name
		ParameterPack atom_group_pack("AtomGroup");
		atom_group_pack.values.insert( ParameterPack::ValuePair("name", standard_name) );
		atom_group_pack.vectors.insert( ParameterPack::VectorPair("selection", tokens) );
		simulation_state_.addAtomGroup(atom_group_pack);
	}


	//----- Save index mappings -----//

	const std::vector<AtomGroup>& atom_groups = simulation_state_.get_atom_groups();
	int num_atom_groups = atom_groups.size();
	op_atom_global_indices_.resize(0);
	atom_group_to_op_indices_.resize(num_atom_groups);
	int counter = 0;
	for ( int k=0; k<num_atom_groups; ++k ) {
		// Combined array of all global indices across all groups (in order)
		const std::vector<int>& global_indices = atom_groups[k].get_global_atom_indices();
		op_atom_global_indices_.insert( op_atom_global_indices_.end(),
		                                global_indices.begin(), global_indices.end() );

		// For each group: mapping from group index to OP index
		int num_atoms_in_group = global_indices.size();
		atom_group_to_op_indices_[k].resize(num_atoms_in_group);
		for ( int i=0; i<num_atoms_in_group; ++i ) {
			atom_group_to_op_indices_[k][i] = counter;
			++counter;
		}
	}

	if ( simulation_state_.debug_mode() ) {
		std::cout << "DONE INITIALIZING ATOM GROUPS (" << LOCATION_IN_SOURCE_STRING << ")" << std::endl;
	}
}


#ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE
void OrderParametersDriver::run()
{

	//----- Basic input checks -----//

	if ( xtc_file_.empty() ) {
		throw std::runtime_error("error in run(): no XTC file provided.");
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
	// TODO check num_atoms vs. Topology, if possible
	int num_atoms;
	std::string xtc_file_copy = xtc_file_;
	int flag = read_xtc_natoms(const_cast<char*>(xtc_file_copy.c_str()), &num_atoms);
	if ( flag != exdrOK ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
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
			err_ss << "  OrderParametersDriver::calculatePvQ: Unable to open xtc file." << "\n";
			throw std::runtime_error( err_ss.str() );
		}
	}

	// Allocate memory for xtc_atom_positions (and derivatives, if requested)
	try { 
		xtc_atom_positions.resize(num_atoms);
	}
	catch (const std::bad_alloc& ba) {
    std::cerr << "  OrderParametersDriver::calculatePvQ: Unable to allocate memory." << "\n"
		          << "  .what() = " << ba.what() << "\n";

		// Close the xtc file (don't leave the file in a bad state!)
		xdrfile_close( xtc_file_ptr );
    throw;
	}

	// Time series
	std::vector<double> t_samples, ntilde_v_samples; 
	std::vector<int>    n_v_samples;

	// For printing biasing forces for analyzed frames
	bool print_bias = false;
	std::ofstream forces_ofs, numerical_forces_ofs, bias_ofs;
	if ( (need_derivatives_ or print_forces_ or print_numerical_forces_ ) and 
	     want_output_ and my_rank_ == master_rank_ 
	) {
		// Header
		std::stringstream header_ss;
		header_ss << "# Biasing forces on OrderParameters atoms\n"
		          << "#   Units: t(ps), box_matrix(nm), virial(kJ/mol), force(kJ/mol/nm)\n"
		          << "#   Each line starts with the global serial number of the affected atom\n"
		          << "N_total= " << num_atoms << "  N_biased= " << op_atom_global_indices_.size() << "\n";

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

		if ( print_forces_ or print_numerical_forces_ ) {
			print_bias = true;
			bias_ofs.open("bias.out");
			bias_ofs << "# t[ps]  u_bias_total[kJ/mol]  u_bias_i(active_biases)\n";
		}
	}


	//----- Run Calculation -----//

	int num_frames = 0;

	// User feedback
	std::string spaces("  ");
	if ( my_rank_ == master_rank_ ) {
		std::cout << spaces << "INDUS calculation initialized." << "\n"
		          << spaces << "Beginning to read XTC file for atom positions." << "\n"
		          << spaces << "Number of atoms: " << "\t" << num_atoms << std::endl;

		// Runtime parallelization settings 
		// (also provides a quick way to check what was successfully included at compile time)
		bool is_mpi_enabled = MpiCommunicator::is_mpi_enabled();
		bool is_omp_enabled = OpenMP::is_enabled();
		if ( is_mpi_enabled or is_omp_enabled ) {
			std::cout << spaces << "Parallelization settings" << std::endl;
			if ( is_mpi_enabled ) {
				std::cout << spaces << spaces << mpi_communicator_.get_size() << " MPI ranks" << std::endl;
			}
			if ( OpenMP::is_enabled() ) {
				std::cout << spaces << spaces << OpenMP::get_max_threads() << " OpenMP threads (max)";
				if ( is_mpi_enabled ) { std::cout << " per MPI rank"; }
				std::cout << std::endl;
			}
		}
	}

	int reading_xtc_file = 1;
	int analyze_frame    = 1;

	while (true) {
		// Master rank reads the XTC file
		if ( my_rank_ == master_rank_ ) {
			read_xtc_timer_.start();
			int xdr_return_code = read_xtc(xtc_file_ptr, num_atoms,
			                               &step, &xdrfile_time, xdrfile_box, 
			                               xtc_atom_positions.data(), &xdrfile_prec);
			read_xtc_timer_.stop();

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
			this->update(xtc_atom_positions, box_matrix, step, time);

			// Calculate q (and/or Ntilde_v)
			this->calculate();

			// Get INDUS order parameters from probe volume
			const int    n_v      = indus_ptr_->get_n_v();
			const double ntilde_v = indus_ptr_->get_ntilde_v();

			// User feedback
			if ( (num_frames % 500) == 0 and my_rank_ == master_rank_ ) {
				std::cout << spaces << "Time[ps] " << time << ": Nv = " << n_v << ", Ntilde_v = " << ntilde_v << std::endl;
			}

			// Store time series of samples
			t_samples.push_back( time );
			n_v_samples.push_back( n_v );
			ntilde_v_samples.push_back( ntilde_v );


			//----- Printing Forces -----//

			// TODO Move to function for handling bias and associated forces
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

			// Bias
			if ( print_bias ) {
				bias_ofs << time << "  " << u_bias_total_;
				int num_ops = u_bias_ops_.size();
				for ( int i=0; i<num_ops; ++i ) {
					if ( bias_ptrs_[i] != nullptr ) {
						bias_ofs << "  " << u_bias_ops_[i];
					}
				}
				bias_ofs << "\n";
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
		if ( bias_ofs.is_open() ) {
			bias_ofs.close();
		}
	}

	// Finalize OPs
	int num_ops = op_ptrs_.size();
	for ( int i=0; i<num_ops; ++i ) {
		if ( op_ptrs_[i] != nullptr ) {
			op_ptrs_[i]->finalize();
		}
	}


	//----- Statistics -----//

	std::vector<double> n_v_samples_d( n_v_samples.begin(), 
	                                   n_v_samples.end() );  // convert to doubles
	double avg_n_v = Statistics::average( n_v_samples_d );
	double var_n_v = Statistics::variance( n_v_samples_d );


	//----- Write output files -----//

	if ( want_output_ and my_rank_ == master_rank_ ) {
		int alpha_kappa     = 2.0;  // multiple of kappa_0 to use for spring constant
		int delta_f_overlap = 4.0;  // Fv(q) for adjacent windows should overlap at this level [in kBT]

		// Prefactor for window spacing
		double delta_factor = 2.0*sqrt(delta_f_overlap/(alpha_kappa + 1.0));

		// Suggested biasing parameters (kappa-values in kBT)
		double kappa_0_n_v    = 1.0/var_n_v;
		double kappa_n_v      = alpha_kappa*kappa_0_n_v;
		double delta_n_v_star = delta_factor*sqrt(var_n_v);

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
		              << "#   Suggested biasing parameters (" << alpha_kappa << "*kappa_0 and "
		                << delta_f_overlap << " kBT overlap rules of thumb)\n"
		              << "#     kappa_n_v      " << kappa_n_v      << " [kBT]\n"
		              << "#     delta_n_v_star " << delta_n_v_star << "\n"
		              << "#" << "\n";

		//----- Write time series of N_v and Ntilde_v to file -----//

		std::ostringstream timeSeriesFileName;

		timeSeriesFileName << "time_samples_indus" << output_file_suffix_ << ".out";
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
#endif // #ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE


// TODO should probably move some things to SimulationState
void OrderParametersDriver::update(
		const RvecArray& atom_positions, const Box& box_matrix, const int step, const double time)
{
	update_timer_.start();

	//----- Update all atom positions -----//

#ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE
	const RvecArray& global_atom_positions = atom_positions;
	updateAtomGroups(global_atom_positions);
#else /* PLUMED MODE */
	const RvecArray& op_atom_positions = atom_positions;
	updateAtomGroupsPlumed(op_atom_positions);
#endif /* ORDER_PARAMETERS_DRIVER_PLUMED_MODE */


	//----- Update other state variables -----//

	// Create local copies that can be changed in standalone mode
	Box    local_box_matrix = box_matrix;
	int    local_step = step;
	double local_time = time;

#ifndef ORDER_PARAMETERS_DRIVER_PLUMED_MODE
	if ( mpi_communicator_.is_mpi_initialized() ) {
		// Share box
		mpi_communicator_.bcast(local_box_matrix, master_rank_);

		// Share step and time
		mpi_communicator_.bcast(local_step, master_rank_);
		mpi_communicator_.bcast(local_time, master_rank_);
	}
#endif /* ORDER_PARAMETERS_DRIVER_PLUMED_MODE */ 

	// Update simulation box (assumes orthorhombic box)
	SimulationBox& simulation_box = simulation_state_.access_simulation_box();
	simulation_box.setLengths( local_box_matrix[X_DIM][X_DIM], 
	                           local_box_matrix[Y_DIM][Y_DIM], 
	                           local_box_matrix[Z_DIM][Z_DIM] );

	simulation_state_.set_step( local_step );
	simulation_state_.set_time( local_time );

	// Ensure that atoms are fully within the simulation box
	// - TODO: Move elsewhere
	std::vector<AtomGroup>& atom_groups = simulation_state_.access_atom_groups();
	int num_groups = atom_groups.size();
	for ( int g=0; g<num_groups; ++g ) {
		auto& group_atom_positions = atom_groups[g].access_atom_positions();
		int num_atoms_in_group = group_atom_positions.size();
		for ( int k=0; k<num_atoms_in_group; ++k ) {
			simulation_box.putInBox( group_atom_positions[k] );
		}
	}

	//----- Update Domain Decomposition -----//

	// TODO update the number of cells in each direction based
	//      on box dimensions?
	// - Put the most cells along the longest box dimension

	// Update the ProbeVolume's bounding box
	probe_volume_ptr_->update();
	probe_volume_ptr_->setBoundingBox();  // TODO move to ProbeVolume::update()/setGeometry()?

	// Set DD region based on the ProbeVolume
	dd_cell_grid_.setCellGridRegion( probe_volume_ptr_->getBoundingBox() );
	//dd_cell_grid_.setCellGridRegion(); // sets DD region to be equal to SimulationBox

	// Find the atoms which are local to this rank
	// FIXME check whether group requires DD and/or which style
	// TODO move elsewhere?
	update_atom_groups_timer_.start();
	for ( int i=0; i<num_groups; ++i ) {
		dd_cell_grid_.findNearbyAtoms(my_rank_, atom_groups[i]);
	}
	update_atom_groups_timer_.stop();

	// Now that AtomGroups have been updated, perform updates on the OrderParameters
	int num_ops = op_ptrs_.size();
	for ( int n=0; n<num_ops; ++n ) {
		op_ptrs_[n]->update();
	}

	/*
	// DEBUG: Uncomment the following lines for debugging domain decomposition
	if ( my_rank_ == master_rank_ ) {
		std::cout << dd_cell_grid_ << std::endl;
	}
	dd_cell_grid_.printAtomGroupDomainDecompositionCheck(std::cout, atom_groups[0]);
	*/

	update_timer_.stop();
}


void OrderParametersDriver::updateAtomGroups(const RvecArray& global_atom_positions)
{
	std::vector<AtomGroup>& atom_groups = simulation_state_.access_atom_groups();
	int num_groups = atom_groups.size();
	for ( int i=0; i<num_groups; ++i ) {
		// Reserve space
		std::vector<Real3>& group_atom_positions = atom_groups[i].access_atom_positions();
		int num_atoms_in_group = atom_groups[i].get_num_atoms();
		group_atom_positions.resize(num_atoms_in_group);

		// Copy positions of group atoms from array of all positions
		if ( my_rank_ == master_rank_ ) {
			const std::vector<int>& global_atom_indices = atom_groups[i].get_global_atom_indices();
			for ( int k=0; k<num_atoms_in_group; ++k ) {
				int global_index = global_atom_indices[k];
				for ( int d=0; d<DIM_; ++d ) {
					group_atom_positions[k][d] = global_atom_positions[global_index][d];
				}
			}
		}

		// Share positions as needed
		if ( mpi_communicator_.is_mpi_initialized() ) {
			mpi_communicator_.bcast(group_atom_positions, master_rank_);
		}
	}
}


void OrderParametersDriver::updateAtomGroupsPlumed(const RvecArray& op_atom_positions)
{
	// When compiled as part of PLUMED, each rank already has the positions of 
	// all OP atoms (no further communication is necessary)
	std::vector<AtomGroup>& atom_groups = simulation_state_.access_atom_groups();
	int num_groups = atom_groups.size();
	for ( int i=0; i<num_groups; ++i ) {
		// Reserve space
		std::vector<Real3>& group_atom_positions = atom_groups[i].access_atom_positions();
		int num_atoms_in_group = atom_groups[i].get_num_atoms();
		group_atom_positions.resize(num_atoms_in_group);
	}

	int op_atom_offset = 0;
	for ( int i=0; i<num_groups; ++i ) {
		std::vector<Real3>& group_atom_positions = atom_groups[i].access_atom_positions();
		int num_atoms_in_group = atom_groups[i].get_num_atoms();

		// Copy positions to AtomGroups
		for ( int k=0; k<num_atoms_in_group; ++k ) {
			for ( int d=0; d<DIM_; ++d ) {
				group_atom_positions[k][d] = op_atom_positions[op_atom_offset + k][d];
			}
		}
		op_atom_offset += num_atoms_in_group;
	}
}


void OrderParametersDriver::calculate()
{
	calculate_timer_.start();

	if ( need_derivatives_ ) {
		derivatives_are_ready_ = false;
	}

	// *** INDUS must go before Steinhardt ***
	// - That way, the probe volume indicator functions that Steinhardt needs
	//   will be available in INDUS
	int num_ops = op_ptrs_.size();
	op_calculate_timer_.start();
	for ( int n=0; n<num_ops; ++n ) {
		op_ptrs_[n]->calculate();
	}
	op_calculate_timer_.stop();

	// TODO Synchronize all at once (concentrate global communication here)
	op_synchronize_timer_.start();
	for ( int n=0; n<num_ops; ++n ) {
		op_ptrs_[n]->synchronize();
	}
	op_synchronize_timer_.stop();

	if ( need_derivatives_ ) {
		// Gather local derivatives of the OPs
		this->collectLocalOrderParameterDerivatives();
	}

	// Compute biasing potentials and their derivatives 
	// (the negatives of the biasing forces)
	this->calculateTotalBias();

	// Optionally, share all derivatives for all OPs
	// - TODO allow this on a per-OP basis
	if ( need_derivatives_ and share_all_op_derivatives_ ) {
		shareAllOrderParameterDerivatives();
	}

	if ( need_derivatives_ ) {
		derivatives_are_ready_ = true;
	}

	calculate_timer_.stop();

	return;
}


void OrderParametersDriver::collectLocalOrderParameterDerivatives()
{
	collect_local_derivs_timer_.start();

	int num_ops = op_ptrs_.size();
	for ( int n=0; n<num_ops; ++n ) {
		// Indices of groups involved
		const std::vector<int>& atom_group_indices = op_ptrs_[n]->get_atom_group_indices();
		int num_atom_groups = atom_group_indices.size();

		// Access derivatives and indices
		const std::vector<std::vector<Real3>>& local_group_derivatives = op_ptrs_[n]->get_local_derivatives();
		const std::vector<std::vector<int>>&   local_group_indices = op_ptrs_[n]->get_local_group_indices();

		// Put together derivatives across all AtomGroups involved
		local_derivatives_ops_[n].resize(0);
		local_op_indices_for_derivatives_[n].resize(0);
		for ( int k=0; k<num_atom_groups; ++k ) {
			int group_index = atom_group_indices[k];

			local_derivatives_ops_[n].insert( 
				local_derivatives_ops_[n].end(),
				local_group_derivatives[k].begin(), local_group_derivatives[k].end() );

			// Map from group indices to OP indices
			int old_size        = local_op_indices_for_derivatives_[n].size();
			int num_new_indices = local_group_indices[k].size();
			local_op_indices_for_derivatives_[n].resize( old_size + num_new_indices );
			for ( int i=0; i<num_new_indices; ++i ) {
				local_op_indices_for_derivatives_[n][old_size + i] = 
					atom_group_to_op_indices_[group_index][ local_group_indices[k][i] ];
			}
		}

		sum_r_cross_op_derivatives_[n] = op_ptrs_[n]->get_sum_r_cross_derivatives();
	}

	collect_local_derivs_timer_.stop();
}


void OrderParametersDriver::shareAllOrderParameterDerivatives()
{
	share_all_derivs_timer_.start();

	int num_ops = op_ptrs_.size();
	if ( static_cast<int>(derivatives_ops_.size()) != num_ops ) {
		derivatives_ops_.resize(num_ops);
	}

	for ( int n=0; n<num_ops; ++n ) {
		if ( mpi_communicator_.is_mpi_initialized() and mpi_communicator_.get_size() > 1 ) {
			// Share local derivatives and corresponding group indices across all ranks
			mpi_communicator_.allgatherv( 
					local_derivatives_ops_[n], 
					// Output
					derivatives_buffer_.buffer, derivatives_buffer_.block_offsets,
					derivatives_buffer_.block_sizes );
			mpi_communicator_.allgatherv( 
					local_op_indices_for_derivatives_[n],
					// Output
					group_indices_buffer_.buffer, group_indices_buffer_.block_offsets,
					group_indices_buffer_.block_sizes );
		}
		else {
			int num_ranks = 1;
			int num_arrays = local_derivatives_ops_[n].size();  // numbers of Real3 elements
			derivatives_buffer_.buffer        = local_derivatives_ops_[n];
			derivatives_buffer_.block_offsets = std::vector<int>(num_ranks, 0);
			derivatives_buffer_.block_sizes   = std::vector<int>(num_ranks, num_arrays);
			group_indices_buffer_.buffer        = local_op_indices_for_derivatives_[n];
			group_indices_buffer_.block_offsets = std::vector<int>(num_ranks, 0);
			group_indices_buffer_.block_sizes   = std::vector<int>(num_ranks, num_arrays);
		}

		// Reset array
		int num_op_atoms = op_atom_global_indices_.size();
		derivatives_ops_[n].assign(num_op_atoms, CommonTypes::zero_3);

		// Add together all derivatives
		// TODO OPENMP
		int num_derivatives = derivatives_buffer_.buffer.size();
		int op_index;
		for ( int i=0; i<num_derivatives; ++i ) {
			op_index = group_indices_buffer_.buffer[i];

			for ( int d=0; d<DIM_; ++d ) {
				derivatives_ops_[n][op_index][d] += derivatives_buffer_.buffer[i][d];
			}
		}
	}

	share_all_derivs_timer_.stop();
}


void OrderParametersDriver::calculateTotalBias()
{
	bias_timer_.start();

	// Prepare output variables
	int num_ops = op_ptrs_.size();
	u_bias_ops_.assign(num_ops, 0.0);
	u_bias_total_ = 0.0;
	if ( need_derivatives_ ) {
		local_derivatives_u_bias_total_.resize(0);
		local_op_indices_u_bias_total_.resize(0);

		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				for ( int n=0; n<num_ops; ++n ) {
					virial_u_bias_ops_[n][a][b] = 0.0;
				}
				virial_u_bias_total_[a][b] = 0.0;
			}
		}
	}

	for ( int n=0; n<num_ops; ++n ) {
		if ( bias_ptrs_[n] != nullptr ) {
			const double value = op_ptrs_[n]->get_value();

			bias_ptrs_[n]->applyBias(
					value, need_derivatives_, local_derivatives_ops_[n], sum_r_cross_op_derivatives_[n],
					// Output
					u_bias_ops_[n], local_derivatives_u_bias_ops_[n], virial_u_bias_ops_[n]
			);
			u_bias_total_ += u_bias_ops_[n];

			if ( need_derivatives_ ) {
				// Combine arrays across multiple OPs for ease of handling
				local_derivatives_u_bias_total_.insert( 
						local_derivatives_u_bias_total_.end(),
						local_derivatives_u_bias_ops_[n].begin(), local_derivatives_u_bias_ops_[n].end() );
				local_op_indices_u_bias_total_.insert( 
						local_op_indices_u_bias_total_.end(),
						local_op_indices_for_derivatives_[n].begin(), local_op_indices_for_derivatives_[n].end() );

				for ( int a=0; a<DIM_; ++a ) {
					for ( int b=0; b<DIM_; ++b ) {
						virial_u_bias_total_[a][b] += virial_u_bias_ops_[n][a][b];
					}
				}
			}
		}
	}

	// Share derivatives across ranks
	if ( need_derivatives_ ) {
		if ( mpi_communicator_.is_mpi_initialized() and mpi_communicator_.get_size() > 1 ) {
			// Share local derivatives and corresponding group indices across all ranks
			// TODO combine into one call with a custom datatype?
			bias_gather_timer_.start();
			mpi_communicator_.allgatherv( 
					local_derivatives_u_bias_total_, 
					// Output
					derivatives_buffer_.buffer, derivatives_buffer_.block_offsets,
					derivatives_buffer_.block_sizes );
			mpi_communicator_.allgatherv( 
					local_op_indices_u_bias_total_,
					// Output
					group_indices_buffer_.buffer, group_indices_buffer_.block_offsets,
					group_indices_buffer_.block_sizes );
			bias_gather_timer_.stop();
		}
		else {
			int num_ranks = 1;
			int num_arrays = local_derivatives_u_bias_total_.size();  // numbers of Real3 elements
			derivatives_buffer_.buffer        = local_derivatives_u_bias_total_;
			derivatives_buffer_.block_offsets = std::vector<int>(num_ranks, 0);
			derivatives_buffer_.block_sizes   = std::vector<int>(num_ranks, num_arrays);
			group_indices_buffer_.buffer        = local_op_indices_u_bias_total_;
			group_indices_buffer_.block_offsets = std::vector<int>(num_ranks, 0);
			group_indices_buffer_.block_sizes   = std::vector<int>(num_ranks, num_arrays);
		}

		// Add together all derivatives
		int num_op_atoms = op_atom_global_indices_.size();
		derivatives_u_bias_total_.assign(num_op_atoms, CommonTypes::zero_3);
		int num_derivatives = derivatives_buffer_.buffer.size();
		int op_index;
		bias_sum_timer_.start();
		for ( int i=0; i<num_derivatives; ++i ) {
			op_index = group_indices_buffer_.buffer[i];

			for ( int d=0; d<DIM_; ++d ) {
				derivatives_u_bias_total_[op_index][d] += derivatives_buffer_.buffer[i][d];
			}
		}
		bias_sum_timer_.stop();

		// TODO Optional to also share all local derivatives for individual OPs across all ranks
	}

	bias_timer_.stop();
}


// Returns a big string with a whole bunch of input variables.
std::string OrderParametersDriver::getInputSummary() const 
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
	   << "#\n";

	// Target group selection
	int target_group_index = indus_ptr_->get_target_atom_group_index();
	const AtomGroup& target_group = simulation_state_.get_atom_group( target_group_index );
	ss << "# Target_Atoms_String = " << target_group.get_selection_string() << "\n";

	ss << "# Gro_File = " << gro_file_ << "\n"
		 << "# Xtc_File = " << xtc_file_ << "\n"
	   << "#   Production_phase\n"  
		 << "#     Start = " << t_min_ << " [ps]\n"   // TODO if at numeric_limits, print "none/unbounded"
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


void OrderParametersDriver::printForcesForFrame(
		std::ofstream& forces_ofs, const double time, const Matrix& box_matrix, 
		const std::vector<Real3>& derivatives_u_bias_total, const Matrix& virial_u_bias_total
) const 
{
	print_forces_timer_.start();

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
	int num_op_atoms = op_atom_global_indices_.size();
	for ( int k=0; k<num_op_atoms; ++k ) {
		forces_ofs << op_atom_global_indices_[k] + 1;
		//for ( int d=0; d<DIM_; ++d ) { forces_ofs << "  " << xtc_atom_positions[k][d]; }
		for ( int d=0; d<DIM_; ++d ) { 
			forces_ofs << "  " << -derivatives_u_bias_total[k][d]; 
		}
		forces_ofs << "\n";
	}

	print_forces_timer_.stop();
}


void OrderParametersDriver::computeNumericalDerivatives()
{
	// Save initial values of key variables
	const double ntilde_v_0        = indus_ptr_->get_ntilde_v();
	const double u_bias_total_0    = u_bias_total_;
	//std::vector<double> u_bias_ops_0 = u_bias_ops_;

	// Allocate memory and set everything to zero
	int num_op_atoms = op_atom_global_indices_.size();
	numerical_derivatives_ntilde_v_.assign( num_op_atoms, CommonTypes::zero_3 );
	numerical_derivatives_u_bias_total_.assign( num_op_atoms, CommonTypes::zero_3 );
	for ( int a=0; a<DIM_; ++a ) {
		numerical_virial_u_bias_total_[a] = CommonTypes::zero_3;
	}

	if ( simulation_state_.debug_mode() ) {
		std::cout << "NUMERICAL DERIVATIVES (" << LOCATION_IN_SOURCE_STRING << ")" << std::endl;
	}

	// Manually perturb each atom position and compute the resulting changes
	// to estimate the derivatives
	std::vector<AtomGroup>& atom_groups = simulation_state_.access_atom_groups();
	int num_atom_groups = atom_groups.size();
	for ( int n=0; n<num_atom_groups; ++n ) {
		if ( simulation_state_.debug_mode() ) {
			std::cout << "  GROUP \"" << atom_groups[n].get_name() << "\" "
			            << "(" << LOCATION_IN_SOURCE_STRING << ")" << std::endl;
		}

		// Access positions array
		std::vector<Real3>& positions = atom_groups[n].access_atom_positions();
		int num_atoms = positions.size();

		for ( int i=0; i<num_atoms; ++i ) {
			// Save actual position
			Real3 x_i_0 = positions[i];

			try {
				int op_atom_index = atom_group_to_op_indices_[n][i];

				for ( int d=0; d<DIM_; ++d ) {
					// Perturb position
					positions[i][d] += delta_x_;

					// Recompute everything
					this->calculate();
					const double ntilde_v_new = indus_ptr_->get_ntilde_v();

					// Accumulate numerical derivatives
					numerical_derivatives_ntilde_v_[op_atom_index][d] += 
							(ntilde_v_new - ntilde_v_0)/delta_x_;
					numerical_derivatives_u_bias_total_[op_atom_index][d] += 
							(u_bias_total_ - u_bias_total_0)/delta_x_;
					for ( int a=0; a<DIM_; ++a ) {
						numerical_virial_u_bias_total_[a][d] += 
								0.5*x_i_0[a]*numerical_derivatives_u_bias_total_[op_atom_index][d];
					}

					// Put atom back where it started
					positions[i][d] = x_i_0[d];
				}
			}
			catch (const std::exception& ex) {
				// Replace original position
				positions[i] = x_i_0;

				// Rethrow
				throw;
			}
		}
	}
}


void OrderParametersDriver::finalize()
{
	printProfilingInfo();
}


void OrderParametersDriver::printProfilingInfo()
{
	// TODO custom file names based on input option

	// Each rank prints a file with suffix ".<myrank>"
	GPTL::print("gptl.log", mpi_communicator_);  // 

	// Summary across all ranks and threads
	GPTL::printSummary("gptl_summary.log", mpi_communicator_);
}
