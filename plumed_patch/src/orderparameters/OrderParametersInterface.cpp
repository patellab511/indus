#include "OrderParametersInterface.h"

namespace PLMD {
namespace colvar {


// Register this action with the Plumed core
PLUMED_REGISTER_ACTION(OrderParametersInterface,"INDUS")


const std::map<int, std::string>
OrderParametersInterface::axis_index_to_label_map_ = {
	{ X_DIM, "x" },
	{ Y_DIM, "y" },
	{ Z_DIM, "z" }
};


void OrderParametersInterface::registerKeywords( Keywords& keys )
{
	// Register the keywords that this action inherits as a derived class of Colvar
  Colvar::registerKeywords(keys);

	// Input keys
  keys.add("compulsory", "INPUTFILE", "order parameters input file");
  keys.add("optional",   "TOL",       "tolerance on the largest norm expected for a derivative (for debugging)");


	//----- Output keys -----//

	// Collective variables
	keys.addOutputComponent("n", "default",
	                        "number of target particles in probe volume");
	keys.addOutputComponent("ntilde", "default",
	                        "coarse-grained number of target particles in probe volume");

	// Biasing potentials
	keys.addOutputComponent("ubias", "default",
	                        "total biasing potential");
	keys.addOutputComponent("ubiasntilde", "default",
	                        "bias on ntilde");

	// Tracking large forces and edge cases
	keys.addOutputComponent("maxnormderiv", "default",
	                        "largest norm across all derivatives (for debugging)");

	// Pressure due to biasing potential
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			std::string p_ab_name = "P" + axis_index_to_label_map_.at(a) + axis_index_to_label_map_.at(b);
			keys.addOutputComponent(p_ab_name, "default",
			                        "pressure [bar] due to biasing potential: component " + p_ab_name);
		}
	}
}


// Constructor
OrderParametersInterface::OrderParametersInterface(const ActionOptions& ao):
	// This line sets up various things in the plumed core which colvars rely on.
	PLUMED_COLVAR_INIT(ao),
	input_file_(""),
#ifdef MPI_ENABLED
	mpi_communicator_( comm.Get_comm() ),
#else
	// Don't try to copy dummy object PLMD::MPI_Comm - use the default constructor instead
	mpi_communicator_(),
#endif
	my_rank_( mpi_communicator_.get_rank() ),
	plumed_simulation_box_( getBox() ),
	max_norm_deriv_tol_(-1.0),
	check_max_norm_deriv_(false)
{
	//----- Parse input -----//

	// Get name of external input file
	parse("INPUTFILE", input_file_);

	// Tolerance on the maximum derivative expected
	parse("TOL", max_norm_deriv_tol_);

	// Check that everything on the input line has been read properly
  checkRead();


	//----- Initialize external objects necessary for calculation -----//

	// OrderParametersDriver object will process its input file
	try {
		driver_ptr_ = std::unique_ptr<OrderParametersDriver>( 
			new OrderParametersDriver(input_file_, mpi_communicator_)
		);

		// Ensure that all derivatives are available on all ranks in case
		// Ntilde is passed directly to another Action like RESTRAINT
		// - TODO: Automatic detection?
		driver_ptr_->set_share_all_derivatives(true);
	}
	catch (const std::exception& e) {
		std::cerr << "OrderParametersInterface: Unable to allocate memory for external "
		          << "objects in constructor.\n";
		std::cerr << "  .what(): " << e.what() << "\n";
		throw;
	}

	// Global indices of targets are parsed by OrderParametersDriver
	// - Convert to PLUMED format
	const std::vector<int>& op_atom_global_indices = driver_ptr_->get_op_atom_global_indices();
	int num_op_atoms = op_atom_global_indices.size();
	op_atom_numbers_.resize(num_op_atoms);
	for ( int i=0; i<num_op_atoms; ++i ) {
		op_atom_numbers_[i].setIndex( op_atom_global_indices[i] );
	}


	//----- Handle PLUMED components -----//

	// Tell the plumed core that we require space to store the value of CV components, that they're
	// not periodic, and that the CV will act on a particular list of atoms.
	// - These functions are inherited from PLMD::ActionWithValue <-- PLMD::Colvar
	if ( driver_ptr_->share_all_derivatives() ) {
		addComponentWithDerivatives("q");
		addComponentWithDerivatives("ntilde");
	}
	else {
		addComponent("q");
		addComponent("ntilde");
	}
	componentIsNotPeriodic("q");
	componentIsNotPeriodic("ntilde");

	addComponentWithDerivatives("ubias");
	componentIsNotPeriodic("ubias");

	addComponent("n");                componentIsNotPeriodic("n");
	addComponent("ubiasntilde");      componentIsNotPeriodic("ubiasntilde");
	addComponent("maxnormderiv");     componentIsNotPeriodic("maxnormderiv");

	// Pressure due to the bias
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			std::string p_ab_name = "P" + axis_index_to_label_map_.at(a) + axis_index_to_label_map_.at(b);
			addComponent(p_ab_name);  componentIsNotPeriodic(p_ab_name);
		}
	}

	//----- Log file -----//

	log.printf("  with external input file %s\n", input_file_.c_str());

	// Record target atoms
	log.printf("  with %d target atoms\n", num_op_atoms);
  log.printf("  list of target atoms (by serial):");
  for (int i=0; i<num_op_atoms; ++i) {
    if ( i % 25 == 0 ) { log << "\n"; }
    log << " " << op_atom_numbers_[i].serial();
  }
  log.printf("\n");

	check_max_norm_deriv_ = ( max_norm_deriv_tol_ >= 0.0 );
	if ( check_max_norm_deriv_ ) {
		log << "  with largest expected norm of a derivative set at tol="
		    << max_norm_deriv_tol_ << "\n";
		log << "    PLUMED will now check for abnormally large derivatives in this action.\n";
	}

	// Parallelization support (TODO: print more detailed info/move to OPDriver function)
	bool is_mpi_enabled        = MpiCommunicator::is_mpi_enabled();
	bool is_openmp_enabled     = OpenMP::is_enabled();
	bool share_all_derivatives = driver_ptr_->share_all_derivatives();
	log << "  Parallelization\n"
	    << "    is_mpi_enabled        = " << is_mpi_enabled        << "\n"
	    << "    is_openmp_enabled     = " << is_openmp_enabled     << "\n"
	    << "    share_all_derivatives = " << share_all_derivatives << "\n";

	// Summary of everything else that was read by Steinhardt
	log << "  BEGIN STEINHARDT INPUT SUMMARY\n"
	    << driver_ptr_->getInputSummary()
	    << "  END STEINHARDT INPUT SUMMARY\n";


	//----- Request Atoms -----//

	// This step must be done after all components have been added with/without derivatives.
	// Otherwise, the necessary arrays for derivatives will not be allocated.
	try {
		requestAtoms(op_atom_numbers_);
	}
	catch (const std::exception& e) {
		std::cerr << "OrderParametersInterface: Unable to request atoms.\n"
		          << "  .what(): " << e.what() << "\n";
		throw;
	}
}


void OrderParametersInterface::calculate()
{
	//--------------------------------------//
	//----- Calculate Order Parameters -----//
	//--------------------------------------//

	// Get variables from the PLUMED core
	double time = getTime();
	int    step = getStep();
	const std::vector<PLMD::Vector>& op_atom_positions = getPositions();
	int num_op_atoms = op_atom_numbers_.size();
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			box_matrix_[a][b] = plumed_simulation_box_[a][b];
		}
	}

	try {
		// Update simulation state
		update_timer_.start();
		driver_ptr_->update(op_atom_positions, box_matrix_, step, time);
		update_timer_.stop();

		// Calculate everything
		calculate_timer_.start();
		driver_ptr_->calculate();
		calculate_timer_.stop();
	}
	catch ( const std::exception& ex ) {
		// Dump coordinates and derivatives for debugging

		std::stringstream err_ss;
		log << "Exception at " << __PRETTY_FUNCTION__
		      << " (" << __FILE__ << ":" << __LINE__ << ")\n"
		    << "  Printing frame at t = " << time << " ps (step " << step << ") for debugging.\n";

		std::stringstream ss;
		ss << "frame_at_exception_" << time;
		std::string file_name = ss.str();
		printFrameWithDerivatives(file_name, "Frame at exception", FrameFileType::xyz);

		// Rethrow
		throw;
	}

	component_output_timer_.start();


	//-------------------------------//
	//----- Output for "ntilde" -----//
	//-------------------------------//

	// Access OrderParameters output variables
	ntilde_v_ = driver_ptr_->get_ntilde_v();
	const std::vector<Real3>& derivatives_ntilde_v = driver_ptr_->get_derivatives_ntilde_v();
	const Matrix& sum_r_cross_dntilde_v_dr         = driver_ptr_->get_sum_r_cross_dntilde_v_dr();

	// "Box derivatives"
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			box_derivatives_ntilde_v_[a][b] = -sum_r_cross_dntilde_v_dr[a][b];
		}
	}

	// Pass output to PLUMED core
	Value* ntilde_ptr = getPntrToComponent("ntilde");

	ntilde_ptr->set(ntilde_v_);
	if ( ntilde_ptr->hasDerivatives() ) {
		for ( int i=0; i<num_op_atoms; i++ ) {
			// PLUMED stores derivatives in a linear array
			int atom_offset = DIM_*i;
			for ( int d=0; d<DIM_; ++d ) {
				ntilde_ptr->addDerivative(atom_offset + d, derivatives_ntilde_v[i][d]);
			}
		}
		setBoxDerivatives(ntilde_ptr, box_derivatives_ntilde_v_);
	}


	//------------------------------//
	//----- Output for "ubias" -----//
	//------------------------------//

	// Access OrderParameters output variables
	u_bias_ = driver_ptr_->get_u_bias_total();
	const std::vector<Real3>& derivatives_u_bias_total
			= driver_ptr_->get_derivatives_u_bias_total();
	const Matrix& virial_u_bias_total
			= driver_ptr_->get_virial_u_bias_total();

	// "Box derivatives"
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			box_derivatives_u_bias_[a][b] = -2.0*virial_u_bias_total[a][b];
		}
	}

	// Pass output to PLUMED core
	Value* ubias_ptr = getPntrToComponent("ubias");
	ubias_ptr->set(u_bias_);
	if ( ubias_ptr->hasDerivatives() ) {
		for ( int i=0; i<num_op_atoms; i++ ) {
			// PLUMED stores derivatives in a linear array
			int atom_offset = DIM_*i;
			for ( int d=0; d<DIM_; ++d ) {
				ubias_ptr->addDerivative(atom_offset + d, derivatives_u_bias_total[i][d]);
			}
		}
		setBoxDerivatives(ubias_ptr, box_derivatives_u_bias_);
	}

	// Convert pressure from kJ/mol/nm^3 to bar
	constexpr double n_Av      = 6.022e23;  // Avogadro's number (1/mol)
	constexpr double convert_p = 1.0e25/n_Av;

	// Pressure due to the bias
	double volume           = driver_ptr_->get_simulation_box().getVolume();
	double virial_prefactor = -2.0/volume;
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			double p_bias_ab = virial_prefactor*virial_u_bias_total[a][b];
			std::string p_ab_name = "P" + axis_index_to_label_map_.at(a) + axis_index_to_label_map_.at(b);
			getPntrToComponent(p_ab_name)->set(convert_p*p_bias_ab);
		}
	}


	//------------------------------------------//
	//----- Checking for large derivatives -----//
	//------------------------------------------//

	if ( check_max_norm_deriv_ and ntilde_ptr->hasDerivatives() ) {
		this->calculateMaxNormDeriv( derivatives_ntilde_v,
		                             max_norm_deriv_, global_index_max_norm_deriv_ );
	}

	// Set output value
	getPntrToComponent("maxnormderiv")->set(max_norm_deriv_);


	//---------------------------------------------//
	//----- Handle output for everything else -----//
	//---------------------------------------------//

	// Integer N_v
	n_v_ = driver_ptr_->get_n_v();
	getPntrToComponent("n")->set(n_v_);

	// Bias on ntilde
	u_bias_ntilde_v_ = driver_ptr_->get_u_bias_ntilde_v();
	getPntrToComponent("ubiasntilde")->set(u_bias_ntilde_v_);

	component_output_timer_.stop();


	//---------------------------------------------//
	//----- Debugging and checking for errors -----//
	//---------------------------------------------//

	// Only check if max_norm_deriv_tol_ is non-negative
	// - The constructor gives it a negative sentinel value so as to not bother checking/printing
	//   if the user never specified a tolerance
	if ( max_norm_deriv_tol_ >= 0.0 and max_norm_deriv_ > max_norm_deriv_tol_ and
	     my_rank_ == 0
	) {
			// Dump coordinates and derivatives for debugging

			// File name
			double time = getTime();
			std::stringstream ss;
			ss << "frame_" << time;
			std::string file_name = ss.str();

			// Message to print at top of file
			ss.str(""); ss.clear();
			ss << "max[norm(deriv)]=" << max_norm_deriv_ << " on atom " << global_index_max_norm_deriv_+1;
			std::string message = ss.str();

			printFrameWithDerivatives(file_name, message, FrameFileType::xyz);
	}
}


void OrderParametersInterface::calculateMaxNormDeriv(
		const std::vector<Real3>& derivatives,
		double& max_norm_deriv, int& global_index_max_norm_deriv) const
{
	// Dummy values
	max_norm_deriv = -1.0;
	global_index_max_norm_deriv = -1;

	// Get derivatives of q (and corresponding global indices of affected atoms)
	const int num_indus_atoms = derivatives.size();

	double max_norm_sq = 0.0, norm_sq, deriv;
	int    op_index = 0;
	try {
		for ( int u=0; u<num_indus_atoms; u++ ) {
			norm_sq = 0.0;
			for ( int d=0; d<DIM_; ++d ) {
				deriv = derivatives[u][d]; // extract for readability
				plumed_massert( std::isfinite(deriv), "A derivative is not finite.\n" );
				norm_sq += deriv*deriv;
			}

			if ( norm_sq > max_norm_sq ) {
				max_norm_sq = norm_sq;
				op_index    = u;
			}
		}

		// Maximum found
		max_norm_deriv = sqrt(max_norm_sq);

		// Convert from Indus index to global index
		const std::vector<int>& global_indices = driver_ptr_->get_op_atom_global_indices();
		global_index_max_norm_deriv = global_indices[op_index];
	}
	catch (const std::exception ex) {
		// Dump coordinates and derivatives for debugging
		std::string file_name = "frame_at_crash";
		std::string message( ex.what() );
		if ( my_rank_ == 0 ) {
			printFrameWithDerivatives(file_name, message, FrameFileType::xyz);
		}

		// Rethrow
		throw;
	}
}


// Print files with atom coordinates and order parameter derivatives
void OrderParametersInterface::printFrameWithDerivatives(
		const std::string& file_name, const std::string& message, const FrameFileType& file_type
) const
{
	double time = getTime();

	std::string atom_name = "ATOM";
	std::string resname = "RES";
	int resnum = 1;
	PLMD::Vector dummy_velocity(0.0,0.0,0.0);

	// Derivatives of the total biasing potential (negatives of the total biasing force on each atom)
	const std::vector<CommonTypes::Real3>& derivatives_u_bias_total
			= driver_ptr_->get_derivatives_u_bias_total();
	const std::vector<int>& op_atom_global_indices
			= driver_ptr_->get_op_atom_global_indices();
	int num_op_atoms = op_atom_global_indices.size();

	// Don't try to print derivatives if they aren't ready (yet)
	bool are_derivatives_ready = driver_ptr_->areDerivativesReady();

	const std::vector<PLMD::Vector>& op_atom_positions = getPositions();

	// File
	std::string file_name_with_ext = file_name;
	if ( file_type == FrameFileType::gro ) {
		file_name_with_ext += ".gro";
	}
	else if ( file_type == FrameFileType::xyz ) {
		file_name_with_ext += ".xyz";
	}
	std::ofstream ofs( file_name_with_ext );

	// Header
	ofs << "PLUMED OrderParameters at t= " << time << " ps with derivs of u_bias_total "
	      << "(message: " << message << ")\n"
	    << num_op_atoms << "\n";

	if ( file_type == FrameFileType::gro ) {
		// All fixed-point numbers
		ofs << std::fixed;
	}

	// Notes
	// - GRO format: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"

	// Print derivatives for OP atoms
	for ( int i=0; i<num_op_atoms; ++i ) {
		// Residue
		if ( file_type == FrameFileType::gro ) {
			// Number (%5d)
			ofs << std::setw(5) << resnum;

			// Name (%-5s, left-justified)
			ofs << std::left << std::setw(5) << resname;
			ofs << std::right;
		}

		// Atom name (%5s)
		ofs << std::setw(5) << atom_name;
	
		// Atom number (%5d)
		ofs << std::setw(7) << op_atom_numbers_[i].serial();

		// Position
		if ( file_type == FrameFileType::gro ) {
			// %8.3f each
			for ( int d=0; d<DIM_; ++d ) {
				ofs << std::setw(8) << std::setprecision(3) << op_atom_positions[i][d];
			}
		}
		else if ( file_type == FrameFileType::xyz ) {
			// Extended precision
			for ( int d=0; d<DIM_; ++d ) {
				ofs << std::setw(15) << std::setprecision(10) << op_atom_positions[i][d];
			}
		}

		// Velocity
		if ( file_type == FrameFileType::gro ) {
			// Velocities in nm/ps (%8.4f for each)
			for ( int d=0; d<DIM_; ++d ) {
				ofs << std::setw(8) << std::setprecision(4) << dummy_velocity[d];
			}
		}

		// Derivatives
		if ( are_derivatives_ready ) {
			// Derivatives of total bias wrt. the atom's position
			for ( int d=0; d<DIM_; ++d ) {
				ofs << std::setw(1) << " ";
				ofs << std::setw(10) << std::setprecision(6) << derivatives_u_bias_total[i][d];
			}
		}

		ofs << "\n";
	}

	// End with box lengths
	for ( int d=0; d<DIM_; ++d ) {
		ofs << "  " << plumed_simulation_box_(d,d);
	}
	ofs << "\n";

	ofs.close();
}


void OrderParametersInterface::runFinalJobs() {
	log << "DEBUG: OrderParametersInterface::runFinalJobs\n";
}


} // namespace colvar
} // namespace PLMD
