
#include "OutputHandler.h"

OutputHandler::OutputHandler( 
	const SimulationState& simulation_state, MpiCommunicator& mpi_communicator,
	const bool debug_mode
):
	simulation_state_(simulation_state),
	mpi_communicator_(mpi_communicator),
	debug_mode_(debug_mode)
{
}


void OutputHandler::writeIndicesForFrame(
	std::ofstream& ofs, const std::vector<int>& indices, const int num_atoms_per_index
)
{
	// First token is the time (in ps)
	ofs << simulation_state_.get_time();

	int num_indices = indices.size();
	for ( int i=0; i<num_indices; ++i ) {
		for ( int j=0; j<num_atoms_per_index; ++j ) {
			ofs << " " << indices[i] + j + 1;
		}
	}
	ofs << "\n";
}


void OutputHandler::make_PvN(
	const std::vector<double>& n_values, const int max_n, 
	std::vector<double>& bins, std::vector<double>& p_v_n,
  double& avg_n, double& var_n
) const
{
	if ( max_n < 0 ) {
		std::stringstream err_ss;
		err_ss << "OrderParametersDriver::make_PvN: ERROR: max(N) can't be less than zero (input: max_n="
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

	int num_samples = 0;
	int n;
	int num_values  = static_cast<int>( n_values.size() );
	for ( int i=0; i<num_values; ++i ) {
		n = static_cast<int>( n_values[i] );
		if ( n >= 0 and n <= max_n ) {
			p_v_n[n] += 1.0;

			++num_samples;
		}
	}

	// Normalize bins (note that bin size = 1)
	double num_samples_d = static_cast<double>( num_samples );
	for ( int b=0; b<num_bins; ++b ) {
		p_v_n[b] /= num_samples_d;
	}

	// Statistics
	avg_n = Statistics::average( n_values );
	var_n = Statistics::variance( n_values );

	return;
}
