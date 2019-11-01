// OutputHandler

#ifndef OUTPUT_HANDLER_H
#define OUTPUT_HANDLER_H

// Standard headers
#include <array>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Project headers
#include "MpiCommunicator.h"
#include "SimulationState.h"
#include "Statistics.h"

class OutputHandler
{
 public:
	OutputHandler(
		const SimulationState& simulation_state,
		MpiCommunicator& mpi_communicator,
		const bool debug_mode = false
	);

	// Use a virtual destructor so that derived classes may implement destruction
	virtual ~OutputHandler() {}

	// Open/close ofstreams
	virtual void openOutputFiles() { return; };
	virtual void closeOutputFiles() { return; };

	virtual void writeOutputForFrame() { return; }

 protected:
	const SimulationState& simulation_state_;
	MpiCommunicator& mpi_communicator_;
	bool debug_mode_ = false;

	// Writes the given list of indices (and "num_atoms_per_index" indices after it)
	// to "ofs" as a single line of serial numbers
	// - For each index k, writes: indices[k]+1, indices[k]+2, ..., indices[k]+num_atoms_per_index
	// - indices (indexed from zero) are written as serials (indexed from 1)
	virtual void writeIndicesForFrame(
		std::ofstream& ofs,
		const std::vector<int>& indices,
		const int num_atoms_per_index = 1
	);

	// Compute P_v(N) from a series of numbers and bin it over the range [0,max_N],
	// where each bin "b" spans the range [b,b+1)
	void make_PvN(
		const std::vector<double>& n_values,
		const int max_N,
		// Output
		std::vector<double>& bins,
		std::vector<double>& pvn,
		double& avg_N,
		double& var_N
	) const;
};

#endif /* OUTPUT_HANDLER_H */
