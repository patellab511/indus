/* ProbeVolume.cpp
 */

#include "ProbeVolume.h"

ProbeVolume::ProbeVolume(ProbeVolumeInputPack& input_pack):
	// Core variables
	simulation_state_( input_pack.simulation_state ),
	simulation_box_( simulation_state_.get_simulation_box() ),
	mpi_communicator_( input_pack.mpi_communicator ),
	need_derivatives_( input_pack.need_derivatives ),
	// Bounding box (for domain decomposition)
	bounding_box_(simulation_box_),
	use_default_bounding_box_(false),
	// Coarse-graining
	sigma_(0.01),
	alpha_c_(0.02),
	// Shells
	width_shell_1_(0.0),
	width_shell_2_(0.0),
	// Flags
	is_dynamic_(false)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Coarse-graining parameters
	input_parameter_pack.readNumber("sigma",   KeyType::Optional, sigma_);
	input_parameter_pack.readNumber("alpha_c", KeyType::Optional, alpha_c_);

	// Flags
	input_parameter_pack.readFlag("is_dynamic", KeyType::Optional, is_dynamic_);

	// TODO Probe volume atoms

	// Bounding box (i.e. region to partition for domain decomposition)
	input_parameter_pack.readFlag("use_default_bounding_box", KeyType::Optional, use_default_bounding_box_);
	if ( use_default_bounding_box_ ) {
		bounding_box_fxn_ = [this]()->BoundingBox { return ProbeVolume::constructBoundingBox(); };
	}
	else {
		bounding_box_fxn_ = [this]()->BoundingBox { return this->constructBoundingBox(); };
	}
}


BoundingBox ProbeVolume::constructBoundingBox() const
{
	// Default implementation: entire simulation box
	return BoundingBox(simulation_box_);
}


void ProbeVolume::setBoundingBox()
{
	bounding_box_ = bounding_box_fxn_();
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume::getInputSummary(const std::string& prepend_string) const
{
	return (prepend_string + "Probe_volume output_summary_not_yet_implemented\n")
	       + this->getSharedAttributesSummary(prepend_string + "  ");
}


// Get summary of attributes shared by all probe volumes
// ( e.g. shell widths, coarse-graining parameters, etc.)
std::string ProbeVolume::getSharedAttributesSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Coarse-graining\n"
		 << prepend_string << "  sigma            = " << sigma_   << " [nm]\n"
	   << prepend_string << "  alpha_c          = " << alpha_c_ << " [nm]\n"
	   << prepend_string << "Shells\n"
	   << prepend_string << "  width_shell_1  = " << width_shell_1_   << " [nm]\n"
	   << prepend_string << "  width_shell_2  = " << width_shell_2_   << " [nm]\n";

	return ss.str();
}
