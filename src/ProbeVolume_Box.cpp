#include "ProbeVolume_Box.h"


// Register this probe volume
namespace ProbeVolumeRegistry {
static const RegisterInFactory<
		ProbeVolume, 
		ProbeVolume_Box, 
		const std::string, 
		ProbeVolume::ProbeVolumeInputPack
	> register_ProbeVolume_Box("box");
}


ProbeVolume_Box::ProbeVolume_Box(
	ProbeVolume::ProbeVolumeInputPack& input_pack
) 
 : ProbeVolume(input_pack)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// x, y, and z bounds (assumes orthorhombic box)
	std::array<double,2> x_range, y_range, z_range;
	input_parameter_pack.readArray("x_range", KeyType::Required, x_range);
	input_parameter_pack.readArray("y_range", KeyType::Required, y_range);
	input_parameter_pack.readArray("z_range", KeyType::Required, z_range);
	box_offset_[X_DIM] = x_range[0];
	box_offset_[Y_DIM] = y_range[0];
	box_offset_[Z_DIM] = z_range[0];
	box_matrix_[X_DIM][X_DIM] = x_range[1] - x_range[0];
	box_matrix_[Y_DIM][Y_DIM] = y_range[1] - y_range[0];
	box_matrix_[Z_DIM][Z_DIM] = z_range[1] - z_range[0];
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( a != b ) {
				box_matrix_[a][b] = 0.0;
			}
		}
	}

	// Coarse-graining parameters
	input_parameter_pack.readNumber("sigma",   KeyType::Optional, sigma_);
	input_parameter_pack.readNumber("alpha_c", KeyType::Optional, alpha_c_);

	// Update state based on what was parsed
	setCoarseGrainingParameters(sigma_, alpha_c_);
	setGeometry(box_offset_, box_matrix_);
}


ProbeVolume_Box::ProbeVolume_Box(
		ProbeVolume::ProbeVolumeInputPack& input_pack,
		const Real3& box_offset, const Box& box_matrix,
		const double sigma, const double alpha_c
) 
 : ProbeVolume(input_pack),
   box_offset_(box_offset),
   box_matrix_(box_matrix)
{
	// Ensure switching function has correct parameters
  sigma_   = sigma;
  alpha_c_ = alpha_c;
	this->setCoarseGrainingParameters(sigma_, alpha_c_);

	this->setGeometry(box_offset_, box_matrix_);
}


void ProbeVolume_Box::setCoarseGrainingParameters(const double sigma, const double alpha_c)
{
	sigma_   = sigma;
	alpha_c_ = alpha_c; 

	// Update other member variables accordingly
	for ( int d=0; d<DIM_; ++d ) {
		switching_functions_shifted_box_[d].setCoarseGrainingParameters(sigma_, alpha_c_);
	}
	this->setGeometry(box_offset_, box_matrix_);
}


void ProbeVolume_Box::setShellParameters(
		const double width_shell_1, const double width_shell_2, const double alpha_c_shells) 
{
	width_shell_1_  = width_shell_1;
	width_shell_2_  = width_shell_2;
	alpha_c_shells_ = alpha_c_shells;

	// Update related member variables
	this->setGeometry(box_offset_, box_matrix_);
}


void ProbeVolume_Box::setGeometry(const Real3& box_offset, const Box& box_matrix)
{
	// Nominal box geometry
	box_offset_ = box_offset;
	box_matrix_ = box_matrix;

	// Assumes orthorhombic box where the nominal box does not cross PBCs,
	// but the shells might
	for ( int d=0; d<DIM_; ++d ) {
		center_[d] = box_offset_[d] + 0.5*box_matrix_[d][d];

		box_half_lengths_[d] = 0.5*box_matrix_[d][d];
		box_half_lengths_eff_[d] = box_half_lengths_[d] + alpha_c_;
		box_shell_1_half_lengths_eff_[d] = box_half_lengths_eff_[d] + 
		                                  (width_shell_1_ + alpha_c_shells_);
		box_shell_2_half_lengths_eff_[d] = box_shell_1_half_lengths_eff_[d] + 
		                                  (width_shell_2_ + alpha_c_shells_);

		// Switching function bounds
		switching_functions_shifted_box_[d].setLimits( -box_half_lengths_[d], box_half_lengths_[d] );
	}

	return;
}


void ProbeVolume_Box::updateUsingSimulationState() 
{
	// Nothing to do for this geometry
	return;
}


// Assign a position to one of the following regions: probe volume, shell 1, shell 2,
// or "other" (irrelevant to local neighbor list)
bool ProbeVolume_Box::isInProbeVolume(const ProbeVolume::Real3& x,
                                      double& h_v, double& htilde_v, Real3& dhtilde_v_dx,
                                      bool& is_in_shell_1, bool& is_in_shell_2) const 
{
	// Compute the minimum image vector from the box's center to the particle, "i"
	ProbeVolume::Real3 dx_center_i;
	double dist_sq;
	simulation_box_.calculateDistance(center_, x, dx_center_i, dist_sq);

	// Default assumption: outside region of interest
	bool is_in_vtilde = false;
	is_in_shell_1 = false;
	is_in_shell_2 = false;

	ProbeVolume::Real3 abs_dx_center_i;
	for ( int d=0; d<DIM_; ++d ) {
		abs_dx_center_i[d] = std::abs( dx_center_i[d] );
	}

	// Assumes orthorhombic box
	if ( abs_dx_center_i[X_DIM] <= box_half_lengths_eff_[X_DIM] and
	     abs_dx_center_i[Y_DIM] <= box_half_lengths_eff_[Y_DIM] and
	     abs_dx_center_i[Z_DIM] <= box_half_lengths_eff_[Z_DIM]
	) {
		// In coarse-grained probe volume
		is_in_vtilde = true;
	}
	else if ( abs_dx_center_i[X_DIM] > box_shell_2_half_lengths_eff_[X_DIM] or
	          abs_dx_center_i[Y_DIM] > box_shell_2_half_lengths_eff_[Y_DIM] or
	          abs_dx_center_i[Z_DIM] > box_shell_2_half_lengths_eff_[Z_DIM]
	) {
		// Completely outside region of interest
	}
	else if ( abs_dx_center_i[X_DIM] <= box_shell_1_half_lengths_eff_[X_DIM] and
	          abs_dx_center_i[Y_DIM] <= box_shell_1_half_lengths_eff_[Y_DIM] and
	          abs_dx_center_i[Z_DIM] <= box_shell_1_half_lengths_eff_[Z_DIM]
	) {
		// Shell 1
		is_in_shell_1 = true;
	}
	else {
		// Must be in shell 2
		is_in_shell_2 = true;
	}

	// Edge cases at the boundary between vtilde and shell 1 
	// - The following assumes an orthorhombic box
	if ( is_in_vtilde or is_in_shell_1 ) {
		// Ensure the atom is fully within the box (XTC positions are sometimes
		// slightly negative)
		Real3 x_in_box = x; 
		simulation_box_.putInBox(x_in_box);

		if ( is_in_shell_1 ) {
			// Check if the atom is really in vtilde, and was placed in shell 1 due to roundoff error
			bool is_at_edge_of_vtilde = true;
			for ( int d=0; d<DIM_; ++d ) {
				if ( not ( x_in_box[d] >= box_offset_[d] and 
				           x_in_box[d] <= box_offset_[d] + box_matrix_[d][d] )
				) {
					is_at_edge_of_vtilde = false;
					break;
				}
			}

			if ( is_at_edge_of_vtilde ) {
				is_in_vtilde  = true;
				is_in_shell_1 = false;
			}
		}
	}

	bool is_in_probe_volume = false;
	if ( is_in_vtilde ) {
		// In coarse-grained probe volume: need to evaluate indicator functions
		// - Assumes orthorhombic box
		htilde_v = 1.0;
		Real3 h_components, htilde_components, deriv_htilde_components;
		for ( int d=0; d<DIM_; ++d ) {
			switching_functions_shifted_box_[d].calculate( 
					dx_center_i[d], need_derivatives_,
			    h_components[d], htilde_components[d], deriv_htilde_components[d] );

			htilde_v *= htilde_components[d];
		}

		// Compute h_v
		if ( h_components[X_DIM] == 1.0 and 
		     h_components[Y_DIM] == 1.0 and 
		     h_components[Z_DIM] == 1.0 
		) {
			h_v = 1.0;
			is_in_probe_volume = true;
		}

		if ( need_derivatives_ ) {
			for ( int d=0; d<DIM_; ++d ) {
				dhtilde_v_dx[d] = deriv_htilde_components[d];

				for ( int f=0; f<DIM_; ++f ) {
					if ( d != f ) {
						dhtilde_v_dx[d] *= htilde_components[f];
					}
				}
			}
		}
	}
	else {
		// Not in either v or vtilde
		h_v      = 0.0;
		htilde_v = 0.0;
		if ( need_derivatives_ ) {
			for ( int d=0; d<DIM_; ++d ) {
				dhtilde_v_dx[d] = 0.0;
			}
		}
	}

	return is_in_probe_volume;
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume_Box::getInputSummary() const
{
	std::stringstream ss;

	ss << "# Probe_volume = box\n";

	ss << "#   x_offset(lower-left corner) = {";
	for ( int d=0; d<DIM_; ++d ) {
		ss << " " << box_offset_[d];
	}
	ss << " }\n";

	ss << "#   x_center = {";
	for ( int d=0; d<DIM_; ++d ) {
		ss << " " << center_[d];
	}
	ss << " }\n";

	ss << "#   box_matrix = {\n";
	for ( int a=0; a<DIM_; ++a ) {
		ss << "#    ";
		for ( int b=0; b<DIM_; ++b ) {
			ss << " " << box_matrix_[a][b];
		}
		ss << "\n";
	}
	ss << "#   }\n";

	ss << getSharedAttributesSummary();

	/*
	ss << " Coarse-grained half-lengths\n";
	ss << "   box: ";
	for ( int d=0; d<DIM_; ++d ) {
		ss << " " << box_half_lengths_eff_[d];
	}
	ss << "\n";
	ss << "   shell_1: ";
	for ( int d=0; d<DIM_; ++d ) {
		ss << " " << box_shell_1_half_lengths_eff_[d];
	}
	ss << "\n";
	ss << "\n";
	ss << "   shell_2: ";
	for ( int d=0; d<DIM_; ++d ) {
		ss << " " << box_shell_2_half_lengths_eff_[d];
	}
	ss << "\n";
	*/

	return ss.str();
}
