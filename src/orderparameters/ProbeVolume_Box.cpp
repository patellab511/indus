#include "ProbeVolume_Box.h"


// Register this probe volume
namespace ProbeVolumeRegistry {
static const Register<ProbeVolume_Box> register_ProbeVolume_Box("box");
}


ProbeVolume_Box::ProbeVolume_Box(
	ProbeVolume::ProbeVolumeInputPack& input_pack
) 
 : ProbeVolume(input_pack)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// x, y, and z bounds (assumes orthorhombic box)
	input_parameter_pack.readArray("x_range", KeyType::Required, axis_ranges_[X_DIM]);
	input_parameter_pack.readArray("y_range", KeyType::Required, axis_ranges_[Y_DIM]);
	input_parameter_pack.readArray("z_range", KeyType::Required, axis_ranges_[Z_DIM]);

	for ( int a=0; a<DIM_; ++a ) {
		box_offset_[a] = axis_ranges_[a][0];
		for ( int b=0; b<DIM_; ++b ) {
			if ( a == b ) {
				box_matrix_[a][a] = axis_ranges_[a][1] - axis_ranges_[a][0];
			}
			else {
				box_matrix_[a][b] = 0.0;
			}
		}
	}

	// Update state based on what was parsed
	setGeometry();
}


ProbeVolume_Box::ProbeVolume_Box(
	ProbeVolume::ProbeVolumeInputPack& input_pack,
	const Real3& box_offset, const Box& box_matrix, const double sigma, const double alpha_c
):
	ProbeVolume(input_pack),
	box_offset_(box_offset),
	box_matrix_(box_matrix)
{
	// Ensure switching function has correct parameters
  sigma_   = sigma;
  alpha_c_ = alpha_c;
	
	setGeometry();
}


void ProbeVolume_Box::setGeometry()
{
	// Assumes orthorhombic box where the nominal box does not cross PBCs,
	// but the shells might
	if ( not is_orthorhombic() ) {
		throw std::runtime_error("Non-orthorhombic probe boxes are not currently supported.");
	}

	// Axis range form
	for ( int d=0; d<DIM_; ++d ) {
		inner_axis_ranges_[d][0] = axis_ranges_[d][0] + alpha_c_;
		inner_axis_ranges_[d][1] = axis_ranges_[d][1] - alpha_c_;

		outer_axis_ranges_[d][0] = axis_ranges_[d][0] - alpha_c_;
		outer_axis_ranges_[d][1] = axis_ranges_[d][1] + alpha_c_;

		center_[d]           = 0.5*(axis_ranges_[d][1] + axis_ranges_[d][0]);
		box_half_lengths_[d] = 0.5*(axis_ranges_[d][1] - axis_ranges_[d][0]);

		// Switching function bounds and coarse-graining
		switching_functions_shifted_box_[d].setBoundaries( 
				-box_half_lengths_[d], box_half_lengths_[d], width_shell_1_, width_shell_2_ );
		switching_functions_shifted_box_[d].setCoarseGrainingParameters(sigma_, alpha_c_);
	}

	// Box offset/matrix form
	for ( int a=0; a<DIM_; ++a ) {
		box_offset_[a] = axis_ranges_[a][0];
		for ( int b=0; b<DIM_; ++b ) {
			if ( a == b ) {
				box_matrix_[a][a] = axis_ranges_[a][1] - axis_ranges_[a][0];
			}
			else {
				box_matrix_[a][b] = 0.0;
			}
		}
	}

	return;
}


void ProbeVolume_Box::setGeometry(const Real3& x_lower, const Real3& x_upper)
{
	for ( int d=0; d<DIM_; ++d ) {
		axis_ranges_[d][0] = x_lower[d];
		axis_ranges_[d][1] = x_upper[d];
	}

	setGeometry();
}


void ProbeVolume_Box::setGeometry(const Real3& box_offset, const Box& box_matrix)
{
	// Nominal box geometry
	box_offset_ = box_offset;
	box_matrix_ = box_matrix;

	// Axis ranges (assuming an orthorhombic box)
	for ( int d=0; d<DIM_; ++d ) {
		axis_ranges_[d][0] = box_offset_[d];
		axis_ranges_[d][1] = box_offset_[d] + box_matrix[d][d];
	}

	setGeometry();
}


BoundingBox ProbeVolume_Box::constructBoundingBox() const
{
	// Assumes orthorhombic box
	if ( not simulation_box_.is_orthorhombic() ) {
		std::stringstream err_ss; 
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
					 << "  Non-orthorhombic simulation boxes are not supported.\n";
		throw std::runtime_error( err_ss.str() );
	}

	// "Extra" distance beyond nominal cutoffs
	double dx_buffer = alpha_c_ + bounding_box_tol_;

	Real3 x_lower, x_upper;
	for ( int d=0; d<DIM_; ++d ) {
		x_lower[d] = axis_ranges_[d][0] - dx_buffer;
		x_upper[d] = axis_ranges_[d][1] + dx_buffer;
	}

	return BoundingBox(x_lower, x_upper, simulation_box_);
}


void ProbeVolume_Box::calculateIndicator(
	const ProbeVolume::Real3& x,
	double& h_v, double& htilde_v, Real3& dhtilde_v_dx, RegionEnum& region
) const
{
	// Compute the minimum image vector from the box's center to the particle's position
	ProbeVolume::Real3 dx_center_i;
	double dist_sq;
	simulation_box_.calculateDistance(center_, x, dx_center_i, dist_sq);

	// Evaluate components of indicator function
	Real3 h_components, htilde_components, deriv_htilde_components;
	std::array<RegionEnum, DIM_> regions;
	for ( int d=0; d<DIM_; ++d ) {
		switching_functions_shifted_box_[d].calculate( 
			dx_center_i[d], need_derivatives_,
			h_components[d], htilde_components[d], deriv_htilde_components[d], regions[d] );
	}

	// Assumes orthorhombic box
	if ( regions[X_DIM] == RegionEnum::Unimportant or
	     regions[Y_DIM] == RegionEnum::Unimportant or
	     regions[Z_DIM] == RegionEnum::Unimportant 
	) {
		// Completely outside region of interest
		region = RegionEnum::Unimportant;
	}
	else if ( regions[X_DIM] == RegionEnum::Vtilde and
	          regions[Y_DIM] == RegionEnum::Vtilde and
	          regions[Z_DIM] == RegionEnum::Vtilde 
	) {
		// In coarse-grained probe volume
		region = RegionEnum::Vtilde;
	}
	else if ( (not (regions[X_DIM] == RegionEnum::Shell_2)) and
	          (not (regions[Y_DIM] == RegionEnum::Shell_2)) and
	          (not (regions[Z_DIM] == RegionEnum::Shell_2))
	) {
		// Shell 1
		region = RegionEnum::Shell_1;
	}
	else {
		// Must be in shell 2
		region = RegionEnum::Shell_2;
	}

	if ( region == RegionEnum::Vtilde ) {
		// Compute htilde_v
		htilde_v = 1.0;
		for ( int d=0; d<DIM_; ++d ) {
			htilde_v *= htilde_components[d];
		}

		// Compute h_v
		if ( h_components[X_DIM] == 1.0 and h_components[Y_DIM] == 1.0 and 
		     h_components[Z_DIM] == 1.0 
		) {
			h_v = 1.0;
		}
		else {
			h_v = 0.0;
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
			dhtilde_v_dx.fill(0.0);
		}
	}

	return;
}


std::string ProbeVolume_Box::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;

	ss << prepend_string << "Probe_volume = box\n";

	if ( is_orthorhombic() ) {
		for ( int d=0; d<DIM_; ++d ) {
			ss << prepend_string << "  " << axis_index_to_label_map_.at(d) << "_range = [ "
			   << axis_ranges_[d][0] << "  " << axis_ranges_[d][1] << " ]\n";
		}
	}
	else {
		ss << prepend_string << "  x_offset(lower-left corner) = {";
		for ( int d=0; d<DIM_; ++d ) {
			ss << " " << box_offset_[d];
		}
		ss << " }\n";

		ss << prepend_string << "  x_center = {";
		for ( int d=0; d<DIM_; ++d ) {
			ss << " " << center_[d];
		}
		ss << " }\n";

		ss << prepend_string << "  box_matrix = {\n";
		for ( int a=0; a<DIM_; ++a ) {
			ss << prepend_string << "   ";
			for ( int b=0; b<DIM_; ++b ) {
				ss << " " << box_matrix_[a][b];
			}
			ss << "\n";
		}
		ss << prepend_string << "  }\n";
	}
	ss << getSharedAttributesSummary(prepend_string + "  ");

	return ss.str();
}
