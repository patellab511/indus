// ProbeVolume_Cylinder.cpp

#include "ProbeVolume_Cylinder.h"

// Register this probe volume
namespace ProbeVolumeRegistry {
static const Register<ProbeVolume_Cylinder> register_ProbeVolume_Cylinder("cylinder");
}


ProbeVolume_Cylinder::ProbeVolume_Cylinder(
	ProbeVolume::ProbeVolumeInputPack& input_pack
): 
	ProbeVolume(input_pack),
	r_max_(1.0),
	zeta_range_{{0.0, 1.0}},
	axis_index_(Z_), 
	axis_label_("z"),
	switching_function_r_(r_max_, sigma_, alpha_c_, width_shell_1_, width_shell_2_, false),
	switching_function_zeta_(zeta_range_[0], zeta_range_[1], sigma_, alpha_c_, 
	                         width_shell_1_, width_shell_2_, false),
	// Variables for when this probe volume is a cylindrical shell with an inner axial range
	has_zeta_range_inner_(false), has_r_min_(false),
	r_min_(-1.0), 
	zeta_range_inner_{{0.2, 0.8}}, 
	switching_function_r_min_(r_min_, sigma_, alpha_c_, width_shell_1_, width_shell_2_, false),
	switching_function_zeta_inner_(zeta_range_inner_[0], zeta_range_inner_[1], sigma_, alpha_c_, 
	                               width_shell_1_, width_shell_2_, false),
	// Hacks to exclude atoms in particular regions
	exclude_atoms_below_base_(false)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Principal axis of the cylinder
	if ( input_parameter_pack.readString("axis", KeyType::Optional, axis_label_) ) {
		StringTools string_tools;
		axis_label_ = string_tools.toLowercase(axis_label_);
		if      ( axis_label_ == "x" ) { axis_index_ = X_; }
		else if ( axis_label_ == "y" ) { axis_index_ = Y_; }
		else if ( axis_label_ == "z" ) { axis_index_ = Z_; }
		else {
			throw std::runtime_error("invalid cylinder axis: " + axis_label_);
		}
	}

	// Cylinder location
	input_parameter_pack.readArray("base", KeyType::Required, x_base_);

	// Cylinder radius (outer)
	bool found_r_max = input_parameter_pack.readNumber("r_max", KeyType::Optional, r_max_);
	if ( not found_r_max ) {
		// Alias for r_max
		bool found_radius = input_parameter_pack.readNumber("radius", KeyType::Required, r_max_);
		if ( not found_radius ) {
			throw std::runtime_error("For probe cylinders, you must specify either 'r_max' or 'radius'.");
		}
	}

	// Limits along principle axis (outer)
	bool found_zeta_range = input_parameter_pack.readArray("axis_range", KeyType::Optional, zeta_range_);
	if ( not found_zeta_range ) {
		// Alternative key
		bool found_height = input_parameter_pack.readNumber("height", KeyType::Optional, h_);
		if ( not found_height ) {
			throw std::runtime_error("For probe cylinders, you must specify either 'height' or 'axis_range'.");
		}
		zeta_range_ = {{ x_base_[axis_index_],  x_base_[axis_index_] + h_ }};
	}
	if ( x_base_[axis_index_] != zeta_range_[0] ) {
		throw std::runtime_error("error: cylinder axis range and location of base don't agree.");
	}

	// Check for inner limits
	has_r_min_ = input_parameter_pack.readNumber("r_min", KeyType::Optional, r_min_);
	has_zeta_range_inner_ = input_parameter_pack.readArray("axis_range_inner", KeyType::Optional,
	                                                       zeta_range_inner_);

	// Misc. options
	input_parameter_pack.readFlag("exclude_atoms_below_base", KeyType::Optional, 
	                              exclude_atoms_below_base_);

	// Update geometry based on what was parsed
	setGeometry();
}


void ProbeVolume_Cylinder::setGeometry()
{
	// Radial switching function
	switching_function_r_.setBoundaries(r_max_, width_shell_1_, width_shell_2_);
	switching_function_r_.setCoarseGrainingParameters(sigma_, alpha_c_);

	// Bounds along cylinder axis
	switching_function_zeta_.setBoundaries(zeta_range_[0], zeta_range_[1], width_shell_1_, width_shell_2_);
	switching_function_zeta_.setCoarseGrainingParameters(sigma_, alpha_c_);

	// Center of the probe cylinder (midpoint between axial boundaries along the principal axis)
	center_ = x_base_;
	h_ = zeta_range_[1] - zeta_range_[0];
	center_[axis_index_] += 0.5*h_;

	// TODO Shells
	// - Currently not supported when there is an inner cylinder
	// - For now, throw an exception if the combination is used
	if ( (width_shell_1_ != 0.0 or width_shell_1_ != 0.0) and (has_r_min_ or has_zeta_range_inner_) ){
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
					 << "  Unsupported input: Probe volume \"shells\" have not been implemented "
		         << "for cylindrical shell probe volumes.\n"
		       << "  This means that this geometry cannot currently be used for Steinhardt parameters.\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Inner radius
	if ( has_r_min_ ) {
		if ( r_min_ > r_max_ ) {
			std::stringstream err_ss;
			err_ss << "Error in " << FANCY_FUNCTION << "\n"
						 << "  Inner radius of cylinder must be smaller than outer radius.\n"
						 << "  Input: r_min = " << r_min_ << ", r_max = " << r_max_ << "\n";
			throw std::runtime_error( err_ss.str() );
		}
		switching_function_r_min_.setBoundaries(r_min_, width_shell_1_, width_shell_2_);
		switching_function_r_min_.setCoarseGrainingParameters(sigma_, alpha_c_);
	}

	// Inner zeta-range
	if ( has_zeta_range_inner_ ) {
		if ( zeta_range_inner_[0] < zeta_range_[0] or zeta_range_inner_[1] > zeta_range_[1] ) {
			std::stringstream err_ss;
			std::string zeta_min = axis_label_ + "min";
			std::string zeta_max = axis_label_ + "max";
			err_ss << "Error in " << FANCY_FUNCTION << "\n"
						 << "  Inner cylinder range must be within outer cylinder range.\n"
						 << "  Input: \n"
			       << "    Inner: " << zeta_min << " = " << zeta_range_inner_[0] << ", " 
			                        << zeta_max << " = " << zeta_range_inner_[1] << "\n"
			       << "    Outer: " << zeta_min << " = " << zeta_range_[0] << ", " 
			                        << zeta_max << " = " << zeta_range_[1] << "\n";
			throw std::runtime_error( err_ss.str() );
		}
		switching_function_zeta_inner_.setBoundaries(
			zeta_range_inner_[0], zeta_range_inner_[1], 
			width_shell_1_, width_shell_2_
		);
		switching_function_zeta_inner_.setCoarseGrainingParameters(sigma_, alpha_c_);
	}

	return;
}


BoundingBox ProbeVolume_Cylinder::constructBoundingBox() const
{
	// "Extra" distance beyond nominal cutoffs
	double dr_buffer = alpha_c_ + bounding_box_tol_;

	Real3 x_lower, x_upper;
	for ( int d=0; d<DIM_; ++d ) {
		if ( d == axis_index_ ) {
			// Axial coordinate
			x_lower[d] = zeta_range_[0] - dr_buffer;
			x_upper[d] = zeta_range_[1] + dr_buffer;
		}
		else {
			// Part of radial coordinate
			x_lower[d] = x_base_[d] - (r_max_ + dr_buffer);
			x_upper[d] = x_base_[d] + (r_max_ + dr_buffer);
		}
	}

	return BoundingBox(x_lower, x_upper, simulation_box_);
}


// Assign a position to one of the following regions: probe volume, shell 1, shell 2,
// or "other" (irrelevant to local neighbor list)
void ProbeVolume_Cylinder::calculateIndicator(
	const ProbeVolume::Real3& x,
	double& h_v, double& htilde_v, Real3& dhtilde_v_dx, RegionEnum& region
) const
{
	// Find the image closest to the center of the (outer) cylinder
	// (TODO: if ghost atoms are implemented, this will not be necessary here: can
	// just compute the difference vector and distance directly, since PBCs are handled elsewhere)
	Real3 x_shifted, shift;
	simulation_box_.calculateShift(x, center_, x_shifted, shift);

	// Minimum image distance along r, where the origin is at
	// the center of the cylinder
	double r_dist = 0.0;
	double dx;
	for ( int d=0; d<DIM_; ++d ) {
		if ( d !=	axis_index_ ) {
			dx = x_shifted[d] - center_[d];
			r_dist += dx*dx;
		}
	}
	r_dist = sqrt(r_dist);

	// Radial component
	double h_r, htilde_r, deriv_htilde_r;
	RegionEnum region_r = RegionEnum::Unimportant;
	switching_function_r_.calculate(r_dist, need_derivatives_,
																	h_r, htilde_r, deriv_htilde_r, region_r);

	// Axial component
	double h_zeta = 0.0, htilde_zeta = 0.0, deriv_htilde_zeta = 0.0;
	RegionEnum region_zeta = RegionEnum::Unimportant;
	switching_function_zeta_.calculate(x_shifted[axis_index_], need_derivatives_, 
	                                   h_zeta, htilde_zeta, deriv_htilde_zeta, region_zeta);

	// Combine components
	htilde_v = htilde_r*htilde_zeta;
	if ( htilde_v > 1.0 ) { htilde_v = 1.0; }

	// Handle regions
	// - Default: outside probe volume
	h_v = 0.0;
	bool is_in_probe_volume = ( h_r == 1.0 and h_zeta == 1.0 );
	region = RegionEnum::Unimportant;
	if ( htilde_v > 0.0 ) {
		region = RegionEnum::Vtilde;
		if ( is_in_probe_volume ) {
			// Particle is in the nominal probe volume, and therefore part of the count for n_v
			h_v = 1.0;
		}

		if ( need_derivatives_ ) {
			// Apply the chain rule carefully, given that the origin for the r-direction
			// is the axis of the cylinder
			if ( r_dist > 0.0 ) {
				for ( int d=0; d<DIM_; ++d ) {
					if ( d != axis_index_ ) {
						dhtilde_v_dx[d] = deriv_htilde_r*(x_shifted[d] - x_base_[d])/r_dist * htilde_zeta;
					}
				}
			}
			else {
				for ( int d=0; d<DIM_; ++d ) {
					if ( d != axis_index_ ) {
						dhtilde_v_dx[d] = 0.0;
					}
				}
			}
			// zeta-direction
			dhtilde_v_dx[axis_index_] = htilde_r*deriv_htilde_zeta;
		}
	}
	else {
		// Not in vtilde
		if ( not ( region_r == RegionEnum::Unimportant or region_zeta == RegionEnum::Unimportant ) ) {
			// Somewhere in the shells
			if ( region_r == RegionEnum::Shell_2 or region_zeta == RegionEnum::Shell_2 ) {
				region = RegionEnum::Shell_2;
			}
			else if ( region_r == RegionEnum::Shell_1 or region_zeta == RegionEnum::Shell_1 ) {
				region = RegionEnum::Shell_1;
			}
		}

		if ( need_derivatives_ ) {
			dhtilde_v_dx.fill(0.0);
		}
	}

	// Inner cylinder
	if ( has_r_min_ or has_zeta_range_inner_ ) {
		// Note: It's important to set the defaults for each inner dimension to 
		// interpret the particle as "in vtilde", in case that dimension is not checked.

		// r_min
		double h_r_min = 1.0, htilde_r_min = 1.0, deriv_htilde_r_min = 0.0;
		RegionEnum region_r_min = RegionEnum::Vtilde;
		if ( has_r_min_ ) {
			switching_function_r_min_.calculate(
				r_dist, need_derivatives_, 
				h_r_min, htilde_r_min, deriv_htilde_r_min, region_r_min
			);
		}

		// zeta_range_inner
		double h_zeta_inner = 1.0, htilde_zeta_inner = 1.0, deriv_htilde_zeta_inner = 0.0;
		RegionEnum region_zeta_inner = RegionEnum::Vtilde;
		if ( has_zeta_range_inner_ ) {
			switching_function_zeta_inner_.calculate(
				x_shifted[axis_index_], need_derivatives_, 
				h_zeta_inner, htilde_zeta_inner, deriv_htilde_zeta_inner, region_zeta_inner
			);
		}

		// Combine r_min and zeta_range_inner to get v_inner
		double htilde_v_inner   = htilde_r_min*htilde_zeta_inner;
		/*
		RegionEnum region_inner = RegionEnum::Unimportant;  // TODO Shells
		if ( htilde_v_inner > 0.0 ) {
			region_inner = RegionEnum::Vtilde;
		}
		*/
		bool is_in_v_inner = ( h_r_min == 1.0 and h_zeta_inner == 1.0 );
		Real3 dhtilde_v_inner_dx;
		if ( need_derivatives_ ) {
			if ( r_dist > 0.0 ) {
				for ( int d=0; d<DIM_; ++d ) {
					if ( d != axis_index_ ) {
						dhtilde_v_inner_dx[d] = deriv_htilde_r_min*(x_shifted[d] - x_base_[d])/r_dist 
						                         *htilde_zeta_inner;
					}
				}
			}
			else {
				for ( int d=0; d<DIM_; ++d ) {
					if ( d != axis_index_ ) {
						dhtilde_v_inner_dx[d] = 0.0;
					}
				}
			}
			// zeta-direction
			dhtilde_v_inner_dx[axis_index_] = htilde_r_min*deriv_htilde_zeta_inner;
		}


		//----- Combine Cylinders -----// 

		// Boolean algebra: v = v_outer AND (NOT(v_inner))

		// Evaluate htilde_v
		double htilde_v_outer = htilde_v;  // save value for use below
		htilde_v = htilde_v_outer*(1.0 - htilde_v_inner);
		if      ( htilde_v > 1.0 ) { htilde_v = 1.0; }  // beware of roundoff error
		else if ( htilde_v < 0.0 ) { htilde_v = 0.0; }

		// Update region (TODO: shells)
		region = RegionEnum::Unimportant;
		if ( htilde_v > 0.0 ) {
			region = RegionEnum::Vtilde;
		}

		// Evaluate h_v
		is_in_probe_volume = is_in_probe_volume and (not is_in_v_inner);
		if ( is_in_probe_volume ) { h_v = 1.0; }
		else                      { h_v = 0.0; }

		if ( need_derivatives_ ) {
			if ( region == RegionEnum::Vtilde ) {
				// Note: at this point, dhtilde_v_dx = dhtilde_v_outer_dx
				for ( int d=0; d<DIM_; ++d ) {
					dhtilde_v_dx[d] = dhtilde_v_dx[d]*(1.0 - htilde_v_inner) 
					                    - htilde_v_outer*dhtilde_v_inner_dx[d];
				}
			}
			else {
				for ( int d=0; d<DIM_; ++d ) {
					dhtilde_v_dx[d] = 0.0;
				}
			}
		}
	} // end section for inner cylinder

	// Hack for excluding atoms below the base of the probe volume
	if ( exclude_atoms_below_base_ and (x[axis_index_] < x_base_[axis_index_]) ) {
		is_in_probe_volume = false;
		h_v = htilde_v = 0.0;
		region = RegionEnum::Unimportant;

		if ( need_derivatives_ ) {
			dhtilde_v_dx.fill(0.0);
		}
	}

	return;
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume_Cylinder::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;

	ss << prepend_string << "Probe_volume = cylinder\n"
	   << prepend_string << "  parallel_to_axis = " << axis_label_ << "\n"
	   << prepend_string << "  center_of_base = {" 
	       << x_base_[0] << ", " << x_base_[1] << ", " << x_base_[2] << "} [nm]\n";
	// Radii
	if ( has_r_min_ ) {
		ss << prepend_string << "  r_min = " << r_min_ << " [nm]\n";
	}
	ss << prepend_string << "  r_max = " << r_max_ << " [nm]\n";
	// Axis ranges
	ss << prepend_string << "  axis_range";
	if ( has_zeta_range_inner_ ) { ss << "(outer)"; }
	ss << " = [ " << zeta_range_[0] << "  " << zeta_range_[1] << "] [nm]\n";
	if ( has_zeta_range_inner_ ) {
		ss << prepend_string << "  axis_range(inner) = [ " << zeta_range_inner_[0] << "  " 
		                                   << zeta_range_inner_[1] << "] [nm]\n";
	}
	ss << getSharedAttributesSummary(prepend_string + "  ");

	if ( exclude_atoms_below_base_ ) {
		ss << prepend_string << "  WARNING: EXCLUDED ATOMS BELOW BASE.\n";
	}
	return ss.str();
}
