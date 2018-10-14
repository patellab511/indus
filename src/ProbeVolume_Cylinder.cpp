/* ProbeVolume_Cylinder.cpp
 *
 * ABOUT: 
 * DEVELOPMENT:
 * - TODO Validate isInProbeVolume further
 */

#include "ProbeVolume_Cylinder.h"

// Register this probe volume
namespace ProbeVolumeRegistry {
static const RegisterInFactory<
		ProbeVolume, 
		ProbeVolume_Cylinder, 
		const std::string, 
		ProbeVolume::ProbeVolumeInputPack
	> register_ProbeVolume_Cylinder("cylinder");
}


ProbeVolume_Cylinder::ProbeVolume_Cylinder(
	ProbeVolume::ProbeVolumeInputPack& input_pack
) 
 : ProbeVolume(input_pack),
   axis_index_(Z_),
   axis_label_("z"),
   exclude_atoms_below_base_(false),
   switching_function_r_(1.0, sigma_, alpha_c_),
   switching_function_zeta_(0.0, 1.0, sigma_, alpha_c_)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Coarse-graining parameters
	input_parameter_pack.readNumber("sigma",   KeyType::Optional, sigma_);
	input_parameter_pack.readNumber("alpha_c", KeyType::Optional, alpha_c_);

	// Cylinder dimensions
	input_parameter_pack.readNumber("radius", KeyType::Required, r_);
	input_parameter_pack.readNumber("height", KeyType::Required, h_);

	// Cylinder location
	input_parameter_pack.readArray("base", KeyType::Required, x_base_);

	// Principle axis of the cylinder
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

	// Misc. options
	input_parameter_pack.readFlag("exclude_atoms_below_base", KeyType::Optional, 
	                              exclude_atoms_below_base_);

	// Update geometry based on what was parsed
	setCoarseGrainingParameters(sigma_, alpha_c_);
	setGeometry(x_base_, r_, h_);
}


void ProbeVolume_Cylinder::setCoarseGrainingParameters(const double sigma, const double alpha_c)
{
	sigma_   = sigma;
	alpha_c_ = alpha_c;

	// Update other member variables accordingly
	switching_function_r_.setCoarseGrainingParameters(sigma_, alpha_c_);
	switching_function_zeta_.setCoarseGrainingParameters(sigma_, alpha_c_);
	this->setGeometry(x_base_, r_, h_);
}


void ProbeVolume_Cylinder::setShellParameters(
		const double width_shell_1, const double width_shell_2, const double alpha_c_shells)
{
	width_shell_1_  = width_shell_1;
	width_shell_2_  = width_shell_2;
	alpha_c_shells_ = alpha_c_shells;

	// Update geometry for new shell widths
	this->setGeometry(x_base_, r_, h_);
}


// Defines all necessary variables for the implementation of a probe cylinder
// which can cross the periodic boundaries
void ProbeVolume_Cylinder::setGeometry(const Real3& x_base, const double radius, 
                                       const double height) 
{
	// Probe volume (nominal dimensions)
	x_base_ = x_base;
	r_ = radius;
	h_ = height; 

	// Radial switching function
	switching_function_r_.set_r_max(r_);

	// Effective probe volume incorporates smoothing over length scale alpha_c
	// on all sides
	r_eff_ = r_ + alpha_c_;
	h_eff_ = h_ + 2.0*alpha_c_;

	// Shell 1
	r_shell_1_ = r_eff_ + (width_shell_1_ + alpha_c_shells_);
	h_shell_1_ = h_eff_ + 2.0*(width_shell_1_ + alpha_c_shells_);

	// Shell 2
	r_shell_2_ = r_shell_1_ + (width_shell_2_ + alpha_c_shells_);
	h_shell_2_ = h_shell_1_ + 2.0*(width_shell_2_ + alpha_c_shells_);

	// Precompute squared radii and half-heights
	h_eff_h_     = 0.5*h_eff_;
	h_shell_1_h_ = 0.5*h_shell_1_;
	h_shell_2_h_ = 0.5*h_shell_2_;

	r_eff_sq_     = r_eff_*r_eff_;
	r_shell_1_sq_ = r_shell_1_*r_shell_1_;
	r_shell_2_sq_ = r_shell_2_*r_shell_2_;

	// Bounds along cylinder axis (before coarse-graining)
	zeta_min_ = x_base_[axis_index_];
	zeta_max_ = zeta_min_ + h_;
	switching_function_zeta_.setLimits(zeta_min_, zeta_max_);

	// Center of the probe cylinder (midpoint along the principal axis)
	center_ = x_base_;
	center_[axis_index_] += 0.5*h_;

	return;
}


void ProbeVolume_Cylinder::updateUsingSimulationState()
{
	// Nothing to do here
	return;
}


// Assign a position to one of the following regions: probe volume, shell 1, shell 2,
// or "other" (irrelevant to local neighbor list)
bool ProbeVolume_Cylinder::isInProbeVolume(
	const ProbeVolume::Real3& x,
	double& h_v, double& htilde_v, Real3& dhtilde_v_dx, 
	bool& is_in_shell_1, bool& is_in_shell_2
) const
{
	// Compute the minimum image vector from the cylinder's center to the particle, 
	ProbeVolume::Real3 dx_center_i;
	double dist_sq;
	simulation_box_.calculateDistance(center_, x, dx_center_i, dist_sq);

	// Minimum image distances along r and zeta, where the origin is at
	// the center of the cylinder
	double axial_distance  = std::abs(dx_center_i[axis_index_]);
	double radial_distance_sq = dist_sq - axial_distance*axial_distance;

	// Default assumption: outside region of interest
	bool is_in_probe_volume = false;
	bool is_in_vtilde = false;
	is_in_shell_1 = false;
	is_in_shell_2 = false;

	// Check if particle is in shells, smoothed probe volume, or otherwise
	if ( radial_distance_sq > r_shell_2_sq_ or axial_distance > h_shell_2_h_ ) {
		// Since most particles won't be in either the probe volume or the shells,
		// it's most efficient to check this case first
		is_in_vtilde = false;
	}
	else if ( radial_distance_sq <= r_eff_sq_ && axial_distance <= h_eff_h_ ) {
		is_in_vtilde = true;
	}
	else if ( radial_distance_sq <= r_shell_1_sq_ && axial_distance <= h_shell_1_h_) {
		is_in_shell_1 = true;
	}
	else {
		is_in_shell_2 = true;
	}

	if ( exclude_atoms_below_base_ ) {
		if ( x[axis_index_] < x_base_[axis_index_] ) {
			// Default assumption: outside region of interest
			is_in_vtilde = false;
			is_in_shell_1 = false;
			is_in_shell_2 = false;
		}
	}

	if ( is_in_vtilde ) {
		// In coarse-grained probe volume: need to evaluate indicator functions 

		// Radial component
		double h_radial, htilde_radial, deriv_htilde_radial;
		double radial_distance = sqrt(radial_distance_sq);
		switching_function_r_.calculate(radial_distance, need_derivatives_,
		                                h_radial, htilde_radial, deriv_htilde_radial);

		// Axial component
		Real3 x_in_box = x;
		simulation_box_.putInBox( x_in_box );  // sometimes, positions are slightly outside the box
		double h_zeta, htilde_zeta, deriv_htilde_zeta;
		switching_function_zeta_.calculate(x_in_box[axis_index_], need_derivatives_, 
		                                   h_zeta, htilde_zeta, deriv_htilde_zeta);

		htilde_v = htilde_radial*htilde_zeta;

		if ( h_radial == 1.0 and h_zeta == 1.0 ) {
			// Particle is in the nominal probe volume, and therefore part of the count for n_v
			h_v = 1.0;
			is_in_probe_volume = true;
		}
		else {
			h_v = 0.0;
			is_in_probe_volume = false;
		}

		if ( need_derivatives_ ) {
			// Apply the chain rule carefully, given that the origin for the r-direction
			// is the axis of the cylinder
			if ( radial_distance > 0.0 ) {
				for ( int d=0; d<DIM_; ++d ) {
					if ( d != axis_index_ ) {
						dhtilde_v_dx[d] = deriv_htilde_radial*(x_in_box[d] - x_base_[d])/radial_distance * htilde_zeta;
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
			dhtilde_v_dx[axis_index_] = htilde_radial*deriv_htilde_zeta;
		}
	}
	else {
		// Not in vtilde
		h_v      = 0.0;
		htilde_v = 0.0;
		if ( need_derivatives_ ) {
			for ( int d=0; d<DIM_; ++d ) { dhtilde_v_dx[d] = 0.0; }
		}
	}

	return is_in_probe_volume;
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume_Cylinder::getInputSummary() const
{
	std::stringstream ss;

	ss << "# Probe_volume = cylinder\n"
	   << "#   parallel_to_axis = " << axis_label_ << "\n"
	   << "#   center_of_base = {" 
	       << x_base_[0] << ", " << x_base_[1] << ", " << x_base_[2] << "} [nm]\n"
	   << "#   radius = " << r_ << " [nm]\n"
	   << "#   height = " << h_ << " [nm]\n"
	   << getSharedAttributesSummary();

	if ( exclude_atoms_below_base_ ) {
		ss << "#   WARNING: EXCLUDED ATOMS BELOW BASE.\n";
	}
	return ss.str();
}
