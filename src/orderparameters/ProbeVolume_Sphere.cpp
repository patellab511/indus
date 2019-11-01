/* ProbeVolume_Sphere.cpp
 *
 * ABOUT: 
 *
 */

#include "ProbeVolume_Sphere.h"

// Register this probe volume
namespace ProbeVolumeRegistry {
static const Register<ProbeVolume_Sphere> register_ProbeVolume_Sphere("sphere");
}


ProbeVolume_Sphere::ProbeVolume_Sphere(ProbeVolume::ProbeVolumeInputPack& input_pack):
	ProbeVolume(input_pack),
	r_min_(-1.0),
	r_max_(1.0),
	switching_function_r_shell_(r_min_, r_max_, sigma_, alpha_c_)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Sphere geometry
	input_parameter_pack.readNumber("r_min", KeyType::Optional, r_min_);
	input_parameter_pack.readNumber("r_max", KeyType::Required, r_max_);
	input_parameter_pack.readArray("center", KeyType::Required, center_);
	setGeometry();
}


ProbeVolume_Sphere::ProbeVolume_Sphere(
	const Real3& center, const double r_min, const double r_max,
	ProbeVolume::ProbeVolumeInputPack& input_pack
):
	ProbeVolume(input_pack),
	center_(center), r_min_(r_min), r_max_(r_max)
{
	setGeometry();
}


void ProbeVolume_Sphere::setGeometry()
{
	// Update switching function
	switching_function_r_shell_.setBoundaries(r_min_, r_max_, width_shell_1_, width_shell_2_);
	switching_function_r_shell_.setCoarseGrainingParameters(sigma_, alpha_c_);
}


void ProbeVolume_Sphere::setGeometry(const Real3& center, const double r_min, const double r_max)
{
	center_ = center;
	r_min_  = r_min;
	r_max_  = r_max;

	setGeometry();
}


void ProbeVolume_Sphere::set_center(const Real3& center) 
{
	center_ = center;
	setGeometry();
}

BoundingBox ProbeVolume_Sphere::constructBoundingBox() const
{
	// *** Assumes orthorhombic box *** //

	double dr_buffer = this->get_rtilde_max() + bounding_box_tol_;
	Real3 x_lower, x_upper;
	for ( int d=0; d<DIM_; ++d ) {
		x_lower[d] = center_[d] - dr_buffer;
		x_upper[d] = center_[d] + dr_buffer;
	}

	return BoundingBox(x_lower, x_upper, simulation_box_);
}



void ProbeVolume_Sphere::calculateIndicator(
	const CommonTypes::Real3& x,
	double& h_v, double& htilde_v, Real3& dhtilde_v_dx, RegionEnum& region
) const
{
	// Apply the minimum image convention to the vector between the sphere's
	// center and the particle's location
	double dist = 0.0, dist_sq;
	Real3  x_center_i;
	simulation_box_.calculateDistance(center_, x, x_center_i, dist_sq);
	dist = sqrt(dist_sq);

	// Evaluate the switching function
	double dhtilde_r_dr;
	switching_function_r_shell_.calculate( dist, need_derivatives_, h_v, htilde_v, dhtilde_r_dr, region );

	if ( need_derivatives_ ) {
		if ( htilde_v > 0.0 and dist > 0.0 ) {
			// Common prefactor: dhtilde(r)/dr * 1/r where r = norm(x_i - x_center) 
			dhtilde_r_dr /= dist;
			for ( int d=0; d<DIM_; ++d ) { dhtilde_v_dx[d] = dhtilde_r_dr*x_center_i[d]; }
		}
		else {
			for ( int d=0; d<DIM_; ++d ) { dhtilde_v_dx[d] = 0.0; }
		}
	}

	return;
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume_Sphere::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Probe_volume sphere\n"
	   << prepend_string << "  center = {" 
	      << center_[0] << ", " << center_[1] << ", " << center_[2] << "} [nm]\n"
	   << prepend_string << "  r_min = " << r_min_ << " [nm]\n"
	   << prepend_string << "  r_max = " << r_max_ << " [nm]\n"
	   << getSharedAttributesSummary(prepend_string + "  ");

	return ss.str();
}
