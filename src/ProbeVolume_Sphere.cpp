/* ProbeVolume_Sphere.cpp
 *
 * ABOUT: 
 *
 */

#include "ProbeVolume_Sphere.h"

// Register this probe volume
namespace ProbeVolumeRegistry {
static const RegisterInFactory<
		ProbeVolume, 
		ProbeVolume_Sphere, 
		const std::string, 
		ProbeVolume::ProbeVolumeInputPack
	> register_ProbeVolume_Sphere("sphere");
}


ProbeVolume_Sphere::ProbeVolume_Sphere(ProbeVolume::ProbeVolumeInputPack& input_pack)
 : ProbeVolume(input_pack),
   r_(1.0),
   switching_function_r_(r_, sigma_, alpha_c_)
{
	const ParameterPack& input_parameter_pack = input_pack.input_parameter_pack;
	using KeyType = ParameterPack::KeyType;

	// Coarse-graining parameters
	input_parameter_pack.readNumber("sigma",   KeyType::Optional, sigma_);
	input_parameter_pack.readNumber("alpha_c", KeyType::Optional, alpha_c_);
	this->setCoarseGrainingParameters(sigma_, alpha_c_);

	// Sphere geometry
	input_parameter_pack.readNumber("radius", KeyType::Required, r_);
	input_parameter_pack.readArray("center", KeyType::Required, center_);
	this->setGeometry(center_, r_);
}


void ProbeVolume_Sphere::setCoarseGrainingParameters(const double sigma, const double alpha_c)
{
	sigma_   = sigma;
	alpha_c_ = alpha_c_;

	// Update other member variables accordingly
	switching_function_r_.setCoarseGrainingParameters(sigma_, alpha_c_);
	this->setGeometry(center_, r_);
}


void ProbeVolume_Sphere::setShellWidths(const double width_shell_1, const double width_shell_2)
{
	width_shell_1_  = width_shell_1;
	width_shell_2_  = width_shell_2;

	// Update probe sphere
	setGeometry(center_, r_);
}


void ProbeVolume_Sphere::setGeometry(const Real3& center, const double r)
{
	// Location of center: { x, y, z }
	center_ = center;

	// Radius
	r_    = r;
	r_sq_ = r*r;
	switching_function_r_.set_r_max(r_);

	// Effective radius, accounting for buffer region
	r_eff_    = r_ + alpha_c_;
	r_eff_sq_ = r_eff_*r_eff_;

	// Shells (with buffer)
	r_shell_1_ = r_eff_     + width_shell_1_;   r_shell_1_sq_ = r_shell_1_*r_shell_1_;
  r_shell_2_ = r_shell_1_ + width_shell_2_;   r_shell_2_sq_ = r_shell_2_*r_shell_2_;
}


void ProbeVolume_Sphere::updateUsingSimulationState() 
{ 
	// Nothing to do for this probe volume geometry
	return; 
}


bool ProbeVolume_Sphere::isInProbeVolume( const CommonTypes::Real3& x,
                                          double& h_v, double& htilde_v, Real3& dhtilde_v_dx,
	                                        bool& is_in_shell_1, bool& is_in_shell_2 ) const 
{
	// Apply the minimum image convention to the vector between the sphere's
	// center and the particle's location
	double dist = 0.0, dist_sq;
	Real3  x_center_i;
	simulation_box_.calculateDistance(center_, x, x_center_i, dist_sq);

	bool   is_in_probe_volume;
	double dhtilde_r_dr;

	is_in_probe_volume = false;
	is_in_shell_1 = false;
	is_in_shell_2 = false;

	// Check for complete exclusion first (true for most particles if the probe volume
	// is significantly smaller than the simulation box)
	if (dist_sq > r_shell_2_sq_) {
		// Far away from probe volume
		h_v      = 0.0;
		htilde_v = 0.0;
	}
	else {
		// Definitely involved with the probe volume in some way
		dist = sqrt(dist_sq);
		switching_function_r_.calculate( dist, need_derivatives_, h_v, htilde_v, dhtilde_r_dr );

		if ( htilde_v > 0.0 ) {
			// Inside at least the coarse-grained probe volume.
			// If h_v = 1.0, it's also inside the nominal probe volume (not coarse-grained)
			if ( h_v == 1.0 ) { 
				is_in_probe_volume = true;  
			}
		}
		else if (dist_sq <= r_shell_1_sq_) {
			is_in_shell_1 = true;
		}
		else { // (dist_sq <= r_shell_2_sq_)
			is_in_shell_2 = true;
		}
	}

	if ( need_derivatives_ ) {
		if ( htilde_v > 0.0 && dist > 0.0 ) {
			// Common prefactor: dhtilde(r)/dr * 1/r where r = norm(x_i - x_center) 
			dhtilde_r_dr /= dist;
			for ( int d=0; d<DIM_; ++d ) { dhtilde_v_dx[d] = dhtilde_r_dr*x_center_i[d]; }
		}
		else {
			for ( int d=0; d<DIM_; ++d ) { dhtilde_v_dx[d] = 0.0; }
		}
	}

	return is_in_probe_volume;
}


// Returns a string with complete, formatted information about the probe volume's
// location and geometry which can be included in output files. 
std::string ProbeVolume_Sphere::getInputSummary() const
{
	std::stringstream ss;
	ss << "# Probe_volume sphere\n"
	   << "#   center = {" 
	      << center_[0] << ", " << center_[1] << ", " << center_[2] << "} [nm]\n"
	   << "#   radius = " << r_ << " [nm]\n"
	   << getSharedAttributesSummary();

	return ss.str();
}
