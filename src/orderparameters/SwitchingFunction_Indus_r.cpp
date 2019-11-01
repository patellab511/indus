// SwitchingFunction_Indus_r
#include "SwitchingFunction_Indus_r.h"


SwitchingFunction_Indus_r::SwitchingFunction_Indus_r( const double r_max, 
                                          const double sigma, const double alpha_c )
	: SwitchingFunction_Indus(sigma, alpha_c),
    r_max_(r_max)
{
	this->set_r_max(r_max_);
}


void SwitchingFunction_Indus_r::set_r_max(const double r_max)
{
	r_max_ = r_max;

	r_inner_ = r_max_ - alpha_c_;
	r_outer_ = r_max_ + alpha_c_;

	/*
	// Trust that users know what they're doing
	if ( r_inner_ < 0.0 ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  The inner boundary of the coarse-grained region is less than zero. "
		         << "  This is not currently supported.\n"
		       << "  Input: r_max = " << r_max_ << ", alpha_c = " << alpha_c_
		         << ", r_inner = " << r_inner_ << "\n";
		throw std::runtime_error( err_ss.str() );
	}
	*/
}


void SwitchingFunction_Indus_r::setCoarseGrainingParameters(const double sigma, const double alpha_c)
{
	setBaseClassParameters(sigma, alpha_c);

	// Update constants associated with r_max
	set_r_max(r_max_);
}


// Compute h(r) (and dh/dr)
void SwitchingFunction_Indus_r::calculate(
	const double r, const bool want_derivative,
	double& h_r, double& htilde_r, double& dhtilde_r_dr
) const
{
	if ( r > r_outer_ ) {
		// Fully outside region
		h_r      = 0.0;
		htilde_r = 0.0;
		if ( want_derivative ) { dhtilde_r_dr = 0.0; }
	}
	else if ( r <= r_inner_ ) {
		// Fully inside region
		h_r      = 1.0;
		htilde_r = 1.0;
		if ( want_derivative ) { dhtilde_r_dr = 0.0; }
	}
	else {
		// In the Gaussian buffer zone
		if ( r <= r_max_ ) { h_r = 1.0; }
		else               { h_r = 0.0; } 

		if ( sigma_ > 0.0 ) {
			double delta_r = r_max_ - r;
			htilde_r = 0.5 + k_1_*erf(delta_r/sqrt_two_sigma_) - k_2_*delta_r;
			if ( want_derivative ) { dhtilde_r_dr = -calcTruncatedGaussian(delta_r); }
		}
		else {
			// No coarse-graining
			// - Note that the derivative at r=r_max is actually undefined
			htilde_r = h_r;
			if ( want_derivative ) { dhtilde_r_dr = 0.0; }
		}
	}
}
