// SwitchingFunction_Indus_x.cpp
#include "SwitchingFunction_Indus_x.h"


SwitchingFunction_Indus_x::SwitchingFunction_Indus_x( 
	const double x_min, const double x_max, const double sigma, const double alpha_c 
):
	SwitchingFunction_Indus(sigma, alpha_c)
{
	this->setLimits(x_min, x_max);
}


void SwitchingFunction_Indus_x::setLimits(const double x_min, const double x_max)
{
	x_min_ = x_min;
	x_max_ = x_max;

	x_lowest_  = x_min_ - alpha_c_;
	x_highest_ = x_max_ + alpha_c_;

	/*
	if ( x_highest_ < x_lowest_ ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  The distance between the boundaries is less than 2*alpha_c. "
		         << "  This is not currently supported.\n"
		       << "  Input: x_min = " << x_min << ", x_max = " 
		         << x_max << ", alpha_c = " << alpha_c_ << "\n";
		throw std::runtime_error( err_ss.str() );
	}
	*/
}


void SwitchingFunction_Indus_x::setCoarseGrainingParameters(const double sigma, const double alpha_c)
{
	setBaseClassParameters(sigma, alpha_c);

	// Update constants associated with x_min and x_max
	setLimits(x_min_, x_max_);
}


// Compute htilde_x(x) (and dhtilde_x/dx) 
void SwitchingFunction_Indus_x::calculate(
	const double x, const bool want_derivative,
	double& h_x, double& htilde_x, double& dhtilde_x_dx
) const
{
	if ( x < x_lowest_ or x > x_highest_ ) {
		// Completely outside region of interest
		h_x      = 0.0;
		htilde_x = 0.0;
		if ( want_derivative ) { dhtilde_x_dx = 0.0; }
	}
	else {
		// Somewhere in region of interest

		if ( x >= x_min_ and x <= x_max_ ) { h_x = 1.0; }
		else                               { h_x = 0.0; }

		if ( sigma_ > 0.0 ) {
			double dx_min = x - x_min_;
			double dx_max = x - x_max_;
			if ( std::abs(dx_min) <= alpha_c_ ) { // Inside left buffer region
				htilde_x = 0.5 + k_1_*erf(dx_min/sqrt_two_sigma_) - k_2_*dx_min;
				if ( want_derivative ) { dhtilde_x_dx = calcTruncatedGaussian(dx_min); }
			}
			else if ( std::abs(dx_max) <= alpha_c_ ) { // Inside right buffer region
				htilde_x = 0.5 + k_1_*erf(-dx_max/sqrt_two_sigma_) + k_2_*dx_max;
				if ( want_derivative ) { dhtilde_x_dx = -calcTruncatedGaussian(dx_max); }
			}
			else { // Completely inside volume of interest
				htilde_x = 1.0;
				if ( want_derivative ) { dhtilde_x_dx = 0.0; }
			}
		}
		else {
			// No coarse-graining
			// - Note: the derivative at x=x_min_ and x_max_ is actually undefined
			htilde_x = h_x;
			if ( want_derivative ) { dhtilde_x_dx = 0.0; }
		}
	}
}
