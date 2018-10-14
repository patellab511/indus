#include "SwitchingFunction.h"


SwitchingFunction::SwitchingFunction(const double sigma, const double alpha_c)
 : sigma_(sigma), alpha_c_(alpha_c)
{
	this->setBaseClassParameters(sigma, alpha_c);
}


void SwitchingFunction::setBaseClassParameters(const double sigma, const double alpha_c)
{
	sigma_  = sigma;
	alpha_c_ = alpha_c;

	two_sigma_sq_   = 2.0*sigma_*sigma_,
	sqrt_two_sigma_ = sqrt(2.0)*sigma_;

	k_   = sqrt(PI_*two_sigma_sq_)*erf(alpha_c_/sqrt_two_sigma_) 
	       - 2.0*alpha_c_*exp(-alpha_c_*alpha_c_/two_sigma_sq_);
	k_1_ = sqrt(0.5*PI_*sigma_*sigma_)/k_;
	k_2_ = exp(-alpha_c_*alpha_c_/two_sigma_sq_)/k_;
}


// Compute phi(alpha), a Gaussian which has been truncated at alpha_c_
// and shifted down so that it's normalized
double SwitchingFunction::calcTruncatedGaussian(const double alpha) const
{
	// Note: if sigma = 0.0, the derivative at alpha=0.0 is actually undefined
	// - Here, I choose to set it to zero
	if ( std::abs(alpha) < alpha_c_ and sigma_ > 0.0 ) {
		return exp(-alpha*alpha/two_sigma_sq_)/k_ - k_2_;
	}
	else {
		return 0.0;
	}
}
