#include "SwitchingFunction_Indus.h"


SwitchingFunction_Indus::SwitchingFunction_Indus(const double sigma, const double alpha_c):
	SwitchingFunction(),
	sigma_(sigma), alpha_c_(alpha_c)
{
	this->setBaseClassParameters(sigma, alpha_c);
}


void SwitchingFunction_Indus::setBaseClassParameters(const double sigma, const double alpha_c)
{
	sigma_  = sigma;
	alpha_c_ = alpha_c;
	if ( sigma_ > 0.0 and alpha_c_ > 0.0 ) {
		two_sigma_sq_   = 2.0*sigma_*sigma_,
		sqrt_two_sigma_ = sqrt(2.0)*sigma_;

		k_   = sqrt(PI_*two_sigma_sq_)*erf(alpha_c_/sqrt_two_sigma_) 
					 - 2.0*alpha_c_*exp(-alpha_c_*alpha_c_/two_sigma_sq_);
		k_1_ = sqrt(0.5*PI_*sigma_*sigma_)/k_;
		k_2_ = exp(-alpha_c_*alpha_c_/two_sigma_sq_)/k_;
	}
	else {
		// Interpret zero/negative values as turning off coarse-graining
		sigma_   = 0.0;
		alpha_c_ = 0.0;

		two_sigma_sq_   = 0.0;
		sqrt_two_sigma_ = 0.0;
		k_   = 0.0;
		k_1_ = 0.0;
		k_2_ = 0.0;
	}
}


// Compute phi(alpha), a Gaussian which has been truncated at alpha_c_
// and shifted down so that it's normalized
double SwitchingFunction_Indus::calcTruncatedGaussian(const double alpha) const
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
