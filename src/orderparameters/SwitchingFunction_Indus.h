/* SwitchingFunction_Indus.h
 *
 * ABOUT: Abstract base class for simple INDUS-type switching functions
 */

#pragma once
#ifndef SWITCHING_FUNCTION_INDUS_H
#define SWITCHING_FUNCTION_INDUS_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "SwitchingFunction.h"

class SwitchingFunction_Indus : public SwitchingFunction
{
 public:
	SwitchingFunction_Indus(
		const double sigma   = 0.01, 
		const double alpha_c = 0.02
	);

	/*
	// Calculate indicator function, and optionally its derivative
	virtual void calculate(
		const double x,
		const bool want_derivative,
		// Output
		double& h_x,          // value of h(x) before coarse-graining
		double& htilde_x,     // htilde_x(x) = coarse-grained h(x)
		double& dhtilde_x_dx  // Derivative of htilde_x(x)
	) const = 0;
	*/

	// Compute phi(alpha), a Gaussian with standard deviation sigma_ and mean 0.0
	// that has been truncated at alpha_c_ and shifted down so that it's normalized
	double calcTruncatedGaussian(const double alpha) const;

	virtual void setCoarseGrainingParameters(
		const double sigma, 
		const double alpha_c
	) = 0;

	void getCoarseGrainingParameters(double& sigma, double& alpha_c) const {
		sigma   = sigma_;
		alpha_c = alpha_c_;
	}

 protected:
	// Coarse-graining parameters
	double sigma_;          // Coarse-grained length scale [nm]
	double alpha_c_;        // Width of Gaussian buffer (usually 2*sigma) [nm]

	// Commonly-used constants that depend only on the 
	// coarse-graining parameters
	double k_, k_1_, k_2_;  // Lumped constants
	double two_sigma_sq_;   // 2*sigma^2
	double sqrt_two_sigma_; // sqrt(2)*sigma

	static constexpr double PI_ = 3.14159265358979323846;

	// Compute the commonly-used constants used above
	// - When updating the coarse-graining parameters in an implementation
	//   of setCoarseGrainingParameters(), this function *must* always be called first!
	void setBaseClassParameters(const double sigma, const double alpha_c);
};

#endif // SWITCHING_FUNCTION_INDUS_H
