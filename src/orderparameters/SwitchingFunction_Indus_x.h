/* SwitchingFunction_Indus_x.h
 *
 * ABOUT: Implements the INDUS-type switching function with two boundaries, h_x(x)
 */

#pragma once
#ifndef SWITCHING_FUNCTION_INDUS_X_H
#define SWITCHING_FUNCTION_INDUS_X_H

// Standard headers
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "SwitchingFunction_Indus.h"

class SwitchingFunction_Indus_x : public SwitchingFunction_Indus
{
 public:
	SwitchingFunction_Indus_x(
		const double x_min   = 0.0,
		const double x_max   = 1.0,
		const double sigma   = 0.01, 
		const double alpha_c = 0.02
	);

	// Compute h_x = h(x), htilde_x, and dhtilde_x/dx
	virtual void calculate(
		const double x, 
		const bool want_derivative,
		// Output
		double& h_x,
		double& htilde_x,
		double& dhtilde_x_dx
	) const override;

	virtual void setLimits(const double x_min, const double x_max);
	void setLimits(const std::array<double,2>& x_limits) {
		setLimits(x_limits[0], x_limits[1]);
	}

	void getLimits(double& x_min, double& x_max) const {
		x_min = x_min_;
		x_max = x_max_;
	}
	std::array<double,2> getLimits() const {
		return {{ x_min_, x_max_ }};
	}
	void getOuterLimits(double& x_lowest, double& x_highest) const {
		x_lowest  = x_lowest_;
		x_highest = x_highest_;
	}
	void getInnerLimits(double& x_min_inner, double& x_max_inner) const {
		x_min_inner = x_min_ + alpha_c_;
		x_max_inner = x_max_ - alpha_c_;
	}

	virtual void setCoarseGrainingParameters(const double sigma, const double alpha_c) override;

 protected:
	double x_min_,    x_max_;      // nominal bounds (before coarse-graining)
	double x_lowest_, x_highest_;  // (outer) coarse-grained bounds
};

#endif // SWITCHING_FUNCTION_INDUS_X_H
