/* SwitchingFunction_x.h
 *
 * ABOUT: Implements the INDUS-type switching function with two boundaries, h_x(x)
 */

#pragma once
#ifndef SWITCHING_FUNCTION_X_H
#define SWITCHING_FUNCTION_X_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "SwitchingFunction.h"

class SwitchingFunction_x : public SwitchingFunction
{
 public:
	SwitchingFunction_x(
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

	void setLimits(const double x_min, const double x_max);

	void getLimits(double& x_min, double& x_max) const {
		x_min = x_min_;
		x_max = x_max_;
	}

	virtual void setCoarseGrainingParameters(const double sigma, const double alpha_c) override;

 private:
	double x_min_,    x_max_;      // nominal bounds (before coarse-graining)
	double x_lowest_, x_highest_;  // coarse-grained bounds
};

#endif // SWITCHING_FUNCTION_X_H
