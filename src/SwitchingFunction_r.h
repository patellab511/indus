/* SwitchingFunction_r.h
 *
 * ABOUT: Implements the INDUS-type switching function with one boundary, h_r(r)
 */

#pragma once
#ifndef SWITCHING_FUNCTION_R_H
#define SWITCHING_FUNCTION_R_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "SwitchingFunction.h"

class SwitchingFunction_r : public SwitchingFunction
{
 public:
	SwitchingFunction_r(
		const double r_max   = 0.35,
		const double sigma   = 0.01, 
		const double alpha_c = 0.02
	);

	// Compute h_r = h(r), htilde_r(r), and dhtilde_r/dr
	virtual void calculate(
		const double r,
		const bool want_derivative,
		// Output
		double& h_r,
		double& htilde_r,
		double& dhtilde_r_dr
	) const override;

	void set_r_max(const double r_max);

	virtual void setCoarseGrainingParameters(const double sigma, const double alpha_c) override;

 private:
	double r_max_;              // nominal boundary (before coarse-graining)
	double r_inner_, r_outer_;  // boundaries of coarse-grained region
};

#endif // SWITCHING_FUNCTION_R_H
