/* SwitchingFunction_Indus_r.h
 *
 * ABOUT: Implements the INDUS-type switching function with one boundary, h_r(r)
 */

#pragma once
#ifndef SWITCHING_FUNCTION_INDUS_R_H
#define SWITCHING_FUNCTION_INDUS_R_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "SwitchingFunction_Indus.h"

class SwitchingFunction_Indus_r : public SwitchingFunction_Indus
{
 public:
	SwitchingFunction_Indus_r(
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

	double get_r_max() const { return r_max_; }

	// Get the radius at which the coarse-grained indicator function becomes zero
	double get_rtilde_max() const { return r_outer_; }

	// r_max +/- alpha_c
	double get_r_outer() const { return r_outer_; }
	double get_r_inner() const { return r_inner_; }

	virtual void setCoarseGrainingParameters(const double sigma, const double alpha_c) override;

 private:
	double r_max_;              // nominal boundary (before coarse-graining)
	double r_inner_, r_outer_;  // boundaries of coarse-grained region
};

#endif // SWITCHING_FUNCTION_INDUS_R_H
