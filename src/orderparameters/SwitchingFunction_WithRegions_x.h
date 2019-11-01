/* SwitchingFunction_WithRegions_x.h
 *
 * ABOUT: Implements a INDUS-type switching function for a region with
 *        two boundaries, where regions are also checked
 */

#pragma once
#ifndef SWITCHING_FUNCTION_WITH_REGIONS_X_H
#define SWITCHING_FUNCTION_WITH_REGIONS_X_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "Region.h"
#include "SwitchingFunction_WithRegions.h"
#include "SwitchingFunction_Indus_x.h"

class SwitchingFunction_WithRegions_x : public SwitchingFunction_WithRegions
{
 public:
	SwitchingFunction_WithRegions_x(
		const double x_min = 0.0, 
		const double x_max = 1.0, 
		const double sigma = 0.01, 
		const double alpha_c = 0.02,
		const double width_shell_1 = 0.0, 
		const double width_shell_2 = 0.0,
		const bool   is_inverted = false
	);

	// Calculate indicator function (and optionally its derivative), and also
	// return which region 'x' is in
	virtual void calculate(
		const double x,
		const bool want_derivative,
		// Output
		double& h_x,          // value of h(x) before coarse-graining
		double& htilde_x,     // htilde_x(x) = coarse-grained h(x)
		double& dhtilde_x_dx, // Derivative of htilde_x(x)
		RegionEnum& region    // Type of region which contains 'x'
	) const override;

	// 'Core' set function that sets all relevant member variables
	void setBoundaries(
		const double x_min, 
		const double x_max,
		const double width_shell_1,
		const double width_shell_2
	);

	// Change x_min/max without changing shells
	void setLimits(const double x_min, const double x_max) {
		setBoundaries(x_min, x_max, width_shell_1_, width_shell_2_);
	}
	// Change shell widths without changing x_min/max
	void setShellWidths(const double width_shell_1, const double width_shell_2) {
		double x_min, x_max;
		switch_x_.getLimits(x_min, x_max);
		setBoundaries(x_min, x_max, width_shell_1, width_shell_2);
	}
	// Update coarse-graining parameters
	void setCoarseGrainingParameters(const double sigma, const double alpha_c) {
		switch_x_.setCoarseGrainingParameters(sigma, alpha_c);
		double x_min, x_max;
		switch_x_.getLimits(x_min, x_max);
		setBoundaries(x_min, x_max, width_shell_1_, width_shell_2_);
	}
	void setInversion(const bool is_inverted) {
		is_inverted_ = is_inverted;
		double x_min, x_max;
		switch_x_.getLimits(x_min, x_max);
		setBoundaries(x_min, x_max, width_shell_1_, width_shell_2_);
	}

	void getShellWidths(double& width_shell_1, double& width_shell_2) const {
		width_shell_1 = width_shell_1_;
		width_shell_2 = width_shell_2_;
	}
	void getLimits(double& x_min, double& x_max) const {
		switch_x_.getLimits(x_min, x_max);
	}
	void getOuterLimits(double& x_min_outer, double& x_max_outer) const {
		switch_x_.getOuterLimits(x_min_outer, x_max_outer);
	}
	void getCoarseGrainingParameters(double& sigma, double& alpha_c) const {
		switch_x_.getCoarseGrainingParameters(sigma, alpha_c);
	}
	bool is_inverted() const { return is_inverted_; }

 protected:
	// Switching function defining the primary region (nominal and coarse-grained)
	SwitchingFunction_Indus_x switch_x_;

	// Size of shells around the coarse-grained primary region
	double width_shell_1_, width_shell_2_;
	bool nontrivial_shells_;  // if true, shells must be checked

	// Furthest limits (upper/lower) for shells
	double x_min_shell_1_, x_max_shell_1_;
	double x_min_shell_2_, x_max_shell_2_;

	// Whether to "invert" the switching function, making it cover
	// x <= x_min and x >= x_max instead of x_min <= x <= x_max
	bool is_inverted_;
};

#endif // SWITCHING_FUNCTION_WITH_REGIONS_X_H
