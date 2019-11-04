/* SwitchingFunction_WithRegions_r.h
 *
 * ABOUT: Implements a INDUS-type switching function for a region with
 *        two boundaries, where regions are also checked
 */

#pragma once
#ifndef SWITCHING_FUNCTION_WITH_REGIONS_R_H
#define SWITCHING_FUNCTION_WITH_REGIONS_R_H

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
#include "SwitchingFunction_Indus_r.h"

class SwitchingFunction_WithRegions_r : public SwitchingFunction_WithRegions
{
 public:
	SwitchingFunction_WithRegions_r(
		const double r     = 1.0,
		const double sigma = 0.01, 
		const double alpha_c = 0.02,
		const double width_shell_1 = 0.0, 
		const double width_shell_2 = 0.0,
		const bool   is_inverted = false
	);

	// Calculate indicator function (and optionally its derivative), and also
	// return which region 'r' is in
	virtual void calculate(
		const double r,
		const bool want_derivative,
		// Output
		double& h_r,          // value of h(r) before coarse-graining
		double& htilde_r,     // htilde_r(r) = coarse-grained h(r)
		double& dhtilde_r_dr, // Derivative of htilde_r(r)
		RegionEnum& region    // Type of region which contains 'r'
	) const override;

	// 'Core' set function that sets all relevant member variables
	void setBoundaries(
		const double r_max,
		const double width_shell_1,
		const double width_shell_2
	);

	// Change r_min/max without changing shells
	void set_r_max(const double r_max) {
		setBoundaries(r_max, width_shell_1_, width_shell_2_);
	}
	// Change shell widths without changing r_min/max
	void setShellWidths(const double width_shell_1, const double width_shell_2) {
		double r_max = switch_r_.get_r_max();
		setBoundaries(r_max, width_shell_1, width_shell_2);
	}
	// Update coarse-graining parameters
	void setCoarseGrainingParameters(const double sigma, const double alpha_c) {
		switch_r_.setCoarseGrainingParameters(sigma, alpha_c);
		double r_max = switch_r_.get_r_max();
		setBoundaries(r_max, width_shell_1_, width_shell_2_);
	}
	void setInversion(const bool is_inverted) {
		is_inverted_ = is_inverted;
		double r_max = switch_r_.get_r_max();
		setBoundaries(r_max, width_shell_1_, width_shell_2_);
	}

	void getShellWidths(double& width_shell_1, double& width_shell_2) const {
		width_shell_1 = width_shell_1_;
		width_shell_2 = width_shell_2_;
	}
	double get_r_max() const { 
		return switch_r_.get_r_max(); 
	}
	void getCoarseGrainingParameters(double& sigma, double& alpha_c) const {
		switch_r_.getCoarseGrainingParameters(sigma, alpha_c);
	}
	bool is_inverted() const { return is_inverted_; }

 protected:
	// Switching function defining the primary region (nominal and coarse-grained)
	SwitchingFunction_Indus_r switch_r_;

	// Size of shells around the coarse-grained primary region
	double width_shell_1_, width_shell_2_;
	bool nontrivial_shells_;  // if true, shells must be checked

	// Furthest limits for shells
	double r_max_shell_1_;
	double r_max_shell_2_;

	// Whether to "invert" the switching function, making it cover
	// r = r_max to r = +infinity instead of 0 to r_max
	bool is_inverted_;
};

#endif // SWITCHING_FUNCTION_WITH_REGIONS_R_H
