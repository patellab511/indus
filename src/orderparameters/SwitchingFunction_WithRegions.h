/* SwitchingFunction_WithRegions.h
 *
 * ABOUT: Abstract interface for switching functions that also consider 
 *        regions of interest (as identified by RegionEnum)
 */

#pragma once
#ifndef SWITCHING_FUNCTION_WITH_REGIONS_H
#define SWITCHING_FUNCTION_WITH_REGIONS_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "Region.h"
#include "SwitchingFunction.h"

class SwitchingFunction_WithRegions : public SwitchingFunction
{
 public:
	using RegionEnum = Region::RegionEnum;

	SwitchingFunction_WithRegions():
		SwitchingFunction()
	{};

	// Calculate indicator function (and optionally its derivative)
	// - Calls the version below which returns a RegionEnum (the RegionEnum is discarded)
	virtual void calculate(
		const double x,
		const bool want_derivative,
		// Output
		double& h_x,          // value of h(x) before coarse-graining
		double& htilde_x,     // htilde_x(x) = coarse-grained h(x)
		double& dhtilde_x_dx  // Derivative of htilde_x(x)
	) const override;

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
	) const = 0;
};

inline
void SwitchingFunction_WithRegions::calculate(
	const double x, const bool want_derivative,
	double& h_x, double& htilde_x, double& dhtilde_x_dx
) const 
{
	RegionEnum tmp;
	this->calculate(x, want_derivative, h_x, htilde_x, dhtilde_x_dx, tmp);
}
#endif // SWITCHING_FUNCTION_WITH_REGIONS_H
