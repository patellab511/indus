/* SwitchingFunction.h
 *
 * ABOUT: Abstract interface for switching functions
 */

#pragma once
#ifndef SWITCHING_FUNCTION_H
#define SWITCHING_FUNCTION_H

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>    // unique_ptr
#include <sstream>
#include <stdexcept>
#include <string>

#include "utils.h"

class SwitchingFunction
{
 public:
	SwitchingFunction() {};
	virtual ~SwitchingFunction() {};

	// Calculate indicator function, and optionally its derivative
	virtual void calculate(
		const double x,
		const bool want_derivative,
		// Output
		double& h_x,          // value of h(x) before coarse-graining
		double& htilde_x,     // htilde_x(x) = coarse-grained h(x)
		double& dhtilde_x_dx  // Derivative of htilde_x(x)
	) const = 0;
};

#endif // SWITCHING_FUNCTION_H
