/* Bias.h
 *
 * ABOUT: Abstract base class for implementing different probe volume geometries
 * NOTES: 
 *   U_bias(x) is the sum of the following terms:
 *
 *     U_harmonic(x)   = 0.5*kappa*(x - x_star)^2 
 *     U_linear        = phi*x + constant
 *     U_left_well(x)  = k_left*(x - x_left)*Theta(x_left - x)
 *     U_right_well(x) = k_right*(x - x_right)*Theta(x - x_right)
 *   
 *  Theta(x) is the unit step function, which is 1 for x>=0 and 0 for x<0
 */

#pragma once

#ifndef BIAS_H
#define BIAS_H

// Standard headers
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

// Project headers
#include "CommonTypes.h"
#include "InputParser.h"
#include "OpenMP.h"
#include "Potential.h"
#include "SimulationState.h"
#include "StringTools.h"

class Ramp;

class Bias
{
 public:
  static const int DIM_  = CommonTypes::DIM_;  // Dimensionality of simulation 
  using rvec      = CommonTypes::Rvec;
  using Real3     = CommonTypes::Real3;
	using Matrix    = CommonTypes::Matrix3x3;
	using RvecArray = CommonTypes::RvecArray;
  using Range     = CommonTypes::Range;

	Bias(
		const ParameterPack& input_parameter_pack,  // contains parameters
		const SimulationState& simulation_state     // need time for time-dependent ramps
	);

	// Calculate u_bias(x), and optionally du_bias(x)/dx
	// - Returns u_bias(x)
	double calculate(
		const double x,
		const bool   want_derivative,  // whether to compute du_bias_dx
		// Output
		double&      u_bias,
		double&      du_bias_dx
	) const;

	// Computes the bias, u_bias(x), and optionally its derivatives, and applies
	// the chain rule to produce the derivatives of the bias wrt. the atom positions
	void applyBias(
		const double              x,                 // order parameter, x
		const bool                want_derivatives,
		const std::vector<Real3>& derivatives_x,     // dx/dr_i for a set of positions, r_i
		const Matrix&             sum_r_cross_dx_dr,
		// Output
		double&             u_bias_x,              // u_bias(x), where x = x(r^N)
		std::vector<Real3>& derivatives_u_bias_x,  // du_bias_(x)/dr_i (length: num_atoms_with_bias)
		Matrix&             virial_u_bias_x        // Contribution to virial due to total bias
	) const;

	// Computes the bias, u_bias(x), and optionally its derivatives, and applies
	// the chain rule to add the contribution from the bias to the given arrays
	// of derivatives of x wrt. atom positions
	void accumulateBias(
		const double              x,                 // order parameter, x = x(r^N) = x({r_i})
		const bool                want_derivatives,  // (of U_bias_total)
		const std::vector<Real3>& derivatives_x,     // dx/dr_i for each position, r_i
		const Matrix&             sum_r_cross_dx_dr,
		// Output
		double&             u_bias_x,                  // U_bias(x)
		double&             u_bias_total,              // U_bias_total = running sum of all biases
		std::vector<Real3>& derivatives_u_bias_total,  // dU_bias_total/dr_i (length: num_atoms_with_bias)
		Matrix&             virial_u_bias_total        // Contribution to virial due to total bias
	) const;

	void setParameters(const ParameterPack& input_parameter_pack);

	//----- Get functions -----//

	// Returns a big string with a summary of the input
	// - Only reports on terms which are not guaranteed to be zero
	//   - e.g. harmonic term of not printed if kappa_ = 0.0
	std::string getInputSummary(
		const std::string& prepend_string = "# " // prepend to each line
	) const;

 private:
	const std::string input_string_;

	const SimulationState& simulation_state_;  // owned by driver

	std::vector<std::unique_ptr<Potential>> potential_ptrs_;

	// If true (and ramp is active), scale the "location" of the bias with time
	// Else scale the bias as a prefactor
	// - e.g. harmonic bias with ramp coefficient c(t)
	//   (a) True:   u_bias = 0.5*kappa*( x - c(t)*xstar )^2
	//   (b) False:  u_bias = c(t) * 0.5*kappa*( x - xstar )^2
	std::unique_ptr<Ramp> bias_ramp_ptr_ = nullptr;
};


// For ramping a bias up/down
// - i.e. scaling a parameter that changes from parameter_range_ramp_[0] 
//   to parameter_range_ramp_[1] over the time interval [t_range_ramp_[0], to t_range_ramp_[1]]
class Ramp
{
 public:
	using Range = CommonTypes::Range;

	Ramp(
		const ParameterPack& input_parameter_pack,  // contains parameters
		const SimulationState& simulation_state     // need time for time-dependent ramps
	);
	Ramp(
		const Range& t_range,
		const Range& parameter_range,
		const SimulationState& simulation_state     // need time for time-dependent ramps
	);
	virtual ~Ramp() {}

	virtual double calculate() const = 0;

	// Get functions
	const Range& get_t_range()         const { return t_range_;         }
	const Range& get_parameter_range() const { return parameter_range_; }

 protected:
	const SimulationState& simulation_state_;  // owned by driver

	Range t_range_, parameter_range_;
};


class LinearRamp : public Ramp
{
 public:
	LinearRamp(
		const ParameterPack& input_parameter_pack,  // contains parameters
		const SimulationState& simulation_state     // need time for time-dependent ramps
	);
	LinearRamp(
		const Range& t_range,
		const Range& parameter_range,
		const SimulationState& simulation_state     // need time for time-dependent ramps
	);

	virtual double calculate() const override;

 protected:
	double slope_;

	double calculate_slope() const;
};

#endif // BIAS_H	
