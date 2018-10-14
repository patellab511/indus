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
 *  Theta(x) is the unit step function which is 1 for x>=0 and 0 for x<0
 */

#pragma once

#ifndef BIAS_H
#define BIAS_H

// Standard headers
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

// Project headers
#include "AtomGroup.h"
#include "CommonTypes.h"
#include "InputParser.h"
#include "SimulationBox.h"
#include "StringTools.h"


class Bias
{
 public:
  static const int DIM_  = CommonTypes::DIM_;  // Dimensionality of simulation 
  using rvec      = CommonTypes::Rvec;
  using Real3     = CommonTypes::Real3;
	using Matrix    = CommonTypes::Matrix3x3;
	using RvecArray = CommonTypes::RvecArray;

	Bias(
		const ParameterPack& input_parameter_pack,  // contains parameters
		const SimulationBox& simulation_box         // need to know about box for virial
	);

	void accumulateBias(
		const double              x,                 // order parameter, x = x(r^N) = x({r_i})
		const bool                want_derivatives,  // (of U_bias_total)
		const std::vector<Real3>& derivatives_x,     // dx/dr_i for each position, r_i
		const Matrix&             sum_r_cross_dx_dr,
		// Output
		double&             u_bias_x,                  // U_bias(x)
		double&             u_bias_total,              // U_bias_total = running sum of all biases
		std::vector<Real3>& derivatives_u_bias_total,  // dU_bias_total/dr_i [num_local_atoms_with_bias x 1]
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

	const SimulationBox& simulation_box_; // owned by driver

	//----- Biasing parameters

	// Harmonic bias
	double x_star_, kappa_;

	// Linear bias
	double phi_, constant_;

	// Left and right one-sided harmonics
	double x_left_,  k_left_;
	double x_right_, k_right_;
};

#endif // BIAS_H	
