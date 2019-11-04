/* NeighborSphere.h
 *
 * ABOUT: 
 */

#pragma once
#ifndef NEIGHBOR_SPHERE_H
#define NEIGHBOR_SPHERE_H

#ifndef INDUS_STANDALONE_MODE
#define NEIGHBOR_SPHERE_PLUMED_MODE
#endif

// Standard headers
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>    // unique_ptr
#include <sstream>
#include <string>
#include <sstream>

// Project headers
#include "CommonTypes.h"
#include "InputParser.h"
#include "StringTools.h"
#include "SwitchingFunction_Indus_r.h"

class NeighborSphere
{
 public:
	// Typedefs
  static const int DIM_  = CommonTypes::DIM_;  // Dimensionality of simulation 
  using Real3 = CommonTypes::Real3;

	NeighborSphere(
		const double radius, 
		const double sigma, 
		const double alpha_c,
		const bool   need_derivatives
	);

	NeighborSphere(
		const ParameterPack& input_parameter_pack, 
		const bool need_derivatives
	);

	// Computes h^[nn(i)](r_j) = h_r(r_i,j) and its coarse-grained version,
	// htilde_nn(r) (and derivatives of htilde w.r.t. position r_j)
	void calculateIndicator(
		const Real3& r_ij,     // Minimim image vector from center "i" to neighbor "j"
		const double dist_ij,
		// Output
		double& h_nn,
		double& htilde_nn,
		Real3& dhtilde_nn_dr_j
	) const;

	// Get/set the radius
	double getRadius() const { return r_; }
	void setRadius(const double r) { 
		r_ = r;
		switching_function_r_.set_r_max(r_);
	}

	// Returns a summary of settings
	std::string getInputSummary(
		const std::string& prefix = std::string("# ")
	) const;

	void setCoarseGrainingParameters(const double sigma, const double alpha_c) {
		sigma_   = sigma;
		alpha_c_ = alpha_c;
		switching_function_r_.setCoarseGrainingParameters(sigma_, alpha_c_);
	}

	void getParameters(double& r, double& sigma, double& alpha_c) const {
		r       = r_;
		sigma   = sigma_;
		alpha_c = alpha_c_;
	}

	// Get the radius at which the coarse-grained indicator function becomes zero
	double get_rtilde_max() const { return switching_function_r_.get_rtilde_max(); }

 private:
	double r_;       // Nominal radius (excluding buffer) [nm]
	double sigma_;   // Coarse-grained length scale [nm]
	double alpha_c_; // Width of Gaussian buffer (usually 2*sigma) [nm]

	SwitchingFunction_Indus_r switching_function_r_;

	bool need_derivatives_;

	static constexpr double PI_ = 3.14159265358979323846;
};

#endif // NEIGHBOR_SPHERE_H
