/* ProbeVolume_Cylinder.h
 *
 * ABOUT:
 *  - Does not work if the cylinder crosses the edges of the box
 * SYNTAX: 
 *   ProbeVolume = {
 *     # Geometry (required)
 *     type  = cylinder
 *     axis  = <x|y|z>    # specify the principle axis of the cylinder
 *     # Basic options (required)
 *     r_max = <r_max>
 *     axis_range = [ <min> <max> ]   # limits along cylinder axis
 *     base = [ <x> <y> <z> ]
 *     # Options to make the cylinder a cylindrical shell (optional)
 *     r_min = <r_min>
 *     axis_range_inner = [ <min> <max> ]   # inner limits along axis
 *     # Coarse-graining
 *     sigma   = <sigma>      # (default: 0.01 nm)
 *     alpha_c = <alpha>      # (default: 0.02 nm)
 *   }
 * DEVELOPMENT:
 *  - "zeta" is shorthand for the axial coordinate (may or may not be z!)
 */

#pragma once
#ifndef PROBE_VOLUME_CYLINDER_H
#define PROBE_VOLUME_CYLINDER_H

// Standard headers
#include <map>

// Core headers
#include "ProbeVolume.h"     // Parent class
#include "GenericFactory.h"  // Register as a ProbeVolume


class ProbeVolume_Cylinder : public ProbeVolume
{
 public:
	using Box   = ProbeVolume::Box;
	using Real3 = ProbeVolume::Real3;
	using Rvec  = ProbeVolume::Rvec;

	ProbeVolume_Cylinder(ProbeVolumeInputPack& input_pack);

	virtual void setGeometry() override;

	virtual void calculateIndicator(
		const Real3& x,
		// Output
		double& h_v,
		double& htilde_v, 
		Real3&  derivatives,
		RegionEnum& region
	) const override;

	std::string getInputSummary(const std::string& prepend_string) const override;

 protected:
	virtual BoundingBox constructBoundingBox() const override;

 private:
	// Cylinder geometry
	Real3  x_base_;   // Position of base: {x,y,z}
	Real3  center_;   // Center point: {x,y,z}
	double r_max_;    // (Outer) radius
	double h_;
	Range       zeta_range_;    // nominal bounds along cylinder axis
	int         axis_index_;    // index of cylinder axis (0=x, 1=y, 2=z)
	std::string axis_label_;

	// Switching functions for cylinder
	SwitchingFunction_WithRegions_r switching_function_r_;  // r_max
	SwitchingFunction_WithRegions_x switching_function_zeta_;

	// Variables for when the probe volume is set up as a cylindrical shell 
	// with an inner zeta-range and/or a lower limit on r)
	bool has_zeta_range_inner_, has_r_min_;
	double r_min_;
	Range zeta_range_inner_;
	SwitchingFunction_WithRegions_r switching_function_r_min_;
	SwitchingFunction_WithRegions_x switching_function_zeta_inner_;

	// Whether to exclude atoms below the probe volume from the shells
	bool exclude_atoms_below_base_;
};

#endif // PROBE_VOLUME_CYLINDER_H
