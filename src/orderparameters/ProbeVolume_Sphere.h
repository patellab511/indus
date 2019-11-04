/* ProbeVolume_Sphere.h
 *
 * ABOUT: 
 * SYNTAX
 *   ProbeVolume = {
 *     type    = sphere
 *     center  = [ <x> <y> <z> ]
 *     r_max   = <r_max>
 *     r_min   = <r_min>      # (default: -1 nm)
 *     sigma   = <sigma>      # (default: 0.01 nm)
 *     alpha_c = <alpha>      # (default: 0.02 nm)
 *   }
 *
 */

#pragma once
#ifndef PROBE_VOLUME_SPHERE_H
#define PROBE_VOLUME_SPHERE_H

// Project headers
#include "ProbeVolume.h"     // Parent class
#include "GenericFactory.h"  // Register as a ProbeVolume

class ProbeVolume_Sphere : public ProbeVolume
{
 public:
	using Real3 = ProbeVolume::Real3;
	using Rvec  = ProbeVolume::Rvec;

	ProbeVolume_Sphere(ProbeVolumeInputPack& input_pack);

	ProbeVolume_Sphere(
		const Real3& center, 
		const double r_min,
		const double r_max,
		ProbeVolumeInputPack& input_pack  // used only to set base class options
	);	

	virtual void setGeometry() override;

	void setGeometry(
		const Real3& center, 
		const double r_min,
		const double r_max
	);

	void set_center(const Real3& center);

	virtual void calculateIndicator(
		const Real3& x, double& h_v, double& htilde_v, Real3& dhtilde_v_dx, RegionEnum& region
	) const override;

	std::string getInputSummary(const std::string& prepend_string) const override;

	// Get the upper radius at which the coarse-grained indicator function becomes zero
	double get_rtilde_max() const {
		double rtilde_min, rtilde_max;
		switching_function_r_shell_.getOuterLimits(rtilde_min, rtilde_max);
		return rtilde_max;
	}

	const Real3& get_center() const {
		return center_;
	}

 protected:
	virtual BoundingBox constructBoundingBox() const override;

 private:
	// Sphere geometry
	Real3  center_;	    // Center of the probe sphere {x,y,z} [nm]
	double r_min_;      // Inner radius of the probe sphere (if a spherical shell) [nm]
	double r_max_;   	  // Outer radius of the probe sphere [nm]

	// Switching function for r (use an 'x' switch to allow spherical shells)
	SwitchingFunction_WithRegions_x switching_function_r_shell_;
};

#endif // PROBE_VOLUME_SPHERE_H	
