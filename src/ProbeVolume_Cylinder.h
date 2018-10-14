/* ProbeVolume_Cylinder.h
 *
 * ABOUT:
 *  - Does not work if the cylinder crosses the edges of the box
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

	virtual void setCoarseGrainingParameters(
			const double sigma, const double alpha_c
	) override;

	virtual void setShellParameters(
			const double width_shell_1, const double width_shell_2, const double alpha_c_shells
	) override;

	void setGeometry(
		const Real3& x_base,  // Position center of base: {x,y,z}
		const double radius,
		const double height
	);

	// Updates the box lengths for ProbeVolume
	virtual void updateUsingSimulationState() override;

	virtual bool isInProbeVolume(
		const Real3& x,
		double& h_v,
		double& htilde_v, 
		Real3& derivatives,
		bool& is_in_shell_1, 
		bool& is_in_shell_2
	) const override;

	// Returns a string with complete, formatted information about the probe volume's
	// location and geometry which can be included in output files. Each line begins
	// with a "#" which signifies that it's a comment line.
	std::string getInputSummary() const override;

 private:
	// Cylinder geometry (nominal and effective)
	// - Effective probe volume includes smoothing over alpha_c
	Real3 x_base_;      // Position of base: {x,y,z}
	Real3 center_;      // Center point: {x,y,z}
	double r_, r_eff_;  // radius
	double h_, h_eff_;  // height
	double zeta_min_, zeta_max_; // Bounds along cylinder axis
	int    axis_index_;         // index of cylinder axis (0=x, 1=y, 2=z)
	std::string axis_label_;

	// Dimensions of shell cylinders
	double r_shell_1_, r_shell_2_;
	double h_shell_1_, h_shell_2_;

	// Precompute squared radii
	double r_eff_sq_, r_shell_1_sq_, r_shell_2_sq_;

	// Half-heights of effective probe volume and shell cylinders
	// - Useful to precompute for checking z-direction
	double h_eff_h_, h_shell_1_h_, h_shell_2_h_;

	// Whether to exclude atoms below the probe volume from the shells
	bool exclude_atoms_below_base_;

	// Switching functions
	SwitchingFunction_r switching_function_r_;
	SwitchingFunction_x switching_function_zeta_;
};

#endif // PROBE_VOLUME_CYLINDER_H
