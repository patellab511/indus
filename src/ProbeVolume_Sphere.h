/* ProbeVolume_Sphere.h
 *
 * ABOUT: 
 * SYNTAX
 *   ProbeVolume = {
 *     type    = sphere
 *     center  = [ <x> <y> <z> ]
 *     radius  = <r>
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

	virtual void setCoarseGrainingParameters(
			const double sigma, const double alpha_c
	) override;

	virtual void setShellWidths(
			const double width_shell_1, const double width_shell_2
	) override;

	void setGeometry(
		const Real3& center, 
		const double r
	);

	// Updates the ProbeVolume for each simulation step
	virtual void updateUsingSimulationState() override;

	virtual bool isInProbeVolume(
			const Real3& x, double& h_v, double& htilde_v, Real3& dhtilde_v_dx, 
			bool& is_in_shell_1, bool& is_in_shell_2
	) const override;

	// Returns a string with complete, formatted information about the probe volume's
	// location and geometry which can be included in output files.
	virtual std::string getInputSummary() const override;

 private:
	// Sphere geometry
	Real3 center_;	  // Center of the probe sphere {x,y,z} [nm]
	double r_;   	    // Radius of the probe sphere [nm]
	double r_sq_;

	// Switching function for r
	SwitchingFunction_r switching_function_r_;

	// Effective radius of the probe volume, including the buffer [nm]
	// - The switching function is zero at: r_eff = r + alpha_c
	double r_eff_, r_eff_sq_;

	// Define shells around the probe volume
	// - This is useful for efficiently locating atoms NEAR the probe volume
	double r_shell_1_, r_shell_1_sq_;
	double r_shell_2_, r_shell_2_sq_;
};

#endif // PROBE_VOLUME_SPHERE_H	
