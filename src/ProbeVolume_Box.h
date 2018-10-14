/* ProbeVolume_Box.h
 *
 * ABOUT:
 *  - Currently only supports orthorhombic boxes
 *
 * DEVELOPMENT:
 */

#pragma once
#ifndef PROBE_VOLUME_BOX_H
#define PROBE_VOLUME_BOX_H

// Standard headers
#include <map>

// Core headers
#include "ProbeVolume.h"     // Parent class
#include "GenericFactory.h"  // Register as a ProbeVolume


class ProbeVolume_Box : public ProbeVolume
{
 public:
	using Box   = ProbeVolume::Box;
	using Real3 = ProbeVolume::Real3;
	using Rvec  = ProbeVolume::Rvec;

	ProbeVolume_Box(ProbeVolumeInputPack& input_pack);

	// Constructor which ignores the ParameterPack in the input pack
	// and instead provides the necessary variables explicitly
	ProbeVolume_Box(
		ProbeVolumeInputPack& input_pack,
		const Real3& box_offset,
		const Box&   box_matrix,
		const double sigma,
		const double alpha_c
	);

	virtual void setCoarseGrainingParameters(
			const double sigma, const double alpha_c
	) override;

	virtual void setShellParameters(
			const double width_shell_1, const double width_shell_2, const double alpha_c_shells
	) override;

	void setGeometry(
		const Real3& box_offset,  // Position of lower-left corner: {x, y, z}
		const Box&   box_matrix   // Size and shape of box
	);

	// Updates the box lengths for ProbeVolume
	virtual void updateUsingSimulationState() override;

	virtual bool isInProbeVolume(
		const   Real3& x,
		// Output
		double& h_v,
		double& htilde_v, 
		Real3&  dhtilde_v_dx,
		bool&   is_in_shell_1, 
		bool&   is_in_shell_2
	) const override;

	// Returns a string with complete, formatted information about the probe volume's
	// location and geometry which can be included in output files. Each line begins
	// with a "#" which signifies that it's a comment line.
	std::string getInputSummary() const override;

	void getGeometry(Real3& box_offset, Box& box_matrix) {
		box_offset = box_offset_;
		box_matrix = box_matrix_;
	}

 private:
	// Box geometry (nominal and effective)
	Real3 box_offset_;  // Position of lower-left corner
	Box   box_matrix_;  // Position of base: {x,y,z}
	Real3 center_;      // Center point: {x,y,z}

	// Half-lengths (assumes an orthorhombic box)
	// - e.g. 
	//     box_half_lengths_     = box_lengths/2.0
	//     box_half_lengths_eff_ = box_lengths/2.0 + alpha_c
	Real3 box_half_lengths_;             // no coarse-graining
	Real3 box_half_lengths_eff_;         // with coarse-graining (alpha_c)
	Real3 box_shell_1_half_lengths_eff_; // with coarse-graining (alpha_c_shells)
	Real3 box_shell_2_half_lengths_eff_; // with coarse-graining (alpha_c_shells)

	// Switching functions for a box shifted to the origin
	// - Each runs from -L_i/2 to +L_i/2
	std::array<SwitchingFunction_x, DIM_> switching_functions_shifted_box_;
};

#endif // PROBE_VOLUME_BOX_H
