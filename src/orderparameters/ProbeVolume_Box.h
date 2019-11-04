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
#include "SwitchingFunction_WithRegions_x.h"

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

	virtual void setGeometry() override;

	void setGeometry(
		const Real3& box_offset,  // Position of lower-left corner: {x, y, z}
		const Box&   box_matrix   // Size and shape of box
	);

	void setGeometry(
		const Real3& x_lower,  // position of lower-left corner
		const Real3& x_upper   // position of upper-right corner
	);

	virtual void calculateIndicator(
		const   Real3& x,
		// Output
		double& h_v,
		double& htilde_v, 
		Real3&  dhtilde_v_dx,
		RegionEnum& region
	) const override;

	// Returns a string with complete, formatted information about the probe volume's
	// location and geometry which can be included in output files. Each line begins
	// with a "#" which signifies that it's a comment line.
	std::string getInputSummary(const std::string& prepend_string) const override;

	void getGeometry(Real3& box_offset, Box& box_matrix) {
		box_offset = box_offset_;
		box_matrix = box_matrix_;
	}

	// Returns whether the probe box is orthrhombic
	bool is_orthorhombic() const {
		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				if ( a != b and box_matrix_[a][b] != 0.0 ) {
					return false;
				}
			}
		}
		return true;
	}

 protected:
	virtual BoundingBox constructBoundingBox() const override;

 private:
	// Nominal ranges for orthorhombic boxes
	std::array<CommonTypes::Range,DIM_> axis_ranges_;
	std::array<CommonTypes::Range,DIM_> inner_axis_ranges_;  // -alpha_c
	std::array<CommonTypes::Range,DIM_> outer_axis_ranges_;  // +alpha_c

	// Box geometry (nominal and effective)
	Real3 box_offset_;  // Position of lower-left corner
	Box   box_matrix_;  // Position of base: {x,y,z}
	Real3 center_;      // Center point: {x,y,z}

	// Half-lengths (assumes an orthorhombic box)
	Real3 box_half_lengths_;

	// Switching functions for a box shifted to the origin
	// - Each runs from -L_i/2 to +L_i/2
	std::array<SwitchingFunction_WithRegions_x, DIM_> switching_functions_shifted_box_;
};

#endif // PROBE_VOLUME_BOX_H
