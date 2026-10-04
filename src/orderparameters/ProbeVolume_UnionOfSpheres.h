/* ProbeVolume_UnionOfSpheres.h
 *
 * ABOUT: Union of K static spheres (or spherical shells) with common radii.
 *
 *   The coarse-grained indicator is the smooth "OR" over the subvolumes,
 *
 *     htilde_v(x) = 1 - prod_k [ 1 - htilde_k(x) ]
 *
 *   so an atom in several overlapping spheres is counted once. Each subvolume
 *   is an ordinary ProbeVolume_Sphere centered at one of the positions listed
 *   in 'reference_positions_file'. The derivative follows from the product rule:
 *
 *     d htilde_v / dx = sum_k [ prod_{l != k} (1 - htilde_l) ] d htilde_k / dx
 *
 * SYNTAX
 *   ProbeVolume = {
 *     type    = union_of_spheres            # alias: union_of_spherical_shells
 *     r_max   = <r_max>
 *     r_min   = <r_min>                     # (default: -1 nm, i.e. solid spheres)
 *     reference_positions_file = <file>     # one "x y z" (nm) per line; '#' starts a comment
 *     sigma   = <sigma>                     # (default: 0.01 nm)
 *     alpha_c = <alpha>                     # (default: 0.02 nm)
 *   }
 *
 * NOTES
 *   - All spheres share r_min, r_max, sigma, alpha_c and the shell widths.
 *   - Every subvolume is checked for every atom (brute force). This is fine for
 *     up to a few hundred spheres; a cell list over the centers can be added
 *     behind the same interface if needed.
 */

#pragma once
#ifndef PROBE_VOLUME_UNION_OF_SPHERES_H
#define PROBE_VOLUME_UNION_OF_SPHERES_H

// Standard headers
#include <memory>
#include <string>
#include <vector>

// Project headers
#include "GenericFactory.h"     // Register as a ProbeVolume
#include "OpenMP.h"
#include "ProbeVolume.h"        // Parent class
#include "ProbeVolume_Sphere.h" // Subvolume geometry

class ProbeVolume_UnionOfSpheres : public ProbeVolume
{
 public:
	using Real3 = ProbeVolume::Real3;

	ProbeVolume_UnionOfSpheres(ProbeVolumeInputPack& input_pack);

	// Pushes coarse-graining parameters, shell widths, radii and centers into the subvolumes
	virtual void setGeometry() override;

	// *** Must be thread-safe: called from inside the OpenMP loop in Indus::calculate() ***
	virtual void calculateIndicator(
		const Real3& x, double& h_v, double& htilde_v, Real3& dhtilde_v_dx, RegionEnum& region
	) const override;

	std::string getInputSummary(const std::string& prepend_string) const override;

	int get_num_subvolumes() const {
		return static_cast<int>( subvolume_ptrs_.size() );
	}

	const std::vector<Real3>& get_sphere_centers() const {
		return sphere_centers_;
	}

	// Largest distance from a sphere center at which htilde_k can be nonzero
	// (r_max + alpha_c; does not include shells)
	double get_rtilde_max() const {
		return rtilde_max_;
	}

 protected:
	// Smallest axis-aligned box around all centers, padded by rtilde_max
	virtual BoundingBox constructBoundingBox() const override;

 private:
	// Reads "x y z" triplets (nm) from a plain-text file; '#' begins a comment
	static std::vector<Real3> readSphereCenters(const std::string& file);

	// Geometry shared by all subvolumes
	double r_min_;  // Inner radius [nm] (negative: solid spheres)
	double r_max_;  // Outer radius [nm]
	double rtilde_max_ = 0.0;

	std::string reference_positions_file_;

	std::vector<Real3> sphere_centers_;

	using SubvolumePtr = std::unique_ptr<ProbeVolume_Sphere>;
	std::vector<SubvolumePtr> subvolume_ptrs_;

	// Per-thread scratch space for the subvolumes that are "active" for the current atom
	// (0 < htilde_k < 1), needed for the derivative loop
	struct SubvolumesBuffer {
		std::vector<double> htilde;      // htilde_k
		std::vector<Real3>  dhtilde_dx;  // d htilde_k / dx

		void clear() {
			htilde.resize(0);
			dhtilde_dx.resize(0);
		}
	};
	mutable std::vector<SubvolumesBuffer> subvolumes_buffers_;

	// Returns the "stronger" of two regions (Vtilde > Shell_1 > Shell_2 > Unimportant)
	static RegionEnum strongerRegion(const RegionEnum a, const RegionEnum b) {
		return ( static_cast<int>(a) <= static_cast<int>(b) ) ? a : b;
	}
};

#endif // PROBE_VOLUME_UNION_OF_SPHERES_H
