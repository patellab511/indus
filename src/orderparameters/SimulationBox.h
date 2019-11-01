// SimulationBox.h

#pragma once
#ifndef SIMULATION_BOX_H
#define SIMULATION_BOX_H

#include "CommonTypes.h"

#include <array>
#include <vector>

class SimulationBox
{
 public:
	static const int DIM_ = CommonTypes::DIM_;
	using Real   = CommonTypes::Real;
	using Real3  = CommonTypes::Real3;
	using Box    = CommonTypes::Box;

	SimulationBox();
	SimulationBox(const Box& box_matrix);
	SimulationBox(double l_x, double l_y, double l_z);
	~SimulationBox();


	//----- Minimum Image Convention/PBCs -----//

	// Calculate the shift vector that brings 'x' closest to 'ref' under PBCs
	void calculateShift(
		const Real3& x, 
		const Real3& ref,
		// Output 
		Real3& x_shifted,  // image of 'x' closest to 'ref'
		Real3& shift       // Vector used to shift x = x_shifted - x
	) const;

	// Calculate the minimum-image vector between 2 points and its squared magnitude
	//  Input:
	//		x1, x2 = positions {x,y,z}
	//  Output:
	//  	x12    = x_{1,2} = x_2 - x_1
	//  	dist_sq = norm(x12)^2
	void calculateDistance(const Real3& x1, const Real3& x2, Real3& x12, double& dist_sq) const;

	// C-array input
	template <typename T>
	void calculateDistance(const T x1[3], const T x2[3], T x12[3], double& dist_sq) const;
	template <typename T>
	void calculateDistance(const T x1[3], const T x2[3], Real3& x12, double& dist_sq) const;

	// Mixed input types
	template <typename T>
	void calculateDistance(const T x1[3], const Real3& x2, T x12[3], double& dist_sq) const;
	template <typename T>
	void calculateDistance(const Real3& x1, const T x2[3], T x12[3], double& dist_sq) const;
	template <typename T>
	void calculateDistance(const Real3& x1, const T x2[3], Real3& x12, double& dist_sq) const;

	// Template for generic rvec-like objects (3 elements) with operator[] access
	template <typename T, typename U, typename V>
	void calculateDistance(const T& x1, const U& x2, V& x12, double& dist_sq) const;

	// Modifies the given position so that it is within the simulation box
	// - For positions that are right outside the simulation box
	// - Only works if position is within one periodic image of the real
	//   simulation box
	// - Assumes an orthorhombic box
	void putInBox(Real3& x) const {
		for ( int d=0; d<DIM_; ++d ) {
			if ( x[d] > lengths_[d] ) {
				x[d] -= lengths_[d];
				if ( x[d] < 0.0 ) {
					x[d] = lengths_[d];  // atom is at the edge of the box to within precision
				}
			}
			else if ( x[d] < 0.0 ) {
				x[d] += lengths_[d];
				if ( x[d] > lengths_[d] ) {
					x[d] = 0.0;  // atom is at the edge of the box to within precision
				}
			}
		}
	}

	Real3 putInBox(const Real3& x) const {
		Real3 x_out = x;
		putInBox(x_out);
		return x_out;
	}

	bool is_in_box(const Real3& x) const {
		// Assume box is orthorhombic
		for ( int d=0; d<DIM_; ++d ) {
			if ( x[d] < 0.0 or x[d] > lengths_[d] ) {
				return false;
			}
		}
		return true;
	}


	//----- Managing box matrix -----//

	// Sets the full box matrix
	void set_box_matrix(const Box& box_matrix);

	// Set the box lengths and half-lengths (orthorhombic box)
	void setLengths(const double l_x, const double l_y, const double l_z);

	// Return the box length(s)
	Real  getLength(const int dim) const { return lengths_[dim]; }
	Real3 getLengths()             const { return lengths_; }

	// Returns the volume of the box
	double getVolume() const { return volume_; }

	const Box& get_box_matrix() const { return box_matrix_; }


	//----- Box Type -----//

	// Returns "true" if the box is orthorhombic; returns false otherwise
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

	bool is_pbc_3d() const {
		return is_pbc_3d_;
	}


	//----- Debugging -----//

	void print(std::ostream& os) const;
	friend std::ostream& operator<<(std::ostream& os, const SimulationBox& simulation_box) {
		simulation_box.print(os);
		return os;
	}

 protected:
	// Whether 3D PBCs are in effect
	bool is_pbc_3d_ = true;

	// Simulation box
	// - Lower triangular matrix (as in GROMACS)
	// - Origin is at the lower-left corner of the box
	Box box_matrix_;
	double volume_;

	// For orthorhombic boxes
	Real3  lengths_;	  // Box side lengths
	Real3  lengths_h_;	// Box side half-lengths
};



//----- Templates -----//

// Calculate the minimum-image vector between 2 points and its squared magnitude
template <typename T>
void SimulationBox::calculateDistance(const T x1[3], const T x2[3], T x12[3], double& dist_sq) const
{
	dist_sq = 0.0;
	for ( int d=0; d<3; ++d ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] > lengths_h_[d]  ) { x12[d] -= lengths_[d]; }
		else if ( x12[d] < -lengths_h_[d] ) { x12[d] += lengths_[d]; }

		dist_sq += x12[d]*x12[d];
	}
	return;
} 

// Calculate the minimum-image vector between 2 points and its squared magnitude
template <typename T>
void SimulationBox::calculateDistance(const T x1[3], const T x2[3], SimulationBox::Real3& x12, 
                             double& dist_sq) const
{
	dist_sq = 0.0;
	for ( int d=0; d<3; ++d ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] > lengths_h_[d]  ) { x12[d] -= lengths_[d]; }
		else if ( x12[d] < -lengths_h_[d] ) { x12[d] += lengths_[d]; }

		dist_sq += x12[d]*x12[d];
	}
	return;
}

template <typename T>
void SimulationBox::calculateDistance(const T x1[3], const SimulationBox::Real3& x2, 
                             T x12[3], double& dist_sq) const
{
	dist_sq = 0.0;
	for ( int d=0; d<3; ++d ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] > lengths_h_[d]  ) { x12[d] -= lengths_[d]; }
		else if ( x12[d] < -lengths_h_[d] ) { x12[d] += lengths_[d]; }

		dist_sq += x12[d]*x12[d];
	}
	return;
} 

template <typename T>
void SimulationBox::calculateDistance(const SimulationBox::Real3& x1, const T x2[3],
                             T x12[3], double& dist_sq) const
{
	dist_sq = 0.0;
	for ( int d=0; d<3; ++d ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] > lengths_h_[d]  ) { x12[d] -= lengths_[d]; }
		else if ( x12[d] < -lengths_h_[d] ) { x12[d] += lengths_[d]; }

		dist_sq += x12[d]*x12[d];
	}
	return;
} 

template <typename T>
void SimulationBox::calculateDistance(const SimulationBox::Real3& x1, const T x2[3], 
                             SimulationBox::Real3& x12, double& dist_sq) const
{
	dist_sq = 0.0;
	for ( int d=0; d<3; ++d ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] > lengths_h_[d]  ) { x12[d] -= lengths_[d]; }
		else if ( x12[d] < -lengths_h_[d] ) { x12[d] += lengths_[d]; }

		dist_sq += x12[d]*x12[d];
	}
	return;
} 

// Template for generic rvec-like objects (3 elements) with operator[] access
template <typename T, typename U, typename V>
void SimulationBox::calculateDistance(const T& x1, const U& x2, V& x12, double& dist_sq) const
{
	dist_sq = 0.0;
	for ( int d=0; d<3; ++d ) {
		x12[d] = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( x12[d] > lengths_h_[d]  ) { x12[d] -= lengths_[d]; }
		else if ( x12[d] < -lengths_h_[d] ) { x12[d] += lengths_[d]; }

		dist_sq += x12[d]*x12[d];
	}
	return;
}

#endif // SIMULATION_BOX_H
