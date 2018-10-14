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
	SimulationBox();
	SimulationBox(double l_x, double l_y, double l_z);
	~SimulationBox();

	static const int DIM_ = CommonTypes::DIM_;
	using Real   = CommonTypes::Real;
	using Real3  = CommonTypes::Real3;
	using Box    = CommonTypes::Box;

	//----- Minimum Image Convention/PBCs -----//

	// Calculate the minimum-image vector between 2 points and its squared magnitude
	//  Input:
	//		x1, x2 = positions {x,y,z}
	//  Output:
	//  	x12    = x_{1,2} = x_2 - x_1
	//  	dist_sq = norm(x12)^2

	void calculateDistance(const std::vector<double>& x1, const std::vector<double>& x2,
	              std::vector<double>&x12, double& dist_sq) const;

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
			}
			else if ( x[d] < 0.0 ) {
				x[d] += lengths_[d];
			}
		}
	}

	Real3 putInBox(const Real3& x) const {
		Real3 x_out = x;
		putInBox(x_out);
		return x_out;
	}

	//----- Managing box lengths -----//

	// Set the box lengths and half-lengths (orthorhombic box)
	void setLengths(const double l_x, const double l_y, const double l_z);

	// Return the box length(s)
	Real  getLength(const int dim) const { return lengths_[dim]; }
	Real3 getLengths() const             { return lengths_; }

	// Returns the volume of the box
	double getVolume() const { return volume_; }

	// Returns the box matix (asumes orthorhombic box)
	void getBoxMatrix(Box& box_matrix) const {
		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				if ( a == b ) {
					box_matrix[a][a] = lengths_[a];
				}
				else {
					box_matrix[a][b] = 0.0;
				}
			}
		}
	}

 protected:
	// Simulation box
	double volume_;
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
