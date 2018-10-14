// SimulationBox.cpp

#include <array>
#include <vector>

#include "SimulationBox.h"

SimulationBox::SimulationBox()
{
	setLengths(0.0, 0.0, 0.0);
}

SimulationBox::SimulationBox(double l_x, double l_y, double l_z)
{
	setLengths(l_x, l_y, l_z);
}


SimulationBox::~SimulationBox() {}


// Calculate the minimum-image vector between 2 points and its squared magnitude
void SimulationBox::calculateDistance(const std::vector<double>& x1, const std::vector<double>& x2,
                             std::vector<double>& x12, double& distSq) const
{
	distSq = 0.0;
	x12.resize(3); // r_1,2 = r_2 - r_1

	double dx;
	for ( int d=0; d<3; ++d ) {
		dx = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( dx > lengths_h_[d]  ) { dx -= lengths_[d]; }
		else if ( dx < -lengths_h_[d] ) { dx += lengths_[d]; }

		distSq += dx*dx;
		x12[d] = dx;
	}
	return;
}


void SimulationBox::calculateDistance(const Real3& x1, const Real3& x2, 
                             Real3& x12, double& distSq) const
{
	distSq = 0.0;

	double dx;
	for ( int d=0; d<3; ++d ) {
		dx = x2[d] - x1[d];
		// Apply minimum image convention
		if      ( dx > lengths_h_[d]  ) { dx -= lengths_[d]; }
		else if ( dx < -lengths_h_[d] ) { dx += lengths_[d]; }

		distSq += dx*dx;
		x12[d] = dx;
	}
	return;
}


// Set the box lengths and half-lengths
void SimulationBox::setLengths(const double l_x, const double l_y, const double l_z)
{
	lengths_[0] = l_x;    lengths_h_[0] = 0.5*l_x;
	lengths_[1] = l_y;    lengths_h_[1] = 0.5*l_y;
	lengths_[2] = l_z;    lengths_h_[2] = 0.5*l_z;

	volume_ = l_x*l_y*l_z;
}

