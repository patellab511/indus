// SimulationBox.cpp
#include "SimulationBox.h"

SimulationBox::SimulationBox()
{
	double l_dummy = 1.0;
	setLengths(l_dummy, l_dummy, l_dummy);
}


SimulationBox::SimulationBox(const Box& box_matrix)
{
	set_box_matrix(box_matrix);
}


SimulationBox::SimulationBox(double l_x, double l_y, double l_z)
{
	setLengths(l_x, l_y, l_z);
}


SimulationBox::~SimulationBox() {}


// Calculate the shift vector that brings 'x' closest to 'ref' under PBCs
// - Assumes orthorhombic box
void SimulationBox::calculateShift(
	const Real3& x, const Real3& ref, Real3& x_shifted, Real3& shift) const
{
	double dx;
	for ( int d=0; d<DIM_; ++d ) {
		// Default: no shift
		x_shifted[d] = x[d];
		shift[d]     = 0.0;

		// Apply minimum image convention
		dx = x[d] - ref[d];
		if ( dx > lengths_h_[d]  ) {
			shift[d]     = -lengths_[d];
			x_shifted[d] += shift[d];
		}
		else if ( dx < -lengths_h_[d] ) {
			shift[d]     =  lengths_[d];
			x_shifted[d] += shift[d];
		}
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


void SimulationBox::set_box_matrix(const Box& box_matrix)
{
	box_matrix_ = box_matrix;

	volume_ = 1.0;
	for ( int d=0; d<DIM_; ++d ) {
		// Diagonal and half-diagonal
		lengths_[d]   = box_matrix_[d][d];
		lengths_h_[d] = 0.5*lengths_[d];

		// Volume of the box (assuming lower triangular matrix)
		volume_ *= lengths_[d];
	}
}


void SimulationBox::setLengths(const double l_x, const double l_y, const double l_z)
{
	// Create full box matrix from diagonal
	Real3 lengths = {{ l_x, l_y, l_z }};
	Box box_matrix;
	for ( int a=0; a<DIM_; ++a ) {
		for ( int b=0; b<DIM_; ++b ) {
			if ( a == b ) {
				box_matrix[a][a] = lengths[a];
			}
			else {
				box_matrix[a][b] = 0.0;
			}
		}
	}

	// Set box matrix and related quantities
	set_box_matrix(box_matrix);
}


void SimulationBox::print(std::ostream& os) const
{
	os << "  box_matrix = {\n";
	for ( int a=0; a<DIM_; ++a ) {
		os << "  ";
		for ( int b=0; b<DIM_; ++b ) {
			if ( b > 0 ) { os << ","; }
			os << " " << box_matrix_[a][b];
		}
		os << "\n";
	}
	os << "  }\n";
}
