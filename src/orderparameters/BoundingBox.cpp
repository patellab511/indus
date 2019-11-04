#include "BoundingBox.h"

void BoundingBox::print(std::ostream& os) const
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

	os << "  x_lower = {";
	for ( int d=0; d<DIM_; ++d ) {
		if ( d > 0 ) { os << ","; }
		os << " " << x_lower_[d];
	}
	os << " }\n";

	os << "  x_upper = {";
	for ( int d=0; d<DIM_; ++d ) {
		if ( d > 0 ) { os << ","; }
		os << " " << x_upper_[d];
	}
	os << " }\n";
}

std::ostream& operator<<(std::ostream& os, const BoundingBox& bounding_box)
{
	bounding_box.print(os);
	return os;
}
