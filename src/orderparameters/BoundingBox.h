// The "bounding box" incorporates any extra distance around the nominal region due to coarse-graining
// - TODO: should it include the shells in its definition?

// TODO: Move "bounding box" treatment from other parts of the code to here
// to cut down on code duplication

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "CommonTypes.h"
#include "SimulationBox.h"

class BoundingBox
{
 public:
	static constexpr int DIM_ = CommonTypes::DIM_;
	using Real3 = CommonTypes::Real3;
	using Box   = CommonTypes::Box;

	BoundingBox(const Real3& x_lower, const Real3& x_upper, const SimulationBox& simulation_box):
		simulation_box_ptr_(&simulation_box)
	{
		setUsingCorners(x_lower, x_upper);
	} 

	BoundingBox(const SimulationBox& simulation_box):
		simulation_box_ptr_(&simulation_box)
	{
		setUsingBoxMatrix( simulation_box.get_box_matrix() );
	} 

	// Assumes the box is orthorhombic
	void setUsingCorners(const Real3& x_lower, const Real3& x_upper)
	{
		x_lower_    = x_lower;
		x_upper_    = x_upper;
		box_matrix_ = makeBoxMatrixUsingCorners(x_lower_, x_upper_);
		keepBoxWithinPbcs();
	}

	void setUsingBoxMatrix(const Box& box_matrix) {
		box_matrix_ = box_matrix;
		for ( int d=0; d<DIM_; ++d ) {
			x_lower_[d] = 0.0;
			x_upper_[d] = box_matrix_[d][d];
		}
		keepBoxWithinPbcs();
	}

	// Access variables
	const Real3& get_x_lower()    const { return x_lower_;    }
	const Real3& get_x_upper()    const { return x_upper_;    }
	const Box&   get_box_matrix() const { return box_matrix_; }

	// Extends the box in all directions by 'dx'
	// - Assumes an orthorhombic box
	void extend(const double dx) {
		for ( int d=0; d<DIM_; ++d ) {
			x_lower_[d] -= dx;
			x_upper_[d] += dx;
		}
		setUsingCorners(x_lower_, x_upper_);
	}

	// Printing
	void print(std::ostream& os) const;
	friend std::ostream& operator<<(std::ostream& os, const BoundingBox& bounding_box);

 private:
	const SimulationBox* simulation_box_ptr_ = nullptr;

	// Upper and lower corners of the box
	Real3 x_lower_, x_upper_;  

	Box box_matrix_;

	// True if the box spans PBCs along the given axis
	std::array<bool,DIM_> spans_pbcs_;

	Box makeBoxMatrixUsingCorners(const Real3& x_lower, const Real3& x_upper) const {
		Box box_matrix;
		for ( int a=0; a<DIM_; ++a ) {
			for ( int b=0; b<DIM_; ++b ) {
				if ( a == b ) {
					box_matrix[a][a] = x_upper[a] - x_lower[a];
				}
				else {
					box_matrix[a][b] = 0.0;
				}
			}
		}
		return box_matrix;
	}

	// If any part of the box crosses PBCs in a particular direction, extend the box
	// so that it spans PBCs along that axis
	void keepBoxWithinPbcs() {
		const Box& simulation_box_matrix = simulation_box_ptr_->get_box_matrix();
		spans_pbcs_.fill(false);
		for ( int d=0; d<DIM_; ++d ) {
			if ( (x_lower_[d] < 0.0) or (x_upper_[d] > simulation_box_matrix[d][d]) ) {
				x_lower_[d]       = 0.0;
				x_upper_[d]       = simulation_box_matrix[d][d];
				box_matrix_[d][d] = simulation_box_matrix[d][d];
				spans_pbcs_[d]    = true;
			}
		}
	}
};

#endif // ifndef BOUNDING_BOX_H
