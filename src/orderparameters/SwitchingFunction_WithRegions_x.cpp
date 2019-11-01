#include "SwitchingFunction_WithRegions_x.h"

SwitchingFunction_WithRegions_x::SwitchingFunction_WithRegions_x(
	const double x_min, const double x_max, const double sigma, const double alpha_c,
	const double width_shell_1, const double width_shell_2, const bool is_inverted
):
	SwitchingFunction_WithRegions(),
	switch_x_(x_min, x_max, sigma, alpha_c),
	is_inverted_(is_inverted)
{
	setBoundaries(x_min, x_max, width_shell_1, width_shell_2);
}


// 'Core' set function that sets all relevant member variables
void SwitchingFunction_WithRegions_x::setBoundaries(
	const double x_min, const double x_max, const double width_shell_1, const double width_shell_2
)
{
	if ( width_shell_1 < 0.0 or width_shell_2 < 0.0 ) {
		std::stringstream err_ss; 
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
					 << "  Shell widths may not be negative. Input:\n"
		       << "    width_shell_1 = " << width_shell_1 << "\n"
		       << "    width_shell_2 = " << width_shell_2 << "\n";
		throw std::runtime_error( err_ss.str() );
	}

	switch_x_.setLimits(x_min, x_max);

	// Shell sizes
	width_shell_1_ = width_shell_1;
	width_shell_2_ = width_shell_2;

	if ( width_shell_1_ > 0.0 or width_shell_2_ > 0.0 ) {
		nontrivial_shells_ = true;
	}
	else {
		nontrivial_shells_ = false;
	}

	// Outer limits of shells
	if ( not is_inverted_ ) {
		double x_lowest, x_highest;
		switch_x_.getOuterLimits(x_lowest, x_highest);

		x_min_shell_1_ = x_lowest  - width_shell_1;
		x_max_shell_1_ = x_highest + width_shell_1;
		x_min_shell_2_ = x_min_shell_1_ - width_shell_2;
		x_max_shell_2_ = x_max_shell_1_ + width_shell_2;
	}
	else {
		double x_min_inner, x_max_inner;
		switch_x_.getInnerLimits(x_min_inner, x_max_inner);

		x_min_shell_1_ = x_min_inner + width_shell_1;
		x_max_shell_1_ = x_max_inner - width_shell_1;
		x_min_shell_2_ = x_min_shell_1_ + width_shell_2;
		x_max_shell_2_ = x_max_shell_1_ - width_shell_2;
	}
}


void SwitchingFunction_WithRegions_x::calculate(
	const double x, const bool want_derivative,
	double& h_x, double& htilde_x, double& dhtilde_x_dx, RegionEnum& region
) const
{
	region = RegionEnum::Unimportant;  // default region

	if ( not is_inverted_ ) {
		if ( nontrivial_shells_ and (x < x_min_shell_2_ or x > x_max_shell_2_) ) {
			// Completely outside region of interest
			// - Consider this case first for performance, assuming that most 
			//   values will be in irrelevant regions
			h_x      = 0.0;
			htilde_x = 0.0;
			if ( want_derivative ) { dhtilde_x_dx = 0.0; }
			return;
		}
		else {
			switch_x_.calculate( x, want_derivative, 
													 h_x, htilde_x, dhtilde_x_dx );

			// Set region
			if ( htilde_x > 0.0 ) {
				region = RegionEnum::Vtilde;
			}
			else if ( nontrivial_shells_ ) {
				if ( x >= x_min_shell_1_ and x <= x_max_shell_1_ ) {
					region = RegionEnum::Shell_1;
				}
				else if ( x >= x_min_shell_2_ and x <= x_max_shell_2_ ) {
					region = RegionEnum::Shell_2;
				}
			}
		}
	}
	else {
		// Inverted switching function
		switch_x_.calculate( x, want_derivative, 
												 h_x, htilde_x, dhtilde_x_dx );
		h_x          = 1.0 - h_x;
		htilde_x     = 1.0 - htilde_x;
		dhtilde_x_dx = -dhtilde_x_dx;

		// Set region
		if ( htilde_x > 0.0 ) {
			region = RegionEnum::Vtilde;
		}
		else if ( nontrivial_shells_ ) {
			h_x      = 0.0;
			htilde_x = 0.0;
			if ( want_derivative ) { dhtilde_x_dx = 0.0; }

			// Note that the following will produce the correct result, even when
			// width_shell_1/2 is large enough to cause overlap between the min/max regions
			if ( x <= x_min_shell_1_ or x >= x_max_shell_1_ ) {
				region = RegionEnum::Shell_1;
			}
			else if ( x <= x_min_shell_2_ or x >= x_max_shell_2_ ) {
				region = RegionEnum::Shell_2;
			}
		}
	}
}
