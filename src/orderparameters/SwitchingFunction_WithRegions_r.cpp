#include "SwitchingFunction_WithRegions_r.h"

SwitchingFunction_WithRegions_r::SwitchingFunction_WithRegions_r(
	const double r_max, const double sigma, const double alpha_c,
	const double width_shell_1, const double width_shell_2, const bool is_inverted
):
	SwitchingFunction_WithRegions(),
	switch_r_(r_max, sigma, alpha_c),
	is_inverted_(is_inverted)
{
	setBoundaries(r_max, width_shell_1, width_shell_2);
}


// 'Core' set function that sets all relevant member variables
void SwitchingFunction_WithRegions_r::setBoundaries(
	const double r_max, const double width_shell_1, const double width_shell_2
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

	switch_r_.set_r_max(r_max);

	// Shell sizes
	width_shell_1_ = width_shell_1;
	width_shell_2_ = width_shell_2;

	if ( width_shell_1_ > 0.0 or width_shell_2_ > 0.0 ) {
		nontrivial_shells_ = true;
	}
	else {
		nontrivial_shells_ = false;
	}

	if ( not is_inverted_ ) {
		// Where the shells begin
		// - This is r_max plus the primary region's coarse-graining over alpha_c
		double r_highest = switch_r_.get_r_outer();

		// Outer limits of shells
		r_max_shell_1_ = r_highest + width_shell_1;
		r_max_shell_2_ = r_max_shell_1_ + width_shell_2;
	}
	else {
		double r_lowest = switch_r_.get_r_inner();
		r_max_shell_1_ = r_lowest - width_shell_1;
		r_max_shell_2_ = r_max_shell_1_ - width_shell_2;
	}
}


void SwitchingFunction_WithRegions_r::calculate(
	const double r, const bool want_derivative,
	double& h_r, double& htilde_r, double& dhtilde_r_dr, RegionEnum& region
) const
{
	region = RegionEnum::Unimportant;  // default region

	if ( not is_inverted_ ) {
		if ( nontrivial_shells_ and r > r_max_shell_2_ ) {
			// Completely outside region of interest
			// - Consider this case first for performance, assuming that most 
			//   values will be in irrelevant regions
			h_r      = 0.0;
			htilde_r = 0.0;
			if ( want_derivative ) { dhtilde_r_dr = 0.0; }
			return;
		}
		else {
			switch_r_.calculate( r, want_derivative, 
													 h_r, htilde_r, dhtilde_r_dr );

			// Set region
			if ( htilde_r > 0.0 ) {
				region = RegionEnum::Vtilde;
			}
			else if ( nontrivial_shells_ ) {
				if ( r <= r_max_shell_1_ ) {
					region = RegionEnum::Shell_1;
				}
				else if ( r <= r_max_shell_2_ ) {
					region = RegionEnum::Shell_2;
				}
			}
		}
	}
	else {
		// Inverted switching function: NOT[htilde_r(r)] = 1 - htilde_r(r)
		if ( r < r_max_shell_2_ ) {
			h_r      = 0.0;
			htilde_r = 0.0;
			if ( want_derivative ) { dhtilde_r_dr = 0.0; }
			return;
		}
		else {
			switch_r_.calculate( r, want_derivative, 
													 h_r, htilde_r, dhtilde_r_dr );
			h_r          = 1.0 - h_r;
			htilde_r     = 1.0 - htilde_r;
			dhtilde_r_dr = -dhtilde_r_dr;

			// Set region
			if ( htilde_r > 0.0 ) {
				region = RegionEnum::Vtilde;
			}
			else if ( nontrivial_shells_ ) {
				if ( r >= r_max_shell_1_ ) {
					region = RegionEnum::Shell_1;
				}
				else if ( r >= r_max_shell_2_ ) {
					region = RegionEnum::Shell_2;
				}
			}
		}
	}
}

