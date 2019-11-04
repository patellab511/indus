// Bins.h - Class for managing histogram bins
// - Currently only supports bins of equal size

#ifndef BINS_H
#define BINS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "InputParser.h"

class Bins
{
 public:
	// For bins that span the range [x_i, x_i + delta_x], this determines
	// the nominal x-value of the bin, x_b
	// - Left:   x_b = x_i
	// - Center: x_b = x_i + 0.5*delta_x
	// - Right:  x_b = x_i + delta_x
	enum class BinStyle { Left, Center, Right };

	Bins();

	Bins(
		const double min, 
		const double max, 
		const int num_bins,
		const BinStyle& bin_style = BinStyle::Left
	);

	Bins(const ParameterPack& input_pack);

	int    get_num_bins() const { return num_bins_; };
	double get_bin_size() const { return bin_size_; };
	//std::vector<double> get_bins() const { return bins_; };
	const std::vector<double>& get_bins() const { return bins_; };

	void set_bins(
		const double min, 
		const double max, 
		const int num_bins,
		const BinStyle& bin_style
	);

	void set_bins(const ParameterPack& input_pack);

	// Returns the index of the appropriate bin for x (>=0),
	// else returns -1
	int find_bin(const double x) const;

	//----- Operators ----//

	const double& operator[](const int b) const {
		return bins_[b];
	}

	/* 
	// TODO Is this copy operation possible, for convenience?
	//friend class std::vector<double>;
	std::vector<double> operator=(const Bins& bins) {
		std::vector<double> bins_tmp = bins.bins_;
		return bins_tmp;
	}
	*/

 private:
	double min_, max_, bin_size_;
	int num_bins_;  // store explicitly for convenience
	BinStyle bin_style_;
	std::vector<double> bins_;

	BinStyle parseBinStyle(const std::string& bin_style_token) const;
};

#endif /* BINS_H */
