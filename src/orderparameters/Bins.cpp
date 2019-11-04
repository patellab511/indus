#include "Bins.h"

Bins::Bins(): 
	min_(0.0), max_(0.0), num_bins_(0), 
	bin_style_(BinStyle::Left) 
{
	set_bins(min_, max_, num_bins_, bin_style_);
}


Bins::Bins(
	const double min, const double max, const int num_bins, const BinStyle& bin_style)
{
	set_bins(min, max, num_bins, bin_style);
}


Bins::Bins(const ParameterPack& input_pack)
 : bin_style_(BinStyle::Left)
{
	set_bins(input_pack);
}



void Bins::set_bins(
	const double min, const double max, const int num_bins, const BinStyle& bin_style)
{
	// Check input
	if ( num_bins < 0 ) {
		throw std::runtime_error("error in Bins::set_bins: the number of bins must be non-negative");
	}
	else if ( min > max ) {
		throw std::runtime_error("error in Bins::set_bins: min must not be greater than max");
	}

	min_ = min;
	max_ = max;
	bin_style_ = bin_style;
	num_bins_ = num_bins;

	bins_.resize(num_bins_);
	if ( num_bins > 0 ) {
		bin_size_ = (max_ - min_)/static_cast<double>(num_bins_);
	}
	else {
		// No bins: nothing more to do
		bin_size_ = 0.0;
		return;
	}

	// Set the bin alignment for the first bin
	bins_[0] = min_;  // default: align Left
	if ( bin_style_ == BinStyle::Center ) {
		bins_[0] += 0.5*bin_size_;
	}
	else if ( bin_style == BinStyle::Right ) {
		bins_[0] += bin_size_;
	}

	// Subsequent bins are evenly spaced
	for ( int b=1; b<num_bins_; ++b ) {
		bins_[b] = bins_[b-1] + bin_size_;
	}
}

void Bins::set_bins(const ParameterPack& input_pack)
{
	using KeyType = ParameterPack::KeyType;

	// Working variables
	std::string token;
	bool found;
	int num_bins;
	BinStyle bin_style = BinStyle::Left;  // default style

	input_pack.readNumber("NumBins", KeyType::Required, num_bins);
	std::array<double,2> x_range;
	input_pack.readArray("BinRange", KeyType::Required, x_range);

	found = input_pack.readString("BinStyle", KeyType::Optional, token);
	if ( found ) {
		bin_style = parseBinStyle(token);
	}

	set_bins(x_range[0], x_range[1], num_bins, bin_style);
}


int Bins::find_bin(const double x) const 
{
	int bin = static_cast<int>(floor( (x - min_)/bin_size_ ));
	if ( (bin >= 0) and (bin <= num_bins_) and (x <= max_) ) {
		if ( bin == num_bins_ ) {
			// x = max_ belongs to the rightmost bin
			--bin;
		}
		return bin;
	}
	else {
		return -1;
	}
}


Bins::BinStyle Bins::parseBinStyle(const std::string& bin_style_token) const
{
	// Case-insensitive parsing
	std::string lowercase_token = bin_style_token;
	std::transform( lowercase_token.begin(), lowercase_token.end(), 
	                lowercase_token.begin(), ::tolower );

	if ( lowercase_token == "left" ) {
		return BinStyle::Left;
	}
	else if ( lowercase_token == "center" ) {
		return BinStyle::Center;
	}
	else if ( lowercase_token == "right" ) {
		return BinStyle::Right;
	}
	else {
		throw std::runtime_error("bin style \"" + bin_style_token + "\" was not recognized");
	}
}
