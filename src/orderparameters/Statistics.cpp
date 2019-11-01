#include "Statistics.h"

// Constructor/destructor
Statistics::Statistics(const int seed)
{
	// Generate series of unsigned ints for RNG (below)
	std::seed_seq seedSequence = { seed };

	// Random number generator: Mersenne Twister
	rng_ptr_ = Statistics::MersenneTwisterPtr( new std::mt19937(seedSequence) );
}



Statistics::~Statistics()
{
	// Nothing to see here;
}



// Compute the sample mean
double Statistics::average(const std::vector<double>& x)
{
	int num_samples = x.size();
	if ( num_samples < 1 ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  Insufficient number of samples\n";
		throw std::runtime_error( err_ss.str() );
	}

	// TODO: This implementation is subject to roundoff error if num_samples is large
	double avg_x = 0.0;
	for ( int i=0; i<num_samples; i++ ) {
		avg_x += x[i];
	}
	avg_x /= static_cast<double>(num_samples);

	return avg_x;
}


double Statistics::variance(const std::vector<double>& x, const int delta_dof)
{
	int num_samples = x.size();
	int num_samples_final = num_samples - delta_dof;  // adjusted number of samples
	if ( num_samples_final < 1 ) {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  Insufficient number of samples\n";
		throw std::runtime_error( err_ss.str() );
	}

	double avg_x = Statistics::average(x);
	double var_x = 0.0;
	double dx_i;
	for ( int i=0; i<num_samples; i++ ) {
		dx_i = x[i] - avg_x;
		var_x += dx_i*dx_i;
	}

	// Normalize
	var_x /= static_cast<double>(num_samples_final);
	return var_x;
}


// Compute the sample variance divided by the sample mean
double Statistics::varianceOverAverage(const std::vector<double>& x)
{
	return ( Statistics::variance(x) )/( Statistics::average(x) );
}


// Compute the sample standard deviation 
double Statistics::std_dev(const std::vector<double>& x)
{
	return sqrt( Statistics::variance(x, 1) );
}


// Compute the sample kurtosis
// - Note that the estimator for non-normal X is consistent but sensitive to outliers
double Statistics::kurtosis(const std::vector<double>& x)
{
	double avg_x = Statistics::average(x);
	double var_x = Statistics::variance(x);

	int    num_samples  = x.size();
	double num_samples_d = static_cast<double>(num_samples);

	// Compute the fourth central moment (biased but consistent estimator)
	double mu_4_x = 0.0;
	for ( int i=0; i<num_samples; i++ )
	{
		mu_4_x += pow(x[i] - avg_x, 4);
	}
	mu_4_x *= (6.0*num_samples_d - 9.0)/(num_samples_d*(num_samples_d - 3.0) + 3.0);

	return mu_4_x/(var_x*var_x);
}


// Constructs a histogram of the data (data bins must be evenly spaced!)
void Statistics::constructHistogram(
		const std::vector<double>& x, const std::vector<double>& bins, 
		const bool isNormalized,
		// Output
		std::vector<double>& histogram, int& numOutliers)
{
	double dx = bins[1] - bins[0];
	double x_min = bins[0] - dx/2.0;

	int    numBins = bins.size();
	int    num_samples = x.size();

	// Place the data into bins
	histogram.resize(numBins, 0.0);
	numOutliers = 0;
	int bin;
	for ( int i=0; i<num_samples; ++i ) {
		bin = static_cast<int>( floor((x[i] - x_min)/dx) );
		if ( (bin >= 0) && (bin < numBins) ) { 
			histogram[bin] += 1.0; 
		}
		else { ++numOutliers; }
	}

	// TODO What if all samples are outliers?!

	// If the histogram is normalized, it's a probability distribution: p(x)
	// - Don't forget to divide by bin size!
	if ( isNormalized ) {
		double sum = 0.0;
		for ( int b=0; b<numBins; ++b ) { sum += histogram[b]; }
		for ( int b=0; b<numBins; ++b ) { histogram[b] /= sum*dx; }
	}
}



// Standard error of the mean assuming i.i.d. samples
double Statistics::std_err_mean(const std::vector<double>& x)
{
	// Standard deviation (unbiased for normal random variables)
	double std_dev = Statistics::std_dev(x);

	double num_samples_d = static_cast<double>( x.size() );
	return std_dev/sqrt(num_samples_d);
}



// Use the moving block bootstrap (MBB) method to estimate the std. error of a point statistic
void Statistics::movingBlockBootstrap_stdError_PointStatistic(
		const std::string& statistic, const std::vector<double>& x, const int num_samplesPerBlock, 
		const int numBootstrapSamples, 
		// Output
		double& bootstrapAverage, double& stdError, double& std_dev_stdError)
{
	// TODO Check block length vs. number of bootstrap samples and number of samples

	// Point estimator
	Statistics::PointEstmator estimator_fxn = this->getPointEstimatorMethod(statistic);

	// Allocate memory
	int 				num_samples = x.size();
	std::vector<double> x_bootstrap(num_samples);
	std::vector<double> bootstrapEstimates(numBootstrapSamples);

	for ( int i=0; i<numBootstrapSamples; ++i ) {
		this->getMovingBlockBootstrapSample(x, num_samplesPerBlock, x_bootstrap);

		bootstrapEstimates[i] = estimator_fxn(x_bootstrap);
	}

	// Average value of the statistic
	bootstrapAverage = Statistics::average(bootstrapEstimates);

	// Calculate the std. error of the statistic (and its std. deviation in turn)
	// - Note: With enough bootstrap samples, the central limit thm. applies
	// - Calculate c_N here (not in std_dev) to avoid duplicate computation
	double c_N = calc_correction_c_N(numBootstrapSamples);
	stdError = c_N*sqrt((Statistics::variance(bootstrapEstimates)));
	std_dev_stdError = sqrt(c_N*c_N - 1.0)*stdError; 
}



// Use the moving block bootstrap (MBB) to estimate the std. errors of a probability histogram
void Statistics::movingBlockBootstrap_stdError_p_x(
		const std::vector<double>& x, const std::vector<double>& bins, const int num_samplesPerBlock, 
		const int numBootstrapSamples, const bool doNegativeLog,
		// Output
		std::vector<double>& bootstrapAverages, std::vector<double>& std_devs) 
{
	// TODO Check block length vs. number of bootstrap samples and number of samples

	// Allocate memory
	int num_samples = x.size();
	int numBins    = bins.size();
	std::vector<double> x_bootstrap(num_samples);
	std::vector<double> p_x_bootstrap(numBins);

	bootstrapAverages.resize(numBins, 0.0);
	std_devs.resize(numBins, 0.0);

	// Number of bootstrap samples which are accepted into each bin
	// - Important for negative-log calculations, in order to track when
	//   a sample yields NaN for a bin and must be excluded
	std::vector<int> num_samplesPerBin(numBins, 0);

	// Perform MBB resampling
	int numOutliers;
	bool isNormalized = true;
	double sample;
	for ( int i=0; i<numBootstrapSamples; ++i ) {
		getMovingBlockBootstrapSample(x, num_samplesPerBlock, 
		                              x_bootstrap);

		constructHistogram(
				x_bootstrap, bins, isNormalized,
				// Output
				p_x_bootstrap, numOutliers);

		/*
		if ( doNegativeLog ) {
			for ( int b=0; b<numBins; ++b ) {
				p_x_bootstrap[b] = -log( p_x_bootstrap[b] );
			}
		}
		*/

		// Record the sample
		for ( int b=0; b<numBins; ++b ) {
			sample = p_x_bootstrap[b];
			if ( doNegativeLog ) {
				sample = -log(sample);
			}

			// Make sure the sample is finite before including it in the data set
			if ( std::isfinite(sample) ) {
				bootstrapAverages[b] += sample;  // sum of p(x)
				std_devs[b] += sample*sample;     // sum of p^2(x)
				++( num_samplesPerBin[b] );
			}
		}
	} // end loop over bootstrap samples

	// Finalize statistics
	double numBootstrapSamples_d;
	for ( int b=0; b<numBins; ++b ) {
		numBootstrapSamples_d = static_cast<double>(num_samplesPerBin[b]);

		bootstrapAverages[b] /= numBootstrapSamples_d; // <p(x)>

		// Std. dev. of bootstrap samples is a proxy for the expected error
		std_devs[b] /= numBootstrapSamples_d;         // <p^2(x)>
		std_devs[b] = sqrt( numBootstrapSamples_d/(numBootstrapSamples_d - 1.0)
		                   * (std_devs[b] - bootstrapAverages[b]*bootstrapAverages[b]) );
	}
}



// Perform block averaging and return the standard error of the mean
void Statistics::blockStatistics(
		const std::vector<double>& x, const int num_samplesPerBlock, 
		int& numBlocks, std::vector<double>& blockAverages, std::vector<double>& blockVariances,
 		double& stdError_avg, double& std_dev_stdError_avg)
{
	// Input checks
	if ( num_samplesPerBlock < 1 )
	{
		std::cerr << "Statistics.blockStatistics: Can't do block averages with less "
				  << "than 1 sample per block "
				  << "(input: numBlocks=" << numBlocks << ")." << "\n";
		exit(1);
	}

	int    num_samples  		   = x.size();
	double num_samples_d 		   = static_cast<double>(num_samples);
	double num_samplesPerBlock_d = static_cast<double>(num_samplesPerBlock);

	numBlocks         = static_cast<int>(num_samples_d/num_samplesPerBlock_d); // round down
	double numBlocks_d = static_cast<double>(numBlocks);

	if ( numBlocks < 2 )
	{
		std::cerr << "Statistics.blockStatistics: Insufficient data for at least "
				  << "two blocks of the same size. Decrease the block size " 
				  << "(input: num_samplesPerBlock=" << num_samplesPerBlock << ")."
				  << "\n";
		exit(1);
	}

	// Compute averages for each block
	blockAverages.assign(numBlocks, 0.0);
	blockVariances.assign(numBlocks, 0.0);
	double avg_x, avg_xSq;
	int    blockStart, blockStop;

	for ( int b=0; b<numBlocks; b++ )
	{
		blockStart = b*num_samplesPerBlock;
		blockStop  = blockStart + num_samplesPerBlock;

		avg_x   = 0.0;
		avg_xSq = 0.0;

		for ( int i=blockStart; i<blockStop; i++ )
		{
			avg_x   += x[i];	   // avg(x)
			avg_xSq += x[i]*x[i];  // avg(x^2)
		}

		avg_x   /= num_samplesPerBlock_d;
		avg_xSq /= num_samplesPerBlock_d;

		blockAverages[b]  = avg_x;
		blockVariances[b] = num_samplesPerBlock_d/(num_samplesPerBlock_d - 1.0)
								*(avg_xSq - avg_x*avg_x);
	}

	// Calculate the std. error (and its std. deviation/std. error)
	// - Note: the sample mean is normally distributed! Both of these estimators are unbiased
	// - Calculate c_N here (not in std_dev) to avoid duplicate computation
	double c_N = calc_correction_c_N(numBlocks);
	stdError_avg = c_N*sqrt((Statistics::variance(blockAverages))/numBlocks_d);
	std_dev_stdError_avg = sqrt(c_N*c_N - 1.0)*stdError_avg; 

	return;
}



// Performs bootstrap on subsamples whose elements are separated by the correlation time
void Statistics::timeCorrelationBootstrap(
		const std::string& statistic, const std::vector<double> x, const int correlationTime, 
		const int bootstrap_sample_size, const int numBootstrapSamples,
		// Output
		double& bootstrapAverage, double& bootstrapStdError)
{
	// Input checks
	int num_samples = x.size();
	if ( correlationTime >= num_samples ) {
		std::cerr << "Statistics.timeCorrelationBootstrap: Can't bootstrap a sample when the "
				  << "correlation time is longer than the data set (input: num_samples = "
				  << num_samples << ", correlation time = " << correlationTime << ")."
				  << "\n";
		exit(1);
	}
	if ( correlationTime == 0 ) { // FIXME Way to handle correlationTime = 0, 1 better?
		std::cerr << "Statistics.timeCorrelationBootstrap: The correlation time for your "
			<< "data is 0 [sample times]. Use a classical bootstrap method instead."
			<< "\n";
		exit(1);
	}

	if ( 10*bootstrap_sample_size < num_samples )
	{
		std::cerr << "Statistics.timeCorrelationBootstrap: The bootstrap sample size will be "
			<< "more than 10 times smaller than the original data set. Results may be "
			<< "inaccurate." 
			<< "\n";
	}

	// Point estimator
	Statistics::PointEstmator estimator_fxn = this->getPointEstimatorMethod(statistic);

	// RNG for selecting subsamples
	std::uniform_int_distribution<int> randomSubsampleOrigin(0, correlationTime-1);

	// Bootstrap estimates
	std::vector<double> bootstrapEstimates(numBootstrapSamples);
	std::vector<double> bootstrap_sample(bootstrap_sample_size);
	int subsampleOrigin;
	for ( int i=0; i<numBootstrapSamples; ++i ) {
		subsampleOrigin = randomSubsampleOrigin(*rng_ptr_);

		this->getBootstrapSampleFromTimeCorrelatedData(
				x, subsampleOrigin, correlationTime, bootstrap_sample_size, bootstrap_sample);

		bootstrapEstimates[i] = estimator_fxn(bootstrap_sample);
	}

	// Bootstrap-derived statistics
	bootstrapAverage  = Statistics::average(bootstrapEstimates);
	bootstrapStdError = Statistics::std_err_mean(bootstrapEstimates);

	return;
}



//----- Private -----//



// (Private) Calculate the correction factor c_N for the std. dev. of normal random variables
double Statistics::calc_correction_c_N(const int num_samples)
{
	if ( num_samples < 100 )
	{
		double num_samples_d = static_cast<double>(num_samples);
		double arg1 = (num_samples_d - 1.0)/2.0;
		double arg2 = num_samples_d/2.0;
		// TODO Check for overflow of lgamma? Still possible with large sample sizes...
		double c_N = exp(lgamma(arg1) - lgamma(arg2) + 0.5*log(arg1));
		return c_N;
	}
	else
	{
		// c_N(100) < 1.005 ==> negligible correction vs. uncertainty in data
		return 1.0;
	}
}



// (Private) Generates a bootstrap sample from the data
void Statistics::getBootstrapSample(
		const std::vector<double>& x, std::vector<double>& x_bootstrap)
{
	int num_samples = x.size();

	if ( static_cast<int>(x_bootstrap.size()) != num_samples )
	{
		x_bootstrap.resize(num_samples);
	}

	// Random int generator
	std::uniform_int_distribution<int> gen_int(0, num_samples-1);

	// Pick samples randomly (with replacement) to generate a bootstrap sample which
	// is as large as the original one
	int j;
	for ( int i=0; i<num_samples; ++i )
	{
		j = gen_int(*rng_ptr_);
		x_bootstrap[i] = x[j];
	}

	return;
}



// (Private) Generates a bootstrap sample of independent data from time-correlated data
void Statistics::getBootstrapSampleFromTimeCorrelatedData(
		const std::vector<double> data, const int subsampleOrigin,
		const int correlationTime, const int bootstrap_sample_size,
		std::vector<double>& bootstrap_sample)
{
	if ( (subsampleOrigin < 0) || (subsampleOrigin >= correlationTime) )
	{
		std::cerr << "Statistics.getBootstrapSampleFromTimeCorrelatedData: "
				  << "the parameter \"subsampleOrigin\" must be >= 0 and < correlationTime "
				  << "(input: subsampleOrigin = " << subsampleOrigin << ", " 
				  << "correlationTime = " << correlationTime << ")." 
				  << "\n";
		exit(1);
	}

	if ( static_cast<int>(bootstrap_sample.size()) != bootstrap_sample_size ) {
		bootstrap_sample.resize(bootstrap_sample_size);
	}

	// RNG for picking a random sample within the designated subsample
	int num_samples = data.size();
	int numIndependentSamplesPerTimeOrigin = num_samples/correlationTime; // rounds down
	std::uniform_int_distribution<int> randomSample(0, numIndependentSamplesPerTimeOrigin - 1);

	// Generate the bootstrap sample
	int j;
	for ( int i=0; i<bootstrap_sample_size; ++i ) {
		j = randomSample(*rng_ptr_);
		bootstrap_sample[i] = data[subsampleOrigin + j*correlationTime];
	}

	return;
}



// (Private) Generate a moving block bootstrap sample from the data
void Statistics::getMovingBlockBootstrapSample(
		const std::vector<double>& x, const int num_samplesPerBlock, 
		std::vector<double>& x_bootstrap)
{
	int num_samples = x.size(); // "n"
	int numBlocks  = num_samples - num_samplesPerBlock + 1; // "N"

	if ( static_cast<int>(x_bootstrap.size()) != num_samples ) {
		x_bootstrap.resize(num_samples);
	}

	double num_samples_d 				       = static_cast<double>(num_samples);
	double num_samplesPerBlock_d 		   = static_cast<double>(num_samplesPerBlock); // "l"
	int    numBlocksPerBootstrapSample = static_cast<int>(num_samples_d/num_samplesPerBlock_d); // "b"

	// Random int generator for choosing blocks
	std::uniform_int_distribution<int> gen_int(0, numBlocks-1);	// Interface to generate ints

	// Pick blocks randomly (with replacement) to generate a bootstrap sample set which
	// is as large as the original sample set
	int blockStart, offset;
	for ( int j=0; j<numBlocksPerBootstrapSample; ++j ) {
		// Block number = index in 'x' where it begins
		blockStart = gen_int(*rng_ptr_);

		// Build the bootstrap sample set
		offset = j*num_samplesPerBlock;	// where to start placing data in bootstrap sample
		for ( int k=0; k<num_samplesPerBlock; ++k ) {
			x_bootstrap[offset + k] = x[blockStart + k];
		}
	}

	return;
}



// (Private) Returns a ptr to the member function which estimates the indicated point statistic
Statistics::PointEstmator Statistics::getPointEstimatorMethod(
		const std::string& statistic)
{
	// TODO Use a map instead
	PointEstmator estimator_fxn;
	if ( (statistic == "average") || (statistic == "avg") )  		
	{ 
		//estimator_fxn = [=](const std::vector<double>& x) { return Statistics::average(x); };
		estimator_fxn = &Statistics::average;
	}
	else if ( (statistic == "variance") || (statistic == "var") ) 
	{ 
		// Bind with default delta_dof (implicit)
		estimator_fxn = [=](const std::vector<double>& x) { return Statistics::variance(x); };
	}
	else if ( (statistic == "std_deviation") 
				|| (statistic == "std deviation")
				|| (statistic == "standard deviation")
				|| (statistic == "std dev")
				|| (statistic == "standardDeviation") 	
				|| (statistic == "std_dev") )
	{
		estimator_fxn = &Statistics::std_dev;
	}
	else if ( statistic == "kurtosis" ) 
	{ 
		estimator_fxn = &Statistics::kurtosis; 
	}
	else if ( (statistic == "variance_over_average") || (statistic == "var_over_avg") )
	{
		estimator_fxn = &Statistics::varianceOverAverage;
	}
	else
	{
		std::cerr << "Statistics.getPointEstimatorMethod: Statistic not recognized/supported "
				  << "(input: statistic = " << statistic << ")." << "\n";
		exit(1);
	}

	return estimator_fxn;
}
