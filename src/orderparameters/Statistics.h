/* Statistics.h
 *
 * ABOUT: Object with statistics routines
 * DEVELOPMENT: TODO
 *  - Swich to throwing exceptions
 *  - Rename variables to be consistent with prevailing style
 */

// Standard headers
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "utils.h"

#ifndef STATISTICS_H
#define STATISTICS_H

class Statistics
{
  public:
	Statistics(const int seed);
	~Statistics();

	static double average(const std::vector<double>& x);

	// Computes the variance of a set of data
	// - Default: delta_dof = 1 (sample variance) --> unbiased estimator
	static double variance(const std::vector<double>& x, const int delta_dof = 1);

	static double varianceOverAverage(const std::vector<double>& x); 

	static double std_dev(const std::vector<double>& x);

	static double kurtosis(const std::vector<double>& x);

	static double covariance(const std::vector<double>& x, const std::vector<double>& y); // TODO

	// Constructs a histogram of the data TODO
	// - Note: data bins must be evenly spaced
	void constructHistogram(
			const std::vector<double>& x,            // Data
			const std::vector<double>& bins,		 // Bin centers
			const bool 				   isNormalized, // Option to normalize the histogram
			// Output
			std::vector<double>& histogram,
			int&                 numOutliers         // # data pts. which didn't fit in the bins
	);

	// Standard error of the mean (assuming i.i.d. samples)
	static double std_err_mean(const std::vector<double>& x);

	// Uses block averaging to compute the standard error and its standard deviation
	void blockStatistics(
			const std::vector<double>& x, 
			const int numSamplesPerBlock,
			// Output
			int& numBlocks,
			std::vector<double>& blockAverages,
			std::vector<double>& blockVariances,
			double& stdError_avg,       // Standard error of the average
			double&	std_dev_stdError_avg // Std. deviation of the above
	);

	// Use moving block bootstrap to estimate the std. error of a point statistic and its std. dev. 
	// - The point statistic should be expressable as an expected value 
	void movingBlockBootstrap_stdError_PointStatistic(
			const std::string&         statistic,             // Point statistic of 'x' to evaluate
			const std::vector<double>& x,                     // Data
			const int 				   numSamplesPerBlock,
			const int				   numBootstrapSamples,
			// Output
			double&                    bootstrapAverage,      // Average across boostrap samples
			double&					   stdError,              // Std. error of the avg. statistic
			double&					   std_dev_stdError        // Std. dev. of the std. error
	);

	// Use moving block bootstrap to estimate the std. error of the probability distribution p(x) TODO
	void movingBlockBootstrap_stdError_p_x(
			const std::vector<double>& x,                    // Data
			const std::vector<double>& bins,                 // Histogram bins
			const int                  numSamplesPerBlock,
			const int                  numBootstrapSamples,
			const bool                 doNegativeLog,        // Do -log(p(x)) instead
			// Output
			std::vector<double>&       bootstrapAverages,    // Averages across boostrap samples
			std::vector<double>&       std_devs               // Std. errors of <p(x)>
	);

	// Perform bootstrap on subsamples whose elements are separated by the correlation time
	// - Subsamples chosen this way contain independent data
	// - Each bootstrap sample should only be as large as the independent subsample
	//   - Make sure you have enough data that a subsample can be used to produce a reliable 
	//     estimate of the target statistic!
	void timeCorrelationBootstrap(
			const std::string&        statistic,			// Statistic of 'x' to evaluate
			const std::vector<double> x,					// Data
			const int				  correlationTime,		// Correlation time for x [samples]
			const int				  bootstrapSampleSize,
			const int				  numBootstrapSamples,
			// Output
			double& bootstrapAverage,	// Average of sample estimates
			double& bootstrapStdError	// Standard error of the bootstrap estimates
	);

	/*
	// Use moving block bootstrap to estimate the std. errors of a probability distribution TODO
	// - Note: Bins must be evenly spaced
	void blockBootstrap_stdError_p_x(
			const std::vector<double>& x,					// Data
			const std::vector<double>& bins,				// Bins for p(x) histogram
			const int 				   numSamplesPerBlock,
			const int 				   numBootstrapSamples,
			// Output
			std::vector<double>& stdErrors, 
			std::vector<double>& std_dev_stdErrors
	);
	*/

	// Moving block bootstrap for std. error of the mean with auto-block size selection TODO
	void blockBootstrap_stdError_avg_auto(
			const std::vector<double>& x,
			const int 				   blockSizeGuess,
			const int				   numSubsamples,
			const int				   numBootstrapSamples,
			// Output
			double&					   stdError_avg,
			double&					   std_dev_stdError_avg
	);

 private:

	//---- Typedefs -----//

	// Member function pointer for flexibility in bootstrapping algorithms
	// - A single bootstrap routine can work for many different point statistics
	// - TODO switch to std::function
	//using MethodPointer_PointEstimator = double (*)(const std::vector<double>&);
	using PointEstmator = std::function<double(const std::vector<double>&)>;

	// Random number generator
	using MersenneTwisterPtr = std::unique_ptr<std::mt19937>;
	MersenneTwisterPtr rng_ptr_;

	// Calculates the correction factor c_N for the std. dev. of normal random variables
	// - Multiplying the sample standard deviation by this factor makes it an unbiased
	//   estimator for normal random variables
	double calc_correction_c_N(const int numSamples);

	// Generates a bootstrap sample from the data
	// - Method assumes independent samples
	void getBootstrapSample(
			const std::vector<double>& x,
			std::vector<double>&       x_bootstrap
	);

	// Generates a bootstrap sample of independent data from time-correlated data
	// - A "subsample of independent data" is defined here as a sequence of points 
	//   which are separated by the correlation time.
	// - The parameter 0 <= subsampleOrigin < correlationTime is the index of the
	//   first point in the subsample which will be used to produce the bootstrap
	//   sample
	void getBootstrapSampleFromTimeCorrelatedData(
			const std::vector<double> data,
			const int				  subsampleOrigin,
			const int				  correlationTime,
			const int				  bootstrapSampleSize,
			// Output
			std::vector<double>& bootstrapSample
	);

	// Generates a moving block bootstrap sample from the data
	void getMovingBlockBootstrapSample(
			const std::vector<double>& x, 
			const int 				   numSamplesPerBlock, // block size
			// Output
			std::vector<double>& 	   x_bootstrap
	);

	// Estimates the optimal block size for computing the std. error of the mean TODO
	int estimateOptimalBlockSize_stdError_avg(
			const std::vector<double>& x,
			const int 				   blockSizeGuess,  // Initial guess for block size
			const int 				   numSubsamples	// Num. of subsamples for MSE estimate
	);

	// Returns a ptr to the member function which estimates the indicated point statistic
	PointEstmator getPointEstimatorMethod(const std::string& statistic);
};

#endif // STATISTICS_H
