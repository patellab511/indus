#include "NeighborSphere.h"

NeighborSphere::NeighborSphere( const double radius, const double sigma, const double alpha_c, 
                                const bool need_derivatives )
	: r_(radius), sigma_(sigma), alpha_c_(alpha_c),
		switching_function_r_(r_, sigma, alpha_c),
	  need_derivatives_(need_derivatives)
{}

NeighborSphere::NeighborSphere( 
		const ParameterPack& input_parameter_pack, const bool need_derivatives 
)
 : r_(0.35), sigma_(0.01), alpha_c_(0.02),
   switching_function_r_(r_, sigma_, alpha_c_),
   need_derivatives_(need_derivatives)
{
	// Parse input
	using KeyType = ParameterPack::KeyType;
	input_parameter_pack.readNumber("sigma",   KeyType::Optional, sigma_);
	input_parameter_pack.readNumber("alpha_c", KeyType::Optional, alpha_c_);
	input_parameter_pack.readNumber("radius",  KeyType::Optional, r_);

	// Update parameters
	switching_function_r_.setCoarseGrainingParameters(sigma_, alpha_c_);
	switching_function_r_.set_r_max(r_);
}


void NeighborSphere::calculateIndicator(
		const NeighborSphere::Real3& r_ij, const double dist_ij,
		double& h_nn, double& htilde_nn, NeighborSphere::Real3& dhtilde_nn_dr_j
) const
{
	double dhtilde_r_dr;
	switching_function_r_.calculate( dist_ij, need_derivatives_,
	                                 h_nn, htilde_nn, dhtilde_r_dr );

	if ( need_derivatives_ ) {
		if ( dist_ij > 0.0 ) {
			dhtilde_r_dr /= dist_ij;  // Reuse variable to precompute common prefactor
			for ( int k=0; k<DIM_; ++k ) { dhtilde_nn_dr_j[k] = dhtilde_r_dr*r_ij[k]; }
		} 
		else {
			// Note: when r_i = r_j, h_nn is ill-defined
			for ( int k=0; k<DIM_; ++k ) { 
				dhtilde_nn_dr_j[k] = 0.0; 
			}
		}
	}
}


// Return a summary of settings
std::string NeighborSphere::getInputSummary(const std::string& prefix) const
{
	std::stringstream ss;
	ss << prefix << "Neighbor_sphere \n"
	   << prefix << "  radius  = " << r_       << " [nm]\n"
	   << prefix << "  sigma   = " << sigma_   << " [nm]\n"
	   << prefix << "  alpha_c = " << alpha_c_ << " [nm]\n";

	return ss.str();
}
