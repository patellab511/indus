#include "Bias.h"

Bias::Bias(const ParameterPack& input_parameter_pack, const SimulationBox& simulation_box)
 : simulation_box_(simulation_box)
{
	setParameters(input_parameter_pack);
}

void Bias::setParameters(const ParameterPack& input_parameter_pack)
{
	using KeyType = ParameterPack::KeyType;

	// Harmonic
	x_star_ = kappa_ = 0.0;
	input_parameter_pack.readNumber("kappa",  KeyType::Optional, kappa_);
	input_parameter_pack.readNumber("x_star", KeyType::Optional, x_star_);

	// Linear
	phi_ = constant_ = 0.0;
	input_parameter_pack.readNumber("phi",      KeyType::Optional, phi_);
	input_parameter_pack.readNumber("constant", KeyType::Optional, constant_);

	// Left one-sided harmonic
	x_left_ = k_left_  = 0.0;
	input_parameter_pack.readNumber("x_left", KeyType::Optional, x_left_);
	input_parameter_pack.readNumber("k_left", KeyType::Optional, k_left_);

	// Right one-sided harmonic
	x_right_ = k_right_ = 0.0;
	input_parameter_pack.readNumber("x_right", KeyType::Optional, x_right_);
	input_parameter_pack.readNumber("k_right", KeyType::Optional, k_right_);
}


void Bias::accumulateBias(
		const double x, const bool want_derivatives, 
		const std::vector<Bias::Real3>& derivatives_x, const Matrix& sum_r_cross_dx_dr,
		// Output
		double& u_bias_x, double& u_bias_total, 
		std::vector<Real3>& derivatives_u_bias_total, Matrix& virial_u_bias_total
) const
{
	u_bias_x = 0.0;
	double du_bias_dx = 0.0;

	// Harmonic
	double dx = x - x_star_;
	u_bias_x += 0.5*kappa_*dx*dx;
	if ( want_derivatives == true ) {
		du_bias_dx += kappa_*dx;
	}

	// Linear
	u_bias_x += phi_*x + constant_;
	if ( want_derivatives == true ) {
		du_bias_dx += phi_;
	}

	// Left one-sided harmonic
	dx = x - x_left_;
	if ( dx < 0.0 ) {
		u_bias_x += 0.5*k_left_*dx*dx;
		if ( want_derivatives == true ) {
			du_bias_dx += k_left_*dx;
		}
	}

	// Right one-sided harmonic
	dx = x - x_right_;
	if ( dx > 0.0 ) {
		u_bias_x += 0.5*k_right_*dx*dx;
		if ( want_derivatives == true ) {
			du_bias_dx += k_right_*dx;
		}
	}

	// Add to running sum of total bias
	u_bias_total += u_bias_x;

	if ( want_derivatives ) {
		// Derivatives wrt. atom positions
		int num_biased_atoms = derivatives_x.size();
		for ( int i=0; i<num_biased_atoms; ++i ) {
			for ( int d=0; d<DIM_; ++d ) {
				derivatives_u_bias_total[i][d] += du_bias_dx*derivatives_x[i][d];
			}
		}

		// Virial
		for ( int a=0; a<DIM_; ++a ) {	
			for ( int b=0; b<DIM_; ++b ) {
				virial_u_bias_total[a][b] += 0.5*du_bias_dx*sum_r_cross_dx_dr[a][b];
			}
		}
	}
}


std::string Bias::getInputSummary(const std::string& prepend_string /*= "# "*/) const
{
	// Don't report parameters for terms which are guaranteed to be zero
	std::stringstream ss;
	ss << prepend_string << "Bias_parameters\n";
	if ( kappa_ != 0.0 ) {
		ss << prepend_string << "  Harmonic\n"
		   << prepend_string << "    x_star = " << x_star_ << " [units of x]\n"
		   << prepend_string << "    kappa  = " << kappa_  << " [kJ/mol]\n";
	}
	if ( phi_ != 0.0 or constant_ != 0.0 ) {
		ss << prepend_string << "  Linear\n"
		   << prepend_string << "    phi       = " << phi_      << " [kJ/mol]\n"
		   << prepend_string << "    constant_ = " << constant_ << " [kJ/mol]\n";
		
	}
	if ( k_left_ != 0.0 ) {
		ss << prepend_string << "  Left_one-sided_harmonic\n"
		   << prepend_string << "    x_left = " << x_left_ << " [units of x]\n"
		   << prepend_string << "    k_left = " << k_left_ << " [kJ/mol]\n";
	}
	if ( k_right_ != 0.0 ) {
		ss << prepend_string << "  Right_one-sided_harmonic\n"
		   << prepend_string << "    x_right = " << x_right_ << " [units of x]\n"
		   << prepend_string << "    k_right = " << k_right_ << " [kJ/mol]\n";
	}

	return ss.str();
}
