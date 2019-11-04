#include "Bias.h"

Bias::Bias(const ParameterPack& input_parameter_pack, const SimulationState& simulation_state)
 : simulation_state_( simulation_state )
{
	setParameters(input_parameter_pack);
}

void Bias::setParameters(const ParameterPack& input_parameter_pack)
{
	using KeyType = ParameterPack::KeyType;

	//----- Construct terms in the bias -----//

	potential_ptrs_.clear();

	// FIXME 
	// This is still kinda cludgy, querying key values to determine terms present
	// However, doing much more would require significant refactoring, and probably
	// a change in input structure
	bool found;
	double param;

	// Harmonic
	found = input_parameter_pack.readNumber("kappa", KeyType::Optional, param);
	if ( found ) {
		auto tmp_ptr = std::unique_ptr<Potential>( 
			new HarmonicPotential(input_parameter_pack, simulation_state_) 
		);
		potential_ptrs_.push_back( std::move(tmp_ptr) );
	}

	// Linear
	found = input_parameter_pack.readNumber("phi", KeyType::Optional, param);
	Range range;
	found = (found or input_parameter_pack.readArray("phi_range", KeyType::Optional, range));
	if ( found ) {
		auto tmp_ptr = std::unique_ptr<Potential>( 
			new LinearPotential(input_parameter_pack, simulation_state_)
		);
		potential_ptrs_.push_back( std::move(tmp_ptr) );
	}

	// Left one-sided harmonic
	found = input_parameter_pack.readNumber("k_left", KeyType::Optional, param);
	if ( found ) {
		auto tmp_ptr = std::unique_ptr<Potential>( 
			new LeftHarmonicPotential(input_parameter_pack, simulation_state_)
		);
		potential_ptrs_.push_back( std::move(tmp_ptr) );
	}

	// Right one-sided harmonic
	found = input_parameter_pack.readNumber("k_right", KeyType::Optional, param);
	if ( found ) {
		auto tmp_ptr = std::unique_ptr<Potential>( 
			new RightHarmonicPotential(input_parameter_pack, simulation_state_)
		);
		potential_ptrs_.push_back( std::move(tmp_ptr) );
	}

	// Ramp (optional)
	// TODO Turn into Ramp input object?
	bool ramp_bias = false;
	input_parameter_pack.readFlag("ramp_bias", KeyType::Optional, ramp_bias);
	if ( ramp_bias ) {
		CommonTypes::Range t_range_tmp;
		input_parameter_pack.readArray("t_range_ramp", KeyType::Required, t_range_tmp);

		bias_ramp_ptr_ = std::unique_ptr<Ramp>( new LinearRamp(input_parameter_pack, simulation_state_) );
	}
}


// Returns u_bias
double Bias::calculate( const double x, const bool want_derivative, 
                        double& u_bias, double& du_bias_dx ) const
{
	// Reset bias
	u_bias = 0.0;
	if ( want_derivative ) {
		du_bias_dx = 0.0;
	}

	// Accumulate terms
	double u_bias_k, du_bias_k_dx;
	int num_terms = potential_ptrs_.size();
	for ( int k=0; k<num_terms; ++k ) {
		potential_ptrs_[k]->evaluate(x, u_bias_k, du_bias_k_dx);

		u_bias += u_bias_k;
		if ( want_derivative ) {
			du_bias_dx += du_bias_k_dx;
		}
	}

	if ( bias_ramp_ptr_ != nullptr ) { 
		// Scale entire bias by coeff(t)
		double coeff = bias_ramp_ptr_->calculate();

		u_bias *= coeff;
		if ( want_derivative ) {
			du_bias_dx *= coeff;
		}
	}

	return u_bias;
}


void Bias::applyBias(
	const double x, const bool want_derivatives,
	const std::vector<Real3>& derivatives_x, const Matrix& sum_r_cross_dx_dr,
	double& u_bias_x, std::vector<Real3>& derivatives_u_bias_x, Matrix& virial_u_bias_x
) const
{
	double du_bias_dx;
	this->calculate(x, want_derivatives, u_bias_x, du_bias_dx);

	if ( want_derivatives ) {
		// Derivatives wrt. atom positions
		int num_biased_atoms = derivatives_x.size();
		derivatives_u_bias_x.resize(num_biased_atoms);
		#pragma omp parallel for schedule(static, 8)
		for ( int i=0; i<num_biased_atoms; ++i ) {
			for ( int d=0; d<DIM_; ++d ) {
				derivatives_u_bias_x[i][d] = du_bias_dx*derivatives_x[i][d];
			}
		}

		// Virial
		for ( int a=0; a<DIM_; ++a ) {	
			for ( int b=0; b<DIM_; ++b ) {
				virial_u_bias_x[a][b] = 0.5*du_bias_dx*sum_r_cross_dx_dr[a][b];
			}
		}
	}
}


void Bias::accumulateBias(
		const double x, const bool want_derivatives, 
		const std::vector<Bias::Real3>& derivatives_x, const Matrix& sum_r_cross_dx_dr,
		// Output
		double& u_bias_x, double& u_bias_total, 
		std::vector<Real3>& derivatives_u_bias_total, Matrix& virial_u_bias_total
) const
{
	double du_bias_dx;
	this->calculate(x, want_derivatives, u_bias_x, du_bias_dx);

	// Add to running sum of total bias
	u_bias_total += u_bias_x;

	if ( want_derivatives ) {
		// Derivatives wrt. atom positions
		int num_biased_atoms = derivatives_x.size();
		#pragma omp parallel for schedule(static, 8)
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
	std::stringstream ss;
	ss << prepend_string << "Bias_parameters\n";

	int num_terms = potential_ptrs_.size();
	std::string indent = "  ";
	std::string prepend_term = prepend_string + indent;
	for ( int k=0; k<num_terms; ++k ) {
		ss << potential_ptrs_[k]->getInputSummary(prepend_term);
	}

	if ( bias_ramp_ptr_ != nullptr ) {
		// TODO Move this section to functions in Ramp and its derived classes?
		const CommonTypes::Range& t_range_ramp     = bias_ramp_ptr_->get_t_range();
		const CommonTypes::Range& coeff_range_ramp = bias_ramp_ptr_->get_parameter_range();

		ss << prepend_string << "  Ramp bias with a changing prefactor\n"
		   << prepend_string << "    type        = linear\n"
		   << prepend_string << "    time_range  = { " << t_range_ramp[0]     << ", " << t_range_ramp[1]     << " } [ps]\n"
		   << prepend_string << "    coeff_range = { " << coeff_range_ramp[0] << ", " << coeff_range_ramp[1] << " }\n";
	}

	return ss.str();
}


Ramp::Ramp(const ParameterPack& input_parameter_pack, const SimulationState& simulation_state):
	simulation_state_(simulation_state),
	t_range_{{0.0, 0.0}},
	parameter_range_{{0.0, 1.0}}
{
	using KeyType = ParameterPack::KeyType;

	input_parameter_pack.readArray("t_range_ramp", KeyType::Required, t_range_);
	if ( t_range_[0] > t_range_[1] ) {
		std::stringstream err_ss;
		err_ss << "Error in Ramp (for Bias): for the ramp-up period, t0 must be less than tf "
		       << "(input: t0 = " << t_range_[0] << " ps, tf = " << t_range_[1]
		       << " ps)\n";
		throw std::runtime_error( err_ss.str() );
	}

	input_parameter_pack.readArray("coeff_range_ramp", KeyType::Optional, parameter_range_);
}


Ramp::Ramp(const Range& t_range, const Range& parameter_range, 
           const SimulationState& simulation_state):
	simulation_state_(simulation_state),
	t_range_(t_range),
	parameter_range_(parameter_range)
{}


LinearRamp::LinearRamp(
	const ParameterPack& input_parameter_pack, const SimulationState& simulation_state
):
	Ramp(input_parameter_pack, simulation_state),
	slope_(calculate_slope())
{}

LinearRamp::LinearRamp(const Range& t_range, const Range& parameter_range, 
                       const SimulationState& simulation_state):
	Ramp(t_range, parameter_range, simulation_state),
	slope_(calculate_slope())
{}

double LinearRamp::calculate_slope() const
{
	return (parameter_range_[1] - parameter_range_[0]) / (t_range_[1] - t_range_[0]);
}

double LinearRamp::calculate() const 
{
	double p;
	const double time = simulation_state_.get_time();
	if ( time < t_range_[0] ) {
		p = parameter_range_[0];
	}
	else if ( time > t_range_[1] ) {
		p = parameter_range_[1];
	}
	else {
		p = slope_*(time - t_range_[0]) + parameter_range_[0];
	}

	return p;
}
