#include "Potential.h"

#include "Bias.h"  // defines the Ramp

//----- Potential base class -----//

Potential::Potential(const ParameterPack& input_pack, const SimulationState& simulation_state):
	simulation_state_(simulation_state)
{}

std::unique_ptr<Ramp> Potential::createRamp(
	const Range& t_range, const Range& parameter_range
) const
{
	// Default: make a linear ramp
	return std::unique_ptr<Ramp>( new LinearRamp(t_range, parameter_range, simulation_state_) );
}

const std::string Potential::t_range_input_key = "t_range_ramp";


//----- Harmonic potential -----//

HarmonicPotential::HarmonicPotential(
	const ParameterPack& input_pack, const SimulationState& simulation_state
):
	Potential(input_pack, simulation_state)
{
	using KeyType = ParameterPack::KeyType;

	input_pack.readNumber("kappa", KeyType::Required, kappa_);

	bool found_range = input_pack.readArray("x_star_range", KeyType::Optional, x_star_range_);
	if ( found_range ) {
		// x*(t)
		input_pack.readArray(t_range_input_key, KeyType::Required, t_range_ramp_);
		ramp_ptr_ = createRamp(t_range_ramp_, x_star_range_);
	}
	else {
		// x* is constant
		input_pack.readNumber("x_star", KeyType::Required, x_star_);
	}
}

void HarmonicPotential::evaluate(
		const double x, double& u_bias, double& du_bias_dx) const
{
	// Location of bias: x*(t)
	double x_star_t = x_star_;
	if ( ramp_ptr_ != nullptr ) {
		x_star_t = ramp_ptr_->calculate();
	}

	double dx = x - x_star_t;
	u_bias     = 0.5*kappa_*dx*dx;
	du_bias_dx = kappa_*dx;
}

std::string HarmonicPotential::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Harmonic\n";
	if ( ramp_ptr_ == nullptr ) {
		ss << prepend_string << "  x_star = " << x_star_ << " [units of x]\n";
	}
	else {
		ss << prepend_string << "  x_star_range = [ " << x_star_range_[0] << "  " << x_star_range_[1] << " ]" << " [units of x]\n";
		ss << prepend_string << "  t_range_ramp = [ " << t_range_ramp_[0] << "  " << t_range_ramp_[1] << " ]" << " [ps]\n";
	}
	
	ss << prepend_string << "  kappa  = " << kappa_  << " [kJ/mol]\n";
	return ss.str();
}


//----- Linear potential -----//

LinearPotential::LinearPotential(
	const ParameterPack& input_pack, const SimulationState& simulation_state
):
	Potential(input_pack, simulation_state)
{
	using KeyType = ParameterPack::KeyType;

	input_pack.readNumber("constant", KeyType::Optional, c_);  // FIXME

	bool found_range = input_pack.readArray("phi_range", KeyType::Optional, phi_range_);
	if ( found_range ) {
		// phi(t)
		input_pack.readArray(t_range_input_key, KeyType::Required, t_range_ramp_);
		ramp_ptr_ = createRamp(t_range_ramp_, phi_range_);
	}
	else {
		// phi is constant
		input_pack.readNumber("phi", KeyType::Required, phi_);
	}
}

void LinearPotential::evaluate(const double x, double& u_bias, double& du_bias_dx) const
{
	// Slope of bias: phi(t)
	double phi_t = phi_;
	if ( ramp_ptr_ != nullptr ) {
		phi_t = ramp_ptr_->calculate();
	}

	u_bias     = phi_t*x + c_;
	du_bias_dx = phi_t;
}

std::string LinearPotential::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Linear\n";
	if ( ramp_ptr_ == nullptr ) {
		ss << prepend_string << "  phi      = " << phi_ << " [kJ/mol]\n";
	}
	else {
		ss << prepend_string << "  phi_range    = [ " << phi_range_[0]    << "  " << phi_range_[1]    << " ]" << " [kJ/mol]\n";
		ss << prepend_string << "  t_range_ramp = [ " << t_range_ramp_[0] << "  " << t_range_ramp_[1] << " ]" << " [ps]\n";
	}
	
	ss << prepend_string << "  constant = " << c_  << " [kJ/mol]\n";
	return ss.str();
}


//----- Left one-sided harmonic potential -----//

LeftHarmonicPotential::LeftHarmonicPotential(
	const ParameterPack& input_pack, const SimulationState& simulation_state
):
	Potential(input_pack, simulation_state)
{
	using KeyType = ParameterPack::KeyType;

	input_pack.readNumber("k_left", KeyType::Required, k_left_);

	bool found_range = input_pack.readArray("x_left_range", KeyType::Optional, x_left_range_);
	if ( found_range ) {
		// x_left(t)
		input_pack.readArray(t_range_input_key, KeyType::Required, t_range_ramp_);
		ramp_ptr_ = createRamp(t_range_ramp_, x_left_range_);
	}
	else {
		// x_left is constant
		input_pack.readNumber("x_left", KeyType::Required, x_left_);
	}
}

void LeftHarmonicPotential::evaluate(const double x, double& u_bias, double& du_bias_dx) const
{
	// Location of bias: x_left(t)
	double x_left_t = x_left_;
	if ( ramp_ptr_ != nullptr ) {
		x_left_t = ramp_ptr_->calculate();
	}

	double dx = x - x_left_t;
	if ( dx < 0.0 ) {
		u_bias     = 0.5*k_left_*dx*dx;
		du_bias_dx = k_left_*dx;
	}
	else {
		u_bias     = 0.0;
		du_bias_dx = 0.0;
	}
}

std::string LeftHarmonicPotential::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Left_one-sided_harmonic\n";
	if ( ramp_ptr_ == nullptr ) {
		ss << prepend_string << "  x_left = " << x_left_ << " [units of x]\n";
	}
	else {
		ss << prepend_string << "  x_left_range = [ " << x_left_range_[0] << "  " << x_left_range_[1] << " ]" << " [units of x]\n";
		ss << prepend_string << "  t_range_ramp = [ " << t_range_ramp_[0] << "  " << t_range_ramp_[1] << " ]" << " [ps]\n";
	}
	
	ss << prepend_string << "  k_left = " << std::setprecision(9) << k_left_  << " [kJ/mol]\n";
	return ss.str();
}

//----- Right one-sided harmonic potential -----//

RightHarmonicPotential::RightHarmonicPotential(
	const ParameterPack& input_pack, const SimulationState& simulation_state
):
	Potential(input_pack, simulation_state)
{
	using KeyType = ParameterPack::KeyType;

	input_pack.readNumber("k_right", KeyType::Required, k_right_);

	bool found_range = input_pack.readArray("x_right_range", KeyType::Optional, x_right_range_);
	if ( found_range ) {
		// x_right(t)
		input_pack.readArray(t_range_input_key, KeyType::Required, t_range_ramp_);
		ramp_ptr_ = createRamp(t_range_ramp_, x_right_range_);
	}
	else {
		// x_right is constant
		input_pack.readNumber("x_right", KeyType::Required, x_right_);
	}
}
void RightHarmonicPotential::evaluate(const double x, double& u_bias, double& du_bias_dx) const
{
	// Location of bias: x_right(t)
	double x_right_t = x_right_;
	if ( ramp_ptr_ != nullptr ) {
		x_right_t = ramp_ptr_->calculate();
	}

	double dx = x - x_right_t;
	if ( dx > 0.0 ) {
		u_bias     = 0.5*k_right_*dx*dx;
		du_bias_dx = k_right_*dx;
	}
	else {
		u_bias     = 0.0;
		du_bias_dx = 0.0;
	}
}
std::string RightHarmonicPotential::getInputSummary(const std::string& prepend_string) const
{
	std::stringstream ss;
	ss << prepend_string << "Right_one-sided_harmonic\n";
	if ( ramp_ptr_ == nullptr ) {
		ss << prepend_string << "  x_right = " << x_right_ << " [units of x]\n";
	}
	else {
		ss << prepend_string << "  x_right_range = [ " << x_right_range_[0] << "  " << x_right_range_[1] << " ]" << " [units of x]\n";
		ss << prepend_string << "  t_range_ramp  = [ " << t_range_ramp_[0]  << "  " << t_range_ramp_[1]  << " ]" << " [ps]\n";
	}
	
	ss << prepend_string << "  k_right = " << k_right_  << " [kJ/mol]\n";
	return ss.str();
}
