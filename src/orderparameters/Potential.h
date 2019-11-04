// Potential: a term in a Bias

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <memory>

#include "CommonTypes.h"
#include "InputParser.h"
#include "SimulationState.h"

class Ramp;

class Potential
{
 public:
	using Range = CommonTypes::Range;

	Potential(
		const ParameterPack&   input_pack,
		const SimulationState& simulation_state
	);

	virtual ~Potential() {};

	// Returns the value of the bias (in kJ/mol)
	virtual void evaluate(
		const double x,
		// Output
		double& u_bias,
		double& du_bias_dx
	) const = 0;

	virtual std::string getInputSummary(const std::string& prepend_string) const { return ""; };

 protected:
	const SimulationState& simulation_state_;  // owned by the driver

	std::unique_ptr<Ramp> ramp_ptr_ = nullptr;
	Range t_range_ramp_;
	static const std::string t_range_input_key;

	virtual std::unique_ptr<Ramp> createRamp(
		const Range& t_range, 
		const Range& parameter_range
	) const;
};


// Harmonic potential
class HarmonicPotential : public Potential
{
 public:
	HarmonicPotential(const ParameterPack& input_pack, const SimulationState& simulation_state);
	virtual void evaluate(const double x, double& u_bias, double& du_bias_dx) const override;
	virtual std::string getInputSummary(const std::string& prepend_string) const override;

 private:
	double x_star_, kappa_;
	Range  x_star_range_;
};


// Linear potential (with a constant)
class LinearPotential : public Potential
{
 public:
	LinearPotential(const ParameterPack& input_pack, const SimulationState& simulation_state);
	virtual void evaluate(const double x, double& u_bias, double& du_bias_dx) const override;
	virtual std::string getInputSummary(const std::string& prepend_string) const override;

 private:
	double phi_ = 0.0;
	double c_   = 0.0;
	Range  phi_range_;
};


// Left (one-sided) harmonic potential
class LeftHarmonicPotential : public Potential
{
 public:
	LeftHarmonicPotential(const ParameterPack& input_pack, const SimulationState& simulation_state);
	virtual void evaluate(const double x, double& u_bias, double& du_bias_dx) const override;
	virtual std::string getInputSummary(const std::string& prepend_string) const override;

 private:
	double x_left_ = 0.0;
	double k_left_ = 0.0;
	Range  x_left_range_;
};


// Right (one-sided) harmonic potential
class RightHarmonicPotential : public Potential
{
 public:
	RightHarmonicPotential(const ParameterPack& input_pack, const SimulationState& simulation_state);
	virtual void evaluate(const double x, double& u_bias, double& du_bias_dx) const override;
	virtual std::string getInputSummary(const std::string& prepend_string) const override;

 private:
	double x_right_ = 0.0;
	double k_right_ = 0.0;
	Range  x_right_range_;
};

#endif // ifndef POTENTIAL_H
