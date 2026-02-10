#include "EqualValueBoundaryConstraintDegraded.h"

#include <algorithm>
#include <cmath>

registerMooseObject("reproductionApp", EqualValueBoundaryConstraintDegraded);

InputParameters
EqualValueBoundaryConstraintDegraded::validParams()
{
  InputParameters params = EqualValueBoundaryConstraint::validParams();
  params.addClassDescription("EqualValueBoundaryConstraint with degradation scaling.");
  params.addCoupledVar("degradation_variable", "Optional variable to degrade penalty");
  params.addParam<Real>("degradation_power", 1.0, "Power for degradation scaling");
  params.addParam<Real>("degradation_a2", 0.0, "Rational degradation parameter a2");
  params.addParam<Real>("degradation_a3", 0.0, "Rational degradation parameter a3");
  params.addParam<Real>("degradation_eta", 0.0, "Residual penalty floor for degradation");
  params.addParam<bool>("use_rational_degradation",
                        false,
                        "Use rational degradation function based on sigma0 field");
  params.addCoupledVar("sigma0_variable", "Optional sigma0 field for rational degradation");
  params.addParam<Real>("youngs_modulus", 0.0, "Young's modulus used to compute a1");
  params.addParam<Real>("fracture_energy", 0.0, "Fracture energy used to compute a1");
  params.addParam<Real>("length_scale", 0.0, "Length scale used to compute a1");
  return params;
}

EqualValueBoundaryConstraintDegraded::EqualValueBoundaryConstraintDegraded(
    const InputParameters & parameters)
  : EqualValueBoundaryConstraint(parameters),
    _degradation_primary(nullptr),
    _degradation_secondary(nullptr),
    _use_degradation(isCoupled("degradation_variable")),
    _degradation_power(getParam<Real>("degradation_power")),
    _degradation_a2(getParam<Real>("degradation_a2")),
    _degradation_a3(getParam<Real>("degradation_a3")),
    _degradation_eta(getParam<Real>("degradation_eta")),
    _use_rational_degradation(getParam<bool>("use_rational_degradation")),
    _sigma0_primary(nullptr),
    _sigma0_secondary(nullptr),
    _youngs_modulus(getParam<Real>("youngs_modulus")),
    _fracture_energy(getParam<Real>("fracture_energy")),
    _length_scale(getParam<Real>("length_scale"))
{
  if (_use_degradation)
  {
    _degradation_primary = &coupledValue("degradation_variable");
    _degradation_secondary = &coupledNeighborValue("degradation_variable");
  }
  if (_use_rational_degradation)
  {
    if (!_use_degradation)
      paramError("degradation_variable", "Rational degradation requires degradation_variable.");
    if (!isCoupled("sigma0_variable"))
      paramError("sigma0_variable", "Rational degradation requires sigma0_variable.");
    if (_youngs_modulus <= 0.0 || _fracture_energy <= 0.0 || _length_scale <= 0.0)
      paramError("youngs_modulus",
                 "youngs_modulus, fracture_energy, and length_scale must be positive.");
    _sigma0_primary = &coupledValue("sigma0_variable");
    _sigma0_secondary = &coupledNeighborValue("sigma0_variable");
  }
}

Real
EqualValueBoundaryConstraintDegraded::computeQpResidual(Moose::ConstraintType type)
{
  switch (type)
  {
    case Moose::Secondary:
      return (_u_secondary[_i] - _u_primary[_j]) * _penalty * degradationFactorSecondary();
    case Moose::Primary:
      return (_u_primary[_j] - _u_secondary[_i]) * _penalty * degradationFactorPrimary();
  }
  return 0.;
}

Real
EqualValueBoundaryConstraintDegraded::computeQpJacobian(Moose::ConstraintJacobianType type)
{
  switch (type)
  {
    case Moose::SecondarySecondary:
      return _penalty * degradationFactorSecondary();
    case Moose::SecondaryPrimary:
      return -_penalty * degradationFactorSecondary();
    case Moose::PrimaryPrimary:
      return _penalty * degradationFactorPrimary();
    case Moose::PrimarySecondary:
      return -_penalty * degradationFactorPrimary();
    default:
      mooseError("Unsupported type");
      break;
  }
  return 0.;
}

Real
EqualValueBoundaryConstraintDegraded::degradationFactorPrimary() const
{
  if (!_use_degradation)
    return 1.0;
  Real d = (*_degradation_primary)[_j];
  d = std::min(1.0, std::max(0.0, d));
  if (_use_rational_degradation)
  {
    const Real sigma0 = (*_sigma0_primary)[_j];
    const Real a1 = 4.0 * _youngs_modulus * _fracture_energy /
                    (sigma0 * sigma0 * _length_scale * 3.14159);
    const Real one_minus_d = 1.0 - d;
    const Real one_minus_d_p = std::pow(one_minus_d, _degradation_power);
    const Real denom = one_minus_d_p +
                       a1 * d * (1.0 + _degradation_a2 * d + _degradation_a3 * d * d);
    const Real g = one_minus_d_p / denom;
    return g * (1.0 - _degradation_eta) + _degradation_eta;
  }
  return std::pow(1.0 - d, _degradation_power) * (1.0 - _degradation_eta) + _degradation_eta;
}

Real
EqualValueBoundaryConstraintDegraded::degradationFactorSecondary() const
{
  if (!_use_degradation)
    return 1.0;
  Real d = (*_degradation_secondary)[_i];
  d = std::min(1.0, std::max(0.0, d));
  if (_use_rational_degradation)
  {
    const Real sigma0 = (*_sigma0_secondary)[_i];
    const Real a1 = 4.0 * _youngs_modulus * _fracture_energy /
                    (sigma0 * sigma0 * _length_scale * 3.14159);
    const Real one_minus_d = 1.0 - d;
    const Real one_minus_d_p = std::pow(one_minus_d, _degradation_power);
    const Real denom = one_minus_d_p +
                       a1 * d * (1.0 + _degradation_a2 * d + _degradation_a3 * d * d);
    const Real g = one_minus_d_p / denom;
    return g * (1.0 - _degradation_eta) + _degradation_eta;
  }
  return std::pow(1.0 - d, _degradation_power) * (1.0 - _degradation_eta) + _degradation_eta;
}
