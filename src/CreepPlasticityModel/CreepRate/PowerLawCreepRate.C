//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PowerLawCreepRate.h"

registerMooseObject("reproductionApp", PowerLawCreepRate);

InputParameters
PowerLawCreepRate::validParams()
{
  InputParameters params = J2CreepPlasticity::validParams();
  params.addClassDescription("Power law creep rate model: dot_epsilon_c_eq = A * sigma_eq^n");
  
  params.addRequiredParam<Real>("A", "Material constant A in power law creep model");
  params.addRequiredParam<Real>("n", "Material constant n in power law creep model");
  
  return params;
}

PowerLawCreepRate::PowerLawCreepRate(const InputParameters & parameters)
  : J2CreepPlasticity(parameters),
    _A(getParam<Real>("A")),
    _n(getParam<Real>("n"))
{
}

ADReal
PowerLawCreepRate::computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // Power law creep rate: dot_epsilon_c_eq = A * sigma_eq^n
  return _A * std::pow(effective_stress, _n);
}

ADReal
PowerLawCreepRate::computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // Derivative of creep rate with respect to effective stress
  // d(dot_epsilon_c_eq)/d(sigma_eq) = A * n * sigma_eq^(n-1)
  if (effective_stress <= 0.0)
    return 0.0;
  
  return _A * _n * std::pow(effective_stress, _n - 1.0);
}

ADReal
PowerLawCreepRate::computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // Power law creep rate does not depend on accumulated creep strain
  // so the derivative is zero
  return 0.0;
} 