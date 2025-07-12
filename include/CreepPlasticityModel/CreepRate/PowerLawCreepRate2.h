//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "J2CreepPlasticity.h"

class PowerLawCreepRate2 : public J2CreepPlasticity
{
public:
  static InputParameters validParams();

  PowerLawCreepRate2(const InputParameters & parameters);

protected:
  /// Compute creep rate using power law: dot_epsilon_c_eq = A * sigma_eq^n
  virtual ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;
  
  /// Compute derivative of creep rate with respect to effective stress
  virtual ADReal computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;
  
  /// Compute derivative of creep rate with respect to effective creep strain (zero for power law)
  virtual ADReal computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;

  /// Material constant A in power law creep model
  const Real _A;
  
  /// Material constant n in power law creep model
  const Real _n;

  /// Material constant Q in power law creep model
  const Real _Q;

  /// Material constant R in power law creep model
  const Real _R;

  /// Temperature in power law creep model
  const ADVariableValue & _T;

}; 