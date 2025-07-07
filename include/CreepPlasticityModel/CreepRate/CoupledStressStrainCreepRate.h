//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "J2CreepPlasticity.h"

class CoupledStressStrainCreepRate : public J2CreepPlasticity
{
public:
  static InputParameters validParams();

  CoupledStressStrainCreepRate(const InputParameters & parameters);

protected:
  /// Compute creep rate using coupled stress-strain model:
  /// dot_epsilon_c_eq = A1*sigma_eq^n1*epsilon_c_eq^m1 + A2*sigma_eq^n2*epsilon_c_eq^m2 + A3*sigma_eq^n3*epsilon_c_eq^m3
  virtual ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;
  
  /// Compute derivative of creep rate with respect to effective stress
  virtual ADReal computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;
  
  /// Compute derivative of creep rate with respect to effective creep strain
  virtual ADReal computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;

  /// Material constants for the first term: A1 * sigma_eq^n1 * epsilon_c_eq^m1
  const Real _A1;
  const Real _n1;
  const Real _m1;
  
  /// Material constants for the second term: A2 * sigma_eq^n2 * epsilon_c_eq^m2
  const Real _A2;
  const Real _n2;
  const Real _m2;
  
  /// Material constants for the third term: A3 * sigma_eq^n3 * epsilon_c_eq^m3
  const Real _A3;
  const Real _n3;
  const Real _m3;
}; 