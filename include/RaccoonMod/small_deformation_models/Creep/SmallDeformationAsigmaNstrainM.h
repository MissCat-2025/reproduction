//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "SmallDeformationJ2PlasticityMod.h"

class SmallDeformationAsigmaNstrainM : public SmallDeformationJ2PlasticityMod
{
  // 在头文件中添加
// private:
//   const Real _stress_unit;
//   const Real _initial_creep_strain;
public:
  static InputParameters validParams();

  SmallDeformationAsigmaNstrainM(const InputParameters & parameters);

protected:
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  // @{ The parameters for the three-term creep model: A * sigma^n * strain^m
  const Real _A1, _n1, _m1;
  const Real _A2, _n2, _m2;
  const Real _A3, _n3, _m3;
  // @}
};