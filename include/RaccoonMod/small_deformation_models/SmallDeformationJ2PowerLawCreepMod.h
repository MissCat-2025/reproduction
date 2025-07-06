//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "SmallDeformationJ2PlasticityMod.h"

class SmallDeformationJ2PowerLawCreepMod : public SmallDeformationJ2PlasticityMod
{
public:
  static InputParameters validParams();

  SmallDeformationJ2PowerLawCreepMod(const InputParameters & parameters);

protected:
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  // @{ The parameters of the power law creep update
  const Real _coefficient;
  const Real _exponent;
  // @}
}; 