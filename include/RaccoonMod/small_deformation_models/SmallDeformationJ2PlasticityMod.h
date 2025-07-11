//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "SmallDeformationPlasticityModelMod.h"

class SmallDeformationJ2PlasticityMod : public SmallDeformationPlasticityModelMod
{
public:
  static InputParameters validParams();

  SmallDeformationJ2PlasticityMod(const InputParameters & parameters);

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) override;

protected:
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ep) override;
};
