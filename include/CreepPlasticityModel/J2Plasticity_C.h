//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "PlasticityModel.h"

class J2Plasticity_C : public PlasticityModel
{
public:
  static InputParameters validParams();

  J2Plasticity_C(const InputParameters & parameters);

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) override;

protected:
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ep) override;
};
