//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "CreepModel.h"
#include "PlasticHardeningModel.h"

class J2Creep_P : public CreepModel
{
public:
  static InputParameters validParams();

  J2Creep_P(const InputParameters & parameters);

  virtual void initialSetup() override;

  virtual void setQp(unsigned int qp) override;

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) override;

  /// Public interface for accessing creep rate calculation (needed by plasticity model)
  virtual ADReal getCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  
  /// Public interface for accessing creep rate stress derivative (needed by plasticity model)
  virtual ADReal getCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  
  /// Public interface for accessing current effective creep strain (needed by plasticity model)
  virtual ADReal getEffectiveCreepStrain() const;

protected:

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ec) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ec) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ec) override;
  //判断非塑性应变下是否满足区分条件
  bool f_no_plastic_strain(const ADReal & effective_stress);
  
  // 蠕变率计算的虚函数接口，由具体的蠕变模型实现
  virtual ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  virtual ADReal computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  virtual ADReal computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);

  /// The plastic hardening model for yield stress calculation
  PlasticHardeningModel * _hardening_model;
};
