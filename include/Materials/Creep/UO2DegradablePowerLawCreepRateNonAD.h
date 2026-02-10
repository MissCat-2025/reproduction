#pragma once

#include "Material.h"
#include "RankTwoTensor.h"

class UO2DegradablePowerLawCreepRateNonAD : public Material
{
public:
  static InputParameters validParams();
  UO2DegradablePowerLawCreepRateNonAD(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const Real _gas_constant;
  const Real _fission_rate;
  const Real _theoretical_density;
  const Real _grain_size;

  const Real _Q3;

  const Real _a1;
  const Real _a2;
  const Real _a3;
  const Real _a6;
  const Real _a7;
  const Real _a8;

  const bool _use_transition_stress;

  const VariableValue & _temperature;
  const VariableValue * _oxygen_ratio;
  const Real _oxygen_ratio_value;

  const MaterialProperty<RankTwoTensor> & _stress_old;

  const MaterialProperty<Real> & _g;
  const bool _use_stress_degradation;

  MaterialProperty<RankTwoTensor> & _creep_rate;
  MaterialProperty<Real> & _effective_creep;
};

