#pragma once

#include "ComputeMultipleInelasticStress.h"

class ElasticityModelNonAD;

class ComputeCreepPlasticityDeformationStressNew : public ComputeMultipleInelasticStress
{
public:
  static InputParameters validParams();

  ComputeCreepPlasticityDeformationStressNew(const InputParameters & parameters);

  void initialSetup() override;

protected:
  void computeQpStress() override;

private:
  const bool _has_elasticity_model;
  ElasticityModelNonAD * _elasticity_model;
  const bool _use_elasticity_model_for_stress;

  const bool _apply_damage_degradation;
  const MaterialPropertyName _g_name;
  const MaterialProperty<Real> * _g;
  const Real _min_degradation;

  const bool _debug;
  const int _debug_qp;
  const unsigned int _debug_step_interval;
};
