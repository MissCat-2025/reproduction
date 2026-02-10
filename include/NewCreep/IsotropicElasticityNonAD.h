#pragma once

#include "ElasticityModelNonAD.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class IsotropicElasticityNonAD : public ElasticityModelNonAD, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  IsotropicElasticityNonAD(const InputParameters & parameters);

  virtual RankTwoTensor computeStress(const RankTwoTensor & elastic_strain) override;
  virtual RankFourTensor computeJacobianMult(const RankTwoTensor & elastic_strain) override;
  void initQpStatefulProperties() override;

protected:
  RankTwoTensor computeStressNone(const RankTwoTensor & strain);
  RankTwoTensor computeStressSpectralDecomposition(const RankTwoTensor & strain);
  RankTwoTensor computeStressVolDevDecomposition(const RankTwoTensor & strain);
  RankTwoTensor computeStressMaxPrincipalDecomposition(const RankTwoTensor & strain);

  void updateHistoryMax(const Real Y_bar);
  Real applyThreshold(const Real psie_active_raw) const;

  const Real _youngs_modulus;
  const Real _poissons_ratio;

  const VariableValue & _phase_field;
  const VariableName _d_name;

  const MaterialPropertyName _g_name;
  const MaterialProperty<Real> & _g;
  const MaterialProperty<Real> & _dg_dd;

  MaterialProperty<Real> & _psie_active;
  MaterialProperty<Real> & _psie;
  MaterialProperty<Real> & _dpsie_dd;

  const MooseEnum _kinematic_assumption;
  const MooseEnum _decomposition;
  const bool _use_threshold;
  const bool _use_history_max;

  const MaterialProperty<Real> & _tensile_strength;

  MaterialProperty<Real> * _history_max;
  const MaterialProperty<Real> * _history_max_old;

  const bool _debug;
  const int _debug_qp;
  const unsigned int _debug_step_interval;

  const bool _degrade_out_of_plane_strain;
};
