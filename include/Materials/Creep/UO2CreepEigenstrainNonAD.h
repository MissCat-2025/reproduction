#pragma once

#include "ComputeEigenstrainBase.h"
#include "RankTwoTensor.h"

class UO2CreepEigenstrainNonAD : public ComputeEigenstrainBase
{
public:
  static InputParameters validParams();
  UO2CreepEigenstrainNonAD(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpEigenstrain() override;

  const MaterialProperty<RankTwoTensor> & _creep_rate;
  const MaterialProperty<Real> & _effective_creep;

  MaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;

  MaterialProperty<Real> & _effective_creep_strain;
  const MaterialProperty<Real> & _effective_creep_strain_old;
};

