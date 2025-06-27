// include/materials/UO2CreepEigenstrain.h
#pragma once
#include "ADComputeEigenstrainBase.h"

class UO2CreepEigenstrain : public ADComputeEigenstrainBase
{
public:
  static InputParameters validParams();
  UO2CreepEigenstrain(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpEigenstrain() override;

  /// 从UO2CreepRate获取蠕变率
  const ADMaterialProperty<RankTwoTensor> & _creep_rate;

  /// 累积的蠕变特征应变 - 旧值不需要AD
  ADMaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;

  // 添加psip_active和相关材料属性
  bool _consider_psip_active;
  ADMaterialProperty<Real> & _psip_active;
  const MaterialProperty<Real> & _psip_active_old;

  // 添加psip_active和相关材料属性
  const ADMaterialProperty<RankTwoTensor> & _stress_deviator;
};