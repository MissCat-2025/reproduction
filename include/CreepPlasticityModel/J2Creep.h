//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "CreepModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class J2Creep : public CreepModel,
                          public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  J2Creep(const InputParameters & parameters);

  virtual void setQp(unsigned int qp) override;

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) override;

protected:

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ec) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ec) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ec) override;

  /// 统一的应力计算函数（根据策略选择使用不同的方法）
  ADRankTwoTensor computeStressUnified(const ADRankTwoTensor & elastic_strain);
  
  // 蠕变率计算的虚函数接口，由具体的蠕变模型实现
  virtual ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  virtual ADReal computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  virtual ADReal computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);

  /// 蠕变输出属性
  /// 有效蠕变应变 (标量)
  ADMaterialProperty<Real> & _effective_creep_strain;
    /// 相场变量名
  const std::string _d_name;
    /// 蠕变能量密度属性名
  const std::string _psic_name;
    /// 总蠕变能量密度 (考虑相场降解)
  ADMaterialProperty<Real> & _psic;
  /// 蠕变能量密度的活跃部分
  ADMaterialProperty<Real> & _psic_active;
  /// 蠕变能量密度的活跃部分的旧值
  const MaterialProperty<Real> & _psic_active_old;
  
  /// 蠕变能量密度对相场的导数
  ADMaterialProperty<Real> & _dpsic_dd;

  /// 应力计算方法选择
  const bool _use_three_shear_modulus;
  
  /// 3倍剪切模量（仅在_use_three_shear_modulus为true时使用）
  ADReal _three_shear_modulus;

  /// 降解函数及其导数
  const ADMaterialProperty<Real> & _gc;
  const ADMaterialProperty<Real> & _dgc_dd;
};