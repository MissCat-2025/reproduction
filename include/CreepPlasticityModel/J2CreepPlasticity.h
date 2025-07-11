//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "CreepModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class PlasticityModel;  // 前向声明
class PlasticHardeningModel;  // 前向声明

class J2CreepPlasticity : public CreepModel,
                          public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  J2CreepPlasticity(const InputParameters & parameters);

  virtual void initialSetup() override;

  virtual void setQp(unsigned int qp) override;

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain) override;

protected:

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ec) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ec) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ec) override;
  //判断非塑性应变下是否满足区分条件（仅在有塑性模型时使用）
  bool f_no_plastic_strain(const ADReal & effective_stress);
  
  /// 求解塑性-蠕变耦合问题（仅在有塑性模型时使用）
  void solveCreepPlasticityCoupled(ADRankTwoTensor & stress, 
                                   ADRankTwoTensor & elastic_strain,
                                   const ADRankTwoTensor & elastic_trial_strain,
                                   const ADReal & effective_stress_trial,
                                   const ADReal & delta_ec_np,
                                   const ADRankTwoTensor & creep_strain_np,
                                   const ADReal & ec_np);
  
  /// 手动更新塑性模型状态变量（仅在有塑性模型时使用）
  void updatePlasticityModelState(const ADReal & delta_ep);
  
  /// 计算弹性修正系数C（根据应力计算方法选择）
  ADReal computeElasticModifier();
  
  /// 统一的应力计算函数（根据策略选择使用不同的方法）
  ADRankTwoTensor computeStressUnified(const ADRankTwoTensor & elastic_strain,bool loop);
  
  // 蠕变率计算的虚函数接口，由具体的蠕变模型实现
  virtual ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  virtual ADReal computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);
  virtual ADReal computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain);

  /// The plasticity model for coupled creep-plasticity computation (可选)
  PlasticityModel * _plasticity_model;
  
  /// The plastic hardening model (仅在有塑性模型时使用)
  PlasticHardeningModel * _hardening_model;

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
  
  /// 是否使用全程3G策略（仅在_use_three_shear_modulus为true时生效）
  const bool _full_three_shear_modulus_strategy;
  
  /// 3倍剪切模量（仅在_use_three_shear_modulus为true时使用）
  ADReal _three_shear_modulus;

  /// 降解函数及其导数
  const ADMaterialProperty<Real> & _gc;
  const ADMaterialProperty<Real> & _dgc_dd;
};