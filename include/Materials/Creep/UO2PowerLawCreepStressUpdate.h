//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADRadialReturnCreepStressUpdateBase.h"

/**
 * 这个类实现了UO2材料的复杂蠕变模型，并考虑相场断裂中的退化函数
 * 蠕变率表达式为：
 * \dot{\varepsilon}^{\mathrm{cr}}=\left[\frac{a_1+a_2\dot{f}}{(a_3+D)g^2}\exp\frac{-Q_1}{RT}+a_7\dot{f}\exp\frac{Q_3}{RT}\right]\sigma+\left(\frac{a_1+a_8\dot{f}}{a_6+D}\exp\frac{Q_2}{RT}\right)\sigma^{4.5}
 */
class UO2PowerLawCreepStressUpdate : public ADRadialReturnCreepStressUpdateBase
{
public:
  static InputParameters validParams();

  UO2PowerLawCreepStressUpdate(const InputParameters & parameters);

  virtual Real computeStrainEnergyRateDensity(
      const ADMaterialProperty<RankTwoTensor> & stress,
      const ADMaterialProperty<RankTwoTensor> & strain_rate) override;

protected:
  virtual void computeStressInitialize(const ADReal & effective_trial_stress,
                                     const ADRankFourTensor & elasticity_tensor) override;

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                               const ADReal & scalar) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override;

  /// 我们还需要在 updateState 里把 trial stress 的偏应力存下来
  virtual void updateState(
      GenericRankTwoTensor<true> & strain_increment,
      GenericRankTwoTensor<true> & inelastic_strain_increment,
      const GenericRankTwoTensor<true> & rotation_increment,
      GenericRankTwoTensor<true> & stress_new,
      const RankTwoTensor & stress_old,
      const GenericRankFourTensor<true> & elasticity_tensor,
      const RankTwoTensor & elastic_strain_old,
      bool compute_full_tangent_operator,
      RankFourTensor & tangent_operator) override;

  ADReal compute_three_shear_modulus_New(const GenericRankFourTensor<true> & elasticity_tensor);
  // 基本物理常数
  const Real _gas_constant;
  
  // 温度相关参数
  const ADVariableValue & _temperature;
  const ADVariableValue & _oxygen_ratio;
  
  // 裂变率和材料特性参数
  const Real _fission_rate;
  const Real _theoretical_density;
  const Real _grain_size;
  
  // 激活能
  ADReal _Q1;
  ADReal _Q2;
  ADReal _Q3;
  
  // 模型常数
  const Real _a1, _a2, _a3, _a6, _a7, _a8;
  
  // 临时变量，用于计算
  ADReal _exp_Q1;
  ADReal _exp_Q2;
  ADReal _exp_Q3;
  // 创蠕变变形
  ADMaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;
  
  // 退化相关参数
  const ADMaterialProperty<Real> & _g;

  // 预计算的常数项
  ADReal _density_term1;   // 晶粒尺寸相关项1
  ADReal _density_term2;   // 晶粒尺寸相关项2
  ADReal _fission_term;    // 裂变率项

    // ——— 新增 —— 用于存储弹性偏应力，后面做谱分解
  ADMaterialProperty<RankTwoTensor> & _deviatoric_trial_stress;
};