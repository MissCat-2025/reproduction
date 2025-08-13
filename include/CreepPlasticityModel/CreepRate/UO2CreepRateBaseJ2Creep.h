//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "J2Creep.h"

/**
 * UO2材料复杂蠕变模型，包含：
 * 1. 热蠕变（线性项和幂律项）
 * 2. 辐照蠕变
 * 3. 转变应力逻辑
 * 4. 时间依赖瞬态蠕变
 * 蠕变率表达式：
 * \dot{\varepsilon}^{\mathrm{cr}}=\left[\frac{a_1+a_2\dot{f}}{(a_3+D)g^2}\exp\frac{-Q_1}{RT}+a_7\dot{f}\exp\frac{-Q_3}{RT}\right]\sigma+\left(\frac{a_1+a_8\dot{f}}{a_6+D}\exp\frac{-Q_2}{RT}\right)\sigma^{4.5}
 * 时间依赖项：\dot{\varepsilon}_T = \dot{\varepsilon}_S[2.5e^{(-1.40\times10^{-6}t)}+1]
 */
class UO2CreepRateBaseJ2Creep : public J2Creep
{
public:
  static InputParameters validParams();

  UO2CreepRateBaseJ2Creep(const InputParameters & parameters);

  virtual void setQp(unsigned int qp) override;

protected:
  /// Compute creep rate using complex UO2 creep model
  virtual ADReal computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;
  
  /// Compute derivative of creep rate with respect to effective stress
  virtual ADReal computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;
  
  /// Compute derivative of creep rate with respect to effective creep strain
  virtual ADReal computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain) override;

  /// 计算稳态蠕变率（不含时间依赖项）
  ADReal computeSteadyStateCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain);

  /// 计算瞬态时间因子
  ADReal computeTransientTimeFactor();

  /// 更新最大应力历史
  void updateMaxStressHistory(const ADReal & effective_stress);

  // 基本物理常数
  const Real _gas_constant;
  
  // 温度相关参数
  const ADVariableValue & _temperature;
  const ADVariableValue & _oxygen_ratio;
  
  // 裂变率和材料特性参数
  const Real _fission_rate;
  const Real _theoretical_density;
  const Real _grain_size;
  
  // 激活能（温度依赖）
  mutable ADReal _Q1;
  mutable ADReal _Q2;
  const Real _Q3;
  
  // 模型常数
  const Real _a1, _a2, _a3, _a5, _a6, _a8;
  
  // 转变应力相关参数
  const bool _use_transition_stress;
  mutable ADReal _sigma_trans;
  
  // 时间依赖瞬态蠕变参数
  const bool _use_transient_creep;
  const Real _transient_decay_constant;  // 1.40e-6
  const Real _transient_multiplier;      // 2.5
  
  // 最大应力历史跟踪
  MaterialProperty<Real> & _max_stress_history;
  const MaterialProperty<Real> & _max_stress_history_old;
  MaterialProperty<Real> & _time_since_max_stress;
  const MaterialProperty<Real> & _time_since_max_stress_old;
  
  // 预计算的常数项
  mutable ADReal _density_term1;   // 晶粒尺寸相关项1
  mutable ADReal _density_term2;   // 晶粒尺寸相关项2
  mutable ADReal _fission_term;    // 裂变率项
  mutable ADReal _exp_Q1;          // 预计算指数项1
  mutable ADReal _exp_Q2;          // 预计算指数项2
  mutable ADReal _exp_Q3;          // 预计算指数项3
}; 