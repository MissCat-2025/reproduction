//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "UO2CreepRate.h"

registerMooseObject("reproductionApp", UO2CreepRate);

InputParameters
UO2CreepRate::validParams()
{
  InputParameters params = J2CreepPlasticity::validParams();
  params.addClassDescription("UO2复杂蠕变模型，包含热蠕变、辐照蠕变、转变应力逻辑和时间依赖瞬态蠕变");

  // 物理常数
  params.addParam<Real>("gas_constant", 8.314, "气体常数 (J/mol-K)");
  
  // 温度相关
  params.addRequiredCoupledVar("temperature", "温度 (K)");
  params.addRequiredCoupledVar("oxygen_ratio", "氧化学计量比");
  
  // 裂变与材料特性参数
  params.addRequiredParam<Real>("fission_rate", "裂变率密度 (fissions/m^3-s)");
  params.addParam<Real>("theoretical_density", 95.0, "理论密度百分比");
  params.addParam<Real>("grain_size", 10.0, "晶粒尺寸 (微米)");
  
  // 模型系数
  params.addParam<Real>("a1", 0.3919, "系数a1");
  params.addParam<Real>("a2", 1.31e-19, "系数a2");
  params.addParam<Real>("a3", 87.7, "系数a3");

  params.addParam<Real>("a5", 2.0391e-25, "系数a5");
  params.addParam<Real>("a6", 90.5, "系数a6");

  params.addParam<Real>("a8", 3.7226e-35, "系数a8");

  
  // 激活能
  params.addParam<Real>("Q3", 21759.0, "激活能Q3 (J/mol)");
  
  // 应力相关选项
  params.addParam<bool>("use_transition_stress", false, "是否使用应力转变点");
  
  // 时间依赖瞬态蠕变参数
  params.addParam<bool>("use_transient_creep", false, "是否使用时间依赖瞬态蠕变");
  params.addParam<Real>("transient_decay_constant", 1.40e-6, "瞬态蠕变衰减常数 (1/s)");
  params.addParam<Real>("transient_multiplier", 2.5, "瞬态蠕变倍数");

  return params;
}

UO2CreepRate::UO2CreepRate(const InputParameters & parameters)
  : J2CreepPlasticity(parameters),
    _gas_constant(getParam<Real>("gas_constant")),
    _temperature(adCoupledValue("temperature")),
    _oxygen_ratio(adCoupledValue("oxygen_ratio")),
    _fission_rate(getParam<Real>("fission_rate")),
    _theoretical_density(getParam<Real>("theoretical_density")),
    _grain_size(getParam<Real>("grain_size")),
    _Q1(0.0),  // 将根据温度和氧化学计量比计算
    _Q2(0.0),  // 将根据温度和氧化学计量比计算
    _Q3(getParam<Real>("Q3")),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _a3(getParam<Real>("a3")),
    _a5(getParam<Real>("a5")),
    _a6(getParam<Real>("a6")),
    _a8(getParam<Real>("a8")),
    _use_transition_stress(getParam<bool>("use_transition_stress")),
    _sigma_trans(0.0),  // 将根据晶粒尺寸计算
    _use_transient_creep(getParam<bool>("use_transient_creep")),
    _transient_decay_constant(getParam<Real>("transient_decay_constant")),
    _transient_multiplier(getParam<Real>("transient_multiplier")),
    _max_stress_history(declareProperty<Real>("max_stress_history")),
    _max_stress_history_old(getMaterialPropertyOld<Real>("max_stress_history")),
    _time_since_max_stress(declareProperty<Real>("time_since_max_stress")),
    _time_since_max_stress_old(getMaterialPropertyOld<Real>("time_since_max_stress")),
    _density_term1(0.0),  // 将在每个时间步计算
    _density_term2(0.0),  // 将在每个时间步计算
    _fission_term(0.0),   // 将在每个时间步计算
    _exp_Q1(0.0),         // 将在每个时间步计算
    _exp_Q2(0.0),         // 将在每个时间步计算
    _exp_Q3(0.0)          // 将在每个时间步计算
{
}

void
UO2CreepRate::setQp(unsigned int qp)
{
  J2CreepPlasticity::setQp(qp);
  
  // 预计算常用值
  const ADReal T = _temperature[_qp];
  const ADReal RT = _gas_constant * T;
  const ADReal inv_RT = 1.0 / RT;
  
  // 使用氧化学计量比计算激活能
  const ADReal x = _oxygen_ratio[_qp];
  const ADReal log_x = std::log10(x);
  const ADReal exp_common = std::exp(-20.0 / log_x - 8.0);
  const ADReal denom = 1.0 / (exp_common + 1.0);
  
  // 计算激活能
  _Q1 = 74829.0 * denom + 301762.0;
  _Q2 = 83143.0 * denom + 469191.0;
  
  // 预计算指数项
  _exp_Q1 = std::exp(-_Q1 * inv_RT);
  _exp_Q2 = std::exp(-_Q2 * inv_RT);
  _exp_Q3 = std::exp(-_Q3 * inv_RT);
  
  // 计算转变应力
  _sigma_trans = 1.6547e7 * std::pow(_grain_size, 0.5714);
  
  // 预计算密度和晶粒尺寸相关项
  _density_term1 = 1.0 / ((_theoretical_density - _a3) * _grain_size * _grain_size);
  _density_term2 = 1.0 / (_theoretical_density - _a6);
  
  // 预计算裂变率项
  _fission_term = _a1 + _a2 * _fission_rate;
}

ADReal
UO2CreepRate::computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 计算稳态蠕变率
  ADReal steady_state_rate = computeSteadyStateCreepRate(effective_stress, effective_creep_strain);
  
  // 如果不使用瞬态蠕变，直接返回稳态蠕变率
  if (!_use_transient_creep)
    return steady_state_rate;
  
  // 更新最大应力历史
  updateMaxStressHistory(effective_stress);
  
  // 计算瞬态时间因子
  ADReal transient_factor = computeTransientTimeFactor();
  
  // 返回考虑瞬态效应的总蠕变率
  return steady_state_rate * transient_factor;
}

ADReal
UO2CreepRate::computeSteadyStateCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 计算各分量蠕变率
  ADReal creep_th1 = 0.0;  // 热蠕变线性项
  ADReal creep_th2 = 0.0;  // 热蠕变幂律项
  ADReal creep_ir = 0.0;   // 辐照蠕变
  
  // 根据是否使用转变应力来计算热蠕变
  if (_use_transition_stress)
  {
    // 使用转变应力逻辑
    if (effective_stress < _sigma_trans)
    {
      // 低应力区域：使用线性项，忽略幂律项
      creep_th1 = _fission_term * _density_term1 * effective_stress * _exp_Q1;
      creep_th2 = 0.0;
    }
    else
    {
      // 高应力区域：线性项使用转变应力，幂律项使用实际应力
      creep_th1 = _fission_term * _density_term1 * _sigma_trans * _exp_Q1;
      creep_th2 = (_a5) * _density_term2 * std::pow(effective_stress, 4.5) * _exp_Q2;
    }
  }
  else
  {
    // 不使用转变应力 - 同时应用两项
    creep_th1 = _fission_term * _density_term1 * effective_stress * _exp_Q1;
    creep_th2 = (_a5) * _density_term2 * std::pow(effective_stress, 4.5) * _exp_Q2;
  }
  
  // 辐照蠕变
  creep_ir = _a8 * _fission_rate * effective_stress * _exp_Q3;
  
  // 总稳态蠕变率
  return creep_th1 + creep_th2 + creep_ir;
}

ADReal
UO2CreepRate::computeTransientTimeFactor()
{
  // 获取自最大应力施加以来的时间
  Real time_since_max = _time_since_max_stress[_qp];
  
  // 计算瞬态时间因子：2.5*exp(-1.40e-6*t) + 1
  Real transient_factor = _transient_multiplier * std::exp(-_transient_decay_constant * time_since_max) + 1.0;
  
  return transient_factor;
}

void
UO2CreepRate::updateMaxStressHistory(const ADReal & effective_stress)
{
  // 获取当前应力的数值
  Real current_stress = MetaPhysicL::raw_value(effective_stress);
  Real max_stress_old = _max_stress_history_old[_qp];
  Real time_since_max_old = _time_since_max_stress_old[_qp];
  
  // 检查是否达到新的最大应力
  if (current_stress > max_stress_old)
  {
    // 更新最大应力历史
    _max_stress_history[_qp] = current_stress;
    // 重置时间计数器
    _time_since_max_stress[_qp] = 0.0;
  }
  else
  {
    // 保持旧的最大应力
    _max_stress_history[_qp] = max_stress_old;
    // 增加时间计数器
    _time_since_max_stress[_qp] = time_since_max_old + _dt;
  }
}

ADReal
UO2CreepRate::computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 计算稳态蠕变率对应力的导数
  ADReal d_creep_th1_d_stress = 0.0;
  ADReal d_creep_th2_d_stress = 0.0;
  ADReal d_creep_ir_d_stress = 0.0;
  
  // 根据是否使用转变应力计算导数
  if (_use_transition_stress)
  {
    if (effective_stress < _sigma_trans)
    {
      // 低应力区域：只有线性项有导数
      d_creep_th1_d_stress = _fission_term * _density_term1 * _exp_Q1;
      d_creep_th2_d_stress = 0.0;
    }
    else
    {
      // 高应力区域：线性项无应力依赖（因为使用转变应力），幂律项有导数
      d_creep_th1_d_stress = 0.0;
      d_creep_th2_d_stress = 4.5 * (_a5) * _density_term2 * std::pow(effective_stress, 3.5) * _exp_Q2;
    }
  }
  else
  {
    // 不使用转变应力：两项都有导数
    d_creep_th1_d_stress = _fission_term * _density_term1 * _exp_Q1;
    d_creep_th2_d_stress = 4.5 * (_a5) * _density_term2 * std::pow(effective_stress, 3.5) * _exp_Q2;
  }
  
  // 辐照蠕变导数
  d_creep_ir_d_stress = _a8 * _fission_rate * _exp_Q3;
  
  // 总稳态蠕变率导数
  ADReal d_steady_state_d_stress = d_creep_th1_d_stress + d_creep_th2_d_stress + d_creep_ir_d_stress;
  
  // 如果不使用瞬态蠕变，直接返回稳态导数
  if (!_use_transient_creep)
    return d_steady_state_d_stress;
  
  // 考虑瞬态时间因子
  ADReal transient_factor = computeTransientTimeFactor();
  
  return d_steady_state_d_stress * transient_factor;
}

ADReal
UO2CreepRate::computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // UO2蠕变率不依赖于累积蠕变应变，所以导数为零
  return 0.0;
} 