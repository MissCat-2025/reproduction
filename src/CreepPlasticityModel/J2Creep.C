//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "J2Creep.h"
#include "ElasticityModel.h"
#include "IsotropicElasticity.h"
#include "libmesh/utility.h"
registerMooseObject("reproductionApp", J2Creep);

InputParameters
J2Creep::validParams()
{
  InputParameters params = CreepModel::validParams();
  params.addClassDescription("J2 creep model. Implements the algorithm "
                             "for separated creep with radial return mapping.");
  //塑性 (可选)
  // 牛顿迭代参数
  params.addParam<int>("max_iter", 1000, "Maximum number of iterations for Newton-Raphson solve");
  params.addParam<Real>("tolerance", 1e-10, "Tolerance for Newton-Raphson convergence");
  
  // 应力计算方法选择
  params.addParam<bool>("use_three_shear_modulus", true, "使用3倍剪切模量计算应力修正（true）还是使用完整弹性模型（false）");
  // 蠕变输出属性参数
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "creep_energy_density",
      "psic",
      "Name of the creep energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "gc", "The degradation function for creep");
  
  return params;
}

J2Creep::J2Creep(const InputParameters & parameters) : 
  CreepModel(parameters),
  DerivativeMaterialPropertyNameInterface(),
  // 初始化蠕变输出属性
  _effective_creep_strain(declareADProperty<Real>("effective_creep_strain")),
  _d_name(getVar("phase_field", 0)->name()),
    // 修复参数获取 - 使用正确的参数名
  _psic_name(prependBaseName("creep_energy_density", true)),
  _psic(declareADProperty<Real>(_psic_name)),
  _psic_active(declareADProperty<Real>(_psic_name + "_active")),
  _psic_active_old(getMaterialPropertyOld<Real>(_psic_name + "_active")),
  _dpsic_dd(declareADProperty<Real>(derivativePropertyName(_psic_name, {_d_name}))),
  
  // 应力计算方法选择
  _use_three_shear_modulus(getParam<bool>("use_three_shear_modulus")),
  _three_shear_modulus(0.0),
  
  // 获取降解函数及其导数
  _gc(getADMaterialProperty<Real>(prependBaseName("degradation_function", true))),
  _dgc_dd(getADMaterialProperty<Real>(derivativePropertyName(prependBaseName("degradation_function", true), {_d_name})))
{
}

void
J2Creep::setQp(unsigned int qp)
{
  // 调用基类的setQp方法
  CreepModel::setQp(qp);
}

void
J2Creep::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain)
{
 
  // Step 1: 计算弹性试验状态（移除旧的蠕变应变）
  ADRankTwoTensor elastic_trial_strain = elastic_strain;
  elastic_trial_strain -= _creep_strain_old[_qp];  // 减去旧的蠕变应变张量
  // Step 2: 计算试验状态
  ADReal delta_ec = 0.0;
  // 计算试验应力,直接使用
  ADRankTwoTensor stress_trial = computeStressUnified(elastic_trial_strain);

  // 计算偏应力和有效应力
  ADRankTwoTensor stress_dev = stress_trial.deviatoric();
  ADReal effective_stress_trial = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 计算流动方向
  if (effective_stress_trial > 1e-15 || effective_stress_trial < -1e-15)
    _Nc[_qp] = 1.5 * stress_dev / effective_stress_trial;
  else
    _Nc[_qp].zero();

  //流动方向
  if (_use_three_shear_modulus)
  {
    // 使用3倍剪切模量：3G
    IsotropicElasticity * iso_elasticity = dynamic_cast<IsotropicElasticity*>(_elasticity_model);
    if (iso_elasticity)
      _three_shear_modulus = iso_elasticity->computeThreeShearModulus();
    else
      mooseError("computeThreeShearModulus() is only available for IsotropicElasticity models");
  }
  else
  {
    _three_shear_modulus =  _elasticity_model->computeStress(_Nc[_qp]).doubleContraction(_Nc[_qp]);
  }


  // 求解非塑性蠕变增量（使用径向返回映射）
  if (effective_stress_trial > 1e-15 || effective_stress_trial < -1e-15)
  {
    returnMappingSolve(effective_stress_trial, delta_ec, _console);
  }
  
  // 更新蠕变应变和有效蠕变应变
  _creep_strain[_qp] = _creep_strain_old[_qp] + delta_ec * _Nc[_qp];
  _ec[_qp] = _ec_old[_qp] + delta_ec;
  
  // 计算弹性应变和应力
  elastic_strain = elastic_trial_strain - delta_ec * _Nc[_qp];
  stress = computeStressUnified(elastic_strain);
  // 输出有效蠕变应变（标量）
  _effective_creep_strain[_qp] = _ec[_qp];
  
  // 计算蠕变能量密度：基于有效应力与有效蠕变应变的积分
  // 对于J2蠕变，蠕变能密度 W_c = ∫ σ_eq * dε_eq^c
  // 在增量形式下，使用梯形积分规则：W_c = W_c_old + 0.5 * (σ_eq_old + σ_eq) * Δε_eq^c
  
  // 计算当前的有效应力
  stress_dev = stress.deviatoric();
  ADReal effective_stress = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 计算蠕变应变增量
  ADReal delta_ec_total = _ec[_qp] - _ec_old[_qp];
  
  // 获取旧的蠕变能量密度
  ADReal creep_energy_old = _psic_active_old[_qp];
  
  // 蠕变能量密度的增量计算
  // 对于小增量，使用简化的矩形积分：W_c += σ_eq * Δε_eq^c
  // 对于更精确的计算，可以使用梯形积分或其他数值积分方法
  ADReal creep_energy_increment = effective_stress * delta_ec_total;
  
  // 确保蠕变能量密度非负（物理要求）
  ADReal creep_energy = creep_energy_old + std::max(0.0, raw_value(creep_energy_increment));
  
  // 设置蠕变能量密度
  _psic_active[_qp] = creep_energy;
  _psic[_qp] = _gc[_qp] * _psic_active[_qp];
  _dpsic_dd[_qp] = _dgc_dd[_qp] * _psic_active[_qp];
}

ADRankTwoTensor
J2Creep::computeStressUnified(const ADRankTwoTensor & elastic_strain)
{
    // 策略3：使用考虑退化的弹性模型
    return _elasticity_model->computeStress(elastic_strain);
}

ADReal
J2Creep::computeResidual(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 根据文档公式:
  // F1(ε_eq^(cr,np)) = ε_eq^(cr,np) - ε_eq0^cr - g(σ_eq^np, ε_eq^(cr,np)) * Δt = 0
  // 其中 σ_eq^np = σ_eq^el - C * (ε_eq^(cr,np) - ε_eq0^cr)
  
  ADReal ec_current = _ec_old[_qp] + delta_ec;
  // 计算非塑性有效应力
  ADReal effective_stress_np = effective_trial_stress - delta_ec*_three_shear_modulus;
  // 计算蠕变率
  ADReal creep_rate = computeCreepRate(effective_stress_np, ec_current);
  // 返回残差
  ADReal residual = creep_rate * _dt - delta_ec;
  
  return residual;
}

ADReal
J2Creep::computeDerivative(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 残差对 delta_ec 的导数
  ADReal ec_current = _ec_old[_qp] + delta_ec;
  // 计算非塑性有效应力
  ADReal effective_stress_np = effective_trial_stress - delta_ec*_three_shear_modulus;
  

  
  // 计算蠕变率的导数
  ADReal dg_dsigma = computeCreepRateStressDerivative(effective_stress_np, ec_current);
  ADReal dg_dec = computeCreepRateStrainDerivative(effective_stress_np, ec_current);
  
  // 返回导数: ∂F/∂(δεc) = 1 - Δt * (∂g/∂σ * (-C) + ∂g/∂εc)
  return  _dt * (dg_dsigma * (-_three_shear_modulus) + dg_dec) - 1.0;
}

Real
J2Creep::computeReferenceResidual(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 计算参考残差（用于收敛判断）
  ADReal residual = computeResidual(effective_trial_stress, delta_ec);
  return raw_value(residual);
}

// 蠕变率计算的默认实现（可以被子类覆盖）
ADReal
J2Creep::computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 默认实现：简单的幂律蠕变 g(σ) = A * σ^n
  // 子类应该覆盖这个函数以实现具体的蠕变模型
  if (effective_stress <= 1e-15)
    return 0.0;
  
  // 临时使用简单的线性关系
  ADReal A = 1.0e-6; // 蠕变系数，应该从输入参数获取
  return A * effective_stress;
}

ADReal
J2Creep::computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 蠕变率对应力的导数
  // 默认实现：对应简单线性关系的导数
  if (effective_stress <= 1e-15)
    return 0.0;
  
  ADReal A = 1.0e-6; // 蠕变系数
  return A;
}

ADReal
J2Creep::computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 蠕变率对蠕变应变的导数
  // 默认实现：对于不依赖蠕变应变的模型，导数为0
  return 0.0;
}