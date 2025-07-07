//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "J2Creep_P.h"
#include "PlasticityModel.h"
#include "ElasticityModel.h"
registerMooseObject("reproductionApp", J2Creep_P);

InputParameters
J2Creep_P::validParams()
{
  InputParameters params = CreepModel::validParams();
  params.addClassDescription("J2 creep model with plasticity check. Implements the algorithm "
                             "for separated creep and plasticity with radial return mapping.");
  params.addRequiredParam<MaterialName>("hardening_model", "Name of the plastic hardening model for yield stress calculation");
  return params;
}

J2Creep_P::J2Creep_P(const InputParameters & parameters) : CreepModel(parameters)
{
}

void
J2Creep_P::initialSetup()
{
  _hardening_model = dynamic_cast<PlasticHardeningModel *>(&getMaterial("hardening_model"));
  if (!_hardening_model)
    paramError("hardening_model",
               "Plastic hardening model " + getParam<MaterialName>("hardening_model") +
                   " is not compatible with " + name());
}

void
J2Creep_P::setQp(unsigned int qp)
{
  // 调用基类的setQp方法
  CreepModel::setQp(qp);
  
  // 设置硬化模型的qp
  if (_hardening_model)
    _hardening_model->setQp(qp);
}

void
J2Creep_P::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain)
{
  // 设置塑性模型的蠕变模型指针
  if (_plasticity_model)
    std::cout << "有塑性模型" << std::endl;

    _plasticity_model->setCreepModel(this);
    
  // Step 1: 计算弹性试验状态（移除旧的蠕变应变和塑性应变）
  ADRankTwoTensor elastic_trial_strain = elastic_strain;
  elastic_trial_strain -= _creep_strain_old[_qp];  // 减去旧的蠕变应变张量
  
  // 如果有塑性模型，也要减去旧的塑性应变张量
  if (_plasticity_model)
    elastic_trial_strain -= _plasticity_model->getPlasticStrainOld();
  
  // Step 2: 计算非塑性试验状态（假设没有塑性增量，只考虑蠕变）
  ADReal delta_ec = 0.0;
  
  // 计算试验应力
  ADRankTwoTensor stress_trial = _elasticity_model->computeStress(elastic_trial_strain);
  
  // 计算偏应力和有效应力
  ADRankTwoTensor stress_dev = stress_trial.deviatoric();
  ADReal effective_stress_trial = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 计算流动方向
  if (effective_stress_trial > 1e-12)
    _Nc[_qp] = 1.5 * stress_dev / effective_stress_trial;
  else
    _Nc[_qp].zero();
  
  // 求解非塑性蠕变增量（使用径向返回映射）
  if (effective_stress_trial > 1e-12)
  {
    returnMappingSolve(effective_stress_trial, delta_ec, _console);
  }
  
  // 更新蠕变应变和有效蠕变应变
  ADRankTwoTensor creep_strain_np = _creep_strain_old[_qp] + delta_ec * _Nc[_qp];
  ADReal ec_np = _ec_old[_qp] + delta_ec;
  
  // 计算非塑性状态下的弹性应变和应力
  ADRankTwoTensor elastic_strain_np = elastic_trial_strain - delta_ec * _Nc[_qp];
  ADRankTwoTensor stress_np = _elasticity_model->computeStress(elastic_strain_np);
  
  // 重新计算偏应力和有效应力
  stress_dev = stress_np.deviatoric();
  ADReal effective_stress_np = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
  
  // 更新流动方向
  if (effective_stress_np > 1e-12)
    _Nc[_qp] = 1.5 * stress_dev / effective_stress_np;
  else
    _Nc[_qp].zero();
  
  // Step 3: 检查屈服条件
  if (f_no_plastic_strain(effective_stress_np))
  {
    // 不屈服，只有蠕变
    _creep_strain[_qp] = creep_strain_np;
    _ec[_qp] = ec_np;
    elastic_strain = elastic_strain_np;
    stress = stress_np;
  }  
  else
  {
    // 屈服，需要同时考虑蠕变和塑性
    // 首先恢复到弹性试验状态
    elastic_strain = elastic_trial_strain;
    stress = stress_trial;
    
    // 调用塑性模型处理塑性变形
    // 注意：这里可能需要迭代求解，因为蠕变和塑性相互影响
    _plasticity_model->updateState(stress, elastic_strain);
    
    // 塑性模型处理后，重新计算蠕变
    // 这里简化处理，实际可能需要更复杂的迭代算法
    elastic_strain -= _creep_strain_old[_qp];
    stress = _elasticity_model->computeStress(elastic_strain);
    
    // 保持蠕变应变不变（这里可以根据需要调整）
    _creep_strain[_qp] = _creep_strain_old[_qp];
    _ec[_qp] = _ec_old[_qp];
  }
}

bool
J2Creep_P::f_no_plastic_strain(const ADReal & effective_stress)
{
  // 获取屈服应力（通过硬化模型）
  if (!_hardening_model)
    return true; // 没有硬化模型，只有蠕变
  
  // 通过硬化模型获取屈服应力
  // 使用上一步的塑性应变进行硬化计算
  Real ep_old_value = 0.0;  // 默认值
  if (_plasticity_model)
    ep_old_value = _plasticity_model->getEffectivePlasticStrainOld();
  
  ADReal yield_stress = _hardening_model->plasticEnergy(ep_old_value, 1);
  
  // 屈服函数: f = σ_eq - σ_Y
  ADReal f = effective_stress - yield_stress;
  
  // 返回 f <= 0 表示不屈服（只有蠕变）
  // 返回 f > 0 表示屈服（需要考虑塑性）
  return f <= 0.0;
}

ADReal
J2Creep_P::computeResidual(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 根据文档公式:
  // F1(ε_eq^(cr,np)) = ε_eq^(cr,np) - ε_eq0^cr - g(σ_eq^np, ε_eq^(cr,np)) * Δt = 0
  // 其中 σ_eq^np = σ_eq^el - C * (ε_eq^(cr,np) - ε_eq0^cr)
  
  ADReal ec_current = _ec_old[_qp] + delta_ec;
  
  // 计算非塑性有效应力
  ADReal C = _elasticity_model->computeStress(_Nc[_qp]).doubleContraction(_Nc[_qp]);
  ADReal effective_stress_np = effective_trial_stress - C * delta_ec;
  
  // 计算蠕变率
  ADReal creep_rate = computeCreepRate(effective_stress_np, ec_current);
  
  // 返回残差
  return delta_ec - creep_rate * _dt;
}

ADReal
J2Creep_P::computeDerivative(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 残差对 delta_ec 的导数
  ADReal ec_current = _ec_old[_qp] + delta_ec;
  
  // 计算非塑性有效应力
  ADReal C = _elasticity_model->computeStress(_Nc[_qp]).doubleContraction(_Nc[_qp]);
  ADReal effective_stress_np = effective_trial_stress - C * delta_ec;
  
  // 计算蠕变率的导数
  ADReal dg_dsigma = computeCreepRateStressDerivative(effective_stress_np, ec_current);
  ADReal dg_dec = computeCreepRateStrainDerivative(effective_stress_np, ec_current);
  
  // 返回导数: ∂F/∂(δεc) = 1 - Δt * (∂g/∂σ * (-C) + ∂g/∂εc)
  return 1.0 - _dt * (dg_dsigma * (-C) + dg_dec);
}

Real
J2Creep_P::computeReferenceResidual(const ADReal & effective_trial_stress, const ADReal & delta_ec)
{
  // 计算参考残差（用于收敛判断）
  ADReal residual = computeResidual(effective_trial_stress, delta_ec);
  return raw_value(residual);
}

// 蠕变率计算的默认实现（可以被子类覆盖）
ADReal
J2Creep_P::computeCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 默认实现：简单的幂律蠕变 g(σ) = A * σ^n
  // 子类应该覆盖这个函数以实现具体的蠕变模型
  if (effective_stress <= 1e-12)
    return 0.0;
  
  // 临时使用简单的线性关系
  ADReal A = 1.0e-6; // 蠕变系数，应该从输入参数获取
  return A * effective_stress;
}

ADReal
J2Creep_P::computeCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 蠕变率对应力的导数
  // 默认实现：对应简单线性关系的导数
  if (effective_stress <= 1e-12)
    return 0.0;
  
  ADReal A = 1.0e-6; // 蠕变系数
  return A;
}

ADReal
J2Creep_P::computeCreepRateStrainDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  // 蠕变率对蠕变应变的导数
  // 默认实现：对于不依赖蠕变应变的模型，导数为0
  return 0.0;
}

ADReal
J2Creep_P::getCreepRate(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  return computeCreepRate(effective_stress, effective_creep_strain);
}

ADReal
J2Creep_P::getCreepRateStressDerivative(const ADReal & effective_stress, const ADReal & effective_creep_strain)
{
  return computeCreepRateStressDerivative(effective_stress, effective_creep_strain);
}

ADReal
J2Creep_P::getEffectiveCreepStrain() const
{
  return _ec[_qp];
}
